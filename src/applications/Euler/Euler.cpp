#include "applications/Euler/Euler.hpp"

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/PatchDataRestartManager.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/mesh/TreeLoadBalancer.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/pdat/CellIterator.h"
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/pdat/FaceIndex.h"
#include "SAMRAI/pdat/FaceVariable.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Utilities.h"

#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>

#ifndef LACKS_SSTREAM
#ifndef included_sstream
#define included_sstream
#include <sstream>
#endif
#else
#ifndef included_strstream
#define included_strstream
#include <strstream>
#endif
#endif

#define EPSILON 1e-40

boost::shared_ptr<tbox::Timer> Euler::t_init = NULL;
boost::shared_ptr<tbox::Timer> Euler::t_compute_dt = NULL;
boost::shared_ptr<tbox::Timer> Euler::t_compute_hyperbolicfluxes = NULL;
boost::shared_ptr<tbox::Timer> Euler::t_advance_steps = NULL;
boost::shared_ptr<tbox::Timer> Euler::t_synchronize_hyperbloicfluxes = NULL;
boost::shared_ptr<tbox::Timer> Euler::t_setphysbcs = NULL;
boost::shared_ptr<tbox::Timer> Euler::t_taggradient = NULL;

Euler::Euler(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const boost::shared_ptr<tbox::Database>& input_db,
    const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry):
        RungeKuttaPatchStrategy(),
        d_object_name(object_name),
        d_dim(dim),
        d_grid_geometry(grid_geometry),
#ifdef HAVE_HDF5
        d_visit_writer(NULL),
#endif
        d_workload_variable(NULL),
        d_use_nonuniform_workload(false),
        d_num_ghosts(hier::IntVector::getZero(d_dim)),
        d_equation_of_state(NULL),
        d_equation_of_state_db(NULL),
        d_conv_flux_reconstructor(NULL),
        d_shock_capturing_scheme_db(NULL),
        d_initial_conditions(NULL),
        d_Euler_boundary_conditions(NULL),
        d_Euler_boundary_conditions_db(NULL),
        d_Euler_boundary_conditions_db_is_from_restart(false),
        d_feature_driven_tagger(NULL),
        d_feature_driven_tagger_db(NULL),
        d_flow_model_manager(NULL),
        d_is_preserving_positivity(false)
{
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(input_db);
    TBOX_ASSERT(grid_geometry);
    
    tbox::RestartManager::getManager()->registerRestartItem(d_object_name, this);
    
    if (!t_init)
    {
        t_init = tbox::TimerManager::getManager()->
            getTimer("Euler::initializeDataOnPatch()");
        t_compute_dt = tbox::TimerManager::getManager()->
            getTimer("Euler::computeStableDtOnPatch()");
        t_compute_hyperbolicfluxes = tbox::TimerManager::getManager()->
            getTimer("Euler::computeHyperbolicFluxesOnPatch()");
        t_advance_steps = tbox::TimerManager::getManager()->
            getTimer("Euler::advanceSingleStep()");
        t_synchronize_hyperbloicfluxes = tbox::TimerManager::getManager()->
            getTimer("Euler::Euler::synchronizeHyperbolicFluxes()");
        t_setphysbcs = tbox::TimerManager::getManager()->
            getTimer("Euler::setPhysicalBoundaryConditions()");
        t_taggradient = tbox::TimerManager::getManager()->
            getTimer("Euler::tagGradientDetectorCells()");
    }
    
    /*
     * Initialize object with data read from given input/restart databases.
     */
    bool is_from_restart = tbox::RestartManager::getManager()->isFromRestart();
    if (is_from_restart)
    {
        getFromRestart();
    }
    getFromInput(input_db, is_from_restart);
    
    /*
     * Initialize d_flow_model_manager.
     */
    d_flow_model_manager.reset(new FlowModelManager(
        "flow model manager",
        d_dim,
        d_grid_geometry,
        d_flow_model_str,
        d_num_species));
    
    /*
     * Get the string describing the equation of state.
     */
    std::string equation_of_state_string;
    if (d_equation_of_state_db->keyExists("equation_of_state"))
    {
        equation_of_state_string = d_equation_of_state_db->getString("equation_of_state");
    }
    else if (d_equation_of_state_db->keyExists("d_equation_of_state"))
    {
        equation_of_state_string = d_equation_of_state_db->getString("d_equation_of_state");
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "No key 'equation_of_state'/'d_equation_of_state' found in data for"
                   << " Equation_of_state."
                   << std::endl);
    }
    
    /*
     * Initialize the flux.
     */
    d_convective_flux = boost::shared_ptr<pdat::FaceVariable<double> > (
        new pdat::FaceVariable<double>(dim, "convective flux", d_flow_model_manager->getNumberOfEquations()));
    
    /*
     * Initialize the source.
     */
    d_source = boost::shared_ptr<pdat::CellVariable<double> > (
        new pdat::CellVariable<double>(dim, "source", d_flow_model_manager->getNumberOfEquations()));
    
    /*
     * Get the string describing the shock capturing scheme.
     */
    std::string shock_capturing_scheme_str;
    if (d_shock_capturing_scheme_db->keyExists("shock_capturing_scheme"))
    {
        shock_capturing_scheme_str = d_shock_capturing_scheme_db->
            getString("shock_capturing_scheme");
    }
    else if (d_shock_capturing_scheme_db->keyExists("d_shock_capturing_scheme"))
    {
        shock_capturing_scheme_str = d_shock_capturing_scheme_db->
            getString("d_shock_capturing_scheme");
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "No key 'shock_capturing_scheme'/'d_shock_capturing_scheme' found in data for"
                   << " Shock_capturing_scheme."
                   << std::endl);
    }
    
    /*
     * Initialize d_equation_of_state.
     */
    d_flow_model_manager->initializeEquationOfState(
        equation_of_state_string,
        d_equation_of_state_db,
        d_equation_of_state);
    
    /*
     * Initialize d_conv_flux_reconstructor.
     */
    d_flow_model_manager->initializeConvectiveFluxReconstructor(
        shock_capturing_scheme_str,
        d_shock_capturing_scheme_db,
        d_conv_flux_reconstructor,
        d_convective_flux,
        d_source);
    
    /*
     * Initialize d_initial_conditions.
     */
    d_flow_model_manager->initializeInitialConditions(
        d_project_name,
        d_initial_conditions);
    
    /*
     * Initialize d_Euler_boundary_conditions.
     */
    d_flow_model_manager->initializeEulerBoundaryConditions(
        d_project_name,
        d_Euler_boundary_conditions_db,
        d_Euler_boundary_conditions_db_is_from_restart,
        d_Euler_boundary_conditions);
    
    /*
     * Initialize d_feature_driven_tagger.
     */
    if (d_feature_driven_tagger_db != nullptr)
    {
        d_flow_model_manager->initializeFeatureDrivenTagger(
            d_feature_driven_tagger_db,
            d_feature_driven_tagger);
    }
    
    /*
     * Get the number of ghost cells needed.
     */
    d_num_ghosts = d_flow_model_manager->
        getNumberOfGhostCells();
}


Euler::~Euler()
{
    t_init.reset();
    t_compute_dt.reset();
    t_compute_hyperbolicfluxes.reset();
    t_advance_steps.reset();
    t_synchronize_hyperbloicfluxes.reset();
    t_setphysbcs.reset();
    t_taggradient.reset();
}


void
Euler::registerModelVariables(
    RungeKuttaLevelIntegrator* integrator)
{
    TBOX_ASSERT(integrator != 0);
    
    /*
     * Register the conservative variables.
     */
    d_flow_model_manager->registerConservativeVariables(integrator);
    
    /*
     * Register the fluxes and sources.
     */
    
    integrator->registerVariable(
        d_convective_flux,
        hier::IntVector::getZero(d_dim),
        RungeKuttaLevelIntegrator::HYP_FLUX,
        d_grid_geometry,
        "CONSERVATIVE_COARSEN",
        "NO_REFINE");
    
    integrator->registerVariable(
        d_source,
        hier::IntVector::getZero(d_dim),
        RungeKuttaLevelIntegrator::SOURCE,
        d_grid_geometry,
        "NO_COARSEN",
        "NO_REFINE");
    
    /*
     * Set the plotting context.
     */
    d_flow_model_manager->setPlotContext(integrator->getPlotContext());
    
    /*
     * Register the plotting quantities.
     */
#ifdef HAVE_HDF5
    if (d_visit_writer)
    {
        d_flow_model_manager->registerPlotQuantities(
            integrator,
            d_visit_writer);
    }
    
    if (!d_visit_writer)
    {
        TBOX_WARNING(d_object_name
                     << ": registerModelVariables()\n"
                     << "VisIt data writer was not registered\n"
                     << "Consequently, no plot data will\n"
                     << "be written."
                     << std::endl);
    }
#endif
}


void
Euler::setupLoadBalancer(
    RungeKuttaLevelIntegrator* integrator,
    mesh::GriddingAlgorithm* gridding_algorithm)
{
    NULL_USE(integrator);
    
    const hier::IntVector& zero_vec = hier::IntVector::getZero(d_dim);
    
    hier::VariableDatabase* vardb = hier::VariableDatabase::getDatabase();
    hier::PatchDataRestartManager* pdrm = hier::PatchDataRestartManager::getManager();
    
    if (d_use_nonuniform_workload && gridding_algorithm)
    {
        boost::shared_ptr<mesh::TreeLoadBalancer> load_balancer(
            boost::dynamic_pointer_cast<mesh::TreeLoadBalancer, mesh::LoadBalanceStrategy>(
                gridding_algorithm->getLoadBalanceStrategy()));
        
        if (load_balancer)
        {
            d_workload_variable.reset(new pdat::CellVariable<double>(
                d_dim,
                "workload_variable",
                1));
            d_workload_data_id = vardb->registerVariableAndContext(
                d_workload_variable,
                vardb->getContext("WORKLOAD"),
                zero_vec);
            load_balancer->setWorkloadPatchDataIndex(d_workload_data_id);
            pdrm->registerPatchDataForRestart(d_workload_data_id);
        }
        else
        {
            TBOX_WARNING(d_object_name << ": "
                                       << "  Unknown load balancer used in gridding algorithm."
                                       << "  Ignoring request for nonuniform load balancing."
                                       << std::endl);
            d_use_nonuniform_workload = false;
        }
    }
    else
    {
        d_use_nonuniform_workload = false;
    }
}


void
Euler::initializeDataOnPatch(
    hier::Patch& patch,
    const double data_time,
    const bool initial_time)
{   
    t_init->start();
    
    d_initial_conditions->initializeDataOnPatch(
        patch,
        data_time,
        initial_time,
        getDataContext());
    
    if (d_use_nonuniform_workload)
    {
        if (!patch.checkAllocated(d_workload_data_id))
        {
            patch.allocatePatchData(d_workload_data_id);
        }
        
        boost::shared_ptr<pdat::CellData<double> > workload_data(
            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                patch.getPatchData(d_workload_data_id)));
        TBOX_ASSERT(workload_data);
        workload_data->fillAll(1.0);
    }

    t_init->stop();
}


double
Euler::computeStableDtOnPatch(
    hier::Patch& patch,
    const bool initial_time,
    const double dt_time)
{
    t_compute_dt->start();
    
    double stable_dt = d_flow_model_manager->
        computeStableDtOnPatch(
            patch,
            initial_time,
            dt_time,
            getDataContext());
    
    t_compute_dt->stop();
    
    return stable_dt;
}


void
Euler::computeHyperbolicFluxesAndSourcesOnPatch(
    hier::Patch& patch,
    const double time,
    const double dt,
    const int RK_step_number)
{
    NULL_USE(time);
    
    t_compute_hyperbolicfluxes->start();
    
    /*
     * Set zero for the source.
     */
    
    boost::shared_ptr<pdat::CellData<double> > source(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_source, getDataContext())));
    
    source->fillAll(0.0);
    
    /*
     * Compute the fluxes and sources.
     */
    
    d_conv_flux_reconstructor->computeConvectiveFluxesAndSources(patch,
        time,
        dt,
        RK_step_number,
        getDataContext());
    
    t_compute_hyperbolicfluxes->stop();
}


void
Euler::advanceSingleStep(
    hier::Patch& patch,
    const double time,
    const double dt,
    const std::vector<double>& alpha,
    const std::vector<double>& beta,
    const std::vector<double>& gamma,
    const std::vector<boost::shared_ptr<hier::VariableContext> >& intermediate_context)
{
    NULL_USE(time);
    NULL_USE(dt);
    
    t_advance_steps->start();
    
    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(patch_geom);
#endif
    
    const double* dx = patch_geom->getDx();
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box dummy_box = patch.getBox();
    const hier::Box interior_box = dummy_box;
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // ghost cells.
    dummy_box.grow(d_num_ghosts);
    const hier::Box ghost_box = dummy_box;
    const hier::IntVector ghostcell_dims = ghost_box.numberCells();
    
    /*
     * Create a vector of pointers to time-dependent variables for the
     * current data context (SCRATCH).
     */
    std::vector<double*> Q = d_flow_model_manager->
        createConservativeVariableVector(
            patch,
            getDataContext(),
            true);
    
    /*
     * Use alpha, beta and gamma values to update the time-dependent solution,
     * flux and source
     */
    
    boost::shared_ptr<pdat::FaceData<double> > convective_flux(
        BOOST_CAST<pdat::FaceData<double>, hier::PatchData>(
        patch.getPatchData(d_convective_flux, getDataContext())));

    boost::shared_ptr<pdat::CellData<double> > source(
        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
            patch.getPatchData(d_source, getDataContext())));
    
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(convective_flux);
    TBOX_ASSERT(source);
    
    TBOX_ASSERT(convective_flux->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
    TBOX_ASSERT(source->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
#endif
    
    int num_coeffs = static_cast<int>(alpha.size());
    
    for (int n = 0; n < num_coeffs; n++)
    {
        boost::shared_ptr<pdat::FaceData<double> > convective_flux_intermediate(
                BOOST_CAST<pdat::FaceData<double>, hier::PatchData>(
                    patch.getPatchData(d_convective_flux, intermediate_context[n])));
                
        boost::shared_ptr<pdat::CellData<double> > source_intermediate(
            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                patch.getPatchData(d_source, intermediate_context[n])));
        
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(convective_flux_intermediate);
        TBOX_ASSERT(source_intermediate);
        
        TBOX_ASSERT(convective_flux_intermediate->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
        TBOX_ASSERT(source_intermediate->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
#endif
        
        /*
         * Create a vector pointers to the time-dependent variables for the
         * current intermediate data context.
         */
        std::vector<double*> Q_intermediate = d_flow_model_manager->
            createConservativeVariableVector(
                patch,
                intermediate_context[n],
                false);
        
        if (d_dim == tbox::Dimension(1))
        {
            if (alpha[n] != 0.0)
            {
                for (int i = 0; i < interior_dims[0]; i++)
                {
                    // Compute indices of time-dependent data.
                    int idx_cell   = i + d_num_ghosts[0];
                    
                    for (int ei = 0; ei < d_flow_model_manager->getNumberOfEquations(); ei++)
                    {
                        Q[ei][idx_cell] += alpha[n]*Q_intermediate[ei][idx_cell];
                    }
                }
            }
            
            if (d_is_preserving_positivity && (n == num_coeffs - 1))
            {
                preservePositivity(Q,
                    convective_flux_intermediate,
                    source_intermediate,
                    interior_dims,
                    ghostcell_dims,
                    dx,
                    dt,
                    beta[n]);
            }
            
            if (beta[n] != 0.0)
            {
                for (int i = 0; i < interior_dims[0]; i++)
                {
                    // Compute indices of time-dependent data, flux and source.
                    int idx_cell   = i + d_num_ghosts[0];
                    int idx_source = i;
                    int idx_flux_x = i + 1;
                    
                    for (int ei = 0; ei < d_flow_model_manager->getNumberOfEquations(); ei++)
                    {
                        double* F_x_intermediate = convective_flux_intermediate->getPointer(0, ei);
                        double* S_intermediate = source_intermediate->getPointer(ei);
                        
                        Q[ei][idx_cell] += beta[n]*
                            (-(F_x_intermediate[idx_flux_x] - F_x_intermediate[idx_flux_x - 1])/dx[0] +
                            S_intermediate[idx_source]);
                    }
                }
            }
                
            if (gamma[n] != 0.0)
            {
                // Accumulate the flux in the x direction.
                for (int i = 0; i < interior_dims[0] + 1; i++)
                {
                    int idx_flux_x = i;
                    
                    for (int ei = 0; ei < d_flow_model_manager->getNumberOfEquations(); ei++)
                    {
                        double* F_x              = convective_flux->getPointer(0, ei);
                        double* F_x_intermediate = convective_flux_intermediate->getPointer(0, ei);
                        
                        F_x[idx_flux_x] += gamma[n]*F_x_intermediate[idx_flux_x];
                    }                        
                }
                
                // Accumulate the source.
                for (int i = 0; i < interior_dims[0]; i++)
                {
                    int idx_cell = i;
                    
                    for (int ei = 0; ei < d_flow_model_manager->getNumberOfEquations(); ei++)
                    {
                        double* S              = source->getPointer(ei);
                        double* S_intermediate = source_intermediate->getPointer(ei);
                        
                        S[idx_cell] += gamma[n]*S_intermediate[idx_cell];
                    }
                }
            } // if (gamma[n] != 0.0)
        } // if (d_dim == tbox::Dimension(1))
        else if (d_dim == tbox::Dimension(2))
        {
            if (alpha[n] != 0.0)
            {
                for (int j = 0; j < interior_dims[1]; j++)
                {
                    for (int i = 0; i < interior_dims[0]; i++)
                    {
                        // Compute indices of time-dependent data.
                        int idx_cell   = (i + d_num_ghosts[0]) +
                            (j + d_num_ghosts[1])*ghostcell_dims[0];
                        
                        for (int ei = 0; ei < d_flow_model_manager->getNumberOfEquations(); ei++)
                        {
                            Q[ei][idx_cell] += alpha[n]*Q_intermediate[ei][idx_cell];
                        }
                    }
                }
            }
            
            if (d_is_preserving_positivity && (n == num_coeffs - 1))
            {
                preservePositivity(Q,
                    convective_flux_intermediate,
                    source_intermediate,
                    interior_dims,
                    ghostcell_dims,
                    dx,
                    dt,
                    beta[n]);
            }
            
            if (beta[n] != 0.0)
            {
                for (int j = 0; j < interior_dims[1]; j++)
                {
                    for (int i = 0; i < interior_dims[0]; i++)
                    {
                        // Compute indices of time-dependent data, flux and source.
                        int idx_cell   = (i + d_num_ghosts[0]) +
                            (j + d_num_ghosts[1])*ghostcell_dims[0];
                        int idx_source = i + j*interior_dims[0];
                        int idx_flux_x = (i + 1) + j*(interior_dims[0] + 1);
                        int idx_flux_y = (j + 1) + i*(interior_dims[1] + 1);
                        
                        for (int ei = 0; ei < d_flow_model_manager->getNumberOfEquations(); ei++)
                        {
                            double* F_x_intermediate = convective_flux_intermediate->getPointer(0, ei);
                            double* F_y_intermediate = convective_flux_intermediate->getPointer(1, ei);
                            double* S_intermediate = source_intermediate->getPointer(ei);
                            
                            Q[ei][idx_cell] += beta[n]*
                                (-(F_x_intermediate[idx_flux_x] - F_x_intermediate[idx_flux_x - 1])/dx[0] -
                                (F_y_intermediate[idx_flux_y] - F_y_intermediate[idx_flux_y - 1])/dx[1] +
                                S_intermediate[idx_source]);
                        }
                    }
                }
            }
            
            if (gamma[n] != 0.0)
            {
                // Accumulate the flux in the x direction.
                for (int j = 0; j < interior_dims[1]; j++)
                {
                    for (int i = 0; i < interior_dims[0] + 1; i++)
                    {
                        int idx_flux_x = i + j*(interior_dims[0] + 1);
                        
                        for (int ei = 0; ei < d_flow_model_manager->getNumberOfEquations(); ei++)
                        {
                            double* F_x              = convective_flux->getPointer(0, ei);
                            double* F_x_intermediate = convective_flux_intermediate->getPointer(0, ei);
                            
                            F_x[idx_flux_x] += gamma[n]*F_x_intermediate[idx_flux_x];
                        }                        
                    }
                }
                
                // Accumulate the flux in the y direction.
                for (int i = 0; i < interior_dims[0]; i++)
                {
                    for (int j = 0; j < interior_dims[1] + 1; j++)
                    {
                        int idx_flux_y = j + i*(interior_dims[1] + 1);
                        
                        for (int ei = 0; ei < d_flow_model_manager->getNumberOfEquations(); ei++)
                        {
                            double* F_y              = convective_flux->getPointer(1, ei);
                            double* F_y_intermediate = convective_flux_intermediate->getPointer(1, ei);
                            
                            F_y[idx_flux_y] += gamma[n]*F_y_intermediate[idx_flux_y];
                        }
                    }
                }
                
                // Accumulate the source.
                for (int j = 0; j < interior_dims[1]; j++)
                {
                    for (int i = 0; i < interior_dims[0]; i++)
                    {
                        int idx_cell = i + j*interior_dims[0];
                        
                        for (int ei = 0; ei < d_flow_model_manager->getNumberOfEquations(); ei++)
                        {
                            double* S              = source->getPointer(ei);
                            double* S_intermediate = source_intermediate->getPointer(ei);
                            
                            S[idx_cell] += gamma[n]*S_intermediate[idx_cell];
                        }
                    }
                }
            } // if (gamma[n] != 0.0)
        } // if (d_dim == tbox::Dimension(2))
        else if (d_dim == tbox::Dimension(3))
        {
            if (alpha[n] != 0.0)
            {
                for (int k = 0; k < interior_dims[2]; k++)
                {
                    for (int j = 0; j < interior_dims[1]; j++)
                    {
                        for (int i = 0; i < interior_dims[0]; i++)
                        {
                            // Compute indices of time-dependent data.
                            int idx_cell   = (i + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            for (int ei = 0; ei < d_flow_model_manager->getNumberOfEquations(); ei++)
                            {
                                Q[ei][idx_cell] += alpha[n]*Q_intermediate[ei][idx_cell];
                            }
                        }
                    }
                }
            }
            
            if (d_is_preserving_positivity && (n == num_coeffs - 1))
            {
                preservePositivity(Q,
                    convective_flux_intermediate,
                    source_intermediate,
                    interior_dims,
                    ghostcell_dims,
                    dx,
                    dt,
                    beta[n]);
            }
            
            if (beta[n] != 0.0)
            {
                for (int k = 0; k < interior_dims[2]; k++)
                {
                    for (int j = 0; j < interior_dims[1]; j++)
                    {
                        for (int i = 0; i < interior_dims[0]; i++)
                        {
                            // Compute indices of time-dependent data, flux and source.
                            int idx_cell   = (i + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            int idx_source = i +
                                j*interior_dims[0] +
                                k*interior_dims[0]*interior_dims[1];
                            
                            int idx_flux_x = (i + 1) +
                                j*(interior_dims[0] + 1) +
                                k*(interior_dims[0] + 1)*interior_dims[1];
                            
                            int idx_flux_y = (j + 1) +
                                k*(interior_dims[1] + 1) +
                                i*(interior_dims[1] + 1)*interior_dims[2];
                            
                            int idx_flux_z = (k + 1) +
                                i*(interior_dims[2] + 1) +
                                j*(interior_dims[2] + 1)*interior_dims[0];
                            
                            for (int ei = 0; ei < d_flow_model_manager->getNumberOfEquations(); ei++)
                            {
                                double* F_x_intermediate = convective_flux_intermediate->getPointer(0, ei);
                                double* F_y_intermediate = convective_flux_intermediate->getPointer(1, ei);
                                double* F_z_intermediate = convective_flux_intermediate->getPointer(2, ei);
                                double* S_intermediate = source_intermediate->getPointer(ei);
                                
                                Q[ei][idx_cell] += beta[n]*
                                    (-(F_x_intermediate[idx_flux_x] - F_x_intermediate[idx_flux_x - 1])/dx[0] -
                                    (F_y_intermediate[idx_flux_y] - F_y_intermediate[idx_flux_y - 1])/dx[1] -
                                    (F_z_intermediate[idx_flux_z] - F_z_intermediate[idx_flux_z - 1])/dx[2] +
                                    S_intermediate[idx_source]);
                            }
                        }
                    }
                }
            }
            
            if (gamma[n] != 0.0)
            {
                // Accumulate the flux in the x direction.
                for (int k = 0; k < interior_dims[2]; k++)
                {
                    for (int j = 0; j < interior_dims[1]; j++)
                    {
                        for (int i = 0; i < interior_dims[0] + 1; i++)
                        {
                            int idx_flux_x = i +
                                j*(interior_dims[0] + 1) +
                                k*(interior_dims[0] + 1)*interior_dims[1];
                            
                            for (int ei = 0; ei < d_flow_model_manager->getNumberOfEquations(); ei++)
                            {
                                double* F_x              = convective_flux->getPointer(0, ei);
                                double* F_x_intermediate = convective_flux_intermediate->getPointer(0, ei);
                                
                                F_x[idx_flux_x] += gamma[n]*F_x_intermediate[idx_flux_x];
                            }                        
                        }
                    }
                }
                
                // Accumulate the flux in the y direction.
                for (int i = 0; i < interior_dims[0]; i++)
                {
                    for (int k = 0; k < interior_dims[2]; k++)
                    {
                        for (int j = 0; j < interior_dims[1] + 1; j++)
                        {
                            int idx_flux_y = j +
                                k*(interior_dims[1] + 1) +
                                i*(interior_dims[1] + 1)*interior_dims[2];
                            
                            for (int ei = 0; ei < d_flow_model_manager->getNumberOfEquations(); ei++)
                            {
                                double* F_y              = convective_flux->getPointer(1, ei);
                                double* F_y_intermediate = convective_flux_intermediate->getPointer(1, ei);
                                
                                F_y[idx_flux_y] += gamma[n]*F_y_intermediate[idx_flux_y];
                            }
                        }
                    }
                }
                
                // Accumulate the flux in the z direction.
                for (int j = 0; j < interior_dims[1]; j++)
                {
                    for (int i = 0; i < interior_dims[0]; i++)
                    {
                        for (int k = 0; k < interior_dims[2]; k++)
                        {
                            int idx_flux_z = k +
                                i*(interior_dims[2] + 1) +
                                j*(interior_dims[2] + 1)*interior_dims[0];
                            
                            for (int ei = 0; ei < d_flow_model_manager->getNumberOfEquations(); ei++)
                            {
                                double* F_z              = convective_flux->getPointer(2, ei);
                                double* F_z_intermediate = convective_flux_intermediate->getPointer(2, ei);
                                
                                F_z[idx_flux_z] += gamma[n]*F_z_intermediate[idx_flux_z];
                            }
                        }
                    }
                }
                
                // Accumulate the source.
                for (int k = 0; k < interior_dims[2]; k++)
                {
                    for (int j = 0; j < interior_dims[1]; j++)
                    {
                        for (int i = 0; i < interior_dims[0]; i++)
                        {
                            int idx_cell = i +
                                j*interior_dims[0] +
                                k*interior_dims[0]*interior_dims[1];
                            
                            for (int ei = 0; ei < d_flow_model_manager->getNumberOfEquations(); ei++)
                            {
                                double* S              = source->getPointer(ei);
                                double* S_intermediate = source_intermediate->getPointer(ei);
                                
                                S[idx_cell] += gamma[n]*S_intermediate[idx_cell];
                            }
                        }
                    }
                }
            } // if (gamma[n] != 0.0)
        } // if (d_dim == tbox::Dimension(3))
        
        // Update the mass fraction/volume fraction of the last species.
        if (beta[n] != 0.0)
        {
            d_flow_model_manager->
                updateConservativeVariableVector(
                    patch,
                    Q);
        }
    }
    
    t_advance_steps->stop();
}


void
Euler::synchronizeHyperbolicFluxes(
    hier::Patch& patch,
    const double time,
    const double dt)
{
    NULL_USE(time);
    NULL_USE(dt);
    
    t_synchronize_hyperbloicfluxes->start();
    
    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(patch_geom);
#endif
    
    const double* dx = patch_geom->getDx();
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box dummy_box = patch.getBox();
    const hier::Box interior_box = dummy_box;
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // ghost cells.
    dummy_box.grow(d_num_ghosts);
    const hier::Box ghost_box = dummy_box;
    const hier::IntVector ghostcell_dims = ghost_box.numberCells();
    
    /*
     * Create a vector of pointers to time-dependent variables for the
     * current data context (SCRATCH).
     */
    std::vector<double*> Q = d_flow_model_manager->
        createConservativeVariableVector(
            patch,
            getDataContext(),
            false);
    
    boost::shared_ptr<pdat::FaceData<double> > convective_flux(
    BOOST_CAST<pdat::FaceData<double>, hier::PatchData>(
        patch.getPatchData(d_convective_flux, getDataContext())));

    boost::shared_ptr<pdat::CellData<double> > source(
        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
            patch.getPatchData(d_source, getDataContext())));
    
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(convective_flux);
    TBOX_ASSERT(source);
    
    TBOX_ASSERT(convective_flux->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
    TBOX_ASSERT(source->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
#endif
    
    if (d_dim == tbox::Dimension(1))
    {
        for (int i = 0; i < interior_dims[0]; i++)
        {
            // Compute indices of time-dependent variables, flux and source.
            int idx_cell   = i + d_num_ghosts[0];
            int idx_source = i;
            int idx_flux_x = i + 1;
            
            for (int ei = 0; ei < d_flow_model_manager->getNumberOfEquations(); ei++)
            {
                double *F_x = convective_flux->getPointer(0, ei);
                double *S   = source->getPointer(ei);
                
                Q[ei][idx_cell] +=
                    (-(F_x[idx_flux_x] - F_x[idx_flux_x - 1])/dx[0] +
                    S[idx_source]);
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        for (int j = 0; j < interior_dims[1]; j++)
        {
            for (int i = 0; i < interior_dims[0]; i++)
            {
                // Compute indices of time-dependent variables, flux and source.
                int idx_cell   = (i + d_num_ghosts[0]) + (j + d_num_ghosts[1])*ghostcell_dims[0];
                int idx_source = i + j*interior_dims[0];
                int idx_flux_x = (i + 1) + j*(interior_dims[0] + 1);
                int idx_flux_y = (j + 1) + i*(interior_dims[1] + 1);
                
                for (int ei = 0; ei < d_flow_model_manager->getNumberOfEquations(); ei++)
                {
                    double *F_x = convective_flux->getPointer(0, ei);
                    double *F_y = convective_flux->getPointer(1, ei);
                    double *S   = source->getPointer(ei);
                    
                    Q[ei][idx_cell] +=
                        (-(F_x[idx_flux_x] - F_x[idx_flux_x - 1])/dx[0] -
                        (F_y[idx_flux_y] - F_y[idx_flux_y - 1])/dx[1] +
                        S[idx_source]);
                }
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        for (int k = 0; k < interior_dims[2]; k++)
        {
            for (int j = 0; j < interior_dims[1]; j++)
            {
                for (int i = 0; i < interior_dims[0]; i++)
                {
                    // Compute indices of time-dependent variables, flux and source.
                    int idx_cell   = (i + d_num_ghosts[0]) +
                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                    
                    int idx_source = i +
                        j*interior_dims[0] +
                        k*interior_dims[0]*interior_dims[1];
                    
                    int idx_flux_x = (i + 1) +
                        j*(interior_dims[0] + 1) +
                        k*(interior_dims[0] + 1)*interior_dims[1];
                    
                    int idx_flux_y = (j + 1) +
                        k*(interior_dims[1] + 1) +
                        i*(interior_dims[1] + 1)*interior_dims[2];
                    
                    int idx_flux_z = (k + 1) +
                        i*(interior_dims[2] + 1) +
                        j*(interior_dims[2] + 1)*interior_dims[0];
                    
                    for (int ei = 0; ei < d_flow_model_manager->getNumberOfEquations(); ei++)
                    {
                        double *F_x = convective_flux->getPointer(0, ei);
                        double *F_y = convective_flux->getPointer(1, ei);
                        double *F_z = convective_flux->getPointer(2, ei);
                        double *S   = source->getPointer(ei);
                        
                        Q[ei][idx_cell] +=
                            (-(F_x[idx_flux_x] - F_x[idx_flux_x - 1])/dx[0] -
                            (F_y[idx_flux_y] - F_y[idx_flux_y - 1])/dx[1] -
                            (F_z[idx_flux_z] - F_z[idx_flux_z - 1])/dx[2] +
                            S[idx_source]);
                    }
                }
            }
        }
    }
    
    // Update the mass fraction/volume fraction of the last species.
    d_flow_model_manager->updateConservativeVariableVector(
        patch,
        Q);
    
    t_synchronize_hyperbloicfluxes->stop();
}


void
Euler::tagGradientDetectorCells(
    hier::Patch& patch,
    const double regrid_time,
    const bool initial_error,
    const int tag_indx,
    const bool uses_richardson_extrapolation_too)
{
    t_taggradient->start();
    
    // Get the tags.
    boost::shared_ptr<pdat::CellData<int> > tags(
        BOOST_CAST<pdat::CellData<int>, hier::PatchData>(
            patch.getPatchData(tag_indx)));
    
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(tags);
    TBOX_ASSERT(tags->getGhostCellWidth() == 0);
#endif
    
    // Initialize values of all tags to zero.
    tags->fillAll(0.0);
    
    // Tag the cells by using d_feature_driven_tagger.
    d_feature_driven_tagger->tagCells(
        patch,
        regrid_time,
        initial_error,
        uses_richardson_extrapolation_too,
        tags,
        getDataContext());
    
    t_taggradient->stop();
}


void
Euler::setPhysicalBoundaryConditions(
    hier::Patch& patch,
    const double fill_time,
    const hier::IntVector& ghost_width_to_fill)
{
    t_setphysbcs->start();
    
    d_Euler_boundary_conditions->setPhysicalBoundaryConditions(
        patch,
        fill_time,
        ghost_width_to_fill,
        getDataContext());
    
    t_setphysbcs->stop();
}


void
Euler::putToRestart(
    const boost::shared_ptr<tbox::Database>& restart_db) const
{
    TBOX_ASSERT(restart_db);
    
    restart_db->putString("d_project_name", d_project_name);
    
    restart_db->putInteger("d_num_species", d_num_species);
    
    // Put the d_flow_model into the restart database.
    d_flow_model_manager->putToRestart(restart_db);
    
    restart_db->putBool("d_is_preserving_positivity", d_is_preserving_positivity);
    
    boost::shared_ptr<tbox::Database> restart_equation_of_state_db =
        restart_db->putDatabase("Equation_of_state");
    
    d_equation_of_state->putToRestart(restart_equation_of_state_db);
    
    boost::shared_ptr<tbox::Database> restart_shock_capturing_scheme_db =
        restart_db->putDatabase("Shock_capturing_scheme");
    
    d_conv_flux_reconstructor->putToRestart(restart_shock_capturing_scheme_db);
    
    restart_db->putIntegerArray("d_num_ghosts", &d_num_ghosts[0], d_dim.getValue());
    
    boost::shared_ptr<tbox::Database> restart_Euler_boundary_conditions_db =
        restart_db->putDatabase("Boundary_data");
    
    d_Euler_boundary_conditions->putToRestart(restart_Euler_boundary_conditions_db);
    
    boost::shared_ptr<tbox::Database> restart_feature_driven_tagger_db =
        restart_db->putDatabase("Feature_driven_tagger");
    
    d_feature_driven_tagger->putToRestart(restart_feature_driven_tagger_db);
}


#ifdef HAVE_HDF5
void
Euler::registerVisItDataWriter(
    const boost::shared_ptr<appu::VisItDataWriter>& viz_writer)
{
    TBOX_ASSERT(viz_writer);
    d_visit_writer = viz_writer;
}
#endif


bool
Euler::packDerivedDataIntoDoubleBuffer(
    double* buffer,
    const hier::Patch& patch,
    const hier::Box& region,
    const std::string& variable_name,
    int depth_id,
    double simulation_time) const
{
    bool data_on_patch = d_flow_model_manager->
        packDerivedDataIntoDoubleBuffer(
                buffer,
                patch,
                region,
                variable_name,
                depth_id,
                simulation_time);
    
    return data_on_patch;
}


void Euler::printClassData(std::ostream& os) const
{
    os << "\nPrint Euler object... " << std::endl;
    os << std::endl;
    
    os << "Euler: this = " << (Euler *)this << std::endl;
    os << "d_object_name = " << d_object_name << std::endl;
    os << "d_project_name = " << d_project_name << std::endl;
    os << "d_dim = " << d_dim.getValue() << std::endl;
    os << "d_grid_geometry = " << d_grid_geometry.get() << std::endl;
    os << "d_num_ghosts = " << d_num_ghosts << std::endl;
    
    // Print all characteristics of d_flow_model.
    d_flow_model_manager->printClassData(os);
    
    os << "d_num_species = " << d_num_species << std::endl;
    os << "d_is_preserving_positivity = " << d_is_preserving_positivity << std::endl;
    
    os << "--------------------------------------------------------------------------------";
    
    /*
     * Print data of d_grid_geometry.
     */
    os << "\nPrint CartesianGridGeometry object..." << std::endl;
    os << std::endl;
    
    os << "CartesianGridGeometry: this = "
       << (geom::CartesianGridGeometry *) &d_grid_geometry
       << std::endl;
    
    os << "d_object_name = "
       << d_grid_geometry->getObjectName()
       << std::endl;
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_up = d_grid_geometry->getXUpper();
    
    os << "d_x_lo = (" << x_lo[0];
    for (int di = 1; di < d_dim.getValue(); di++)
    {
        os << "," << x_lo[di];
    }
    os << ")" << std::endl;
    
    os << "d_x_up = (" << x_up[0];
    for (int di = 1; di < d_dim.getValue(); di++)
    {
        os << "," << x_up[di];
    }
    os << ")" << std::endl;
    
    const double* dx = d_grid_geometry->getDx();
    
    os << "d_dx = (" << dx[0];
    for (int di = 1; di < d_dim.getValue(); di++)
    {
        os << "," << dx[di];
    }
    os << ")" << std::endl;
    
    const hier::BoxContainer domain_box = d_grid_geometry->getPhysicalDomain();
    
    os << "d_domain_box = ";
    domain_box.print(os);
    
    os << "d_periodic_shift = " << d_grid_geometry->getPeriodicShift(hier::IntVector::getOne(d_dim));
    os << std::endl;
    
    os << "--------------------------------------------------------------------------------";
    
    /*
     * Print data of d_equation_of_state object
     */
    
    d_equation_of_state->printClassData(os);
    os << "--------------------------------------------------------------------------------";
    
    /*
     * Print data of d_conv_flux_reconstructor.
     */
    
    d_conv_flux_reconstructor->printClassData(os);
    os << "--------------------------------------------------------------------------------";
    
    /*
     * Print Refinement data
     */
    
    /*
     * Print data of d_Euler_boundary_conditions.
     */
    
    d_Euler_boundary_conditions->printClassData(os);
}


void
Euler::printDataStatistics(
    std::ostream& os,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy) const
{
    d_flow_model_manager->
        printDataStatistics(
            os,
            patch_hierarchy);
}


void
Euler::getFromInput(
    const boost::shared_ptr<tbox::Database>& input_db,
    bool is_from_restart)
{
    /*
     * Note: if we are restarting, then we only allow nonuniform
     * workload to be used if nonuniform workload was used originally.
     */
    if (!is_from_restart)
    {
        d_use_nonuniform_workload = input_db->
            getBoolWithDefault(
                "use_nonuniform_workload",
                d_use_nonuniform_workload);
    }
    else
    {
        if (d_use_nonuniform_workload)
        {
            d_use_nonuniform_workload = input_db->getBool("use_nonuniform_workload");
        }
    }
    
    if (!is_from_restart)
    {
        if (input_db->keyExists("project_name"))
        {
            d_project_name = input_db->getString("project_name");
        }
        else
        {
            d_project_name = "Unnamed";
        }
        
        if (input_db->keyExists("num_species"))
        {
            d_num_species = input_db->getInteger("num_species");
            
            if (d_num_species <= 0)
            {
                TBOX_ERROR(d_object_name
                           << ": "
                           << "Non-positive number of species is specified."
                           << " Number of species should be positive."
                           << std::endl);
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "Key data 'num_species' not found in input database."
                       << " Number of species is unknown."
                       << std::endl);
        }
        
        /*
         * Initialize the flow model.
         */
        if (input_db->keyExists("flow_model"))
        {
            d_flow_model_str = input_db->getString("flow_model");
        }
        else
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "Key data 'flow model' not found in input database."
                       << " Compressible flow model is unknown."
                       << std::endl);            
        }
        
        if (input_db->keyExists("preserving_positivity"))
        {
            d_is_preserving_positivity =
                input_db-> getBool("preserving_positivity");
        }
        else
        {
            d_is_preserving_positivity = false;
        }
        
        /*
         * Get the database of the equation of state.
         */
        if (input_db->keyExists("Equation_of_state"))
        {
            d_equation_of_state_db =
                input_db->getDatabase("Equation_of_state");
        }
        else
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "Key data 'Equation_of_state' not found in input database."
                       << std::endl);
        }
        
        /*
         * Get the database of the convective flux reconstructor.
         */
        if (input_db->keyExists("Shock_capturing_scheme"))
        {
            d_shock_capturing_scheme_db = input_db->getDatabase("Shock_capturing_scheme");
        }
        else
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "Key data 'Shock_capturing_scheme' not found in input database."
                       << std::endl);
        }
        
        if (input_db->keyExists("Feature_driven_tagger"))
        {
            d_feature_driven_tagger_db =
                input_db->getDatabase("Feature_driven_tagger");
        }
    }
    
    /*
     * Get the boundary conditions from the input database.
     */
    
    const hier::IntVector &one_vec = hier::IntVector::getOne(d_dim);
    hier::IntVector periodic = d_grid_geometry->getPeriodicShift(one_vec);
    int num_per_dirs = 0;
    for (int di = 0; di < d_dim.getValue(); di++)
    {
        if (periodic(di))
        {
            num_per_dirs++;
        }
    }
    
    if (num_per_dirs < d_dim.getValue())
    {
        if (input_db->keyExists("Boundary_data"))
        {
            d_Euler_boundary_conditions_db = input_db->getDatabase(
                "Boundary_data");
            
            d_Euler_boundary_conditions_db_is_from_restart = false;
        }
        else
        {
            if (!is_from_restart)
            {
                TBOX_ERROR(d_object_name
                           << ": "
                           << "Key data 'Boundary_data' not found in input database."
                           << std::endl);
            }
        }
    }
}


void Euler::getFromRestart()
{
    boost::shared_ptr<tbox::Database> root_db(tbox::RestartManager::getManager()->getRootDatabase());
    
    if (!root_db->isDatabase(d_object_name))
    {
        TBOX_ERROR("Restart database corresponding to "
                   << d_object_name
                   << " not found in restart file."
                   << std::endl);
    }
    
    boost::shared_ptr<tbox::Database> db(root_db->getDatabase(d_object_name));
    
    d_project_name = db->getString("d_project_name");
    
    d_num_species = db->getInteger("d_num_species");
    
    d_flow_model_str = db->getString("d_flow_model");
    
    d_is_preserving_positivity = db->getBool("d_is_preserving_positivity");
    
    d_equation_of_state_db = db->getDatabase("Equation_of_state");
    
    d_shock_capturing_scheme_db = db->getDatabase("Shock_capturing_scheme");
    
    int* tmp_num_ghosts = &d_num_ghosts[0];
    db->getIntegerArray("d_num_ghosts", tmp_num_ghosts, d_dim.getValue());
    
    d_Euler_boundary_conditions_db = db->getDatabase("Boundary_data");
    
    d_Euler_boundary_conditions_db_is_from_restart = true;
    
    if (db->keyExists("Feature_driven_tagger"))
    {
        d_feature_driven_tagger_db = db->getDatabase("Feature_driven_tagger");
    }
}


void
Euler::preservePositivity(
    std::vector<double*>& Q,
    const boost::shared_ptr<pdat::FaceData<double> >& convective_flux_intermediate,
    const boost::shared_ptr<pdat::CellData<double> >& source_intermediate,
    const hier::IntVector interior_dims,
    const hier::IntVector ghostcell_dims,
    const double* const dx,
    const double& dt,
    const double& beta)
{
    NULL_USE(Q);
    NULL_USE(convective_flux_intermediate);
    NULL_USE(source_intermediate);
    NULL_USE(interior_dims);
    NULL_USE(ghostcell_dims);
    NULL_USE(dx);
    NULL_USE(dt);
    NULL_USE(beta);
}
