#include "apps/Euler/Euler.hpp"

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
#include "SAMRAI/pdat/SideIndex.h"
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
#include <limits>
#include <fstream>

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

boost::shared_ptr<tbox::Timer> Euler::t_init;
boost::shared_ptr<tbox::Timer> Euler::t_compute_dt;
boost::shared_ptr<tbox::Timer> Euler::t_compute_fluxes_sources;
boost::shared_ptr<tbox::Timer> Euler::t_advance_step;
boost::shared_ptr<tbox::Timer> Euler::t_synchronize_fluxes;
boost::shared_ptr<tbox::Timer> Euler::t_setphysbcs;
boost::shared_ptr<tbox::Timer> Euler::t_tagvalue;
boost::shared_ptr<tbox::Timer> Euler::t_taggradient;
boost::shared_ptr<tbox::Timer> Euler::t_tagmultiresolution;

Euler::Euler(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const boost::shared_ptr<tbox::Database>& input_db,
    const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
    const std::string& stat_dump_filename):
        RungeKuttaPatchStrategy(),
        d_object_name(object_name),
        d_dim(dim),
        d_grid_geometry(grid_geometry),
        d_stat_dump_filename(stat_dump_filename),
        d_use_nonuniform_workload(false),
        d_Euler_boundary_conditions_db_is_from_restart(false)
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
        t_compute_fluxes_sources = tbox::TimerManager::getManager()->
            getTimer("Euler::computeHyperbolicFluxesOnPatch()");
        t_advance_step = tbox::TimerManager::getManager()->
            getTimer("Euler::advanceSingleStepOnPatch()");
        t_synchronize_fluxes = tbox::TimerManager::getManager()->
            getTimer("Euler::synchronizeHyperbolicFluxes()");
        t_setphysbcs = tbox::TimerManager::getManager()->
            getTimer("Euler::setPhysicalBoundaryConditions()");
        t_tagvalue = tbox::TimerManager::getManager()->
            getTimer("Euler::tagValueDetectorCells()");
        t_taggradient = tbox::TimerManager::getManager()->
            getTimer("Euler::tagGradientDetectorCells()");
        t_tagmultiresolution = tbox::TimerManager::getManager()->
            getTimer("Euler::tagMultiresolutionDetectorCells()");
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
     * Initialize d_flow_model_manager and get the flow model object.
     */
    
    d_flow_model_manager.reset(new FlowModelManager(
        "d_flow_model_manager",
        d_dim,
        d_grid_geometry,
        d_num_species,
        d_flow_model_db,
        d_flow_model_str));
    
    d_flow_model = d_flow_model_manager->getFlowModel();
    
    /*
     * Initialize d_convective_flux_reconstructor_manager and get the convective flux reconstructor object.
     */
    
    d_convective_flux_reconstructor_manager.reset(new ConvectiveFluxReconstructorManager(
        "d_convective_flux_reconstructor_manager",
        d_dim,
        d_grid_geometry,
        d_flow_model->getNumberOfEquations(),
        d_num_species,
        d_flow_model,
        d_convective_flux_reconstructor_db,
        d_convective_flux_reconstructor_str));
    
    d_convective_flux_reconstructor = d_convective_flux_reconstructor_manager->getConvectiveFluxReconstructor();
    
    /*
     * Initialize d_Euler_initial_conditions.
     */
    
    d_Euler_initial_conditions.reset(new EulerInitialConditions(
        "d_Euler_initial_conditions",
        d_project_name,
        d_dim,
        d_grid_geometry,
        d_flow_model_manager->getFlowModelType(),
        d_num_species));
    
    /*
     * Initialize d_Euler_boundary_conditions.
     */
    
    d_Euler_boundary_conditions.reset(new EulerBoundaryConditions(
        "d_Euler_boundary_conditions",
        d_project_name,
        d_dim,
        d_grid_geometry,
        d_num_species,
        d_flow_model_manager->getFlowModelType(),
        d_flow_model,
        d_Euler_boundary_conditions_db,
        d_Euler_boundary_conditions_db_is_from_restart));
    
    /*
     * Initialize d_value_tagger.
     */
    
    if (d_value_tagger_db != nullptr)
    {
        d_value_tagger.reset(new ValueTagger(
            "d_value_tagger",
            d_dim,
            d_grid_geometry,
            d_num_species,
            d_flow_model,
            d_value_tagger_db));
    }
    
    /*
     * Initialize d_gradient_tagger.
     */
    
    if (d_gradient_tagger_db != nullptr)
    {
        d_gradient_tagger.reset(new GradientTagger(
            "d_gradient_tagger",
            d_dim,
            d_grid_geometry,
            d_num_species,
            d_flow_model,
            d_gradient_tagger_db));
    }
    
    /*
     * Initialize d_multiresolution_tagger.
     */
    
    if (d_multiresolution_tagger_db != nullptr)
    {
        d_multiresolution_tagger.reset(new MultiresolutionTagger(
            "d_multiresolution_tagger",
            d_dim,
            d_grid_geometry,
            d_num_species,
            d_flow_model,
            d_multiresolution_tagger_db));
    }
    
    /*
     * Initialize the side variable of convective flux.
     */
    
    d_variable_convective_flux = boost::shared_ptr<pdat::SideVariable<double> > (
        new pdat::SideVariable<double>(dim, "convective flux", d_flow_model->getNumberOfEquations()));
    
    /*
     * Initialize the cell variable of source.
     */
    
    d_variable_source = boost::shared_ptr<pdat::CellVariable<double> > (
        new pdat::CellVariable<double>(dim, "source", d_flow_model->getNumberOfEquations()));
    
    if ((!d_stat_dump_filename.empty()))
    {
        d_flow_model->setupStatisticsUtilities();
        
        const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
        if (mpi.getRank() == 0)
        {
            std::ofstream f_out;
            f_out.open(d_stat_dump_filename.c_str(), std::ios::app);
            
            if (!f_out.is_open())
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "Failed to open file to output statistics!"
                    << std::endl);
            }
            
            f_out << "TIME                 ";
            f_out.close();
        }
        
        boost::shared_ptr<FlowModelStatisticsUtilities> flow_model_statistics_utilities =
            d_flow_model->getFlowModelStatisticsUtilities();
        
        flow_model_statistics_utilities->outputStatisticalQuantitiesNames(
            d_stat_dump_filename);
        
        if (mpi.getRank() == 0)
        {
            std::ofstream f_out;
            f_out.open(d_stat_dump_filename.c_str(), std::ios::app);
            
            if (!f_out.is_open())
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "Failed to open file to output statistics!"
                    << std::endl);
            }
            
            f_out << std::endl;
            f_out.close();
        }
    }
}


Euler::~Euler()
{
    t_init.reset();
    t_compute_dt.reset();
    t_compute_fluxes_sources.reset();
    t_advance_step.reset();
    t_synchronize_fluxes.reset();
    t_setphysbcs.reset();
    t_tagvalue.reset();
    t_taggradient.reset();
    t_tagmultiresolution.reset();
}


void
Euler::registerModelVariables(
    RungeKuttaLevelIntegrator* integrator)
{
    TBOX_ASSERT(integrator != 0);
    
    /*
     * Determine the numbers of ghost cells needed.
     */
    
    hier::IntVector num_ghosts_intermediate = hier::IntVector::getZero(d_dim);
    
    num_ghosts_intermediate = hier::IntVector::max(
        num_ghosts_intermediate,
        d_convective_flux_reconstructor->getConvectiveFluxNumberOfGhostCells());
    
    hier::IntVector num_ghosts = num_ghosts_intermediate;
    
    if (d_value_tagger != nullptr)
    {
        num_ghosts = hier::IntVector::max(
            num_ghosts,
            d_value_tagger->getValueTaggerNumberOfGhostCells());
    }
    
    if (d_gradient_tagger != nullptr)
    {
        num_ghosts = hier::IntVector::max(
            num_ghosts,
            d_gradient_tagger->getGradientTaggerNumberOfGhostCells());
    }
    
    if (d_multiresolution_tagger != nullptr)
    {
        num_ghosts = hier::IntVector::max(
            num_ghosts,
            d_multiresolution_tagger->getMultiresolutionTaggerNumberOfGhostCells());
    }
    
    /*
     * Register the conservative variables of d_flow_model.
     */
    d_flow_model->registerConservativeVariables(
        integrator,
        num_ghosts,
        num_ghosts_intermediate);
    
    /*
     * Register the fluxes and sources.
     */
    
    integrator->registerVariable(
        d_variable_convective_flux,
        hier::IntVector::getZero(d_dim),
        hier::IntVector::getZero(d_dim),
        RungeKuttaLevelIntegrator::FLUX,
        d_grid_geometry,
        "CONSERVATIVE_COARSEN",
        "NO_REFINE");
    
    integrator->registerVariable(
        d_variable_source,
        hier::IntVector::getZero(d_dim),
        hier::IntVector::getZero(d_dim),
        RungeKuttaLevelIntegrator::SOURCE,
        d_grid_geometry,
        "NO_COARSEN",
        "NO_REFINE");
    
    /*
     * Register the temporary variables used in refinement taggers.
     */
    
    if (d_value_tagger != nullptr)
    {
        d_value_tagger->registerValueTaggerVariables(integrator);
    }
    
    if (d_gradient_tagger != nullptr)
    {
        d_gradient_tagger->registerGradientTaggerVariables(integrator);
    }
    
    if (d_multiresolution_tagger != nullptr)
    {
        d_multiresolution_tagger->registerMultiresolutionTaggerVariables(integrator);
    }
    
    /*
     * Register the temporary variables used in statistics utilities.
     */
    
    d_flow_model->setupStatisticsUtilities();
    
    boost::shared_ptr<FlowModelStatisticsUtilities> flow_model_statistics_utilities =
        d_flow_model->getFlowModelStatisticsUtilities();
    
    flow_model_statistics_utilities->registerVariables(
        integrator,
        num_ghosts);
    
    /*
     * Set the plotting context.
     */
    setPlotContext(integrator->getPlotContext());
    
    d_flow_model->setPlotContext(integrator->getPlotContext());
    
    /*
     * Register the plotting quantities.
     */
#ifdef HAVE_HDF5
    if (d_visit_writer)
    {
        d_flow_model->registerPlotQuantities(
            d_visit_writer);
        
        if (d_value_tagger != nullptr)
        {
            d_value_tagger->registerPlotQuantities(
                d_visit_writer,
                integrator->getPlotContext());
        }
        
        if (d_gradient_tagger != nullptr)
        {
            d_gradient_tagger->registerPlotQuantities(
                d_visit_writer,
                integrator->getPlotContext());
        }
        
        if (d_multiresolution_tagger != nullptr)
        {
            d_multiresolution_tagger->registerPlotQuantities(
                d_visit_writer,
                integrator->getPlotContext());
        }
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
    
    d_flow_model->registerPatchWithDataContext(patch, getDataContext());
    
    std::vector<boost::shared_ptr<pdat::CellData<double> > > conservative_var_data =
        d_flow_model->getGlobalCellDataConservativeVariables();
    
    d_Euler_initial_conditions->initializeDataOnPatch(
        patch,
        conservative_var_data,
        data_time,
        initial_time);
    
    d_flow_model->unregisterPatch();
    
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
    
    double stable_dt;
    
    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(patch_geom);
#endif
    
    const double* dx = patch_geom->getDx();
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    double stable_spectral_radius = 0.0;
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the dimension and grid spacing.
         */
        
        const int interior_dim_0 = interior_dims[0];
        
        const double dx_0 = dx[0];
        
        /*
         * Register the patch and maximum wave speed in the flow model and compute the corresponding cell data.
         */
        
        d_flow_model->registerPatchWithDataContext(patch, getDataContext());
        
        hier::IntVector num_ghosts = d_flow_model->getNumberOfGhostCells();
        
        std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
        num_subghosts_of_data.insert(
            std::pair<std::string, hier::IntVector>(
                "MAX_WAVE_SPEED_X", num_ghosts));
        
        d_flow_model->registerDerivedCellVariable(num_subghosts_of_data);
        
        d_flow_model->computeGlobalDerivedCellData();
        
        /*
         * Get the pointer to the maximum wave speed inside the flow model.
         * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
         */
        
        boost::shared_ptr<pdat::CellData<double> > max_wave_speed_x =
            d_flow_model->getGlobalCellData("MAX_WAVE_SPEED_X");
        
        hier::IntVector num_subghosts_max_wave_speed_x = max_wave_speed_x->getGhostCellWidth();
        
        TBOX_ASSERT(num_subghosts_max_wave_speed_x == num_ghosts);
        
        const int num_ghosts_0 = num_ghosts[0];
        
        double* max_lambda_x = max_wave_speed_x->getPointer(0);
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd reduction(max:stable_spectral_radius)
#endif
        for (int i = -num_ghosts_0;
             i < interior_dim_0 + num_ghosts_0;
             i++)
        {
            // Compute the linear index.
            const int idx = i + num_ghosts_0;
            
            const double spectral_radius = max_lambda_x[idx]/dx_0;
            stable_spectral_radius = fmax(stable_spectral_radius, spectral_radius);
        }
        
        /*
         * Unregister the patch and data of all registered derived cell variables in the flow model.
         */
        
        d_flow_model->unregisterPatch();
        
    }
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the dimensions and grid spacings.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        const double dx_0 = dx[0];
        const double dx_1 = dx[1];
        
        /*
         * Register the patch and maximum wave speeds in the flow model and compute the corresponding cell data.
         */
        
        d_flow_model->registerPatchWithDataContext(patch, getDataContext());
        
        hier::IntVector num_ghosts = d_flow_model->getNumberOfGhostCells();
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts);
        const hier::IntVector ghostcell_dims = ghost_box.numberCells();
        
        std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
        num_subghosts_of_data.insert(
            std::pair<std::string, hier::IntVector>(
                "MAX_WAVE_SPEED_X", num_ghosts));
        num_subghosts_of_data.insert(
            std::pair<std::string, hier::IntVector>(
                "MAX_WAVE_SPEED_Y", num_ghosts));
        
        d_flow_model->registerDerivedCellVariable(num_subghosts_of_data);
        
        d_flow_model->computeGlobalDerivedCellData();
        
        /*
         * Get the pointers to the maximum wave speeds inside the flow model.
         * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
         */
        
        boost::shared_ptr<pdat::CellData<double> > max_wave_speed_x =
            d_flow_model->getGlobalCellData("MAX_WAVE_SPEED_X");
        
        boost::shared_ptr<pdat::CellData<double> > max_wave_speed_y =
            d_flow_model->getGlobalCellData("MAX_WAVE_SPEED_Y");
        
        hier::IntVector num_subghosts_max_wave_speed_x = max_wave_speed_x->getGhostCellWidth();
        hier::IntVector num_subghosts_max_wave_speed_y = max_wave_speed_y->getGhostCellWidth();
        
        TBOX_ASSERT(num_subghosts_max_wave_speed_x == num_ghosts);
        TBOX_ASSERT(num_subghosts_max_wave_speed_y == num_ghosts);
        
        const int num_ghosts_0 = num_ghosts[0];
        const int num_ghosts_1 = num_ghosts[1];
        const int ghostcell_dim_0 = ghostcell_dims[0];
        
        double* max_lambda_x = max_wave_speed_x->getPointer(0);
        double* max_lambda_y = max_wave_speed_y->getPointer(0);
        
        for (int j = -num_ghosts_1;
             j < interior_dim_1 + num_ghosts_1;
             j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd reduction(max:stable_spectral_radius)
#endif
            for (int i = -num_ghosts_0;
                 i < interior_dim_0 + num_ghosts_0;
                 i++)
            {
                // Compute the linear indices.
                const int idx = (i + num_ghosts_0) +
                    (j + num_ghosts_1)*ghostcell_dim_0;
                
                const double spectral_radius = max_lambda_x[idx]/dx_0 +
                    max_lambda_y[idx]/dx_1;
                
                stable_spectral_radius = fmax(stable_spectral_radius, spectral_radius);
            }
        }
        
        /*
         * Unregister the patch and data of all registered derived cell variables in the flow model.
         */
        
        d_flow_model->unregisterPatch();
        
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the dimensions and grid spacings.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        const int interior_dim_2 = interior_dims[2];
        
        const double dx_0 = dx[0];
        const double dx_1 = dx[1];
        const double dx_2 = dx[2];
        
        /*
         * Register the patch and maximum wave speeds in the flow model and compute the corresponding
         * cell data.
         */
        
        d_flow_model->registerPatchWithDataContext(patch, getDataContext());
        
        hier::IntVector num_ghosts = d_flow_model->getNumberOfGhostCells();
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts);
        const hier::IntVector ghostcell_dims = ghost_box.numberCells();
        
        std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
        num_subghosts_of_data.insert(
            std::pair<std::string, hier::IntVector>(
                "MAX_WAVE_SPEED_X", num_ghosts));
        num_subghosts_of_data.insert(
            std::pair<std::string, hier::IntVector>(
                "MAX_WAVE_SPEED_Y", num_ghosts));
        num_subghosts_of_data.insert(
            std::pair<std::string, hier::IntVector>(
                "MAX_WAVE_SPEED_Z", num_ghosts));
        
        d_flow_model->registerDerivedCellVariable(num_subghosts_of_data);
        
        d_flow_model->computeGlobalDerivedCellData();
        
        /*
         * Get the pointers to the maximum wave speeds inside the flow model.
         * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
         */
        
        boost::shared_ptr<pdat::CellData<double> > max_wave_speed_x =
            d_flow_model->getGlobalCellData("MAX_WAVE_SPEED_X");
        
        boost::shared_ptr<pdat::CellData<double> > max_wave_speed_y =
            d_flow_model->getGlobalCellData("MAX_WAVE_SPEED_Y");
        
        boost::shared_ptr<pdat::CellData<double> > max_wave_speed_z =
            d_flow_model->getGlobalCellData("MAX_WAVE_SPEED_Z");
        
        hier::IntVector num_subghosts_max_wave_speed_x = max_wave_speed_x->getGhostCellWidth();
        hier::IntVector num_subghosts_max_wave_speed_y = max_wave_speed_y->getGhostCellWidth();
        hier::IntVector num_subghosts_max_wave_speed_z = max_wave_speed_z->getGhostCellWidth();
        
        TBOX_ASSERT(num_subghosts_max_wave_speed_x == num_ghosts);
        TBOX_ASSERT(num_subghosts_max_wave_speed_y == num_ghosts);
        TBOX_ASSERT(num_subghosts_max_wave_speed_z == num_ghosts);
        
        const int num_ghosts_0 = num_ghosts[0];
        const int num_ghosts_1 = num_ghosts[1];
        const int num_ghosts_2 = num_ghosts[2];
        const int ghostcell_dim_0 = ghostcell_dims[0];
        const int ghostcell_dim_1 = ghostcell_dims[1];
        
        double* max_lambda_x = max_wave_speed_x->getPointer(0);
        double* max_lambda_y = max_wave_speed_y->getPointer(0);
        double* max_lambda_z = max_wave_speed_z->getPointer(0);
        
        for (int k = -num_ghosts_2;
             k < interior_dim_2 + num_ghosts_2;
             k++)
        {
            for (int j = -num_ghosts_1;
                 j < interior_dim_1 + num_ghosts_1;
                 j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd reduction(max:stable_spectral_radius)
#endif
                for (int i = -num_ghosts_0;
                     i < interior_dim_0 + num_ghosts_0;
                     i++)
                {
                    // Compute the linear indices.
                    const int idx = (i + num_ghosts_0) +
                        (j + num_ghosts_1)*ghostcell_dim_0 +
                        (k + num_ghosts_2)*ghostcell_dim_0*
                            ghostcell_dim_1;
                    
                    const double spectral_radius = max_lambda_x[idx]/dx_0 +
                        max_lambda_y[idx]/dx_1 +
                        max_lambda_z[idx]/dx_2;
                    
                    stable_spectral_radius = fmax(stable_spectral_radius, spectral_radius);
                }
            }
        }
        
        /*
         * Unregister the patch and data of all registered derived cell variables in the flow model.
         */
        
        d_flow_model->unregisterPatch();
    }
    
    stable_dt = 1.0/stable_spectral_radius;
    
    t_compute_dt->stop();
    
    return stable_dt;
}


void
Euler::computeFluxesAndSourcesOnPatch(
    hier::Patch& patch,
    const double time,
    const double dt,
    const int RK_step_number,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
    t_compute_fluxes_sources->start();
    
    /*
     * Set zero for the source.
     */
    
    if (data_context)
    {
        boost::shared_ptr<pdat::CellData<double> > data_source(
            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                patch.getPatchData(d_variable_source, data_context)));
        
        data_source->fillAll(0.0);
    }
    else
    {
        boost::shared_ptr<pdat::CellData<double> > data_source(
            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                patch.getPatchData(d_variable_source, getDataContext())));
        
        data_source->fillAll(0.0);
    }
    
    /*
     * Compute the convective flux and source due to splitting of convective term.
     */
    
    if (data_context)
    {
        d_convective_flux_reconstructor->
            computeConvectiveFluxAndSourceOnPatch(
                patch,
                d_variable_convective_flux,
                d_variable_source,
                data_context,
                time,
                dt,
                RK_step_number);
    }
    else
    {
        d_convective_flux_reconstructor->
            computeConvectiveFluxAndSourceOnPatch(
                patch,
                d_variable_convective_flux,
                d_variable_source,
                getDataContext(),
                time,
                dt,
                RK_step_number);
    }
    
    t_compute_fluxes_sources->stop();
}


void
Euler::advanceSingleStepOnPatch(
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
    
    t_advance_step->start();
    
    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(patch_geom);
#endif
    
    const double* dx = patch_geom->getDx();
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    /*
     * Get a vector of pointers to the conservative variables for the current data context (SCRATCH).
     * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
     */
    
    d_flow_model->registerPatchWithDataContext(patch, getDataContext());
    
    std::vector<boost::shared_ptr<pdat::CellData<double> > > conservative_variables =
        d_flow_model->getGlobalCellDataConservativeVariables();
    
    std::vector<hier::IntVector> num_ghosts_conservative_var;
    num_ghosts_conservative_var.reserve(d_flow_model->getNumberOfEquations());
    
    std::vector<hier::IntVector> ghostcell_dims_conservative_var;
    ghostcell_dims_conservative_var.reserve(d_flow_model->getNumberOfEquations());
    
    std::vector<double*> Q;
    Q.reserve(d_flow_model->getNumberOfEquations());
    
    int count_eqn = 0;
    
    for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
    {
        int depth = conservative_variables[vi]->getDepth();
        
        for (int di = 0; di < depth; di++)
        {
            // If the last element of the conservative variable vector is not in the system of
            // equations, ignore it.
            if (count_eqn >= d_flow_model->getNumberOfEquations())
                break;
            
            Q.push_back(conservative_variables[vi]->getPointer(di));
            num_ghosts_conservative_var.push_back(conservative_variables[vi]->getGhostCellWidth());
            ghostcell_dims_conservative_var.push_back(conservative_variables[vi]->getGhostBox().numberCells());
            
            count_eqn++;
        }
    }
    
    d_flow_model->fillZeroGlobalCellDataConservativeVariables();
    
    // Unregister the patch.
    d_flow_model->unregisterPatch();
    
    /*
     * Use alpha, beta and gamma values to update the time-dependent solution, flux and source.
     */
    
    boost::shared_ptr<pdat::SideData<double> > convective_flux(
        BOOST_CAST<pdat::SideData<double>, hier::PatchData>(
            patch.getPatchData(d_variable_convective_flux, getDataContext())));
    
    boost::shared_ptr<pdat::CellData<double> > source(
        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
            patch.getPatchData(d_variable_source, getDataContext())));
    
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(convective_flux);
    TBOX_ASSERT(source);
    
    TBOX_ASSERT(convective_flux->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
    TBOX_ASSERT(source->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
#endif
    
    int num_coeffs = static_cast<int>(alpha.size());
    
    for (int n = 0; n < num_coeffs; n++)
    {
        boost::shared_ptr<pdat::SideData<double> > convective_flux_intermediate(
            BOOST_CAST<pdat::SideData<double>, hier::PatchData>(
                    patch.getPatchData(d_variable_convective_flux, intermediate_context[n])));
        
        boost::shared_ptr<pdat::CellData<double> > source_intermediate(
            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                patch.getPatchData(d_variable_source, intermediate_context[n])));
        
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
        
        d_flow_model->registerPatchWithDataContext(patch, intermediate_context[n]);
        
        std::vector<boost::shared_ptr<pdat::CellData<double> > > conservative_variables_intermediate =
            d_flow_model->getGlobalCellDataConservativeVariables();
        
        std::vector<hier::IntVector> num_ghosts_conservative_var_intermediate;
        num_ghosts_conservative_var_intermediate.reserve(d_flow_model->getNumberOfEquations());
        
        std::vector<hier::IntVector> ghostcell_dims_conservative_var_intermediate;
        ghostcell_dims_conservative_var_intermediate.reserve(d_flow_model->getNumberOfEquations());
        
        std::vector<double*> Q_intermediate;
        Q_intermediate.reserve(d_flow_model->getNumberOfEquations());
        
        count_eqn = 0;
        
        for (int vi = 0; vi < static_cast<int>(conservative_variables_intermediate.size()); vi++)
        {
            int depth = conservative_variables_intermediate[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                // If the last element of the conservative variable vector is not in the system of
                // equations, ignore it.
                if (count_eqn >= d_flow_model->getNumberOfEquations())
                    break;
                
                Q_intermediate.push_back(conservative_variables_intermediate[vi]->getPointer(di));
                num_ghosts_conservative_var_intermediate.push_back(
                    conservative_variables_intermediate[vi]->getGhostCellWidth());
                ghostcell_dims_conservative_var_intermediate.push_back(
                    conservative_variables_intermediate[vi]->getGhostBox().numberCells());
                
                count_eqn++;
            }
        }
        
        // Unregister the patch.
        d_flow_model->unregisterPatch();
        
        if (d_dim == tbox::Dimension(1))
        {
            /*
             * Get the dimension and grid spacing.
             */
            
            const int interior_dim_0 = interior_dims[0];
            
            const double dx_0 = dx[0];
            
            if (alpha[n] != 0.0)
            {
                for (int ei = 0; ei < d_flow_model->getNumberOfEquations(); ei++)
                {
                    const int num_ghosts_0_conservative_var = num_ghosts_conservative_var[ei][0];
                    const int num_ghosts_0_conservative_var_intermediate =
                        num_ghosts_conservative_var_intermediate[ei][0];
                    
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute linear indices of conservative variable data.
                        const int idx = i + num_ghosts_0_conservative_var;
                        const int idx_intermediate = i + num_ghosts_0_conservative_var_intermediate;
                        
                        Q[ei][idx] += alpha[n]*Q_intermediate[ei][idx_intermediate];
                    }
                }
            }
            
            if (beta[n] != 0.0)
            {
                for (int ei = 0; ei < d_flow_model->getNumberOfEquations(); ei++)
                {
                    double* F_x_intermediate = convective_flux_intermediate->getPointer(0, ei);
                    double* S_intermediate = source_intermediate->getPointer(ei);
                    
                    const int num_ghosts_0_conservative_var = num_ghosts_conservative_var[ei][0];
                    
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute linear indices.
                        const int idx = i + num_ghosts_0_conservative_var;
                        const int idx_source = i;
                        const int idx_flux_x = i + 1;
                        
                        Q[ei][idx] += beta[n]*
                            (-(F_x_intermediate[idx_flux_x] - F_x_intermediate[idx_flux_x - 1])/dx_0 +
                             S_intermediate[idx_source]);
                    }
                }
            }
            
            if (gamma[n] != 0.0)
            {
                // Accumulate the flux in the x direction.
                for (int ei = 0; ei < d_flow_model->getNumberOfEquations(); ei++)
                {
                    double* F_x = convective_flux->getPointer(0, ei);
                    double* F_x_intermediate = convective_flux_intermediate->getPointer(0, ei);
                    
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0 + 1; i++)
                    {
                        // Compute linear index.
                        const int idx_flux_x = i;
                        
                        F_x[idx_flux_x] += gamma[n]*F_x_intermediate[idx_flux_x];
                    }                        
                }
                
                // Accumulate the source.
                for (int ei = 0; ei < d_flow_model->getNumberOfEquations(); ei++)
                {
                    double* S = source->getPointer(ei);
                    double* S_intermediate = source_intermediate->getPointer(ei);
                    
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute linear index.
                        const int idx = i;
                        
                        S[idx] += gamma[n]*S_intermediate[idx];
                    }
                }
            } // if (gamma[n] != 0.0)
        } // if (d_dim == tbox::Dimension(1))
        else if (d_dim == tbox::Dimension(2))
        {
            /*
             * Get the dimensions and grid spacings.
             */
            
            const int interior_dim_0 = interior_dims[0];
            const int interior_dim_1 = interior_dims[1];
            
            const double dx_0 = dx[0];
            const double dx_1 = dx[1];
            
            if (alpha[n] != 0.0)
            {
                for (int ei = 0; ei < d_flow_model->getNumberOfEquations(); ei++)
                {
                    const int num_ghosts_0_conservative_var = num_ghosts_conservative_var[ei][0];
                    const int num_ghosts_1_conservative_var = num_ghosts_conservative_var[ei][1];
                    const int ghostcell_dim_0_conservative_var = ghostcell_dims_conservative_var[ei][0];
                    
                    const int num_ghosts_0_conservative_var_intermediate =
                        num_ghosts_conservative_var_intermediate[ei][0];
                    const int num_ghosts_1_conservative_var_intermediate =
                        num_ghosts_conservative_var_intermediate[ei][1];
                    const int ghostcell_dim_0_conservative_var_intermediate =
                        ghostcell_dims_conservative_var_intermediate[ei][0];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute linear indices of conservative data.
                            const int idx = (i + num_ghosts_0_conservative_var) +
                                (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var;
                            
                            const int idx_intermediate = (i + num_ghosts_0_conservative_var_intermediate) +
                                (j + num_ghosts_1_conservative_var_intermediate)*
                                    ghostcell_dim_0_conservative_var_intermediate;
                            
                            Q[ei][idx] += alpha[n]*Q_intermediate[ei][idx_intermediate];
                        }
                    }
                }
            }
            
            if (beta[n] != 0.0)
            {
                for (int ei = 0; ei < d_flow_model->getNumberOfEquations(); ei++)
                {
                    double* F_x_intermediate = convective_flux_intermediate->getPointer(0, ei);
                    double* F_y_intermediate = convective_flux_intermediate->getPointer(1, ei);
                    double* S_intermediate = source_intermediate->getPointer(ei);
                    
                    const int num_ghosts_0_conservative_var = num_ghosts_conservative_var[ei][0];
                    const int num_ghosts_1_conservative_var = num_ghosts_conservative_var[ei][1];
                    const int ghostcell_dim_0_conservative_var = ghostcell_dims_conservative_var[ei][0];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute linear indices.
                            const int idx = (i + num_ghosts_0_conservative_var) +
                                (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var;
                            
                            const int idx_flux_x_L = i +
                                j*(interior_dim_0 + 1);
                            
                            const int idx_flux_x_R = (i + 1) +
                                j*(interior_dim_0 + 1);
                            
                            const int idx_flux_y_B = i +
                                j*interior_dim_0;
                            
                            const int idx_flux_y_T = i +
                                (j + 1)*interior_dim_0;
                            
                            const int idx_source = i +
                                j*interior_dim_0;
                            
                            Q[ei][idx] += beta[n]*
                                (-(F_x_intermediate[idx_flux_x_R] - F_x_intermediate[idx_flux_x_L])/dx_0 -
                                  (F_y_intermediate[idx_flux_y_T] - F_y_intermediate[idx_flux_y_B])/dx_1 +
                                  S_intermediate[idx_source]);
                        }
                    }
                }
            }
            
            if (gamma[n] != 0.0)
            {
                // Accumulate the flux in the x direction.
                for (int ei = 0; ei < d_flow_model->getNumberOfEquations(); ei++)
                {
                    double* F_x = convective_flux->getPointer(0, ei);
                    double* F_x_intermediate = convective_flux_intermediate->getPointer(0, ei);
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = 0; i < interior_dim_0 + 1; i++)
                        {
                            // Compute linear index.
                            const int idx_flux_x = i +
                                j*(interior_dim_0 + 1);
                            
                            F_x[idx_flux_x] += gamma[n]*F_x_intermediate[idx_flux_x];
                        }                        
                    }
                }
                
                // Accumulate the flux in the y direction.
                for (int ei = 0; ei < d_flow_model->getNumberOfEquations(); ei++)
                {
                    double* F_y = convective_flux->getPointer(1, ei);
                    double* F_y_intermediate = convective_flux_intermediate->getPointer(1, ei);
                    
                    for (int j = 0; j < interior_dim_1 + 1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute linear index.
                            const int idx_flux_y = i +
                                j*interior_dim_0;
                            
                            F_y[idx_flux_y] += gamma[n]*F_y_intermediate[idx_flux_y];
                        }
                    }
                }
                
                // Accumulate the source.
                for (int ei = 0; ei < d_flow_model->getNumberOfEquations(); ei++)
                {
                    double* S = source->getPointer(ei);
                    double* S_intermediate = source_intermediate->getPointer(ei);
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute linear index.
                            const int idx = i +
                                j*interior_dim_0;
                            
                            S[idx] += gamma[n]*S_intermediate[idx];
                        }
                    }
                }
            } // if (gamma[n] != 0.0)
        } // if (d_dim == tbox::Dimension(2))
        else if (d_dim == tbox::Dimension(3))
        {
            /*
             * Get the dimensions and grid spacings.
             */
            
            const int interior_dim_0 = interior_dims[0];
            const int interior_dim_1 = interior_dims[1];
            const int interior_dim_2 = interior_dims[2];
            
            const double dx_0 = dx[0];
            const double dx_1 = dx[1];
            const double dx_2 = dx[2];
            
            if (alpha[n] != 0.0)
            {
                for (int ei = 0; ei < d_flow_model->getNumberOfEquations(); ei++)
                {
                    const int num_ghosts_0_conservative_var = num_ghosts_conservative_var[ei][0];
                    const int num_ghosts_1_conservative_var = num_ghosts_conservative_var[ei][1];
                    const int num_ghosts_2_conservative_var = num_ghosts_conservative_var[ei][2];
                    const int ghostcell_dim_0_conservative_var = ghostcell_dims_conservative_var[ei][0];
                    const int ghostcell_dim_1_conservative_var = ghostcell_dims_conservative_var[ei][1];
                    
                    const int num_ghosts_0_conservative_var_intermediate =
                        num_ghosts_conservative_var_intermediate[ei][0];
                    const int num_ghosts_1_conservative_var_intermediate =
                        num_ghosts_conservative_var_intermediate[ei][1];
                    const int num_ghosts_2_conservative_var_intermediate =
                        num_ghosts_conservative_var_intermediate[ei][2];
                    const int ghostcell_dim_0_conservative_var_intermediate =
                        ghostcell_dims_conservative_var_intermediate[ei][0];
                    const int ghostcell_dim_1_conservative_var_intermediate =
                        ghostcell_dims_conservative_var_intermediate[ei][1];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                // Compute linear indices of conservative variable data.
                                const int idx = (i + num_ghosts_0_conservative_var) +
                                    (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                                    (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                                        ghostcell_dim_1_conservative_var;
                                
                                const int idx_intermediate = (i + num_ghosts_0_conservative_var_intermediate) +
                                    (j + num_ghosts_1_conservative_var_intermediate)*
                                        ghostcell_dim_0_conservative_var_intermediate +
                                    (k + num_ghosts_2_conservative_var_intermediate)*
                                        ghostcell_dim_0_conservative_var_intermediate*
                                            ghostcell_dim_1_conservative_var_intermediate;
                                
                                Q[ei][idx] += alpha[n]*Q_intermediate[ei][idx_intermediate];
                            }
                        }
                    }
                }
            }
            
            if (beta[n] != 0.0)
            {
                for (int ei = 0; ei < d_flow_model->getNumberOfEquations(); ei++)
                {
                    double* F_x_intermediate = convective_flux_intermediate->getPointer(0, ei);
                    double* F_y_intermediate = convective_flux_intermediate->getPointer(1, ei);
                    double* F_z_intermediate = convective_flux_intermediate->getPointer(2, ei);
                    double* S_intermediate = source_intermediate->getPointer(ei);
                    
                    const int num_ghosts_0_conservative_var = num_ghosts_conservative_var[ei][0];
                    const int num_ghosts_1_conservative_var = num_ghosts_conservative_var[ei][1];
                    const int num_ghosts_2_conservative_var = num_ghosts_conservative_var[ei][2];
                    const int ghostcell_dim_0_conservative_var = ghostcell_dims_conservative_var[ei][0];
                    const int ghostcell_dim_1_conservative_var = ghostcell_dims_conservative_var[ei][1];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                // Compute linear indices.
                                const int idx = (i + num_ghosts_0_conservative_var) +
                                    (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                                    (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                                        ghostcell_dim_1_conservative_var;
                                
                                const int idx_flux_x_L = i +
                                    j*(interior_dim_0 + 1) +
                                    k*(interior_dim_0 + 1)*interior_dim_1;
                                
                                const int idx_flux_x_R = (i + 1) +
                                    j*(interior_dim_0 + 1) +
                                    k*(interior_dim_0 + 1)*interior_dim_1;
                                
                                const int idx_flux_y_B = i +
                                    j*interior_dim_0 +
                                    k*interior_dim_0*(interior_dim_1 + 1);
                                
                                const int idx_flux_y_T = i +
                                    (j + 1)*interior_dim_0 +
                                    k*interior_dim_0*(interior_dim_1 + 1);
                                
                                const int idx_flux_z_B = i +
                                    j*interior_dim_0 +
                                    k*interior_dim_0*interior_dim_1;
                                
                                const int idx_flux_z_F = i +
                                    j*interior_dim_0 +
                                    (k + 1)*interior_dim_0*interior_dim_1;
                                
                                const int idx_source = i +
                                    j*interior_dim_0 +
                                    k*interior_dim_0*interior_dim_1;
                                
                                Q[ei][idx] += beta[n]*
                                    (-(F_x_intermediate[idx_flux_x_R] - F_x_intermediate[idx_flux_x_L])/dx_0 -
                                      (F_y_intermediate[idx_flux_y_T] - F_y_intermediate[idx_flux_y_B])/dx_1 -
                                      (F_z_intermediate[idx_flux_z_F] - F_z_intermediate[idx_flux_z_B])/dx_2 +
                                      S_intermediate[idx_source]);
                            }
                        }
                    }
                }
            }
            
            if (gamma[n] != 0.0)
            {
                // Accumulate the flux in the x direction.
                for (int ei = 0; ei < d_flow_model->getNumberOfEquations(); ei++)
                {
                    double* F_x = convective_flux->getPointer(0, ei);
                    double* F_x_intermediate = convective_flux_intermediate->getPointer(0, ei);
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = 0; i < interior_dim_0 + 1; i++)
                            {
                                // Compute linear index.
                                const int idx_flux_x = i +
                                    j*(interior_dim_0 + 1) +
                                    k*(interior_dim_0 + 1)*interior_dim_1;
                                
                                F_x[idx_flux_x] += gamma[n]*F_x_intermediate[idx_flux_x];
                            }                        
                        }
                    }
                }
                
                // Accumulate the flux in the y direction.
                for (int ei = 0; ei < d_flow_model->getNumberOfEquations(); ei++)
                {
                    double* F_y = convective_flux->getPointer(1, ei);
                    double* F_y_intermediate = convective_flux_intermediate->getPointer(1, ei);
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1 + 1; j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                // Compute linear index.
                                const int idx_flux_y = i +
                                    j*interior_dim_0 +
                                    k*interior_dim_0*(interior_dim_1 + 1);
                                
                                F_y[idx_flux_y] += gamma[n]*F_y_intermediate[idx_flux_y];
                            }
                        }
                    }
                }
                
                // Accumulate the flux in the z direction.
                for (int ei = 0; ei < d_flow_model->getNumberOfEquations(); ei++)
                {
                    double* F_z = convective_flux->getPointer(2, ei);
                    double* F_z_intermediate = convective_flux_intermediate->getPointer(2, ei);
                    
                    for (int k = 0; k < interior_dim_2 + 1; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                // Compute linear index.
                                const int idx_flux_z = i +
                                    j*interior_dim_0 +
                                    k*interior_dim_0*interior_dim_1;
                                
                                F_z[idx_flux_z] += gamma[n]*F_z_intermediate[idx_flux_z];
                            }
                        }
                    }
                }
                
                // Accumulate the source.
                for (int ei = 0; ei < d_flow_model->getNumberOfEquations(); ei++)
                {
                    double* S = source->getPointer(ei);
                    double* S_intermediate = source_intermediate->getPointer(ei);
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                // Compute linear index.
                                const int idx = i +
                                    j*interior_dim_0 +
                                    k*interior_dim_0*interior_dim_1;
                                
                                S[idx] += gamma[n]*S_intermediate[idx];
                            }
                        }
                    }
                }
            } // if (gamma[n] != 0.0)
        } // if (d_dim == tbox::Dimension(3))
        
        /*
         * Update the conservative variables.
         */
        
        if (beta[n] != 0.0)
        {
            d_flow_model->registerPatchWithDataContext(patch, getDataContext());
            
            d_flow_model->updateGlobalCellDataConservativeVariables();
            
            d_flow_model->unregisterPatch();
        }
    }
    
    t_advance_step->stop();
}


void
Euler::synchronizeFluxes(
    hier::Patch& patch,
    const double time,
    const double dt)
{
    NULL_USE(time);
    NULL_USE(dt);
    
    t_synchronize_fluxes->start();
    
    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(patch_geom);
#endif
    
    const double* dx = patch_geom->getDx();
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    /*
     * Get a vector of pointers to the conservative variables for the current data context (SCRATCH).
     * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
     */
    
    d_flow_model->registerPatchWithDataContext(patch, getDataContext());
    
    std::vector<boost::shared_ptr<pdat::CellData<double> > > conservative_variables =
        d_flow_model->getGlobalCellDataConservativeVariables();
    
    std::vector<hier::IntVector> num_ghosts_conservative_var;
    num_ghosts_conservative_var.reserve(d_flow_model->getNumberOfEquations());
    
    std::vector<hier::IntVector> ghostcell_dims_conservative_var;
    ghostcell_dims_conservative_var.reserve(d_flow_model->getNumberOfEquations());
    
    std::vector<double*> Q;
    Q.reserve(d_flow_model->getNumberOfEquations());
    
    int count_eqn = 0;
    
    for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
    {
        int depth = conservative_variables[vi]->getDepth();
        
        for (int di = 0; di < depth; di++)
        {
            // If the last element of the conservative variable vector is not in the system of equations, ignore it.
            if (count_eqn >= d_flow_model->getNumberOfEquations())
                break;
            
            Q.push_back(conservative_variables[vi]->getPointer(di));
            num_ghosts_conservative_var.push_back(conservative_variables[vi]->getGhostCellWidth());
            ghostcell_dims_conservative_var.push_back(conservative_variables[vi]->getGhostBox().numberCells());
            
            count_eqn++;
        }
    }
    
    // Unregister the patch.
    d_flow_model->unregisterPatch();
    
    boost::shared_ptr<pdat::SideData<double> > convective_flux(
        BOOST_CAST<pdat::SideData<double>, hier::PatchData>(
            patch.getPatchData(d_variable_convective_flux, getDataContext())));
    
    boost::shared_ptr<pdat::CellData<double> > source(
        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
            patch.getPatchData(d_variable_source, getDataContext())));
    
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(convective_flux);
    TBOX_ASSERT(source);
    
    TBOX_ASSERT(convective_flux->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
    TBOX_ASSERT(source->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
#endif
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the dimension and grid spacing.
         */
        
        const int interior_dim_0 = interior_dims[0];
        
        const double dx_0 = dx[0];
        
        for (int ei = 0; ei < d_flow_model->getNumberOfEquations(); ei++)
        {
            double *F_x = convective_flux->getPointer(0, ei);
            double *S = source->getPointer(ei);
            
            const int num_ghosts_0_conservative_var = num_ghosts_conservative_var[ei][0];
            
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute linear indices.
                const int idx = i + num_ghosts_0_conservative_var;
                const int idx_flux_x_L = i;
                const int idx_flux_x_R = i + 1;
                const int idx_source = i;
                
                Q[ei][idx] += (-(F_x[idx_flux_x_R] - F_x[idx_flux_x_L])/dx_0 +
                                S[idx_source]);
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the dimensions and grid spacings.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        const double dx_0 = dx[0];
        const double dx_1 = dx[1];
        
        for (int ei = 0; ei < d_flow_model->getNumberOfEquations(); ei++)
        {
            double *F_x = convective_flux->getPointer(0, ei);
            double *F_y = convective_flux->getPointer(1, ei);
            double *S = source->getPointer(ei);
            
            const int num_ghosts_0_conservative_var = num_ghosts_conservative_var[ei][0];
            const int num_ghosts_1_conservative_var = num_ghosts_conservative_var[ei][1];
            const int ghostcell_dim_0_conservative_var = ghostcell_dims_conservative_var[ei][0];
            
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute linear indices.
                    const int idx = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var;
                    
                    const int idx_flux_x_L = i +
                        j*(interior_dim_0 + 1);
                    
                    const int idx_flux_x_R = (i + 1) +
                        j*(interior_dim_0 + 1);
                    
                    const int idx_flux_y_B = i +
                        j*interior_dim_0;
                    
                    const int idx_flux_y_T = i +
                        (j + 1)*interior_dim_0;
                    
                    const int idx_source = i +
                        j*interior_dim_0;
                    
                    Q[ei][idx] += (-(F_x[idx_flux_x_R] - F_x[idx_flux_x_L])/dx_0 -
                                    (F_y[idx_flux_y_T] - F_y[idx_flux_y_B])/dx_1 +
                                    S[idx_source]);
                }
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the dimensions and grid spacings.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        const int interior_dim_2 = interior_dims[2];
        
        const double dx_0 = dx[0];
        const double dx_1 = dx[1];
        const double dx_2 = dx[2];
        
        for (int ei = 0; ei < d_flow_model->getNumberOfEquations(); ei++)
        {
            double *F_x = convective_flux->getPointer(0, ei);
            double *F_y = convective_flux->getPointer(1, ei);
            double *F_z = convective_flux->getPointer(2, ei);
            double *S = source->getPointer(ei);
            
            const int num_ghosts_0_conservative_var = num_ghosts_conservative_var[ei][0];
            const int num_ghosts_1_conservative_var = num_ghosts_conservative_var[ei][1];
            const int num_ghosts_2_conservative_var = num_ghosts_conservative_var[ei][2];
            const int ghostcell_dim_0_conservative_var = ghostcell_dims_conservative_var[ei][0];
            const int ghostcell_dim_1_conservative_var = ghostcell_dims_conservative_var[ei][1];
            
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute linear indices.
                        const int idx = (i + num_ghosts_0_conservative_var) +
                            (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                            (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                                ghostcell_dim_1_conservative_var;
                        
                        const int idx_flux_x_L = i +
                            j*(interior_dim_0 + 1) +
                            k*(interior_dim_0 + 1)*interior_dim_1;
                        
                        const int idx_flux_x_R = (i + 1) +
                            j*(interior_dim_0 + 1) +
                            k*(interior_dim_0 + 1)*interior_dim_1;
                        
                        const int idx_flux_y_B = i +
                            j*interior_dim_0 +
                            k*interior_dim_0*(interior_dim_1 + 1);
                        
                        const int idx_flux_y_T = i +
                            (j + 1)*interior_dim_0 +
                            k*interior_dim_0*(interior_dim_1 + 1);
                        
                        const int idx_flux_z_B = i +
                            j*interior_dim_0 +
                            k*interior_dim_0*interior_dim_1;
                        
                        const int idx_flux_z_F = i +
                            j*interior_dim_0 +
                            (k + 1)*interior_dim_0*interior_dim_1;
                        
                        const int idx_source = i +
                            j*interior_dim_0 +
                            k*interior_dim_0*interior_dim_1;
                        
                        Q[ei][idx] += (-(F_x[idx_flux_x_R] - F_x[idx_flux_x_L])/dx_0 -
                                        (F_y[idx_flux_y_T] - F_y[idx_flux_y_B])/dx_1 -
                                        (F_z[idx_flux_z_F] - F_z[idx_flux_z_B])/dx_2 +
                                        S[idx_source]);
                    }
                }
            }
        }
    }
    
    /*
     * Update the conservative variables.
     */
    
    d_flow_model->registerPatchWithDataContext(patch, getDataContext());
    
    d_flow_model->updateGlobalCellDataConservativeVariables();
    
    d_flow_model->unregisterPatch();
    
    t_synchronize_fluxes->stop();
}


/*
 * Preprocess before tagging cells using value detector.
 */
void
Euler::preprocessTagCellsValueDetector(
   const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
   const int level_number,
   const double regrid_time,
   const bool initial_error,
   const bool uses_gradient_detector_too,
   const bool uses_multiresolution_detector_too,
   const bool uses_integral_detector_too,
   const bool uses_richardson_extrapolation_too)
{
    NULL_USE(regrid_time);
    NULL_USE(initial_error);
    NULL_USE(uses_gradient_detector_too);
    NULL_USE(uses_multiresolution_detector_too);
    NULL_USE(uses_integral_detector_too);
    NULL_USE(uses_richardson_extrapolation_too);
    
    if (d_value_tagger != nullptr)
    {
        boost::shared_ptr<hier::PatchLevel> level(
            patch_hierarchy->getPatchLevel(level_number));
        
        for (hier::PatchLevel::iterator ip(level->begin());
             ip != level->end();
             ip++)
        {
            const boost::shared_ptr<hier::Patch>& patch = *ip;
            
            d_value_tagger->computeValueTaggerValuesOnPatch(
                *patch,
                getDataContext());
        }
        
        d_value_tagger->getValueStatistics(
            patch_hierarchy,
            level_number,
            getDataContext());
    }
}


/*
 * Tag cells for refinement using value detector.
 */
void
Euler::tagCellsOnPatchValueDetector(
    hier::Patch& patch,
    const double regrid_time,
    const bool initial_error,
    const int tag_indx,
    const bool uses_gradient_detector_too,
    const bool uses_multiresolution_detector_too,
    const bool uses_integral_detector_too,
    const bool uses_richardson_extrapolation_too)
{
    NULL_USE(regrid_time);
    NULL_USE(initial_error);
    
    t_tagvalue->start();
    
    // Get the tags.
    boost::shared_ptr<pdat::CellData<int> > tags(
        BOOST_CAST<pdat::CellData<int>, hier::PatchData>(
            patch.getPatchData(tag_indx)));
    
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(tags);
    TBOX_ASSERT(tags->getGhostCellWidth() == 0);
#endif
    
    // Initialize values of all tags to zero.
    if ((!uses_richardson_extrapolation_too) &&
        (!uses_integral_detector_too) &&
        (!uses_multiresolution_detector_too) &&
        (!uses_gradient_detector_too)) 
    {
        tags->fillAll(0);
    }
    
    // Tag the cells by using d_value_tagger.
    if (d_value_tagger != nullptr)
    {
        d_value_tagger->tagCellsOnPatch(
            patch,
            tags,
            getDataContext());
    }
    
    t_tagvalue->stop();
}


/*
 * Preprocess before tagging cells using gradient detector.
 */
void
Euler::preprocessTagCellsGradientDetector(
   const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
   const int level_number,
   const double regrid_time,
   const bool initial_error,
   const bool uses_value_detector_too,
   const bool uses_multiresolution_detector_too,
   const bool uses_integral_detector_too,
   const bool uses_richardson_extrapolation_too)
{
    NULL_USE(regrid_time);
    NULL_USE(initial_error);
    NULL_USE(uses_value_detector_too);
    NULL_USE(uses_multiresolution_detector_too);
    NULL_USE(uses_integral_detector_too);
    NULL_USE(uses_richardson_extrapolation_too);
    
    if (d_gradient_tagger != nullptr)
    {
        boost::shared_ptr<hier::PatchLevel> level(
            patch_hierarchy->getPatchLevel(level_number));
        
        for (hier::PatchLevel::iterator ip(level->begin());
             ip != level->end();
             ip++)
        {
            const boost::shared_ptr<hier::Patch>& patch = *ip;
            
            d_gradient_tagger->computeGradientSensorValuesOnPatch(
                *patch,
                getDataContext());
        }
        
        d_gradient_tagger->getSensorValueStatistics(
            patch_hierarchy,
            level_number,
            getDataContext());
    }
}


/*
 * Tag cells for refinement using gradient detector.
 */
void
Euler::tagCellsOnPatchGradientDetector(
    hier::Patch& patch,
    const double regrid_time,
    const bool initial_error,
    const int tag_indx,
    const bool uses_value_detector_too,
    const bool uses_multiresolution_detector_too,
    const bool uses_integral_detector_too,
    const bool uses_richardson_extrapolation_too)
{
    NULL_USE(regrid_time);
    NULL_USE(initial_error);
    NULL_USE(uses_value_detector_too);
    
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
    if ((!uses_richardson_extrapolation_too) &&
        (!uses_integral_detector_too) &&
        (!uses_multiresolution_detector_too)) 
    {
        tags->fillAll(0);
    }
    
    // Tag the cells by using d_gradient_tagger.
    if (d_gradient_tagger != nullptr)
    {
        d_gradient_tagger->tagCellsOnPatch(
            patch,
            tags,
            getDataContext());
    }
    
    t_taggradient->stop();
}


/*
 * Preprocess before tagging cells using multiresolution detector.
 */
void
Euler::preprocessTagCellsMultiresolutionDetector(
   const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
   const int level_number,
   const double regrid_time,
   const bool initial_error,
   const bool uses_value_detector_too,
   const bool uses_gradient_detector_too,
   const bool uses_integral_detector_too,
   const bool uses_richardson_extrapolation_too)
{
    NULL_USE(regrid_time);
    NULL_USE(initial_error);
    NULL_USE(uses_value_detector_too);
    NULL_USE(uses_gradient_detector_too);
    NULL_USE(uses_integral_detector_too);
    NULL_USE(uses_richardson_extrapolation_too);
    
    if (d_multiresolution_tagger != nullptr)
    {
        boost::shared_ptr<hier::PatchLevel> level(
            patch_hierarchy->getPatchLevel(level_number));
        
        for (hier::PatchLevel::iterator ip(level->begin());
             ip != level->end();
             ip++)
        {
            const boost::shared_ptr<hier::Patch>& patch = *ip;
            
            d_multiresolution_tagger->computeMultiresolutionSensorValuesOnPatch(
                *patch,
                getDataContext());
        }
        
        d_multiresolution_tagger->getSensorValueStatistics(
            patch_hierarchy,
            level_number,
            getDataContext());
    }
}


/*
 * Tag cells for refinement using multiresolution detector.
 */
void
Euler::tagCellsOnPatchMultiresolutionDetector(
    hier::Patch& patch,
    const double regrid_time,
    const bool initial_error,
    const int tag_indx,
    const bool uses_value_detector_too,
    const bool uses_gradient_detector_too,
    const bool uses_integral_detector_too,
    const bool uses_richardson_extrapolation_too)
{
    NULL_USE(regrid_time);
    NULL_USE(initial_error);
    NULL_USE(uses_value_detector_too);
    NULL_USE(uses_gradient_detector_too);
    
    t_tagmultiresolution->start();
    
    // Get the tags.
    boost::shared_ptr<pdat::CellData<int> > tags(
        BOOST_CAST<pdat::CellData<int>, hier::PatchData>(
            patch.getPatchData(tag_indx)));
    
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(tags);
    TBOX_ASSERT(tags->getGhostCellWidth() == 0);
#endif
    
    // Initialize values of all tags to zero.
    if ((!uses_richardson_extrapolation_too) &&
        (!uses_integral_detector_too))
    {
        tags->fillAll(0);
    }
    
    // Tag the cells by using d_multiresolution_tagger.
    if (d_multiresolution_tagger != nullptr)
    {
        d_multiresolution_tagger->tagCellsOnPatch(
            patch,
            tags,
            getDataContext());
    }
    
    t_tagmultiresolution->stop();
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
    
    restart_db->putString("d_flow_model_str", d_flow_model_str);
    
    // Put the properties of d_flow_model into the restart database.
    boost::shared_ptr<tbox::Database> restart_flow_model_db =
        restart_db->putDatabase("d_flow_model_db");
    d_flow_model->putToRestart(restart_flow_model_db);
    
    restart_db->putString("d_convective_flux_reconstructor_str", d_convective_flux_reconstructor_str);
    
    // Put the properties of d_convective_flux_reconstructor into the restart database.
    boost::shared_ptr<tbox::Database> restart_convective_flux_reconstructor_db =
        restart_db->putDatabase("d_convective_flux_reconstructor_db");
    d_convective_flux_reconstructor->putToRestart(restart_convective_flux_reconstructor_db);
    
    boost::shared_ptr<tbox::Database> restart_Euler_boundary_conditions_db =
        restart_db->putDatabase("d_Euler_boundary_conditions_db");
    
    d_Euler_boundary_conditions->putToRestart(restart_Euler_boundary_conditions_db);
    
    if (d_value_tagger != nullptr)
    {
        boost::shared_ptr<tbox::Database> restart_value_tagger_db =
            restart_db->putDatabase("d_value_tagger_db");
        
        d_value_tagger->putToRestart(restart_value_tagger_db);
    }
    
    if (d_gradient_tagger != nullptr)
    {
        boost::shared_ptr<tbox::Database> restart_gradient_tagger_db =
            restart_db->putDatabase("d_gradient_tagger_db");
        
        d_gradient_tagger->putToRestart(restart_gradient_tagger_db);
    }
    
    if (d_multiresolution_tagger != nullptr)
    {
        boost::shared_ptr<tbox::Database> restart_multiresolution_tagger_db =
            restart_db->putDatabase("d_multiresolution_tagger_db");
        
        d_multiresolution_tagger->putToRestart(restart_multiresolution_tagger_db);
    }
}


#ifdef HAVE_HDF5
void
Euler::registerVisItDataWriter(
    const boost::shared_ptr<ExtendedVisItDataWriter>& viz_writer)
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
    bool data_on_patch = d_flow_model->
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
    
    // Print all characteristics of d_flow_model.
    d_flow_model_manager->printClassData(os);
    
    os << "d_num_species = " << d_num_species << std::endl;
    
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
     * Print data of d_convective_flux_reconstructor.
     */
    
    d_convective_flux_reconstructor->printClassData(os);
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
Euler::printErrorStatistics(
    std::ostream& os,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy) const
{
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    math::HierarchyCellDataOpsReal<double> cell_double_operator(patch_hierarchy, 0, 0);
    
    std::vector<std::string> variable_names = d_flow_model->getNamesOfConservativeVariables();
    
    std::vector<boost::shared_ptr<pdat::CellVariable<double> > > variables =
        d_flow_model->getConservativeVariables();
    
   if (d_project_name == "2D single-species convergence test")
   {
        for (int vi = 0; vi < static_cast<int>(variables.size()); vi++)
        {
            if (variable_names[vi] == "density")
            {
                boost::shared_ptr<hier::PatchLevel> level_root(
                    patch_hierarchy->getPatchLevel(0));
                
                double error_sum_local = 0.0;
                double error_squared_sum_local = 0.0;
                double error_max_local = 0.0;
                
                for (hier::PatchLevel::iterator ip(level_root->begin());
                     ip != level_root->end();
                     ip++)
                {
                    const boost::shared_ptr<hier::Patch>& patch = *ip;
                    
                    // Get the dimensions of box that covers the interior of Patch.
                    hier::Box patch_box = patch->getBox();
                    const hier::IntVector patch_dims = patch_box.numberCells();
                    
                    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                            patch->getPatchGeometry()));
                    
                    const double* const dx = patch_geom->getDx();
                    const double* const patch_xlo = patch_geom->getXLower();
                    
                    boost::shared_ptr<pdat::CellData<double> > density(
                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                            patch->getPatchData(variables[vi], d_plot_context)));
                    
                    double* rho = density->getPointer(0);
                    
                    for (int j = 0; j < patch_dims[1]; j++)
                    {
                        for (int i = 0; i < patch_dims[0]; i++)
                        {
                            // Compute index into linear data array.
                            int idx_cell = i + j*patch_dims[0];
                            
                            // Compute the coordinates.
                            double x[2];
                            x[0] = patch_xlo[0] + (i + 0.5)*dx[0];
                            x[1] = patch_xlo[1] + (j + 0.5)*dx[1];
                            
                            const double rho_exact   = 1.0 + 0.5*sin(M_PI*(x[0] + x[1]));
                            const double error = fabs(rho_exact - rho[idx_cell]);
                            
                            error_sum_local += dx[0]*dx[0]*error;
                            error_squared_sum_local += dx[0]*dx[0]*error*error;
                            error_max_local = fmax(error, error_max_local);
                        }
                    }
                }
                
                double error_sum_global = 0.0;
                double error_squared_sum_global = 0.0;
                double error_max_global = 0.0;
                
                mpi.Allreduce(
                    &error_sum_local,
                    &error_sum_global,
                    1,
                    MPI_DOUBLE,
                    MPI_SUM);
                
                mpi.Allreduce(
                    &error_squared_sum_local,
                    &error_squared_sum_global,
                    1,
                    MPI_DOUBLE,
                    MPI_SUM);
                
                mpi.Allreduce(
                    &error_max_local,
                    &error_max_global,
                    1,
                    MPI_DOUBLE,
                    MPI_MAX);
                
                
                os.precision(17);
                os << "L1_error: " << std::scientific << error_sum_global << std::endl;
                os << "L2_error: " << std::scientific << sqrt(error_squared_sum_global) << std::endl;
                os << "Linf_error: " << std::scientific << error_max_global << std::endl;
            }
        }
    }
    else if (d_project_name == "2D multi-species convergence test")
    {
        for (int vi = 0; vi < static_cast<int>(variables.size()); vi++)
        {
            if (variable_names[vi] == "volume fraction")
            {
                boost::shared_ptr<hier::PatchLevel> level_root(
                patch_hierarchy->getPatchLevel(0));
                
                double error_sum_local = 0.0;
                double error_squared_sum_local = 0.0;
                double error_max_local = 0.0;
                
                for (hier::PatchLevel::iterator ip(level_root->begin());
                     ip != level_root->end();
                     ip++)
                {
                    const boost::shared_ptr<hier::Patch>& patch = *ip;
                    
                    // Get the dimensions of box that covers the interior of Patch.
                    hier::Box patch_box = patch->getBox();
                    const hier::IntVector patch_dims = patch_box.numberCells();
                    
                    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                            patch->getPatchGeometry()));
                    
                    const double* const dx = patch_geom->getDx();
                    const double* const patch_xlo = patch_geom->getXLower();
                    
                    boost::shared_ptr<pdat::CellData<double> > volume_fraction(
                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                            patch->getPatchData(variables[vi], d_plot_context)));
                    
                    double* Z_1 = volume_fraction->getPointer(0);
                    
                    for (int j = 0; j < patch_dims[1]; j++)
                    {
                        for (int i = 0; i < patch_dims[0]; i++)
                        {
                            // Compute index into linear data array.
                            int idx_cell = i + j*patch_dims[0];
                            
                            // Compute the coordinates.
                            double x[2];
                            x[0] = patch_xlo[0] + (i + 0.5)*dx[0];
                            x[1] = patch_xlo[1] + (j + 0.5)*dx[1];
                            
                            const double Z_1_exact   = 0.5 + 0.25*sin(M_PI*(x[0] + x[1]));
                            const double error = fabs(Z_1_exact - Z_1[idx_cell]);
                            
                            error_sum_local += dx[0]*dx[0]*error;
                            error_squared_sum_local += dx[0]*dx[0]*error*error;
                            error_max_local = fmax(error, error_max_local);
                        }
                    }        
                }
                
                double error_sum_global = 0.0;
                double error_squared_sum_global = 0.0;
                double error_max_global = 0.0;
                
                mpi.Allreduce(
                    &error_sum_local,
                    &error_sum_global,
                    1,
                    MPI_DOUBLE,
                    MPI_SUM);
                
                mpi.Allreduce(
                    &error_squared_sum_local,
                    &error_squared_sum_global,
                    1,
                    MPI_DOUBLE,
                    MPI_SUM);
                
                mpi.Allreduce(
                    &error_max_local,
                    &error_max_global,
                    1,
                    MPI_DOUBLE,
                    MPI_MAX);
                
                
                os.precision(17);
                os << "L1_error: " << std::scientific << error_sum_global << std::endl;
                os << "L2_error: " << std::scientific << sqrt(error_squared_sum_global) << std::endl;
                os << "Linf_error: " << std::scientific << error_max_global << std::endl;
            }
        }
    }
}


void
Euler::printDataStatistics(
    std::ostream& os,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy) const
{
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    math::HierarchyCellDataOpsReal<double> cell_double_operator(patch_hierarchy, 0, 0);
    
    hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();
    
    std::vector<std::string> variable_names = d_flow_model->getNamesOfConservativeVariables();
    
    std::vector<boost::shared_ptr<pdat::CellVariable<double> > > variables =
        d_flow_model->getConservativeVariables();
    
    for (int vi = 0; vi < static_cast<int>(variables.size()); vi++)
    {
        const int var_depth = variables[vi]->getDepth();
        const int var_id = variable_db->mapVariableAndContextToIndex(
            variables[vi],
            d_plot_context);
        
        double var_max_local = cell_double_operator.max(var_id);
        double var_min_local = cell_double_operator.min(var_id);
        
        double var_max_global = 0.0;
        double var_min_global = 0.0;
        
        mpi.Allreduce(
            &var_max_local,
            &var_max_global,
            1,
            MPI_DOUBLE,
            MPI_MAX);
        
        mpi.Allreduce(
            &var_min_local,
            &var_min_global,
            1,
            MPI_DOUBLE,
            MPI_MAX);
        
        if (var_depth > 1)
        {
            os << "Max/min " << variable_names[vi] << " component: " << var_max_global << "/" << var_min_global << std::endl;
        }
        else
        {
            os << "Max/min " << variable_names[vi] << ": " << var_max_global << "/" << var_min_global << std::endl;
        }
    }
}


/**
 * Compute variables for computing the statistics of data.
 */
void
Euler::computeStatisticsVariables(
   const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy)
{
    d_flow_model->setupStatisticsUtilities();
    
    boost::shared_ptr<FlowModelStatisticsUtilities> flow_model_statistics_utilities =
        d_flow_model->getFlowModelStatisticsUtilities();
    
    flow_model_statistics_utilities->computeVariables(
        patch_hierarchy,
        getDataContext());
}


/**
 * Output the statistics of data.
 */
void
Euler::outputDataStatistics(
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const double output_time)
{
    if ((!d_stat_dump_filename.empty()))
    {
        const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
        
        if (mpi.getRank() == 0)
        {
            std::ofstream f_out;
            f_out.open(d_stat_dump_filename.c_str(), std::ios::app);
            
            if (!f_out.is_open())
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "Failed to open file to output statistics!"
                    << std::endl);
            }
            
            f_out << std::fixed << std::setprecision(std::numeric_limits<double>::digits10) << output_time;
            f_out.close();
        }
        
        d_flow_model->setupStatisticsUtilities();
        
        boost::shared_ptr<FlowModelStatisticsUtilities> flow_model_statistics_utilities =
            d_flow_model->getFlowModelStatisticsUtilities();
        
        flow_model_statistics_utilities->outputStatisticalQuantities(
            d_stat_dump_filename,
            patch_hierarchy,
            getDataContext(),
            output_time);
        
        if (mpi.getRank() == 0)
        {
            std::ofstream f_out;
            f_out.open(d_stat_dump_filename.c_str(), std::ios::app);
            
            if (!f_out.is_open())
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "Failed to open file to output statistics!"
                    << std::endl);
            }
            
            f_out << std::endl;
            f_out.close();
        }
    }
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
         * Get the flow model.
         */
        if (input_db->keyExists("flow_model"))
        {
            d_flow_model_str = input_db->getString("flow_model");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Key data 'flow_model' not found in input database."
                << " Compressible flow model is unknown."
                << std::endl);            
        }
        
        /*
         * Get the database of the flow model.
         */
        if (input_db->keyExists("Flow_model"))
        {
            d_flow_model_db = input_db->getDatabase("Flow_model");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Key data 'Flow_model' not found in input database."
                << std::endl);
        }
        
        /*
         * Get the convective flux reconstructor.
         */
        if (input_db->keyExists("convective_flux_reconstructor"))
        {
            d_convective_flux_reconstructor_str = input_db->getString("convective_flux_reconstructor");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Key data 'convective_flux_reconstructor' not found in input database."
                << std::endl);            
        }
        
        /*
         * Get the database of the convective flux reconstructor.
         */
        if (input_db->keyExists("Convective_flux_reconstructor"))
        {
            d_convective_flux_reconstructor_db = input_db->getDatabase("Convective_flux_reconstructor");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Key data 'Convective_flux_reconstructor' not found in input database."
                << std::endl);
        }
        
        if (input_db->keyExists("Value_tagger"))
        {
            d_value_tagger_db =
                input_db->getDatabase("Value_tagger");
        }
        
        if (input_db->keyExists("Gradient_tagger"))
        {
            d_gradient_tagger_db =
                input_db->getDatabase("Gradient_tagger");
        }
        
        if (input_db->keyExists("Multiresolution_tagger"))
        {
            d_multiresolution_tagger_db =
                input_db->getDatabase("Multiresolution_tagger");
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
    
    d_flow_model_str = db->getString("d_flow_model_str");
    
    d_flow_model_db = db->getDatabase("d_flow_model_db");
    
    d_convective_flux_reconstructor_str = db->getString("d_convective_flux_reconstructor_str");
    
    d_convective_flux_reconstructor_db = db->getDatabase("d_convective_flux_reconstructor_db");
    
    d_Euler_boundary_conditions_db = db->getDatabase("d_Euler_boundary_conditions_db");
    
    d_Euler_boundary_conditions_db_is_from_restart = true;
    
    if (db->keyExists("d_value_tagger_db"))
    {
        d_value_tagger_db = db->getDatabase("d_value_tagger_db");
    }
    
    if (db->keyExists("d_gradient_tagger_db"))
    {
        d_gradient_tagger_db = db->getDatabase("d_gradient_tagger_db");
    }
    
    if (db->keyExists("d_multiresolution_tagger_db"))
    {
        d_multiresolution_tagger_db = db->getDatabase("d_multiresolution_tagger_db");
    }
}
