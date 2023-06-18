#include "apps/Navier-Stokes/NavierStokes.hpp"

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

HAMERS_SHARED_PTR<tbox::Timer> NavierStokes::t_init;
HAMERS_SHARED_PTR<tbox::Timer> NavierStokes::t_compute_dt;
HAMERS_SHARED_PTR<tbox::Timer> NavierStokes::t_compute_fluxes_sources;
HAMERS_SHARED_PTR<tbox::Timer> NavierStokes::t_advance_step;
HAMERS_SHARED_PTR<tbox::Timer> NavierStokes::t_synchronize_fluxes;
HAMERS_SHARED_PTR<tbox::Timer> NavierStokes::t_setphysbcs;
HAMERS_SHARED_PTR<tbox::Timer> NavierStokes::t_tagvalue;
HAMERS_SHARED_PTR<tbox::Timer> NavierStokes::t_tagimmersedbdry;
HAMERS_SHARED_PTR<tbox::Timer> NavierStokes::t_taggradient;
HAMERS_SHARED_PTR<tbox::Timer> NavierStokes::t_tagmultiresolution;

NavierStokes::NavierStokes(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const HAMERS_SHARED_PTR<tbox::Database>& input_db,
    const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
    const std::string& stat_dump_filename):
        RungeKuttaPatchStrategy(),
        d_object_name(object_name),
        d_dim(dim),
        d_grid_geometry(grid_geometry),
        d_monitoring_stat_dump_filename("monitoring_stats.txt"),
        d_stat_dump_filename(stat_dump_filename),
        d_use_nonuniform_workload(false),
        d_use_conservative_form_diffusive_flux(true),
        d_Navier_Stokes_boundary_conditions_db_is_from_restart(false),
        d_use_immersed_boundaries(false)
{
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(input_db);
    TBOX_ASSERT(grid_geometry);
    
    tbox::RestartManager::getManager()->registerRestartItem(d_object_name, this);
    
    if (!t_init)
    {
        t_init = tbox::TimerManager::getManager()->
            getTimer("NavierStokes::initializeDataOnPatch()");
        t_compute_dt = tbox::TimerManager::getManager()->
            getTimer("NavierStokes::computeStableDtOnPatch()");
        t_compute_fluxes_sources = tbox::TimerManager::getManager()->
            getTimer("NavierStokes::computeHyperbolicFluxesOnPatch()");
        t_advance_step = tbox::TimerManager::getManager()->
            getTimer("NavierStokes::advanceSingleStepOnPatch()");
        t_synchronize_fluxes = tbox::TimerManager::getManager()->
            getTimer("NavierStokes::synchronizeHyperbolicFluxes()");
        t_setphysbcs = tbox::TimerManager::getManager()->
            getTimer("NavierStokes::setPhysicalBoundaryConditions()");
        t_tagimmersedbdry = tbox::TimerManager::getManager()->
            getTimer("NavierStokes::tagCellsOnPatchImmersedBdryDetector()");
        t_tagvalue = tbox::TimerManager::getManager()->
            getTimer("NavierStokes::tagCellsOnPatchValueDetector()");
        t_taggradient = tbox::TimerManager::getManager()->
            getTimer("NavierStokes::tagCellsOnPatchGradientDetector()");
        t_tagmultiresolution = tbox::TimerManager::getManager()->
            getTimer("NavierStokes::tagCellsOnPatchMultiresolutionDetector()");
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
    
    if (d_use_ghost_cell_immersed_boundary_method)
    {
        d_use_immersed_boundaries = true;
    }
    
    /*
     * Initialize d_flow_model_manager and get the flow model object.
     */
    
    d_flow_model_manager.reset(new FlowModelManager(
        "d_flow_model_manager",
        d_project_name,
        d_dim,
        d_grid_geometry,
        d_num_species,
        d_flow_model_db,
        d_flow_model_str));
    
    d_flow_model = d_flow_model_manager->getFlowModel();
    
    /*
     * Initialize the immersed boundaries and the flow model immersed boundary method objects.
     */
    
    if (d_use_immersed_boundaries)
    {
        d_immersed_boundaries.reset(new ImmersedBoundaries(
            "d_immersed_boundaries",
            d_project_name,
            d_dim,
            d_grid_geometry));
        
        d_flow_model->initializeImmersedBoundaryMethod(
            d_immersed_boundaries,
            d_immersed_boundary_method_db);
    }
    
    /*
     * Initialize d_convective_flux_reconstructor_manager and get the convective flux reconstructor object.
     */
    
    d_convective_flux_reconstructor_manager.reset(new ConvectiveFluxReconstructorManager(
        "d_convective_flux_reconstructor_manager",
        d_dim,
        d_grid_geometry,
        d_flow_model->getNumberOfEquations(),
        d_flow_model_manager->getFlowModelType(),
        d_flow_model,
        d_convective_flux_reconstructor_db,
        d_convective_flux_reconstructor_str));
    
    d_convective_flux_reconstructor = d_convective_flux_reconstructor_manager->getConvectiveFluxReconstructor();
    
    if (d_use_conservative_form_diffusive_flux)
    {
        /*
         * Initialize d_diffusive_flux_reconstructor_manager and get the diffusive flux reconstructor object.
         */
        
        d_diffusive_flux_reconstructor_manager.reset(new DiffusiveFluxReconstructorManager(
            "d_diffusive_flux_reconstructor_manager",
            d_dim,
            d_grid_geometry,
            d_flow_model->getNumberOfEquations(),
            d_flow_model,
            d_diffusive_flux_reconstructor_db,
            d_diffusive_flux_reconstructor_str));
        
        d_diffusive_flux_reconstructor = d_diffusive_flux_reconstructor_manager->getDiffusiveFluxReconstructor();
    }
    else
    {
        /*
         * Initialize d_non-conservative_diffusive_flux_divergence_operator_manager and get the
         * non-conservative diffusive flux divergence operator object.
         */
        
        d_nonconservative_diffusive_flux_divergence_operator_manager.reset(
            new NonconservativeDiffusiveFluxDivergenceOperatorManager(
                "d_nonconservative_diffusive_flux_divergence_operator_manager",
                d_dim,
                d_grid_geometry,
                d_flow_model->getNumberOfEquations(),
                d_flow_model,
                d_nonconservative_diffusive_flux_divergence_operator_db,
                d_nonconservative_diffusive_flux_divergence_operator_str));
        
        d_nonconservative_diffusive_flux_divergence_operator =
            d_nonconservative_diffusive_flux_divergence_operator_manager->
                getNonconservativeDiffusiveFluxDivergenceOperator();
    }
    
    /*
     * Initialize d_Navier_Stokes_initial_conditions.
     */
    
    d_Navier_Stokes_initial_conditions.reset(new NavierStokesInitialConditions(
        "d_Navier_Stokes_initial_conditions",
        d_project_name,
        d_dim,
        d_grid_geometry,
        d_flow_model_manager->getFlowModelType(),
        d_flow_model,
        d_Navier_Stokes_initial_conditions_db));
    
    /*
     * Initialize d_Navier_Stokes_boundary_conditions.
     */
    
    d_Navier_Stokes_boundary_conditions.reset(new NavierStokesBoundaryConditions(
        "d_Navier_Stokes_boundary_conditions",
        d_project_name,
        d_dim,
        d_grid_geometry,
        d_flow_model_manager->getFlowModelType(),
        d_flow_model,
        d_Navier_Stokes_boundary_conditions_db,
        d_Navier_Stokes_boundary_conditions_db_is_from_restart));
    
    /*
     * Initialize d_Navier_Stokes_error_statistics.
     */
    
    d_Navier_Stokes_error_statistics.reset(new NavierStokesErrorStatistics(
        "d_Navier_Stokes_error_statistics",
        d_project_name,
        d_dim,
        d_grid_geometry,
        d_flow_model_manager->getFlowModelType(),
        d_flow_model));
    
    /*
     * Initialize d_immersed_boundary_tagger.
     */
    
    if (d_use_immersed_boundaries)
    {
        hier::IntVector num_cells_buffer_required = hier::IntVector::getZero(d_dim);
        num_cells_buffer_required = hier::IntVector::max(
            num_cells_buffer_required,
            d_convective_flux_reconstructor->getConvectiveFluxNumberOfGhostCells());
        
        hier::IntVector num_ghosts_diffusive_flux = hier::IntVector::getZero(d_dim);
        
        if (d_use_conservative_form_diffusive_flux)
        {
            num_ghosts_diffusive_flux = d_diffusive_flux_reconstructor->getDiffusiveFluxNumberOfGhostCells();
        }
        else
        {
            num_ghosts_diffusive_flux = d_nonconservative_diffusive_flux_divergence_operator->
                getNonconservativeDiffusiveFluxDivergenceOperatorNumberOfGhostCells();
        }
        
        bool use_transverse_derivatives_bc = d_Navier_Stokes_boundary_conditions->useTransverseDerivativesBoundaryConditions();
        
        // Increase the number of ghosts for computing transverse derivatives at corner ghost cells for diffusive/viscous flux.
        // Euler assumes corner ghost cells are not used.
        if (use_transverse_derivatives_bc)
        {
            hier::IntVector num_ghosts_transverse_derivatives_bc = d_Navier_Stokes_boundary_conditions->
                getBoundaryConditionsTransverseDerivativesNumberOfGhostCells();
            
            num_ghosts_diffusive_flux = num_ghosts_diffusive_flux + num_ghosts_transverse_derivatives_bc;
        }
        
        num_cells_buffer_required = hier::IntVector::max(num_cells_buffer_required, num_ghosts_diffusive_flux);
        
        if (hier::IntVector::getOne(d_dim)*d_immersed_boundary_tagger_num_cells_buffer >= hier::IntVector::getZero(d_dim))
        {
            if (hier::IntVector::getOne(d_dim)*d_immersed_boundary_tagger_num_cells_buffer < num_cells_buffer_required)
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "d_immersed_boundary_tagger_num_cells_buffer is smaller than the required number!"
                    << std::endl);
            }
            else
            {
                num_cells_buffer_required = hier::IntVector::getOne(d_dim)*d_immersed_boundary_tagger_num_cells_buffer;
            }
        }
        
        d_immersed_boundary_tagger.reset(new ImmersedBoundaryTagger(
            "d_immersed_boundary_tagger",
            d_dim,
            d_grid_geometry,
            d_flow_model,
            num_cells_buffer_required));
    }
    else
    {
        d_immersed_boundary_tagger = nullptr;
    }
    
    /*
     * Initialize d_value_tagger.
     */
    
    if (d_value_tagger_db != nullptr)
    {
        d_value_tagger.reset(new ValueTagger(
            "d_value_tagger",
            d_dim,
            d_grid_geometry,
            d_flow_model,
            d_value_tagger_db));
    }
    else
    {
        d_value_tagger = nullptr;
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
            d_flow_model,
            d_gradient_tagger_db));
    }
    else
    {
        d_gradient_tagger = nullptr;
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
            d_flow_model,
            d_multiresolution_tagger_db));
    }
    else
    {
        d_multiresolution_tagger = nullptr;
    }
    
    /*
     * Initialize the side variable of convective flux.
     */
    
    d_variable_convective_flux = HAMERS_SHARED_PTR<pdat::SideVariable<double> > (
        new pdat::SideVariable<double>(dim, "convective flux", d_flow_model->getNumberOfEquations()));
    
    if (d_use_conservative_form_diffusive_flux)
    {
        /*
         * Initialize the side variable of diffusive flux.
         */
        
        d_variable_diffusive_flux = HAMERS_SHARED_PTR<pdat::SideVariable<double> > (
            new pdat::SideVariable<double>(dim, "diffusive flux", d_flow_model->getNumberOfEquations()));
    }
    else
    {
        /*
         * Initialize the cell variable of diffusive flux divergence.
         */
        
        d_variable_diffusive_flux_divergence = HAMERS_SHARED_PTR<pdat::CellVariable<double> > (
            new pdat::CellVariable<double>(dim, "diffusive flux divergence", d_flow_model->getNumberOfEquations()));
    }
    
    /*
     * Initialize the cell variable of source.
     */
    
    d_variable_source = HAMERS_SHARED_PTR<pdat::CellVariable<double> > (
        new pdat::CellVariable<double>(dim, "source", d_flow_model->getNumberOfEquations()));
}


NavierStokes::~NavierStokes()
{
    t_init.reset();
    t_compute_dt.reset();
    t_compute_fluxes_sources.reset();
    t_advance_step.reset();
    t_synchronize_fluxes.reset();
    t_setphysbcs.reset();
    t_tagimmersedbdry.reset();
    t_tagvalue.reset();
    t_taggradient.reset();
    t_tagmultiresolution.reset();
}


void
NavierStokes::registerModelVariables(
    RungeKuttaLevelIntegrator* integrator)
{
    TBOX_ASSERT(integrator != 0);
    
    /*
     * Determine the number of ghost cells needed.
     */
    
    hier::IntVector num_ghosts_intermediate = hier::IntVector::getZero(d_dim);
    
    num_ghosts_intermediate = hier::IntVector::max(
        num_ghosts_intermediate,
        d_convective_flux_reconstructor->getConvectiveFluxNumberOfGhostCells());
    
    hier::IntVector num_ghosts_diffusive_flux = hier::IntVector::getZero(d_dim);
    
    if (d_use_conservative_form_diffusive_flux)
    {
        num_ghosts_diffusive_flux = d_diffusive_flux_reconstructor->getDiffusiveFluxNumberOfGhostCells();
    }
    else
    {
        num_ghosts_diffusive_flux = d_nonconservative_diffusive_flux_divergence_operator->
            getNonconservativeDiffusiveFluxDivergenceOperatorNumberOfGhostCells();
    }
    
    bool use_transverse_derivatives_bc = d_Navier_Stokes_boundary_conditions->useTransverseDerivativesBoundaryConditions();
    
    // Increase the number of ghosts for computing transverse derivatives at corner ghost cells for diffusive/viscous flux.
    // Euler assumes corner ghost cells are not used.
    if (use_transverse_derivatives_bc)
    {
        hier::IntVector num_ghosts_transverse_derivatives_bc = d_Navier_Stokes_boundary_conditions->
            getBoundaryConditionsTransverseDerivativesNumberOfGhostCells();
        
        num_ghosts_diffusive_flux = num_ghosts_diffusive_flux + num_ghosts_transverse_derivatives_bc;
    }
    
    num_ghosts_intermediate = hier::IntVector::max(num_ghosts_intermediate, num_ghosts_diffusive_flux);
    
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
     * Set the number of immersed boundary ghost cells of d_immersed_boundaries.
     */
    
    if (d_use_immersed_boundaries)
    {
        d_immersed_boundaries->setNumberOfImmersedBoundaryGhosts(num_ghosts_intermediate);
        
        HAMERS_SHARED_PTR<FlowModelImmersedBoundaryMethod> flow_model_immersed_boundary_method =
            d_flow_model->getFlowModelImmersedBoundaryMethod();
        
        hier::IntVector num_ghosts_IB = flow_model_immersed_boundary_method->
            getImmersedBoundaryMethodAdditionalNumberOfGhostCells();
        
        num_ghosts_intermediate = num_ghosts_intermediate + num_ghosts_IB;
        num_ghosts              = num_ghosts              + num_ghosts_IB;
    }
    
    /*
     * Register the conservative variables of d_flow_model.
     */
    
    d_flow_model->registerConservativeVariables(
        integrator,
        num_ghosts,
        num_ghosts_intermediate);
    
    /*
     * Register the variables of flow model immersed boundary method.
     */
    
    if (d_use_immersed_boundaries)
    {
        HAMERS_SHARED_PTR<FlowModelImmersedBoundaryMethod> flow_model_immersed_boundary_method =
            d_flow_model->getFlowModelImmersedBoundaryMethod();
        
        flow_model_immersed_boundary_method->registerImmersedBoundaryMethodVariables(
            integrator,
            num_ghosts,
            num_ghosts_intermediate);
    }
    
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
    
    if (d_use_conservative_form_diffusive_flux)
    {
        integrator->registerVariable(
            d_variable_diffusive_flux,
            hier::IntVector::getZero(d_dim),
            hier::IntVector::getZero(d_dim),
            RungeKuttaLevelIntegrator::FLUX,
            d_grid_geometry,
            "CONSERVATIVE_COARSEN",
            "NO_REFINE");
    }
    else
    {
        integrator->registerVariable(
            d_variable_diffusive_flux_divergence,
            hier::IntVector::getZero(d_dim),
            hier::IntVector::getZero(d_dim),
            RungeKuttaLevelIntegrator::SOURCE,
            d_grid_geometry,
            "NO_COARSEN",
            "NO_REFINE");
    }
    
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
    
    HAMERS_SHARED_PTR<FlowModelStatisticsUtilities> flow_model_statistics_utilities =
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
        
        if (d_use_immersed_boundaries)
        {
            HAMERS_SHARED_PTR<FlowModelImmersedBoundaryMethod> flow_model_immersed_boundary_method =
                d_flow_model->getFlowModelImmersedBoundaryMethod();
            
            flow_model_immersed_boundary_method->registerPlotQuantities(
                d_visit_writer,
                integrator->getPlotContext());
        }
        
        
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
NavierStokes::setupLoadBalancer(
    RungeKuttaLevelIntegrator* integrator,
    mesh::GriddingAlgorithm* gridding_algorithm)
{
    NULL_USE(integrator);
    
    const hier::IntVector& zero_vec = hier::IntVector::getZero(d_dim);
    
    hier::VariableDatabase* vardb = hier::VariableDatabase::getDatabase();
    hier::PatchDataRestartManager* pdrm = hier::PatchDataRestartManager::getManager();
    
    if (d_use_nonuniform_workload && gridding_algorithm)
    {
        HAMERS_SHARED_PTR<mesh::TreeLoadBalancer> load_balancer(
            HAMERS_DYNAMIC_POINTER_CAST<mesh::TreeLoadBalancer, mesh::LoadBalanceStrategy>(
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
NavierStokes::initializeDataOnPatch(
    hier::Patch& patch,
    const double data_time,
    const bool initial_time)
{   
    t_init->start();
    
    d_flow_model->registerPatchWithDataContext(patch, getDataContext());
    
    std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > conservative_var_data =
        d_flow_model->getCellDataOfConservativeVariables();
    
    d_Navier_Stokes_initial_conditions->initializeDataOnPatch(
        patch,
        conservative_var_data,
        data_time,
        initial_time);
    
    if (d_use_immersed_boundaries)
    {
        /*
         * Initialize the immersed boundary method variables.
         */
        
        d_flow_model->setupImmersedBoundaryMethod();
        
        HAMERS_SHARED_PTR<FlowModelImmersedBoundaryMethod> flow_model_immersed_boundary_method =
            d_flow_model->getFlowModelImmersedBoundaryMethod();
        
        const hier::Box empty_box = hier::Box::getEmptyBox(d_dim);
        
        flow_model_immersed_boundary_method->setImmersedBoundaryMethodVariables(
            empty_box,
            data_time,
            initial_time,
            getDataContext());
        
        flow_model_immersed_boundary_method->setConservativeVariablesCellDataImmersedBoundaryGhosts(
            empty_box,
            data_time,
            initial_time,
            getDataContext());
    }
    
    d_flow_model->unregisterPatch();
    
    if (d_use_nonuniform_workload)
    {
        if (!patch.checkAllocated(d_workload_data_id))
        {
            patch.allocatePatchData(d_workload_data_id);
        }
        
        HAMERS_SHARED_PTR<pdat::CellData<double> > workload_data(
            HAMERS_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                patch.getPatchData(d_workload_data_id)));
        TBOX_ASSERT(workload_data);
        workload_data->fillAll(1.0);
    }

    t_init->stop();
}


std::vector<double>
NavierStokes::computeSpectralRadiusesAndStableDtOnPatch(
    hier::Patch& patch,
    const bool initial_time,
    const double dt_time)
{
    t_compute_dt->start();
    
    std::vector<double> spectral_radiuses_and_dt;
    
    const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
        HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(patch_geom);
#endif
    
    const double* dx = patch_geom->getDx();
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    /*
     * Register the patch, maximum wave speed and maximum diffusivity in the flow model and
     * compute the corresponding cell data.
     */
    
    d_flow_model->registerPatchWithDataContext(patch, getDataContext());
    
    hier::IntVector num_ghosts = d_flow_model->getNumberOfGhostCells();
    
    hier::Box ghost_box = interior_box;
    ghost_box.grow(num_ghosts);
    const hier::IntVector ghostcell_dims = ghost_box.numberCells();
    
    d_flow_model->setupSourceUtilities();
    
    HAMERS_SHARED_PTR<FlowModelSourceUtilities> source_utilities =
        d_flow_model->getFlowModelSourceUtilities();
    
    std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
    
    if (d_dim == tbox::Dimension(1))
    {
        num_subghosts_of_data.insert(
            std::pair<std::string, hier::IntVector>(
                "MAX_WAVE_SPEED_X", num_ghosts));
        num_subghosts_of_data.insert(
            std::pair<std::string, hier::IntVector>(
                "MAX_DIFFUSIVITY", num_ghosts));
    }
    else if (d_dim == tbox::Dimension(2))
    {
        num_subghosts_of_data.insert(
            std::pair<std::string, hier::IntVector>(
                "MAX_WAVE_SPEED_X", num_ghosts));
        num_subghosts_of_data.insert(
            std::pair<std::string, hier::IntVector>(
                "MAX_WAVE_SPEED_Y", num_ghosts));
        num_subghosts_of_data.insert(
            std::pair<std::string, hier::IntVector>(
                "MAX_DIFFUSIVITY", num_ghosts));
    }
    else if (d_dim == tbox::Dimension(3))
    {
        num_subghosts_of_data.insert(
            std::pair<std::string, hier::IntVector>(
                "MAX_WAVE_SPEED_X", num_ghosts));
        num_subghosts_of_data.insert(
            std::pair<std::string, hier::IntVector>(
                "MAX_WAVE_SPEED_Y", num_ghosts));
        num_subghosts_of_data.insert(
            std::pair<std::string, hier::IntVector>(
                "MAX_WAVE_SPEED_Z", num_ghosts));
        num_subghosts_of_data.insert(
            std::pair<std::string, hier::IntVector>(
                "MAX_DIFFUSIVITY", num_ghosts));
    }
    
    d_flow_model->registerDerivedVariables(num_subghosts_of_data);
    
    if (source_utilities->hasSourceTerms())
    {
        source_utilities->registerDerivedVariablesForSourceTermsStableDt(hier::IntVector::getZero(d_dim));
    }
    
    d_flow_model->allocateMemoryForDerivedCellData();
    
    if (source_utilities->hasSourceTerms())
    {
        source_utilities->allocateMemoryForDerivedCellData();
    }
    
    d_flow_model->computeDerivedCellData();
    
    if (source_utilities->hasSourceTerms())
    {
        source_utilities->computeDerivedCellData();
    }
    
    /*
     * Get the pointer to the cell data of the immersed boundary mask.
     * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
     */
    
    int* IB_mask = nullptr;
    HAMERS_SHARED_PTR<pdat::CellData<int> > IB_mask_cell_data;
    hier::IntVector num_ghosts_IB_mask(d_dim);
    hier::IntVector ghostcell_dims_IB_mask(d_dim);
    const int fluid = int(IB_MASK::FLUID);
    
    if (d_use_immersed_boundaries)
    {
        d_flow_model->setupImmersedBoundaryMethod();
        
        HAMERS_SHARED_PTR<FlowModelImmersedBoundaryMethod> flow_model_immersed_boundary_method =
            d_flow_model->getFlowModelImmersedBoundaryMethod();
        
        IB_mask_cell_data = flow_model_immersed_boundary_method->
            getCellDataOfImmersedBoundaryMask(getDataContext());
        
        IB_mask                = IB_mask_cell_data->getPointer(0);
        num_ghosts_IB_mask     = IB_mask_cell_data->getGhostCellWidth();
        ghostcell_dims_IB_mask = IB_mask_cell_data->getGhostBox().numberCells();
    }
    
    if (d_dim == tbox::Dimension(1))
    {
        spectral_radiuses_and_dt.resize(2, double(0));
        
        /*
         * Get the dimension and grid spacing.
         */
        
        const int interior_dim_0 = interior_dims[0];
        
        const double dx_0 = dx[0];
        
        /*
         * Get the pointers to the maximum wave speed and maximum diffusivity inside the flow model.
         * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
         */
        
        HAMERS_SHARED_PTR<pdat::CellData<double> > max_wave_speed_x = d_flow_model->getCellData("MAX_WAVE_SPEED_X");
        HAMERS_SHARED_PTR<pdat::CellData<double> > max_diffusivity = d_flow_model->getCellData("MAX_DIFFUSIVITY");
        
        hier::IntVector num_subghosts_max_wave_speed_x = max_wave_speed_x->getGhostCellWidth();
        hier::IntVector num_subghosts_max_diffusivity = max_diffusivity->getGhostCellWidth();
        
        TBOX_ASSERT(num_subghosts_max_wave_speed_x == num_ghosts);
        TBOX_ASSERT(num_subghosts_max_diffusivity == num_ghosts);
        
        const int num_ghosts_0 = num_ghosts[0];
        
        double* max_lambda_x = max_wave_speed_x->getPointer(0);
        double* max_D = max_diffusivity->getPointer(0);
        
        double spectral_radiuses_and_dt_0 = double(0);
        double spectral_radiuses_and_dt_1 = double(0);
        
        if (d_use_immersed_boundaries)
        {
            const int num_ghosts_0_IB_mask = num_ghosts_IB_mask[0];
            
            // HAMERS_PRAGMA_VEC("omp simd reduction(max: spectral_radiuses_and_dt_0, spectral_radiuses_and_dt_1")
            for (int i = -num_ghosts_0;
                 i < interior_dim_0 + num_ghosts_0;
                 i++)
            {
                // Compute the linear indices.
                const int idx = i + num_ghosts_0;
                const int idx_IB_mask = i + num_ghosts_0_IB_mask;
                
                if (IB_mask[idx_IB_mask] == fluid)
                {
                    const double spectral_radius_acoustic_x = max_lambda_x[idx]/dx_0;
                    
                    spectral_radiuses_and_dt_0 = fmax(spectral_radiuses_and_dt_0, spectral_radius_acoustic_x);
                    
                    spectral_radiuses_and_dt_1 = fmax(spectral_radiuses_and_dt_1, spectral_radius_acoustic_x);
                }
            }
        }
        else
        {
            // HAMERS_PRAGMA_VEC("omp simd reduction(max: spectral_radiuses_and_dt_0, spectral_radiuses_and_dt_1")
            for (int i = -num_ghosts_0;
                 i < interior_dim_0 + num_ghosts_0;
                 i++)
            {
                // Compute the linear index.
                const int idx = i + num_ghosts_0;
                
                const double spectral_radius_acoustic_x = max_lambda_x[idx]/dx_0;
                
                spectral_radiuses_and_dt_0 = fmax(spectral_radiuses_and_dt_0, spectral_radius_acoustic_x);
                
                spectral_radiuses_and_dt_1 = fmax(spectral_radiuses_and_dt_1, spectral_radius_acoustic_x);
            }
        }
        
        spectral_radiuses_and_dt[0] = spectral_radiuses_and_dt_0;
        spectral_radiuses_and_dt[1] = spectral_radiuses_and_dt_1;
        
        double spectral_radius_tmp = double(0);
        
        if (d_use_immersed_boundaries)
        {
            const int num_ghosts_0_IB_mask = num_ghosts_IB_mask[0];
            
            // HAMERS_PRAGMA_VEC("omp simd reduction(max: spectral_radius_tmp")
            for (int i = -num_ghosts_0;
                 i < interior_dim_0 + num_ghosts_0;
                 i++)
            {
                // Compute the linear indices.
                const int idx = i + num_ghosts_0;
                const int idx_IB_mask = i + num_ghosts_0_IB_mask;
                
                if (IB_mask[idx_IB_mask] == fluid)
                {
                    const double spectral_radius_diffusive = double(2)*max_D[idx]/(dx_0*dx_0);
                    
                    spectral_radius_tmp = fmax(spectral_radius_tmp, spectral_radius_diffusive);
                }
            }
        }
        else
        {
            // HAMERS_PRAGMA_VEC("omp simd reduction(max: spectral_radius_tmp")
            for (int i = -num_ghosts_0;
                 i < interior_dim_0 + num_ghosts_0;
                 i++)
            {
                // Compute the linear index.
                const int idx = i + num_ghosts_0;
                
                const double spectral_radius_diffusive = double(2)*max_D[idx]/(dx_0*dx_0);
                
                spectral_radius_tmp = fmax(spectral_radius_tmp, spectral_radius_diffusive);
            }
        }
        
        spectral_radiuses_and_dt[1] = fmax(spectral_radius_tmp, spectral_radiuses_and_dt[1]);
        
        spectral_radiuses_and_dt[1] = double(1)/(spectral_radiuses_and_dt[1] + HAMERS_EPSILON);
        
        if (source_utilities->hasSourceTerms())
        {
            const double stable_dt_source = source_utilities->getStableDtOnPatch();
            
            spectral_radiuses_and_dt[1] = fmin(stable_dt_source, spectral_radiuses_and_dt[1]);
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        spectral_radiuses_and_dt.resize(3, double(0));
        
        /*
         * Get the dimensions and grid spacings.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        const double dx_0 = dx[0];
        const double dx_1 = dx[1];
        
        /*
         * Get the pointers to the maximum wave speeds and maximum diffusivity inside the flow model.
         * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
         */
        
        HAMERS_SHARED_PTR<pdat::CellData<double> > max_wave_speed_x = d_flow_model->getCellData("MAX_WAVE_SPEED_X");
        HAMERS_SHARED_PTR<pdat::CellData<double> > max_wave_speed_y = d_flow_model->getCellData("MAX_WAVE_SPEED_Y");
        HAMERS_SHARED_PTR<pdat::CellData<double> > max_diffusivity = d_flow_model->getCellData("MAX_DIFFUSIVITY");
        
        hier::IntVector num_subghosts_max_wave_speed_x = max_wave_speed_x->getGhostCellWidth();
        hier::IntVector num_subghosts_max_wave_speed_y = max_wave_speed_y->getGhostCellWidth();
        hier::IntVector num_subghosts_max_diffusivity = max_diffusivity->getGhostCellWidth();
        
        TBOX_ASSERT(num_subghosts_max_wave_speed_x == num_ghosts);
        TBOX_ASSERT(num_subghosts_max_wave_speed_y == num_ghosts);
        TBOX_ASSERT(num_subghosts_max_diffusivity == num_ghosts);
        
        const int num_ghosts_0 = num_ghosts[0];
        const int num_ghosts_1 = num_ghosts[1];
        const int ghostcell_dim_0 = ghostcell_dims[0];
        
        double* max_lambda_x = max_wave_speed_x->getPointer(0);
        double* max_lambda_y = max_wave_speed_y->getPointer(0);
        double* max_D = max_diffusivity->getPointer(0);
        
        double spectral_radiuses_and_dt_0 = double(0);
        double spectral_radiuses_and_dt_1 = double(0);
        double spectral_radiuses_and_dt_2 = double(0);
        
        if (d_use_immersed_boundaries)
        {
            const int num_ghosts_0_IB_mask = num_ghosts_IB_mask[0];
            const int num_ghosts_1_IB_mask = num_ghosts_IB_mask[1];
            const int ghostcell_dim_0_IB_mask = ghostcell_dims_IB_mask[0];
            
            for (int j = -num_ghosts_1;
                 j < interior_dim_1 + num_ghosts_1;
                 j++)
            {
                // HAMERS_PRAGMA_VEC("omp simd reduction(max: spectral_radiuses_and_dt_0, spectral_radiuses_and_dt_1, spectral_radiuses_and_dt_2")
                for (int i = -num_ghosts_0;
                     i < interior_dim_0 + num_ghosts_0;
                     i++)
                {
                    // Compute the linear indices.
                    const int idx = (i + num_ghosts_0) +
                        (j + num_ghosts_1)*ghostcell_dim_0;
                    
                    const int idx_IB_mask = (i + num_ghosts_0_IB_mask) +
                        (j + num_ghosts_1_IB_mask)*ghostcell_dim_0_IB_mask;
                    
                    if (IB_mask[idx_IB_mask] == fluid)
                    {
                        const double spectral_radius_acoustic_x = max_lambda_x[idx]/dx_0;
                        const double spectral_radius_acoustic_y = max_lambda_y[idx]/dx_1;
                        
                        spectral_radiuses_and_dt_0 = fmax(spectral_radiuses_and_dt_0, spectral_radius_acoustic_x);
                        spectral_radiuses_and_dt_1 = fmax(spectral_radiuses_and_dt_1, spectral_radius_acoustic_y);
                        
                        spectral_radiuses_and_dt_2 = fmax(spectral_radiuses_and_dt_2,
                            spectral_radius_acoustic_x + spectral_radius_acoustic_y);
                    }
                }
            }
        }
        else
        {
            for (int j = -num_ghosts_1;
                 j < interior_dim_1 + num_ghosts_1;
                 j++)
            {
                // HAMERS_PRAGMA_VEC("omp simd reduction(max: spectral_radiuses_and_dt_0, spectral_radiuses_and_dt_1, spectral_radiuses_and_dt_2")
                for (int i = -num_ghosts_0;
                     i < interior_dim_0 + num_ghosts_0;
                     i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0) +
                        (j + num_ghosts_1)*ghostcell_dim_0;
                    
                    const double spectral_radius_acoustic_x = max_lambda_x[idx]/dx_0;
                    const double spectral_radius_acoustic_y = max_lambda_y[idx]/dx_1;
                    
                    spectral_radiuses_and_dt_0 = fmax(spectral_radiuses_and_dt_0, spectral_radius_acoustic_x);
                    spectral_radiuses_and_dt_1 = fmax(spectral_radiuses_and_dt_1, spectral_radius_acoustic_y);
                    
                    spectral_radiuses_and_dt_2 = fmax(spectral_radiuses_and_dt_2,
                        spectral_radius_acoustic_x + spectral_radius_acoustic_y);
                }
            }
        }
        
        spectral_radiuses_and_dt[0] = spectral_radiuses_and_dt_0;
        spectral_radiuses_and_dt[1] = spectral_radiuses_and_dt_1;
        spectral_radiuses_and_dt[2] = spectral_radiuses_and_dt_2;
        
        double spectral_radius_tmp = double(0);
        
        if (d_use_immersed_boundaries)
        {
            const int num_ghosts_0_IB_mask = num_ghosts_IB_mask[0];
            const int num_ghosts_1_IB_mask = num_ghosts_IB_mask[1];
            const int ghostcell_dim_0_IB_mask = ghostcell_dims_IB_mask[0];
            
            for (int j = -num_ghosts_1;
                 j < interior_dim_1 + num_ghosts_1;
                 j++)
            {
                // HAMERS_PRAGMA_VEC("omp simd reduction(max: spectral_radius_tmp")
                for (int i = -num_ghosts_0;
                     i < interior_dim_0 + num_ghosts_0;
                     i++)
                {
                    // Compute the linear indices.
                    const int idx = (i + num_ghosts_0) +
                        (j + num_ghosts_1)*ghostcell_dim_0;
                    
                    const int idx_IB_mask = (i + num_ghosts_0_IB_mask) +
                        (j + num_ghosts_1_IB_mask)*ghostcell_dim_0_IB_mask;
                    
                    if (IB_mask[idx_IB_mask] == fluid)
                    {
                        const double spectral_radius_diffusive = double(2)*fmax(
                            max_D[idx]/(dx_0*dx_0),
                            max_D[idx]/(dx_1*dx_1));
                        
                        spectral_radius_tmp = fmax(spectral_radius_tmp, spectral_radius_diffusive);
                    }
                }
            }
        }
        else
        {
            for (int j = -num_ghosts_1;
                 j < interior_dim_1 + num_ghosts_1;
                 j++)
            {
                // HAMERS_PRAGMA_VEC("omp simd reduction(max: spectral_radius_tmp")
                for (int i = -num_ghosts_0;
                     i < interior_dim_0 + num_ghosts_0;
                     i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0) +
                        (j + num_ghosts_1)*ghostcell_dim_0;
                    
                    const double spectral_radius_diffusive = double(2)*fmax(
                        max_D[idx]/(dx_0*dx_0),
                        max_D[idx]/(dx_1*dx_1));
                    
                    spectral_radius_tmp = fmax(spectral_radius_tmp, spectral_radius_diffusive);
                }
            }
        }
        
        spectral_radiuses_and_dt[2] = fmax(spectral_radius_tmp, spectral_radiuses_and_dt[2]);
        
        spectral_radiuses_and_dt[2] = double(1)/(spectral_radiuses_and_dt[2] + HAMERS_EPSILON);
        
        if (source_utilities->hasSourceTerms())
        {
            const double stable_dt_source = source_utilities->getStableDtOnPatch();
            
            spectral_radiuses_and_dt[2] = fmin(stable_dt_source, spectral_radiuses_and_dt[2]);
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        spectral_radiuses_and_dt.resize(4, double(0));
        
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
         * Get the pointers to the maximum wave speeds and maximum diffusivity inside the flow model.
         * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
         */
        
        HAMERS_SHARED_PTR<pdat::CellData<double> > max_wave_speed_x = d_flow_model->getCellData("MAX_WAVE_SPEED_X");
        HAMERS_SHARED_PTR<pdat::CellData<double> > max_wave_speed_y = d_flow_model->getCellData("MAX_WAVE_SPEED_Y");
        HAMERS_SHARED_PTR<pdat::CellData<double> > max_wave_speed_z = d_flow_model->getCellData("MAX_WAVE_SPEED_Z");
        HAMERS_SHARED_PTR<pdat::CellData<double> > max_diffusivity = d_flow_model->getCellData("MAX_DIFFUSIVITY");
        
        hier::IntVector num_subghosts_max_wave_speed_x = max_wave_speed_x->getGhostCellWidth();
        hier::IntVector num_subghosts_max_wave_speed_y = max_wave_speed_y->getGhostCellWidth();
        hier::IntVector num_subghosts_max_wave_speed_z = max_wave_speed_z->getGhostCellWidth();
        hier::IntVector num_subghosts_max_diffusivity = max_diffusivity->getGhostCellWidth();
        
        TBOX_ASSERT(num_subghosts_max_wave_speed_x == num_ghosts);
        TBOX_ASSERT(num_subghosts_max_wave_speed_y == num_ghosts);
        TBOX_ASSERT(num_subghosts_max_wave_speed_z == num_ghosts);
        TBOX_ASSERT(num_subghosts_max_diffusivity == num_ghosts);
        
        const int num_ghosts_0 = num_ghosts[0];
        const int num_ghosts_1 = num_ghosts[1];
        const int num_ghosts_2 = num_ghosts[2];
        const int ghostcell_dim_0 = ghostcell_dims[0];
        const int ghostcell_dim_1 = ghostcell_dims[1];
        
        double* max_lambda_x = max_wave_speed_x->getPointer(0);
        double* max_lambda_y = max_wave_speed_y->getPointer(0);
        double* max_lambda_z = max_wave_speed_z->getPointer(0);
        double* max_D = max_diffusivity->getPointer(0);
        
        double spectral_radiuses_and_dt_0 = double(0);
        double spectral_radiuses_and_dt_1 = double(0);
        double spectral_radiuses_and_dt_2 = double(0);
        double spectral_radiuses_and_dt_3 = double(0);
        
        if (d_use_immersed_boundaries)
        {
            const int num_ghosts_0_IB_mask = num_ghosts_IB_mask[0];
            const int num_ghosts_1_IB_mask = num_ghosts_IB_mask[1];
            const int num_ghosts_2_IB_mask = num_ghosts_IB_mask[2];
            const int ghostcell_dim_0_IB_mask = ghostcell_dims_IB_mask[0];
            const int ghostcell_dim_1_IB_mask = ghostcell_dims_IB_mask[1];
            
            for (int k = -num_ghosts_2;
                 k < interior_dim_2 + num_ghosts_2;
                 k++)
            {
                for (int j = -num_ghosts_1;
                     j < interior_dim_1 + num_ghosts_1;
                     j++)
                {
                    // HAMERS_PRAGMA_VEC("omp simd reduction(max: spectral_radiuses_and_dt_0, spectral_radiuses_and_dt_1, spectral_radiuses_and_dt_2, spectral_radiuses_and_dt_3")
                    for (int i = -num_ghosts_0;
                         i < interior_dim_0 + num_ghosts_0;
                         i++)
                    {
                        // Compute the linear indices.
                        const int idx = (i + num_ghosts_0) +
                            (j + num_ghosts_1)*ghostcell_dim_0 +
                            (k + num_ghosts_2)*ghostcell_dim_0*
                                ghostcell_dim_1;
                        
                        const int idx_IB_mask = (i + num_ghosts_0_IB_mask) +
                            (j + num_ghosts_1_IB_mask)*ghostcell_dim_0_IB_mask +
                            (k + num_ghosts_2_IB_mask)*ghostcell_dim_0_IB_mask*
                                ghostcell_dim_1_IB_mask;
                        
                        if (IB_mask[idx_IB_mask] == fluid)
                        {
                            const double spectral_radius_acoustic_x = max_lambda_x[idx]/dx_0;
                            const double spectral_radius_acoustic_y = max_lambda_y[idx]/dx_1;
                            const double spectral_radius_acoustic_z = max_lambda_z[idx]/dx_2;
                            
                            spectral_radiuses_and_dt_0 = fmax(spectral_radiuses_and_dt_0, spectral_radius_acoustic_x);
                            spectral_radiuses_and_dt_1 = fmax(spectral_radiuses_and_dt_1, spectral_radius_acoustic_y);
                            spectral_radiuses_and_dt_2 = fmax(spectral_radiuses_and_dt_2, spectral_radius_acoustic_z);
                        
                            spectral_radiuses_and_dt_3 = fmax(spectral_radiuses_and_dt_3,
                                spectral_radius_acoustic_x + spectral_radius_acoustic_y + spectral_radius_acoustic_z);
                        }
                    }
                }
            }
        }
        else
        {
            for (int k = -num_ghosts_2;
                 k < interior_dim_2 + num_ghosts_2;
                 k++)
            {
                for (int j = -num_ghosts_1;
                     j < interior_dim_1 + num_ghosts_1;
                     j++)
                {
                    // HAMERS_PRAGMA_VEC("omp simd reduction(max: spectral_radiuses_and_dt_0, spectral_radiuses_and_dt_1, spectral_radiuses_and_dt_2, spectral_radiuses_and_dt_3")
                    for (int i = -num_ghosts_0;
                         i < interior_dim_0 + num_ghosts_0;
                         i++)
                    {
                        // Compute the linear indices.
                        const int idx = (i + num_ghosts_0) +
                            (j + num_ghosts_1)*ghostcell_dim_0 +
                            (k + num_ghosts_2)*ghostcell_dim_0*
                                ghostcell_dim_1;
                        
                        const double spectral_radius_acoustic_x = max_lambda_x[idx]/dx_0;
                        const double spectral_radius_acoustic_y = max_lambda_y[idx]/dx_1;
                        const double spectral_radius_acoustic_z = max_lambda_z[idx]/dx_2;
                        
                        spectral_radiuses_and_dt_0 = fmax(spectral_radiuses_and_dt_0, spectral_radius_acoustic_x);
                        spectral_radiuses_and_dt_1 = fmax(spectral_radiuses_and_dt_1, spectral_radius_acoustic_y);
                        spectral_radiuses_and_dt_2 = fmax(spectral_radiuses_and_dt_2, spectral_radius_acoustic_z);
                    
                        spectral_radiuses_and_dt_3 = fmax(spectral_radiuses_and_dt_3,
                            spectral_radius_acoustic_x + spectral_radius_acoustic_y + spectral_radius_acoustic_z);
                    }
                }
            }
        }
        
        spectral_radiuses_and_dt[0] = spectral_radiuses_and_dt_0;
        spectral_radiuses_and_dt[1] = spectral_radiuses_and_dt_1;
        spectral_radiuses_and_dt[2] = spectral_radiuses_and_dt_2;
        spectral_radiuses_and_dt[3] = spectral_radiuses_and_dt_3;
        
        double spectral_radius_tmp = double(0);
        
        if (d_use_immersed_boundaries)
        {
            const int num_ghosts_0_IB_mask = num_ghosts_IB_mask[0];
            const int num_ghosts_1_IB_mask = num_ghosts_IB_mask[1];
            const int num_ghosts_2_IB_mask = num_ghosts_IB_mask[2];
            const int ghostcell_dim_0_IB_mask = ghostcell_dims_IB_mask[0];
            const int ghostcell_dim_1_IB_mask = ghostcell_dims_IB_mask[1];
            
            for (int k = -num_ghosts_2;
                 k < interior_dim_2 + num_ghosts_2;
                 k++)
            {
                for (int j = -num_ghosts_1;
                     j < interior_dim_1 + num_ghosts_1;
                     j++)
                {
                    // HAMERS_PRAGMA_VEC("omp simd reduction(max: spectral_radius_tmp")
                    for (int i = -num_ghosts_0;
                         i < interior_dim_0 + num_ghosts_0;
                         i++)
                    {
                        // Compute the linear indices.
                        const int idx = (i + num_ghosts_0) +
                            (j + num_ghosts_1)*ghostcell_dim_0 +
                            (k + num_ghosts_2)*ghostcell_dim_0*
                                ghostcell_dim_1;
                        
                        const int idx_IB_mask = (i + num_ghosts_0_IB_mask) +
                            (j + num_ghosts_1_IB_mask)*ghostcell_dim_0_IB_mask +
                            (k + num_ghosts_2_IB_mask)*ghostcell_dim_0_IB_mask*
                                ghostcell_dim_1_IB_mask;
                        
                        if (IB_mask[idx_IB_mask] == fluid)
                        {
                            const double spectral_radius_diffusive = double(2)*fmax(
                                max_D[idx]/(dx_0*dx_0),
                                fmax(max_D[idx]/(dx_1*dx_1),
                                    max_D[idx]/(dx_2*dx_2)));
                            
                            spectral_radius_tmp = fmax(spectral_radius_tmp, spectral_radius_diffusive);
                        }
                    }
                }
            }
        }
        else
        {
            for (int k = -num_ghosts_2;
                 k < interior_dim_2 + num_ghosts_2;
                 k++)
            {
                for (int j = -num_ghosts_1;
                     j < interior_dim_1 + num_ghosts_1;
                     j++)
                {
                    // HAMERS_PRAGMA_VEC("max: spectral_radius_tmp")
                    for (int i = -num_ghosts_0;
                         i < interior_dim_0 + num_ghosts_0;
                         i++)
                    {
                        // Compute the linear index.
                        const int idx = (i + num_ghosts_0) +
                            (j + num_ghosts_1)*ghostcell_dim_0 +
                            (k + num_ghosts_2)*ghostcell_dim_0*
                                ghostcell_dim_1;
                        
                        const double spectral_radius_diffusive = double(2)*fmax(
                            max_D[idx]/(dx_0*dx_0),
                            fmax(max_D[idx]/(dx_1*dx_1),
                                max_D[idx]/(dx_2*dx_2)));
                        
                        spectral_radius_tmp = fmax(spectral_radius_tmp, spectral_radius_diffusive);
                    }
                }
            }
        }
        
        spectral_radiuses_and_dt[3] = fmax(spectral_radius_tmp, spectral_radiuses_and_dt[3]);
        
        spectral_radiuses_and_dt[3] = double(1)/(spectral_radiuses_and_dt[3] + HAMERS_EPSILON);;
        
        if (source_utilities->hasSourceTerms())
        {
            const double stable_dt_source = source_utilities->getStableDtOnPatch();
            
            spectral_radiuses_and_dt[3] = fmin(stable_dt_source, spectral_radiuses_and_dt[3]);
        }
    }
    
    /*
     * Unregister the patch and data of all registered derived cell variables in the flow model.
     */
    
    d_flow_model->unregisterPatch();
    
    t_compute_dt->stop();
    
    return spectral_radiuses_and_dt;
}


/**
 * Set the immersed boundary ghost cells.
 */
void
NavierStokes::setImmersedBoundaryGhostCells(
    hier::Patch& patch,
    const double time,
    const int RK_step_number,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    if (!d_use_immersed_boundaries || !d_use_ghost_cell_immersed_boundary_method)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "d_use_immersed_boundaries or d_use_ghost_cell_immersed_boundary_method are set to false!"
            << std::endl);
    }
    
    d_flow_model->registerPatchWithDataContext(patch, data_context);
    
    /*
     * Initialize the immersed boundary method variables.
     */
    
    d_flow_model->setupImmersedBoundaryMethod();
    
    HAMERS_SHARED_PTR<FlowModelImmersedBoundaryMethod> flow_model_immersed_boundary_method =
        d_flow_model->getFlowModelImmersedBoundaryMethod();
    
    const hier::Box empty_box = hier::Box::getEmptyBox(d_dim);
    
    if (RK_step_number == 0 || !d_use_static_immersed_boundaries)
    {
        flow_model_immersed_boundary_method->setImmersedBoundaryMethodVariables(
            empty_box,
            time,
            false,
            getDataContext());
    }
    
    // Compute the immersed boundary ghost cells here.
    
    flow_model_immersed_boundary_method->setConservativeVariablesCellDataImmersedBoundaryGhosts(
        empty_box,
        time,
        false,
        getDataContext());
    
    d_flow_model->unregisterPatch();
}


void
NavierStokes::computeFluxesAndSourcesOnPatch(
    hier::Patch& patch,
    const double time,
    const double dt,
    const int RK_step_number,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    t_compute_fluxes_sources->start();
    
    /*
     * Set zero for the source.
     */
    
    if (data_context)
    {
        HAMERS_SHARED_PTR<pdat::CellData<double> > data_source(
            HAMERS_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                patch.getPatchData(d_variable_source, data_context)));
        
        data_source->fillAll(0.0);
    }
    else
    {
        HAMERS_SHARED_PTR<pdat::CellData<double> > data_source(
            HAMERS_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                patch.getPatchData(d_variable_source, getDataContext())));
        
        data_source->fillAll(0.0);
    }
    
    /*
     * Compute the convective flux, source due to splitting of convective term and diffusive flux.
     */
    
    if (data_context)
    {
        d_convective_flux_reconstructor->computeConvectiveFluxAndSourceOnPatch(
            patch,
            d_variable_convective_flux,
            d_variable_source,
            data_context,
            time,
            dt,
            RK_step_number);
        
        if (d_use_conservative_form_diffusive_flux)
        {
            d_diffusive_flux_reconstructor->computeDiffusiveFluxOnPatch(
                patch,
                d_variable_diffusive_flux,
                data_context,
                time,
                dt,
                RK_step_number);
        }
        else
        {
            d_nonconservative_diffusive_flux_divergence_operator->
                computeNonconservativeDiffusiveFluxDivergenceOnPatch(
                    patch,
                    d_variable_diffusive_flux_divergence,
                    data_context,
                    time,
                    dt,
                    RK_step_number);
        }
    }
    else
    {
        d_convective_flux_reconstructor->computeConvectiveFluxAndSourceOnPatch(
            patch,
            d_variable_convective_flux,
            d_variable_source,
            getDataContext(),
            time,
            dt,
            RK_step_number);
        
        if (d_use_conservative_form_diffusive_flux)
        {
            d_diffusive_flux_reconstructor->computeDiffusiveFluxOnPatch(
                patch,
                d_variable_diffusive_flux,
                getDataContext(),
                time,
                dt,
                RK_step_number);
        }
        else
        {
            d_nonconservative_diffusive_flux_divergence_operator->
                computeNonconservativeDiffusiveFluxDivergenceOnPatch(
                    patch,
                    d_variable_diffusive_flux_divergence,
                    getDataContext(),
                    time,
                    dt,
                    RK_step_number);
        }
    }
    
    /*
     * Compute the source terms.
     */
    
    d_flow_model->setupSourceUtilities();
    
    HAMERS_SHARED_PTR<FlowModelSourceUtilities> source_utilities =
        d_flow_model->getFlowModelSourceUtilities();
    
    if (source_utilities->hasSourceTerms())
    {
        if (data_context)
        {
            d_flow_model->registerPatchWithDataContext(patch, data_context);
        }
        else
        {
            d_flow_model->registerPatchWithDataContext(patch, getDataContext());
        }
        
        source_utilities->registerDerivedVariablesForSourceTerms(hier::IntVector::getZero(d_dim));
        
        source_utilities->allocateMemoryForDerivedCellData();
        
        source_utilities->computeDerivedCellData();
        
        source_utilities->computeSourceTermsOnPatch(
            d_variable_source,
            time,
            dt,
            RK_step_number);
        
        d_flow_model->unregisterPatch();
    }
    
    t_compute_fluxes_sources->stop();
}


void
NavierStokes::advanceSingleStepOnPatch(
    hier::Patch& patch,
    const double time,
    const double dt,
    const std::vector<double>& alpha,
    const std::vector<double>& beta,
    const std::vector<double>& gamma,
    const std::vector<HAMERS_SHARED_PTR<hier::VariableContext> >& intermediate_context)
{
    NULL_USE(time);
    NULL_USE(dt);
    
    t_advance_step->start();
    
    const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
        HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
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
    
    std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > conservative_variables =
        d_flow_model->getCellDataOfConservativeVariables();
    
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
    
    /*
     * Get the pointer to the cell data of the immersed boundary mask.
     * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
     */
    
    int* IB_mask = nullptr;
    HAMERS_SHARED_PTR<pdat::CellData<int> > IB_mask_cell_data;
    hier::IntVector num_ghosts_IB_mask(d_dim);
    hier::IntVector ghostcell_dims_IB_mask(d_dim);
    const int fluid = int(IB_MASK::FLUID);
    
    if (d_use_immersed_boundaries)
    {
        d_flow_model->setupImmersedBoundaryMethod();
        
        HAMERS_SHARED_PTR<FlowModelImmersedBoundaryMethod> flow_model_immersed_boundary_method =
            d_flow_model->getFlowModelImmersedBoundaryMethod();
        
        IB_mask_cell_data = flow_model_immersed_boundary_method->
            getCellDataOfImmersedBoundaryMask(getDataContext());
        
        IB_mask                = IB_mask_cell_data->getPointer(0);
        num_ghosts_IB_mask     = IB_mask_cell_data->getGhostCellWidth();
        ghostcell_dims_IB_mask = IB_mask_cell_data->getGhostBox().numberCells();
    }
    
    if (d_use_immersed_boundaries)
    {
        d_flow_model->fillCellDataOfConservativeVariablesWithZero(
            IB_mask_cell_data,
            fluid);
    }
    else
    {
        d_flow_model->fillCellDataOfConservativeVariablesWithZero();
    }
    
    // Unregister the patch.
    d_flow_model->unregisterPatch();
    
    /*
     * Use alpha, beta and gamma values to update the time-dependent solution,
     * flux and source
     */
    
    HAMERS_SHARED_PTR<pdat::SideData<double> > convective_flux(
        HAMERS_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
            patch.getPatchData(d_variable_convective_flux, getDataContext())));
    
    HAMERS_SHARED_PTR<pdat::SideData<double> > diffusive_flux;
    HAMERS_SHARED_PTR<pdat::CellData<double> > diffusive_flux_divergence;
    
    if (d_use_conservative_form_diffusive_flux)
    {
        diffusive_flux =
            HAMERS_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                patch.getPatchData(d_variable_diffusive_flux, getDataContext()));
    }
    else
    {
        diffusive_flux_divergence =
            HAMERS_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                patch.getPatchData(d_variable_diffusive_flux_divergence, getDataContext()));
    }
    
    HAMERS_SHARED_PTR<pdat::CellData<double> > source(
        HAMERS_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
            patch.getPatchData(d_variable_source, getDataContext())));
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(convective_flux);
    if (d_use_conservative_form_diffusive_flux)
    {
        TBOX_ASSERT(diffusive_flux);
    }
    else
    {
        TBOX_ASSERT(diffusive_flux_divergence);
    }
    TBOX_ASSERT(source);
    
    TBOX_ASSERT(convective_flux->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
    if (d_use_conservative_form_diffusive_flux)
    {
        TBOX_ASSERT(diffusive_flux->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
    }
    else
    {
        TBOX_ASSERT(diffusive_flux_divergence->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
    }
    TBOX_ASSERT(source->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
#endif
    
    int num_coeffs = static_cast<int>(alpha.size());
    
    for (int n = 0; n < num_coeffs; n++)
    {
        HAMERS_SHARED_PTR<pdat::SideData<double> > convective_flux_intermediate(
            HAMERS_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                    patch.getPatchData(d_variable_convective_flux, intermediate_context[n])));
        
        HAMERS_SHARED_PTR<pdat::SideData<double> > diffusive_flux_intermediate;
        HAMERS_SHARED_PTR<pdat::CellData<double> > diffusive_flux_divergence_intermediate;
        
        if (d_use_conservative_form_diffusive_flux)
        {
            diffusive_flux_intermediate =
                HAMERS_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                        patch.getPatchData(d_variable_diffusive_flux, intermediate_context[n]));
        }
        else
        {
            diffusive_flux_divergence_intermediate =
                HAMERS_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_variable_diffusive_flux_divergence, intermediate_context[n]));
        }
        
        HAMERS_SHARED_PTR<pdat::CellData<double> > source_intermediate(
            HAMERS_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                patch.getPatchData(d_variable_source, intermediate_context[n])));
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(convective_flux_intermediate);
        if (d_use_conservative_form_diffusive_flux)
        {
            TBOX_ASSERT(diffusive_flux_intermediate);
        }
        else
        {
            TBOX_ASSERT(diffusive_flux_divergence_intermediate);
        }
        TBOX_ASSERT(source_intermediate);
        
        TBOX_ASSERT(convective_flux_intermediate->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
        if (d_use_conservative_form_diffusive_flux)
        {
            TBOX_ASSERT(diffusive_flux_intermediate->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
        }
        else
        {
            TBOX_ASSERT(diffusive_flux_divergence_intermediate->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
        }
        TBOX_ASSERT(source_intermediate->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
#endif
        
        /*
         * Create a vector pointers to the time-dependent variables for the
         * current intermediate data context.
         */
        
        d_flow_model->registerPatchWithDataContext(patch, intermediate_context[n]);
        
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > conservative_variables_intermediate =
            d_flow_model->getCellDataOfConservativeVariables();
        
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
                    
                    if (d_use_immersed_boundaries)
                    {
                        const int num_ghosts_0_IB_mask = num_ghosts_IB_mask[0];
                        
                        HAMERS_PRAGMA_SIMD
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute linear indices of conservative variable data and immersed boundary mask.
                            const int idx = i + num_ghosts_0_conservative_var;
                            const int idx_intermediate = i + num_ghosts_0_conservative_var_intermediate;
                            const int idx_IB_mask = i + num_ghosts_0_IB_mask;
                            
                            if (IB_mask[idx_IB_mask] == fluid)
                            {
                                Q[ei][idx] += alpha[n]*Q_intermediate[ei][idx_intermediate];
                            }
                            else
                            {
                                Q[ei][idx] = Q_intermediate[ei][idx_intermediate];
                            }
                        }
                    }
                    else
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute linear indices of conservative variable data.
                            const int idx = i + num_ghosts_0_conservative_var;
                            const int idx_intermediate = i + num_ghosts_0_conservative_var_intermediate;
                            
                            Q[ei][idx] += alpha[n]*Q_intermediate[ei][idx_intermediate];
                        }
                    }
                }
            }
            
            if (beta[n] != 0.0)
            {
                if (d_use_conservative_form_diffusive_flux)
                {
                    for (int ei = 0; ei < d_flow_model->getNumberOfEquations(); ei++)
                    {
                        double* F_c_x_intermediate = convective_flux_intermediate->getPointer(0, ei);
                        double* F_d_x_intermediate = diffusive_flux_intermediate->getPointer(0, ei);
                        double* S_intermediate = source_intermediate->getPointer(ei);
                        
                        const int num_ghosts_0_conservative_var = num_ghosts_conservative_var[ei][0];
                        
                        if (d_use_immersed_boundaries)
                        {
                            const int num_ghosts_0_IB_mask = num_ghosts_IB_mask[0];
                            
                            HAMERS_PRAGMA_SIMD
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                // Compute linear indices.
                                const int idx = i + num_ghosts_0_conservative_var;
                                const int idx_flux_x = i + 1;
                                const int idx_source = i;
                                const int idx_IB_mask = i + num_ghosts_0_IB_mask;
                                
                                if (IB_mask[idx_IB_mask] == fluid)
                                {
                                    Q[ei][idx] += beta[n]*
                                        (-(F_c_x_intermediate[idx_flux_x] - F_c_x_intermediate[idx_flux_x - 1] +
                                           F_d_x_intermediate[idx_flux_x] - F_d_x_intermediate[idx_flux_x - 1])/dx_0 +
                                         S_intermediate[idx_source]);
                                }
                            }
                        }
                        else
                        {
                            HAMERS_PRAGMA_SIMD
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                // Compute linear indices.
                                const int idx = i + num_ghosts_0_conservative_var;
                                const int idx_flux_x = i + 1;
                                const int idx_source = i;
                                
                                Q[ei][idx] += beta[n]*
                                    (-(F_c_x_intermediate[idx_flux_x] - F_c_x_intermediate[idx_flux_x - 1] +
                                       F_d_x_intermediate[idx_flux_x] - F_d_x_intermediate[idx_flux_x - 1])/dx_0 +
                                     S_intermediate[idx_source]);
                            }
                        }
                    }
                }
                else
                {
                    for (int ei = 0; ei < d_flow_model->getNumberOfEquations(); ei++)
                    {
                        double* F_c_x_intermediate = convective_flux_intermediate->getPointer(0, ei);
                        double* nabla_F_d_intermediate = diffusive_flux_divergence_intermediate->getPointer(ei);
                        double* S_intermediate = source_intermediate->getPointer(ei);
                        
                        const int num_ghosts_0_conservative_var = num_ghosts_conservative_var[ei][0];
                        
                        if (d_use_immersed_boundaries)
                        {
                            const int num_ghosts_0_IB_mask = num_ghosts_IB_mask[0];
                            
                            HAMERS_PRAGMA_SIMD
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                // Compute linear indices.
                                const int idx = i + num_ghosts_0_conservative_var;
                                const int idx_flux_x = i + 1;
                                const int idx_cell = i;
                                const int idx_IB_mask = i + num_ghosts_0_IB_mask;
                                
                                if (IB_mask[idx_IB_mask] == fluid)
                                {
                                    Q[ei][idx] += beta[n]*
                                        (-(F_c_x_intermediate[idx_flux_x] - F_c_x_intermediate[idx_flux_x - 1])/dx_0 -
                                         nabla_F_d_intermediate[idx_cell] +
                                         S_intermediate[idx_cell]);
                                }
                            }
                        }
                        else
                        {
                            HAMERS_PRAGMA_SIMD
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                // Compute linear indices.
                                const int idx = i + num_ghosts_0_conservative_var;
                                const int idx_flux_x = i + 1;
                                const int idx_cell = i;
                                
                                Q[ei][idx] += beta[n]*
                                    (-(F_c_x_intermediate[idx_flux_x] - F_c_x_intermediate[idx_flux_x - 1])/dx_0 -
                                     nabla_F_d_intermediate[idx_cell] +
                                     S_intermediate[idx_cell]);
                            }
                        }
                    }
                }
            }
            
            if (gamma[n] != 0.0)
            {
                // Accumulate the flux in the x direction.
                if (d_use_conservative_form_diffusive_flux)
                {
                    for (int ei = 0; ei < d_flow_model->getNumberOfEquations(); ei++)
                    {
                        double* F_c_x = convective_flux->getPointer(0, ei);
                        double* F_c_x_intermediate = convective_flux_intermediate->getPointer(0, ei);
                        
                        double* F_d_x = diffusive_flux->getPointer(0, ei);
                        double* F_d_x_intermediate = diffusive_flux_intermediate->getPointer(0, ei);
                        
                        HAMERS_PRAGMA_SIMD
                        for (int i = 0; i < interior_dim_0 + 1; i++)
                        {
                            // Compute linear index.
                            const int idx_flux_x = i;
                            
                            F_c_x[idx_flux_x] += gamma[n]*F_c_x_intermediate[idx_flux_x];
                            F_d_x[idx_flux_x] += gamma[n]*F_d_x_intermediate[idx_flux_x];
                        }
                    }
                }
                else
                {
                    for (int ei = 0; ei < d_flow_model->getNumberOfEquations(); ei++)
                    {
                        double* F_c_x = convective_flux->getPointer(0, ei);
                        double* F_c_x_intermediate = convective_flux_intermediate->getPointer(0, ei);
                        
                        HAMERS_PRAGMA_SIMD
                        for (int i = 0; i < interior_dim_0 + 1; i++)
                        {
                            // Compute linear index.
                            const int idx_flux_x = i;
                            
                            F_c_x[idx_flux_x] += gamma[n]*F_c_x_intermediate[idx_flux_x];
                        }
                    }
                }
                
                // Accumulate the source terms.
                if (d_use_conservative_form_diffusive_flux)
                {
                    for (int ei = 0; ei < d_flow_model->getNumberOfEquations(); ei++)
                    {
                        double* S = source->getPointer(ei);
                        double* S_intermediate = source_intermediate->getPointer(ei);
                        
                        HAMERS_PRAGMA_SIMD
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute linear index.
                            const int idx = i;
                            
                            S[idx] += gamma[n]*S_intermediate[idx];
                        }
                    }
                }
                else
                {
                    for (int ei = 0; ei < d_flow_model->getNumberOfEquations(); ei++)
                    {
                        double* nabla_F_d = diffusive_flux_divergence->getPointer(ei);
                        double* nabla_F_d_intermediate = diffusive_flux_divergence_intermediate->getPointer(ei);
                        
                        double* S = source->getPointer(ei);
                        double* S_intermediate = source_intermediate->getPointer(ei);
                        
                        HAMERS_PRAGMA_SIMD
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute linear index.
                            const int idx = i;
                            
                            nabla_F_d[idx] += gamma[n]*nabla_F_d_intermediate[idx];
                            S[idx] += gamma[n]*S_intermediate[idx];
                        }
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
                    
                    if (d_use_immersed_boundaries)
                    {
                        const int num_ghosts_0_IB_mask = num_ghosts_IB_mask[0];
                        const int num_ghosts_1_IB_mask = num_ghosts_IB_mask[1];
                        const int ghostcell_dim_0_IB_mask = ghostcell_dims_IB_mask[0];
                        
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            HAMERS_PRAGMA_SIMD
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                // Compute linear indices of conservative data and immersed boundary mask.
                                const int idx = (i + num_ghosts_0_conservative_var) +
                                    (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var;
                                
                                const int idx_intermediate = (i + num_ghosts_0_conservative_var_intermediate) +
                                    (j + num_ghosts_1_conservative_var_intermediate)*
                                        ghostcell_dim_0_conservative_var_intermediate;
                                
                                const int idx_IB_mask = (i + num_ghosts_0_IB_mask) +
                                    (j + num_ghosts_1_IB_mask)*ghostcell_dim_0_IB_mask;
                                
                                if (IB_mask[idx_IB_mask] == fluid)
                                {
                                    Q[ei][idx] += alpha[n]*Q_intermediate[ei][idx_intermediate];
                                }
                                else
                                {
                                    Q[ei][idx] = Q_intermediate[ei][idx_intermediate];
                                }
                            }
                        }
                    }
                    else
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            HAMERS_PRAGMA_SIMD
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
            }
            
            if (beta[n] != 0.0)
            {
                if (d_use_conservative_form_diffusive_flux)
                {
                    for (int ei = 0; ei < d_flow_model->getNumberOfEquations(); ei++)
                    {
                        double* F_c_x_intermediate = convective_flux_intermediate->getPointer(0, ei);
                        double* F_c_y_intermediate = convective_flux_intermediate->getPointer(1, ei);
                        double* F_d_x_intermediate = diffusive_flux_intermediate->getPointer(0, ei);
                        double* F_d_y_intermediate = diffusive_flux_intermediate->getPointer(1, ei);
                        double* S_intermediate = source_intermediate->getPointer(ei);
                        
                        const int num_ghosts_0_conservative_var = num_ghosts_conservative_var[ei][0];
                        const int num_ghosts_1_conservative_var = num_ghosts_conservative_var[ei][1];
                        const int ghostcell_dim_0_conservative_var = ghostcell_dims_conservative_var[ei][0];
                        
                        if (d_use_immersed_boundaries)
                        {
                            const int num_ghosts_0_IB_mask = num_ghosts_IB_mask[0];
                            const int num_ghosts_1_IB_mask = num_ghosts_IB_mask[1];
                            const int ghostcell_dim_0_IB_mask = ghostcell_dims_IB_mask[0];
                            
                            for (int j = 0; j < interior_dim_1; j++)
                            {
                                HAMERS_PRAGMA_SIMD
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
                                    
                                    const int idx_IB_mask = (i + num_ghosts_0_IB_mask) +
                                        (j + num_ghosts_1_IB_mask)*ghostcell_dim_0_IB_mask;
                                    
                                    if (IB_mask[idx_IB_mask] == fluid)
                                    {
                                        Q[ei][idx] += beta[n]*
                                            (-(F_c_x_intermediate[idx_flux_x_R] - F_c_x_intermediate[idx_flux_x_L] +
                                               F_d_x_intermediate[idx_flux_x_R] - F_d_x_intermediate[idx_flux_x_L])/dx_0 -
                                              (F_c_y_intermediate[idx_flux_y_T] - F_c_y_intermediate[idx_flux_y_B] +
                                               F_d_y_intermediate[idx_flux_y_T] - F_d_y_intermediate[idx_flux_y_B])/dx_1 +
                                              S_intermediate[idx_source]);
                                    }
                                }
                            }
                        }
                        else
                        {
                            for (int j = 0; j < interior_dim_1; j++)
                            {
                                HAMERS_PRAGMA_SIMD
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
                                        (-(F_c_x_intermediate[idx_flux_x_R] - F_c_x_intermediate[idx_flux_x_L] +
                                           F_d_x_intermediate[idx_flux_x_R] - F_d_x_intermediate[idx_flux_x_L])/dx_0 -
                                          (F_c_y_intermediate[idx_flux_y_T] - F_c_y_intermediate[idx_flux_y_B] +
                                           F_d_y_intermediate[idx_flux_y_T] - F_d_y_intermediate[idx_flux_y_B])/dx_1 +
                                          S_intermediate[idx_source]);
                                }
                            }
                        }
                    }
                }
                else
                {
                    for (int ei = 0; ei < d_flow_model->getNumberOfEquations(); ei++)
                    {
                        double* F_c_x_intermediate = convective_flux_intermediate->getPointer(0, ei);
                        double* F_c_y_intermediate = convective_flux_intermediate->getPointer(1, ei);
                        double* nabla_F_d_intermediate = diffusive_flux_divergence_intermediate->getPointer(ei);
                        double* S_intermediate = source_intermediate->getPointer(ei);
                        
                        const int num_ghosts_0_conservative_var = num_ghosts_conservative_var[ei][0];
                        const int num_ghosts_1_conservative_var = num_ghosts_conservative_var[ei][1];
                        const int ghostcell_dim_0_conservative_var = ghostcell_dims_conservative_var[ei][0];
                        
                        if (d_use_immersed_boundaries)
                        {
                            const int num_ghosts_0_IB_mask = num_ghosts_IB_mask[0];
                            const int num_ghosts_1_IB_mask = num_ghosts_IB_mask[1];
                            const int ghostcell_dim_0_IB_mask = ghostcell_dims_IB_mask[0];
                            
                            for (int j = 0; j < interior_dim_1; j++)
                            {
                                HAMERS_PRAGMA_SIMD
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
                                    
                                    const int idx_cell = i +
                                        j*interior_dim_0;
                                    
                                    const int idx_IB_mask = (i + num_ghosts_0_IB_mask) +
                                        (j + num_ghosts_1_IB_mask)*ghostcell_dim_0_IB_mask;
                                    
                                    if (IB_mask[idx_IB_mask] == fluid)
                                    {
                                        Q[ei][idx] += beta[n]*
                                            (-(F_c_x_intermediate[idx_flux_x_R] - F_c_x_intermediate[idx_flux_x_L])/dx_0 -
                                              (F_c_y_intermediate[idx_flux_y_T] - F_c_y_intermediate[idx_flux_y_B])/dx_1 -
                                              nabla_F_d_intermediate[idx_cell] +
                                              S_intermediate[idx_cell]);
                                    }
                                }
                            }
                        }
                        else
                        {
                            for (int j = 0; j < interior_dim_1; j++)
                            {
                                HAMERS_PRAGMA_SIMD
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
                                    
                                    const int idx_cell = i +
                                        j*interior_dim_0;
                                    
                                    Q[ei][idx] += beta[n]*
                                        (-(F_c_x_intermediate[idx_flux_x_R] - F_c_x_intermediate[idx_flux_x_L])/dx_0 -
                                          (F_c_y_intermediate[idx_flux_y_T] - F_c_y_intermediate[idx_flux_y_B])/dx_1 -
                                          nabla_F_d_intermediate[idx_cell] +
                                          S_intermediate[idx_cell]);
                                }
                            }
                        }
                    }
                }
            }
            
            if (gamma[n] != 0.0)
            {
                // Accumulate the flux in the x direction.
                if (d_use_conservative_form_diffusive_flux)
                {
                    for (int ei = 0; ei < d_flow_model->getNumberOfEquations(); ei++)
                    {
                        double* F_c_x = convective_flux->getPointer(0, ei);
                        double* F_c_x_intermediate = convective_flux_intermediate->getPointer(0, ei);
                        
                        double* F_d_x = diffusive_flux->getPointer(0, ei);
                        double* F_d_x_intermediate = diffusive_flux_intermediate->getPointer(0, ei);
                        
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            HAMERS_PRAGMA_SIMD
                            for (int i = 0; i < interior_dim_0 + 1; i++)
                            {
                                // Compute linear index.
                                const int idx_flux_x = i +
                                    j*(interior_dim_0 + 1);
                                
                                F_c_x[idx_flux_x] += gamma[n]*F_c_x_intermediate[idx_flux_x];
                                F_d_x[idx_flux_x] += gamma[n]*F_d_x_intermediate[idx_flux_x];
                            }                        
                        }
                    }
                    
                    // Accumulate the flux in the y direction.
                    for (int ei = 0; ei < d_flow_model->getNumberOfEquations(); ei++)
                    {
                        double* F_c_y = convective_flux->getPointer(1, ei);
                        double* F_c_y_intermediate = convective_flux_intermediate->getPointer(1, ei);
                        
                        double* F_d_y = diffusive_flux->getPointer(1, ei);
                        double* F_d_y_intermediate = diffusive_flux_intermediate->getPointer(1, ei);
                        
                        for (int j = 0; j < interior_dim_1 + 1; j++)
                        {
                            HAMERS_PRAGMA_SIMD
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                // Compute linear index.
                                const int idx_flux_y = i +
                                    j*interior_dim_0;
                                
                                F_c_y[idx_flux_y] += gamma[n]*F_c_y_intermediate[idx_flux_y];
                                F_d_y[idx_flux_y] += gamma[n]*F_d_y_intermediate[idx_flux_y];
                            }
                        }
                    }
                }
                else
                {
                    for (int ei = 0; ei < d_flow_model->getNumberOfEquations(); ei++)
                    {
                        double* F_c_x = convective_flux->getPointer(0, ei);
                        double* F_c_x_intermediate = convective_flux_intermediate->getPointer(0, ei);
                        
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            HAMERS_PRAGMA_SIMD
                            for (int i = 0; i < interior_dim_0 + 1; i++)
                            {
                                // Compute linear index.
                                const int idx_flux_x = i +
                                    j*(interior_dim_0 + 1);
                                
                                F_c_x[idx_flux_x] += gamma[n]*F_c_x_intermediate[idx_flux_x];
                            }                        
                        }
                    }
                    
                    // Accumulate the flux in the y direction.
                    for (int ei = 0; ei < d_flow_model->getNumberOfEquations(); ei++)
                    {
                        double* F_c_y = convective_flux->getPointer(1, ei);
                        double* F_c_y_intermediate = convective_flux_intermediate->getPointer(1, ei);
                        
                        for (int j = 0; j < interior_dim_1 + 1; j++)
                        {
                            HAMERS_PRAGMA_SIMD
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                // Compute linear index.
                                const int idx_flux_y = i +
                                    j*interior_dim_0;
                                
                                F_c_y[idx_flux_y] += gamma[n]*F_c_y_intermediate[idx_flux_y];
                            }
                        }
                    }
                }
                
                // Accumulate the source.
                if (d_use_conservative_form_diffusive_flux)
                {
                    for (int ei = 0; ei < d_flow_model->getNumberOfEquations(); ei++)
                    {
                        double* S = source->getPointer(ei);
                        double* S_intermediate = source_intermediate->getPointer(ei);
                        
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            HAMERS_PRAGMA_SIMD
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                // Compute linear index.
                                const int idx = i +
                                    j*interior_dim_0;
                                
                                S[idx] += gamma[n]*S_intermediate[idx];
                            }
                        }
                    }
                }
                else
                {
                    for (int ei = 0; ei < d_flow_model->getNumberOfEquations(); ei++)
                    {
                        double* nabla_F_d = diffusive_flux_divergence->getPointer(ei);
                        double* nabla_F_d_intermediate = diffusive_flux_divergence_intermediate->getPointer(ei);
                        
                        double* S = source->getPointer(ei);
                        double* S_intermediate = source_intermediate->getPointer(ei);
                        
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            HAMERS_PRAGMA_SIMD
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                // Compute linear index.
                                const int idx = i +
                                    j*interior_dim_0;
                                
                                nabla_F_d[idx] += gamma[n]*nabla_F_d_intermediate[idx];
                                S[idx] += gamma[n]*S_intermediate[idx];
                            }
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
                    
                    if (d_use_immersed_boundaries)
                    {
                        const int num_ghosts_0_IB_mask = num_ghosts_IB_mask[0];
                        const int num_ghosts_1_IB_mask = num_ghosts_IB_mask[1];
                        const int num_ghosts_2_IB_mask = num_ghosts_IB_mask[2];
                        const int ghostcell_dim_0_IB_mask = ghostcell_dims_IB_mask[0];
                        const int ghostcell_dim_1_IB_mask = ghostcell_dims_IB_mask[1];
                        
                        for (int k = 0; k < interior_dim_2; k++)
                        {
                            for (int j = 0; j < interior_dim_1; j++)
                            {
                                HAMERS_PRAGMA_SIMD
                                for (int i = 0; i < interior_dim_0; i++)
                                {
                                    // Compute linear indices of conservative variable data and immersed boundary mask.
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
                                    
                                    const int idx_IB_mask = (i + num_ghosts_0_IB_mask) +
                                        (j + num_ghosts_1_IB_mask)*ghostcell_dim_0_IB_mask +
                                        (k + num_ghosts_2_IB_mask)*ghostcell_dim_0_IB_mask*
                                            ghostcell_dim_1_IB_mask;
                                    
                                    if (IB_mask[idx_IB_mask] == fluid)
                                    {
                                        Q[ei][idx] += alpha[n]*Q_intermediate[ei][idx_intermediate];
                                    }
                                    else
                                    {
                                        Q[ei][idx] = Q_intermediate[ei][idx_intermediate];
                                    }
                                }
                            }
                        }
                    }
                    else
                    {
                        for (int k = 0; k < interior_dim_2; k++)
                        {
                            for (int j = 0; j < interior_dim_1; j++)
                            {
                                HAMERS_PRAGMA_SIMD
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
            }
            
            if (beta[n] != 0.0)
            {
                if (d_use_conservative_form_diffusive_flux)
                {
                    for (int ei = 0; ei < d_flow_model->getNumberOfEquations(); ei++)
                    {
                        double* F_c_x_intermediate = convective_flux_intermediate->getPointer(0, ei);
                        double* F_c_y_intermediate = convective_flux_intermediate->getPointer(1, ei);
                        double* F_c_z_intermediate = convective_flux_intermediate->getPointer(2, ei);
                        double* F_d_x_intermediate = diffusive_flux_intermediate->getPointer(0, ei);
                        double* F_d_y_intermediate = diffusive_flux_intermediate->getPointer(1, ei);
                        double* F_d_z_intermediate = diffusive_flux_intermediate->getPointer(2, ei);
                        double* S_intermediate = source_intermediate->getPointer(ei);
                        
                        const int num_ghosts_0_conservative_var = num_ghosts_conservative_var[ei][0];
                        const int num_ghosts_1_conservative_var = num_ghosts_conservative_var[ei][1];
                        const int num_ghosts_2_conservative_var = num_ghosts_conservative_var[ei][2];
                        const int ghostcell_dim_0_conservative_var = ghostcell_dims_conservative_var[ei][0];
                        const int ghostcell_dim_1_conservative_var = ghostcell_dims_conservative_var[ei][1];
                        
                        if (d_use_immersed_boundaries)
                        {
                            const int num_ghosts_0_IB_mask = num_ghosts_IB_mask[0];
                            const int num_ghosts_1_IB_mask = num_ghosts_IB_mask[1];
                            const int num_ghosts_2_IB_mask = num_ghosts_IB_mask[2];
                            const int ghostcell_dim_0_IB_mask = ghostcell_dims_IB_mask[0];
                            const int ghostcell_dim_1_IB_mask = ghostcell_dims_IB_mask[1];
                            
                            for (int k = 0; k < interior_dim_2; k++)
                            {
                                for (int j = 0; j < interior_dim_1; j++)
                                {
                                    HAMERS_PRAGMA_SIMD
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
                                        
                                        const int idx_IB_mask = (i + num_ghosts_0_IB_mask) +
                                            (j + num_ghosts_1_IB_mask)*ghostcell_dim_0_IB_mask +
                                            (k + num_ghosts_2_IB_mask)*ghostcell_dim_0_IB_mask*
                                                ghostcell_dim_1_IB_mask;
                                        
                                        if (IB_mask[idx_IB_mask] == fluid)
                                        {
                                            Q[ei][idx] += beta[n]*
                                                (-(F_c_x_intermediate[idx_flux_x_R] - F_c_x_intermediate[idx_flux_x_L] +
                                                   F_d_x_intermediate[idx_flux_x_R] - F_d_x_intermediate[idx_flux_x_L])/dx_0 -
                                                  (F_c_y_intermediate[idx_flux_y_T] - F_c_y_intermediate[idx_flux_y_B] +
                                                   F_d_y_intermediate[idx_flux_y_T] - F_d_y_intermediate[idx_flux_y_B])/dx_1 -
                                                  (F_c_z_intermediate[idx_flux_z_F] - F_c_z_intermediate[idx_flux_z_B] +
                                                   F_d_z_intermediate[idx_flux_z_F] - F_d_z_intermediate[idx_flux_z_B])/dx_2 +
                                                  S_intermediate[idx_source]);
                                        }
                                    }
                                }
                            }
                        }
                        else
                        {
                            for (int k = 0; k < interior_dim_2; k++)
                            {
                                for (int j = 0; j < interior_dim_1; j++)
                                {
                                    HAMERS_PRAGMA_SIMD
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
                                            (-(F_c_x_intermediate[idx_flux_x_R] - F_c_x_intermediate[idx_flux_x_L] +
                                               F_d_x_intermediate[idx_flux_x_R] - F_d_x_intermediate[idx_flux_x_L])/dx_0 -
                                              (F_c_y_intermediate[idx_flux_y_T] - F_c_y_intermediate[idx_flux_y_B] +
                                               F_d_y_intermediate[idx_flux_y_T] - F_d_y_intermediate[idx_flux_y_B])/dx_1 -
                                              (F_c_z_intermediate[idx_flux_z_F] - F_c_z_intermediate[idx_flux_z_B] +
                                               F_d_z_intermediate[idx_flux_z_F] - F_d_z_intermediate[idx_flux_z_B])/dx_2 +
                                              S_intermediate[idx_source]);
                                    }
                                }
                            }
                        }
                    }
                }
                else
                {
                    for (int ei = 0; ei < d_flow_model->getNumberOfEquations(); ei++)
                    {
                        double* F_c_x_intermediate = convective_flux_intermediate->getPointer(0, ei);
                        double* F_c_y_intermediate = convective_flux_intermediate->getPointer(1, ei);
                        double* F_c_z_intermediate = convective_flux_intermediate->getPointer(2, ei);
                        double* nabla_F_d_intermediate = diffusive_flux_divergence_intermediate->getPointer(ei);
                        double* S_intermediate = source_intermediate->getPointer(ei);
                        
                        const int num_ghosts_0_conservative_var = num_ghosts_conservative_var[ei][0];
                        const int num_ghosts_1_conservative_var = num_ghosts_conservative_var[ei][1];
                        const int num_ghosts_2_conservative_var = num_ghosts_conservative_var[ei][2];
                        const int ghostcell_dim_0_conservative_var = ghostcell_dims_conservative_var[ei][0];
                        const int ghostcell_dim_1_conservative_var = ghostcell_dims_conservative_var[ei][1];
                        
                        if (d_use_immersed_boundaries)
                        {
                            const int num_ghosts_0_IB_mask = num_ghosts_IB_mask[0];
                            const int num_ghosts_1_IB_mask = num_ghosts_IB_mask[1];
                            const int num_ghosts_2_IB_mask = num_ghosts_IB_mask[2];
                            const int ghostcell_dim_0_IB_mask = ghostcell_dims_IB_mask[0];
                            const int ghostcell_dim_1_IB_mask = ghostcell_dims_IB_mask[1];
                            
                            for (int k = 0; k < interior_dim_2; k++)
                            {
                                for (int j = 0; j < interior_dim_1; j++)
                                {
                                    HAMERS_PRAGMA_SIMD
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
                                        
                                        const int idx_cell = i +
                                            j*interior_dim_0 +
                                            k*interior_dim_0*interior_dim_1;
                                        
                                        const int idx_IB_mask = (i + num_ghosts_0_IB_mask) +
                                            (j + num_ghosts_1_IB_mask)*ghostcell_dim_0_IB_mask +
                                            (k + num_ghosts_2_IB_mask)*ghostcell_dim_0_IB_mask*
                                                ghostcell_dim_1_IB_mask;
                                        
                                        if (IB_mask[idx_IB_mask] == fluid)
                                        {
                                            Q[ei][idx] += beta[n]*
                                                (-(F_c_x_intermediate[idx_flux_x_R] - F_c_x_intermediate[idx_flux_x_L])/dx_0 -
                                                  (F_c_y_intermediate[idx_flux_y_T] - F_c_y_intermediate[idx_flux_y_B])/dx_1 -
                                                  (F_c_z_intermediate[idx_flux_z_F] - F_c_z_intermediate[idx_flux_z_B])/dx_2 -
                                                  nabla_F_d_intermediate[idx_cell] +
                                                  S_intermediate[idx_cell]);
                                        }
                                    }
                                }
                            }
                        }
                        else
                        {
                            for (int k = 0; k < interior_dim_2; k++)
                            {
                                for (int j = 0; j < interior_dim_1; j++)
                                {
                                    HAMERS_PRAGMA_SIMD
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
                                        
                                        const int idx_cell = i +
                                            j*interior_dim_0 +
                                            k*interior_dim_0*interior_dim_1;
                                        
                                        Q[ei][idx] += beta[n]*
                                            (-(F_c_x_intermediate[idx_flux_x_R] - F_c_x_intermediate[idx_flux_x_L])/dx_0 -
                                              (F_c_y_intermediate[idx_flux_y_T] - F_c_y_intermediate[idx_flux_y_B])/dx_1 -
                                              (F_c_z_intermediate[idx_flux_z_F] - F_c_z_intermediate[idx_flux_z_B])/dx_2 -
                                              nabla_F_d_intermediate[idx_cell] +
                                              S_intermediate[idx_cell]);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            
            if (gamma[n] != 0.0)
            {
                // Accumulate the flux in the x direction.
                if (d_use_conservative_form_diffusive_flux)
                {
                    for (int ei = 0; ei < d_flow_model->getNumberOfEquations(); ei++)
                    {
                        double* F_c_x = convective_flux->getPointer(0, ei);
                        double* F_c_x_intermediate = convective_flux_intermediate->getPointer(0, ei);
                        
                        double* F_d_x = diffusive_flux->getPointer(0, ei);
                        double* F_d_x_intermediate = diffusive_flux_intermediate->getPointer(0, ei);
                        
                        for (int k = 0; k < interior_dim_2; k++)
                        {
                            for (int j = 0; j < interior_dim_1; j++)
                            {
                                HAMERS_PRAGMA_SIMD
                                for (int i = 0; i < interior_dim_0 + 1; i++)
                                {
                                    // Compute linear index.
                                    const int idx_flux_x = i +
                                        j*(interior_dim_0 + 1) +
                                        k*(interior_dim_0 + 1)*interior_dim_1;
                                    
                                    F_c_x[idx_flux_x] += gamma[n]*F_c_x_intermediate[idx_flux_x];
                                    F_d_x[idx_flux_x] += gamma[n]*F_d_x_intermediate[idx_flux_x];
                                }                        
                            }
                        }
                    }
                    
                    // Accumulate the flux in the y direction.
                    for (int ei = 0; ei < d_flow_model->getNumberOfEquations(); ei++)
                    {
                        double* F_c_y = convective_flux->getPointer(1, ei);
                        double* F_c_y_intermediate = convective_flux_intermediate->getPointer(1, ei);
                        
                        double* F_d_y = diffusive_flux->getPointer(1, ei);
                        double* F_d_y_intermediate = diffusive_flux_intermediate->getPointer(1, ei);
                        
                        for (int k = 0; k < interior_dim_2; k++)
                        {
                            for (int j = 0; j < interior_dim_1 + 1; j++)
                            {
                                HAMERS_PRAGMA_SIMD
                                for (int i = 0; i < interior_dim_0; i++)
                                {
                                    // Compute linear index.
                                    const int idx_flux_y = i +
                                        j*interior_dim_0 +
                                        k*interior_dim_0*(interior_dim_1 + 1);
                                    
                                    F_c_y[idx_flux_y] += gamma[n]*F_c_y_intermediate[idx_flux_y];
                                    F_d_y[idx_flux_y] += gamma[n]*F_d_y_intermediate[idx_flux_y];
                                }
                            }
                        }
                    }
                    
                    // Accumulate the flux in the z direction.
                    for (int ei = 0; ei < d_flow_model->getNumberOfEquations(); ei++)
                    {
                        double* F_c_z = convective_flux->getPointer(2, ei);
                        double* F_c_z_intermediate = convective_flux_intermediate->getPointer(2, ei);
                        
                        double* F_d_z = diffusive_flux->getPointer(2, ei);
                        double* F_d_z_intermediate = diffusive_flux_intermediate->getPointer(2, ei);
                        
                        for (int k = 0; k < interior_dim_2 + 1; k++)
                        {
                            for (int j = 0; j < interior_dim_1; j++)
                            {
                                HAMERS_PRAGMA_SIMD
                                for (int i = 0; i < interior_dim_0; i++)
                                {
                                    // Compute linear index.
                                    const int idx_flux_z = i +
                                        j*interior_dim_0 +
                                        k*interior_dim_0*interior_dim_1;
                                    
                                    F_c_z[idx_flux_z] += gamma[n]*F_c_z_intermediate[idx_flux_z];
                                    F_d_z[idx_flux_z] += gamma[n]*F_d_z_intermediate[idx_flux_z];
                                }
                            }
                        }
                    }
                }
                else
                {
                    for (int ei = 0; ei < d_flow_model->getNumberOfEquations(); ei++)
                    {
                        double* F_c_x = convective_flux->getPointer(0, ei);
                        double* F_c_x_intermediate = convective_flux_intermediate->getPointer(0, ei);
                        
                        for (int k = 0; k < interior_dim_2; k++)
                        {
                            for (int j = 0; j < interior_dim_1; j++)
                            {
                                HAMERS_PRAGMA_SIMD
                                for (int i = 0; i < interior_dim_0 + 1; i++)
                                {
                                    // Compute linear index.
                                    const int idx_flux_x = i +
                                        j*(interior_dim_0 + 1) +
                                        k*(interior_dim_0 + 1)*interior_dim_1;
                                    
                                    F_c_x[idx_flux_x] += gamma[n]*F_c_x_intermediate[idx_flux_x];
                                }                        
                            }
                        }
                    }
                    
                    // Accumulate the flux in the y direction.
                    for (int ei = 0; ei < d_flow_model->getNumberOfEquations(); ei++)
                    {
                        double* F_c_y = convective_flux->getPointer(1, ei);
                        double* F_c_y_intermediate = convective_flux_intermediate->getPointer(1, ei);
                        
                        for (int k = 0; k < interior_dim_2; k++)
                        {
                            for (int j = 0; j < interior_dim_1 + 1; j++)
                            {
                                HAMERS_PRAGMA_SIMD
                                for (int i = 0; i < interior_dim_0; i++)
                                {
                                    // Compute linear index.
                                    const int idx_flux_y = i +
                                        j*interior_dim_0 +
                                        k*interior_dim_0*(interior_dim_1 + 1);
                                    
                                    F_c_y[idx_flux_y] += gamma[n]*F_c_y_intermediate[idx_flux_y];
                                }
                            }
                        }
                    }
                    
                    // Accumulate the flux in the z direction.
                    for (int ei = 0; ei < d_flow_model->getNumberOfEquations(); ei++)
                    {
                        double* F_c_z = convective_flux->getPointer(2, ei);
                        double* F_c_z_intermediate = convective_flux_intermediate->getPointer(2, ei);
                        
                        for (int k = 0; k < interior_dim_2 + 1; k++)
                        {
                            for (int j = 0; j < interior_dim_1; j++)
                            {
                                HAMERS_PRAGMA_SIMD
                                for (int i = 0; i < interior_dim_0; i++)
                                {
                                    // Compute linear index.
                                    const int idx_flux_z = i +
                                        j*interior_dim_0 +
                                        k*interior_dim_0*interior_dim_1;
                                    
                                    F_c_z[idx_flux_z] += gamma[n]*F_c_z_intermediate[idx_flux_z];
                                }
                            }
                        }
                    }
                }
                
                // Accumulate the source.
                if (d_use_conservative_form_diffusive_flux)
                {
                    for (int ei = 0; ei < d_flow_model->getNumberOfEquations(); ei++)
                    {
                        double* S = source->getPointer(ei);
                        double* S_intermediate = source_intermediate->getPointer(ei);
                        
                        for (int k = 0; k < interior_dim_2; k++)
                        {
                            for (int j = 0; j < interior_dim_1; j++)
                            {
                                HAMERS_PRAGMA_SIMD
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
                }
                else
                {
                    for (int ei = 0; ei < d_flow_model->getNumberOfEquations(); ei++)
                    {
                        double* nabla_F_d = diffusive_flux_divergence->getPointer(ei);
                        double* nabla_F_d_intermediate = diffusive_flux_divergence_intermediate->getPointer(ei);
                        
                        double* S = source->getPointer(ei);
                        double* S_intermediate = source_intermediate->getPointer(ei);
                        
                        for (int k = 0; k < interior_dim_2; k++)
                        {
                            for (int j = 0; j < interior_dim_1; j++)
                            {
                                HAMERS_PRAGMA_SIMD
                                for (int i = 0; i < interior_dim_0; i++)
                                {
                                    // Compute linear index.
                                    const int idx = i +
                                        j*interior_dim_0 +
                                        k*interior_dim_0*interior_dim_1;
                                    
                                    nabla_F_d[idx] += gamma[n]*nabla_F_d_intermediate[idx];
                                    S[idx] += gamma[n]*S_intermediate[idx];
                                }
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
            
            if (d_use_immersed_boundaries)
            {
                d_flow_model->updateCellDataOfConservativeVariables(
                    IB_mask_cell_data,
                    fluid);
            }
            else
            {
                d_flow_model->updateCellDataOfConservativeVariables();
            }
            
            d_flow_model->unregisterPatch();
        }
    }
    
    t_advance_step->stop();
}


void
NavierStokes::synchronizeFluxes(
    hier::Patch& patch,
    const double time,
    const double dt)
{
    NULL_USE(time);
    NULL_USE(dt);
    
    t_synchronize_fluxes->start();
    
    const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
        HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
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
    
    std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > conservative_variables =
        d_flow_model->getCellDataOfConservativeVariables();
    
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
    
    /*
     * Get the pointer to the cell data of the immersed boundary mask.
     * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
     */
    
    int* IB_mask = nullptr;
    HAMERS_SHARED_PTR<pdat::CellData<int> > IB_mask_cell_data;
    hier::IntVector num_ghosts_IB_mask(d_dim);
    hier::IntVector ghostcell_dims_IB_mask(d_dim);
    const int fluid = int(IB_MASK::FLUID);
    
    if (d_use_immersed_boundaries)
    {
        d_flow_model->setupImmersedBoundaryMethod();
        
        HAMERS_SHARED_PTR<FlowModelImmersedBoundaryMethod> flow_model_immersed_boundary_method =
            d_flow_model->getFlowModelImmersedBoundaryMethod();
        
        IB_mask_cell_data = flow_model_immersed_boundary_method->
            getCellDataOfImmersedBoundaryMask(getDataContext());
        
        IB_mask                = IB_mask_cell_data->getPointer(0);
        num_ghosts_IB_mask     = IB_mask_cell_data->getGhostCellWidth();
        ghostcell_dims_IB_mask = IB_mask_cell_data->getGhostBox().numberCells();
    }
    
    // Unregister the patch.
    d_flow_model->unregisterPatch();
    
    HAMERS_SHARED_PTR<pdat::SideData<double> > convective_flux(
        HAMERS_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
            patch.getPatchData(d_variable_convective_flux, getDataContext())));
    
    HAMERS_SHARED_PTR<pdat::SideData<double> > diffusive_flux;
    HAMERS_SHARED_PTR<pdat::CellData<double> > diffusive_flux_divergence;
    
    if (d_use_conservative_form_diffusive_flux)
    {
        diffusive_flux =
            HAMERS_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                patch.getPatchData(d_variable_diffusive_flux, getDataContext()));
    }
    else
    {
        diffusive_flux_divergence =
            HAMERS_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                patch.getPatchData(d_variable_diffusive_flux_divergence, getDataContext()));
    }
    
    HAMERS_SHARED_PTR<pdat::CellData<double> > source(
        HAMERS_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
            patch.getPatchData(d_variable_source, getDataContext())));
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(convective_flux);
    if (d_use_conservative_form_diffusive_flux)
    {
        TBOX_ASSERT(diffusive_flux);
    }
    else
    {
        TBOX_ASSERT(diffusive_flux_divergence);
    }
    TBOX_ASSERT(source);
    
    TBOX_ASSERT(convective_flux->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
    if (d_use_conservative_form_diffusive_flux)
    {
        TBOX_ASSERT(diffusive_flux->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
    }
    else
    {
        TBOX_ASSERT(diffusive_flux_divergence->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
    }
    TBOX_ASSERT(source->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
#endif
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the dimension and grid spacing.
         */
        
        const int interior_dim_0 = interior_dims[0];
        
        const double dx_0 = dx[0];
        
        if (d_use_conservative_form_diffusive_flux)
        {
            for (int ei = 0; ei < d_flow_model->getNumberOfEquations(); ei++)
            {
                double* F_c_x = convective_flux->getPointer(0, ei);
                double* F_d_x = diffusive_flux->getPointer(0, ei);
                double* S = source->getPointer(ei);
                
                const int num_ghosts_0_conservative_var = num_ghosts_conservative_var[ei][0];
                
                if (d_use_immersed_boundaries)
                {
                    const int num_ghosts_0_IB_mask = num_ghosts_IB_mask[0];
                    
                    HAMERS_PRAGMA_SIMD
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute linear indices.
                        const int idx = i + num_ghosts_0_conservative_var;
                        const int idx_flux_x_L = i;
                        const int idx_flux_x_R = i + 1;
                        const int idx_source = i;
                        const int idx_IB_mask = i + num_ghosts_0_IB_mask;
                        
                        if (IB_mask[idx_IB_mask] == fluid)
                        {
                            Q[ei][idx] += (-(F_c_x[idx_flux_x_R] - F_c_x[idx_flux_x_L] +
                                             F_d_x[idx_flux_x_R] - F_d_x[idx_flux_x_L])/dx_0 +
                                            S[idx_source]);
                        }
                    }
                }
                else
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute linear indices.
                        const int idx = i + num_ghosts_0_conservative_var;
                        const int idx_flux_x_L = i;
                        const int idx_flux_x_R = i + 1;
                        const int idx_source = i;
                        
                        Q[ei][idx] += (-(F_c_x[idx_flux_x_R] - F_c_x[idx_flux_x_L] +
                                         F_d_x[idx_flux_x_R] - F_d_x[idx_flux_x_L])/dx_0 +
                                        S[idx_source]);
                    }
                }
            }
        }
        else
        {
            for (int ei = 0; ei < d_flow_model->getNumberOfEquations(); ei++)
            {
                double* F_c_x = convective_flux->getPointer(0, ei);
                double* nabla_F_d = diffusive_flux_divergence->getPointer(ei);
                double* S = source->getPointer(ei);
                
                const int num_ghosts_0_conservative_var = num_ghosts_conservative_var[ei][0];
                
                if (d_use_immersed_boundaries)
                {
                    const int num_ghosts_0_IB_mask = num_ghosts_IB_mask[0];
                    
                    HAMERS_PRAGMA_SIMD
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute linear indices.
                        const int idx = i + num_ghosts_0_conservative_var;
                        const int idx_flux_x_L = i;
                        const int idx_flux_x_R = i + 1;
                        const int idx_cell = i;
                        const int idx_IB_mask = i + num_ghosts_0_IB_mask;
                        
                        if (IB_mask[idx_IB_mask] == fluid)
                        {
                            Q[ei][idx] += (-(F_c_x[idx_flux_x_R] - F_c_x[idx_flux_x_L])/dx_0 -
                                            nabla_F_d[idx_cell] +
                                            S[idx_cell]);
                        }
                    }
                }
                else
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute linear indices.
                        const int idx = i + num_ghosts_0_conservative_var;
                        const int idx_flux_x_L = i;
                        const int idx_flux_x_R = i + 1;
                        const int idx_cell = i;
                        
                        Q[ei][idx] += (-(F_c_x[idx_flux_x_R] - F_c_x[idx_flux_x_L])/dx_0 -
                                        nabla_F_d[idx_cell] +
                                        S[idx_cell]);
                    }
                }
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
        
        if (d_use_conservative_form_diffusive_flux)
        {
            for (int ei = 0; ei < d_flow_model->getNumberOfEquations(); ei++)
            {
                double* F_c_x = convective_flux->getPointer(0, ei);
                double* F_c_y = convective_flux->getPointer(1, ei);
                double* F_d_x = diffusive_flux->getPointer(0, ei);
                double* F_d_y = diffusive_flux->getPointer(1, ei);
                double* S = source->getPointer(ei);
                
                const int num_ghosts_0_conservative_var = num_ghosts_conservative_var[ei][0];
                const int num_ghosts_1_conservative_var = num_ghosts_conservative_var[ei][1];
                const int ghostcell_dim_0_conservative_var = ghostcell_dims_conservative_var[ei][0];
                
                if (d_use_immersed_boundaries)
                {
                    const int num_ghosts_0_IB_mask = num_ghosts_IB_mask[0];
                    const int num_ghosts_1_IB_mask = num_ghosts_IB_mask[1];
                    const int ghostcell_dim_0_IB_mask = ghostcell_dims_IB_mask[0];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
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
                            
                            const int idx_IB_mask = (i + num_ghosts_0_IB_mask) +
                                (j + num_ghosts_1_IB_mask)*ghostcell_dim_0_IB_mask;
                            
                            if (IB_mask[idx_IB_mask] == fluid)
                            {
                                Q[ei][idx] += (-(F_c_x[idx_flux_x_R] - F_c_x[idx_flux_x_L] +
                                                 F_d_x[idx_flux_x_R] - F_d_x[idx_flux_x_L])/dx_0 -
                                                (F_c_y[idx_flux_y_T] - F_c_y[idx_flux_y_B] +
                                                 F_d_y[idx_flux_y_T] - F_d_y[idx_flux_y_B])/dx_1 +
                                                S[idx_source]);
                            }
                        }
                    }
                }
                else
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
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
                            
                            Q[ei][idx] += (-(F_c_x[idx_flux_x_R] - F_c_x[idx_flux_x_L] +
                                             F_d_x[idx_flux_x_R] - F_d_x[idx_flux_x_L])/dx_0 -
                                            (F_c_y[idx_flux_y_T] - F_c_y[idx_flux_y_B] +
                                             F_d_y[idx_flux_y_T] - F_d_y[idx_flux_y_B])/dx_1 +
                                            S[idx_source]);
                        }
                    }
                }
            }
        }
        else
        {
            for (int ei = 0; ei < d_flow_model->getNumberOfEquations(); ei++)
            {
                double* F_c_x = convective_flux->getPointer(0, ei);
                double* F_c_y = convective_flux->getPointer(1, ei);
                double* nabla_F_d = diffusive_flux_divergence->getPointer(ei);
                double* S = source->getPointer(ei);
                
                const int num_ghosts_0_conservative_var = num_ghosts_conservative_var[ei][0];
                const int num_ghosts_1_conservative_var = num_ghosts_conservative_var[ei][1];
                const int ghostcell_dim_0_conservative_var = ghostcell_dims_conservative_var[ei][0];
                
                if (d_use_immersed_boundaries)
                {
                    const int num_ghosts_0_IB_mask = num_ghosts_IB_mask[0];
                    const int num_ghosts_1_IB_mask = num_ghosts_IB_mask[1];
                    const int ghostcell_dim_0_IB_mask = ghostcell_dims_IB_mask[0];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
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
                            
                            const int idx_cell = i +
                                j*interior_dim_0;
                            
                            const int idx_IB_mask = (i + num_ghosts_0_IB_mask) +
                                (j + num_ghosts_1_IB_mask)*ghostcell_dim_0_IB_mask;
                            
                            if (IB_mask[idx_IB_mask] == fluid)
                            {
                                Q[ei][idx] += (-(F_c_x[idx_flux_x_R] - F_c_x[idx_flux_x_L])/dx_0 -
                                                (F_c_y[idx_flux_y_T] - F_c_y[idx_flux_y_B])/dx_1 -
                                                nabla_F_d[idx_cell] +
                                                S[idx_cell]);
                            }   
                        }
                    }
                }
                else
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
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
                            
                            const int idx_cell = i +
                                j*interior_dim_0;
                            
                            Q[ei][idx] += (-(F_c_x[idx_flux_x_R] - F_c_x[idx_flux_x_L])/dx_0 -
                                            (F_c_y[idx_flux_y_T] - F_c_y[idx_flux_y_B])/dx_1 -
                                            nabla_F_d[idx_cell] +
                                            S[idx_cell]);
                        }
                    }
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
        
        if (d_use_conservative_form_diffusive_flux)
        {
            for (int ei = 0; ei < d_flow_model->getNumberOfEquations(); ei++)
            {
                double* F_c_x = convective_flux->getPointer(0, ei);
                double* F_c_y = convective_flux->getPointer(1, ei);
                double* F_c_z = convective_flux->getPointer(2, ei);
                double* F_d_x = diffusive_flux->getPointer(0, ei);
                double* F_d_y = diffusive_flux->getPointer(1, ei);
                double* F_d_z = diffusive_flux->getPointer(2, ei);
                double* S = source->getPointer(ei);
                
                const int num_ghosts_0_conservative_var = num_ghosts_conservative_var[ei][0];
                const int num_ghosts_1_conservative_var = num_ghosts_conservative_var[ei][1];
                const int num_ghosts_2_conservative_var = num_ghosts_conservative_var[ei][2];
                const int ghostcell_dim_0_conservative_var = ghostcell_dims_conservative_var[ei][0];
                const int ghostcell_dim_1_conservative_var = ghostcell_dims_conservative_var[ei][1];
                
                if (d_use_immersed_boundaries)
                {
                    const int num_ghosts_0_IB_mask = num_ghosts_IB_mask[0];
                    const int num_ghosts_1_IB_mask = num_ghosts_IB_mask[1];
                    const int num_ghosts_2_IB_mask = num_ghosts_IB_mask[2];
                    const int ghostcell_dim_0_IB_mask = ghostcell_dims_IB_mask[0];
                    const int ghostcell_dim_1_IB_mask = ghostcell_dims_IB_mask[1];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            HAMERS_PRAGMA_SIMD
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
                                
                                const int idx_IB_mask = (i + num_ghosts_0_IB_mask) +
                                    (j + num_ghosts_1_IB_mask)*ghostcell_dim_0_IB_mask +
                                    (k + num_ghosts_2_IB_mask)*ghostcell_dim_0_IB_mask*
                                        ghostcell_dim_1_IB_mask;
                                
                                if (IB_mask[idx_IB_mask] == fluid)
                                {
                                    Q[ei][idx] += (-(F_c_x[idx_flux_x_R] - F_c_x[idx_flux_x_L] +
                                                     F_d_x[idx_flux_x_R] - F_d_x[idx_flux_x_L])/dx_0 -
                                                    (F_c_y[idx_flux_y_T] - F_c_y[idx_flux_y_B] +
                                                     F_d_y[idx_flux_y_T] - F_d_y[idx_flux_y_B])/dx_1 -
                                                    (F_c_z[idx_flux_z_F] - F_c_z[idx_flux_z_B] +
                                                     F_d_z[idx_flux_z_F] - F_d_z[idx_flux_z_B])/dx_2 +
                                                    S[idx_source]);
                                }
                            }
                        }
                    }
                }
                else
                {
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            HAMERS_PRAGMA_SIMD
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
                                
                                Q[ei][idx] += (-(F_c_x[idx_flux_x_R] - F_c_x[idx_flux_x_L] +
                                                 F_d_x[idx_flux_x_R] - F_d_x[idx_flux_x_L])/dx_0 -
                                                (F_c_y[idx_flux_y_T] - F_c_y[idx_flux_y_B] +
                                                 F_d_y[idx_flux_y_T] - F_d_y[idx_flux_y_B])/dx_1 -
                                                (F_c_z[idx_flux_z_F] - F_c_z[idx_flux_z_B] +
                                                 F_d_z[idx_flux_z_F] - F_d_z[idx_flux_z_B])/dx_2 +
                                                S[idx_source]);
                            }
                        }
                    }
                }
            }
        }
        else
        {
            for (int ei = 0; ei < d_flow_model->getNumberOfEquations(); ei++)
            {
                double* F_c_x = convective_flux->getPointer(0, ei);
                double* F_c_y = convective_flux->getPointer(1, ei);
                double* F_c_z = convective_flux->getPointer(2, ei);
                double* nabla_F_d = diffusive_flux_divergence->getPointer(ei);
                double* S = source->getPointer(ei);
                
                const int num_ghosts_0_conservative_var = num_ghosts_conservative_var[ei][0];
                const int num_ghosts_1_conservative_var = num_ghosts_conservative_var[ei][1];
                const int num_ghosts_2_conservative_var = num_ghosts_conservative_var[ei][2];
                const int ghostcell_dim_0_conservative_var = ghostcell_dims_conservative_var[ei][0];
                const int ghostcell_dim_1_conservative_var = ghostcell_dims_conservative_var[ei][1];
                
                if (d_use_immersed_boundaries)
                {
                    const int num_ghosts_0_IB_mask = num_ghosts_IB_mask[0];
                    const int num_ghosts_1_IB_mask = num_ghosts_IB_mask[1];
                    const int num_ghosts_2_IB_mask = num_ghosts_IB_mask[2];
                    const int ghostcell_dim_0_IB_mask = ghostcell_dims_IB_mask[0];
                    const int ghostcell_dim_1_IB_mask = ghostcell_dims_IB_mask[1];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            HAMERS_PRAGMA_SIMD
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
                                
                                const int idx_cell = i +
                                    j*interior_dim_0 +
                                    k*interior_dim_0*interior_dim_1;
                                
                                const int idx_IB_mask = (i + num_ghosts_0_IB_mask) +
                                    (j + num_ghosts_1_IB_mask)*ghostcell_dim_0_IB_mask +
                                    (k + num_ghosts_2_IB_mask)*ghostcell_dim_0_IB_mask*
                                        ghostcell_dim_1_IB_mask;
                                
                                if (IB_mask[idx_IB_mask] == fluid)
                                {
                                    Q[ei][idx] += (-(F_c_x[idx_flux_x_R] - F_c_x[idx_flux_x_L])/dx_0 -
                                                    (F_c_y[idx_flux_y_T] - F_c_y[idx_flux_y_B])/dx_1 -
                                                    (F_c_z[idx_flux_z_F] - F_c_z[idx_flux_z_B])/dx_2 -
                                                    nabla_F_d[idx_cell] +
                                                    S[idx_cell]);
                                }
                            }
                        }
                    }
                }
                else
                {
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            HAMERS_PRAGMA_SIMD
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
                                
                                const int idx_cell = i +
                                    j*interior_dim_0 +
                                    k*interior_dim_0*interior_dim_1;
                                
                                Q[ei][idx] += (-(F_c_x[idx_flux_x_R] - F_c_x[idx_flux_x_L])/dx_0 -
                                                (F_c_y[idx_flux_y_T] - F_c_y[idx_flux_y_B])/dx_1 -
                                                (F_c_z[idx_flux_z_F] - F_c_z[idx_flux_z_B])/dx_2 -
                                                nabla_F_d[idx_cell] +
                                                S[idx_cell]);
                            }
                        }
                    }
                }
            }
        }
    }
    
    /*
     * Update the conservative variables.
     */
    
    d_flow_model->registerPatchWithDataContext(patch, getDataContext());
    
    if (d_use_immersed_boundaries)
    {
        d_flow_model->updateCellDataOfConservativeVariables(
            IB_mask_cell_data,
            fluid);
    }
    else
    {
        d_flow_model->updateCellDataOfConservativeVariables();
    }
    
    d_flow_model->unregisterPatch();
    
    t_synchronize_fluxes->stop();
}


/*
 * Tag cells for refinement using immersed boundary detector.
 */
void
NavierStokes::tagCellsOnPatchImmersedBdryDetector(
    hier::Patch& patch,
    const double regrid_time,
    const bool initial_error,
    const int tag_index,
    const bool uses_refine_regions_too,
    const bool uses_value_detector_too,
    const bool uses_gradient_detector_too,
    const bool uses_multiresolution_detector_too,
    const bool uses_integral_detector_too,
    const bool uses_richardson_extrapolation_too)
{
    NULL_USE(regrid_time);
    NULL_USE(initial_error);
    NULL_USE(uses_refine_regions_too);
    
    t_tagimmersedbdry->start();
    
    // Get the tags.
    HAMERS_SHARED_PTR<pdat::CellData<int> > tags(
        HAMERS_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
            patch.getPatchData(tag_index)));
    
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(tags);
    TBOX_ASSERT(tags->getGhostCellWidth() == 0);
#endif
    
    // Initialize values of all tags to zero.
    if ((!uses_richardson_extrapolation_too) &&
        (!uses_integral_detector_too) &&
        (!uses_multiresolution_detector_too) &&
        (!uses_gradient_detector_too) &&
        (!uses_value_detector_too))
    {
        tags->fillAll(0);
    }
    
    // Tag the cells by using d_immersed_boundary_tagger.
    if (d_immersed_boundary_tagger != nullptr)
    {
        d_flow_model->registerPatchWithDataContext(patch, getDataContext());
        
        /*
         * Compute the immersed boundary method variables.
         */
        
        d_flow_model->setupImmersedBoundaryMethod();
        
        HAMERS_SHARED_PTR<FlowModelImmersedBoundaryMethod> flow_model_immersed_boundary_method =
            d_flow_model->getFlowModelImmersedBoundaryMethod();
        
        const hier::Box empty_box = hier::Box::getEmptyBox(d_dim);
        
        flow_model_immersed_boundary_method->setImmersedBoundaryMethodVariables(
            empty_box,
            regrid_time,
            false,
            getDataContext());
        
        d_flow_model->unregisterPatch();
        
        /*
         * Tag the cells using immersed boundary method variables.
         */
        
        d_immersed_boundary_tagger->tagCellsOnPatch(
            patch,
            tags,
            getDataContext());
    }
    
    t_tagimmersedbdry->stop();
}


/*
 * Preprocess before tagging cells using value detector.
 */
void
NavierStokes::preprocessTagCellsValueDetector(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const int level_number,
    const double regrid_time,
    const bool initial_error,
    const bool uses_refine_regions_too,
    const bool uses_immersed_bdry_detector_too,
    const bool uses_gradient_detector_too,
    const bool uses_multiresolution_detector_too,
    const bool uses_integral_detector_too,
    const bool uses_richardson_extrapolation_too)
{
    NULL_USE(regrid_time);
    NULL_USE(initial_error);
    NULL_USE(uses_refine_regions_too);
    NULL_USE(uses_immersed_bdry_detector_too);
    NULL_USE(uses_gradient_detector_too);
    NULL_USE(uses_multiresolution_detector_too);
    NULL_USE(uses_integral_detector_too);
    NULL_USE(uses_richardson_extrapolation_too);
    
    if (d_value_tagger != nullptr)
    {
        HAMERS_SHARED_PTR<hier::PatchLevel> level(
            patch_hierarchy->getPatchLevel(level_number));
        
        for (hier::PatchLevel::iterator ip(level->begin());
             ip != level->end();
             ip++)
        {
            const HAMERS_SHARED_PTR<hier::Patch>& patch = *ip;
            
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
NavierStokes::tagCellsOnPatchValueDetector(
    hier::Patch& patch,
    const double regrid_time,
    const bool initial_error,
    const int tag_index,
    const bool uses_refine_regions_too,
    const bool uses_immersed_bdry_detector_too,
    const bool uses_gradient_detector_too,
    const bool uses_multiresolution_detector_too,
    const bool uses_integral_detector_too,
    const bool uses_richardson_extrapolation_too)
{
    NULL_USE(regrid_time);
    NULL_USE(initial_error);
    NULL_USE(uses_refine_regions_too);
    NULL_USE(uses_immersed_bdry_detector_too);
    
    t_tagvalue->start();
    
    // Get the tags.
    HAMERS_SHARED_PTR<pdat::CellData<int> > tags(
        HAMERS_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
            patch.getPatchData(tag_index)));
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
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
NavierStokes::preprocessTagCellsGradientDetector(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const int level_number,
    const double regrid_time,
    const bool initial_error,
    const bool uses_refine_regions_too,
    const bool uses_immersed_bdry_detector_too,
    const bool uses_value_detector_too,
    const bool uses_multiresolution_detector_too,
    const bool uses_integral_detector_too,
    const bool uses_richardson_extrapolation_too)
{
    NULL_USE(regrid_time);
    NULL_USE(initial_error);
    NULL_USE(uses_refine_regions_too);
    NULL_USE(uses_immersed_bdry_detector_too);
    NULL_USE(uses_value_detector_too);
    NULL_USE(uses_multiresolution_detector_too);
    NULL_USE(uses_integral_detector_too);
    NULL_USE(uses_richardson_extrapolation_too);
    
    if (d_gradient_tagger != nullptr)
    {
        HAMERS_SHARED_PTR<hier::PatchLevel> level(
            patch_hierarchy->getPatchLevel(level_number));
        
        for (hier::PatchLevel::iterator ip(level->begin());
             ip != level->end();
             ip++)
        {
            const HAMERS_SHARED_PTR<hier::Patch>& patch = *ip;
            
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
NavierStokes::tagCellsOnPatchGradientDetector(
    hier::Patch& patch,
    const double regrid_time,
    const bool initial_error,
    const int tag_index,
    const bool uses_refine_regions_too,
    const bool uses_immersed_bdry_detector_too,
    const bool uses_value_detector_too,
    const bool uses_multiresolution_detector_too,
    const bool uses_integral_detector_too,
    const bool uses_richardson_extrapolation_too)
{
    NULL_USE(regrid_time);
    NULL_USE(initial_error);
    NULL_USE(uses_refine_regions_too);
    NULL_USE(uses_immersed_bdry_detector_too);
    NULL_USE(uses_value_detector_too);
    
    t_taggradient->start();
    
    // Get the tags.
    HAMERS_SHARED_PTR<pdat::CellData<int> > tags(
        HAMERS_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
            patch.getPatchData(tag_index)));
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
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
NavierStokes::preprocessTagCellsMultiresolutionDetector(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const int level_number,
    const double regrid_time,
    const bool initial_error,
    const bool uses_refine_regions_too,
    const bool uses_immersed_bdry_detector_too,
    const bool uses_value_detector_too,
    const bool uses_gradient_detector_too,
    const bool uses_integral_detector_too,
    const bool uses_richardson_extrapolation_too)
{
    NULL_USE(regrid_time);
    NULL_USE(initial_error);
    NULL_USE(uses_refine_regions_too);
    NULL_USE(uses_immersed_bdry_detector_too);
    NULL_USE(uses_value_detector_too);
    NULL_USE(uses_gradient_detector_too);
    NULL_USE(uses_integral_detector_too);
    NULL_USE(uses_richardson_extrapolation_too);
    
    if (d_multiresolution_tagger != nullptr)
    {
        HAMERS_SHARED_PTR<hier::PatchLevel> level(
            patch_hierarchy->getPatchLevel(level_number));
        
        for (hier::PatchLevel::iterator ip(level->begin());
             ip != level->end();
             ip++)
        {
            const HAMERS_SHARED_PTR<hier::Patch>& patch = *ip;
            
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
NavierStokes::tagCellsOnPatchMultiresolutionDetector(
    hier::Patch& patch,
    const double regrid_time,
    const bool initial_error,
    const int tag_index,
    const bool uses_refine_regions_too,
    const bool uses_immersed_bdry_detector_too,
    const bool uses_value_detector_too,
    const bool uses_gradient_detector_too,
    const bool uses_integral_detector_too,
    const bool uses_richardson_extrapolation_too)
{
    NULL_USE(regrid_time);
    NULL_USE(initial_error);
    NULL_USE(uses_refine_regions_too);
    NULL_USE(uses_immersed_bdry_detector_too);
    NULL_USE(uses_value_detector_too);
    NULL_USE(uses_gradient_detector_too);
    
    t_tagmultiresolution->start();
    
    // Get the tags.
    HAMERS_SHARED_PTR<pdat::CellData<int> > tags(
        HAMERS_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
            patch.getPatchData(tag_index)));
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
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
NavierStokes::setPhysicalBoundaryConditions(
    hier::Patch& patch,
    const double fill_time,
    const hier::IntVector& ghost_width_to_fill)
{
    t_setphysbcs->start();
    
    d_Navier_Stokes_boundary_conditions->setPhysicalBoundaryConditions(
        patch,
        fill_time,
        ghost_width_to_fill,
        getDataContext());
    
    t_setphysbcs->stop();
}


void
NavierStokes::putToRestart(
    const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const
{
    TBOX_ASSERT(restart_db);
    
    restart_db->putString("d_project_name", d_project_name);
    
    restart_db->putInteger("d_num_species", d_num_species);
    
    restart_db->putString("d_flow_model_str", d_flow_model_str);
    
    // Put the properties of d_flow_model into the restart database.
    HAMERS_SHARED_PTR<tbox::Database> restart_flow_model_db =
        restart_db->putDatabase("d_flow_model_db");
    d_flow_model->putToRestart(restart_flow_model_db);
    
    restart_db->putString("d_convective_flux_reconstructor_str", d_convective_flux_reconstructor_str);
    
    // Put the properties of d_convective_flux_reconstructor into the restart database.
    HAMERS_SHARED_PTR<tbox::Database> restart_convective_flux_reconstructor_db =
        restart_db->putDatabase("d_convective_flux_reconstructor_db");
    d_convective_flux_reconstructor->putToRestart(restart_convective_flux_reconstructor_db);
    
    restart_db->putBool("d_use_conservative_form_diffusive_flux", d_use_conservative_form_diffusive_flux);
    
    if (d_use_conservative_form_diffusive_flux)
    {
        restart_db->putString("d_diffusive_flux_reconstructor_str", d_diffusive_flux_reconstructor_str);
        
        // Put the properties of d_diffusive_flux_reconstructor into the restart database.
        HAMERS_SHARED_PTR<tbox::Database> restart_diffusive_flux_reconstructor_db =
            restart_db->putDatabase("d_diffusive_flux_reconstructor_db");
        d_diffusive_flux_reconstructor->putToRestart(restart_diffusive_flux_reconstructor_db);
    }
    else
    {
        restart_db->putString("d_nonconservative_diffusive_flux_divergence_operator_str",
            d_nonconservative_diffusive_flux_divergence_operator_str);
        
        // Put the properties of d_nonconservative_diffusive_flux_divergence_operator into the restart database.
        HAMERS_SHARED_PTR<tbox::Database> restart_nonconservative_diffusive_flux_divergence_operator_db =
            restart_db->putDatabase("d_nonconservative_diffusive_flux_divergence_operator_db");
        d_nonconservative_diffusive_flux_divergence_operator->putToRestart(
            restart_nonconservative_diffusive_flux_divergence_operator_db);
    }
    
    HAMERS_SHARED_PTR<tbox::Database> restart_Navier_Stokes_boundary_conditions_db =
        restart_db->putDatabase("d_Navier_Stokes_boundary_conditions_db");
    
    d_Navier_Stokes_boundary_conditions->putToRestart(restart_Navier_Stokes_boundary_conditions_db);
    
    restart_db->putBool("d_use_ghost_cell_immersed_boundary_method", d_use_ghost_cell_immersed_boundary_method);
    if (d_use_ghost_cell_immersed_boundary_method)
    {
        HAMERS_SHARED_PTR<tbox::Database> immersed_boundary_method_db =
            restart_db->putDatabase("immersed_boundary_method_db");
        
        HAMERS_SHARED_PTR<FlowModelImmersedBoundaryMethod> flow_model_immersed_boundary_method =
            d_flow_model->getFlowModelImmersedBoundaryMethod();
        
        flow_model_immersed_boundary_method->putToRestart(immersed_boundary_method_db);
    }
    
    restart_db->putInteger("d_immersed_boundary_tagger_num_cells_buffer", d_immersed_boundary_tagger_num_cells_buffer);
    
    if (d_value_tagger != nullptr)
    {
        HAMERS_SHARED_PTR<tbox::Database> restart_value_tagger_db =
            restart_db->putDatabase("d_value_tagger_db");
        
        d_value_tagger->putToRestart(restart_value_tagger_db);
    }
    
    if (d_gradient_tagger != nullptr)
    {
        HAMERS_SHARED_PTR<tbox::Database> restart_gradient_tagger_db =
            restart_db->putDatabase("d_gradient_tagger_db");
        
        d_gradient_tagger->putToRestart(restart_gradient_tagger_db);
    }
    
    if (d_multiresolution_tagger != nullptr)
    {
        HAMERS_SHARED_PTR<tbox::Database> restart_multiresolution_tagger_db =
            restart_db->putDatabase("d_multiresolution_tagger_db");
        
        d_multiresolution_tagger->putToRestart(restart_multiresolution_tagger_db);
    }
}


#ifdef HAVE_HDF5
void
NavierStokes::registerVisItDataWriter(
    const HAMERS_SHARED_PTR<ExtendedVisItDataWriter>& viz_writer)
{
    TBOX_ASSERT(viz_writer);
    d_visit_writer = viz_writer;
}
#endif


bool
NavierStokes::packDerivedDataIntoDoubleBuffer(
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


void NavierStokes::printClassData(std::ostream& os) const
{
    os << "\nPrint NavierStokes object... " << std::endl;
    os << std::endl;
    
    os << "NavierStokes: this = " << (NavierStokes *)this << std::endl;
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
    
    if (d_use_conservative_form_diffusive_flux)
    {
        /*
         * Print data of d_diffusive_flux_reconstructor.
         */
        
        d_diffusive_flux_reconstructor->printClassData(os);
    }
    else
    {
        /*
         * Print data of d_nonconservative_diffusive_flux_divergence_operator.
         */
        
        d_nonconservative_diffusive_flux_divergence_operator->printClassData(os);
    }
    os << "--------------------------------------------------------------------------------";
    
    /*
     * Print Refinement data
     */
    
    /*
     * Print data of d_Navier_Stokes_boundary_conditions.
     */
    
    d_Navier_Stokes_boundary_conditions->printClassData(os);
}


void
NavierStokes::printErrorStatistics(
    std::ostream& os,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const double time) const
{
    d_Navier_Stokes_error_statistics->printErrorStatistics(os, patch_hierarchy, d_plot_context, time);
}


void
NavierStokes::computeAndOutputMonitoringDataStatistics(
    std::ostream& os,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const int step_num,
    const double time) const
{
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    math::HierarchyCellDataOpsReal<double> cell_double_operator(patch_hierarchy, 0, 0);
    
    hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();
    
    std::vector<std::string> variable_names = d_flow_model->getNamesOfConservativeVariables();
    
    std::vector<HAMERS_SHARED_PTR<pdat::CellVariable<double> > > variables =
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
    
    d_flow_model->setupMonitoringStatisticsUtilities();
    
    HAMERS_SHARED_PTR<FlowModelMonitoringStatisticsUtilities> monitoring_statistics_utilities =
        d_flow_model->getFlowModelMonitoringStatisticsUtilities();
    
    if (monitoring_statistics_utilities->hasMonitoringStatistics())
    {
        if (monitoring_statistics_utilities->isStepToOutputMonitoringStatistics(step_num))
        {
            monitoring_statistics_utilities->computeMonitoringStatistics(
                patch_hierarchy,
                d_plot_context,
                step_num,
                time);
            
            monitoring_statistics_utilities->outputMonitoringStatistics(
                os,
                d_monitoring_stat_dump_filename,
                step_num,
                time);
        }
    }
}


/**
 * Compute variables for computing the statistics of data.
 */
void
NavierStokes::computeStatisticsVariables(
   const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy)
{
    d_flow_model->setupStatisticsUtilities();
    
    HAMERS_SHARED_PTR<FlowModelStatisticsUtilities> flow_model_statistics_utilities =
        d_flow_model->getFlowModelStatisticsUtilities();
    
    flow_model_statistics_utilities->computeVariables(
        patch_hierarchy,
        getDataContext());
}


/**
 * Filter variables for computing the statistics of data.
 */
void
NavierStokes::filterStatisticsVariables(
    const int level,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy)
{
    d_flow_model->setupStatisticsUtilities();
    
     HAMERS_SHARED_PTR<FlowModelStatisticsUtilities> flow_model_statistics_utilities =
        d_flow_model->getFlowModelStatisticsUtilities();
    
    flow_model_statistics_utilities->filterVariables(
        level,
        patch_hierarchy,
        getDataContext());
}


/**
 * Output the header of monitoring statistics.
 */
void
NavierStokes::outputHeaderMonitoringStatistics()
{
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    d_flow_model->setupMonitoringStatisticsUtilities();
    
    HAMERS_SHARED_PTR<FlowModelMonitoringStatisticsUtilities> monitoring_statistics_utilities =
        d_flow_model->getFlowModelMonitoringStatisticsUtilities();
    
    if (monitoring_statistics_utilities->hasMonitoringStatistics())
    {
        if (mpi.getRank() == 0)
        {
            std::ofstream f_out;
            f_out.open(d_monitoring_stat_dump_filename.c_str(), std::ios::app);
            
            if (!f_out.is_open())
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "Failed to open file to output statistics!"
                    << std::endl);
            }
            
            f_out << "# TIME               ";
            f_out.close();
        }
        
        monitoring_statistics_utilities->outputMonitoringStatisticalQuantitiesNames(
            d_monitoring_stat_dump_filename);
        
        if (mpi.getRank() == 0)
        {
            std::ofstream f_out;
            f_out.open(d_monitoring_stat_dump_filename.c_str(), std::ios::app);
            
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


/**
 * Output the header of statistics.
 */
void
NavierStokes::outputHeaderStatistics()
{
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
            
            f_out << "# TIME               ";
            f_out.close();
        }
        
        HAMERS_SHARED_PTR<FlowModelStatisticsUtilities> flow_model_statistics_utilities =
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


/**
 * Compute the statistics of data.
 */
void
NavierStokes::computeDataStatistics(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const double statistics_data_time)
{
    d_flow_model->setupStatisticsUtilities();
    
    HAMERS_SHARED_PTR<FlowModelStatisticsUtilities> flow_model_statistics_utilities =
        d_flow_model->getFlowModelStatisticsUtilities();
    
    flow_model_statistics_utilities->computeStatisticalQuantities(
        patch_hierarchy,
        getDataContext(),
        statistics_data_time);
}


/**
 * Output the statistics of data.
 */
void
NavierStokes::outputDataStatistics(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
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
            
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << output_time;
            f_out.close();
        }
        
        d_flow_model->setupStatisticsUtilities();
        
        HAMERS_SHARED_PTR<FlowModelStatisticsUtilities> flow_model_statistics_utilities =
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


/**
 * Get object of storing ensemble statistics.
 */
HAMERS_SHARED_PTR<EnsembleStatistics>
NavierStokes::getEnsembleStatistics()
{
    d_flow_model->setupStatisticsUtilities();
    
    HAMERS_SHARED_PTR<FlowModelStatisticsUtilities> flow_model_statistics_utilities =
        d_flow_model->getFlowModelStatisticsUtilities();
    
    return flow_model_statistics_utilities->getEnsembleStatistics();
}


/**
 * Set object of storing ensemble statistics.
 */
void
NavierStokes::setEnsembleStatistics(const HAMERS_SHARED_PTR<EnsembleStatistics> ensemble_statistics)
{
    d_flow_model->setupStatisticsUtilities();
    
    HAMERS_SHARED_PTR<FlowModelStatisticsUtilities> flow_model_statistics_utilities =
        d_flow_model->getFlowModelStatisticsUtilities();
    
    flow_model_statistics_utilities->setEnsembleStatistics(ensemble_statistics);
}


void
NavierStokes::getFromInput(
    const HAMERS_SHARED_PTR<tbox::Database>& input_db,
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
    
    const bool read_on_restart = input_db->getBoolWithDefault("read_on_restart", false);
    
    if (!is_from_restart || read_on_restart)
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
        
        d_use_conservative_form_diffusive_flux = input_db->getBoolWithDefault(
            "use_conservative_form_diffusive_flux", true);
        
        if (d_use_conservative_form_diffusive_flux)
        {
            /*
             * Get the diffusive flux reconstructor.
             */
            if (input_db->keyExists("diffusive_flux_reconstructor"))
            {
                d_diffusive_flux_reconstructor_str = input_db->getString("diffusive_flux_reconstructor");
            }
            else
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "Key data 'diffusive_flux_reconstructor' not found in input database."
                    << std::endl);
            }
            
            /*
             * Get the database of the diffusive flux reconstructor.
             */
            if (input_db->keyExists("Diffusive_flux_reconstructor"))
            {
                d_diffusive_flux_reconstructor_db = input_db->getDatabase("Diffusive_flux_reconstructor");
            }
            else
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "Key data 'Diffusive_flux_reconstructor' not found in input database."
                    << std::endl);
            }
        }
        else
        {
            /*
             * Get the non-conservative diffusive flux divergence operator.
             */
            if (input_db->keyExists("nonconservative_diffusive_flux_divergence_operator"))
            {
                d_nonconservative_diffusive_flux_divergence_operator_str = input_db->getString(
                    "nonconservative_diffusive_flux_divergence_operator");
            }
            else
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "Key data 'nonconservative_diffusive_flux_divergence_operator' not found in input database."
                    << std::endl);
            }
            
            /*
             * Get the database of the non-conservative diffusive flux divergence operator.
             */
            if (input_db->keyExists("Nonconservative_diffusive_flux_divergence_operator"))
            {
                d_diffusive_flux_reconstructor_db = input_db->getDatabase(
                    "Nonconservative_diffusive_flux_divergence_operator");
            }
            else
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "Key data 'Nonconservative_diffusive_flux_divergence_operator' not found in input database."
                    << std::endl);
            }
        }
        
        if (input_db->keyExists("immersed_boundary_tagger_num_cells_buffer"))
        {
            d_immersed_boundary_tagger_num_cells_buffer =
                input_db->getInteger("immersed_boundary_tagger_num_cells_buffer");
        }
        else
        {
            d_immersed_boundary_tagger_num_cells_buffer = -1;
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
     * Get the initial conditions from the input database.
     */
    if (input_db->keyExists("Initial_conditions"))
    {
        d_Navier_Stokes_initial_conditions_db = input_db->getDatabase(
            "Initial_conditions");
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
            d_Navier_Stokes_boundary_conditions_db = input_db->getDatabase(
                "Boundary_data");
            
            d_Navier_Stokes_boundary_conditions_db_is_from_restart = false;
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
    
    /*
     * Get whether to use immersed boundaries and the database for the immersed boundary method from the input database.
     */
    
    if (input_db->keyExists("use_ghost_cell_immersed_boundary_method"))
    {
        d_use_ghost_cell_immersed_boundary_method = input_db->getBool("use_ghost_cell_immersed_boundary_method");
        if (d_use_ghost_cell_immersed_boundary_method)
        {
            d_immersed_boundary_method_db = input_db->getDatabase("Immersed_boundary_method");
        }
    }
}


void NavierStokes::getFromRestart()
{
    HAMERS_SHARED_PTR<tbox::Database> root_db(tbox::RestartManager::getManager()->getRootDatabase());
    
    if (!root_db->isDatabase(d_object_name))
    {
        TBOX_ERROR("Restart database corresponding to "
                   << d_object_name
                   << " not found in restart file."
                   << std::endl);
    }
    
    HAMERS_SHARED_PTR<tbox::Database> db(root_db->getDatabase(d_object_name));
    
    d_project_name = db->getString("d_project_name");
    
    d_num_species = db->getInteger("d_num_species");
    
    d_flow_model_str = db->getString("d_flow_model_str");
    
    d_flow_model_db = db->getDatabase("d_flow_model_db");
    
    d_convective_flux_reconstructor_str = db->getString("d_convective_flux_reconstructor_str");
    
    d_convective_flux_reconstructor_db = db->getDatabase("d_convective_flux_reconstructor_db");
    
    d_use_conservative_form_diffusive_flux = db->getBool("d_use_conservative_form_diffusive_flux");
    
    if (d_use_conservative_form_diffusive_flux)
    {
        d_diffusive_flux_reconstructor_str = db->getString("d_diffusive_flux_reconstructor_str");
        d_diffusive_flux_reconstructor_db = db->getDatabase("d_diffusive_flux_reconstructor_db");
    }
    else
    {
        d_nonconservative_diffusive_flux_divergence_operator_str =
            db->getString("d_nonconservative_diffusive_flux_divergence_operator_str");
        d_nonconservative_diffusive_flux_divergence_operator_db =
            db->getDatabase("d_nonconservative_diffusive_flux_divergence_operator_db");
    }
    
    d_Navier_Stokes_boundary_conditions_db = db->getDatabase("d_Navier_Stokes_boundary_conditions_db");
    
    d_Navier_Stokes_boundary_conditions_db_is_from_restart = true;
    
    d_use_ghost_cell_immersed_boundary_method = db->getBool("d_use_ghost_cell_immersed_boundary_method");
    if (d_use_ghost_cell_immersed_boundary_method)
    {
        d_immersed_boundary_method_db = db->getDatabase("d_immersed_boundary_method_db");
    }
    
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
