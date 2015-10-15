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

//integer constants for boundary conditions
#define CHECK_BDRY_DATA (0)
#include "SAMRAI/appu/CartesianBoundaryDefines.h"

//integer constant for debugging improperly set boundary data
#define BOGUS_BDRY_DATA (-9999)

// routines for managing boundary data
#include "SAMRAI/appu/CartesianBoundaryUtilities2.h"
#include "SAMRAI/appu/CartesianBoundaryUtilities3.h"

#include "flow_model/convective_flux_reconstructor/ConvectiveFluxReconstructorLLF.hpp"
#include "flow_model/convective_flux_reconstructor/ConvectiveFluxReconstructorFirstOrderHLLC.hpp"
#include "flow_model/convective_flux_reconstructor/ConvectiveFluxReconstructorWENO-JS5-LLF.hpp"
#include "flow_model/convective_flux_reconstructor/ConvectiveFluxReconstructorWENO-CU6-M2-LLF.hpp"
#include "flow_model/convective_flux_reconstructor/ConvectiveFluxReconstructorWCNS-JS5-HLLC-HLL.hpp"
#include "flow_model/convective_flux_reconstructor/ConvectiveFluxReconstructorWCNS-CU6-M2-HLLC-HLL.hpp"
#include "flow_model/convective_flux_reconstructor/ConvectiveFluxReconstructorWCNS-HW6-HLLC-HLL.hpp"
#include "flow_model/convective_flux_reconstructor/ConvectiveFluxReconstructorWCNS-HW6-LD-HLLC-HLL.hpp"

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
    boost::shared_ptr<tbox::Database> input_db,
    boost::shared_ptr<geom::CartesianGridGeometry> grid_geom):
        RungeKuttaPatchStrategy(),
        d_object_name(object_name),
        d_dim(dim),
        d_grid_geometry(grid_geom),
#ifdef HAVE_HDF5
        d_visit_writer(NULL),
#endif
        d_plot_context(NULL),
        d_workload_variable(NULL),
        d_use_nonuniform_workload(false),
        d_num_ghosts(hier::IntVector::getZero(d_dim)),
        d_equation_of_state(NULL),
        d_equation_of_state_db(NULL),
        d_conv_flux_reconstructor(NULL),
        d_shock_capturing_scheme_db(NULL),
        d_density(NULL),
        d_partial_density(NULL),
        d_momentum(NULL),
        d_total_energy(NULL),
        d_mass_fraction(NULL),
        d_volume_fraction(NULL),
        d_convective_flux(NULL),
        d_source(NULL),
        d_is_preserving_positivity(false)
{
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(input_db);
    TBOX_ASSERT(grid_geom);
    
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
     * Initialize the d_equation_of_state.
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
    
    if (equation_of_state_string == "IDEAL_GAS")
    {
        switch (d_flow_model)
        {
            case SINGLE_SPECIES:
                d_equation_of_state.reset(new EquationOfStateIdealGas(
                    "ideal gas",
                    d_dim,
                    d_num_species,
                    d_equation_of_state_db,
                    NO_ASSUMPTION));
                break;
            case FOUR_EQN_SHYUE:
                d_equation_of_state.reset(new EquationOfStateIdealGas(
                    "ideal gas",
                    d_dim,
                    d_num_species,
                    d_equation_of_state_db,
                    ISOTHERMAL));
                break;
            case FIVE_EQN_ALLAIRE:
                d_equation_of_state.reset(new EquationOfStateIdealGas(
                    "ideal gas",
                    d_dim,
                    d_num_species,
                    d_equation_of_state_db,
                    ISOBARIC));
                break;
            default:
                TBOX_ERROR(d_object_name
                           << ": "
                           << "Unknown d_flow_model."
                           << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "Unknown equation_of_state string = "
                   << equation_of_state_string
                   << " found in data for Equation_of_state."
                   << std::endl);  
    }
    
    /*
     * Initialize the time-independent variables.
     */
    switch (d_flow_model)
    {
        case SINGLE_SPECIES:
            d_density = boost::shared_ptr<pdat::CellVariable<double> > (
                new pdat::CellVariable<double>(dim, "density", 1));
            break;
        case FOUR_EQN_SHYUE:
            d_density = boost::shared_ptr<pdat::CellVariable<double> > (
                new pdat::CellVariable<double>(dim, "density", 1));
            break;
        case FIVE_EQN_ALLAIRE:
            d_partial_density = boost::shared_ptr<pdat::CellVariable<double> > (
                new pdat::CellVariable<double>(dim, "partial density", d_num_species));
            break;
        default:
            TBOX_ERROR(d_object_name
                       << ": "
                       << "Unknown d_flow_model."
                       << std::endl);
    }
    
    d_momentum = boost::shared_ptr<pdat::CellVariable<double> > (
        new pdat::CellVariable<double>(dim, "momentum", d_dim.getValue()));
    
    d_total_energy = boost::shared_ptr<pdat::CellVariable<double> > (
        new pdat::CellVariable<double>(dim, "total energy", 1));
    
    switch (d_flow_model)
    {
        case SINGLE_SPECIES:
            break;
        case FOUR_EQN_SHYUE:
            d_mass_fraction = boost::shared_ptr<pdat::CellVariable<double> > (
                new pdat::CellVariable<double>(dim, "mass fraction", d_num_species));
            break;
        case FIVE_EQN_ALLAIRE:
            d_volume_fraction = boost::shared_ptr<pdat::CellVariable<double> > (
                new pdat::CellVariable<double>(dim, "volume fraction", d_num_species));
            break;
        default:
            TBOX_ERROR(d_object_name
                       << ": "
                       << "Unknown d_flow_model."
                       << std::endl);
    }
    
    /*
     * Initialize the flux.
     */
    d_convective_flux = boost::shared_ptr<pdat::FaceVariable<double> > (
        new pdat::FaceVariable<double>(dim, "convective flux", d_num_eqn));
    
    /*
     * Initialize the source.
     */
    d_source = boost::shared_ptr<pdat::CellVariable<double> > (
        new pdat::CellVariable<double>(dim, "source", d_num_eqn));
    
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
     * Initialize d_conv_flux_reconstructor.
     */
    if (shock_capturing_scheme_str == "LLF")
    {
        d_conv_flux_reconstructor.reset(new ConvectiveFluxReconstructorLLF(
            "LLF",
            d_dim,
            d_grid_geometry,
            d_flow_model,
            d_num_eqn,
            d_num_species,
            d_equation_of_state,
            d_shock_capturing_scheme_db));
    }
    else if (shock_capturing_scheme_str == "FIRST_ORDER_HLLC")
    {
        d_conv_flux_reconstructor.reset(new ConvectiveFluxReconstructorFirstOrderHLLC(
            "first order HLLC",
            d_dim,
            d_grid_geometry,
            d_flow_model,
            d_num_eqn,
            d_num_species,
            d_equation_of_state,
            d_shock_capturing_scheme_db));
    }
    else if (shock_capturing_scheme_str == "WENO_JS5_LLF")
    {
        d_conv_flux_reconstructor.reset(new ConvectiveFluxReconstructorWENO_JS5_LLF(
            "WENO-JS5-LLF",
            d_dim,
            d_grid_geometry,
            d_flow_model,
            d_num_eqn,
            d_num_species,
            d_equation_of_state,
            d_shock_capturing_scheme_db));
    }
    else if (shock_capturing_scheme_str == "WENO_CU6_M2_LLF")
    {
        d_conv_flux_reconstructor.reset(new ConvectiveFluxReconstructorWENO_CU6_M2_LLF(
            "WENO-CU6-M2-LLF",
            d_dim,
            d_grid_geometry,
            d_flow_model,
            d_num_eqn,
            d_num_species,
            d_equation_of_state,
            d_shock_capturing_scheme_db));
    }
    else if (shock_capturing_scheme_str == "WCNS_JS5_HLLC_HLL")
    {
        d_conv_flux_reconstructor.reset(new ConvectiveFluxReconstructorWCNS_JS5_HLLC_HLL(
            "WCNS-JS5-HLLC-HLL",
            d_dim,
            d_grid_geometry,
            d_flow_model,
            d_num_eqn,
            d_num_species,
            d_equation_of_state,
            d_shock_capturing_scheme_db));
    }
    else if (shock_capturing_scheme_str == "WCNS_CU6_M2_HLLC_HLL")
    {
        d_conv_flux_reconstructor.reset(new ConvectiveFluxReconstructorWCNS_CU6_M2_HLLC_HLL(
            "WCNS-CU6-M2-HLLC-HLL",
            d_dim,
            d_grid_geometry,
            d_flow_model,
            d_num_eqn,
            d_num_species,
            d_equation_of_state,
            d_shock_capturing_scheme_db));
    }
    else if (shock_capturing_scheme_str == "WCNS_HW6_HLLC_HLL")
    {
        d_conv_flux_reconstructor.reset(new ConvectiveFluxReconstructorWCNS_HW6_HLLC_HLL(
            "WCNS-HW6-HLLC-HLL",
            d_dim,
            d_grid_geometry,
            d_flow_model,
            d_num_eqn,
            d_num_species,
            d_equation_of_state,
            d_shock_capturing_scheme_db));
    }
    else if (shock_capturing_scheme_str == "WCNS_HW6_LD_HLLC_HLL")
    {
        d_conv_flux_reconstructor.reset(new ConvectiveFluxReconstructorWCNS_HW6_LD_HLLC_HLL(
            "WCNS-HW6-LD-HLLC-HLL",
            d_dim,
            d_grid_geometry,
            d_flow_model,
            d_num_eqn,
            d_num_species,
            d_equation_of_state,
            d_shock_capturing_scheme_db));
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "Unknown shock_capturing_scheme string = "
                   << shock_capturing_scheme_str
                   << " found in input."
                   << std::endl);        
    }
    
    /*
     * Initialize the number of ghost cells needed.
     */
    d_num_ghosts = d_conv_flux_reconstructor->
        getConvectiveFluxNumberOfGhostCells();
    
    /*
     * Initialize the number of ghost cells and boost::shared_ptr of the variables
     * in the d_conv_flux_reconstructor.
     */
    switch (d_flow_model)
    {
        case SINGLE_SPECIES:
        {
            d_conv_flux_reconstructor->setVariablesForSingleSpecies(
                d_num_ghosts,
                d_density,
                d_momentum,
                d_total_energy,
                d_convective_flux,
                d_source);
            
            break;
        }
        case FOUR_EQN_SHYUE:
        {
            d_conv_flux_reconstructor->setVariablesForFourEqnShyue(
                d_num_ghosts,
                d_density,
                d_momentum,
                d_total_energy,
                d_mass_fraction,
                d_convective_flux,
                d_source);
            
            break;
        }
        case FIVE_EQN_ALLAIRE:
        {
            d_conv_flux_reconstructor->setVariablesForFiveEqnAllaire(
                d_num_ghosts,
                d_partial_density,
                d_momentum,
                d_total_energy,
                d_volume_fraction,
                d_convective_flux,
                d_source);
            
            break;
        }
        default:
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "Unknown d_flow_model."
                       << std::endl);
        }
    }
    
    /*
     * Postprocess boundary data from input/restart values.
     */
    if (d_dim == tbox::Dimension(1))
    {
        // NOT YET IMPLEMENTED
    }
    else if (d_dim == tbox::Dimension(2))
    {
        for (int i = 0; i < NUM_2D_EDGES; i++)
        {
            d_scalar_bdry_edge_conds[i] = d_master_bdry_edge_conds[i];
            d_vector_bdry_edge_conds[i] = d_master_bdry_edge_conds[i];
            
            if (d_master_bdry_edge_conds[i] == BdryCond::REFLECT)
            {
                d_scalar_bdry_edge_conds[i] = BdryCond::FLOW;
            }
        }
     
        for (int i = 0; i < NUM_2D_NODES; i++)
        {
            d_scalar_bdry_node_conds[i] = d_master_bdry_node_conds[i];
            d_vector_bdry_node_conds[i] = d_master_bdry_node_conds[i];
      
            if (d_master_bdry_node_conds[i] == BdryCond::XREFLECT)
            {
                d_scalar_bdry_node_conds[i] = BdryCond::XFLOW;
            }
            if (d_master_bdry_node_conds[i] == BdryCond::YREFLECT)
            {
                d_scalar_bdry_node_conds[i] = BdryCond::YFLOW;
            }
      
            if (d_master_bdry_node_conds[i] != BOGUS_BDRY_DATA)
            {
                d_node_bdry_edge[i] =
                    appu::CartesianBoundaryUtilities2::getEdgeLocationForNodeBdry(
                        i, d_master_bdry_node_conds[i]);
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        for (int i = 0; i < NUM_3D_FACES; i++)
        {
            d_scalar_bdry_face_conds[i] = d_master_bdry_face_conds[i];
            d_vector_bdry_face_conds[i] = d_master_bdry_face_conds[i];
      
            if (d_master_bdry_face_conds[i] == BdryCond::REFLECT)
            {
                d_scalar_bdry_face_conds[i] = BdryCond::FLOW;
            }
        }
        
        for (int i = 0; i < NUM_3D_EDGES; i++)
        {
            d_scalar_bdry_edge_conds[i] = d_master_bdry_edge_conds[i];
            d_vector_bdry_edge_conds[i] = d_master_bdry_edge_conds[i];
      
            if (d_master_bdry_edge_conds[i] == BdryCond::XREFLECT)
            {
                d_scalar_bdry_edge_conds[i] = BdryCond::XFLOW;
            }
            if (d_master_bdry_edge_conds[i] == BdryCond::YREFLECT)
            {
                d_scalar_bdry_edge_conds[i] = BdryCond::YFLOW;
            }
            if (d_master_bdry_edge_conds[i] == BdryCond::ZREFLECT)
            {
                d_scalar_bdry_edge_conds[i] = BdryCond::ZFLOW;
            }
      
            if (d_master_bdry_edge_conds[i] != BOGUS_BDRY_DATA)
            {
                d_edge_bdry_face[i] =
                    appu::CartesianBoundaryUtilities3::getFaceLocationForEdgeBdry(
                        i, d_master_bdry_edge_conds[i]);
            }
        }
    
        for (int i = 0; i < NUM_3D_NODES; i++)
        {
            d_scalar_bdry_node_conds[i] = d_master_bdry_node_conds[i];
            d_vector_bdry_node_conds[i] = d_master_bdry_node_conds[i];
      
            if (d_master_bdry_node_conds[i] == BdryCond::XREFLECT)
            {
                d_scalar_bdry_node_conds[i] = BdryCond::XFLOW;
            }
            if (d_master_bdry_node_conds[i] == BdryCond::YREFLECT)
            {
                d_scalar_bdry_node_conds[i] = BdryCond::YFLOW;
            }
            if (d_master_bdry_node_conds[i] == BdryCond::ZREFLECT)
            {
                d_scalar_bdry_node_conds[i] = BdryCond::ZFLOW;
            }
      
            if (d_master_bdry_node_conds[i] != BOGUS_BDRY_DATA)
            {
                d_node_bdry_face[i] =
                    appu::CartesianBoundaryUtilities3::getFaceLocationForNodeBdry(
                        i, d_master_bdry_node_conds[i]);
            }
        }
    }
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
    
    // Register the time-dependent variables.
    switch (d_flow_model)
    {
        case SINGLE_SPECIES:
        {
            integrator->registerVariable(
                d_density,
                d_num_ghosts,
                RungeKuttaLevelIntegrator::TIME_DEP,
                d_grid_geometry,
                "CONSERVATIVE_COARSEN",
                "CONSERVATIVE_LINEAR_REFINE");
            
            integrator->registerVariable(
                d_momentum,
                d_num_ghosts,
                RungeKuttaLevelIntegrator::TIME_DEP,
                d_grid_geometry,
                "CONSERVATIVE_COARSEN",
                "CONSERVATIVE_LINEAR_REFINE");
            
            integrator->registerVariable(
                d_total_energy,
                d_num_ghosts,
                RungeKuttaLevelIntegrator::TIME_DEP,
                d_grid_geometry,
                "CONSERVATIVE_COARSEN",
                "CONSERVATIVE_LINEAR_REFINE");
            
            break;
        }
        case FOUR_EQN_SHYUE:
        {
            integrator->registerVariable(
                d_density,
                d_num_ghosts,
                RungeKuttaLevelIntegrator::TIME_DEP,
                d_grid_geometry,
                "CONSERVATIVE_COARSEN",
                "CONSERVATIVE_LINEAR_REFINE");
            
            integrator->registerVariable(
                d_momentum,
                d_num_ghosts,
                RungeKuttaLevelIntegrator::TIME_DEP,
                d_grid_geometry,
                "CONSERVATIVE_COARSEN",
                "CONSERVATIVE_LINEAR_REFINE");
            
            integrator->registerVariable(
                d_total_energy,
                d_num_ghosts,
                RungeKuttaLevelIntegrator::TIME_DEP,
                d_grid_geometry,
                "CONSERVATIVE_COARSEN",
                "CONSERVATIVE_LINEAR_REFINE");
            
            integrator->registerVariable(
                d_mass_fraction,
                d_num_ghosts,
                RungeKuttaLevelIntegrator::TIME_DEP,
                d_grid_geometry,
                "CONSERVATIVE_COARSEN",
                "CONSERVATIVE_LINEAR_REFINE");
            
            break;
        }
        case FIVE_EQN_ALLAIRE:
        {
            integrator->registerVariable(
                d_partial_density,
                d_num_ghosts,
                RungeKuttaLevelIntegrator::TIME_DEP,
                d_grid_geometry,
                "CONSERVATIVE_COARSEN",
                "CONSERVATIVE_LINEAR_REFINE");
            
            integrator->registerVariable(
                d_momentum,
                d_num_ghosts,
                RungeKuttaLevelIntegrator::TIME_DEP,
                d_grid_geometry,
                "CONSERVATIVE_COARSEN",
                "CONSERVATIVE_LINEAR_REFINE");
            
            integrator->registerVariable(
                d_total_energy,
                d_num_ghosts,
                RungeKuttaLevelIntegrator::TIME_DEP,
                d_grid_geometry,
                "CONSERVATIVE_COARSEN",
                "CONSERVATIVE_LINEAR_REFINE");
            
            integrator->registerVariable(
                d_volume_fraction,
                d_num_ghosts,
                RungeKuttaLevelIntegrator::TIME_DEP,
                d_grid_geometry,
                "CONSERVATIVE_COARSEN",
                "CONSERVATIVE_LINEAR_REFINE");
            
            break;
        }
        default:
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "Unknown d_flow_model."
                       << std::endl);
        }
    }
    
    // Register the fluxes and sources.
    
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
    
    hier::VariableDatabase* vardb = hier::VariableDatabase::getDatabase();
    
    d_plot_context = integrator->getPlotContext();

#ifdef HAVE_HDF5
    // Register the plotting quantities.
    if (d_visit_writer)
    {
        switch (d_flow_model)
        {
            case SINGLE_SPECIES:
            {
                d_visit_writer->registerPlotQuantity(
                    "density",
                    "SCALAR",
                    vardb->mapVariableAndContextToIndex(
                       d_density,
                       d_plot_context));
                
                d_visit_writer->registerPlotQuantity(
                    "momentum",
                    "VECTOR",
                    vardb->mapVariableAndContextToIndex(
                       d_momentum,
                       d_plot_context));
                
                d_visit_writer->registerPlotQuantity(
                    "total energy",
                    "SCALAR",
                    vardb->mapVariableAndContextToIndex(
                       d_total_energy,
                       d_plot_context));
                
                d_visit_writer->registerDerivedPlotQuantity("pressure",
                    "SCALAR",
                    this);
                
                d_visit_writer->registerDerivedPlotQuantity("sound speed",
                    "SCALAR",
                    this);
                
                d_visit_writer->registerDerivedPlotQuantity("velocity",
                    "VECTOR",
                    this);
                
                break;
            }
            case FOUR_EQN_SHYUE:
            {
                d_visit_writer->registerPlotQuantity(
                    "density",
                    "SCALAR",
                    vardb->mapVariableAndContextToIndex(
                       d_density,
                       d_plot_context));
                
                d_visit_writer->registerPlotQuantity(
                    "momentum",
                    "VECTOR",
                    vardb->mapVariableAndContextToIndex(
                       d_momentum,
                       d_plot_context));
                
                d_visit_writer->registerPlotQuantity(
                    "total energy",
                    "SCALAR",
                    vardb->mapVariableAndContextToIndex(
                       d_total_energy,
                       d_plot_context));
                
                for (int si = 0; si < d_num_species; si++)
                {
                    std::string mass_fraction_name =
                        "mass fraction " + tbox::Utilities::intToString(si);
                        
                    d_visit_writer->registerPlotQuantity(
                        mass_fraction_name,
                        "SCALAR",
                        vardb->mapVariableAndContextToIndex(
                           d_mass_fraction,
                           d_plot_context),
                        si);
                }
                
                d_visit_writer->registerDerivedPlotQuantity("pressure",
                    "SCALAR",
                    this);
                
                d_visit_writer->registerDerivedPlotQuantity("sound speed",
                    "SCALAR",
                    this);
                
                d_visit_writer->registerDerivedPlotQuantity("velocity",
                    "VECTOR",
                    this);
                
                break;
            }
            case FIVE_EQN_ALLAIRE:
            {
                for (int si = 0; si < d_num_species; si++)
                {
                    std::string partial_density_name =
                        "partial density " + tbox::Utilities::intToString(si);
                    
                    d_visit_writer->registerPlotQuantity(
                        partial_density_name,
                        "SCALAR",
                        vardb->mapVariableAndContextToIndex(
                            d_partial_density,
                            d_plot_context),
                        si);
                }
                
                d_visit_writer->registerPlotQuantity(
                    "momentum",
                    "VECTOR",
                    vardb->mapVariableAndContextToIndex(
                       d_momentum,
                       d_plot_context));
                
                d_visit_writer->registerPlotQuantity(
                    "total energy",
                    "SCALAR",
                    vardb->mapVariableAndContextToIndex(
                       d_total_energy,
                       d_plot_context));
                
                for (int si = 0; si < d_num_species; si++)
                {
                    std::string volume_fraction_name =
                        "volume fraction " + tbox::Utilities::intToString(si);
                        
                    d_visit_writer->registerPlotQuantity(
                        volume_fraction_name,
                        "SCALAR",
                        vardb->mapVariableAndContextToIndex(
                           d_volume_fraction,
                           d_plot_context),
                        si);
                }
                
                d_visit_writer->registerDerivedPlotQuantity("pressure",
                    "SCALAR",
                    this);
                
                d_visit_writer->registerDerivedPlotQuantity("sound speed",
                    "SCALAR",
                    this);
                
                d_visit_writer->registerDerivedPlotQuantity("velocity",
                    "VECTOR",
                    this);
                
                d_visit_writer->registerDerivedPlotQuantity("density",
                    "SCALAR",
                    this);
                
                for (int si = 0; si < d_num_species; si++)
                {
                    std::string mass_fraction_name =
                        "mass fraction " + tbox::Utilities::intToString(si);
                        
                    d_visit_writer->registerDerivedPlotQuantity(mass_fraction_name,
                                                                "SCALAR",
                                                                this);
                }
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                           << ": "
                           << "Unknown d_flow_model."
                           << std::endl);
            }
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
    NULL_USE(data_time);
    
    t_init->start();
    
    if (initial_time)
    {
        const double* const domain_xlo = d_grid_geometry->getXLower();
        const double* const domain_xhi = d_grid_geometry->getXUpper();
        
        const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
            BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                patch.getPatchGeometry()));
        
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(patch_geom);
#endif
        
        const double* const dx = patch_geom->getDx();
        const double* const patch_xlo = patch_geom->getXLower();
        
        // Get the dimensions of box that covers the interior of Patch.
        hier::Box patch_box = patch.getBox();
        const hier::IntVector patch_dims = patch_box.numberCells();
        
        switch (d_flow_model)
        {
            case SINGLE_SPECIES:
            {
                boost::shared_ptr<pdat::CellData<double> > density(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_density, getDataContext())));
                
                boost::shared_ptr<pdat::CellData<double> > momentum(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_momentum, getDataContext())));
                
                boost::shared_ptr<pdat::CellData<double> > total_energy(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_total_energy, getDataContext())));
                
#ifdef DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(density);
                TBOX_ASSERT(momentum);
                TBOX_ASSERT(total_energy);
#endif
                
                if (d_dim == tbox::Dimension(1))
                {
                    // NOT YET IMPLEMENTED
                }
                else if (d_dim == tbox::Dimension(2))
                {
                    if (d_project_name == "2D double-Mach reflection")
                    {
                        // Initialize data for a 2D double-Mach reflection problem.
                        double* rho   = density->getPointer(0);
                        double* rho_u = momentum->getPointer(0);
                        double* rho_v = momentum->getPointer(1);
                        double* E     = total_energy->getPointer(0);
                        
                        const double x_0 = 1.0/6.0;
                        
                        const double gamma = 1.4;
                        
                        const double rho_post_shock = 8.0;
                        const double u_post_shock   = 8.25*cos(M_PI/6.0);
                        const double v_post_shock   = -8.25*sin(M_PI/6.0);
                        const double p_post_shock   = 116.5;
                        
                        const double rho_pre_shock  = 1.4;
                        const double u_pre_shock    = 0.0;
                        const double v_pre_shock    = 0.0;
                        const double p_pre_shock    = 1.0;
                        
                        const double rho_u_post_shock = rho_post_shock*u_post_shock;
                        const double rho_v_post_shock = rho_post_shock*v_post_shock;
                        
                        const double rho_u_pre_shock  = rho_pre_shock*u_pre_shock;
                        const double rho_v_pre_shock  = rho_pre_shock*v_pre_shock;
                        
                        const double E_pre_shock = p_pre_shock/(gamma - 1.0) +
                            0.5*rho_pre_shock*(u_pre_shock*u_pre_shock + v_pre_shock*v_pre_shock);
                        
                        const double E_post_shock = p_post_shock/(gamma - 1.0) +
                            0.5*rho_post_shock*(u_post_shock*u_post_shock + v_post_shock*v_post_shock);
                            
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
                                
                                if (x[0] < x_0 + x[1]*sqrt(1.0/3.0) )
                                {
                                    rho[idx_cell] = rho_post_shock;
                                    rho_u[idx_cell] = rho_u_post_shock;
                                    rho_v[idx_cell] = rho_v_post_shock;
                                    E[idx_cell] = E_post_shock;
                                }
                                else
                                {
                                    rho[idx_cell] = rho_pre_shock;
                                    rho_u[idx_cell] = rho_u_pre_shock;
                                    rho_v[idx_cell] = rho_v_pre_shock;
                                    E[idx_cell] = E_pre_shock;
                                }
                            }
                        }
                    }
                    else if (d_project_name == "2D wedge flow")
                    {
                        // Initialize data for a 2D density wave advection problem.
                        double* rho   = density->getPointer(0);
                        double* rho_u = momentum->getPointer(0);
                        double* rho_v = momentum->getPointer(1);
                        double* E     = total_energy->getPointer(0);
                        
                        const double gamma = 1.4;
                        const double R     = 287.058;
                        
                        const double p_inf = 1e5;
                        const double T_inf = 300;
                        const double M_inf = 2;
                        const double theta = 10.0/180*M_PI;
                        
                        const double U_inf   = M_inf * sqrt(gamma*R*T_inf);
                        const double rho_inf = p_inf/R/T_inf;
                        
                        const double u_inf = U_inf*cos(theta);
                        const double v_inf = -U_inf*sin(theta);
                        
                        for (int j = 0; j < patch_dims[1]; j++)
                        {
                            for (int i = 0; i < patch_dims[0]; i++)
                            {
                                // Compute index into linear data array.
                                int idx_cell = i + j*patch_dims[0];
                                
                                rho[idx_cell]   = rho_inf;
                                rho_u[idx_cell] = rho_inf*u_inf;
                                rho_v[idx_cell] = rho_inf*v_inf;
                                E[idx_cell]     = p_inf/(gamma - 1.0) + 0.5*rho_inf*(
                                    u_inf*u_inf + v_inf*v_inf);
                            }
                        }
                    }
                    else if (d_project_name == "Advection of 2D density wave")
                    {
                        // Initialize data for a 2D density wave advection problem.
                        double* rho   = density->getPointer(0);
                        double* rho_u = momentum->getPointer(0);
                        double* rho_v = momentum->getPointer(1);
                        double* E     = total_energy->getPointer(0);
                        
                        const double x_a = 1.0/3*(domain_xlo[0] + domain_xhi[0]);
                        const double x_b = 2.0/3*(domain_xlo[0] + domain_xhi[0]);
                        
                        const double y_a = 1.0/3*(domain_xlo[1] + domain_xhi[1]);
                        const double y_b = 2.0/3*(domain_xlo[1] + domain_xhi[1]);
                        
                        double gamma = d_equation_of_state->getSpeciesThermodynamicProperty(
                            "gamma",
                            0);
                        
                        // Initial conditions inside the square.
                        double rho_i = 10.0;
                        double u_i   = 1.0;
                        double v_i   = 1.0;
                        double p_i   = 1.0;
                        
                        // Initial conditions outside the square.
                        double rho_o = 1.0;
                        double u_o   = 1.0;
                        double v_o   = 1.0;
                        double p_o   = 1.0;
                        
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
                                
                                if ((x[0] >= x_a) &&
                                    (x[0] <= x_b) &&
                                    (x[1] >= y_a) &&
                                    (x[1] <= y_b))
                                {
                                    rho[idx_cell]   = rho_i;
                                    rho_u[idx_cell] = rho_i*u_i;
                                    rho_v[idx_cell] = rho_i*v_i;
                                    E[idx_cell]     = p_i/(gamma - 1.0) + 0.5*rho_i*(u_i*u_i +
                                        v_i*v_i);
                                }
                                else
                                {
                                    rho[idx_cell]   = rho_o;
                                    rho_u[idx_cell] = rho_o*u_o;
                                    rho_v[idx_cell] = rho_o*v_o;
                                    E[idx_cell]     = p_o/(gamma - 1.0) + 0.5*rho_o*(u_o*u_o +
                                        v_o*v_o);
                                }
                            }
                        }
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name
                                   << ": "
                                   << "Cannot initialize data for unknown problem."
                                   << std::endl);
                    }
                }
                else if (d_dim == tbox::Dimension(3))
                {
                    if (d_project_name == "3D inviscid Taylor-Green vortex")
                    {
                        // Initialize data for a 3D inviscid Taylor-Green problem.
                        double* rho   = density->getPointer(0);
                        double* rho_u = momentum->getPointer(0);
                        double* rho_v = momentum->getPointer(1);
                        double* rho_w = momentum->getPointer(2);
                        double* E     = total_energy->getPointer(0);
                        
                        double gamma = d_equation_of_state->getSpeciesThermodynamicProperty(
                            "gamma",
                            0);
                        
                        for (int k = 0; k < patch_dims[2]; k++)
                        {
                            for (int j = 0; j < patch_dims[1]; j++)
                            {
                                for (int i = 0; i < patch_dims[0]; i++)
                                {
                                    // Compute index into linear data array.
                                    int idx_cell = i +
                                        j*patch_dims[0] +
                                        k*patch_dims[0]*patch_dims[1];
                                    
                                    // Compute the coordinates.
                                    double x[3];
                                    x[0] = patch_xlo[0] + (i + 0.5)*dx[0];
                                    x[1] = patch_xlo[1] + (j + 0.5)*dx[1];
                                    x[2] = patch_xlo[2] + (k + 0.5)*dx[2];
                                    
                                    // Compute density, velocities and pressure.
                                    const double rho_i = 1.0;
                                    const double u_i = sin(x[0])*cos(x[1])*cos(x[2]);
                                    const double v_i = -cos(x[0])*sin(x[1])*cos(x[2]);
                                    const double w_i = 0.0;
                                    const double p_i = 100.0 +
                                        ((cos(2*x[2]) + 2.0)*(cos(2*x[0]) + cos(2*x[1])) - 2.0)/
                                            16.0;
                                    
                                    rho[idx_cell]   = rho_i;
                                    rho_u[idx_cell] = rho_i*u_i;
                                    rho_v[idx_cell] = rho_i*v_i;
                                    rho_w[idx_cell] = rho_i*w_i;
                                    E[idx_cell]     = p_i/(gamma - 1.0) + 0.5*rho_i*(u_i*u_i +
                                        v_i*v_i + w_i*w_i);
                                }
                            }
                        }
                    }
                    else if (d_project_name == "Advection of 3D density wave")
                    {
                        // Initialize data for a 3D density wave advection problem.
                        double* rho   = density->getPointer(0);
                        double* rho_u = momentum->getPointer(0);
                        double* rho_v = momentum->getPointer(1);
                        double* rho_w = momentum->getPointer(2);
                        double* E     = total_energy->getPointer(0);
                        
                        const double x_a = 1.0/3*(domain_xlo[0] + domain_xhi[0]);
                        const double x_b = 2.0/3*(domain_xlo[0] + domain_xhi[0]);
                        
                        const double y_a = 1.0/3*(domain_xlo[1] + domain_xhi[1]);
                        const double y_b = 2.0/3*(domain_xlo[1] + domain_xhi[1]);
                        
                        const double z_a = 1.0/3*(domain_xlo[2] + domain_xhi[2]);
                        const double z_b = 2.0/3*(domain_xlo[2] + domain_xhi[2]);
                        
                        double gamma = d_equation_of_state->getSpeciesThermodynamicProperty(
                            "gamma",
                            0);
                        
                        // Initial conditions inside the cube.
                        double rho_i = 10.0;
                        double u_i   = 1.0;
                        double v_i   = 1.0;
                        double w_i   = 1.0;
                        double p_i   = 1.0;
                        
                        // Initial conditions outside the cube.
                        double rho_o = 1.0;
                        double u_o   = 1.0;
                        double v_o   = 1.0;
                        double w_o   = 1.0;
                        double p_o   = 1.0;
                        
                        for (int k = 0; k < patch_dims[2]; k++)
                        {
                            for (int j = 0; j < patch_dims[1]; j++)
                            {
                                for (int i = 0; i < patch_dims[0]; i++)
                                {
                                    // Compute index into linear data array.
                                    int idx_cell = i +
                                        j*patch_dims[0] +
                                        k*patch_dims[0]*patch_dims[1];
                                    
                                    // Compute the coordinates.
                                    double x[3];
                                    x[0] = patch_xlo[0] + (i + 0.5)*dx[0];
                                    x[1] = patch_xlo[1] + (j + 0.5)*dx[1];
                                    x[2] = patch_xlo[2] + (k + 0.5)*dx[2];
                                    
                                    if ((x[0] >= x_a) &&
                                        (x[0] <= x_b) &&
                                        (x[1] >= y_a) &&
                                        (x[1] <= y_b) &&
                                        (x[2] >= z_a) &&
                                        (x[2] <= z_b))
                                    {
                                        rho[idx_cell]   = rho_i;
                                        rho_u[idx_cell] = rho_i*u_i;
                                        rho_v[idx_cell] = rho_i*v_i;
                                        rho_w[idx_cell] = rho_i*w_i;
                                        E[idx_cell]     = p_i/(gamma - 1.0) + 0.5*rho_i*(u_i*u_i +
                                            v_i*v_i + w_i*w_i);
                                    }
                                    else
                                    {
                                        rho[idx_cell]   = rho_o;
                                        rho_u[idx_cell] = rho_o*u_o;
                                        rho_v[idx_cell] = rho_o*v_o;
                                        rho_w[idx_cell] = rho_o*w_o;
                                        E[idx_cell]     = p_o/(gamma - 1.0) + 0.5*rho_o*(u_o*u_o +
                                            v_o*v_o + w_i*w_i);
                                    }
                                }
                            }
                        }
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name
                                   << ": "
                                   << "Cannot initialize data for unknown problem."
                                   << std::endl);
                    }
                }
                
                break;
            }
            case FOUR_EQN_SHYUE:
            {
                boost::shared_ptr<pdat::CellData<double> > density(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_density, getDataContext())));
                
                boost::shared_ptr<pdat::CellData<double> > momentum(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_momentum, getDataContext())));
                
                boost::shared_ptr<pdat::CellData<double> > total_energy(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_total_energy, getDataContext())));
                
                boost::shared_ptr<pdat::CellData<double> > mass_fraction(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_mass_fraction, getDataContext())));
                
#ifdef DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(density);
                TBOX_ASSERT(momentum);
                TBOX_ASSERT(total_energy);
                TBOX_ASSERT(mass_fraction);
#endif
                
                if (d_dim == tbox::Dimension(1))
                {
                    // NOT YET IMPLEMENTED
                }
                else if (d_dim == tbox::Dimension(2))
                {
                    if (d_project_name == "Advection of 2D material interface")
                    {
                        // Initialize data for a 2D material interface advection problem.
                        
                        if (d_num_species != 2)
                        {
                            TBOX_ERROR(d_object_name
                                       << ": "
                                       << "Please provide only two-species for multi-species simulation testing."
                                       << std::endl);
                        }
                        
                        double* rho       = density->getPointer(0);
                        double* rho_u     = momentum->getPointer(0);
                        double* rho_v     = momentum->getPointer(1);
                        double* E         = total_energy->getPointer(0);
                        double* Y_1       = mass_fraction->getPointer(0);
                        double* Y_2       = mass_fraction->getPointer(1);
                        
                        const double x_a = 1.0/3*(domain_xlo[0] + domain_xhi[0]);
                        const double x_b = 2.0/3*(domain_xlo[0] + domain_xhi[0]);
                        
                        const double y_a = 1.0/3*(domain_xlo[1] + domain_xhi[1]);
                        const double y_b = 2.0/3*(domain_xlo[1] + domain_xhi[1]);
                        
                        // Material initial conditions.
                        double gamma_m  = d_equation_of_state->
                            getSpeciesThermodynamicProperty(
                                "gamma",
                                0);
                        double rho_m    = 10.0;
                        double u_m      = 0.5;
                        double v_m      = 0.5;
                        double p_m      = 1.0/1.4;
                        
                        // Ambient initial conditions.
                        double gamma_a  = d_equation_of_state->
                            getSpeciesThermodynamicProperty(
                                "gamma",
                                1);
                        double rho_a    = 1.0;
                        double u_a      = 0.5;
                        double v_a      = 0.5;
                        double p_a      = 1.0/1.4;
                        
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
                                
                                if ((x[0] >= x_a) &&
                                    (x[0] <= x_b) &&
                                    (x[1] >= y_a) &&
                                    (x[1] <= y_b))
                                {
                                    rho[idx_cell]     = rho_m;
                                    rho_u[idx_cell]   = rho_m*u_m;
                                    rho_v[idx_cell]   = rho_m*v_m;
                                    E[idx_cell]       = p_m/(gamma_m - 1.0) +
                                        0.5*rho_m*(u_m*u_m + v_m*v_m);
                                    Y_1[idx_cell]     = 1.0;
                                    Y_2[idx_cell]     = 0.0;
                                }
                                else
                                {
                                    rho[idx_cell]     = rho_a;
                                    rho_u[idx_cell]   = rho_a*u_a;
                                    rho_v[idx_cell]   = rho_a*v_a;
                                    E[idx_cell]       = p_a/(gamma_a - 1.0) +
                                        0.5*rho_a*(u_a*u_a + v_a*v_a);
                                    Y_1[idx_cell]     = 0.0;
                                    Y_2[idx_cell]     = 1.0;
                                }
                            }
                        }
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name
                                   << ": "
                                   << "Cannot initialize data for unknown problem."
                                   << std::endl);
                    }
                }
                else if (d_dim == tbox::Dimension(3))
                {
                    if (d_project_name == "Advection of 3D material interface")
                    {
                        // Initialize data for a 3D material interface advection problem.
                        
                        if (d_num_species != 2)
                        {
                            TBOX_ERROR(d_object_name
                                       << ": "
                                       << "Please provide only two-species for multi-species simulation testing."
                                       << std::endl);
                        }
                        
                        double* rho       = density->getPointer(0);
                        double* rho_u     = momentum->getPointer(0);
                        double* rho_v     = momentum->getPointer(1);
                        double* rho_w     = momentum->getPointer(2);
                        double* E         = total_energy->getPointer(0);
                        double* Y_1       = mass_fraction->getPointer(0);
                        double* Y_2       = mass_fraction->getPointer(1);
                        
                        const double x_a = 1.0/3*(domain_xlo[0] + domain_xhi[0]);
                        const double x_b = 2.0/3*(domain_xlo[0] + domain_xhi[0]);
                        
                        const double y_a = 1.0/3*(domain_xlo[1] + domain_xhi[1]);
                        const double y_b = 2.0/3*(domain_xlo[1] + domain_xhi[1]);
                        
                        const double z_a = 1.0/3*(domain_xlo[2] + domain_xhi[2]);
                        const double z_b = 2.0/3*(domain_xlo[2] + domain_xhi[2]);
                        
                        // Material initial conditions.
                        double gamma_m  = d_equation_of_state->
                            getSpeciesThermodynamicProperty(
                                "gamma",
                                0);
                        double rho_m    = 10.0;
                        double u_m      = 0.5;
                        double v_m      = 0.5;
                        double w_m      = 0.5;
                        double p_m      = 1.0/1.4;
                        
                        // Ambient initial conditions.
                        double gamma_a  = d_equation_of_state->
                            getSpeciesThermodynamicProperty(
                                "gamma",
                                1);
                        double rho_a    = 1.0;
                        double u_a      = 0.5;
                        double v_a      = 0.5;
                        double w_a      = 0.5;
                        double p_a      = 1.0/1.4;
                        
                        for (int k = 0; k < patch_dims[2]; k++)
                        {
                            for (int j = 0; j < patch_dims[1]; j++)
                            {
                                for (int i = 0; i < patch_dims[0]; i++)
                                {
                                    // Compute index into linear data array.
                                    int idx_cell = i +
                                        j*patch_dims[0] +
                                        k*patch_dims[0]*patch_dims[1];
                                    
                                    // Compute the coordinates.
                                    double x[3];
                                    x[0] = patch_xlo[0] + (i + 0.5)*dx[0];
                                    x[1] = patch_xlo[1] + (j + 0.5)*dx[1];
                                    x[2] = patch_xlo[2] + (k + 0.5)*dx[2];
                                    
                                    if ((x[0] >= x_a) &&
                                        (x[0] <= x_b) &&
                                        (x[1] >= y_a) &&
                                        (x[1] <= y_b) &&
                                        (x[2] >= z_a) &&
                                        (x[2] <= z_b))
                                    {
                                        rho[idx_cell]     = rho_m;
                                        rho_u[idx_cell]   = rho_m*u_m;
                                        rho_v[idx_cell]   = rho_m*v_m;
                                        rho_w[idx_cell]   = rho_m*w_m;
                                        E[idx_cell]       = p_m/(gamma_m - 1.0) +
                                            0.5*rho_m*(u_m*u_m + v_m*v_m + w_m*w_m);
                                        Y_1[idx_cell]     = 1.0;
                                        Y_2[idx_cell]     = 0.0;
                                    }
                                    else
                                    {
                                        rho[idx_cell]     = rho_a;
                                        rho_u[idx_cell]   = rho_a*u_a;
                                        rho_v[idx_cell]   = rho_a*v_a;
                                        rho_w[idx_cell]   = rho_a*w_a;
                                        E[idx_cell]       = p_a/(gamma_a - 1.0) +
                                            0.5*rho_a*(u_a*u_a + v_a*v_a + w_a*w_a);
                                        Y_1[idx_cell]     = 0.0;
                                        Y_2[idx_cell]     = 1.0;
                                    }
                                }
                            }
                        }
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name
                                   << ": "
                                   << "Cannot initialize data for unknown problem."
                                   << std::endl);
                    }
                }
                
                break;
            }
            case FIVE_EQN_ALLAIRE:
            {
                boost::shared_ptr<pdat::CellData<double> > partial_density(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_partial_density, getDataContext())));
                
                boost::shared_ptr<pdat::CellData<double> > momentum(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_momentum, getDataContext())));
                
                boost::shared_ptr<pdat::CellData<double> > total_energy(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_total_energy, getDataContext())));
                
                boost::shared_ptr<pdat::CellData<double> > volume_fraction(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_volume_fraction, getDataContext())));
                
#ifdef DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(partial_density);
                TBOX_ASSERT(momentum);
                TBOX_ASSERT(total_energy);
                TBOX_ASSERT(volume_fraction);
#endif
                
                if (d_dim == tbox::Dimension(1))
                {
                    // NOT YET IMPLEMENTED
                }
                else if (d_dim == tbox::Dimension(2))
                {
                    if (d_project_name == "2D shock-bubble interaction")
                    {
                        if (d_num_species != 2)
                        {
                            TBOX_ERROR(d_object_name
                                       << ": "
                                       << "Please provide only two-species for the 2D shock-bubble interaction simulation."
                                       << std::endl);
                        }
                        
                        const double D = 1.0;
                        
                        double* Z_rho_1   = partial_density->getPointer(0);
                        double* Z_rho_2   = partial_density->getPointer(1);
                        double* rho_u     = momentum->getPointer(0);
                        double* rho_v     = momentum->getPointer(1);
                        double* E         = total_energy->getPointer(0);
                        double* Z_1       = volume_fraction->getPointer(0);
                        double* Z_2       = volume_fraction->getPointer(1);
                        
                        // species 0: He
                        // species 1: air
                        const double gamma_0 = d_equation_of_state->
                            getSpeciesThermodynamicProperty(
                                "gamma",
                                0);
                        
                        const double gamma_1 = d_equation_of_state->
                            getSpeciesThermodynamicProperty(
                                "gamma",
                                1);
                        
                        NULL_USE(gamma_0);
                        
                        // He, pre-shock condition.
                        const double rho_He = 0.1819;
                        const double u_He   = 0.0;
                        const double v_He   = 0.0;
                        const double p_He   = 1.0/1.4;
                        const double Z_He   = 1.0;
                        
                        // air, pre-shock condition.
                        const double rho_pre = 1.0;
                        const double u_pre   = 0.0;
                        const double v_pre   = 0.0;
                        const double p_pre   = 1.0/1.4;
                        const double Z_pre   = 0.0;
                        
                        // air, post-shock condition.
                        const double rho_post = 1.3764;
                        const double u_post   = -0.3336;
                        const double v_post   = 0.0;
                        const double p_post   = 1.5698/1.4;
                        const double Z_post   = 0.0;
                        
                        // Compute the characteristic length of the initial interface thickness.
                        const double epsilon_i = (3.0)*sqrt(dx[0]*dx[1]);
                        
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
                                
                                if (x[0] > 4.5*D)
                                {
                                    Z_rho_1[idx_cell] = 0.0;
                                    Z_rho_2[idx_cell] = rho_post;
                                    rho_u[idx_cell]   = rho_post*u_post;
                                    rho_v[idx_cell]   = rho_post*v_post;
                                    E[idx_cell]       = p_post/(gamma_1 - 1.0) +
                                        0.5*rho_post*(u_post*u_post + v_post*v_post);
                                    Z_1[idx_cell]     = Z_post;
                                    Z_2[idx_cell]     = 1.0 - Z_post;
                                }
                                else
                                {
                                    // Compute the distance from the initial material interface.
                                    const double dR = sqrt(pow(x[0] - 3.5, 2) + x[1]*x[1]) - 0.5*D;
                                    
                                    const double f_sm = 0.5*(1.0 + erf(dR/epsilon_i));
                                    
                                    // Smooth the primitive quantity.
                                    const double Z_rho_1_i = rho_He*(1.0 - f_sm);
                                    const double Z_rho_2_i = rho_pre*f_sm;
                                    const double u_i       = u_He*(1.0 - f_sm) + u_pre*f_sm;
                                    const double v_i       = v_He*(1.0 - f_sm) + v_pre*f_sm;
                                    const double p_i       = p_He*(1.0 - f_sm) + p_pre*f_sm;
                                    const double Z_1_i     = Z_He*(1.0 - f_sm) + Z_pre*f_sm;
                                    
                                    const double rho_i = Z_rho_1_i + Z_rho_2_i;
                                    const double Z_2_i = 1.0 - Z_1_i;
                                    
                                    std::vector<const double*> Z_ptr;
                                    Z_ptr.push_back(&Z_1_i);
                                    Z_ptr.push_back(&Z_2_i);
                                    
                                    const double gamma = d_equation_of_state->
                                        getMixtureThermodynamicPropertyWithVolumeFraction(
                                            "gamma",
                                            Z_ptr);
                                    
                                    Z_rho_1[idx_cell] = Z_rho_1_i;
                                    Z_rho_2[idx_cell] = Z_rho_2_i;
                                    rho_u[idx_cell]   = rho_i*u_i;
                                    rho_v[idx_cell]   = rho_i*v_i;
                                    E[idx_cell]       = p_i/(gamma - 1.0) +
                                        0.5*rho_i*(u_i*u_i + v_i*v_i);
                                    Z_1[idx_cell]     = Z_1_i;
                                    Z_2[idx_cell]     = Z_2_i;
                                }
                            }
                        }
                    }
                    /*
                    if (d_project_name == "2D shock-bubble interaction")
                    {
                        if (d_num_species != 2)
                        {
                            TBOX_ERROR(d_object_name
                                       << ": "
                                       << "Please provide only two-species for the 2D shock-bubble interaction simulation."
                                       << std::endl);
                        }
                        
                        const double D = 1.0;
                        
                        double* Z_rho_1   = partial_density->getPointer(0);
                        double* Z_rho_2   = partial_density->getPointer(1);
                        double* rho_u     = momentum->getPointer(0);
                        double* rho_v     = momentum->getPointer(1);
                        double* E         = total_energy->getPointer(0);
                        double* Z_1       = volume_fraction->getPointer(0);
                        double* Z_2       = volume_fraction->getPointer(1);
                        
                        // species 0: He
                        // species 1: air
                        const double gamma_0 = d_equation_of_state->
                            getSpeciesThermodynamicProperty(
                                "gamma",
                                0);
                        
                        const double gamma_1 = d_equation_of_state->
                            getSpeciesThermodynamicProperty(
                                "gamma",
                                1);
                        
                        // He, pre-shock condition.
                        const double rho_He = 0.1819;
                        const double u_He   = 0.0;
                        const double v_He   = 0.0;
                        const double p_He   = 1.0/1.4;
                        const double Z_He   = 1.0;
                        
                        // air, pre-shock condition.
                        const double rho_pre = 1.0;
                        const double u_pre   = 0.0;
                        const double v_pre   = 0.0;
                        const double p_pre   = 1.0/1.4;
                        const double Z_pre   = 0.0;
                        
                        // air, post-shock condition.
                        const double rho_post = 1.3764;
                        const double u_post   = -0.3336;
                        const double v_post   = 0.0;
                        const double p_post   = 1.5698/1.4;
                        const double Z_post   = 0.0;
                        
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
                                
                                if (x[0] > 4.5*D)
                                {
                                    Z_rho_1[idx_cell] = 0.0;
                                    Z_rho_2[idx_cell] = rho_post;
                                    rho_u[idx_cell]   = rho_post*u_post;
                                    rho_v[idx_cell]   = rho_post*v_post;
                                    E[idx_cell]       = p_post/(gamma_1 - 1.0) +
                                        0.5*rho_post*(u_post*u_post + v_post*v_post);
                                    Z_1[idx_cell]     = Z_post;
                                    Z_2[idx_cell]     = 1.0 - Z_post;
                                }
                                else if (sqrt(pow(x[0] - 3.5, 2) + x[1]*x[1]) < 0.5*D)
                                {
                                    Z_rho_1[idx_cell] = rho_He;
                                    Z_rho_2[idx_cell] = 0.0;
                                    rho_u[idx_cell]   = rho_He*u_He;
                                    rho_v[idx_cell]   = rho_He*v_He;
                                    E[idx_cell]       = p_He/(gamma_0 - 1.0) +
                                        0.5*rho_He*(u_He*u_He + v_He*v_He);
                                    Z_1[idx_cell]     = Z_He;
                                    Z_2[idx_cell]     = 1.0 - Z_He;
                                }
                                else
                                {
                                    Z_rho_1[idx_cell] = 0.0;
                                    Z_rho_2[idx_cell] = rho_pre;
                                    rho_u[idx_cell]   = rho_pre*u_pre;
                                    rho_v[idx_cell]   = rho_pre*v_pre;
                                    E[idx_cell]       = p_pre/(gamma_1 - 1.0) +
                                        0.5*rho_pre*(u_pre*u_pre + v_pre*v_pre);
                                    Z_1[idx_cell]     = Z_pre;
                                    Z_2[idx_cell]     = 1.0 - Z_pre;
                                }
                            }
                        }
                    }
                    */
                    else if (d_project_name == "2D Richtmyer-Meshkov instability")
                    {
                        // Initialize data for a 2D Richtmyer-Meshkov instability.
                        
                        if (d_num_species != 2)
                        {
                            TBOX_ERROR(d_object_name
                                       << ": "
                                       << "Please provide only two-species for multi-species simulation testing."
                                       << std::endl);
                        }
                        
                        const double D = 1.0;
                        
                        double* Z_rho_1   = partial_density->getPointer(0);
                        double* Z_rho_2   = partial_density->getPointer(1);
                        double* rho_u     = momentum->getPointer(0);
                        double* rho_v     = momentum->getPointer(1);
                        double* E         = total_energy->getPointer(0);
                        double* Z_1       = volume_fraction->getPointer(0);
                        double* Z_2       = volume_fraction->getPointer(1);
                        
                        // species 0: SF6
                        // species 1: air
                        const double gamma_0 = d_equation_of_state->
                            getSpeciesThermodynamicProperty(
                                "gamma",
                                0);
                        
                        const double gamma_1 = d_equation_of_state->
                            getSpeciesThermodynamicProperty(
                                "gamma",
                                1);
                        
                        NULL_USE(gamma_0);
                        
                        // SF6, pre-shock condition.
                        const double rho_SF6 = 5.04;
                        const double u_SF6   = 1.24;
                        const double v_SF6   = 0.0;
                        const double p_SF6   = 1.0/1.4;
                        const double Z_SF6   = 1.0;
                        
                        // air, pre-shock condition.
                        const double rho_pre = 1.0;
                        const double u_pre   = 1.24;
                        const double v_pre   = 0.0;
                        const double p_pre   = 1.0/1.4;
                        const double Z_pre   = 0.0;
                        
                        // air, post-shock condition.
                        const double rho_post = 1.4112;
                        const double u_post   = 0.8787;
                        const double v_post   = 0.0;
                        const double p_post   = 1.6272/1.4;
                        const double Z_post   = 0.0;
                        
                        // Compute the characteristic length of the initial interface thickness.
                        const double epsilon_i = (6.0)*sqrt(dx[0]*dx[1]);
                        
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
                                
                                if (x[0] > 0.7*D)
                                {
                                    Z_rho_1[idx_cell] = 0.0;
                                    Z_rho_2[idx_cell] = rho_post;
                                    rho_u[idx_cell]   = rho_post*u_post;
                                    rho_v[idx_cell]   = rho_post*v_post;
                                    E[idx_cell]       = p_post/(gamma_1 - 1.0) +
                                        0.5*rho_post*(u_post*u_post + v_post*v_post);
                                    Z_1[idx_cell]     = Z_post;
                                    Z_2[idx_cell]     = 1.0 - Z_post;
                                }
                                else
                                {
                                    // Compute the distance from the initial material interface.
                                    const double dR = x[0] - (0.4 - 0.1*sin(2*M_PI*(x[1]/D + 0.25)))*D;
                                    
                                    const double f_sm = 0.5*(1.0 + erf(dR/epsilon_i));
                                    
                                    // Smooth the primitive quantity.
                                    const double Z_rho_1_i = rho_SF6*(1.0 - f_sm);
                                    const double Z_rho_2_i = rho_pre*f_sm;
                                    const double u_i       = u_SF6*(1.0 - f_sm) + u_pre*f_sm;
                                    const double v_i       = v_SF6*(1.0 - f_sm) + v_pre*f_sm;
                                    const double p_i       = p_SF6*(1.0 - f_sm) + p_pre*f_sm;
                                    const double Z_1_i     = Z_SF6*(1.0 - f_sm) + Z_pre*f_sm;
                                    
                                    const double rho_i = Z_rho_1_i + Z_rho_2_i;
                                    const double Z_2_i = 1.0 - Z_1_i;
                                    
                                    std::vector<const double*> Z_ptr;
                                    Z_ptr.push_back(&Z_1_i);
                                    Z_ptr.push_back(&Z_2_i);
                                    
                                    const double gamma = d_equation_of_state->
                                        getMixtureThermodynamicPropertyWithVolumeFraction(
                                            "gamma",
                                            Z_ptr);
                                    
                                    Z_rho_1[idx_cell] = Z_rho_1_i;
                                    Z_rho_2[idx_cell] = Z_rho_2_i;
                                    rho_u[idx_cell]   = rho_i*u_i;
                                    rho_v[idx_cell]   = rho_i*v_i;
                                    E[idx_cell]       = p_i/(gamma - 1.0) +
                                        0.5*rho_i*(u_i*u_i + v_i*v_i);
                                    Z_1[idx_cell]     = Z_1_i;
                                    Z_2[idx_cell]     = Z_2_i;
                                }
                            }
                        }
                    }
                    /*
                    else if (d_project_name == "2D Richtmyer-Meshkov instability")
                    {
                        // Initialize data for a 2D Richtmyer-Meshkov instability.
                        
                        if (d_num_species != 2)
                        {
                            TBOX_ERROR(d_object_name
                                       << ": "
                                       << "Please provide only two-species for multi-species simulation testing."
                                       << std::endl);
                        }
                        
                        const double D = 1.0;
                        
                        double* Z_rho_1   = partial_density->getPointer(0);
                        double* Z_rho_2   = partial_density->getPointer(1);
                        double* rho_u     = momentum->getPointer(0);
                        double* rho_v     = momentum->getPointer(1);
                        double* E         = total_energy->getPointer(0);
                        double* Z_1       = volume_fraction->getPointer(0);
                        double* Z_2       = volume_fraction->getPointer(1);
                        
                        // species 0: SF6
                        // species 1: air
                        const double gamma_0 = d_equation_of_state->
                            getSpeciesThermodynamicProperty(
                                "gamma",
                                0);
                        
                        const double gamma_1 = d_equation_of_state->
                            getSpeciesThermodynamicProperty(
                                "gamma",
                                1);
                        
                        // SF6, pre-shock condition.
                        const double rho_SF6 = 5.04;
                        const double u_SF6   = 1.24;
                        const double v_SF6   = 0.0;
                        const double p_SF6   = 1.0/1.4;
                        const double Z_SF6   = 1.0;
                        
                        // air, pre-shock condition.
                        const double rho_pre = 1.0;
                        const double u_pre   = 1.24;
                        const double v_pre   = 0.0;
                        const double p_pre   = 1.0/1.4;
                        const double Z_pre   = 0.0;
                        
                        // air, post-shock condition.
                        const double rho_post = 1.4112;
                        const double u_post   = 0.8787;
                        const double v_post   = 0.0;
                        const double p_post   = 1.6272/1.4;
                        const double Z_post   = 0.0;
                        
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
                                
                                if (x[0] > 0.7*D)
                                {
                                    Z_rho_1[idx_cell] = 0.0;
                                    Z_rho_2[idx_cell] = rho_post;
                                    rho_u[idx_cell]   = rho_post*u_post;
                                    rho_v[idx_cell]   = rho_post*v_post;
                                    E[idx_cell]       = p_post/(gamma_1 - 1.0) +
                                        0.5*rho_post*(u_post*u_post + v_post*v_post);
                                    Z_1[idx_cell]     = Z_post;
                                    Z_2[idx_cell]     = 1.0 - Z_post;
                                }
                                else if (x[0] < (0.4 - 0.1*sin(2*M_PI*(x[1]/D + 0.25)))*D)
                                {
                                    Z_rho_1[idx_cell] = rho_SF6;
                                    Z_rho_2[idx_cell] = 0.0;
                                    rho_u[idx_cell]   = rho_SF6*u_SF6;
                                    rho_v[idx_cell]   = rho_SF6*v_SF6;
                                    E[idx_cell]       = p_SF6/(gamma_0 - 1.0) +
                                        0.5*rho_SF6*(u_SF6*u_SF6 + v_SF6*v_SF6);
                                    Z_1[idx_cell]     = Z_SF6;
                                    Z_2[idx_cell]     = 1.0 - Z_SF6;
                                }
                                else
                                {
                                    Z_rho_1[idx_cell] = 0.0;
                                    Z_rho_2[idx_cell] = rho_pre;
                                    rho_u[idx_cell]   = rho_pre*u_pre;
                                    rho_v[idx_cell]   = rho_pre*v_pre;
                                    E[idx_cell]       = p_pre/(gamma_1 - 1.0) +
                                        0.5*rho_pre*(u_pre*u_pre + v_pre*v_pre);
                                    Z_1[idx_cell]     = Z_pre;
                                    Z_2[idx_cell]     = 1.0 - Z_pre;
                                }
                            }
                        }
                    }
                    */
                    /*
                    else if (d_project_name == "2D Richtmyer-Meshkov instability")
                    {
                        // Initialize data for a 2D Richtmyer-Meshkov instability.
                        
                        if (d_num_species != 2)
                        {
                            TBOX_ERROR(d_object_name
                                       << ": "
                                       << "Please provide only two-species for multi-species simulation testing."
                                       << std::endl);
                        }
                        
                        const double D = 1.0;
                        
                        double* Z_rho_1   = partial_density->getPointer(0);
                        double* Z_rho_2   = partial_density->getPointer(1);
                        double* rho_u     = momentum->getPointer(0);
                        double* rho_v     = momentum->getPointer(1);
                        double* E         = total_energy->getPointer(0);
                        double* Z_1       = volume_fraction->getPointer(0);
                        double* Z_2       = volume_fraction->getPointer(1);
                        
                        // species 0: SF6
                        // species 1: air
                        const double gamma_0 = d_equation_of_state->
                            getSpeciesThermodynamicProperty(
                                "gamma",
                                0);
                        
                        const double gamma_1 = d_equation_of_state->
                            getSpeciesThermodynamicProperty(
                                "gamma",
                                1);
                        
                        // SF6, pre-shock condition.
                        const double rho_SF6 = 5.04;
                        const double u_SF6   = 0.0;
                        const double v_SF6   = 0.0;
                        const double p_SF6   = 1.0/1.4;
                        const double Z_SF6   = 1.0;
                        
                        // air, pre-shock condition.
                        const double rho_pre = 1.0;
                        const double u_pre   = 0.0;
                        const double v_pre   = 0.0;
                        const double p_pre   = 1.0/1.4;
                        const double Z_pre   = 0.0;
                        
                        // air, post-shock condition.
                        const double rho_post = 1.4112;
                        const double u_post   = -0.3613;
                        const double v_post   = 0.0;
                        const double p_post   = 1.6272/1.4;
                        const double Z_post   = 0.0;
                        
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
                                
                                if (x[0] > 3.2*D)
                                {
                                    Z_rho_1[idx_cell] = 0.0;
                                    Z_rho_2[idx_cell] = rho_post;
                                    rho_u[idx_cell]   = rho_post*u_post;
                                    rho_v[idx_cell]   = rho_post*v_post;
                                    E[idx_cell]       = p_post/(gamma_1 - 1.0) +
                                        0.5*rho_post*(u_post*u_post + v_post*v_post);
                                    Z_1[idx_cell]     = Z_post;
                                    Z_2[idx_cell]     = 1.0 - Z_post;
                                }
                                else if (x[0] < (2.9 - 0.1*sin(2*M_PI*(x[1]/D + 0.25)))*D)
                                {
                                    Z_rho_1[idx_cell] = rho_SF6;
                                    Z_rho_2[idx_cell] = 0.0;
                                    rho_u[idx_cell]   = rho_SF6*u_SF6;
                                    rho_v[idx_cell]   = rho_SF6*v_SF6;
                                    E[idx_cell]       = p_SF6/(gamma_0 - 1.0) +
                                        0.5*rho_SF6*(u_SF6*u_SF6 + v_SF6*v_SF6);
                                    Z_1[idx_cell]     = Z_SF6;
                                    Z_2[idx_cell]     = 1.0 - Z_SF6;
                                }
                                else
                                {
                                    Z_rho_1[idx_cell] = 0.0;
                                    Z_rho_2[idx_cell] = rho_pre;
                                    rho_u[idx_cell]   = rho_pre*u_pre;
                                    rho_v[idx_cell]   = rho_pre*v_pre;
                                    E[idx_cell]       = p_pre/(gamma_1 - 1.0) +
                                        0.5*rho_pre*(u_pre*u_pre + v_pre*v_pre);
                                    Z_1[idx_cell]     = Z_pre;
                                    Z_2[idx_cell]     = 1.0 - Z_pre;
                                }
                            }
                        }
                    }
                    */
                    else if (d_project_name == "Advection of 2D material interface")
                    {
                        // Initialize data for a 2D material interface advection problem.
                        
                        if (d_num_species != 2)
                        {
                            TBOX_ERROR(d_object_name
                                       << ": "
                                       << "Please provide only two-species for multi-species simulation testing."
                                       << std::endl);
                        }
                        
                        double* Z_rho_1   = partial_density->getPointer(0);
                        double* Z_rho_2   = partial_density->getPointer(1);
                        double* rho_u     = momentum->getPointer(0);
                        double* rho_v     = momentum->getPointer(1);
                        double* E         = total_energy->getPointer(0);
                        double* Z_1       = volume_fraction->getPointer(0);
                        double* Z_2       = volume_fraction->getPointer(1);
                        
                        const double x_a = 1.0/3*(domain_xlo[0] + domain_xhi[0]);
                        const double x_b = 2.0/3*(domain_xlo[0] + domain_xhi[0]);
                        
                        const double y_a = 1.0/3*(domain_xlo[1] + domain_xhi[1]);
                        const double y_b = 2.0/3*(domain_xlo[1] + domain_xhi[1]);
                        
                        // material initial conditions.
                        double gamma_m  = d_equation_of_state->
                            getSpeciesThermodynamicProperty(
                                "gamma",
                                0);
                        double rho_m    = 10.0;
                        double u_m      = 0.5;
                        double v_m      = 0.5;
                        double p_m      = 1.0/1.4;
                        
                        // ambient initial conditions.
                        double gamma_a  = d_equation_of_state->
                            getSpeciesThermodynamicProperty(
                                "gamma",
                                1);
                        double rho_a    = 1.0;
                        double u_a      = 0.5;
                        double v_a      = 0.5;
                        double p_a      = 1.0/1.4;
                        
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
                                
                                if ((x[0] >= x_a) &&
                                    (x[0] <= x_b) &&
                                    (x[1] >= y_a) &&
                                    (x[1] <= y_b))
                                {
                                    Z_rho_1[idx_cell] = rho_m;
                                    Z_rho_2[idx_cell] = 0.0;
                                    rho_u[idx_cell]   = rho_m*u_m;
                                    rho_v[idx_cell]   = rho_m*v_m;
                                    E[idx_cell]       = p_m/(gamma_m - 1.0) +
                                        0.5*rho_m*(u_m*u_m + v_m*v_m);
                                    Z_1[idx_cell]     = 1.0;
                                    Z_2[idx_cell]     = 0.0;
                                }
                                else
                                {
                                    Z_rho_1[idx_cell] = 0.0;
                                    Z_rho_2[idx_cell] = rho_a;
                                    rho_u[idx_cell]   = rho_a*u_a;
                                    rho_v[idx_cell]   = rho_a*v_a;
                                    E[idx_cell]       = p_a/(gamma_a - 1.0) +
                                        0.5*rho_a*(u_a*u_a + v_a*v_a);
                                    Z_1[idx_cell]     = 0.0;
                                    Z_2[idx_cell]     = 1.0;
                                }
                            }
                        }
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name
                                   << ": "
                                   << "Cannot initialize data for unknown problem."
                                   << std::endl);
                    }
                }
                else if (d_dim == tbox::Dimension(3))
                {
                    if (d_project_name == "Advection of 3D material interface")
                    {
                        // Initialize data for a 3D material interface advection problem.
                        
                        if (d_num_species != 2)
                        {
                            TBOX_ERROR(d_object_name
                                       << ": "
                                       << "Please provide only two-species for multi-species simulation testing."
                                       << std::endl);
                        }
                        
                        double* Z_rho_1   = partial_density->getPointer(0);
                        double* Z_rho_2   = partial_density->getPointer(1);
                        double* rho_u     = momentum->getPointer(0);
                        double* rho_v     = momentum->getPointer(1);
                        double* rho_w     = momentum->getPointer(2);
                        double* E         = total_energy->getPointer(0);
                        double* Z_1       = volume_fraction->getPointer(0);
                        double* Z_2       = volume_fraction->getPointer(1);
                        
                        const double x_a = 1.0/3*(domain_xlo[0] + domain_xhi[0]);
                        const double x_b = 2.0/3*(domain_xlo[0] + domain_xhi[0]);
                        
                        const double y_a = 1.0/3*(domain_xlo[1] + domain_xhi[1]);
                        const double y_b = 2.0/3*(domain_xlo[1] + domain_xhi[1]);
                        
                        const double z_a = 1.0/3*(domain_xlo[2] + domain_xhi[2]);
                        const double z_b = 2.0/3*(domain_xlo[2] + domain_xhi[2]);
                        
                        // material initial conditions.
                        double gamma_m  = d_equation_of_state->
                            getSpeciesThermodynamicProperty(
                                "gamma",
                                0);
                        double rho_m    = 10.0;
                        double u_m      = 0.5;
                        double v_m      = 0.5;
                        double w_m      = 0.5;
                        double p_m      = 1.0/1.4;
                        
                        // ambient initial conditions.
                        double gamma_a  = d_equation_of_state->
                            getSpeciesThermodynamicProperty(
                                "gamma",
                                1);
                        double rho_a    = 1.0;
                        double u_a      = 0.5;
                        double v_a      = 0.5;
                        double w_a      = 0.5;
                        double p_a      = 1.0/1.4;
                        
                        for (int k = 0; k < patch_dims[2]; k++)
                        {
                            for (int j = 0; j < patch_dims[1]; j++)
                            {
                                for (int i = 0; i < patch_dims[0]; i++)
                                {
                                    // Compute index into linear data array.
                                    int idx_cell = i +
                                        j*patch_dims[0] +
                                        k*patch_dims[0]*patch_dims[1];
                                    
                                    // Compute the coordinates.
                                    double x[3];
                                    x[0] = patch_xlo[0] + (i + 0.5)*dx[0];
                                    x[1] = patch_xlo[1] + (j + 0.5)*dx[1];
                                    x[2] = patch_xlo[2] + (k + 0.5)*dx[2];
                                    
                                    if ((x[0] >= x_a) &&
                                        (x[0] <= x_b) &&
                                        (x[1] >= y_a) &&
                                        (x[1] <= y_b) &&
                                        (x[2] >= z_a) &&
                                        (x[2] <= z_b))
                                    {
                                        Z_rho_1[idx_cell] = rho_m;
                                        Z_rho_2[idx_cell] = 0.0;
                                        rho_u[idx_cell]   = rho_m*u_m;
                                        rho_v[idx_cell]   = rho_m*v_m;
                                        rho_w[idx_cell]   = rho_m*w_m;
                                        E[idx_cell]       = p_m/(gamma_m - 1.0) +
                                            0.5*rho_m*(u_m*u_m + v_m*v_m + w_m*w_m);
                                        Z_1[idx_cell]     = 1.0;
                                        Z_2[idx_cell]     = 0.0;
                                    }
                                    else
                                    {
                                        Z_rho_1[idx_cell] = 0.0;
                                        Z_rho_2[idx_cell] = rho_a;
                                        rho_u[idx_cell]   = rho_a*u_a;
                                        rho_v[idx_cell]   = rho_a*v_a;
                                        rho_w[idx_cell]   = rho_a*w_a;
                                        E[idx_cell]       = p_a/(gamma_a - 1.0) +
                                            0.5*rho_a*(u_a*u_a + v_a*v_a + w_a*w_a);
                                        Z_1[idx_cell]     = 0.0;
                                        Z_2[idx_cell]     = 1.0;
                                    }
                                }
                            }
                        }
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name
                                   << ": "
                                   << "Cannot initialize data for unknown problem."
                                   << std::endl);
                    }
                }
                
                break;
            }
        }
    }
    
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
    NULL_USE(initial_time);
    NULL_USE(dt_time);
    
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
    hier::Box dummy_box = patch.getBox();
    const hier::Box interior_box = dummy_box;
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // ghost cells.
    dummy_box.grow(d_num_ghosts);
    const hier::Box ghost_box = dummy_box;
    const hier::IntVector ghostcell_dims = ghost_box.numberCells();
    
    double stable_spectral_radius = 0.0;
    
    switch (d_flow_model)
    {
        case SINGLE_SPECIES:
        {
            boost::shared_ptr<pdat::CellData<double> > density(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_density, getDataContext())));
            
            boost::shared_ptr<pdat::CellData<double> > momentum(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_momentum, getDataContext())));
            
            boost::shared_ptr<pdat::CellData<double> > total_energy(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_total_energy, getDataContext())));
    
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(density);
            TBOX_ASSERT(momentum);
            TBOX_ASSERT(total_energy);
            
            TBOX_ASSERT(density->getGhostCellWidth() == d_num_ghosts);
            TBOX_ASSERT(momentum->getGhostCellWidth() == d_num_ghosts);
            TBOX_ASSERT(total_energy->getGhostCellWidth() == d_num_ghosts);
#endif
            
            if (d_dim == tbox::Dimension(1))
            {
                // Get the pointer of time-dependent variables.
                double* rho   = density->getPointer(0);
                double* rho_u = momentum->getPointer(0);
                double* E     = total_energy->getPointer(0);
                
                for (int i = 0; i < interior_dims[0]; i++)
                {
                    // Compute index of cell into linear data array.
                    const int idx_cell = i + d_num_ghosts[0];
                    
                    const double u = rho_u[idx_cell]/rho[idx_cell];
                    
                    std::vector<const double*> m_ptr;
                    m_ptr.push_back(&(rho_u[idx_cell]));
                    
                    const double c = d_equation_of_state->getSoundSpeed(
                        &(rho[idx_cell]),
                        m_ptr,
                        &(E[idx_cell]));
                    
                    const double spectral_radius = (fabs(u) + c)/dx[0];
                    stable_spectral_radius = fmax(stable_spectral_radius, spectral_radius);
                }
            }
            else if (d_dim == tbox::Dimension(2))
            {
                // Get the pointer of time-dependent variables.
                double* rho   = density->getPointer(0);
                double* rho_u = momentum->getPointer(0);
                double* rho_v = momentum->getPointer(1);
                double* E     = total_energy->getPointer(0);
                
                for (int j = 0; j < interior_dims[1]; j++)
                {
                    for (int i = 0; i < interior_dims[0]; i++)
                    {
                        // Compute index of cell into linear data array.
                        const int idx_cell = (i + d_num_ghosts[0]) +
                            (j + d_num_ghosts[1])*ghostcell_dims[0];
                        
                        const double u = rho_u[idx_cell]/rho[idx_cell];
                        const double v = rho_v[idx_cell]/rho[idx_cell];
                        
                        std::vector<const double*> m_ptr;
                        m_ptr.push_back(&(rho_u[idx_cell]));
                        m_ptr.push_back(&(rho_v[idx_cell]));
                        
                        const double c = d_equation_of_state->getSoundSpeed(
                            &(rho[idx_cell]),
                            m_ptr,
                            &(E[idx_cell]));
                        
                        const double spectral_radius = (fabs(u) + c)/dx[0] + (fabs(v) + c)/dx[1];
                        stable_spectral_radius = fmax(stable_spectral_radius, spectral_radius);
                    }
                }
            }
            else if (d_dim == tbox::Dimension(3))
            {
                // Get the pointer of time-dependent variables.
                double* rho   = density->getPointer(0);
                double* rho_u = momentum->getPointer(0);
                double* rho_v = momentum->getPointer(1);
                double* rho_w = momentum->getPointer(2);
                double* E     = total_energy->getPointer(0);
                
                for (int k = 0; k < interior_dims[2]; k++)
                {
                    for (int j = 0; j < interior_dims[1]; j++)
                    {
                        for (int i = 0; i < interior_dims[0]; i++)
                        {
                            // Compute index of cell into linear data array.
                            const int idx_cell = (i + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            const double u = rho_u[idx_cell]/rho[idx_cell];
                            const double v = rho_v[idx_cell]/rho[idx_cell];
                            const double w = rho_w[idx_cell]/rho[idx_cell];
                            
                            std::vector<const double*> m_ptr;
                            m_ptr.push_back(&(rho_u[idx_cell]));
                            m_ptr.push_back(&(rho_v[idx_cell]));
                            m_ptr.push_back(&(rho_w[idx_cell]));
                            
                            const double c = d_equation_of_state->getSoundSpeed(
                                &(rho[idx_cell]),
                                m_ptr,
                                &(E[idx_cell]));
                            
                            const double spectral_radius = (fabs(u) + c)/dx[0] + (fabs(v) + c)/dx[1] +
                                (fabs(w) + c)/dx[2];
                            
                            stable_spectral_radius = fmax(stable_spectral_radius, spectral_radius);
                        }
                    }
                }
            }
            
            break;
        }
        case FOUR_EQN_SHYUE:
        {
            boost::shared_ptr<pdat::CellData<double> > density(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_density, getDataContext())));
            
            boost::shared_ptr<pdat::CellData<double> > momentum(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_momentum, getDataContext())));
            
            boost::shared_ptr<pdat::CellData<double> > total_energy(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_total_energy, getDataContext())));
            
            boost::shared_ptr<pdat::CellData<double> > mass_fraction(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_mass_fraction, getDataContext())));
            
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(density);
            TBOX_ASSERT(momentum);
            TBOX_ASSERT(total_energy);
            TBOX_ASSERT(mass_fraction);
            
            TBOX_ASSERT(density->getGhostCellWidth() == d_num_ghosts);
            TBOX_ASSERT(momentum->getGhostCellWidth() == d_num_ghosts);
            TBOX_ASSERT(total_energy->getGhostCellWidth() == d_num_ghosts);
            TBOX_ASSERT(mass_fraction->getGhostCellWidth() == d_num_ghosts);
#endif
            
            if (d_dim == tbox::Dimension(1))
            {
                // Get the pointer of time-dependent variables.
                double* rho   = density->getPointer(0);
                double* rho_u = momentum->getPointer(0);
                double* E     = total_energy->getPointer(0);
                std::vector<double*> Y;
                for (int si = 0; si < d_num_species; si++)
                {
                    Y.push_back(mass_fraction->getPointer(si));
                }
                
                for (int i = 0; i < interior_dims[0]; i++)
                {
                    // Compute index of cell into linear data array.
                    const int idx_cell = i + d_num_ghosts[0];
                    
                    const double u = rho_u[idx_cell]/rho[idx_cell];
                    
                    std::vector<const double*> m_ptr;
                    m_ptr.push_back(&(rho_u[idx_cell]));
                    
                    std::vector<const double*> Y_ptr;
                    for (int si = 0; si < d_num_species; si++)
                    {
                        Y_ptr.push_back(&(Y[si][idx_cell]));
                    }
                    
                    const double c = d_equation_of_state->getSoundSpeedWithMassFraction(
                        &(rho[idx_cell]),
                        m_ptr,
                        &(E[idx_cell]),
                        Y_ptr);
                    
                    const double spectral_radius = (fabs(u) + c)/dx[0];
                    
                    stable_spectral_radius = fmax(stable_spectral_radius, spectral_radius);
                }
            }
            else if (d_dim == tbox::Dimension(2))
            {
                // Get the pointer of time-dependent variables.
                double* rho   = density->getPointer(0);
                double* rho_u = momentum->getPointer(0);
                double* rho_v = momentum->getPointer(1);
                double* E     = total_energy->getPointer(0);
                std::vector<double*> Y;
                for (int si = 0; si < d_num_species; si++)
                {
                    Y.push_back(mass_fraction->getPointer(si));
                }
                
                for (int j = 0; j < interior_dims[1]; j++)
                {
                    for (int i = 0; i < interior_dims[0]; i++)
                    {
                        // Compute index of cell into linear data array.
                        const int idx_cell = (i + d_num_ghosts[0]) +
                            (j + d_num_ghosts[1])*ghostcell_dims[0];
                        
                        const double u = rho_u[idx_cell]/rho[idx_cell];
                        const double v = rho_v[idx_cell]/rho[idx_cell];
                        
                        std::vector<const double*> m_ptr;
                        m_ptr.push_back(&(rho_u[idx_cell]));
                        m_ptr.push_back(&(rho_v[idx_cell]));
                        
                        std::vector<const double*> Y_ptr;
                        for (int si = 0; si < d_num_species; si++)
                        {
                            Y_ptr.push_back(&(Y[si][idx_cell]));
                        }
                        
                        const double c = d_equation_of_state->getSoundSpeedWithMassFraction(
                            &(rho[idx_cell]),
                            m_ptr,
                            &(E[idx_cell]),
                            Y_ptr);
                        
                        const double spectral_radius = (fabs(u) + c)/dx[0] + (fabs(v) + c)/dx[1];
                        
                        stable_spectral_radius = fmax(stable_spectral_radius, spectral_radius);
                    }
                }
            }
            else if (d_dim == tbox::Dimension(3))
            {
                // Get the pointer of time-dependent variables.
                double* rho   = density->getPointer(0);
                double* rho_u = momentum->getPointer(0);
                double* rho_v = momentum->getPointer(1);
                double* rho_w = momentum->getPointer(2);
                double* E     = total_energy->getPointer(0);
                std::vector<double*> Y;
                for (int si = 0; si < d_num_species; si++)
                {
                    Y.push_back(mass_fraction->getPointer(si));
                }
                
                for (int k = 0; k < interior_dims[2]; k++)
                {
                    for (int j = 0; j < interior_dims[1]; j++)
                    {
                        for (int i = 0; i < interior_dims[0]; i++)
                        {
                            // Compute index of cell into linear data array.
                            const int idx_cell = (i + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            const double u = rho_u[idx_cell]/rho[idx_cell];
                            const double v = rho_v[idx_cell]/rho[idx_cell];
                            const double w = rho_w[idx_cell]/rho[idx_cell];
                            
                            std::vector<const double*> m_ptr;
                            m_ptr.push_back(&(rho_u[idx_cell]));
                            m_ptr.push_back(&(rho_v[idx_cell]));
                            m_ptr.push_back(&(rho_w[idx_cell]));
                            
                            std::vector<const double*> Y_ptr;
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_ptr.push_back(&(Y[si][idx_cell]));
                            }
                            
                            const double c = d_equation_of_state->getSoundSpeedWithMassFraction(
                                &(rho[idx_cell]),
                                m_ptr,
                                &(E[idx_cell]),
                                Y_ptr);
                            
                            const double spectral_radius = (fabs(u) + c)/dx[0] + (fabs(v) + c)/dx[1] +
                                (fabs(w) + c)/dx[2];
                            
                            stable_spectral_radius = fmax(stable_spectral_radius, spectral_radius);
                        }
                    }
                }
            }
            
            break;
        }
        case FIVE_EQN_ALLAIRE:
        {
            boost::shared_ptr<pdat::CellData<double> > partial_density(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_partial_density, getDataContext())));
            
            boost::shared_ptr<pdat::CellData<double> > momentum(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_momentum, getDataContext())));
            
            boost::shared_ptr<pdat::CellData<double> > total_energy(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_total_energy, getDataContext())));
            
            boost::shared_ptr<pdat::CellData<double> > volume_fraction(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_volume_fraction, getDataContext())));
            
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(partial_density);
            TBOX_ASSERT(momentum);
            TBOX_ASSERT(total_energy);
            TBOX_ASSERT(volume_fraction);
            
            TBOX_ASSERT(partial_density->getGhostCellWidth() == d_num_ghosts);
            TBOX_ASSERT(momentum->getGhostCellWidth() == d_num_ghosts);
            TBOX_ASSERT(total_energy->getGhostCellWidth() == d_num_ghosts);
            TBOX_ASSERT(volume_fraction->getGhostCellWidth() == d_num_ghosts);
#endif
            
            if (d_dim == tbox::Dimension(1))
            {
                // Get the pointer of time-dependent variables.
                std::vector<double*> Z_rho;
                for (int si = 0; si < d_num_species; si++)
                {
                    Z_rho.push_back(partial_density->getPointer(si));
                }
                double* rho_u = momentum->getPointer(0);
                double* E     = total_energy->getPointer(0);
                std::vector<double*> Z;
                for (int si = 0; si < d_num_species; si++)
                {
                    Z.push_back(volume_fraction->getPointer(si));
                }
                
                for (int i = 0; i < interior_dims[0]; i++)
                {
                    // Compute index of cell into linear data array.
                    const int idx_cell = (i + d_num_ghosts[0]);
                    
                    std::vector<const double*> partial_density_idx_cell;
                    for (int si = 0; si < d_num_species; si++)
                    {
                        partial_density_idx_cell.push_back(&(Z_rho[si][idx_cell]));
                    }
                    
                    const double rho = d_equation_of_state->getTotalDensity(
                        partial_density_idx_cell);
                    
                    const double u = rho_u[idx_cell]/rho;
                    
                    std::vector<const double*> m_ptr;
                    m_ptr.push_back(&(rho_u[idx_cell]));
                    
                    std::vector<const double*> Z_ptr;
                    for (int si = 0; si < d_num_species; si++)
                    {
                        Z_ptr.push_back(&(Z[si][idx_cell]));
                    }
                    
                    const double c = d_equation_of_state->getSoundSpeedWithVolumeFraction(
                        &rho,
                        m_ptr,
                        &(E[idx_cell]),
                        Z_ptr);
                    
                    const double spectral_radius = (fabs(u) + c)/dx[0];
                    
                    stable_spectral_radius = fmax(stable_spectral_radius, spectral_radius);
                }
            }
            else if (d_dim == tbox::Dimension(2))
            {
                // Get the pointer of time-dependent variables.
                std::vector<double*> Z_rho;
                for (int si = 0; si < d_num_species; si++)
                {
                    Z_rho.push_back(partial_density->getPointer(si));
                }
                double* rho_u = momentum->getPointer(0);
                double* rho_v = momentum->getPointer(1);
                double* E     = total_energy->getPointer(0);
                std::vector<double*> Z;
                for (int si = 0; si < d_num_species; si++)
                {
                    Z.push_back(volume_fraction->getPointer(si));
                }
                
                for (int j = 0; j < interior_dims[1]; j++)
                {
                    for (int i = 0; i < interior_dims[0]; i++)
                    {
                        // Compute index of cell into linear data array.
                        const int idx_cell = (i + d_num_ghosts[0]) +
                            (j + d_num_ghosts[1])*ghostcell_dims[0];
                        
                        std::vector<const double*> partial_density_idx_cell;
                        for (int si = 0; si < d_num_species; si++)
                        {
                            partial_density_idx_cell.push_back(&(Z_rho[si][idx_cell]));
                        }
                        
                        const double rho = d_equation_of_state->getTotalDensity(
                            partial_density_idx_cell);
                        
                        const double u = rho_u[idx_cell]/rho;
                        const double v = rho_v[idx_cell]/rho;
                        
                        std::vector<const double*> m_ptr;
                        m_ptr.push_back(&(rho_u[idx_cell]));
                        m_ptr.push_back(&(rho_v[idx_cell]));
                        
                        std::vector<const double*> Z_ptr;
                        for (int si = 0; si < d_num_species; si++)
                        {
                            Z_ptr.push_back(&(Z[si][idx_cell]));
                        }
                        
                        const double c = d_equation_of_state->getSoundSpeedWithVolumeFraction(
                            &rho,
                            m_ptr,
                            &(E[idx_cell]),
                            Z_ptr);
                        
                        const double spectral_radius = (fabs(u) + c)/dx[0] + (fabs(v) + c)/dx[1];
                        
                        stable_spectral_radius = fmax(stable_spectral_radius, spectral_radius);
                    }
                }
            }
            else if (d_dim == tbox::Dimension(3))
            {
                // Get the pointer of time-dependent variables.
                std::vector<double*> Z_rho;
                for (int si = 0; si < d_num_species; si++)
                {
                    Z_rho.push_back(partial_density->getPointer(si));
                }
                double* rho_u = momentum->getPointer(0);
                double* rho_v = momentum->getPointer(1);
                double* rho_w = momentum->getPointer(2);
                double* E     = total_energy->getPointer(0);
                std::vector<double*> Z;
                for (int si = 0; si < d_num_species; si++)
                {
                    Z.push_back(volume_fraction->getPointer(si));
                }
                
                for (int k = 0; k < interior_dims[2]; k++)
                {
                    for (int j = 0; j < interior_dims[1]; j++)
                    {
                        for (int i = 0; i < interior_dims[0]; i++)
                        {
                            // Compute index of cell into linear data array.
                            const int idx_cell = (i + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            std::vector<const double*> partial_density_idx_cell;
                            for (int si = 0; si < d_num_species; si++)
                            {
                                partial_density_idx_cell.push_back(&(Z_rho[si][idx_cell]));
                            }
                            
                            const double rho = d_equation_of_state->getTotalDensity(
                                partial_density_idx_cell);
                            
                            const double u = rho_u[idx_cell]/rho;
                            const double v = rho_v[idx_cell]/rho;
                            const double w = rho_w[idx_cell]/rho;
                            
                            std::vector<const double*> m_ptr;
                            m_ptr.push_back(&(rho_u[idx_cell]));
                            m_ptr.push_back(&(rho_v[idx_cell]));
                            m_ptr.push_back(&(rho_w[idx_cell]));
                            
                            std::vector<const double*> Z_ptr;
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Z_ptr.push_back(&(Z[si][idx_cell]));
                            }
                            
                            const double c = d_equation_of_state->getSoundSpeedWithVolumeFraction(
                                &rho,
                                m_ptr,
                                &(E[idx_cell]),
                                Z_ptr);
                            
                            const double spectral_radius = (fabs(u) + c)/dx[0] + (fabs(v) + c)/dx[1] +
                                (fabs(w) + c)/dx[2];
                            
                            stable_spectral_radius = fmax(stable_spectral_radius, spectral_radius);
                        }
                    }
                }
            }
            
            break;
        }
        default:
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "Unknown d_flow_model."
                       << std::endl);
        }
    }
    
    stable_dt = 1.0/stable_spectral_radius;
    
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
    
    std::vector<double*> Q;
    
    switch (d_flow_model)
    {
        case SINGLE_SPECIES:
        {
            boost::shared_ptr<pdat::CellData<double> > density(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_density, getDataContext())));
            
            boost::shared_ptr<pdat::CellData<double> > momentum(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_momentum, getDataContext())));
            
            boost::shared_ptr<pdat::CellData<double> > total_energy(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_total_energy, getDataContext())));
            
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(density);
            TBOX_ASSERT(momentum);
            TBOX_ASSERT(total_energy);
            
            TBOX_ASSERT(density->getGhostCellWidth() == d_num_ghosts);
            TBOX_ASSERT(momentum->getGhostCellWidth() == d_num_ghosts);
            TBOX_ASSERT(total_energy->getGhostCellWidth() == d_num_ghosts);
#endif
            
            // Initialize all time-dependent data within the interior box with zero values
            density->fillAll(0.0, interior_box);
            momentum->fillAll(0.0, interior_box);
            total_energy->fillAll(0.0, interior_box);
            
            Q.push_back(density->getPointer(0));
            for (int di = 0; di < d_dim.getValue(); di++)
            {
                Q.push_back(momentum->getPointer(di));
            }
            Q.push_back(total_energy->getPointer(0));
            
            break;
        }
        case FOUR_EQN_SHYUE:
        {    
            boost::shared_ptr<pdat::CellData<double> > density(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_density, getDataContext())));
            
            boost::shared_ptr<pdat::CellData<double> > momentum(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_momentum, getDataContext())));
            
            boost::shared_ptr<pdat::CellData<double> > total_energy(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_total_energy, getDataContext())));
            
            boost::shared_ptr<pdat::CellData<double> > mass_fraction(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_mass_fraction, getDataContext())));
            
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(density);
            TBOX_ASSERT(momentum);
            TBOX_ASSERT(total_energy);
            TBOX_ASSERT(mass_fraction);
            
            TBOX_ASSERT(density->getGhostCellWidth() == d_num_ghosts);
            TBOX_ASSERT(momentum->getGhostCellWidth() == d_num_ghosts);
            TBOX_ASSERT(total_energy->getGhostCellWidth() == d_num_ghosts);
            TBOX_ASSERT(mass_fraction->getGhostCellWidth() == d_num_ghosts);
#endif
            
            // Initialize all time-dependent data within the interior box with zero values
            density->fillAll(0.0, interior_box);
            momentum->fillAll(0.0, interior_box);
            total_energy->fillAll(0.0, interior_box);
            mass_fraction->fillAll(0.0, interior_box);
            
            Q.push_back(density->getPointer(0));
            for (int di = 0; di < d_dim.getValue(); di++)
            {
                Q.push_back(momentum->getPointer(di));
            }
            Q.push_back(total_energy->getPointer(0));
            for (int si = 0; si < d_num_species; si++)
            {
                Q.push_back(mass_fraction->getPointer(si));
            }
            
            break;
        }
        case FIVE_EQN_ALLAIRE:
        {
            boost::shared_ptr<pdat::CellData<double> > partial_density(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_partial_density, getDataContext())));
            
            boost::shared_ptr<pdat::CellData<double> > momentum(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_momentum, getDataContext())));
            
            boost::shared_ptr<pdat::CellData<double> > total_energy(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_total_energy, getDataContext())));
            
            boost::shared_ptr<pdat::CellData<double> > volume_fraction(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_volume_fraction, getDataContext())));
                        
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(partial_density);
            TBOX_ASSERT(momentum);
            TBOX_ASSERT(total_energy);
            TBOX_ASSERT(volume_fraction);
            
            TBOX_ASSERT(partial_density->getGhostCellWidth() == d_num_ghosts);
            TBOX_ASSERT(momentum->getGhostCellWidth() == d_num_ghosts);
            TBOX_ASSERT(total_energy->getGhostCellWidth() == d_num_ghosts);
            TBOX_ASSERT(volume_fraction->getGhostCellWidth() == d_num_ghosts);
#endif
            
            // Initialize all time-dependent data within the interior box with zero values
            partial_density->fillAll(0.0, interior_box);
            momentum->fillAll(0.0, interior_box);
            total_energy->fillAll(0.0, interior_box);
            volume_fraction->fillAll(0.0, interior_box);
            
            for (int si = 0; si < d_num_species; si++)
            {
                Q.push_back(partial_density->getPointer(si));
            }
            for (int di = 0; di < d_dim.getValue(); di++)
            {
                Q.push_back(momentum->getPointer(di));
            }
            Q.push_back(total_energy->getPointer(0));
            for (int si = 0; si < d_num_species; si++)
            {
                Q.push_back(volume_fraction->getPointer(si));
            }
            
            break;
        }
        default:
            TBOX_ERROR(d_object_name
                       << ": "
                       << "Unknown d_flow_model."
                       << std::endl);
    }
    
    /*
     * Use alpah, beta and gamma values to update the time-dependent solution,
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
        
        std::vector<double*> Q_intermediate;
        
        switch (d_flow_model)
        {
            case SINGLE_SPECIES:
            {
                boost::shared_ptr<pdat::CellData<double> > density_intermediate(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_density, intermediate_context[n])));
                
                boost::shared_ptr<pdat::CellData<double> > momentum_intermediate(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_momentum, intermediate_context[n])));
                
                boost::shared_ptr<pdat::CellData<double> > total_energy_intermediate(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_total_energy, intermediate_context[n])));
                
#ifdef DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(density_intermediate);
                TBOX_ASSERT(momentum_intermediate);
                TBOX_ASSERT(total_energy_intermediate);
                
                TBOX_ASSERT(density_intermediate->getGhostCellWidth() == d_num_ghosts);
                TBOX_ASSERT(momentum_intermediate->getGhostCellWidth() == d_num_ghosts);
                TBOX_ASSERT(total_energy_intermediate->getGhostCellWidth() == d_num_ghosts);
#endif
                
                Q_intermediate.push_back(density_intermediate->getPointer(0));
                for (int di = 0; di < d_dim.getValue(); di++)
                {
                    Q_intermediate.push_back(momentum_intermediate->getPointer(di));
                }
                Q_intermediate.push_back(total_energy_intermediate->getPointer(0));
                
                break;
            }
            case FOUR_EQN_SHYUE:
            {
                boost::shared_ptr<pdat::CellData<double> > density_intermediate(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_density, intermediate_context[n])));
                
                boost::shared_ptr<pdat::CellData<double> > momentum_intermediate(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_momentum, intermediate_context[n])));
                
                boost::shared_ptr<pdat::CellData<double> > total_energy_intermediate(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_total_energy, intermediate_context[n])));
                
                boost::shared_ptr<pdat::CellData<double> > mass_fraction_intermediate(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_mass_fraction, intermediate_context[n])));
                
#ifdef DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(density_intermediate);
                TBOX_ASSERT(momentum_intermediate);
                TBOX_ASSERT(total_energy_intermediate);
                TBOX_ASSERT(mass_fraction_intermediate);
                
                TBOX_ASSERT(density_intermediate->getGhostCellWidth() == d_num_ghosts);
                TBOX_ASSERT(momentum_intermediate->getGhostCellWidth() == d_num_ghosts);
                TBOX_ASSERT(total_energy_intermediate->getGhostCellWidth() == d_num_ghosts);
                TBOX_ASSERT(mass_fraction_intermediate->getGhostCellWidth() == d_num_ghosts);
#endif
                
                Q_intermediate.push_back(density_intermediate->getPointer(0));
                for (int di = 0; di < d_dim.getValue(); di++)
                {
                    Q_intermediate.push_back(momentum_intermediate->getPointer(di));
                }
                Q_intermediate.push_back(total_energy_intermediate->getPointer(0));
                for (int si = 0; si < d_num_species; si++)
                {
                    Q_intermediate.push_back(mass_fraction_intermediate->getPointer(si));
                }
                
                break;
            }
            case FIVE_EQN_ALLAIRE:
            {
                boost::shared_ptr<pdat::CellData<double> > partial_density_intermediate(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_partial_density, intermediate_context[n])));
                
                boost::shared_ptr<pdat::CellData<double> > momentum_intermediate(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_momentum, intermediate_context[n])));
                
                boost::shared_ptr<pdat::CellData<double> > total_energy_intermediate(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_total_energy, intermediate_context[n])));
                
                boost::shared_ptr<pdat::CellData<double> > volume_fraction_intermediate(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_volume_fraction, intermediate_context[n])));
                
#ifdef DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(partial_density_intermediate);
                TBOX_ASSERT(momentum_intermediate);
                TBOX_ASSERT(total_energy_intermediate);
                TBOX_ASSERT(volume_fraction_intermediate);
                
                TBOX_ASSERT(partial_density_intermediate->getGhostCellWidth() == d_num_ghosts);
                TBOX_ASSERT(momentum_intermediate->getGhostCellWidth() == d_num_ghosts);
                TBOX_ASSERT(total_energy_intermediate->getGhostCellWidth() == d_num_ghosts);
                TBOX_ASSERT(volume_fraction_intermediate->getGhostCellWidth() == d_num_ghosts);
#endif
                
                for (int si = 0; si < d_num_species; si++)
                {
                    Q_intermediate.push_back(partial_density_intermediate->getPointer(si));
                }
                for (int di = 0; di < d_dim.getValue(); di++)
                {
                    Q_intermediate.push_back(momentum_intermediate->getPointer(di));
                }
                Q_intermediate.push_back(total_energy_intermediate->getPointer(0));
                for (int si = 0; si < d_num_species; si++)
                {
                    Q_intermediate.push_back(volume_fraction_intermediate->getPointer(si));
                }
                
                break;
            }
            default:
            TBOX_ERROR(d_object_name
                       << ": "
                       << "Unknown d_flow_model."
                       << std::endl);
        }
        
        if (d_dim == tbox::Dimension(1))
        {
            if (alpha[n] != 0.0)
            {
                for (int i = 0; i < interior_dims[0]; i++)
                {
                    // Compute indices of time-dependent data.
                    int idx_cell   = i + d_num_ghosts[0];
                    
                    for (int ei = 0; ei < d_num_eqn; ei++)
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
                    
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        double* F_x_intermediate = convective_flux_intermediate->getPointer(0, ei);
                        double* S_intermediate = source_intermediate->getPointer(ei);
                        
                        Q[ei][idx_cell] += beta[n]*
                            (-(F_x_intermediate[idx_flux_x] - F_x_intermediate[idx_flux_x - 1])/dx[0] +
                            S_intermediate[idx_source]);
                    }
                    
                    // Update the mass fraction/volume fraction of the last species.
                    switch (d_flow_model)
                    {
                        case SINGLE_SPECIES:
                            break;
                        case FOUR_EQN_SHYUE:
                        {
                            Q[d_num_eqn][idx_cell] = 1.0;
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                Q[d_num_eqn][idx_cell] -= Q[d_num_eqn - 1 - si][idx_cell];
                            }
                            break;
                        }
                        case FIVE_EQN_ALLAIRE:
                        {
                            Q[d_num_eqn][idx_cell] = 1.0;
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                Q[d_num_eqn][idx_cell] -= Q[d_num_eqn - 1 - si][idx_cell];
                            }
                            break;
                        }
                        default:
                            TBOX_ERROR(d_object_name
                                       << ": "
                                       << "Unknown d_flow_model."
                                       << std::endl);
                    }
                }
            }
                
            if (gamma[n] != 0.0)
            {
                // Accumulate the flux in the x direction.
                for (int i = 0; i < interior_dims[0] + 1; i++)
                {
                    int idx_flux_x = i;
                    
                    for (int ei = 0; ei < d_num_eqn; ei++)
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
                    
                    for (int ei = 0; ei < d_num_eqn; ei++)
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
                        
                        for (int ei = 0; ei < d_num_eqn; ei++)
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
                        
                        for (int ei = 0; ei < d_num_eqn; ei++)
                        {
                            double* F_x_intermediate = convective_flux_intermediate->getPointer(0, ei);
                            double* F_y_intermediate = convective_flux_intermediate->getPointer(1, ei);
                            double* S_intermediate = source_intermediate->getPointer(ei);
                            
                            Q[ei][idx_cell] += beta[n]*
                                (-(F_x_intermediate[idx_flux_x] - F_x_intermediate[idx_flux_x - 1])/dx[0] -
                                (F_y_intermediate[idx_flux_y] - F_y_intermediate[idx_flux_y - 1])/dx[1] +
                                S_intermediate[idx_source]);
                        }
                        
                        // Update the mass fraction/volume fraction of the last species.
                        switch (d_flow_model)
                        {
                            case SINGLE_SPECIES:
                                break;
                            case FOUR_EQN_SHYUE:
                            {
                                Q[d_num_eqn][idx_cell] = 1.0;
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    Q[d_num_eqn][idx_cell] -= Q[d_num_eqn - 1 - si][idx_cell];
                                }
                                break;
                            }
                            case FIVE_EQN_ALLAIRE:
                            {
                                Q[d_num_eqn][idx_cell] = 1.0;
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    Q[d_num_eqn][idx_cell] -= Q[d_num_eqn - 1 - si][idx_cell];
                                }
                                break;
                            }
                            default:
                                TBOX_ERROR(d_object_name
                                           << ": "
                                           << "Unknown d_flow_model."
                                           << std::endl);
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
                        
                        for (int ei = 0; ei < d_num_eqn; ei++)
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
                        
                        for (int ei = 0; ei < d_num_eqn; ei++)
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
                        
                        for (int ei = 0; ei < d_num_eqn; ei++)
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
                            
                            for (int ei = 0; ei < d_num_eqn; ei++)
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
                            
                            for (int ei = 0; ei < d_num_eqn; ei++)
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
                            
                            // Update the mass fraction/volume fraction of the last species.
                            switch (d_flow_model)
                            {
                                case SINGLE_SPECIES:
                                    break;
                                case FOUR_EQN_SHYUE:
                                {
                                    Q[d_num_eqn][idx_cell] = 1.0;
                                    for (int si = 0; si < d_num_species - 1; si++)
                                    {
                                        Q[d_num_eqn][idx_cell] -= Q[d_num_eqn - 1 - si][idx_cell];
                                    }
                                    break;
                                }
                                case FIVE_EQN_ALLAIRE:
                                {
                                    Q[d_num_eqn][idx_cell] = 1.0;
                                    for (int si = 0; si < d_num_species - 1; si++)
                                    {
                                        Q[d_num_eqn][idx_cell] -= Q[d_num_eqn - 1 - si][idx_cell];
                                    }
                                    break;
                                }
                                default:
                                    TBOX_ERROR(d_object_name
                                               << ": "
                                               << "Unknown d_flow_model."
                                               << std::endl);
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
                            
                            for (int ei = 0; ei < d_num_eqn; ei++)
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
                            
                            for (int ei = 0; ei < d_num_eqn; ei++)
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
                            
                            for (int ei = 0; ei < d_num_eqn; ei++)
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
                            
                            for (int ei = 0; ei < d_num_eqn; ei++)
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
    
    std::vector<double*> Q;
    
    switch (d_flow_model)
    {
        case SINGLE_SPECIES:
        {
            boost::shared_ptr<pdat::CellData<double> > density(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_density, getDataContext())));
            
            boost::shared_ptr<pdat::CellData<double> > momentum(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_momentum, getDataContext())));
            
            boost::shared_ptr<pdat::CellData<double> > total_energy(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_total_energy, getDataContext())));
            
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(density);
            TBOX_ASSERT(momentum);
            TBOX_ASSERT(total_energy);
            
            TBOX_ASSERT(density->getGhostCellWidth() == d_num_ghosts);
            TBOX_ASSERT(momentum->getGhostCellWidth() == d_num_ghosts);
            TBOX_ASSERT(total_energy->getGhostCellWidth() == d_num_ghosts);
#endif
            Q.push_back(density->getPointer(0));
            for (int di = 0; di < d_dim.getValue(); di++)
            {
                Q.push_back(momentum->getPointer(di));
            }
            Q.push_back(total_energy->getPointer(0));
            
            break;
        }
        case FOUR_EQN_SHYUE:
        {    
            boost::shared_ptr<pdat::CellData<double> > density(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_density, getDataContext())));
            
            boost::shared_ptr<pdat::CellData<double> > momentum(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_momentum, getDataContext())));
            
            boost::shared_ptr<pdat::CellData<double> > total_energy(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_total_energy, getDataContext())));
            
            boost::shared_ptr<pdat::CellData<double> > mass_fraction(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_mass_fraction, getDataContext())));
            
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(density);
            TBOX_ASSERT(momentum);
            TBOX_ASSERT(total_energy);
            TBOX_ASSERT(mass_fraction);
            
            TBOX_ASSERT(density->getGhostCellWidth() == d_num_ghosts);
            TBOX_ASSERT(momentum->getGhostCellWidth() == d_num_ghosts);
            TBOX_ASSERT(total_energy->getGhostCellWidth() == d_num_ghosts);
            TBOX_ASSERT(mass_fraction->getGhostCellWidth() == d_num_ghosts);
#endif
            
            Q.push_back(density->getPointer(0));
            for (int di = 0; di < d_dim.getValue(); di++)
            {
                Q.push_back(momentum->getPointer(di));
            }
            Q.push_back(total_energy->getPointer(0));
            for (int si = 0; si < d_num_species; si++)
            {
                Q.push_back(mass_fraction->getPointer(si));
            }
            
            break;
        }
        case FIVE_EQN_ALLAIRE:
        {
            boost::shared_ptr<pdat::CellData<double> > partial_density(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_partial_density, getDataContext())));
            
            boost::shared_ptr<pdat::CellData<double> > momentum(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_momentum, getDataContext())));
            
            boost::shared_ptr<pdat::CellData<double> > total_energy(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_total_energy, getDataContext())));
            
            boost::shared_ptr<pdat::CellData<double> > volume_fraction(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_volume_fraction, getDataContext())));
                        
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(partial_density);
            TBOX_ASSERT(momentum);
            TBOX_ASSERT(total_energy);
            TBOX_ASSERT(volume_fraction);
            
            TBOX_ASSERT(partial_density->getGhostCellWidth() == d_num_ghosts);
            TBOX_ASSERT(momentum->getGhostCellWidth() == d_num_ghosts);
            TBOX_ASSERT(total_energy->getGhostCellWidth() == d_num_ghosts);
            TBOX_ASSERT(volume_fraction->getGhostCellWidth() == d_num_ghosts);
#endif
            for (int si = 0; si < d_num_species; si++)
            {
                Q.push_back(partial_density->getPointer(si));
            }
            for (int di = 0; di < d_dim.getValue(); di++)
            {
                Q.push_back(momentum->getPointer(di));
            }
            Q.push_back(total_energy->getPointer(0));
            for (int si = 0; si < d_num_species; si++)
            {
                Q.push_back(volume_fraction->getPointer(si));
            }
            
            break;
        }
        default:
            TBOX_ERROR(d_object_name
                       << ": "
                       << "Unknown d_flow_model."
                       << std::endl);
    }
    
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
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                double *F_x = convective_flux->getPointer(0, ei);
                double *S   = source->getPointer(ei);
                
                Q[ei][idx_cell] +=
                    (-(F_x[idx_flux_x] - F_x[idx_flux_x - 1])/dx[0] +
                    S[idx_source]);
                
                // Update the mass fraction/volume fraction of the last species.
                switch (d_flow_model)
                {
                    case SINGLE_SPECIES:
                        break;
                    case FOUR_EQN_SHYUE:
                        Q[d_num_eqn][idx_cell] = 1.0;
                        for (int si = 0; si < d_num_species - 1; si++)
                        {
                            Q[d_num_eqn][idx_cell] -= Q[d_num_eqn - 1 - si][idx_cell];
                        }
                        break;
                    case FIVE_EQN_ALLAIRE:
                        Q[d_num_eqn][idx_cell] = 1.0;
                        for (int si = 0; si < d_num_species - 1; si++)
                        {
                            Q[d_num_eqn][idx_cell] -= Q[d_num_eqn - 1 - si][idx_cell];
                        }
                        break;
                    default:
                        TBOX_ERROR(d_object_name
                                   << ": "
                                   << "Unknown d_flow_model."
                                   << std::endl);
                }
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
                
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    double *F_x = convective_flux->getPointer(0, ei);
                    double *F_y = convective_flux->getPointer(1, ei);
                    double *S   = source->getPointer(ei);
                    
                    Q[ei][idx_cell] +=
                        (-(F_x[idx_flux_x] - F_x[idx_flux_x - 1])/dx[0] -
                        (F_y[idx_flux_y] - F_y[idx_flux_y - 1])/dx[1] +
                        S[idx_source]);
                    
                    // Update the mass fraction/volume fraction of the last species.
                    switch (d_flow_model)
                    {
                        case SINGLE_SPECIES:
                            break;
                        case FOUR_EQN_SHYUE:
                            Q[d_num_eqn][idx_cell] = 1.0;
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                Q[d_num_eqn][idx_cell] -= Q[d_num_eqn - 1 - si][idx_cell];
                            }
                            break;
                        case FIVE_EQN_ALLAIRE:
                            Q[d_num_eqn][idx_cell] = 1.0;
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                Q[d_num_eqn][idx_cell] -= Q[d_num_eqn - 1 - si][idx_cell];
                            }
                            break;
                        default:
                            TBOX_ERROR(d_object_name
                                       << ": "
                                       << "Unknown d_flow_model."
                                       << std::endl);
                    }
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
                    
                    for (int ei = 0; ei < d_num_eqn; ei++)
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
                    
                    // Update the mass fraction/volume fraction of the last species.
                    switch (d_flow_model)
                    {
                        case SINGLE_SPECIES:
                            break;
                        case FOUR_EQN_SHYUE:
                            Q[d_num_eqn][idx_cell] = 1.0;
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                Q[d_num_eqn][idx_cell] -= Q[d_num_eqn - 1 - si][idx_cell];
                            }
                            break;
                        case FIVE_EQN_ALLAIRE:
                            Q[d_num_eqn][idx_cell] = 1.0;
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                Q[d_num_eqn][idx_cell] -= Q[d_num_eqn - 1 - si][idx_cell];
                            }
                            break;
                        default:
                            TBOX_ERROR(d_object_name
                                       << ": "
                                       << "Unknown d_flow_model."
                                       << std::endl);
                    }
                }
            }
        }
    }

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
    
    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(patch_geom);
#endif
    
    const double* dx = patch_geom->getDx();
    
    NULL_USE(dx);
    
    boost::shared_ptr<pdat::CellData<int> > tags(
        BOOST_CAST<pdat::CellData<int>, hier::PatchData>(
            patch.getPatchData(tag_indx)));

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(tags);
    TBOX_ASSERT(tags->getGhostCellWidth() == 0);
#endif    
    
    boost::shared_ptr<pdat::CellData<double> > density(
        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
            patch.getPatchData(d_density, getDataContext())));
    
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(density);    
    TBOX_ASSERT(density->getGhostCellWidth() == d_num_ghosts);
#endif
    
    // Get the dimensions of box that covers the interior of Patch.
    hier::Box dummy_box = patch.getBox();
    const hier::Box interior_box = dummy_box;
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of Patch plus
    // ghost cells.
    dummy_box.grow(d_num_ghosts);
    const hier::Box ghost_box = dummy_box;
    const hier::IntVector ghostcell_dims = ghost_box.numberCells();
    
    if (d_dim == tbox::Dimension(1))
    {
        // NOT YET IMPLEMENTED
    }
    else if (d_dim == tbox::Dimension(2))
    {
        for (int ncrit = 0;
             ncrit < static_cast<int>(d_refinement_criteria.size());
             ncrit++)
        {
            std::string ref = d_refinement_criteria[ncrit];
            
            // Get the pointer of the tags
            int* tag_ptr  = tags->getPointer();
            
            if (ref == "DENSITY_SHOCK")
            {
                // Get the pointer of conservative variables.
                double* rho   = density->getPointer(0);
                
                for (int j = 0; j < interior_dims[1]; j++)
                {
                    for (int i = 0; i < interior_dims[0]; i++)
                    {
                        // Compute index into linear data array.
                        // Note: the data is stored in Fortran order.
                        int idx_wghost_x[3];
                        int idx_wghost_y[3];
                        
                        idx_wghost_x[0] = ((i - 1) + d_num_ghosts[0]) + (j + d_num_ghosts[1])*ghostcell_dims[0];
                        idx_wghost_x[1] = ( i      + d_num_ghosts[0]) + (j + d_num_ghosts[1])*ghostcell_dims[0];
                        idx_wghost_x[2] = ((i + 1) + d_num_ghosts[0]) + (j + d_num_ghosts[1])*ghostcell_dims[0];
                        
                        idx_wghost_y[0] = (i + d_num_ghosts[0]) + ((j - 1) + d_num_ghosts[1])*ghostcell_dims[0];
                        idx_wghost_y[1] = (i + d_num_ghosts[0]) + ( j      + d_num_ghosts[1])*ghostcell_dims[0];
                        idx_wghost_y[2] = (i + d_num_ghosts[0]) + ((j + 1) + d_num_ghosts[1])*ghostcell_dims[0];
                        
                        int idx_nghost = i + j*interior_dims[0];
                        
                        // Compute the gradient of density
                        double detector_rho = sqrt(pow(fabs(rho[idx_wghost_x[0]] - 2*rho[idx_wghost_x[1]] + rho[idx_wghost_x[2]])/
                                                   (fabs(rho[idx_wghost_x[1]] - rho[idx_wghost_x[0]])
                                                   + fabs(rho[idx_wghost_x[2]] - rho[idx_wghost_x[1]]) + 1.0e-40), 2.0)
                                                   + pow(fabs(rho[idx_wghost_y[0]] - 2*rho[idx_wghost_y[1]] + rho[idx_wghost_y[2]])/
                                                   (fabs(rho[idx_wghost_y[1]] - rho[idx_wghost_y[0]])
                                                   + fabs(rho[idx_wghost_y[2]] - rho[idx_wghost_y[1]]) + 1.0e-40), 2.0))/sqrt(2.0);
                        
                        NULL_USE(detector_rho);
                        
                        if (rho[idx_wghost_x[1]] > 6.5 && rho[idx_wghost_x[1]] < 8.5)
                        {
                            tag_ptr[idx_nghost] = 1;
                        }
                        else
                        {
                            tag_ptr[idx_nghost] = 0;
                        }
                    }
                }
            }
            else if (ref == "PRESSURE_SHOCK")
            {
                // NOT YET IMPLEMENTED
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        // NOT YET IMPLEMENTED
    }
    
    t_taggradient->stop();
}


void
Euler::setPhysicalBoundaryConditions(
    hier::Patch& patch,
    const double fill_time,
    const hier::IntVector& ghost_width_to_fill)
{
    NULL_USE(fill_time);
    t_setphysbcs->start();
    
    switch (d_flow_model)
    {
        case SINGLE_SPECIES:
        {
            boost::shared_ptr<pdat::CellData<double> > density(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_density, getDataContext())));
            
            boost::shared_ptr<pdat::CellData<double> > momentum(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_momentum, getDataContext())));
            
            boost::shared_ptr<pdat::CellData<double> > total_energy(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_total_energy, getDataContext())));
            
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(density);
            TBOX_ASSERT(momentum);
            TBOX_ASSERT(total_energy);
            
            TBOX_ASSERT(density->getGhostCellWidth() == d_num_ghosts);
            TBOX_ASSERT(momentum->getGhostCellWidth() == d_num_ghosts);
            TBOX_ASSERT(total_energy->getGhostCellWidth() == d_num_ghosts);
#endif
            
            if (d_dim == tbox::Dimension(1))
            {
                // NOT YET IMPLEMENTED
            }
            else if (d_dim == tbox::Dimension(2))
            {
                /*
                 * Set boundary conditions for cells corresponding to patch edges.
                 */
                
                appu::CartesianBoundaryUtilities2::fillEdgeBoundaryData(
                    "density",
                    density,
                    patch,
                    ghost_width_to_fill,
                    d_scalar_bdry_edge_conds,
                    d_bdry_edge_density);
                
                appu::CartesianBoundaryUtilities2::fillEdgeBoundaryData(
                    "momentum",
                    momentum,
                    patch,
                    ghost_width_to_fill,
                    d_vector_bdry_edge_conds,
                    d_bdry_edge_momentum);
                
                appu::CartesianBoundaryUtilities2::fillEdgeBoundaryData(
                    "total energy",
                    total_energy,
                    patch,
                    ghost_width_to_fill,
                    d_scalar_bdry_edge_conds,
                    d_bdry_edge_total_energy);

/*                
#ifdef DEBUG_CHECK_ASSERTIONS
#if CHECK_BDRY_DATA
                checkBoundaryData(Bdry::EDGE2D, patch, ghost_width_to_fill,
                    tmp_edge_scalar_bcond, tmp_edge_vector_bcond);
#endif
#endif
*/
                
                /*
                 *  Set boundary conditions for cells corresponding to patch nodes.
                 */
                
                appu::CartesianBoundaryUtilities2::fillNodeBoundaryData(
                    "density",
                    density,
                    patch,
                    ghost_width_to_fill,
                    d_scalar_bdry_node_conds,
                    d_bdry_edge_density);
                
                appu::CartesianBoundaryUtilities2::fillNodeBoundaryData(
                    "momentum",
                    momentum,
                    patch,
                    ghost_width_to_fill,
                    d_vector_bdry_node_conds,
                    d_bdry_edge_momentum);
                
                appu::CartesianBoundaryUtilities2::fillNodeBoundaryData(
                    "total energy",
                    total_energy,
                    patch,
                    ghost_width_to_fill,
                    d_scalar_bdry_node_conds,
                    d_bdry_edge_total_energy);
                
/*
#ifdef DEBUG_CHECK_ASSERTIONS
#if CHECK_BDRY_DATA
                checkBoundaryData(Bdry::NODE2D, patch, ghost_width_to_fill,
                    d_scalar_bdry_node_conds, d_vector_bdry_node_conds);
#endif
#endif
*/
            }
            else if (d_dim == tbox::Dimension(3))
            {
                /*
                 *  Set boundary conditions for cells corresponding to patch faces.
                 */
                
                appu::CartesianBoundaryUtilities3::fillFaceBoundaryData(
                    "density",
                    density,
                    patch,
                    ghost_width_to_fill,
                    d_scalar_bdry_face_conds,
                    d_bdry_face_density);
                
                appu::CartesianBoundaryUtilities3::fillFaceBoundaryData(
                    "momentum",
                    momentum,
                    patch,
                    ghost_width_to_fill,
                    d_vector_bdry_face_conds,
                    d_bdry_face_momentum);
                
                appu::CartesianBoundaryUtilities3::fillFaceBoundaryData(
                    "total energy",
                    total_energy,
                    patch,
                    ghost_width_to_fill,
                    d_scalar_bdry_face_conds,
                    d_bdry_face_total_energy);
                
/*
#ifdef DEBUG_CHECK_ASSERTIONS
#if CHECK_BDRY_DATA
                checkBoundaryData(Bdry::FACE3D, patch, ghost_width_to_fill,
                    d_scalar_bdry_face_conds, d_vector_bdry_face_conds);
#endif
#endif
*/
                
                /*
                 *  Set boundary conditions for cells corresponding to patch edges.
                 */
                
                appu::CartesianBoundaryUtilities3::fillEdgeBoundaryData(
                    "density",
                    density,
                    patch,
                    ghost_width_to_fill,
                    d_scalar_bdry_edge_conds,
                    d_bdry_face_density);
                
                appu::CartesianBoundaryUtilities3::fillEdgeBoundaryData(
                    "momentum",
                    momentum,
                    patch,
                    ghost_width_to_fill,
                    d_vector_bdry_edge_conds,
                    d_bdry_face_momentum);
                
                appu::CartesianBoundaryUtilities3::fillEdgeBoundaryData(
                    "total energy",
                    total_energy,
                    patch,
                    ghost_width_to_fill,
                    d_scalar_bdry_edge_conds,
                    d_bdry_face_total_energy);

/*
#ifdef DEBUG_CHECK_ASSERTIONS
#if CHECK_BDRY_DATA
      checkBoundaryData(Bdry::EDGE3D, patch, ghost_width_to_fill,
         d_scalar_bdry_edge_conds, d_vector_bdry_edge_conds);
#endif
#endif
*/
                
                /*
                 *  Set boundary conditions for cells corresponding to patch nodes.
                 */
                
                appu::CartesianBoundaryUtilities3::fillNodeBoundaryData(
                    "density",
                    density,
                    patch,
                    ghost_width_to_fill,
                    d_scalar_bdry_node_conds,
                    d_bdry_face_density);
                
                appu::CartesianBoundaryUtilities3::fillNodeBoundaryData(
                    "momentum",
                    momentum,
                    patch,
                    ghost_width_to_fill,
                    d_vector_bdry_node_conds,
                    d_bdry_face_momentum);
                
                appu::CartesianBoundaryUtilities3::fillNodeBoundaryData(
                    "total energy",
                    total_energy,
                    patch,
                    ghost_width_to_fill,
                    d_scalar_bdry_node_conds,
                    d_bdry_face_total_energy);
                
/*
#ifdef DEBUG_CHECK_ASSERTIONS
#if CHECK_BDRY_DATA
      checkBoundaryData(Bdry::NODE3D, patch, ghost_width_to_fill,
         d_scalar_bdry_node_conds, d_scalar_bdry_node_conds);
#endif
#endif
*/
            }
            
            break;
        }
        case FOUR_EQN_SHYUE:
        {
            boost::shared_ptr<pdat::CellData<double> > density(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_density, getDataContext())));
            
            boost::shared_ptr<pdat::CellData<double> > momentum(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_momentum, getDataContext())));
            
            boost::shared_ptr<pdat::CellData<double> > total_energy(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_total_energy, getDataContext())));
            
            boost::shared_ptr<pdat::CellData<double> > mass_fraction(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_mass_fraction, getDataContext())));
            
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(density);
            TBOX_ASSERT(momentum);
            TBOX_ASSERT(total_energy);
            TBOX_ASSERT(mass_fraction);
            
            TBOX_ASSERT(density->getGhostCellWidth() == d_num_ghosts);
            TBOX_ASSERT(momentum->getGhostCellWidth() == d_num_ghosts);
            TBOX_ASSERT(total_energy->getGhostCellWidth() == d_num_ghosts);
            TBOX_ASSERT(mass_fraction->getGhostCellWidth() == d_num_ghosts);
#endif
            
            if (d_dim == tbox::Dimension(1))
            {
                // NOT YET IMPLEMENTED
            }
            else if (d_dim == tbox::Dimension(2))
            {
                /*
                 * Set boundary conditions for cells corresponding to patch edges.
                 */
                
                appu::CartesianBoundaryUtilities2::fillEdgeBoundaryData(
                    "density",
                    density,
                    patch,
                    ghost_width_to_fill,
                    d_scalar_bdry_edge_conds,
                    d_bdry_edge_density);
                
                appu::CartesianBoundaryUtilities2::fillEdgeBoundaryData(
                    "momentum",
                    momentum,
                    patch,
                    ghost_width_to_fill,
                    d_vector_bdry_edge_conds,
                    d_bdry_edge_momentum);
                
                appu::CartesianBoundaryUtilities2::fillEdgeBoundaryData(
                    "total energy",
                    total_energy,
                    patch,
                    ghost_width_to_fill,
                    d_scalar_bdry_edge_conds,
                    d_bdry_edge_total_energy);
                
                appu::CartesianBoundaryUtilities2::fillEdgeBoundaryData(
                    "mass fraction",
                    mass_fraction,
                    patch,
                    ghost_width_to_fill,
                    d_scalar_bdry_edge_conds,
                    d_bdry_edge_mass_fraction);

/*                
#ifdef DEBUG_CHECK_ASSERTIONS
#if CHECK_BDRY_DATA
                checkBoundaryData(Bdry::EDGE2D, patch, ghost_width_to_fill,
                    tmp_edge_scalar_bcond, tmp_edge_vector_bcond);
#endif
#endif
*/
                
                /*
                 *  Set boundary conditions for cells corresponding to patch nodes.
                 */
                
                appu::CartesianBoundaryUtilities2::fillNodeBoundaryData(
                    "density",
                    density,
                    patch,
                    ghost_width_to_fill,
                    d_scalar_bdry_node_conds,
                    d_bdry_edge_density);
                
                appu::CartesianBoundaryUtilities2::fillNodeBoundaryData(
                    "momentum",
                    momentum,
                    patch,
                    ghost_width_to_fill,
                    d_vector_bdry_node_conds,
                    d_bdry_edge_momentum);
                
                appu::CartesianBoundaryUtilities2::fillNodeBoundaryData(
                    "total energy",
                    total_energy,
                    patch,
                    ghost_width_to_fill,
                    d_scalar_bdry_node_conds,
                    d_bdry_edge_total_energy);
                
                appu::CartesianBoundaryUtilities2::fillNodeBoundaryData(
                    "mass fraction",
                    mass_fraction,
                    patch,
                    ghost_width_to_fill,
                    d_scalar_bdry_node_conds,
                    d_bdry_edge_mass_fraction);
                
/*
#ifdef DEBUG_CHECK_ASSERTIONS
#if CHECK_BDRY_DATA
                checkBoundaryData(Bdry::NODE2D, patch, ghost_width_to_fill,
                    d_scalar_bdry_node_conds, d_vector_bdry_node_conds);
#endif
#endif
*/
            }
            else if (d_dim == tbox::Dimension(3))
            {
                /*
                 *  Set boundary conditions for cells corresponding to patch faces.
                 */
                
                appu::CartesianBoundaryUtilities3::fillFaceBoundaryData(
                    "density",
                    density,
                    patch,
                    ghost_width_to_fill,
                    d_scalar_bdry_face_conds,
                    d_bdry_face_density);
                
                appu::CartesianBoundaryUtilities3::fillFaceBoundaryData(
                    "momentum",
                    momentum,
                    patch,
                    ghost_width_to_fill,
                    d_vector_bdry_face_conds,
                    d_bdry_face_momentum);
                
                appu::CartesianBoundaryUtilities3::fillFaceBoundaryData(
                    "total energy",
                    total_energy,
                    patch,
                    ghost_width_to_fill,
                    d_scalar_bdry_face_conds,
                    d_bdry_face_total_energy);
                
                appu::CartesianBoundaryUtilities3::fillFaceBoundaryData(
                    "mass fraction",
                    mass_fraction,
                    patch,
                    ghost_width_to_fill,
                    d_scalar_bdry_face_conds,
                    d_bdry_face_mass_fraction);
                
/*
#ifdef DEBUG_CHECK_ASSERTIONS
#if CHECK_BDRY_DATA
                checkBoundaryData(Bdry::FACE3D, patch, ghost_width_to_fill,
                    d_scalar_bdry_face_conds, d_vector_bdry_face_conds);
#endif
#endif
*/
                
                /*
                 *  Set boundary conditions for cells corresponding to patch edges.
                 */
                
                appu::CartesianBoundaryUtilities3::fillEdgeBoundaryData(
                    "density",
                    density,
                    patch,
                    ghost_width_to_fill,
                    d_scalar_bdry_edge_conds,
                    d_bdry_face_density);
                
                appu::CartesianBoundaryUtilities3::fillEdgeBoundaryData(
                    "momentum",
                    momentum,
                    patch,
                    ghost_width_to_fill,
                    d_vector_bdry_edge_conds,
                    d_bdry_face_momentum);
                
                appu::CartesianBoundaryUtilities3::fillEdgeBoundaryData(
                    "total energy",
                    total_energy,
                    patch,
                    ghost_width_to_fill,
                    d_scalar_bdry_edge_conds,
                    d_bdry_face_total_energy);
                
                appu::CartesianBoundaryUtilities3::fillEdgeBoundaryData(
                    "mass fraction",
                    mass_fraction,
                    patch,
                    ghost_width_to_fill,
                    d_scalar_bdry_edge_conds,
                    d_bdry_face_mass_fraction);

/*
#ifdef DEBUG_CHECK_ASSERTIONS
#if CHECK_BDRY_DATA
      checkBoundaryData(Bdry::EDGE3D, patch, ghost_width_to_fill,
         d_scalar_bdry_edge_conds, d_vector_bdry_edge_conds);
#endif
#endif
*/
                
                /*
                 *  Set boundary conditions for cells corresponding to patch nodes.
                 */
                
                appu::CartesianBoundaryUtilities3::fillNodeBoundaryData(
                    "density",
                    density,
                    patch,
                    ghost_width_to_fill,
                    d_scalar_bdry_node_conds,
                    d_bdry_face_density);
                
                appu::CartesianBoundaryUtilities3::fillNodeBoundaryData(
                    "momentum",
                    momentum,
                    patch,
                    ghost_width_to_fill,
                    d_vector_bdry_node_conds,
                    d_bdry_face_momentum);
                
                appu::CartesianBoundaryUtilities3::fillNodeBoundaryData(
                    "total energy",
                    total_energy,
                    patch,
                    ghost_width_to_fill,
                    d_scalar_bdry_node_conds,
                    d_bdry_face_total_energy);
                
                appu::CartesianBoundaryUtilities3::fillNodeBoundaryData(
                    "mass fraction",
                    mass_fraction,
                    patch,
                    ghost_width_to_fill,
                    d_scalar_bdry_node_conds,
                    d_bdry_face_mass_fraction);
                
/*
#ifdef DEBUG_CHECK_ASSERTIONS
#if CHECK_BDRY_DATA
      checkBoundaryData(Bdry::NODE3D, patch, ghost_width_to_fill,
         d_scalar_bdry_node_conds, d_scalar_bdry_node_conds);
#endif
#endif
*/
            }
            
            break;
        }
        case FIVE_EQN_ALLAIRE:
        {
            boost::shared_ptr<pdat::CellData<double> > partial_density(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_partial_density, getDataContext())));
            
            boost::shared_ptr<pdat::CellData<double> > momentum(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_momentum, getDataContext())));
            
            boost::shared_ptr<pdat::CellData<double> > total_energy(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_total_energy, getDataContext())));
            
            boost::shared_ptr<pdat::CellData<double> > volume_fraction(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_volume_fraction, getDataContext())));
            
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(partial_density);
            TBOX_ASSERT(momentum);
            TBOX_ASSERT(total_energy);
            TBOX_ASSERT(volume_fraction);
            
            TBOX_ASSERT(partial_density->getGhostCellWidth() == d_num_ghosts);
            TBOX_ASSERT(momentum->getGhostCellWidth() == d_num_ghosts);
            TBOX_ASSERT(total_energy->getGhostCellWidth() == d_num_ghosts);
            TBOX_ASSERT(volume_fraction->getGhostCellWidth() == d_num_ghosts);
#endif
            
            if (d_dim == tbox::Dimension(1))
            {
                // NOT YET IMPLEMENTED
            }
            else if (d_dim == tbox::Dimension(2))
            {
                /*
                 * Set boundary conditions for cells corresponding to patch edges.
                 */
                
                appu::CartesianBoundaryUtilities2::fillEdgeBoundaryData(
                    "partial density",
                    partial_density,
                    patch,
                    ghost_width_to_fill,
                    d_scalar_bdry_edge_conds,
                    d_bdry_edge_partial_density);
                
                appu::CartesianBoundaryUtilities2::fillEdgeBoundaryData(
                    "momentum",
                    momentum,
                    patch,
                    ghost_width_to_fill,
                    d_vector_bdry_edge_conds,
                    d_bdry_edge_momentum);
                
                appu::CartesianBoundaryUtilities2::fillEdgeBoundaryData(
                    "total energy",
                    total_energy,
                    patch,
                    ghost_width_to_fill,
                    d_scalar_bdry_edge_conds,
                    d_bdry_edge_total_energy);
                
                appu::CartesianBoundaryUtilities2::fillEdgeBoundaryData(
                    "volume fraction",
                    volume_fraction,
                    patch,
                    ghost_width_to_fill,
                    d_scalar_bdry_edge_conds,
                    d_bdry_edge_volume_fraction);

/*                
#ifdef DEBUG_CHECK_ASSERTIONS
#if CHECK_BDRY_DATA
                checkBoundaryData(Bdry::EDGE2D, patch, ghost_width_to_fill,
                    tmp_edge_scalar_bcond, tmp_edge_vector_bcond);
#endif
#endif
*/
                
                /*
                 *  Set boundary conditions for cells corresponding to patch nodes.
                 */
                
                appu::CartesianBoundaryUtilities2::fillNodeBoundaryData(
                    "partial density",
                    partial_density,
                    patch,
                    ghost_width_to_fill,
                    d_scalar_bdry_node_conds,
                    d_bdry_edge_partial_density);
                
                appu::CartesianBoundaryUtilities2::fillNodeBoundaryData(
                    "momentum",
                    momentum,
                    patch,
                    ghost_width_to_fill,
                    d_vector_bdry_node_conds,
                    d_bdry_edge_momentum);
                
                appu::CartesianBoundaryUtilities2::fillNodeBoundaryData(
                    "total energy",
                    total_energy,
                    patch,
                    ghost_width_to_fill,
                    d_scalar_bdry_node_conds,
                    d_bdry_edge_total_energy);
                
                appu::CartesianBoundaryUtilities2::fillNodeBoundaryData(
                    "volume fraction",
                    volume_fraction,
                    patch,
                    ghost_width_to_fill,
                    d_scalar_bdry_node_conds,
                    d_bdry_edge_volume_fraction);
                
/*
#ifdef DEBUG_CHECK_ASSERTIONS
#if CHECK_BDRY_DATA
                checkBoundaryData(Bdry::NODE2D, patch, ghost_width_to_fill,
                    d_scalar_bdry_node_conds, d_vector_bdry_node_conds);
#endif
#endif
*/
            }
            else if (d_dim == tbox::Dimension(3))
            {
                /*
                 *  Set boundary conditions for cells corresponding to patch faces.
                 */
                
                appu::CartesianBoundaryUtilities3::fillFaceBoundaryData(
                    "partial density",
                    partial_density,
                    patch,
                    ghost_width_to_fill,
                    d_scalar_bdry_face_conds,
                    d_bdry_face_partial_density);
                
                appu::CartesianBoundaryUtilities3::fillFaceBoundaryData(
                    "momentum",
                    momentum,
                    patch,
                    ghost_width_to_fill,
                    d_vector_bdry_face_conds,
                    d_bdry_face_momentum);
                
                appu::CartesianBoundaryUtilities3::fillFaceBoundaryData(
                    "total energy",
                    total_energy,
                    patch,
                    ghost_width_to_fill,
                    d_scalar_bdry_face_conds,
                    d_bdry_face_total_energy);
                
                appu::CartesianBoundaryUtilities3::fillFaceBoundaryData(
                    "volume fraction",
                    volume_fraction,
                    patch,
                    ghost_width_to_fill,
                    d_scalar_bdry_face_conds,
                    d_bdry_face_volume_fraction);
                
/*
#ifdef DEBUG_CHECK_ASSERTIONS
#if CHECK_BDRY_DATA
                checkBoundaryData(Bdry::FACE3D, patch, ghost_width_to_fill,
                    d_scalar_bdry_face_conds, d_vector_bdry_face_conds);
#endif
#endif
*/
                
                /*
                 *  Set boundary conditions for cells corresponding to patch edges.
                 */
                
                appu::CartesianBoundaryUtilities3::fillEdgeBoundaryData(
                    "partial density",
                    partial_density,
                    patch,
                    ghost_width_to_fill,
                    d_scalar_bdry_edge_conds,
                    d_bdry_face_partial_density);
                
                appu::CartesianBoundaryUtilities3::fillEdgeBoundaryData(
                    "momentum",
                    momentum,
                    patch,
                    ghost_width_to_fill,
                    d_vector_bdry_edge_conds,
                    d_bdry_face_momentum);
                
                appu::CartesianBoundaryUtilities3::fillEdgeBoundaryData(
                    "total energy",
                    total_energy,
                    patch,
                    ghost_width_to_fill,
                    d_scalar_bdry_edge_conds,
                    d_bdry_face_total_energy);
                
                appu::CartesianBoundaryUtilities3::fillEdgeBoundaryData(
                    "volume fraction",
                    volume_fraction,
                    patch,
                    ghost_width_to_fill,
                    d_scalar_bdry_edge_conds,
                    d_bdry_face_volume_fraction);

/*
#ifdef DEBUG_CHECK_ASSERTIONS
#if CHECK_BDRY_DATA
      checkBoundaryData(Bdry::EDGE3D, patch, ghost_width_to_fill,
         d_scalar_bdry_edge_conds, d_vector_bdry_edge_conds);
#endif
#endif
*/
                
                /*
                 *  Set boundary conditions for cells corresponding to patch nodes.
                 */
                
                appu::CartesianBoundaryUtilities3::fillNodeBoundaryData(
                    "partial density",
                    partial_density,
                    patch,
                    ghost_width_to_fill,
                    d_scalar_bdry_node_conds,
                    d_bdry_face_partial_density);
                
                appu::CartesianBoundaryUtilities3::fillNodeBoundaryData(
                    "momentum",
                    momentum,
                    patch,
                    ghost_width_to_fill,
                    d_vector_bdry_node_conds,
                    d_bdry_face_momentum);
                
                appu::CartesianBoundaryUtilities3::fillNodeBoundaryData(
                    "total energy",
                    total_energy,
                    patch,
                    ghost_width_to_fill,
                    d_scalar_bdry_node_conds,
                    d_bdry_face_total_energy);
                
                appu::CartesianBoundaryUtilities3::fillNodeBoundaryData(
                    "volume fraction",
                    volume_fraction,
                    patch,
                    ghost_width_to_fill,
                    d_scalar_bdry_node_conds,
                    d_bdry_face_volume_fraction);
                
/*
#ifdef DEBUG_CHECK_ASSERTIONS
#if CHECK_BDRY_DATA
      checkBoundaryData(Bdry::NODE3D, patch, ghost_width_to_fill,
         d_scalar_bdry_node_conds, d_scalar_bdry_node_conds);
#endif
#endif
*/
            }
            
            break;
        }
        default:
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "Unknown d_flow_model."
                       << std::endl);
        }
    }
    
    if (d_project_name == "2D double-Mach reflection")
    {
        TBOX_ASSERT(d_dim == tbox::Dimension(2));
        TBOX_ASSERT(d_flow_model == SINGLE_SPECIES);
        
        // Get the dimensions of box that covers the interior of patch.
        hier::Box dummy_box = patch.getBox();
        const hier::Box interior_box = dummy_box;
        const hier::IntVector interior_dims = interior_box.numberCells();
        
        // Get the dimensions of box that covers interior of patch plus
        // ghost cells.
        dummy_box.grow(d_num_ghosts);
        const hier::Box ghost_box = dummy_box;
        const hier::IntVector ghostcell_dims = ghost_box.numberCells();
        
        const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
            BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                patch.getPatchGeometry()));
        
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(patch_geom);
#endif
        
        if (patch_geom->getTouchesRegularBoundary(1, 0) ||
            patch_geom->getTouchesRegularBoundary(1, 1))
        {
            const double* const dx = patch_geom->getDx();
            const double* const patch_xlo = patch_geom->getXLower();
            
            boost::shared_ptr<pdat::CellData<double> > density(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_density, getDataContext())));
            
            boost::shared_ptr<pdat::CellData<double> > momentum(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_momentum, getDataContext())));
            
            boost::shared_ptr<pdat::CellData<double> > total_energy(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_total_energy, getDataContext())));
            
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(density);
            TBOX_ASSERT(momentum);
            TBOX_ASSERT(total_energy);
            
            TBOX_ASSERT(density->getGhostCellWidth() == d_num_ghosts);
            TBOX_ASSERT(momentum->getGhostCellWidth() == d_num_ghosts);
            TBOX_ASSERT(total_energy->getGhostCellWidth() == d_num_ghosts);
#endif
            
            double* rho   = density->getPointer(0);
            double* rho_u = momentum->getPointer(0);
            double* rho_v = momentum->getPointer(1);
            double* E     = total_energy->getPointer(0);
            
            const double x_0 = 1.0/6.0;
            
            const double gamma = 1.4;
            
            const double rho_post_shock = 8.0;
            const double u_post_shock = 8.25*cos(M_PI/6.0);
            const double v_post_shock = -8.25*sin(M_PI/6.0);
            const double p_post_shock = 116.5;
            
            const double rho_pre_shock = 1.4;
            const double u_pre_shock = 0.0;
            const double v_pre_shock = 0.0;
            const double p_pre_shock = 1.0;
            
            const double rho_u_post_shock = rho_post_shock*u_post_shock;
            const double rho_v_post_shock = rho_post_shock*v_post_shock;
            
            const double rho_u_pre_shock = rho_pre_shock*u_pre_shock;
            const double rho_v_pre_shock = rho_pre_shock*v_pre_shock;
            
            const double E_pre_shock = p_pre_shock/(gamma - 1.0) +
                0.5*rho_pre_shock*(u_pre_shock*u_pre_shock + v_pre_shock*v_pre_shock);
            
            const double E_post_shock = p_post_shock/(gamma - 1.0) +
                0.5*rho_post_shock*(u_post_shock*u_post_shock + v_post_shock*v_post_shock);
            
            /*
             * Update the bottom boundary conditions.
             */
            
            if (patch_geom->getTouchesRegularBoundary(1, 0))
            {
                for (int i = 0; i < interior_dims[0]; i++)
                {
                    for (int j = -ghost_width_to_fill[1];
                         j < 0;
                         j++)
                    {
                        const int idx_cell = (i + d_num_ghosts[0]) +
                            (j + d_num_ghosts[1])*ghostcell_dims[0];
                        
                        // Compute the coordinates.
                        double x[2];
                        x[0] = patch_xlo[0] + (i + 0.5)*dx[0];
                        x[1] = patch_xlo[1] + (j + 0.5)*dx[1];
                        
                        if (x[0] < x_0)
                        {
                            rho[idx_cell] = rho_post_shock;
                            rho_u[idx_cell] = rho_u_post_shock;
                            rho_v[idx_cell] = rho_v_post_shock;
                            E[idx_cell] = E_post_shock;
                        }
                        else
                        {
                            const int idx_mirror_cell = (i + d_num_ghosts[0]) +
                                (-j + d_num_ghosts[1] - 1)*ghostcell_dims[0];
                            
                            rho[idx_cell] = rho[idx_mirror_cell];
                            rho_u[idx_cell] = rho_u[idx_mirror_cell];
                            rho_v[idx_cell] = -rho_v[idx_mirror_cell];
                            E[idx_cell] = E[idx_mirror_cell];
                        }
                    }
                }
            }
            
            /*
             * Update the top boundary conditions.
             */
            
            if (patch_geom->getTouchesRegularBoundary(1, 1))
            {
                const double x_s = x_0 + (1.0 + 20.0*fill_time)/sqrt(3.0);
                
                for (int i = 0; i < interior_dims[0]; i++)
                {
                    for (int j = interior_dims[1];
                         j < interior_dims[1] + ghost_width_to_fill[1];
                         j++)
                    {
                        const int idx_cell = (i + d_num_ghosts[0]) +
                            (j + d_num_ghosts[1])*ghostcell_dims[0];
                        
                        // Compute the coordinates.
                        double x[2];
                        x[0] = patch_xlo[0] + (i + 0.5)*dx[0];
                        x[1] = patch_xlo[1] + (j + 0.5)*dx[1];
                        
                        if (x[0] >= x_s)
                        {
                            rho[idx_cell] = rho_pre_shock;
                            rho_u[idx_cell] = rho_u_pre_shock;
                            rho_v[idx_cell] = rho_v_pre_shock;
                            E[idx_cell] = E_pre_shock;
                        }
                        else
                        {
                            rho[idx_cell] = rho_post_shock;
                            rho_u[idx_cell] = rho_u_post_shock;
                            rho_v[idx_cell] = rho_v_post_shock;
                            E[idx_cell] = E_post_shock;
                        }
                    }
                }
            }
        }
    }
    
    t_setphysbcs->stop();
}


void
Euler::putToRestart(
    const boost::shared_ptr<tbox::Database>& restart_db) const
{
    TBOX_ASSERT(restart_db);
    
    restart_db->putString("d_project_name", d_project_name);
    
    restart_db->putInteger("d_num_species", d_num_species);
    
    switch(d_flow_model)
    {
        case SINGLE_SPECIES:
            restart_db->putString("d_flow_model", "SINGLE_SPECIES");
            break;
        case FOUR_EQN_SHYUE:
            restart_db->putString("d_flow_model", "FOUR_EQN_SHYUE");
            break;
        case FIVE_EQN_ALLAIRE:
            restart_db->putString("d_flow_model", "FIVE_EQN_ALLAIRE");
            break;
        default:
            TBOX_ERROR(d_object_name
                       << ": "
                       << "Unknown d_flow_model."
                       << std::endl);
    }
    
    restart_db->putBool("d_is_preserving_positivity", d_is_preserving_positivity);
    
    boost::shared_ptr<tbox::Database> restart_equation_of_state_db =
        restart_db->putDatabase("Equation_of_state");
    
    d_equation_of_state->putToRestart(restart_equation_of_state_db);
    
    boost::shared_ptr<tbox::Database> restart_shock_capturing_scheme_db =
        restart_db->putDatabase("Shock_capturing_scheme");
    
    d_conv_flux_reconstructor->putToRestart(restart_shock_capturing_scheme_db);
    
    restart_db->putIntegerArray("d_num_ghosts", &d_num_ghosts[0], d_dim.getValue());
    
    restart_db->putIntegerVector("d_master_bdry_node_conds",
                                 d_master_bdry_node_conds);
    
    if (d_dim == tbox::Dimension(1))
    {
        switch (d_flow_model)
        {
            case SINGLE_SPECIES:
            {
                restart_db->putDoubleVector("d_bdry_node_density",
                                            d_bdry_node_density);
                
                restart_db->putDoubleVector("d_bdry_node_momentum",
                                            d_bdry_node_momentum);
                
                restart_db->putDoubleVector("d_bdry_node_total_energy",
                                            d_bdry_node_total_energy);
                
                break;
            }
            case FOUR_EQN_SHYUE:
            {
                restart_db->putDoubleVector("d_bdry_node_density",
                                            d_bdry_node_density);
                
                restart_db->putDoubleVector("d_bdry_node_momentum",
                                            d_bdry_node_momentum);
                
                restart_db->putDoubleVector("d_bdry_node_total_energy",
                                            d_bdry_node_total_energy);
                
                restart_db->putDoubleVector("d_bdry_node_mass_fraction",
                                            d_bdry_node_mass_fraction);
                
                break;
            }
            case FIVE_EQN_ALLAIRE:
            {
                restart_db->putDoubleVector("d_bdry_node_partial_density",
                                            d_bdry_node_partial_density);
                
                restart_db->putDoubleVector("d_bdry_node_momentum",
                                            d_bdry_node_momentum);
                
                restart_db->putDoubleVector("d_bdry_node_total_energy",
                                            d_bdry_node_total_energy);
                
                restart_db->putDoubleVector("d_bdry_node_volume_fraction",
                                            d_bdry_node_volume_fraction);
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                           << ": "
                           << "Unknown d_flow_model."
                           << std::endl);
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        restart_db->putIntegerVector("d_master_bdry_edge_conds",
                                     d_master_bdry_edge_conds);
        
        switch (d_flow_model)
        {
            case SINGLE_SPECIES:
            {
                restart_db->putDoubleVector("d_bdry_edge_density",
                                            d_bdry_edge_density);
                
                restart_db->putDoubleVector("d_bdry_edge_momentum",
                                            d_bdry_edge_momentum);
                
                restart_db->putDoubleVector("d_bdry_edge_total_energy",
                                            d_bdry_edge_total_energy);
                
                break;
            }
            case FOUR_EQN_SHYUE:
            {
                restart_db->putDoubleVector("d_bdry_edge_density",
                                            d_bdry_edge_density);
                
                restart_db->putDoubleVector("d_bdry_edge_momentum",
                                            d_bdry_edge_momentum);
                
                restart_db->putDoubleVector("d_bdry_edge_total_energy",
                                            d_bdry_edge_total_energy);
                
                restart_db->putDoubleVector("d_bdry_edge_mass_fraction",
                                            d_bdry_edge_mass_fraction);
                
                break;
            }
            case FIVE_EQN_ALLAIRE:
            {
                restart_db->putDoubleVector("d_bdry_edge_partial_density",
                                            d_bdry_edge_partial_density);
                
                restart_db->putDoubleVector("d_bdry_edge_momentum",
                                            d_bdry_edge_momentum);
                
                restart_db->putDoubleVector("d_bdry_edge_total_energy",
                                            d_bdry_edge_total_energy);
                
                restart_db->putDoubleVector("d_bdry_edge_volume_fraction",
                                            d_bdry_edge_volume_fraction);
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                           << ": "
                           << "Unknown d_flow_model."
                           << std::endl);
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        restart_db->putIntegerVector("d_master_bdry_edge_conds",
                                     d_master_bdry_edge_conds);
        
        restart_db->putIntegerVector("d_master_bdry_face_conds",
                                     d_master_bdry_face_conds);
        
        switch (d_flow_model)
        {
            case SINGLE_SPECIES:
            {
                restart_db->putDoubleVector("d_bdry_face_density",
                                            d_bdry_face_density);
                
                restart_db->putDoubleVector("d_bdry_face_momentum",
                                            d_bdry_face_momentum);
                
                restart_db->putDoubleVector("d_bdry_face_total_energy",
                                            d_bdry_face_total_energy);
                
                break;
            }
            case FOUR_EQN_SHYUE:
            {
                restart_db->putDoubleVector("d_bdry_face_density",
                                            d_bdry_face_density);
                
                restart_db->putDoubleVector("d_bdry_face_momentum",
                                            d_bdry_face_momentum);
                
                restart_db->putDoubleVector("d_bdry_face_total_energy",
                                            d_bdry_face_total_energy);
                
                restart_db->putDoubleVector("d_bdry_face_mass_fraction",
                                            d_bdry_face_mass_fraction);
                
                break;
            }
            case FIVE_EQN_ALLAIRE:
            {
                restart_db->putDoubleVector("d_bdry_face_partial_density",
                                            d_bdry_face_partial_density);
                
                restart_db->putDoubleVector("d_bdry_face_momentum",
                                            d_bdry_face_momentum);
                
                restart_db->putDoubleVector("d_bdry_face_total_energy",
                                            d_bdry_face_total_energy);
                
                restart_db->putDoubleVector("d_bdry_face_volume_fraction",
                                            d_bdry_face_volume_fraction);
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                           << ": "
                           << "Unknown d_flow_model."
                           << std::endl);
            }
        }
    }
    
    if (d_refinement_criteria.size() > 0)
    {
        restart_db->putStringVector("d_refinement_criteria", d_refinement_criteria);
    }
    for (int i = 0; i < static_cast<int>(d_refinement_criteria.size()); i++)
    {
        if (d_refinement_criteria[i] == "DENSITY_SHOCK")
        {
            restart_db->putDoubleVector("d_density_shock_tol",
                                        d_density_shock_tol);
        }
        else if(d_refinement_criteria[i] == "PRESSURE_SHOCK")
        {
            restart_db->putDoubleVector("d_pressure_shock_tol",
                                        d_pressure_shock_tol);
        }
    }
}


void
Euler::readDirichletBoundaryDataEntry(
    const boost::shared_ptr<tbox::Database>& db,
    std::string& db_name,
    int bdry_location_index)
{
    TBOX_ASSERT(db);
    TBOX_ASSERT(!db_name.empty());
    
    if (d_dim == tbox::Dimension(1))
    {
        // NOT YET IMPLEMENTED
    }
    else if (d_dim == tbox::Dimension(2))
    {
        switch (d_flow_model)
        {
            case SINGLE_SPECIES:
            {
                readStateDataEntryForSingleSpecies(
                    db,
                    db_name,
                    bdry_location_index,
                    d_bdry_edge_density,
                    d_bdry_edge_momentum,
                    d_bdry_edge_total_energy);
                
                break;
            }
            case FOUR_EQN_SHYUE:
            {
                readStateDataEntryForFourEqnShyue(
                    db,
                    db_name,
                    bdry_location_index,
                    d_bdry_edge_density,
                    d_bdry_edge_momentum,
                    d_bdry_edge_total_energy,
                    d_bdry_edge_mass_fraction);
                
                break;
            }
            case FIVE_EQN_ALLAIRE:
            {
                readStateDataEntryForFiveEqnAllaire(
                    db,
                    db_name,
                    bdry_location_index,
                    d_bdry_edge_partial_density,
                    d_bdry_edge_momentum,
                    d_bdry_edge_total_energy,
                    d_bdry_edge_volume_fraction);
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                           << ": "
                           << "Unknown d_flow_model."
                           << std::endl);
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        switch (d_flow_model)
        {
            case SINGLE_SPECIES:
            {
                readStateDataEntryForSingleSpecies(
                    db,
                    db_name,
                    bdry_location_index,
                    d_bdry_face_density,
                    d_bdry_face_momentum,
                    d_bdry_face_total_energy);
                
                break;
            }
            case FOUR_EQN_SHYUE:
            {
                readStateDataEntryForFourEqnShyue(
                    db,
                    db_name,
                    bdry_location_index,
                    d_bdry_face_density,
                    d_bdry_face_momentum,
                    d_bdry_face_total_energy,
                    d_bdry_face_mass_fraction);
                
                break;
            }
            case FIVE_EQN_ALLAIRE:
            {
                readStateDataEntryForFiveEqnAllaire(
                    db,
                    db_name,
                    bdry_location_index,
                    d_bdry_face_partial_density,
                    d_bdry_face_momentum,
                    d_bdry_face_total_energy,
                    d_bdry_face_volume_fraction);
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                           << ": "
                           << "Unknown d_flow_model."
                           << std::endl);
            }
        }
    }
}


void
Euler::readNeumannBoundaryDataEntry(
    const boost::shared_ptr<tbox::Database>& db,
    std::string& db_name,
    int bdry_location_index)
{
    NULL_USE(db);
    NULL_USE(db_name);
    NULL_USE(bdry_location_index);
}


#ifdef HAVE_HDF5
void
Euler::registerVisItDataWriter(
    boost::shared_ptr<appu::VisItDataWriter> viz_writer)
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
    NULL_USE(simulation_time);
    
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT((region * patch.getBox()).isSpatiallyEqual(region));
#endif
    
    bool data_on_patch = false;
    
    // Get the dimensions of the region.
    const hier::IntVector region_dims = region.numberCells();
    
    if (variable_name == "pressure")
    {
        switch (d_flow_model)
        {
            case SINGLE_SPECIES:
            {
                boost::shared_ptr<pdat::CellData<double> > density(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_density, d_plot_context)));
                
                boost::shared_ptr<pdat::CellData<double> > momentum(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_momentum, d_plot_context)));
                
                boost::shared_ptr<pdat::CellData<double> > total_energy(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_total_energy, d_plot_context)));
                
#ifdef DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(density);
                TBOX_ASSERT(momentum);
                TBOX_ASSERT(total_energy);
                TBOX_ASSERT(density->getGhostBox().isSpatiallyEqual(patch.getBox()));
                TBOX_ASSERT(momentum->getGhostBox().isSpatiallyEqual(patch.getBox()));
                TBOX_ASSERT(total_energy->getGhostBox().isSpatiallyEqual(patch.getBox()));
#endif
                
                // Get the dimensions of box that covers the data.
                const hier::Box data_box = density->getGhostBox();
                const hier::IntVector data_dims = data_box.numberCells();
                
                // Get the arrays of time-dependent variables
                const double* const rho   = density->getPointer();
                const double* const rho_u = momentum->getPointer(0);
                const double* const rho_v = d_dim > tbox::Dimension(1) ? momentum->getPointer(1) : NULL;
                const double* const rho_w = d_dim > tbox::Dimension(2) ? momentum->getPointer(2) : NULL;
                const double* const E     = total_energy->getPointer(0);
                
                size_t offset_data = data_box.offset(region.lower());
                
                if (d_dim == tbox::Dimension(1))
                {
                    for (int i = 0; i < region_dims[0]; i++)
                    {
                        size_t idx_data = offset_data + i;
                        
                        size_t idx_region = i;
                        
                        std::vector<const double*> m_ptr;
                        m_ptr.push_back(&rho_u[idx_data]);
                        
                        buffer[idx_region] = d_equation_of_state->getPressure(
                            &rho[idx_data],
                            m_ptr,
                            &E[idx_data]);
                    }
                }
                else if (d_dim == tbox::Dimension(2))
                {
                    for (int j = 0; j < region_dims[1]; j++)
                    {
                        for (int i = 0; i < region_dims[0]; i++)
                        {
                            size_t idx_data = offset_data + i +
                                j*data_dims[0];
                            
                            size_t idx_region = i +
                                j*region_dims[0];
                            
                            std::vector<const double*> m_ptr;
                            m_ptr.push_back(&rho_u[idx_data]);
                            m_ptr.push_back(&rho_v[idx_data]);
                            
                            buffer[idx_region] = d_equation_of_state->getPressure(
                                &rho[idx_data],
                                m_ptr,
                                &E[idx_data]);
                        }
                    }
                }
                else if (d_dim == tbox::Dimension(3))
                {
                    for (int k = 0; k < region_dims[2]; k++)
                    {
                        for (int j = 0; j < region_dims[1]; j++)
                        {
                            for (int i = 0; i < region_dims[0]; i++)
                            {
                                size_t idx_data = offset_data + i +
                                    j*data_dims[0] +
                                    k*data_dims[0]*data_dims[1];
                                
                                size_t idx_region = i +
                                    j*region_dims[0] +
                                    k*region_dims[0]*region_dims[1];
                                
                                std::vector<const double*> m_ptr;
                                m_ptr.push_back(&rho_u[idx_data]);
                                m_ptr.push_back(&rho_v[idx_data]);
                                m_ptr.push_back(&rho_w[idx_data]);
                                
                                buffer[idx_region] = d_equation_of_state->getPressure(
                                    &rho[idx_data],
                                    m_ptr,
                                    &E[idx_data]);
                            }
                        }
                    }
                }
                
                data_on_patch = true;
                
                break;
            }
            case FOUR_EQN_SHYUE:
            {
                boost::shared_ptr<pdat::CellData<double> > density(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_density, d_plot_context)));
                
                boost::shared_ptr<pdat::CellData<double> > momentum(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_momentum, d_plot_context)));
                
                boost::shared_ptr<pdat::CellData<double> > total_energy(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_total_energy, d_plot_context)));
                
                boost::shared_ptr<pdat::CellData<double> > mass_fraction(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_mass_fraction, d_plot_context)));
                
#ifdef DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(density);
                TBOX_ASSERT(momentum);
                TBOX_ASSERT(total_energy);
                TBOX_ASSERT(mass_fraction);
                TBOX_ASSERT(density->getGhostBox().isSpatiallyEqual(patch.getBox()));
                TBOX_ASSERT(momentum->getGhostBox().isSpatiallyEqual(patch.getBox()));
                TBOX_ASSERT(total_energy->getGhostBox().isSpatiallyEqual(patch.getBox()));
                TBOX_ASSERT(mass_fraction->getGhostBox().isSpatiallyEqual(patch.getBox()));
#endif
                
                // Get the dimensions of box that covers the data.
                const hier::Box data_box = density->getGhostBox();
                const hier::IntVector data_dims = data_box.numberCells();
                
                // Get the arrays of time-dependent variables
                const double* const rho   = density->getPointer();
                const double* const rho_u = momentum->getPointer(0);
                const double* const rho_v = d_dim > tbox::Dimension(1) ? momentum->getPointer(1) : NULL;
                const double* const rho_w = d_dim > tbox::Dimension(2) ? momentum->getPointer(2) : NULL;
                const double* const E     = total_energy->getPointer(0);
                std::vector<const double*> Y;
                for (int si = 0; si < d_num_species - 1; si++)
                {
                    Y.push_back(mass_fraction->getPointer(si));
                }
                
                size_t offset_data = data_box.offset(region.lower());
                
                if (d_dim == tbox::Dimension(1))
                {
                    for (int i = 0; i < region_dims[0]; i++)
                    {
                        size_t idx_data = offset_data + i;
                        
                        size_t idx_region = i;
                        
                        std::vector<const double*> m_ptr;
                        m_ptr.push_back(&rho_u[idx_data]);
                        
                        std::vector<const double*> Y_ptr;
                        for (int si = 0; si < d_num_species - 1; si++)
                        {
                            Y_ptr.push_back(&Y[si][idx_data]);
                        }
                        
                        buffer[idx_region] = d_equation_of_state->getPressureWithMassFraction(
                            &rho[idx_data],
                            m_ptr,
                            &E[idx_data],
                            Y_ptr);
                    }
                }
                else if (d_dim == tbox::Dimension(2))
                {
                    for (int j = 0; j < region_dims[1]; j++)
                    {
                        for (int i = 0; i < region_dims[0]; i++)
                        {
                            size_t idx_data = offset_data + i +
                                j*data_dims[0];
                            
                            size_t idx_region = i +
                                j*region_dims[0];
                            
                            std::vector<const double*> m_ptr;
                            m_ptr.push_back(&rho_u[idx_data]);
                            m_ptr.push_back(&rho_v[idx_data]);
                            
                            std::vector<const double*> Y_ptr;
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                Y_ptr.push_back(&Y[si][idx_data]);
                            }
                            
                            buffer[idx_region] = d_equation_of_state->getPressureWithMassFraction(
                                &rho[idx_data],
                                m_ptr,
                                &E[idx_data],
                                Y_ptr);
                            }
                    }
                }
                else if (d_dim == tbox::Dimension(3))
                {
                    for (int k = 0; k < region_dims[2]; k++)
                    {
                        for (int j = 0; j < region_dims[1]; j++)
                        {
                            for (int i = 0; i < region_dims[0]; i++)
                            {
                                size_t idx_data = offset_data + i +
                                    j*data_dims[0] +
                                    k*data_dims[0]*data_dims[1];
                                
                                size_t idx_region = i +
                                    j*region_dims[0] +
                                    k*region_dims[0]*region_dims[1];
                                
                                std::vector<const double*> m_ptr;
                                m_ptr.push_back(&rho_u[idx_data]);
                                m_ptr.push_back(&rho_v[idx_data]);
                                m_ptr.push_back(&rho_w[idx_data]);
                                
                                std::vector<const double*> Y_ptr;
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    Y_ptr.push_back(&Y[si][idx_data]);
                                }
                                
                                buffer[idx_region] = d_equation_of_state->getPressureWithMassFraction(
                                    &rho[idx_data],
                                    m_ptr,
                                    &E[idx_data],
                                    Y_ptr);
                            }
                        }
                    }
                }
                
                data_on_patch = true;
                
                break;
            }
            case FIVE_EQN_ALLAIRE:
            {
                boost::shared_ptr<pdat::CellData<double> > partial_density(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_partial_density, d_plot_context)));
                
                boost::shared_ptr<pdat::CellData<double> > momentum(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_momentum, d_plot_context)));
                
                boost::shared_ptr<pdat::CellData<double> > total_energy(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_total_energy, d_plot_context)));
                
                boost::shared_ptr<pdat::CellData<double> > volume_fraction(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_volume_fraction, d_plot_context)));
                
#ifdef DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(partial_density);
                TBOX_ASSERT(momentum);
                TBOX_ASSERT(total_energy);
                TBOX_ASSERT(volume_fraction);
                TBOX_ASSERT(partial_density->getGhostBox().isSpatiallyEqual(patch.getBox()));
                TBOX_ASSERT(momentum->getGhostBox().isSpatiallyEqual(patch.getBox()));
                TBOX_ASSERT(total_energy->getGhostBox().isSpatiallyEqual(patch.getBox()));
                TBOX_ASSERT(volume_fraction->getGhostBox().isSpatiallyEqual(patch.getBox()));
#endif
                
                // Get the dimensions of box that covers the data.
                const hier::Box data_box = partial_density->getGhostBox();
                const hier::IntVector data_dims = data_box.numberCells();
                
                // Get the arrays of time-dependent variables
                std::vector<const double*> Z_rho;
                for (int si = 0; si < d_num_species; si++)
                {
                    Z_rho.push_back(partial_density->getPointer(si));
                }
                const double* const rho_u = momentum->getPointer(0);
                const double* const rho_v = d_dim > tbox::Dimension(1) ? momentum->getPointer(1) : NULL;
                const double* const rho_w = d_dim > tbox::Dimension(2) ? momentum->getPointer(2) : NULL;
                const double* const E     = total_energy->getPointer(0);
                std::vector<const double*> Z;
                for (int si = 0; si < d_num_species - 1; si++)
                {
                    Z.push_back(volume_fraction->getPointer(si));
                }
                
                size_t offset_data = data_box.offset(region.lower());
                
                if (d_dim == tbox::Dimension(1))
                {
                    for (int i = 0; i < region_dims[0]; i++)
                    {
                        size_t idx_data = offset_data + i;
                        
                        size_t idx_region = i;
                        
                        std::vector<const double*> Z_rho_ptr;
                        for (int si = 0; si < d_num_species; si++)
                        {
                            Z_rho_ptr.push_back(&Z_rho[si][idx_data]);
                        }
                        
                        std::vector<const double*> m_ptr;
                        m_ptr.push_back(&rho_u[idx_data]);
                        
                        std::vector<const double*> Z_ptr;
                        for (int si = 0; si < d_num_species - 1; si++)
                        {
                            Z_ptr.push_back(&Z[si][idx_data]);
                        }
                        
                        buffer[idx_region] = d_equation_of_state->getPressureWithVolumeFraction(
                            Z_rho_ptr,
                            m_ptr,
                            &E[idx_data],
                            Z_ptr);
                    }
                }
                else if (d_dim == tbox::Dimension(2))
                {
                    for (int j = 0; j < region_dims[1]; j++)
                    {
                        for (int i = 0; i < region_dims[0]; i++)
                        {
                            size_t idx_data = offset_data + i +
                                j*data_dims[0];
                            
                            size_t idx_region = i +
                                j*region_dims[0];
                            
                            std::vector<const double*> Z_rho_ptr;
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Z_rho_ptr.push_back(&Z_rho[si][idx_data]);
                            }
                            
                            std::vector<const double*> m_ptr;
                            m_ptr.push_back(&rho_u[idx_data]);
                            m_ptr.push_back(&rho_v[idx_data]);
                            
                            std::vector<const double*> Z_ptr;
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                Z_ptr.push_back(&Z[si][idx_data]);
                            }
                            
                            buffer[idx_region] = d_equation_of_state->getPressureWithVolumeFraction(
                                Z_rho_ptr,
                                m_ptr,
                                &E[idx_data],
                                Z_ptr);
                        }
                    }
                }
                else if (d_dim == tbox::Dimension(3))
                {
                    for (int k = 0; k < region_dims[2]; k++)
                    {
                        for (int j = 0; j < region_dims[1]; j++)
                        {
                            for (int i = 0; i < region_dims[0]; i++)
                            {
                                size_t idx_data = offset_data + i +
                                    j*data_dims[0] +
                                    k*data_dims[0]*data_dims[1];
                                
                                size_t idx_region = i +
                                    j*region_dims[0] +
                                    k*region_dims[0]*region_dims[1];
                                
                                std::vector<const double*> Z_rho_ptr;
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Z_rho_ptr.push_back(&Z_rho[si][idx_data]);
                                }
                                
                                std::vector<const double*> m_ptr;
                                m_ptr.push_back(&rho_u[idx_data]);
                                m_ptr.push_back(&rho_v[idx_data]);
                                m_ptr.push_back(&rho_w[idx_data]);
                                
                                std::vector<const double*> Z_ptr;
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    Z_ptr.push_back(&Z[si][idx_data]);
                                }
                                
                                buffer[idx_region] = d_equation_of_state->getPressureWithVolumeFraction(
                                    Z_rho_ptr,
                                    m_ptr,
                                    &E[idx_data],
                                    Z_ptr);
                            }
                        }
                    }
                }
                
                data_on_patch = true;
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                           << ": "
                           << "Unknown d_flow_model."
                           << std::endl);
            }
        }
    }
    else if (variable_name == "sound speed")
    {
        switch (d_flow_model)
        {
            case SINGLE_SPECIES:
            {
                boost::shared_ptr<pdat::CellData<double> > density(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_density, d_plot_context)));
                
                boost::shared_ptr<pdat::CellData<double> > momentum(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_momentum, d_plot_context)));
                
                boost::shared_ptr<pdat::CellData<double> > total_energy(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_total_energy, d_plot_context)));
                
#ifdef DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(density);
                TBOX_ASSERT(momentum);
                TBOX_ASSERT(total_energy);
                TBOX_ASSERT(density->getGhostBox().isSpatiallyEqual(patch.getBox()));
                TBOX_ASSERT(momentum->getGhostBox().isSpatiallyEqual(patch.getBox()));
                TBOX_ASSERT(total_energy->getGhostBox().isSpatiallyEqual(patch.getBox()));
#endif
                
                // Get the dimensions of box that covers the data.
                const hier::Box data_box = density->getGhostBox();
                const hier::IntVector data_dims = data_box.numberCells();
                
                // Get the arrays of time-dependent variables
                const double* const rho   = density->getPointer();
                const double* const rho_u = momentum->getPointer(0);
                const double* const rho_v = d_dim > tbox::Dimension(1) ? momentum->getPointer(1) : NULL;
                const double* const rho_w = d_dim > tbox::Dimension(2) ? momentum->getPointer(2) : NULL;
                const double* const E     = total_energy->getPointer(0);
                
                size_t offset_data = data_box.offset(region.lower());
                
                if (d_dim == tbox::Dimension(1))
                {
                    for (int i = 0; i < region_dims[0]; i++)
                    {
                        size_t idx_data = offset_data + i;
                        
                        size_t idx_region = i;
                        
                        std::vector<const double*> m_ptr;
                        m_ptr.push_back(&rho_u[idx_data]);
                        
                        buffer[idx_region] = d_equation_of_state->getSoundSpeed(
                            &rho[idx_data],
                            m_ptr,
                            &E[idx_data]);
                    }
                }
                else if (d_dim == tbox::Dimension(2))
                {
                    for (int j = 0; j < region_dims[1]; j++)
                    {
                        for (int i = 0; i < region_dims[0]; i++)
                        {
                            size_t idx_data = offset_data + i +
                                j*data_dims[0];
                            
                            size_t idx_region = i +
                                j*region_dims[0];
                            
                            std::vector<const double*> m_ptr;
                            m_ptr.push_back(&rho_u[idx_data]);
                            m_ptr.push_back(&rho_v[idx_data]);
                            
                            buffer[idx_region] = d_equation_of_state->getSoundSpeed(
                                &rho[idx_data],
                                m_ptr,
                                &E[idx_data]);
                        }
                    }
                }
                else if (d_dim == tbox::Dimension(3))
                {
                    for (int k = 0; k < region_dims[2]; k++)
                    {
                        for (int j = 0; j < region_dims[1]; j++)
                        {
                            for (int i = 0; i < region_dims[0]; i++)
                            {
                                size_t idx_data = offset_data + i +
                                    j*data_dims[0] +
                                    k*data_dims[0]*data_dims[1];
                                
                                size_t idx_region = i +
                                    j*region_dims[0] +
                                    k*region_dims[0]*region_dims[1];
                                
                                std::vector<const double*> m_ptr;
                                m_ptr.push_back(&rho_u[idx_data]);
                                m_ptr.push_back(&rho_v[idx_data]);
                                m_ptr.push_back(&rho_w[idx_data]);
                                
                                buffer[idx_region] = d_equation_of_state->getSoundSpeed(
                                    &rho[idx_data],
                                    m_ptr,
                                    &E[idx_data]);
                            }
                        }
                    }
                }
                
                data_on_patch = true;
                
                break;
            }
            case FOUR_EQN_SHYUE:
            {
                boost::shared_ptr<pdat::CellData<double> > density(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_density, d_plot_context)));
                
                boost::shared_ptr<pdat::CellData<double> > momentum(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_momentum, d_plot_context)));
                
                boost::shared_ptr<pdat::CellData<double> > total_energy(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_total_energy, d_plot_context)));
                
                boost::shared_ptr<pdat::CellData<double> > mass_fraction(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_mass_fraction, d_plot_context)));
                
#ifdef DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(density);
                TBOX_ASSERT(momentum);
                TBOX_ASSERT(total_energy);
                TBOX_ASSERT(mass_fraction);
                TBOX_ASSERT(density->getGhostBox().isSpatiallyEqual(patch.getBox()));
                TBOX_ASSERT(momentum->getGhostBox().isSpatiallyEqual(patch.getBox()));
                TBOX_ASSERT(total_energy->getGhostBox().isSpatiallyEqual(patch.getBox()));
                TBOX_ASSERT(mass_fraction->getGhostBox().isSpatiallyEqual(patch.getBox()));
#endif
                
                // Get the dimensions of box that covers the data.
                const hier::Box data_box = density->getGhostBox();
                const hier::IntVector data_dims = data_box.numberCells();
                
                // Get the arrays of time-dependent variables
                const double* const rho   = density->getPointer();
                const double* const rho_u = momentum->getPointer(0);
                const double* const rho_v = d_dim > tbox::Dimension(1) ? momentum->getPointer(1) : NULL;
                const double* const rho_w = d_dim > tbox::Dimension(2) ? momentum->getPointer(2) : NULL;
                const double* const E     = total_energy->getPointer(0);
                std::vector<const double*> Y;
                for (int si = 0; si < d_num_species - 1; si++)
                {
                    Y.push_back(mass_fraction->getPointer(si));
                }
                
                size_t offset_data = data_box.offset(region.lower());
                
                if (d_dim == tbox::Dimension(1))
                {
                    for (int i = 0; i < region_dims[0]; i++)
                    {
                        size_t idx_data = offset_data + i;
                        
                        size_t idx_region = i;
                        
                        std::vector<const double*> m_ptr;
                        m_ptr.push_back(&rho_u[idx_data]);
                        
                        std::vector<const double*> Y_ptr;
                        for (int si = 0; si < d_num_species - 1; si++)
                        {
                            Y_ptr.push_back(&Y[si][idx_data]);
                        }
                        
                        buffer[idx_region] = d_equation_of_state->getSoundSpeedWithMassFraction(
                            &rho[idx_data],
                            m_ptr,
                            &E[idx_data],
                            Y_ptr);
                    }
                }
                else if (d_dim == tbox::Dimension(2))
                {
                    for (int j = 0; j < region_dims[1]; j++)
                    {
                        for (int i = 0; i < region_dims[0]; i++)
                        {
                            size_t idx_data = offset_data + i +
                                j*data_dims[0];
                            
                            size_t idx_region = i +
                                j*region_dims[0];
                            
                            std::vector<const double*> m_ptr;
                            m_ptr.push_back(&rho_u[idx_data]);
                            m_ptr.push_back(&rho_v[idx_data]);
                            
                            std::vector<const double*> Y_ptr;
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                Y_ptr.push_back(&Y[si][idx_data]);
                            }
                            
                            buffer[idx_region] = d_equation_of_state->getSoundSpeedWithMassFraction(
                                &rho[idx_data],
                                m_ptr,
                                &E[idx_data],
                                Y_ptr);
                            }
                    }
                }
                else if (d_dim == tbox::Dimension(3))
                {
                    for (int k = 0; k < region_dims[2]; k++)
                    {
                        for (int j = 0; j < region_dims[1]; j++)
                        {
                            for (int i = 0; i < region_dims[0]; i++)
                            {
                                size_t idx_data = offset_data + i +
                                    j*data_dims[0] +
                                    k*data_dims[0]*data_dims[1];
                                
                                size_t idx_region = i +
                                    j*region_dims[0] +
                                    k*region_dims[0]*region_dims[1];
                                
                                std::vector<const double*> m_ptr;
                                m_ptr.push_back(&rho_u[idx_data]);
                                m_ptr.push_back(&rho_v[idx_data]);
                                m_ptr.push_back(&rho_w[idx_data]);
                                
                                std::vector<const double*> Y_ptr;
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    Y_ptr.push_back(&Y[si][idx_data]);
                                }
                                
                                buffer[idx_region] = d_equation_of_state->getSoundSpeedWithMassFraction(
                                    &rho[idx_data],
                                    m_ptr,
                                    &E[idx_data],
                                    Y_ptr);
                            }
                        }
                    }
                }
                
                data_on_patch = true;
                
                break;
            }
            case FIVE_EQN_ALLAIRE:
            {
                boost::shared_ptr<pdat::CellData<double> > partial_density(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_partial_density, d_plot_context)));
                
                boost::shared_ptr<pdat::CellData<double> > momentum(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_momentum, d_plot_context)));
                
                boost::shared_ptr<pdat::CellData<double> > total_energy(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_total_energy, d_plot_context)));
                
                boost::shared_ptr<pdat::CellData<double> > volume_fraction(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_volume_fraction, d_plot_context)));
                
#ifdef DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(partial_density);
                TBOX_ASSERT(momentum);
                TBOX_ASSERT(total_energy);
                TBOX_ASSERT(volume_fraction);
                TBOX_ASSERT(partial_density->getGhostBox().isSpatiallyEqual(patch.getBox()));
                TBOX_ASSERT(momentum->getGhostBox().isSpatiallyEqual(patch.getBox()));
                TBOX_ASSERT(total_energy->getGhostBox().isSpatiallyEqual(patch.getBox()));
                TBOX_ASSERT(volume_fraction->getGhostBox().isSpatiallyEqual(patch.getBox()));
#endif
                
                // Get the dimensions of box that covers the data.
                const hier::Box data_box = partial_density->getGhostBox();
                const hier::IntVector data_dims = data_box.numberCells();
                
                // Get the arrays of time-dependent variables
                std::vector<const double*> Z_rho;
                for (int si = 0; si < d_num_species; si++)
                {
                    Z_rho.push_back(partial_density->getPointer(si));
                }
                const double* const rho_u = momentum->getPointer(0);
                const double* const rho_v = d_dim > tbox::Dimension(1) ? momentum->getPointer(1) : NULL;
                const double* const rho_w = d_dim > tbox::Dimension(2) ? momentum->getPointer(2) : NULL;
                const double* const E     = total_energy->getPointer(0);
                std::vector<const double*> Z;
                for (int si = 0; si < d_num_species - 1; si++)
                {
                    Z.push_back(volume_fraction->getPointer(si));
                }
                
                size_t offset_data = data_box.offset(region.lower());
                
                if (d_dim == tbox::Dimension(1))
                {
                    for (int i = 0; i < region_dims[0]; i++)
                    {
                        size_t idx_data = offset_data + i;
                        
                        size_t idx_region = i;
                        
                        std::vector<const double*> Z_rho_ptr;
                        for (int si = 0; si < d_num_species; si++)
                        {
                            Z_rho_ptr.push_back(&Z_rho[si][idx_data]);
                        }
                        
                        std::vector<const double*> m_ptr;
                        m_ptr.push_back(&rho_u[idx_data]);
                        
                        std::vector<const double*> Z_ptr;
                        for (int si = 0; si < d_num_species - 1; si++)
                        {
                            Z_ptr.push_back(&Z[si][idx_data]);
                        }
                        
                        buffer[idx_region] = d_equation_of_state->getSoundSpeedWithVolumeFraction(
                            Z_rho_ptr,
                            m_ptr,
                            &E[idx_data],
                            Z_ptr);
                    }
                }
                else if (d_dim == tbox::Dimension(2))
                {
                    for (int j = 0; j < region_dims[1]; j++)
                    {
                        for (int i = 0; i < region_dims[0]; i++)
                        {
                            size_t idx_data = offset_data + i +
                                j*data_dims[0];
                            
                            size_t idx_region = i +
                                j*region_dims[0];
                            
                            std::vector<const double*> Z_rho_ptr;
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Z_rho_ptr.push_back(&Z_rho[si][idx_data]);
                            }
                            
                            std::vector<const double*> m_ptr;
                            m_ptr.push_back(&rho_u[idx_data]);
                            m_ptr.push_back(&rho_v[idx_data]);
                            
                            std::vector<const double*> Z_ptr;
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                Z_ptr.push_back(&Z[si][idx_data]);
                            }
                            
                            buffer[idx_region] = d_equation_of_state->getSoundSpeedWithVolumeFraction(
                                Z_rho_ptr,
                                m_ptr,
                                &E[idx_data],
                                Z_ptr);
                        }
                    }
                }
                else if (d_dim == tbox::Dimension(3))
                {
                    for (int k = 0; k < region_dims[2]; k++)
                    {
                        for (int j = 0; j < region_dims[1]; j++)
                        {
                            for (int i = 0; i < region_dims[0]; i++)
                            {
                                size_t idx_data = offset_data + i +
                                    j*data_dims[0] +
                                    k*data_dims[0]*data_dims[1];
                                
                                size_t idx_region = i +
                                    j*region_dims[0] +
                                    k*region_dims[0]*region_dims[1];
                                
                                std::vector<const double*> Z_rho_ptr;
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Z_rho_ptr.push_back(&Z_rho[si][idx_data]);
                                }
                                
                                std::vector<const double*> m_ptr;
                                m_ptr.push_back(&rho_u[idx_data]);
                                m_ptr.push_back(&rho_v[idx_data]);
                                m_ptr.push_back(&rho_w[idx_data]);
                                
                                std::vector<const double*> Z_ptr;
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    Z_ptr.push_back(&Z[si][idx_data]);
                                }
                                
                                buffer[idx_region] = d_equation_of_state->getSoundSpeedWithVolumeFraction(
                                    Z_rho_ptr,
                                    m_ptr,
                                    &E[idx_data],
                                    Z_ptr);
                            }
                        }
                    }
                }
                
                data_on_patch = true;
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                           << ": "
                           << "Unknown d_flow_model."
                           << std::endl);
            }
        }
    }
    else if (variable_name == "velocity")
    {
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(depth_id < d_dim.getValue());
#endif
        switch (d_flow_model)
        {
            case SINGLE_SPECIES:
            {
                boost::shared_ptr<pdat::CellData<double> > density(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_density, d_plot_context)));
                
                boost::shared_ptr<pdat::CellData<double> > momentum(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_momentum, d_plot_context)));
                
#ifdef DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(density);
                TBOX_ASSERT(momentum);
                TBOX_ASSERT(density->getGhostBox().isSpatiallyEqual(patch.getBox()));
                TBOX_ASSERT(momentum->getGhostBox().isSpatiallyEqual(patch.getBox()));
#endif
                
                // Get the dimensions of box that covers the data.
                const hier::Box data_box = density->getGhostBox();
                const hier::IntVector data_dims = data_box.numberCells();
                
                // Get the arrays of time-dependent variables
                const double* const rho   = density->getPointer();
                const double* const m     = momentum->getPointer(depth_id);
                
                size_t offset_data = data_box.offset(region.lower());
                
                if (d_dim == tbox::Dimension(1))
                {
                    for (int i = 0; i < region_dims[0]; i++)
                    {
                        size_t idx_data = offset_data + i;
                        
                        size_t idx_region = i;
                        
                        buffer[idx_region] = m[idx_data]/rho[idx_data];
                    }
                }
                else if (d_dim == tbox::Dimension(2))
                {
                    for (int j = 0; j < region_dims[1]; j++)
                    {
                        for (int i = 0; i < region_dims[0]; i++)
                        {
                            size_t idx_data = offset_data + i +
                                j*data_dims[0];
                            
                            size_t idx_region = i +
                                j*region_dims[0];
                            
                            buffer[idx_region] = m[idx_data]/rho[idx_data];
                        }
                    }
                }
                else if (d_dim == tbox::Dimension(3))
                {
                    for (int k = 0; k < region_dims[2]; k++)
                    {
                        for (int j = 0; j < region_dims[1]; j++)
                        {
                            for (int i = 0; i < region_dims[0]; i++)
                            {
                                size_t idx_data = offset_data + i +
                                    j*data_dims[0] +
                                    k*data_dims[0]*data_dims[1];
                                
                                size_t idx_region = i +
                                    j*region_dims[0] +
                                    k*region_dims[0]*region_dims[1];
                                
                                buffer[idx_region] = m[idx_data]/rho[idx_data];
                            }
                        }
                    }
                }
                
                data_on_patch = true;
                
                break;
            }
            case FOUR_EQN_SHYUE:
            {
                boost::shared_ptr<pdat::CellData<double> > density(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_density, d_plot_context)));
                
                boost::shared_ptr<pdat::CellData<double> > momentum(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_momentum, d_plot_context)));
                
#ifdef DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(density);
                TBOX_ASSERT(momentum);
                TBOX_ASSERT(density->getGhostBox().isSpatiallyEqual(patch.getBox()));
                TBOX_ASSERT(momentum->getGhostBox().isSpatiallyEqual(patch.getBox()));
#endif
                
                // Get the dimensions of box that covers the data.
                const hier::Box data_box = density->getGhostBox();
                const hier::IntVector data_dims = data_box.numberCells();
                
                // Get the arrays of time-dependent variables
                const double* const rho   = density->getPointer();
                const double* const m     = momentum->getPointer(depth_id);
                
                size_t offset_data = data_box.offset(region.lower());
                
                if (d_dim == tbox::Dimension(1))
                {
                    for (int i = 0; i < region_dims[0]; i++)
                    {
                        size_t idx_data = offset_data + i;
                        
                        size_t idx_region = i;
                        
                        buffer[idx_region] = m[idx_data]/rho[idx_data];
                    }
                }
                else if (d_dim == tbox::Dimension(2))
                {
                    for (int j = 0; j < region_dims[1]; j++)
                    {
                        for (int i = 0; i < region_dims[0]; i++)
                        {
                            size_t idx_data = offset_data + i +
                                j*data_dims[0];
                            
                            size_t idx_region = i +
                                j*region_dims[0];
                            
                            buffer[idx_region] = m[idx_data]/rho[idx_data];
                        }
                    }
                }
                else if (d_dim == tbox::Dimension(3))
                {
                    for (int k = 0; k < region_dims[2]; k++)
                    {
                        for (int j = 0; j < region_dims[1]; j++)
                        {
                            for (int i = 0; i < region_dims[0]; i++)
                            {
                                size_t idx_data = offset_data + i +
                                    j*data_dims[0] +
                                    k*data_dims[0]*data_dims[1];
                                
                                size_t idx_region = i +
                                    j*region_dims[0] +
                                    k*region_dims[0]*region_dims[1];
                                
                                buffer[idx_region] = m[idx_data]/rho[idx_data];
                            }
                        }
                    }
                }
                
                data_on_patch = true;
                
                break;
            }
            case FIVE_EQN_ALLAIRE:
            {
                boost::shared_ptr<pdat::CellData<double> > partial_density(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_partial_density, d_plot_context)));
                
                boost::shared_ptr<pdat::CellData<double> > momentum(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_momentum, d_plot_context)));
                
#ifdef DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(partial_density);
                TBOX_ASSERT(momentum);
                TBOX_ASSERT(partial_density->getGhostBox().isSpatiallyEqual(patch.getBox()));
                TBOX_ASSERT(momentum->getGhostBox().isSpatiallyEqual(patch.getBox()));
#endif
                
                // Get the dimensions of box that covers the data.
                const hier::Box data_box = partial_density->getGhostBox();
                const hier::IntVector data_dims = data_box.numberCells();
                
                // Get the arrays of time-dependent variables
                std::vector<const double*> Z_rho;
                for (int si = 0; si < d_num_species; si++)
                {
                    Z_rho.push_back(partial_density->getPointer(si));
                }
                const double* const m = momentum->getPointer(depth_id);
                
                size_t offset_data = data_box.offset(region.lower());
                
                if (d_dim == tbox::Dimension(1))
                {
                    for (int i = 0; i < region_dims[0]; i++)
                    {
                        size_t idx_data = offset_data + i;
                        
                        size_t idx_region = i;
                        
                        std::vector<const double*> Z_rho_ptr;
                        for (int si = 0; si < d_num_species; si++)
                        {
                            Z_rho_ptr.push_back(&Z_rho[si][idx_data]);
                        }
                        
                        double rho = d_equation_of_state->getTotalDensity(Z_rho_ptr);
                        
                        buffer[idx_region] = m[idx_data]/rho;
                    }
                }
                else if (d_dim == tbox::Dimension(2))
                {
                    for (int j = 0; j < region_dims[1]; j++)
                    {
                        for (int i = 0; i < region_dims[0]; i++)
                        {
                            size_t idx_data = offset_data + i +
                                j*data_dims[0];
                            
                            size_t idx_region = i +
                                j*region_dims[0];
                            
                            std::vector<const double*> Z_rho_ptr;
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Z_rho_ptr.push_back(&Z_rho[si][idx_data]);
                            }
                            
                            double rho = d_equation_of_state->getTotalDensity(Z_rho_ptr);
                            
                            buffer[idx_region] = m[idx_data]/rho;
                        }
                    }
                }
                else if (d_dim == tbox::Dimension(3))
                {
                    for (int k = 0; k < region_dims[2]; k++)
                    {
                        for (int j = 0; j < region_dims[1]; j++)
                        {
                            for (int i = 0; i < region_dims[0]; i++)
                            {
                                size_t idx_data = offset_data + i +
                                    j*data_dims[0] +
                                    k*data_dims[0]*data_dims[1];
                                
                                size_t idx_region = i +
                                    j*region_dims[0] +
                                    k*region_dims[0]*region_dims[1];
                                
                                std::vector<const double*> Z_rho_ptr;
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Z_rho_ptr.push_back(&Z_rho[si][idx_data]);
                                }
                                
                                double rho = d_equation_of_state->getTotalDensity(Z_rho_ptr);
                                
                                buffer[idx_region] = m[idx_data]/rho;
                            }
                        }
                    }
                }
                
                data_on_patch = true;
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                           << ": "
                           << "Unknown d_flow_model."
                           << std::endl);
            }
        }
    }
    else if (variable_name == "density")
    {
        switch (d_flow_model)
        {
            case SINGLE_SPECIES:
            {
                TBOX_ERROR("Euler::packDerivedDataIntoDoubleBuffer()"
                   << "\n    'Density' is already registered."
                   << std::endl);
                
                data_on_patch = false;
                
                break;
            }
            case FOUR_EQN_SHYUE:
            {
                TBOX_ERROR("Euler::packDerivedDataIntoDoubleBuffer()"
                   << "\n    'Density' is already registered."
                   << std::endl);
                
                data_on_patch = false;
                
                break;
            }
            case FIVE_EQN_ALLAIRE:
            {
                boost::shared_ptr<pdat::CellData<double> > partial_density(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_partial_density, d_plot_context)));
                
#ifdef DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(partial_density);
                TBOX_ASSERT(partial_density->getGhostBox().isSpatiallyEqual(patch.getBox()));
#endif
                
                // Get the dimensions of box that covers the data.
                const hier::Box data_box = partial_density->getGhostBox();
                const hier::IntVector data_dims = data_box.numberCells();
                
                // Get the arrays of time-dependent variables
                std::vector<const double*> Z_rho;
                for (int si = 0; si < d_num_species; si++)
                {
                    Z_rho.push_back(partial_density->getPointer(si));
                }
                
                size_t offset_data = data_box.offset(region.lower());
                
                if (d_dim == tbox::Dimension(1))
                {
                    for (int i = 0; i < region_dims[0]; i++)
                    {
                        size_t idx_data = offset_data + i;
                        
                        size_t idx_region = i;
                        
                        std::vector<const double*> Z_rho_ptr;
                        for (int si = 0; si < d_num_species; si++)
                        {
                            Z_rho_ptr.push_back(&Z_rho[si][idx_data]);
                        }
                        
                        buffer[idx_region] = d_equation_of_state->getTotalDensity(Z_rho_ptr);
                    }
                }
                else if (d_dim == tbox::Dimension(2))
                {
                    for (int j = 0; j < region_dims[1]; j++)
                    {
                        for (int i = 0; i < region_dims[0]; i++)
                        {
                            size_t idx_data = offset_data + i +
                                j*data_dims[0];
                            
                            size_t idx_region = i +
                                j*region_dims[0];
                            
                            std::vector<const double*> Z_rho_ptr;
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Z_rho_ptr.push_back(&Z_rho[si][idx_data]);
                            }
                            
                            buffer[idx_region] = d_equation_of_state->getTotalDensity(Z_rho_ptr);
                        }
                    }
                }
                else if (d_dim == tbox::Dimension(3))
                {
                    for (int k = 0; k < region_dims[2]; k++)
                    {
                        for (int j = 0; j < region_dims[1]; j++)
                        {
                            for (int i = 0; i < region_dims[0]; i++)
                            {
                                size_t idx_data = offset_data + i +
                                    j*data_dims[0] +
                                    k*data_dims[0]*data_dims[1];
                                
                                size_t idx_region = i +
                                    j*region_dims[0] +
                                    k*region_dims[0]*region_dims[1];
                                
                                std::vector<const double*> Z_rho_ptr;
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Z_rho_ptr.push_back(&Z_rho[si][idx_data]);
                                }
                                
                                buffer[idx_region] = d_equation_of_state->getTotalDensity(Z_rho_ptr);
                            }
                        }
                    }
                }
                
                data_on_patch = true;
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                           << ": "
                           << "Unknown d_flow_model."
                           << std::endl);
            }
        }
    }
    else if (variable_name.find("mass fraction") != std::string::npos)
    {
        switch (d_flow_model)
        {
            case SINGLE_SPECIES:
            {
                TBOX_ERROR("Euler::packDerivedDataIntoDoubleBuffer()"
                   << "\n    'Mass fraction' of single-species cannot be registered."
                   << std::endl);
                
                data_on_patch = false;
                
                break;
            }
            case FOUR_EQN_SHYUE:
            {
                TBOX_ERROR("Euler::packDerivedDataIntoDoubleBuffer()"
                   << "\n    'Mass fraction' is already registered."
                   << std::endl);
                
                data_on_patch = false;
                
                break;
            }
            case FIVE_EQN_ALLAIRE:
            {
                int species_idx = std::stoi(variable_name.substr(14));
                
                boost::shared_ptr<pdat::CellData<double> > partial_density(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_partial_density, d_plot_context)));
                
#ifdef DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(partial_density);
                TBOX_ASSERT(partial_density->getGhostBox().isSpatiallyEqual(patch.getBox()));
#endif
                
                // Get the dimensions of box that covers the data.
                const hier::Box data_box = partial_density->getGhostBox();
                const hier::IntVector data_dims = data_box.numberCells();
                
                // Get the arrays of time-dependent variables
                std::vector<const double*> Z_rho;
                for (int si = 0; si < d_num_species; si++)
                {
                    Z_rho.push_back(partial_density->getPointer(si));
                }
                
                size_t offset_data = data_box.offset(region.lower());
                
                if (d_dim == tbox::Dimension(1))
                {
                    for (int i = 0; i < region_dims[0]; i++)
                    {
                        size_t idx_data = offset_data + i;
                        
                        size_t idx_region = i;
                        
                        std::vector<const double*> Z_rho_ptr;
                        for (int si = 0; si < d_num_species; si++)
                        {
                            Z_rho_ptr.push_back(&Z_rho[si][idx_data]);
                        }
                        
                        double rho = d_equation_of_state->getTotalDensity(Z_rho_ptr);
                        
                        buffer[idx_region] = Z_rho[species_idx][idx_data]/rho;
                    }
                }
                else if (d_dim == tbox::Dimension(2))
                {
                    for (int j = 0; j < region_dims[1]; j++)
                    {
                        for (int i = 0; i < region_dims[0]; i++)
                        {
                            size_t idx_data = offset_data + i +
                                j*data_dims[0];
                            
                            size_t idx_region = i +
                                j*region_dims[0];
                            
                            std::vector<const double*> Z_rho_ptr;
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Z_rho_ptr.push_back(&Z_rho[si][idx_data]);
                            }
                            
                            double rho = d_equation_of_state->getTotalDensity(Z_rho_ptr);
                            
                            buffer[idx_region] = Z_rho[species_idx][idx_data]/rho;
                        }
                    }
                }
                else if (d_dim == tbox::Dimension(3))
                {
                    for (int k = 0; k < region_dims[2]; k++)
                    {
                        for (int j = 0; j < region_dims[1]; j++)
                        {
                            for (int i = 0; i < region_dims[0]; i++)
                            {
                                size_t idx_data = offset_data + i +
                                    j*data_dims[0] +
                                    k*data_dims[0]*data_dims[1];
                                
                                size_t idx_region = i +
                                    j*region_dims[0] +
                                    k*region_dims[0]*region_dims[1];
                                
                                std::vector<const double*> Z_rho_ptr;
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Z_rho_ptr.push_back(&Z_rho[si][idx_data]);
                                }
                                
                                double rho = d_equation_of_state->getTotalDensity(Z_rho_ptr);
                                
                                buffer[idx_region] = Z_rho[species_idx][idx_data]/rho;
                            }
                        }
                    }
                }
                
                data_on_patch = true;
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                           << ": "
                           << "Unknown d_flow_model."
                           << std::endl);
            }
        }
    }
    else
    {
        TBOX_ERROR("Euler::packDerivedDataIntoDoubleBuffer()"
                   << "\n    unknown variable_name "
                   << variable_name
                   << std::endl);
    }
    
    return data_on_patch;
}


void
Euler::boundaryReset(
    hier::Patch& patch,
    pdat::FaceData<double>& traced_left,
    pdat::FaceData<double>& traced_right) const
{
    
}


void Euler::printClassData(std::ostream& os) const
{
    os << "\nEuler::printClassData..." << std::endl;
    
    os << std::endl;
    
    os << "Euler: this = " << (Euler *)this << std::endl;
    os << "d_object_name = " << d_object_name << std::endl;
    os << "d_project_name = " << d_project_name << std::endl;
    os << "d_dim = " << d_dim.getValue() << std::endl;
    os << "d_grid_geometry = " << d_grid_geometry.get() << std::endl;
    os << "d_num_ghosts = " << d_num_ghosts << std::endl;
    
    switch (d_flow_model)
    {
        case SINGLE_SPECIES:
        {
            os << "d_flow_model = SINGLE_SPECIES" << std::endl;
            break;
        }
        case FOUR_EQN_SHYUE:
        {
            os << "d_flow_model = FOUR_EQN_SHYUE" << std::endl;
            break;
        }
        case FIVE_EQN_ALLAIRE:
        {
            os << "d_flow_model = FIVE_EQN_ALLAIRE" << std::endl;
            break;
        }
        default:
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "Unknown d_flow_model."
                       << std::endl);
        }
    }
    
    os << "d_num_eqn = " << d_num_eqn << std::endl;
    os << "d_num_species = " << d_num_species << std::endl;
    
    os << "d_is_preserving_positivity = " << d_is_preserving_positivity << std::endl;
    
    /*
     * Print data of d_grid_geometry.
     */
    os << "\nGrid geometry data:" << std::endl;
    d_grid_geometry->printClassData(os);
    os << std::endl;
    os << "End of BaseGridGeometry::printClassData" << std::endl;
    
    /*
     * Print data of d_equation_of_state object
     */
    
    os << "\nEquation of state data:" << std::endl;
    
    d_equation_of_state->printClassData(os);
    
    /*
     * Print data of d_conv_flux_reconstructor.
     */
    
    os << "\nConvective flux reconstructor data:" << std::endl;
    
    d_conv_flux_reconstructor->printClassData(os);
    
    /*
     * Print Refinement data
     */
    
    /*
     * Print boundary condition data.
     */
    
    os << "\nBoundary condition data:" << std::endl;
    
    os << std::endl;
    
    if (d_dim == tbox::Dimension(1))
    {
        // NOT YET IMPLEMENTED
    }
    else if (d_dim == tbox::Dimension(2))
    {
        for (int j = 0; j < static_cast<int>(d_master_bdry_node_conds.size()); j++)
        {
             os << "d_master_bdry_node_conds["
                << j
                << "] = "
                << d_master_bdry_node_conds[j]
                << std::endl;
             os << "d_scalar_bdry_node_conds["
                << j
                << "] = "
                << d_scalar_bdry_node_conds[j]
                << std::endl;
             os << "d_vector_bdry_node_conds["
                << j
                << "] = "
                << d_vector_bdry_node_conds[j]
                << std::endl;
             os << "d_node_bdry_edge["
                << j
                << "] = "
                << d_node_bdry_edge[j]
                << std::endl;
        }
        
        os << std::endl;
        
        for (int j = 0; j < static_cast<int>(d_master_bdry_edge_conds.size()); j++)
        {
            os << "d_master_bdry_edge_conds["
               << j
               << "] = "
               << d_master_bdry_edge_conds[j]
               << std::endl;
            
            os << "d_scalar_bdry_edge_conds["
               << j
               << "] = "
               << d_scalar_bdry_edge_conds[j]
               << std::endl;
            
            os << "d_vector_bdry_edge_conds["
               << j
               << "] = "
               << d_vector_bdry_edge_conds[j]
               << std::endl;
            
            if (d_master_bdry_edge_conds[j] == BdryCond::DIRICHLET)
            {
                switch (d_flow_model)
                {
                    case SINGLE_SPECIES:
                    {
                        os << "d_bdry_edge_density["
                           << j
                           << "] = "
                           << d_bdry_edge_density[j]
                           << std::endl;
                        
                        os << "d_bdry_edge_momentum["
                           << j << "] = "
                           << d_bdry_edge_momentum[j*d_dim.getValue() + 0]
                           << " , "
                           << d_bdry_edge_momentum[j*d_dim.getValue() + 1]
                           << std::endl;
                        
                        os << "d_bdry_edge_total_energy["
                           << j
                           << "] = "
                           << d_bdry_edge_total_energy[j]
                           << std::endl;
                        
                        break;
                    }
                    case FOUR_EQN_SHYUE:
                    {
                        os << "d_bdry_edge_density["
                           << j
                           << "] = "
                           << d_bdry_edge_density[j]
                           << std::endl;
                        
                        os << "d_bdry_edge_momentum["
                           << j
                           << "] = "
                           << d_bdry_edge_momentum[j*d_dim.getValue() + 0]
                           << " , "
                           << d_bdry_edge_momentum[j*d_dim.getValue() + 1]
                           << std::endl;
                        
                        os << "d_bdry_edge_total_energy["
                           << j
                           << "] = "
                           << d_bdry_edge_total_energy[j]
                           << std::endl;
                        
                        os << "d_bdry_edge_mass_fraction["
                           << j
                           << "] = "
                           << d_bdry_edge_mass_fraction[j*d_num_species + 0];
                        for (int si = 1; si < d_num_species - 1; si++)
                        {
                            os << " , "
                               << d_bdry_edge_mass_fraction[j*d_num_species + si];
                        }
                        os << std::endl;
                        
                        break;
                    }
                    case FIVE_EQN_ALLAIRE:
                    {
                        os << "d_bdry_edge_partial_density["
                           << j
                           << "] = "
                           << d_bdry_edge_partial_density[j*d_num_species + 0];
                        for (int si = 1; si < d_num_species; si++)
                        {
                            os << " , "
                               << d_bdry_edge_partial_density[j*d_num_species + si];
                        }
                        os << std::endl;
                        
                        os << "d_bdry_edge_momentum["
                           << j << "] = "
                           << d_bdry_edge_momentum[j*d_dim.getValue() + 0]
                           << " , "
                           << d_bdry_edge_momentum[j*d_dim.getValue() + 1]
                           << std::endl;
                        
                        os << "d_bdry_edge_total_energy["
                           << j
                           << "] = "
                           << d_bdry_edge_total_energy[j]
                           << std::endl;
                        
                        os << "d_bdry_edge_volume_fraction["
                           << j
                           << "] = "
                           << d_bdry_edge_volume_fraction[j*d_num_species + 0];
                        for (int si = 1; si < d_num_species - 1; si++)
                        {
                            os << " , "
                               << d_bdry_edge_volume_fraction[j*d_num_species + si];
                        }
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                                   << ": "
                                   << "Unknown d_flow_model."
                                   << std::endl);
                    }
                }
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        for (int j = 0; j < static_cast<int>(d_master_bdry_node_conds.size()); j++)
        {
            os << "d_master_bdry_node_conds["
               << j
               << "] = "
               << d_master_bdry_node_conds[j]
               << std::endl;
            os << "d_scalar_bdry_node_conds["
               << j
               << "] = "
               << d_scalar_bdry_node_conds[j]
               << std::endl;
            os << "d_vector_bdry_node_conds["
               << j
               << "] = "
               << d_vector_bdry_node_conds[j]
               << std::endl;
            os << "d_node_bdry_face["
               << j
               << "] = "
               << d_node_bdry_face[j]
               << std::endl;
        }
        
        os << std::endl;
        
        for (int j = 0; j < static_cast<int>(d_master_bdry_edge_conds.size()); j++)
        {
            os << "d_master_bdry_edge_conds["
               << j
               << "] = "
               << d_master_bdry_edge_conds[j]
               << std::endl;
            os << "d_scalar_bdry_edge_conds["
               << j
               << "] = "
               << d_scalar_bdry_edge_conds[j]
               << std::endl;
            os << "d_vector_bdry_edge_conds["
               << j
               << "] = "
               << d_vector_bdry_edge_conds[j]
               << std::endl;
            os << "d_edge_bdry_face["
               << j
               << "] = "
               << d_edge_bdry_face[j]
               << std::endl;
        }
        
        os << std::endl;
        
        
        for (int j = 0; j < static_cast<int>(d_master_bdry_face_conds.size()); j++)
        {
            os << "d_master_bdry_face_conds["
               << j
               << "] = "
               << d_master_bdry_face_conds[j]
               << std::endl;
            os << "d_scalar_bdry_face_conds["
               << j
               << "] = "
               << d_scalar_bdry_face_conds[j]
               << std::endl;
            os << "d_vector_bdry_face_conds["
               << j
               << "] = "
               << d_vector_bdry_face_conds[j]
               << std::endl;
            
            if (d_master_bdry_face_conds[j] == BdryCond::DIRICHLET)
            {
                switch (d_flow_model)
                {
                    case SINGLE_SPECIES:
                    {
                        os << "d_bdry_face_density["
                           << j
                           << "] = "
                           << d_bdry_face_density[j]
                           << std::endl;
                        
                        os << "d_bdry_face_momentum["
                           << j << "] = "
                           << d_bdry_face_momentum[j*d_dim.getValue() + 0]
                           << " , "
                           << d_bdry_face_momentum[j*d_dim.getValue() + 1]
                           << " , "
                           << d_bdry_face_momentum[j*d_dim.getValue() + 2]
                           << std::endl;
                        
                        os << "d_bdry_face_total_energy["
                           << j
                           << "] = "
                           << d_bdry_face_total_energy[j]
                           << std::endl;
                        
                        break;
                    }
                    case FOUR_EQN_SHYUE:
                    {
                        os << "d_bdry_face_density["
                           << j
                           << "] = "
                           << d_bdry_face_density[j]
                           << std::endl;
                        
                        os << "d_bdry_face_momentum["
                           << j
                           << "] = "
                           << d_bdry_face_momentum[j*d_dim.getValue() + 0]
                           << " , "
                           << d_bdry_face_momentum[j*d_dim.getValue() + 1]
                           << " , "
                           << d_bdry_face_momentum[j*d_dim.getValue() + 2]
                           << std::endl;
                        
                        os << "d_bdry_face_total_energy["
                           << j
                           << "] = "
                           << d_bdry_face_total_energy[j]
                           << std::endl;
                        
                        os << "d_bdry_face_mass_fraction["
                           << j
                           << "] = "
                           << d_bdry_face_mass_fraction[j*d_num_species + 0];
                        for (int si = 1; si < d_num_species - 1; si++)
                        {
                            os << " , "
                               << d_bdry_face_mass_fraction[j*d_num_species + si];
                        }
                        os << std::endl;
                        
                        break;
                    }
                    case FIVE_EQN_ALLAIRE:
                    {
                        os << "d_bdry_face_partial_density["
                           << j
                           << "] = "
                           << d_bdry_face_partial_density[j*d_num_species + 0];
                        for (int si = 1; si < d_num_species; si++)
                        {
                            os << " , "
                               << d_bdry_face_partial_density[j*d_num_species + si];
                        }
                        os << std::endl;
                        
                        os << "d_bdry_face_momentum["
                           << j << "] = "
                           << d_bdry_face_momentum[j*d_dim.getValue() + 0]
                           << " , "
                           << d_bdry_face_momentum[j*d_dim.getValue() + 1]
                           << " , "
                           << d_bdry_face_momentum[j*d_dim.getValue() + 2]
                           << std::endl;
                        
                        os << "d_bdry_face_total_energy["
                           << j
                           << "] = "
                           << d_bdry_face_total_energy[j]
                           << std::endl;
                        
                        os << "d_bdry_face_volume_fraction["
                           << j
                           << "] = "
                           << d_bdry_face_volume_fraction[j*d_num_species + 0];
                        for (int si = 1; si < d_num_species - 1; si++)
                        {
                            os << " , "
                               << d_bdry_face_volume_fraction[j*d_num_species + si];
                        }
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                                   << ": "
                                   << "Unknown d_flow_model."
                                   << std::endl);
                    }
                }
            }
        }
    }
    
    os << std::endl;
    
    os << "End of Euler::printClassData" << std::endl;
}


void
Euler::printDataStatistics(std::ostream& os,
    const boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy) const
{
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    SAMRAI::math::HierarchyCellDataOpsReal<double> cell_double_operator(patch_hierarchy, 0, 0);
    
    hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();
    
    switch (d_flow_model)
    {
        case SINGLE_SPECIES:
        {
            const int rho_id = variable_db->mapVariableAndContextToIndex(
                d_density,
                d_plot_context);
            
            const int m_id = variable_db->mapVariableAndContextToIndex(
                d_momentum,
                d_plot_context);
            
            const int E_id = variable_db->mapVariableAndContextToIndex(
                d_total_energy,
                d_plot_context);
            
            double rho_max_local = cell_double_operator.max(rho_id);
            double rho_min_local = cell_double_operator.min(rho_id);
            
            double m_max_local = cell_double_operator.max(m_id);
            double m_min_local = cell_double_operator.min(m_id);
            
            double E_max_local = cell_double_operator.max(E_id);
            double E_min_local = cell_double_operator.min(E_id);
            
            double rho_max_global = 0.0;
            double rho_min_global = 0.0;
            double m_max_global = 0.0;
            double m_min_global = 0.0;
            double E_max_global = 0.0;
            double E_min_global = 0.0;
            
            mpi.Allreduce(
                &rho_max_local,
                &rho_max_global,
                1,
                MPI_DOUBLE,
                MPI_MAX);
            
            mpi.Allreduce(
                &rho_min_local,
                &rho_min_global,
                1,
                MPI_DOUBLE,
                MPI_MAX);
            
            mpi.Allreduce(
                &m_max_local,
                &m_max_global,
                1,
                MPI_DOUBLE,
                MPI_MAX);
            
            mpi.Allreduce(
                &m_min_local,
                &m_min_global,
                1,
                MPI_DOUBLE,
                MPI_MAX);
            
            mpi.Allreduce(
                &E_max_local,
                &E_max_global,
                1,
                MPI_DOUBLE,
                MPI_MAX);
            
            mpi.Allreduce(
                &E_min_local,
                &E_min_global,
                1,
                MPI_DOUBLE,
                MPI_MAX);
            
            os << "Max/min density: " << rho_max_global << "/" << rho_min_global << std::endl;
            os << "Max/min momentum component: " << m_max_global << "/" << m_min_global << std::endl;
            os << "Max/min total energy: " << E_max_global << "/" << E_min_global << std::endl;
            
            break;
        }
        case FOUR_EQN_SHYUE:
        {
            const int rho_id = variable_db->mapVariableAndContextToIndex(
                d_density,
                d_plot_context);
            
            const int m_id = variable_db->mapVariableAndContextToIndex(
                d_momentum,
                d_plot_context);
            
            const int E_id = variable_db->mapVariableAndContextToIndex(
                d_total_energy,
                d_plot_context);
            
            const int Y_id = variable_db->mapVariableAndContextToIndex(
                d_mass_fraction,
                d_plot_context);
            
            double rho_max_local = cell_double_operator.max(rho_id);
            double rho_min_local = cell_double_operator.min(rho_id);
            
            double m_max_local = cell_double_operator.max(m_id);
            double m_min_local = cell_double_operator.min(m_id);
            
            double E_max_local = cell_double_operator.max(E_id);
            double E_min_local = cell_double_operator.min(E_id);
            
            double Y_max_local = cell_double_operator.max(Y_id);
            double Y_min_local = cell_double_operator.min(Y_id);
            
            double rho_max_global = 0.0;
            double rho_min_global = 0.0;
            double m_max_global = 0.0;
            double m_min_global = 0.0;
            double E_max_global = 0.0;
            double E_min_global = 0.0;
            double Y_max_global = 0.0;
            double Y_min_global = 0.0;
            
            mpi.Allreduce(
                &rho_max_local,
                &rho_max_global,
                1,
                MPI_DOUBLE,
                MPI_MAX);
            
            mpi.Allreduce(
                &rho_min_local,
                &rho_min_global,
                1,
                MPI_DOUBLE,
                MPI_MAX);
            
            mpi.Allreduce(
                &m_max_local,
                &m_max_global,
                1,
                MPI_DOUBLE,
                MPI_MAX);
            
            mpi.Allreduce(
                &m_min_local,
                &m_min_global,
                1,
                MPI_DOUBLE,
                MPI_MAX);
            
            mpi.Allreduce(
                &E_max_local,
                &E_max_global,
                1,
                MPI_DOUBLE,
                MPI_MAX);
            
            mpi.Allreduce(
                &E_min_local,
                &E_min_global,
                1,
                MPI_DOUBLE,
                MPI_MAX);
            
            mpi.Allreduce(
                &Y_max_local,
                &Y_max_global,
                1,
                MPI_DOUBLE,
                MPI_MAX);
            
            mpi.Allreduce(
                &Y_min_local,
                &Y_min_global,
                1,
                MPI_DOUBLE,
                MPI_MAX);
            
            os << "Max/min density: " << rho_max_global << "/" << rho_min_global << std::endl;
            os << "Max/min momentum component: " << m_max_global << "/" << m_min_global << std::endl;
            os << "Max/min total energy: " << E_max_global << "/" << E_min_global << std::endl;
            os << "Max/min mass fraction component: " << Y_max_global << "/" << Y_min_global << std::endl;
            
            break;
        }
        case FIVE_EQN_ALLAIRE:
        {
            const int Z_rho_id = variable_db->mapVariableAndContextToIndex(
                d_partial_density,
                d_plot_context);
            
            const int m_id = variable_db->mapVariableAndContextToIndex(
                d_momentum,
                d_plot_context);
            
            const int E_id = variable_db->mapVariableAndContextToIndex(
                d_total_energy,
                d_plot_context);
            
            const int Z_id = variable_db->mapVariableAndContextToIndex(
                d_volume_fraction,
                d_plot_context);
            
            double Z_rho_max_local = cell_double_operator.max(Z_rho_id);
            double Z_rho_min_local = cell_double_operator.min(Z_rho_id);
            
            double m_max_local = cell_double_operator.max(m_id);
            double m_min_local = cell_double_operator.min(m_id);
            
            double E_max_local = cell_double_operator.max(E_id);
            double E_min_local = cell_double_operator.min(E_id);
            
            double Z_max_local = cell_double_operator.max(Z_id);
            double Z_min_local = cell_double_operator.min(Z_id);
            
            double Z_rho_max_global = 0.0;
            double Z_rho_min_global = 0.0;
            double m_max_global = 0.0;
            double m_min_global = 0.0;
            double E_max_global = 0.0;
            double E_min_global = 0.0;
            double Z_max_global = 0.0;
            double Z_min_global = 0.0;
            
            mpi.Allreduce(
                &Z_rho_max_local,
                &Z_rho_max_global,
                1,
                MPI_DOUBLE,
                MPI_MAX);
            
            mpi.Allreduce(
                &Z_rho_min_local,
                &Z_rho_min_global,
                1,
                MPI_DOUBLE,
                MPI_MAX);
            
            mpi.Allreduce(
                &m_max_local,
                &m_max_global,
                1,
                MPI_DOUBLE,
                MPI_MAX);
            
            mpi.Allreduce(
                &m_min_local,
                &m_min_global,
                1,
                MPI_DOUBLE,
                MPI_MAX);
            
            mpi.Allreduce(
                &E_max_local,
                &E_max_global,
                1,
                MPI_DOUBLE,
                MPI_MAX);
            
            mpi.Allreduce(
                &E_min_local,
                &E_min_global,
                1,
                MPI_DOUBLE,
                MPI_MAX);
            
            mpi.Allreduce(
                &Z_max_local,
                &Z_max_global,
                1,
                MPI_DOUBLE,
                MPI_MAX);
            
            mpi.Allreduce(
                &Z_min_local,
                &Z_min_global,
                1,
                MPI_DOUBLE,
                MPI_MAX);
            
            os << "Max/min partial density component: " << Z_rho_max_global << "/" << Z_rho_min_global << std::endl;
            os << "Max/min momentum component: " << m_max_global << "/" << m_min_global << std::endl;
            os << "Max/min total energy: " << E_max_global << "/" << E_min_global << std::endl;
            os << "Max/min volume fraction component: " << Z_max_global << "/" << Z_min_global << std::endl;
            
            break;
        }
        default:
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "Unknown d_flow_model."
                       << std::endl);
        }
    }
}


void
Euler::getFromInput(
    boost::shared_ptr<tbox::Database> input_db,
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
                       << "Key data 'num_species' not found in input."
                       << " Number of species is unknown."
                       << std::endl);
        }
        
        /*
         * Initialize the flow model.
         */
        if (input_db->keyExists("flow_model"))
        {
            std::string flow_model_str = input_db->getString("flow_model");
            
            if (flow_model_str == "SINGLE_SPECIES")
            {
                d_flow_model = SINGLE_SPECIES;
                d_num_eqn = 2 + d_dim.getValue();
            }
            else if (flow_model_str == "FOUR_EQN_SHYUE")
            {
                d_flow_model = FOUR_EQN_SHYUE;
                d_num_eqn = 1 + d_dim.getValue() + d_num_species;
            }
            else if (flow_model_str == "FIVE_EQN_ALLAIRE")
            {
                d_flow_model = FIVE_EQN_ALLAIRE;
                d_num_eqn = d_dim.getValue() + 2*d_num_species;
            }
            else
            {
                TBOX_ERROR(d_object_name
                           << ": "
                           << "Unknown flow_model string = "
                           << flow_model_str
                           << " found in input."
                           << std::endl);        
            }
            
            if (d_num_species > 1 && d_flow_model == SINGLE_SPECIES)
            {
                TBOX_ERROR(d_object_name
                           << ": "
                           << "Number of species = "
                           << d_num_species
                           << " shouldn't use single-species model."
                           << std::endl); 
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "Key data 'flow model' not found in input."
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
                       << "Key data 'Equation_of_state' not found in input."
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
                       << "Key data 'Shock_capturing_scheme' not found in input."
                       << std::endl);
        }
    }
    
    /*
     * Defaults for boundary conditions. Set to bogus values
     * for error checking.
     */
    setDefaultBoundaryConditions();
    
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
            boost::shared_ptr<tbox::Database> bdry_db = input_db->getDatabase(
                "Boundary_data");
            
            if (d_dim == tbox::Dimension(1))
            {
                // NOT YET IMPLEMENTED
            }
            if (d_dim == tbox::Dimension(2))
            {
                appu::CartesianBoundaryUtilities2::getFromInput(
                    this,
                    bdry_db,
                    d_master_bdry_edge_conds,
                    d_master_bdry_node_conds,
                    periodic);
            }
            else if (d_dim == tbox::Dimension(3))
            {
                appu::CartesianBoundaryUtilities3::getFromInput(
                    this,
                    bdry_db,
                    d_master_bdry_face_conds,
                    d_master_bdry_edge_conds,
                    d_master_bdry_node_conds,
                    periodic);
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "Key data 'Boundary_data' not found in input. "
                       << std::endl);
        }
    }
    
    if (input_db->keyExists("Refinement_data"))
    {
        boost::shared_ptr<tbox::Database> refine_db(
            input_db->getDatabase("Refinement_data"));
        std::vector<std::string> refinement_keys = refine_db->getAllKeys();
        int num_keys = static_cast<int>(refinement_keys.size());
        
        if (refine_db->keyExists("refine_criteria"))
        {
            d_refinement_criteria = refine_db->getStringVector("refine_criteria");
        }
        else
        {
            TBOX_WARNING(
                d_object_name << ": "
                              << "No key 'refine_criteria' found in data for"
                              << " RefinementData. No refinement will occur."
                              << std::endl);
        }
        
        std::vector<std::string> ref_keys_defined(num_keys);
        int def_key_cnt = 0;
        boost::shared_ptr<tbox::Database> error_db;
        for (int i = 0; i < num_keys; i++)
        {
            std::string error_key = refinement_keys[i];
            error_db.reset();
            
            if (!(error_key == "refine_criteria"))
            {
                if (!(error_key == "DENSITY_SHOCK" ||
                      error_key == "PRESSURE_SHOCK"))
                {
                    TBOX_ERROR(d_object_name
                               << ": "
                               << "Unknown refinement criteria: "
                               << error_key
                               << "\nin input."
                               << std::endl);
                }
                else
                {
                    error_db = refine_db->getDatabase(error_key);
                    ref_keys_defined[def_key_cnt] = error_key;
                    ++def_key_cnt;
                }
                
                if (error_db && error_key == "DENSITY_SHOCK")
                {                    
                    if (error_db->keyExists("shock_tol"))
                    {
                        d_density_shock_tol = error_db->getDoubleVector("shock_tol");
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name
                                   << ": "
                                   << "No key 'shock_tol' found in data for "
                                   << error_key
                                   << "."
                                   << std::endl);
                    }
                }
                
                if (error_db && error_key == "PRESSURE_SHOCK")
                {        
                    if (error_db->keyExists("shock_tol"))
                    {
                        d_pressure_shock_tol = error_db->getDoubleVector("shock_tol");
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name
                                   << ": "
                                   << "No key 'shock_tol' found in data for "
                                   << error_key
                                   << "."
                                   << std::endl);
                    }
                }
            }
        } // loop over refine criteria
        
        /*
         * Check that input is found for each string identifier in key list.
         */
        for (int k0 = 0;
             k0 < static_cast<int>(d_refinement_criteria.size());
             k0++)
        {
            std::string use_key = d_refinement_criteria[k0];
            bool key_found = false;
            for (int k1 = 0; k1 < def_key_cnt; k1++)
            {
                std::string def_key = ref_keys_defined[k1];
                if (def_key == use_key)
                {
                    key_found = true;
                }
            }
            
            if (!key_found)
            {
                TBOX_ERROR(d_object_name
                           << ": "
                           << "No input found for specified refine criteria: "
                           << d_refinement_criteria[k0]
                           << "."
                           << std::endl);
            }
        }
    } // if "Refinement_data" db entry exists
    
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
    
    std::string flow_model_str = db->getString("d_flow_model");
    if (flow_model_str == "SINGLE_SPECIES")
    {
        d_flow_model = SINGLE_SPECIES;
        d_num_eqn = 2 + d_dim.getValue();
    }
    else if (flow_model_str == "FOUR_EQN_SHYUE")
    {
        d_flow_model = FOUR_EQN_SHYUE;
        d_num_eqn = 1 + d_dim.getValue() + d_num_species;
    }
    else if (flow_model_str == "FIVE_EQN_ALLAIRE")
    {
        d_flow_model = FIVE_EQN_ALLAIRE;
        d_num_eqn = d_dim.getValue() + 2*d_num_species;
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "Unknown d_flow_model string = "
                   << flow_model_str
                   << " found in restart file."
                   << std::endl);        
    }
    
    d_is_preserving_positivity = db->getBool("d_is_preserving_positivity");
    
    d_equation_of_state_db = db->getDatabase("Equation_of_state");
    
    d_shock_capturing_scheme_db = db->getDatabase("Shock_capturing_scheme");
    
    int* tmp_num_ghosts = &d_num_ghosts[0];
    db->getIntegerArray("d_num_ghosts", tmp_num_ghosts, d_dim.getValue());
    
    /*
     * Defaults for boundary conditions. Set to bogus values
     * for error checking.
     */
    setDefaultBoundaryConditions();
    
    d_master_bdry_node_conds = db->getIntegerVector("d_master_bdry_node_conds");
    
    if (d_dim == tbox::Dimension(1))
    {
        switch (d_flow_model)
        {
            case SINGLE_SPECIES:
            {
                d_bdry_node_density = db->getDoubleVector("d_bdry_node_density");
                d_bdry_node_momentum = db->getDoubleVector("d_bdry_node_momentum");
                d_bdry_node_total_energy = db->getDoubleVector("d_bdry_node_total_energy");
                
                break;
            }
            case FOUR_EQN_SHYUE:
            {
                d_bdry_node_density = db->getDoubleVector("d_bdry_node_density");
                d_bdry_node_momentum = db->getDoubleVector("d_bdry_node_momentum");
                d_bdry_node_total_energy = db->getDoubleVector("d_bdry_node_total_energy");
                d_bdry_node_mass_fraction = db->getDoubleVector("d_bdry_node_mass_fraction");
                
                break;
            }
            case FIVE_EQN_ALLAIRE:
            {
                d_bdry_node_partial_density = db->getDoubleVector("d_bdry_node_partial_density");
                d_bdry_node_momentum = db->getDoubleVector("d_bdry_node_momentum");
                d_bdry_node_total_energy = db->getDoubleVector("d_bdry_node_total_energy");
                d_bdry_node_volume_fraction = db->getDoubleVector("d_bdry_node_volume_fraction");
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                           << ": "
                           << "Unknown d_flow_model."
                           << std::endl);
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        d_master_bdry_edge_conds = db->getIntegerVector("d_master_bdry_edge_conds");
        
        switch (d_flow_model)
        {
            case SINGLE_SPECIES:
            {
                d_bdry_edge_density = db->getDoubleVector("d_bdry_edge_density");
                d_bdry_edge_momentum = db->getDoubleVector("d_bdry_edge_momentum");
                d_bdry_edge_total_energy = db->getDoubleVector("d_bdry_edge_total_energy");
                
                break;
            }
            case FOUR_EQN_SHYUE:
            {
                d_bdry_edge_density = db->getDoubleVector("d_bdry_edge_density");
                d_bdry_edge_momentum = db->getDoubleVector("d_bdry_edge_momentum");
                d_bdry_edge_total_energy = db->getDoubleVector("d_bdry_edge_total_energy");
                d_bdry_edge_mass_fraction = db->getDoubleVector("d_bdry_edge_mass_fraction");
                
                break;
            }
            case FIVE_EQN_ALLAIRE:
            {
                d_bdry_edge_partial_density = db->getDoubleVector("d_bdry_edge_partial_density");
                d_bdry_edge_momentum = db->getDoubleVector("d_bdry_edge_momentum");
                d_bdry_edge_total_energy = db->getDoubleVector("d_bdry_edge_total_energy");
                d_bdry_edge_volume_fraction = db->getDoubleVector("d_bdry_edge_volume_fraction");
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                           << ": "
                           << "Unknown d_flow_model."
                           << std::endl);
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        d_master_bdry_edge_conds = db->getIntegerVector("d_master_bdry_edge_conds");
        d_master_bdry_face_conds = db->getIntegerVector("d_master_bdry_face_conds");
        
        switch (d_flow_model)
        {
            case SINGLE_SPECIES:
            {
                d_bdry_face_density = db->getDoubleVector("d_bdry_face_density");
                d_bdry_face_momentum = db->getDoubleVector("d_bdry_face_momentum");
                d_bdry_face_total_energy = db->getDoubleVector("d_bdry_face_total_energy");
                
                break;
            }
            case FOUR_EQN_SHYUE:
            {
                d_bdry_face_density = db->getDoubleVector("d_bdry_face_density");
                d_bdry_face_momentum = db->getDoubleVector("d_bdry_face_momentum");
                d_bdry_face_total_energy = db->getDoubleVector("d_bdry_face_total_energy");
                d_bdry_face_mass_fraction = db->getDoubleVector("d_bdry_face_mass_fraction");
                
                break;
            }
            case FIVE_EQN_ALLAIRE:
            {
                d_bdry_face_partial_density = db->getDoubleVector("d_bdry_face_partial_density");
                d_bdry_face_momentum = db->getDoubleVector("d_bdry_face_momentum");
                d_bdry_face_total_energy = db->getDoubleVector("d_bdry_face_total_energy");
                d_bdry_face_volume_fraction = db->getDoubleVector("d_bdry_face_volume_fraction");
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                           << ": "
                           << "Unknown d_flow_model."
                           << std::endl);
            }
        }
    }
    
    for (int i = 0; i < static_cast<int>(d_refinement_criteria.size()); i++)
    {
        if (d_refinement_criteria[i] == "DENSITY_SHOCK")
        {
            d_density_shock_tol = db->getDoubleVector("d_density_shock_tol");
    
        }
        else if (d_refinement_criteria[i] == "PRESSURE_SHOCK")
        {
            d_pressure_shock_tol = db->getDoubleVector("d_pressure_shock_tol");
        }    
    }
}


void
Euler::readStateDataEntryForSingleSpecies(
    boost::shared_ptr<tbox::Database> db,
    const std::string& db_name,
    int array_indx,
    std::vector<double>& density,
    std::vector<double>& momentum,
    std::vector<double>& total_energy)
{
    TBOX_ASSERT(db);
    TBOX_ASSERT(!db_name.empty());
    TBOX_ASSERT(array_indx >= 0);
    TBOX_ASSERT(static_cast<int>(density.size()) > array_indx);
    TBOX_ASSERT(static_cast<int>(momentum.size()) > array_indx*d_dim.getValue());
    TBOX_ASSERT(static_cast<int>(total_energy.size()) > array_indx);
    
    if (db->keyExists("density"))
    {
        density[array_indx] = db->getDouble("density");
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "'density' entry missing from "
                   << db_name
                   << " input database."
                   << std::endl);
    }
    
    if (db->keyExists("momentum"))
    {
        std::vector<double> tmp_m = db->getDoubleVector("momentum");
        if (static_cast<int>(tmp_m.size()) < d_dim.getValue())
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "Insufficient number of 'momentum' values"
                       << " given in "
                       << db_name
                       << " input database."
                       << std::endl);
        }
        for (int di = 0; di < d_dim.getValue(); di++)
        {
            momentum[array_indx*d_dim.getValue() + di] = tmp_m[di];
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "'momentum' entry missing from "
                   << db_name
                   << " input database."
                   << std::endl);
    }
   
    if (db->keyExists("total_energy"))
    {
        total_energy[array_indx] = db->getDouble("total_energy");
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "'total_energy' entry missing from "
                   << db_name
                   << " input database."
                   << std::endl);
    }
}


void
Euler::readStateDataEntryForFourEqnShyue(
    boost::shared_ptr<tbox::Database> db,
    const std::string& db_name,
    int array_indx,
    std::vector<double>& density,
    std::vector<double>& momentum,
    std::vector<double>& total_energy,
    std::vector<double>& mass_fraction)
{
    TBOX_ASSERT(db);
    TBOX_ASSERT(!db_name.empty());
    TBOX_ASSERT(array_indx >= 0);
    TBOX_ASSERT(static_cast<int>(density.size()) > array_indx);
    TBOX_ASSERT(static_cast<int>(momentum.size()) > array_indx*d_dim.getValue());
    TBOX_ASSERT(static_cast<int>(total_energy.size()) > array_indx);
    TBOX_ASSERT(static_cast<int>(mass_fraction.size()) > array_indx*d_num_species);
    
    if (db->keyExists("density"))
    {
        density[array_indx] = db->getDouble("density");
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "'density' entry missing from "
                   << db_name
                   << " input database."
                   << std::endl);
    }
    
    if (db->keyExists("momentum"))
    {
        std::vector<double> tmp_m = db->getDoubleVector("momentum");
        if (static_cast<int>(tmp_m.size()) < d_dim.getValue())
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "Insufficient number of 'momentum' values"
                       << " given in "
                       << db_name
                       << " input database."
                       << std::endl);
        }
        for (int di = 0; di < d_dim.getValue(); di++)
        {
            momentum[array_indx*d_dim.getValue() + di] = tmp_m[di];
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "'momentum' entry missing from "
                   << db_name
                   << " input database."
                   << std::endl);
    }
   
    if (db->keyExists("total_energy"))
    {
        total_energy[array_indx] = db->getDouble("total_energy");
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "'total_energy' entry missing from "
                   << db_name
                   << " input database."
                   << std::endl);
    }
    
    if (db->keyExists("mass_fraction"))
    {
        std::vector<double> tmp_Y = db->getDoubleVector("mass_fraction");
        if (static_cast<int>(tmp_Y.size()) < d_num_species - 1)
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "Insufficient number of 'mass_fraction' values"
                       << " given in "
                       << db_name
                       << " input database."
                       << std::endl);
        }
        double Y_last = 1.0;
        for (int si = 0; si < d_num_species - 1; si++)
        {
            mass_fraction[array_indx*d_num_species + si] = tmp_Y[si];
            Y_last -= tmp_Y[si];
        }
        mass_fraction[(array_indx + 1)*d_num_species - 1] = Y_last;
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "'mass_fraction' entry missing from "
                   << db_name
                   << " input database."
                   << std::endl);
    }
}


void
Euler::readStateDataEntryForFiveEqnAllaire(
    boost::shared_ptr<tbox::Database> db,
    const std::string& db_name,
    int array_indx,
    std::vector<double>& partial_density,
    std::vector<double>& momentum,
    std::vector<double>& total_energy,
    std::vector<double>& volume_fraction)
{
    TBOX_ASSERT(db);
    TBOX_ASSERT(!db_name.empty());
    TBOX_ASSERT(array_indx >= 0);
    TBOX_ASSERT(static_cast<int>(partial_density.size()) > array_indx*d_num_species);
    TBOX_ASSERT(static_cast<int>(momentum.size()) > array_indx*d_dim.getValue());
    TBOX_ASSERT(static_cast<int>(total_energy.size()) > array_indx);
    TBOX_ASSERT(static_cast<int>(volume_fraction.size()) > array_indx*d_num_species);
    
    if (db->keyExists("partial_density"))
    {
        std::vector<double> tmp_Z_rho = db->getDoubleVector("partial_density");
        if (static_cast<int>(tmp_Z_rho.size()) < d_num_species)
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "Insufficient number of 'partial_density' values"
                       << " given in "
                       << db_name
                       << " input database."
                       << std::endl);
        }
        for (int si = 0; si < d_num_species; si++)
        {
            partial_density[array_indx*d_num_species + si] = tmp_Z_rho[si];
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "'partial_density' entry missing from "
                   << db_name
                   << " input database."
                   << std::endl);
    }
    
    if (db->keyExists("momentum"))
    {
        std::vector<double> tmp_m = db->getDoubleVector("momentum");
        if (static_cast<int>(tmp_m.size()) < d_dim.getValue())
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "Insufficient number of 'momentum' values"
                       << " given in "
                       << db_name
                       << " input database."
                       << std::endl);
        }
        for (int di = 0; di < d_dim.getValue(); di++)
        {
            momentum[array_indx*d_dim.getValue() + di] = tmp_m[di];
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "'momentum' entry missing from "
                   << db_name
                   << " input database."
                   << std::endl);
    }
   
    if (db->keyExists("total_energy"))
    {
        total_energy[array_indx] = db->getDouble("total_energy");
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "'total_energy' entry missing from "
                   << db_name
                   << " input database."
                   << std::endl);
    }
    
    if (db->keyExists("volume_fraction"))
    {
        std::vector<double> tmp_Z = db->getDoubleVector("volume_fraction");
        if (static_cast<int>(tmp_Z.size()) < d_num_species - 1)
        {
            TBOX_ERROR(d_object_name
                       << ": "
                       << "Insufficient number of 'volume_fraction' values"
                       << " given in "
                       << db_name
                       << " input database."
                       << std::endl);
        }
        double Z_last = 1.0;
        for (int si = 0; si < d_num_species - 1; si++)
        {
            volume_fraction[array_indx*d_num_species + si] = tmp_Z[si];
            Z_last -= tmp_Z[si];
        }
        volume_fraction[(array_indx + 1)*d_num_species - 1] = Z_last;
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "'volume_fraction' entry missing from "
                   << db_name
                   << " input database."
                   << std::endl);
    }
}

void
Euler::setDefaultBoundaryConditions()
{
    if (d_dim == tbox::Dimension(1))
    {
        d_master_bdry_node_conds.resize(NUM_1D_NODES);
        d_scalar_bdry_node_conds.resize(NUM_1D_NODES);
        d_vector_bdry_node_conds.resize(NUM_1D_NODES);
        for (int ni = 0; ni < NUM_1D_NODES; ni++)
        {
            d_master_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
            d_scalar_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
            d_vector_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
        }
        
        switch (d_flow_model)
        {
            case SINGLE_SPECIES:
            {
                d_bdry_node_density.resize(NUM_1D_NODES);
                d_bdry_node_momentum.resize(NUM_1D_NODES*d_dim.getValue());
                d_bdry_node_total_energy.resize(NUM_1D_NODES);
                
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_node_density);
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_node_momentum);
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_node_total_energy);
                
                break;
            }
            case FOUR_EQN_SHYUE:
            {
                d_bdry_node_density.resize(NUM_1D_NODES);
                d_bdry_node_momentum.resize(NUM_1D_NODES*d_dim.getValue());
                d_bdry_node_total_energy.resize(NUM_1D_NODES);
                d_bdry_node_mass_fraction.resize(NUM_1D_NODES*d_num_species);
                
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_node_density);
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_node_momentum);
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_node_total_energy);
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_node_mass_fraction);
                
                break;
            }
            case FIVE_EQN_ALLAIRE:
            {
                d_bdry_node_partial_density.resize(NUM_1D_NODES*d_num_species);
                d_bdry_node_momentum.resize(NUM_1D_NODES*d_dim.getValue());
                d_bdry_node_total_energy.resize(NUM_1D_NODES);
                d_bdry_node_volume_fraction.resize(NUM_1D_NODES*d_num_species);
                
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_node_partial_density);
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_node_momentum);
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_node_total_energy);
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_node_volume_fraction);
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                           << ": "
                           << "Unknown d_flow_model."
                           << std::endl);
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        d_master_bdry_edge_conds.resize(NUM_2D_EDGES);
        d_scalar_bdry_edge_conds.resize(NUM_2D_EDGES);
        d_vector_bdry_edge_conds.resize(NUM_2D_EDGES);
        for (int ei = 0; ei < NUM_2D_EDGES; ei++)
        {
            d_master_bdry_edge_conds[ei] = BOGUS_BDRY_DATA;
            d_scalar_bdry_edge_conds[ei] = BOGUS_BDRY_DATA;
            d_vector_bdry_edge_conds[ei] = BOGUS_BDRY_DATA;
        }
        
        d_master_bdry_node_conds.resize(NUM_2D_NODES);
        d_scalar_bdry_node_conds.resize(NUM_2D_NODES);
        d_vector_bdry_node_conds.resize(NUM_2D_NODES);
        d_node_bdry_edge.resize(NUM_2D_NODES);
        for (int ni = 0; ni < NUM_2D_NODES; ni++)
        {
            d_master_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
            d_scalar_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
            d_vector_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
            d_node_bdry_edge[ni] = BOGUS_BDRY_DATA;
        }
        
        switch (d_flow_model)
        {
            case SINGLE_SPECIES:
            {
                d_bdry_edge_density.resize(NUM_2D_EDGES);
                d_bdry_edge_momentum.resize(NUM_2D_EDGES*d_dim.getValue());
                d_bdry_edge_total_energy.resize(NUM_2D_EDGES);
                
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_edge_density);
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_edge_momentum);
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_edge_total_energy);
                
                break;
            }
            case FOUR_EQN_SHYUE:
            {
                d_bdry_edge_density.resize(NUM_2D_EDGES);
                d_bdry_edge_momentum.resize(NUM_2D_EDGES*d_dim.getValue());
                d_bdry_edge_total_energy.resize(NUM_2D_EDGES);
                d_bdry_edge_mass_fraction.resize(NUM_2D_NODES*d_num_species);
                
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_edge_density);
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_edge_momentum);
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_edge_total_energy);
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_edge_mass_fraction);
                
                break;
            }
            case FIVE_EQN_ALLAIRE:
            {
                d_bdry_edge_partial_density.resize(NUM_2D_EDGES*d_num_species);
                d_bdry_edge_momentum.resize(NUM_2D_EDGES*d_dim.getValue());
                d_bdry_edge_total_energy.resize(NUM_2D_EDGES);
                d_bdry_edge_volume_fraction.resize(NUM_2D_EDGES*d_num_species);
                
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_edge_density);
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_edge_momentum);
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_edge_total_energy);
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_edge_volume_fraction);
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                           << ": "
                           << "Unknown d_flow_model."
                           << std::endl);
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        d_master_bdry_face_conds.resize(NUM_3D_FACES);
        d_scalar_bdry_face_conds.resize(NUM_3D_FACES);
        d_vector_bdry_face_conds.resize(NUM_3D_FACES);
        for (int fi = 0; fi < NUM_3D_FACES; fi++) {
           d_master_bdry_face_conds[fi] = BOGUS_BDRY_DATA;
           d_scalar_bdry_face_conds[fi] = BOGUS_BDRY_DATA;
           d_vector_bdry_face_conds[fi] = BOGUS_BDRY_DATA;
        }
     
        d_master_bdry_edge_conds.resize(NUM_3D_EDGES);
        d_scalar_bdry_edge_conds.resize(NUM_3D_EDGES);
        d_vector_bdry_edge_conds.resize(NUM_3D_EDGES);
        d_edge_bdry_face.resize(NUM_3D_EDGES);
        for (int ei = 0; ei < NUM_3D_EDGES; ei++)
        {
            d_master_bdry_edge_conds[ei] = BOGUS_BDRY_DATA;
            d_scalar_bdry_edge_conds[ei] = BOGUS_BDRY_DATA;
            d_vector_bdry_edge_conds[ei] = BOGUS_BDRY_DATA;
            d_edge_bdry_face[ei] = BOGUS_BDRY_DATA;
        }
     
        d_master_bdry_node_conds.resize(NUM_3D_NODES);
        d_scalar_bdry_node_conds.resize(NUM_3D_NODES);
        d_vector_bdry_node_conds.resize(NUM_3D_NODES);
        d_node_bdry_face.resize(NUM_3D_NODES);
        for (int ni = 0; ni < NUM_3D_NODES; ni++)
        {
            d_master_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
            d_scalar_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
            d_vector_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
            d_node_bdry_face[ni] = BOGUS_BDRY_DATA;
        }
        
        switch (d_flow_model)
        {
            case SINGLE_SPECIES:
            {
                d_bdry_face_density.resize(NUM_3D_FACES);
                d_bdry_face_momentum.resize(NUM_3D_FACES*d_dim.getValue());
                d_bdry_face_total_energy.resize(NUM_3D_FACES);
                
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_face_density);
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_face_momentum);
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_face_total_energy);
                
                break;
            }
            case FOUR_EQN_SHYUE:
            {
                d_bdry_face_density.resize(NUM_3D_FACES);
                d_bdry_face_momentum.resize(NUM_3D_FACES*d_dim.getValue());
                d_bdry_face_total_energy.resize(NUM_3D_FACES);
                d_bdry_face_mass_fraction.resize(NUM_3D_FACES*d_num_species);
                
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_face_density);
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_face_momentum);
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_face_total_energy);
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_face_mass_fraction);
                
                break;
            }
            case FIVE_EQN_ALLAIRE:
            {
                d_bdry_face_partial_density.resize(NUM_3D_FACES*d_num_species);
                d_bdry_face_momentum.resize(NUM_3D_FACES*d_dim.getValue());
                d_bdry_face_total_energy.resize(NUM_3D_FACES);
                d_bdry_face_volume_fraction.resize(NUM_3D_FACES*d_num_species);
                
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_edge_partial_density);
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_edge_momentum);
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_edge_total_energy);
                tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_edge_volume_fraction);
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                           << ": "
                           << "Unknown d_flow_model."
                           << std::endl);
            }
        }
    }
}


void
Euler::preservePositivity(
    std::vector<double*>& Q,
    boost::shared_ptr<pdat::FaceData<double> >& convective_flux_intermediate,
    boost::shared_ptr<pdat::CellData<double> >& source_intermediate,
    const hier::IntVector interior_dims,
    const hier::IntVector ghostcell_dims,
    const double* const dx,
    const double& dt,
    const double& beta)
{
    /*
     * Set the lower and upper bounds for mass and volume fractions.
     */
    const double Y_bnd_lo = -0.1;
    const double Y_bnd_up = 1.1;
    
    const double Z_bnd_lo = -0.1;
    const double Z_bnd_up = 1.1;
    
    NULL_USE(Z_bnd_lo);
    NULL_USE(Z_bnd_up);
    
    TBOX_ASSERT(beta != 0.0);
    
    std::vector<double> Q_tmp;
    Q_tmp.resize(d_num_eqn);
    
    if (d_dim == tbox::Dimension(1))
    {
        for (int i = 0; i < interior_dims[0]; i++)
        {
            // Compute indices of time-dependent data, flux and source.
            const int idx_cell   = i + d_num_ghosts[0];
            const int idx_source = i;
            const int idx_face_x = i + 1;
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                double* F_x_intermediate = convective_flux_intermediate->getPointer(0, ei);
                double* S_intermediate = source_intermediate->getPointer(ei);
                
                Q_tmp[ei] = Q[ei][idx_cell] + beta*
                    (-(F_x_intermediate[idx_face_x] - F_x_intermediate[idx_face_x - 1])/dx[0] +
                    S_intermediate[idx_source]);
            }
            
            double rho = 0;
            double p   = 0;
            
            /*
             * Check if the density or pressure are negative and if the mass or volume fractions
             * are within the bound. If they are, switch to first order hyperbolic flux.
             */
            
            switch (d_flow_model)
            {
                case SINGLE_SPECIES:
                {
                    std::vector<const double*> m_ptr;
                    for (int di = 0; di < d_dim.getValue(); di++)
                    {
                        m_ptr.push_back(&Q_tmp[1 + di]);
                    }
                    
                    rho = Q_tmp[0];
                    const double& E = Q_tmp[1 + d_dim.getValue()];
                    
                    p = d_equation_of_state->getPressure(
                        &rho,
                        m_ptr,
                        &E);
                    
                    break;
                }
                case FOUR_EQN_SHYUE:
                {
                    std::vector<const double*> m_ptr;
                    for (int di = 0; di < d_dim.getValue(); di++)
                    {
                        m_ptr.push_back(&Q_tmp[1 + di]);
                    }
                    
                    std::vector<const double*> Y_ptr;
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Y_ptr.push_back(&Q_tmp[2 + d_dim.getValue() + si]);
                    }
                    
                    rho = Q_tmp[0];
                    const double& E = Q_tmp[1 + d_dim.getValue()];
                    
                    p = d_equation_of_state->getPressureWithMassFraction(
                        &rho,
                        m_ptr,
                        &E,
                        Y_ptr);
                    
                    break;
                }
                case FIVE_EQN_ALLAIRE:
                {
                    std::vector<const double*> Z_rho_ptr;
                    for (int si = 0; si < d_num_species; si++)
                    {
                        Z_rho_ptr.push_back(&Q_tmp[si]);
                    }
                    
                    std::vector<const double*> m_ptr;
                    for (int di = 0; di < d_dim.getValue(); di++)
                    {
                        m_ptr.push_back(&Q_tmp[d_num_species + di]);
                    }
                    
                    std::vector<const double*> Z_ptr;
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Z_ptr.push_back(&Q_tmp[1 + d_num_species + d_dim.getValue() + si]);
                    }
                    
                    rho = d_equation_of_state->getTotalDensity(Z_rho_ptr);
                    const double& E = Q_tmp[d_num_species + d_dim.getValue()];
                    
                    p = d_equation_of_state->getPressureWithVolumeFraction(
                        &rho,
                        m_ptr,
                        &E,
                        Z_ptr);
                    
                    break;
                }
                default:
                    TBOX_ERROR(d_object_name
                               << ": "
                               << "Unknown d_flow_model."
                               << std::endl);
            }
            
            bool is_first_order = false;
            
            // Check if the density or pressure are negative.
            
            if (rho < 0 || p < 0)
            {
                is_first_order = true;
            }
            
            // Check if the mass or volume fractions are within the bound.
            
            switch (d_flow_model)
            {
                case SINGLE_SPECIES:
                {
                    break;
                }
                case FOUR_EQN_SHYUE:
                {
                    double Y_last = 1.0;
                    
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        double Y = Q_tmp[2 + d_dim.getValue() + si];
                        
                        Y_last -= Y;
                        
                        if (Y < Y_bnd_lo || Y > Y_bnd_up)
                        {
                            is_first_order = true;
                        }
                    }
                    
                    if (Y_last < Y_bnd_lo || Y_last > Y_bnd_up)
                    {
                        is_first_order = true;
                    }
                    
                    break;
                }
                case FIVE_EQN_ALLAIRE:
                {
                    for (int si = 0; si < d_num_species; si++)
                    {
                        double Y = Q_tmp[si]/rho;
                        
                        if (Y < Y_bnd_lo || Y > Y_bnd_up)
                        {
                            is_first_order = true;
                            break;
                        }
                    }
                    
                    /*
                    double Z_last = 1.0;
                    
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        double Z = Q_tmp[1 + d_num_species + d_dim.getValue() + si];
                        
                        Z_last -= Z;
                        
                        if (Z < Z_bnd_lo || Z > Z_bnd_up)
                        {
                            is_first_order = true;
                        }
                    }
                    
                    if (Z_last < Z_bnd_lo || Z_last > Z_bnd_up)
                    {
                        is_first_order = true;
                    }
                    */
                    
                    break;
                }
                default:
                        TBOX_ERROR(d_object_name
                                   << ": "
                                   << "Unknown d_flow_model."
                                   << std::endl);
            }
            
            // Change the fluxes to first order fluxes if it is needed.
            
            if (is_first_order)
            {
                switch (d_flow_model)
                {
                    case SINGLE_SPECIES:
                    {
                        const double& rho_cell = Q[0][idx_cell];
                        
                        const double& E_cell = Q[1 + d_dim.getValue()][idx_cell];
                        
                        const double u_cell = Q[1][idx_cell]/rho_cell;
                        
                        std::vector<const double*> m_ptr_cell;
                        for (int di = 0; di < d_dim.getValue(); di++)
                        {
                            m_ptr_cell.push_back(&Q[1 + di][idx_cell]);
                        }
                        
                        const double p_cell = d_equation_of_state->getPressure(
                            &rho_cell,
                            m_ptr_cell,
                            &E_cell);
                        
                        const double c_cell = d_equation_of_state->getSoundSpeedWithPressure(
                            &rho_cell,
                            &p_cell);
                        
                        const double sp_x_cell = fabs(u_cell) + c_cell;
                        
                        // Correct the fluxes in the x direction.
                        {
                            // Compute indices of left and right cells.
                            const int idx_cell_L   = i - 1 + d_num_ghosts[0];
                            const int idx_cell_R   = i + 1 + d_num_ghosts[0];
                            
                            const double& rho_L = Q[0][idx_cell_L];
                            const double& rho_R = Q[0][idx_cell_R];
                            
                            const double& E_L = Q[1 + d_dim.getValue()][idx_cell_L];
                            const double& E_R = Q[1 + d_dim.getValue()][idx_cell_R];
                            
                            const double u_L = Q[1][idx_cell_L]/rho_L;
                            const double u_R = Q[1][idx_cell_R]/rho_R;
                            
                            std::vector<const double*> m_ptr_L;
                            for (int di = 0; di < d_dim.getValue(); di++)
                            {
                                m_ptr_L.push_back(&Q[1 + di][idx_cell_L]);
                            }
                            
                            std::vector<const double*> m_ptr_R;
                            for (int di = 0; di < d_dim.getValue(); di++)
                            {
                                m_ptr_R.push_back(&Q[1 + di][idx_cell_R]);
                            }
                            
                            const double p_L = d_equation_of_state->getPressure(
                                &rho_L,
                                m_ptr_L,
                                &E_L);
                            
                            const double p_R = d_equation_of_state->getPressure(
                                &rho_R,
                                m_ptr_R,
                                &E_R);
                            
                            const double c_L = d_equation_of_state->getSoundSpeedWithPressure(
                                &rho_L,
                                &p_L);
                            
                            const double c_R = d_equation_of_state->getSoundSpeedWithPressure(
                                &rho_R,
                                &p_R);
                            
                            const double sp_x_L = fabs(u_L) + c_L;
                            const double sp_x_R = fabs(u_R) + c_R;
                            
                            std::vector<double> F_x_cell;
                            std::vector<double> F_x_L;
                            std::vector<double> F_x_R;
                            
                            F_x_cell.push_back(rho_cell*u_cell);
                            F_x_cell.push_back(rho_cell*u_cell*u_cell + p_cell);
                            F_x_cell.push_back(u_cell*(E_cell + p_cell));
                            
                            F_x_L.push_back(rho_L*u_L);
                            F_x_L.push_back(rho_L*u_L*u_L + p_L);
                            F_x_L.push_back(u_L*(E_L + p_L));
                            
                            F_x_R.push_back(rho_R*u_R);
                            F_x_R.push_back(rho_R*u_R*u_R + p_R);
                            F_x_R.push_back(u_R*(E_R + p_R));
                            
                            const double alpha_x_L = fmax(sp_x_L, sp_x_cell);
                            const double alpha_x_R = fmax(sp_x_R, sp_x_cell);
                            
                            for (int ei = 0; ei < d_num_eqn; ei++)
                            {
                                double* F_x = convective_flux_intermediate->getPointer(0, ei);
                                F_x[idx_face_x - 1] = 0.5*dt*(F_x_L[ei] + F_x_cell[ei] -
                                    alpha_x_L*(Q[ei][idx_cell] - Q[ei][idx_cell_L]));
                                
                                F_x[idx_face_x]     = 0.5*dt*(F_x_cell[ei] + F_x_R[ei] -
                                    alpha_x_R*(Q[ei][idx_cell_R] - Q[ei][idx_cell]));
                            }
                        }
                        
                        break;
                    }
                    case FOUR_EQN_SHYUE:
                    {
                        const double& rho_cell = Q[0][idx_cell];
                        
                        const double& E_cell = Q[1 + d_dim.getValue()][idx_cell];
                        
                        const double u_cell = Q[1][idx_cell]/rho_cell;
                        
                        std::vector<const double*> m_ptr_cell;
                        for (int di = 0; di < d_dim.getValue(); di++)
                        {
                            m_ptr_cell.push_back(&Q[1 + di][idx_cell]);
                        }
                        
                        std::vector<const double*> Y_ptr_cell;
                        for (int si = 0; si < d_num_species; si++)
                        {
                            Y_ptr_cell.push_back(&Q[2 + d_dim.getValue() + si][idx_cell]);
                        }
                        
                        const double p_cell = d_equation_of_state->getPressureWithMassFraction(
                            &rho_cell,
                            m_ptr_cell,
                            &E_cell,
                            Y_ptr_cell);
                        
                        const double c_cell = d_equation_of_state->getSoundSpeedWithMassFractionAndPressure(
                            &rho_cell,
                            Y_ptr_cell,
                            &p_cell);
                        
                        const double sp_x_cell = fabs(u_cell) + c_cell;
                        
                        // Correct the fluxes in the x direction.
                        {
                            // Compute indices of left and right cells.
                            const int idx_cell_L   = i - 1 + d_num_ghosts[0];
                            const int idx_cell_R   = i + 1 + d_num_ghosts[0];
                            
                            const double& rho_L = Q[0][idx_cell_L];
                            const double& rho_R = Q[0][idx_cell_R];
                            
                            const double& E_L = Q[1 + d_dim.getValue()][idx_cell_L];
                            const double& E_R = Q[1 + d_dim.getValue()][idx_cell_R];
                            
                            const double u_L = Q[1][idx_cell_L]/rho_L;
                            const double u_R = Q[1][idx_cell_R]/rho_R;
                            
                            std::vector<const double*> m_ptr_L;
                            for (int di = 0; di < d_dim.getValue(); di++)
                            {
                                m_ptr_L.push_back(&Q[1 + di][idx_cell_L]);
                            }
                            
                            std::vector<const double*> m_ptr_R;
                            for (int di = 0; di < d_dim.getValue(); di++)
                            {
                                m_ptr_R.push_back(&Q[1 + di][idx_cell_R]);
                            }
                            
                            std::vector<const double*> Y_ptr_L;
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_ptr_L.push_back(&Q[2 + d_dim.getValue() + si][idx_cell_L]);
                            }
                            
                            std::vector<const double*> Y_ptr_R;
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_ptr_R.push_back(&Q[2 + d_dim.getValue() + si][idx_cell_R]);
                            }
                            
                            const double p_L = d_equation_of_state->getPressureWithMassFraction(
                                &rho_L,
                                m_ptr_L,
                                &E_L,
                                Y_ptr_L);
                            
                            const double p_R = d_equation_of_state->getPressureWithMassFraction(
                                &rho_R,
                                m_ptr_R,
                                &E_R,
                                Y_ptr_R);
                            
                            const double c_L = d_equation_of_state->getSoundSpeedWithMassFractionAndPressure(
                                &rho_L,
                                Y_ptr_L,
                                &p_L);
                            
                            const double c_R = d_equation_of_state->getSoundSpeedWithMassFractionAndPressure(
                                &rho_R,
                                Y_ptr_R,
                                &p_R);
                            
                            const double sp_x_L = fabs(u_L) + c_L;
                            const double sp_x_R = fabs(u_R) + c_R;
                            
                            std::vector<double> F_x_cell;
                            std::vector<double> F_x_L;
                            std::vector<double> F_x_R;
                            
                            F_x_cell.push_back(rho_cell*u_cell);
                            F_x_cell.push_back(rho_cell*u_cell*u_cell + p_cell);
                            F_x_cell.push_back(u_cell*(E_cell + p_cell));
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                F_x_cell.push_back(u_cell*(*Y_ptr_cell[si]));
                            }
                            
                            F_x_L.push_back(rho_L*u_L);
                            F_x_L.push_back(rho_L*u_L*u_L + p_L);
                            F_x_L.push_back(u_L*(E_L + p_L));
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                F_x_L.push_back(u_L*(*Y_ptr_L[si]));
                            }
                            
                            F_x_R.push_back(rho_R*u_R);
                            F_x_R.push_back(rho_R*u_R*u_R + p_R);
                            F_x_R.push_back(u_R*(E_R + p_R));
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                F_x_R.push_back(u_R*(*Y_ptr_R[si]));
                            }
                            
                            const double alpha_x_L = fmax(sp_x_L, sp_x_cell);
                            const double alpha_x_R = fmax(sp_x_R, sp_x_cell);
                            
                            for (int ei = 0; ei < d_num_eqn; ei++)
                            {
                                double* F_x = convective_flux_intermediate->getPointer(0, ei);
                                F_x[idx_face_x - 1] = 0.5*dt*(F_x_L[ei] + F_x_cell[ei] -
                                    alpha_x_L*(Q[ei][idx_cell] - Q[ei][idx_cell_L]));
                                
                                F_x[idx_face_x]     = 0.5*dt*(F_x_cell[ei] + F_x_R[ei] -
                                    alpha_x_R*(Q[ei][idx_cell_R] - Q[ei][idx_cell]));
                            }
                        }
                        
                        break;
                    }
                    case FIVE_EQN_ALLAIRE:
                    {
                        std::vector<const double*> Z_rho_ptr_cell;
                        for (int si = 0; si < d_num_species; si++)
                        {
                            Z_rho_ptr_cell.push_back(&Q[si][idx_cell]);
                        }
                        
                        const double rho_cell = d_equation_of_state->getTotalDensity(Z_rho_ptr_cell);
                        
                        const double& E_cell = Q[d_num_species + d_dim.getValue()][idx_cell];
                        
                        const double u_cell = Q[d_num_species][idx_cell]/rho_cell;
                        
                        std::vector<const double*> m_ptr_cell;
                        for (int di = 0; di < d_dim.getValue(); di++)
                        {
                            m_ptr_cell.push_back(&Q[d_num_species + di][idx_cell]);
                        }
                        
                        std::vector<const double*> Z_ptr_cell;
                        for (int si = 0; si < d_num_species; si++)
                        {
                            Z_ptr_cell.push_back(&Q[1 + d_num_species + d_dim.getValue() + si][idx_cell]);
                        }
                        
                        const double p_cell = d_equation_of_state->getPressureWithVolumeFraction(
                            &rho_cell,
                            m_ptr_cell,
                            &E_cell,
                            Z_ptr_cell);
                        
                        const double c_cell = d_equation_of_state->getSoundSpeedWithVolumeFractionAndPressure(
                            &rho_cell,
                            Z_ptr_cell,
                            &p_cell);
                        
                        const double sp_x_cell = fabs(u_cell) + c_cell;
                        
                        // Correct the fluxes in the x direction.
                        {
                            // Compute indices of left and right cells.
                            const int idx_cell_L   = i - 1 + d_num_ghosts[0];
                            const int idx_cell_R   = i + 1 + d_num_ghosts[0];
                            
                            std::vector<const double*> Z_rho_ptr_L;
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Z_rho_ptr_L.push_back(&Q[si][idx_cell_L]);
                            }
                            
                            std::vector<const double*> Z_rho_ptr_R;
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Z_rho_ptr_R.push_back(&Q[si][idx_cell_R]);
                            }
                            
                            const double& rho_L = d_equation_of_state->getTotalDensity(Z_rho_ptr_L);
                            const double& rho_R = d_equation_of_state->getTotalDensity(Z_rho_ptr_R);
                            
                            const double& E_L = Q[d_num_species + d_dim.getValue()][idx_cell_L];
                            const double& E_R = Q[d_num_species + d_dim.getValue()][idx_cell_R];
                            
                            const double u_L = Q[d_num_species][idx_cell_L]/rho_L;
                            const double u_R = Q[d_num_species][idx_cell_R]/rho_R;
                            
                            std::vector<const double*> m_ptr_L;
                            for (int di = 0; di < d_dim.getValue(); di++)
                            {
                                m_ptr_L.push_back(&Q[d_num_species + di][idx_cell_L]);
                            }
                            
                            std::vector<const double*> m_ptr_R;
                            for (int di = 0; di < d_dim.getValue(); di++)
                            {
                                m_ptr_R.push_back(&Q[d_num_species + di][idx_cell_R]);
                            }
                            
                            std::vector<const double*> Z_ptr_L;
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Z_ptr_L.push_back(&Q[1 + d_num_species + d_dim.getValue() + si][idx_cell_L]);
                            }
                            
                            std::vector<const double*> Z_ptr_R;
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Z_ptr_R.push_back(&Q[1 + d_num_species + d_dim.getValue() + si][idx_cell_R]);
                            }
                            
                            const double p_L = d_equation_of_state->getPressureWithVolumeFraction(
                                &rho_L,
                                m_ptr_L,
                                &E_L,
                                Z_ptr_L);
                            
                            const double p_R = d_equation_of_state->getPressureWithVolumeFraction(
                                &rho_R,
                                m_ptr_R,
                                &E_R,
                                Z_ptr_R);
                            
                            const double c_L = d_equation_of_state->getSoundSpeedWithVolumeFractionAndPressure(
                                &rho_L,
                                Z_ptr_L,
                                &p_L);
                            
                            const double c_R = d_equation_of_state->getSoundSpeedWithVolumeFractionAndPressure(
                                &rho_R,
                                Z_ptr_R,
                                &p_R);
                            
                            const double sp_x_L = fabs(u_L) + c_L;
                            const double sp_x_R = fabs(u_R) + c_R;
                            
                            std::vector<double> F_x_cell;
                            std::vector<double> F_x_L;
                            std::vector<double> F_x_R;
                            
                            for (int si = 0; si < d_num_species; si++)
                            {
                                F_x_cell.push_back(u_cell*(*Z_rho_ptr_cell[si]));
                            }
                            F_x_cell.push_back(rho_cell*u_cell*u_cell + p_cell);
                            F_x_cell.push_back(u_cell*(E_cell + p_cell));
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                F_x_cell.push_back(u_cell*(*Z_ptr_cell[si]));
                            }
                            
                            for (int si = 0; si < d_num_species; si++)
                            {
                                F_x_L.push_back(u_L*(*Z_rho_ptr_L[si]));
                            }
                            F_x_L.push_back(rho_L*u_L*u_L + p_L);
                            F_x_L.push_back(u_L*(E_L + p_L));
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                F_x_L.push_back(u_L*(*Z_ptr_L[si]));
                            }
                            
                            for (int si = 0; si < d_num_species; si++)
                            {
                                F_x_R.push_back(u_R*(*Z_rho_ptr_R[si]));
                            }
                            F_x_R.push_back(rho_R*u_R*u_R + p_R);
                            F_x_R.push_back(u_R*(E_R + p_R));
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                F_x_R.push_back(u_R*(*Z_ptr_R[si]));
                            }
                            
                            const double alpha_x_L = fmax(sp_x_L, sp_x_cell);
                            const double alpha_x_R = fmax(sp_x_R, sp_x_cell);
                            
                            for (int ei = 0; ei < d_num_eqn; ei++)
                            {
                                double* F_x = convective_flux_intermediate->getPointer(0, ei);
                                F_x[idx_face_x - 1] = 0.5*dt*(F_x_L[ei] + F_x_cell[ei] -
                                    alpha_x_L*(Q[ei][idx_cell] - Q[ei][idx_cell_L]));
                                
                                F_x[idx_face_x]     = 0.5*dt*(F_x_cell[ei] + F_x_R[ei] -
                                    alpha_x_R*(Q[ei][idx_cell_R] - Q[ei][idx_cell]));
                            }
                        }
                        
                        break;
                    }
                    default:
                        TBOX_ERROR(d_object_name
                                   << ": "
                                   << "Unknown d_flow_model."
                                   << std::endl);
                }
            }
        }
    } // if (d_dim == tbox::Dimension(1))
    else if (d_dim == tbox::Dimension(2))
    {
        for (int j = 0; j < interior_dims[1]; j++)
        {
            for (int i = 0; i < interior_dims[0]; i++)
            {
                // Compute indices of time-dependent data, flux and source.
                const int idx_cell   = (i + d_num_ghosts[0]) +
                    (j + d_num_ghosts[1])*ghostcell_dims[0];
                const int idx_source = i + j*interior_dims[0];
                const int idx_face_x = (i + 1) + j*(interior_dims[0] + 1);
                const int idx_face_y = (j + 1) + i*(interior_dims[1] + 1);
                
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    double* F_x_intermediate = convective_flux_intermediate->getPointer(0, ei);
                    double* F_y_intermediate = convective_flux_intermediate->getPointer(1, ei);
                    double* S_intermediate = source_intermediate->getPointer(ei);
                    
                    Q_tmp[ei] = Q[ei][idx_cell] + beta*
                        (-(F_x_intermediate[idx_face_x] - F_x_intermediate[idx_face_x - 1])/dx[0] -
                        (F_y_intermediate[idx_face_y] - F_y_intermediate[idx_face_y - 1])/dx[1] +
                        S_intermediate[idx_source]);
                }
                
                double rho = 0;
                double p   = 0;
                
                /*
                 * Check if the density or pressure are negative and if the mass or volume fractions
                 * are within the bound. If they are, switch to first order hyperbolic flux.
                 */
                
                switch (d_flow_model)
                {
                    case SINGLE_SPECIES:
                    {
                        std::vector<const double*> m_ptr;
                        for (int di = 0; di < d_dim.getValue(); di++)
                        {
                            m_ptr.push_back(&Q_tmp[1 + di]);
                        }
                        
                        rho = Q_tmp[0];
                        const double& E = Q_tmp[1 + d_dim.getValue()];
                        
                        p = d_equation_of_state->getPressure(
                            &rho,
                            m_ptr,
                            &E);
                        
                        break;
                    }
                    case FOUR_EQN_SHYUE:
                    {
                        std::vector<const double*> m_ptr;
                        for (int di = 0; di < d_dim.getValue(); di++)
                        {
                            m_ptr.push_back(&Q_tmp[1 + di]);
                        }
                        
                        std::vector<const double*> Y_ptr;
                        for (int si = 0; si < d_num_species - 1; si++)
                        {
                            Y_ptr.push_back(&Q_tmp[2 + d_dim.getValue() + si]);
                        }
                        
                        rho = Q_tmp[0];
                        const double& E = Q_tmp[1 + d_dim.getValue()];
                        
                        p = d_equation_of_state->getPressureWithMassFraction(
                            &rho,
                            m_ptr,
                            &E,
                            Y_ptr);
                        
                        break;
                    }
                    case FIVE_EQN_ALLAIRE:
                    {
                        std::vector<const double*> Z_rho_ptr;
                        for (int si = 0; si < d_num_species; si++)
                        {
                            Z_rho_ptr.push_back(&Q_tmp[si]);
                        }
                        
                        std::vector<const double*> m_ptr;
                        for (int di = 0; di < d_dim.getValue(); di++)
                        {
                            m_ptr.push_back(&Q_tmp[d_num_species + di]);
                        }
                        
                        std::vector<const double*> Z_ptr;
                        for (int si = 0; si < d_num_species - 1; si++)
                        {
                            Z_ptr.push_back(&Q_tmp[1 + d_num_species + d_dim.getValue() + si]);
                        }
                        
                        rho = d_equation_of_state->getTotalDensity(Z_rho_ptr);
                        const double& E = Q_tmp[d_num_species + d_dim.getValue()];
                        
                        p = d_equation_of_state->getPressureWithVolumeFraction(
                            &rho,
                            m_ptr,
                            &E,
                            Z_ptr);
                        
                        break;
                    }
                    default:
                        TBOX_ERROR(d_object_name
                                   << ": "
                                   << "Unknown d_flow_model."
                                   << std::endl);
                }
                
                bool is_first_order = false;
                
                // Check if the density or pressure are negative.
                
                if (rho < 0 || p < 0)
                {
                    is_first_order = true;
                }
                
                // Check if the mass or volume fractions are within the bound.
                
                switch (d_flow_model)
                {
                    case SINGLE_SPECIES:
                    {
                        break;
                    }
                    case FOUR_EQN_SHYUE:
                    {
                        double Y_last = 1.0;
                        
                        for (int si = 0; si < d_num_species - 1; si++)
                        {
                            double Y = Q_tmp[2 + d_dim.getValue() + si];
                            
                            Y_last -= Y;
                            
                            if (Y < Y_bnd_lo || Y > Y_bnd_up)
                            {
                                is_first_order = true;
                            }
                        }
                        
                        if (Y_last < Y_bnd_lo || Y_last > Y_bnd_up)
                        {
                            is_first_order = true;
                        }
                        
                        break;
                    }
                    case FIVE_EQN_ALLAIRE:
                    {
                        for (int si = 0; si < d_num_species; si++)
                        {
                            double Y = Q_tmp[si]/rho;
                            
                            if (Y < Y_bnd_lo || Y > Y_bnd_up)
                            {
                                is_first_order = true;
                                break;
                            }
                        }
                        
                        /*
                        double Z_last = 1.0;
                        
                        for (int si = 0; si < d_num_species - 1; si++)
                        {
                            double Z = Q_tmp[1 + d_num_species + d_dim.getValue() + si];
                            
                            Z_last -= Z;
                            
                            if (Z < Z_bnd_lo || Z > Z_bnd_up)
                            {
                                is_first_order = true;
                            }
                        }
                        
                        if (Z_last < Z_bnd_lo || Z_last > Z_bnd_up)
                        {
                            is_first_order = true;
                        }
                        */
                        
                        break;
                    }
                    default:
                            TBOX_ERROR(d_object_name
                                       << ": "
                                       << "Unknown d_flow_model."
                                       << std::endl);
                }
                
                // Change the fluxes to first order fluxes if it is needed.
                
                if (is_first_order)
                {
                    switch (d_flow_model)
                    {
                        case SINGLE_SPECIES:
                        {
                            const double& rho_cell = Q[0][idx_cell];
                            
                            const double& E_cell = Q[1 + d_dim.getValue()][idx_cell];
                            
                            const double u_cell = Q[1][idx_cell]/rho_cell;
                            
                            const double v_cell = Q[2][idx_cell]/rho_cell;
                            
                            std::vector<const double*> m_ptr_cell;
                            for (int di = 0; di < d_dim.getValue(); di++)
                            {
                                m_ptr_cell.push_back(&Q[1 + di][idx_cell]);
                            }
                            
                            const double p_cell = d_equation_of_state->getPressure(
                                &rho_cell,
                                m_ptr_cell,
                                &E_cell);
                            
                            const double c_cell = d_equation_of_state->getSoundSpeedWithPressure(
                                &rho_cell,
                                &p_cell);
                            
                            const double sp_x_cell = fabs(u_cell) + c_cell;
                            const double sp_y_cell = fabs(v_cell) + c_cell;
                            
                            // Correct the fluxes in the x direction.
                            {
                                // Compute indices of left and right cells.
                                int idx_cell_L   = (i - 1 + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                int idx_cell_R   = (i + 1 + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                const double& rho_L = Q[0][idx_cell_L];
                                const double& rho_R = Q[0][idx_cell_R];
                                
                                const double& E_L = Q[1 + d_dim.getValue()][idx_cell_L];
                                const double& E_R = Q[1 + d_dim.getValue()][idx_cell_R];
                                
                                const double u_L = Q[1][idx_cell_L]/rho_L;
                                const double u_R = Q[1][idx_cell_R]/rho_R;
                                
                                const double v_L = Q[2][idx_cell_L]/rho_L;
                                const double v_R = Q[2][idx_cell_R]/rho_R;
                                
                                std::vector<const double*> m_ptr_L;
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    m_ptr_L.push_back(&Q[1 + di][idx_cell_L]);
                                }
                                
                                std::vector<const double*> m_ptr_R;
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    m_ptr_R.push_back(&Q[1 + di][idx_cell_R]);
                                }
                                
                                const double p_L = d_equation_of_state->getPressure(
                                    &rho_L,
                                    m_ptr_L,
                                    &E_L);
                                
                                const double p_R = d_equation_of_state->getPressure(
                                    &rho_R,
                                    m_ptr_R,
                                    &E_R);
                                
                                const double c_L = d_equation_of_state->getSoundSpeedWithPressure(
                                    &rho_L,
                                    &p_L);
                                
                                const double c_R = d_equation_of_state->getSoundSpeedWithPressure(
                                    &rho_R,
                                    &p_R);
                                
                                const double sp_x_L = fabs(u_L) + c_L;
                                const double sp_x_R = fabs(u_R) + c_R;
                                
                                std::vector<double> F_x_cell;
                                std::vector<double> F_x_L;
                                std::vector<double> F_x_R;
                                
                                F_x_cell.push_back(rho_cell*u_cell);
                                F_x_cell.push_back(rho_cell*u_cell*u_cell + p_cell);
                                F_x_cell.push_back(rho_cell*u_cell*v_cell);
                                F_x_cell.push_back(u_cell*(E_cell + p_cell));
                                
                                F_x_L.push_back(rho_L*u_L);
                                F_x_L.push_back(rho_L*u_L*u_L + p_L);
                                F_x_L.push_back(rho_L*u_L*v_L);
                                F_x_L.push_back(u_L*(E_L + p_L));
                                
                                F_x_R.push_back(rho_R*u_R);
                                F_x_R.push_back(rho_R*u_R*u_R + p_R);
                                F_x_R.push_back(rho_R*u_R*v_R);
                                F_x_R.push_back(u_R*(E_R + p_R));
                                
                                const double alpha_x_L = fmax(sp_x_L, sp_x_cell);
                                const double alpha_x_R = fmax(sp_x_R, sp_x_cell);
                                
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    double* F_x = convective_flux_intermediate->getPointer(0, ei);
                                    F_x[idx_face_x - 1] = 0.5*dt*(F_x_L[ei] + F_x_cell[ei] -
                                        alpha_x_L*(Q[ei][idx_cell] - Q[ei][idx_cell_L]));
                                    
                                    F_x[idx_face_x]     = 0.5*dt*(F_x_cell[ei] + F_x_R[ei] -
                                        alpha_x_R*(Q[ei][idx_cell_R] - Q[ei][idx_cell]));
                                }
                            }
                            
                            // Correct the fluxes in the y direction.
                            {
                                int idx_cell_B   = (i + d_num_ghosts[0]) +
                                    (j - 1 + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                int idx_cell_T   = (i + d_num_ghosts[0]) +
                                    (j + 1 + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                const double& rho_B = Q[0][idx_cell_B];
                                const double& rho_T = Q[0][idx_cell_T];
                                
                                const double& E_B = Q[1 + d_dim.getValue()][idx_cell_B];
                                const double& E_T = Q[1 + d_dim.getValue()][idx_cell_T];
                                
                                const double u_B = Q[1][idx_cell_B]/rho_B;
                                const double u_T = Q[1][idx_cell_T]/rho_T;
                                
                                const double v_B = Q[2][idx_cell_B]/rho_B;
                                const double v_T = Q[2][idx_cell_T]/rho_T;
                                
                                std::vector<const double*> m_ptr_B;
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    m_ptr_B.push_back(&Q[1 + di][idx_cell_B]);
                                }
                                
                                std::vector<const double*> m_ptr_T;
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    m_ptr_T.push_back(&Q[1 + di][idx_cell_T]);
                                }
                                
                                const double p_B = d_equation_of_state->getPressure(
                                    &rho_B,
                                    m_ptr_B,
                                    &E_B);
                                
                                const double p_T = d_equation_of_state->getPressure(
                                    &rho_T,
                                    m_ptr_T,
                                    &E_T);
                                
                                const double c_B = d_equation_of_state->getSoundSpeedWithPressure(
                                    &rho_B,
                                    &p_B);
                                
                                const double c_T = d_equation_of_state->getSoundSpeedWithPressure(
                                    &rho_T,
                                    &p_T);
                                
                                const double sp_y_B = fabs(v_B) + c_B;
                                const double sp_y_T = fabs(v_T) + c_T;
                                
                                std::vector<double> F_y_cell;
                                std::vector<double> F_y_B;
                                std::vector<double> F_y_T;
                                
                                F_y_cell.push_back(rho_cell*v_cell);
                                F_y_cell.push_back(rho_cell*v_cell*u_cell);
                                F_y_cell.push_back(rho_cell*v_cell*v_cell + p_cell);
                                F_y_cell.push_back(v_cell*(E_cell + p_cell));
                                
                                F_y_B.push_back(rho_B*v_B);
                                F_y_B.push_back(rho_B*v_B*u_B);
                                F_y_B.push_back(rho_B*v_B*v_B + p_B);
                                F_y_B.push_back(v_B*(E_B + p_B));
                                
                                F_y_T.push_back(rho_T*v_T);
                                F_y_T.push_back(rho_T*v_T*u_T);
                                F_y_T.push_back(rho_T*v_T*v_T + p_T);
                                F_y_T.push_back(v_T*(E_T + p_T));
                                
                                const double alpha_y_B = fmax(sp_y_B, sp_y_cell);
                                const double alpha_y_T = fmax(sp_y_T, sp_y_cell);
                                
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    double* F_y = convective_flux_intermediate->getPointer(1, ei);
                                    F_y[idx_face_y - 1] = 0.5*dt*(F_y_B[ei] + F_y_cell[ei] -
                                        alpha_y_B*(Q[ei][idx_cell] - Q[ei][idx_cell_B]));
                                    
                                    F_y[idx_face_y]     = 0.5*dt*(F_y_cell[ei] + F_y_T[ei] -
                                        alpha_y_T*(Q[ei][idx_cell_T] - Q[ei][idx_cell]));
                                }
                            }
                            
                            break;
                        }
                        case FOUR_EQN_SHYUE:
                        {
                            const double& rho_cell = Q[0][idx_cell];
                            
                            const double& E_cell = Q[1 + d_dim.getValue()][idx_cell];
                            
                            const double u_cell = Q[1][idx_cell]/rho_cell;
                            
                            const double v_cell = Q[2][idx_cell]/rho_cell;
                            
                            std::vector<const double*> m_ptr_cell;
                            for (int di = 0; di < d_dim.getValue(); di++)
                            {
                                m_ptr_cell.push_back(&Q[1 + di][idx_cell]);
                            }
                            
                            std::vector<const double*> Y_ptr_cell;
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_ptr_cell.push_back(&Q[2 + d_dim.getValue() + si][idx_cell]);
                            }
                            
                            const double p_cell = d_equation_of_state->getPressureWithMassFraction(
                                &rho_cell,
                                m_ptr_cell,
                                &E_cell,
                                Y_ptr_cell);
                            
                            const double c_cell = d_equation_of_state->getSoundSpeedWithMassFractionAndPressure(
                                &rho_cell,
                                Y_ptr_cell,
                                &p_cell);
                            
                            const double sp_x_cell = fabs(u_cell) + c_cell;
                            const double sp_y_cell = fabs(v_cell) + c_cell;
                            
                            // Correct the fluxes in the x direction.
                            {
                                // Compute indices of left and right cells.
                                int idx_cell_L   = (i - 1 + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                int idx_cell_R   = (i + 1 + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                const double& rho_L = Q[0][idx_cell_L];
                                const double& rho_R = Q[0][idx_cell_R];
                                
                                const double& E_L = Q[1 + d_dim.getValue()][idx_cell_L];
                                const double& E_R = Q[1 + d_dim.getValue()][idx_cell_R];
                                
                                const double u_L = Q[1][idx_cell_L]/rho_L;
                                const double u_R = Q[1][idx_cell_R]/rho_R;
                                
                                const double v_L = Q[2][idx_cell_L]/rho_L;
                                const double v_R = Q[2][idx_cell_R]/rho_R;
                                
                                std::vector<const double*> m_ptr_L;
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    m_ptr_L.push_back(&Q[1 + di][idx_cell_L]);
                                }
                                
                                std::vector<const double*> m_ptr_R;
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    m_ptr_R.push_back(&Q[1 + di][idx_cell_R]);
                                }
                                
                                std::vector<const double*> Y_ptr_L;
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_ptr_L.push_back(&Q[2 + d_dim.getValue() + si][idx_cell_L]);
                                }
                                
                                std::vector<const double*> Y_ptr_R;
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_ptr_R.push_back(&Q[2 + d_dim.getValue() + si][idx_cell_R]);
                                }
                                
                                const double p_L = d_equation_of_state->getPressureWithMassFraction(
                                    &rho_L,
                                    m_ptr_L,
                                    &E_L,
                                    Y_ptr_L);
                                
                                const double p_R = d_equation_of_state->getPressureWithMassFraction(
                                    &rho_R,
                                    m_ptr_R,
                                    &E_R,
                                    Y_ptr_R);
                                
                                const double c_L = d_equation_of_state->getSoundSpeedWithMassFractionAndPressure(
                                    &rho_L,
                                    Y_ptr_L,
                                    &p_L);
                                
                                const double c_R = d_equation_of_state->getSoundSpeedWithMassFractionAndPressure(
                                    &rho_R,
                                    Y_ptr_R,
                                    &p_R);
                                
                                const double sp_x_L = fabs(u_L) + c_L;
                                const double sp_x_R = fabs(u_R) + c_R;
                                
                                std::vector<double> F_x_cell;
                                std::vector<double> F_x_L;
                                std::vector<double> F_x_R;
                                
                                F_x_cell.push_back(rho_cell*u_cell);
                                F_x_cell.push_back(rho_cell*u_cell*u_cell + p_cell);
                                F_x_cell.push_back(rho_cell*u_cell*v_cell);
                                F_x_cell.push_back(u_cell*(E_cell + p_cell));
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    F_x_cell.push_back(u_cell*(*Y_ptr_cell[si]));
                                }
                                
                                F_x_L.push_back(rho_L*u_L);
                                F_x_L.push_back(rho_L*u_L*u_L + p_L);
                                F_x_L.push_back(rho_L*u_L*v_L);
                                F_x_L.push_back(u_L*(E_L + p_L));
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    F_x_L.push_back(u_L*(*Y_ptr_L[si]));
                                }
                                
                                F_x_R.push_back(rho_R*u_R);
                                F_x_R.push_back(rho_R*u_R*u_R + p_R);
                                F_x_R.push_back(rho_R*u_R*v_R);
                                F_x_R.push_back(u_R*(E_R + p_R));
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    F_x_R.push_back(u_R*(*Y_ptr_R[si]));
                                }
                                
                                const double alpha_x_L = fmax(sp_x_L, sp_x_cell);
                                const double alpha_x_R = fmax(sp_x_R, sp_x_cell);
                                
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    double* F_x = convective_flux_intermediate->getPointer(0, ei);
                                    F_x[idx_face_x - 1] = 0.5*dt*(F_x_L[ei] + F_x_cell[ei] -
                                        alpha_x_L*(Q[ei][idx_cell] - Q[ei][idx_cell_L]));
                                    
                                    F_x[idx_face_x]     = 0.5*dt*(F_x_cell[ei] + F_x_R[ei] -
                                        alpha_x_R*(Q[ei][idx_cell_R] - Q[ei][idx_cell]));
                                }
                            }
                            
                            // Correct the fluxes in the x direction.
                            {
                                int idx_cell_B   = (i + d_num_ghosts[0]) +
                                    (j - 1 + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                int idx_cell_T   = (i + d_num_ghosts[0]) +
                                    (j + 1 + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                const double& rho_B = Q[0][idx_cell_B];
                                const double& rho_T = Q[0][idx_cell_T];
                                
                                const double& E_B = Q[1 + d_dim.getValue()][idx_cell_B];
                                const double& E_T = Q[1 + d_dim.getValue()][idx_cell_T];
                                
                                const double u_B = Q[1][idx_cell_B]/rho_B;
                                const double u_T = Q[1][idx_cell_T]/rho_T;
                                
                                const double v_B = Q[2][idx_cell_B]/rho_B;
                                const double v_T = Q[2][idx_cell_T]/rho_T;
                                
                                std::vector<const double*> m_ptr_B;
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    m_ptr_B.push_back(&Q[1 + di][idx_cell_B]);
                                }
                                
                                std::vector<const double*> m_ptr_T;
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    m_ptr_T.push_back(&Q[1 + di][idx_cell_T]);
                                }
                                
                                std::vector<const double*> Y_ptr_B;
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_ptr_B.push_back(&Q[2 + d_dim.getValue() + si][idx_cell_B]);
                                }
                                
                                std::vector<const double*> Y_ptr_T;
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_ptr_T.push_back(&Q[2 + d_dim.getValue() + si][idx_cell_T]);
                                }
                                
                                const double p_B = d_equation_of_state->getPressureWithMassFraction(
                                    &rho_B,
                                    m_ptr_B,
                                    &E_B,
                                    Y_ptr_B);
                                
                                const double p_T = d_equation_of_state->getPressureWithMassFraction(
                                    &rho_T,
                                    m_ptr_T,
                                    &E_T,
                                    Y_ptr_T);
                                
                                const double c_B = d_equation_of_state->getSoundSpeedWithMassFractionAndPressure(
                                    &rho_B,
                                    Y_ptr_B,
                                    &p_B);
                                
                                const double c_T = d_equation_of_state->getSoundSpeedWithMassFractionAndPressure(
                                    &rho_T,
                                    Y_ptr_T,
                                    &p_T);
                                
                                const double sp_y_B = fabs(v_B) + c_B;
                                const double sp_y_T = fabs(v_T) + c_T;
                                
                                std::vector<double> F_y_cell;
                                std::vector<double> F_y_B;
                                std::vector<double> F_y_T;
                                
                                F_y_cell.push_back(rho_cell*v_cell);
                                F_y_cell.push_back(rho_cell*v_cell*u_cell);
                                F_y_cell.push_back(rho_cell*v_cell*v_cell + p_cell);
                                F_y_cell.push_back(v_cell*(E_cell + p_cell));
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    F_y_cell.push_back(v_cell*(*Y_ptr_cell[si]));
                                }
                                
                                F_y_B.push_back(rho_B*v_B);
                                F_y_B.push_back(rho_B*v_B*u_B);
                                F_y_B.push_back(rho_B*v_B*v_B + p_B);
                                F_y_B.push_back(v_B*(E_B + p_B));
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    F_y_B.push_back(v_B*(*Y_ptr_B[si]));
                                }
                                
                                F_y_T.push_back(rho_T*v_T);
                                F_y_T.push_back(rho_T*v_T*u_T);
                                F_y_T.push_back(rho_T*v_T*v_T + p_T);
                                F_y_T.push_back(v_T*(E_T + p_T));
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    F_y_T.push_back(v_T*(*Y_ptr_T[si]));
                                }
                                
                                const double alpha_y_B = fmax(sp_y_B, sp_y_cell);
                                const double alpha_y_T = fmax(sp_y_T, sp_y_cell);
                                
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    double* F_y = convective_flux_intermediate->getPointer(1, ei);
                                    F_y[idx_face_y - 1] = 0.5*dt*(F_y_B[ei] + F_y_cell[ei] -
                                        alpha_y_B*(Q[ei][idx_cell] - Q[ei][idx_cell_B]));
                                    
                                    F_y[idx_face_y]     = 0.5*dt*(F_y_cell[ei] + F_y_T[ei] -
                                        alpha_y_T*(Q[ei][idx_cell_T] - Q[ei][idx_cell]));
                                }
                            }
                            
                            break;
                        }
                        case FIVE_EQN_ALLAIRE:
                        {
                            std::vector<const double*> Z_rho_ptr_cell;
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Z_rho_ptr_cell.push_back(&Q[si][idx_cell]);
                            }
                            
                            const double rho_cell = d_equation_of_state->getTotalDensity(Z_rho_ptr_cell);
                            
                            const double& E_cell = Q[d_num_species + d_dim.getValue()][idx_cell];
                            
                            const double u_cell = Q[d_num_species][idx_cell]/rho_cell;
                            
                            const double v_cell = Q[d_num_species + 1][idx_cell]/rho_cell;
                            
                            std::vector<const double*> m_ptr_cell;
                            for (int di = 0; di < d_dim.getValue(); di++)
                            {
                                m_ptr_cell.push_back(&Q[d_num_species + di][idx_cell]);
                            }
                            
                            std::vector<const double*> Z_ptr_cell;
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Z_ptr_cell.push_back(&Q[1 + d_num_species + d_dim.getValue() + si][idx_cell]);
                            }
                            
                            const double p_cell = d_equation_of_state->getPressureWithVolumeFraction(
                                &rho_cell,
                                m_ptr_cell,
                                &E_cell,
                                Z_ptr_cell);
                            
                            const double c_cell = d_equation_of_state->getSoundSpeedWithVolumeFractionAndPressure(
                                &rho_cell,
                                Z_ptr_cell,
                                &p_cell);
                            
                            const double sp_x_cell = fabs(u_cell) + c_cell;
                            const double sp_y_cell = fabs(v_cell) + c_cell;
                            
                            // Correct the fluxes in the x direction.
                            {
                                // Compute indices of left and right cells.
                                int idx_cell_L   = (i - 1 + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                int idx_cell_R   = (i + 1 + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                std::vector<const double*> Z_rho_ptr_L;
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Z_rho_ptr_L.push_back(&Q[si][idx_cell_L]);
                                }
                                
                                std::vector<const double*> Z_rho_ptr_R;
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Z_rho_ptr_R.push_back(&Q[si][idx_cell_R]);
                                }
                                
                                const double& rho_L = d_equation_of_state->getTotalDensity(Z_rho_ptr_L);
                                const double& rho_R = d_equation_of_state->getTotalDensity(Z_rho_ptr_R);
                                
                                const double& E_L = Q[d_num_species + d_dim.getValue()][idx_cell_L];
                                const double& E_R = Q[d_num_species + d_dim.getValue()][idx_cell_R];
                                
                                const double u_L = Q[d_num_species][idx_cell_L]/rho_L;
                                const double u_R = Q[d_num_species][idx_cell_R]/rho_R;
                                
                                const double v_L = Q[d_num_species + 1][idx_cell_L]/rho_L;
                                const double v_R = Q[d_num_species + 1][idx_cell_R]/rho_R;
                                
                                std::vector<const double*> m_ptr_L;
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    m_ptr_L.push_back(&Q[d_num_species + di][idx_cell_L]);
                                }
                                
                                std::vector<const double*> m_ptr_R;
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    m_ptr_R.push_back(&Q[d_num_species + di][idx_cell_R]);
                                }
                                
                                std::vector<const double*> Z_ptr_L;
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Z_ptr_L.push_back(&Q[1 + d_num_species + d_dim.getValue() + si][idx_cell_L]);
                                }
                                
                                std::vector<const double*> Z_ptr_R;
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Z_ptr_R.push_back(&Q[1 + d_num_species + d_dim.getValue() + si][idx_cell_R]);
                                }
                                
                                const double p_L = d_equation_of_state->getPressureWithVolumeFraction(
                                    &rho_L,
                                    m_ptr_L,
                                    &E_L,
                                    Z_ptr_L);
                                
                                const double p_R = d_equation_of_state->getPressureWithVolumeFraction(
                                    &rho_R,
                                    m_ptr_R,
                                    &E_R,
                                    Z_ptr_R);
                                
                                const double c_L = d_equation_of_state->getSoundSpeedWithVolumeFractionAndPressure(
                                    &rho_L,
                                    Z_ptr_L,
                                    &p_L);
                                
                                const double c_R = d_equation_of_state->getSoundSpeedWithVolumeFractionAndPressure(
                                    &rho_R,
                                    Z_ptr_R,
                                    &p_R);
                                
                                const double sp_x_L = fabs(u_L) + c_L;
                                const double sp_x_R = fabs(u_R) + c_R;
                                
                                std::vector<double> F_x_cell;
                                std::vector<double> F_x_L;
                                std::vector<double> F_x_R;
                                
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    F_x_cell.push_back(u_cell*(*Z_rho_ptr_cell[si]));
                                }
                                F_x_cell.push_back(rho_cell*u_cell*u_cell + p_cell);
                                F_x_cell.push_back(rho_cell*u_cell*v_cell);
                                F_x_cell.push_back(u_cell*(E_cell + p_cell));
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    F_x_cell.push_back(u_cell*(*Z_ptr_cell[si]));
                                }
                                
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    F_x_L.push_back(u_L*(*Z_rho_ptr_L[si]));
                                }
                                F_x_L.push_back(rho_L*u_L*u_L + p_L);
                                F_x_L.push_back(rho_L*u_L*v_L);
                                F_x_L.push_back(u_L*(E_L + p_L));
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    F_x_L.push_back(u_L*(*Z_ptr_L[si]));
                                }
                                
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    F_x_R.push_back(u_R*(*Z_rho_ptr_R[si]));
                                }
                                F_x_R.push_back(rho_R*u_R*u_R + p_R);
                                F_x_R.push_back(rho_R*u_R*v_R);
                                F_x_R.push_back(u_R*(E_R + p_R));
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    F_x_R.push_back(u_R*(*Z_ptr_R[si]));
                                }
                                
                                const double alpha_x_L = fmax(sp_x_L, sp_x_cell);
                                const double alpha_x_R = fmax(sp_x_R, sp_x_cell);
                                
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    double* F_x = convective_flux_intermediate->getPointer(0, ei);
                                    F_x[idx_face_x - 1] = 0.5*dt*(F_x_L[ei] + F_x_cell[ei] -
                                        alpha_x_L*(Q[ei][idx_cell] - Q[ei][idx_cell_L]));
                                    
                                    F_x[idx_face_x]     = 0.5*dt*(F_x_cell[ei] + F_x_R[ei] -
                                        alpha_x_R*(Q[ei][idx_cell_R] - Q[ei][idx_cell]));
                                }
                            }
                            
                            // Correct the fluxes in the y direction.
                            {
                                int idx_cell_B   = (i + d_num_ghosts[0]) +
                                    (j - 1 + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                int idx_cell_T   = (i + d_num_ghosts[0]) +
                                    (j + 1 + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                std::vector<const double*> Z_rho_ptr_B;
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Z_rho_ptr_B.push_back(&Q[si][idx_cell_B]);
                                }
                                
                                std::vector<const double*> Z_rho_ptr_T;
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Z_rho_ptr_T.push_back(&Q[si][idx_cell_T]);
                                }
                                
                                const double& rho_B = d_equation_of_state->getTotalDensity(Z_rho_ptr_B);
                                const double& rho_T = d_equation_of_state->getTotalDensity(Z_rho_ptr_T);
                                
                                const double& E_B = Q[d_num_species + d_dim.getValue()][idx_cell_B];
                                const double& E_T = Q[d_num_species + d_dim.getValue()][idx_cell_T];
                                
                                const double u_B = Q[d_num_species][idx_cell_B]/rho_B;
                                const double u_T = Q[d_num_species][idx_cell_T]/rho_T;
                                
                                const double v_B = Q[d_num_species + 1][idx_cell_B]/rho_B;
                                const double v_T = Q[d_num_species + 1][idx_cell_T]/rho_T;
                                
                                std::vector<const double*> m_ptr_B;
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    m_ptr_B.push_back(&Q[d_num_species + di][idx_cell_B]);
                                }
                                
                                std::vector<const double*> m_ptr_T;
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    m_ptr_T.push_back(&Q[d_num_species + di][idx_cell_T]);
                                }
                                
                                std::vector<const double*> Z_ptr_B;
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Z_ptr_B.push_back(&Q[1 + d_num_species + d_dim.getValue() + si][idx_cell_B]);
                                }
                                
                                std::vector<const double*> Z_ptr_T;
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Z_ptr_T.push_back(&Q[1 + d_num_species + d_dim.getValue() + si][idx_cell_T]);
                                }
                                
                                const double p_B = d_equation_of_state->getPressureWithVolumeFraction(
                                    &rho_B,
                                    m_ptr_B,
                                    &E_B,
                                    Z_ptr_B);
                                
                                const double p_T = d_equation_of_state->getPressureWithVolumeFraction(
                                    &rho_T,
                                    m_ptr_T,
                                    &E_T,
                                    Z_ptr_T);
                                
                                const double c_B = d_equation_of_state->getSoundSpeedWithVolumeFractionAndPressure(
                                    &rho_B,
                                    Z_ptr_B,
                                    &p_B);
                                
                                const double c_T = d_equation_of_state->getSoundSpeedWithVolumeFractionAndPressure(
                                    &rho_T,
                                    Z_ptr_T,
                                    &p_T);
                                
                                const double sp_y_B = fabs(v_B) + c_B;
                                const double sp_y_T = fabs(v_T) + c_T;
                                
                                std::vector<double> F_y_cell;
                                std::vector<double> F_y_B;
                                std::vector<double> F_y_T;
                                
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    F_y_cell.push_back(v_cell*(*Z_rho_ptr_cell[si]));
                                }
                                F_y_cell.push_back(rho_cell*v_cell*u_cell);
                                F_y_cell.push_back(rho_cell*v_cell*v_cell + p_cell);
                                F_y_cell.push_back(v_cell*(E_cell + p_cell));
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    F_y_cell.push_back(v_cell*(*Z_ptr_cell[si]));
                                }
                                
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    F_y_B.push_back(v_B*(*Z_rho_ptr_B[si]));
                                }
                                F_y_B.push_back(rho_B*v_B*u_B);
                                F_y_B.push_back(rho_B*v_B*v_B + p_B);
                                F_y_B.push_back(v_B*(E_B + p_B));
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    F_y_B.push_back(v_B*(*Z_ptr_B[si]));
                                }
                                
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    F_y_T.push_back(v_T*(*Z_rho_ptr_T[si]));
                                }
                                F_y_T.push_back(rho_T*v_T*u_T);
                                F_y_T.push_back(rho_T*v_T*v_T + p_T);
                                F_y_T.push_back(v_T*(E_T + p_T));
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    F_y_T.push_back(v_T*(*Z_ptr_T[si]));
                                }
                                
                                const double alpha_y_B = fmax(sp_y_B, sp_y_cell);
                                const double alpha_y_T = fmax(sp_y_T, sp_y_cell);
                                
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    double* F_y = convective_flux_intermediate->getPointer(1, ei);
                                    F_y[idx_face_y - 1] = 0.5*dt*(F_y_B[ei] + F_y_cell[ei] -
                                        alpha_y_B*(Q[ei][idx_cell] - Q[ei][idx_cell_B]));
                                    
                                    F_y[idx_face_y]     = 0.5*dt*(F_y_cell[ei] + F_y_T[ei] -
                                        alpha_y_T*(Q[ei][idx_cell_T] - Q[ei][idx_cell]));
                                }
                            }
                            
                            break;
                        }
                        default:
                            TBOX_ERROR(d_object_name
                                       << ": "
                                       << "Unknown d_flow_model."
                                       << std::endl);
                    }
                }
            }
        }
    } // if (d_dim == tbox::Dimension(2))
    else if (d_dim == tbox::Dimension(3))
    {
        for (int k = 0; k < interior_dims[2]; k++)
        {
            for (int j = 0; j < interior_dims[1]; j++)
            {
                for (int i = 0; i < interior_dims[0]; i++)
                {
                    // Compute indices of time-dependent data, flux and source.
                    const int idx_cell = (i + d_num_ghosts[0]) +
                        (j - 1 + d_num_ghosts[1])*ghostcell_dims[0] +
                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                    
                    const int idx_source = i +
                        j*interior_dims[0] +
                        k*interior_dims[0]*interior_dims[1];
                    
                    const int idx_face_x = (i + 1) +
                        j*(interior_dims[0] + 1) +
                        k*(interior_dims[0] + 1)*interior_dims[1];
                    
                    const int idx_face_y = (j + 1) +
                        k*(interior_dims[1] + 1) +
                        i*(interior_dims[1] + 1)*interior_dims[2];
                    
                    const int idx_face_z = (k + 1) +
                        i*(interior_dims[2] + 1) +
                        j*(interior_dims[2] + 1)*interior_dims[0];
                    
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        double* F_x_intermediate = convective_flux_intermediate->getPointer(0, ei);
                        double* F_y_intermediate = convective_flux_intermediate->getPointer(1, ei);
                        double* F_z_intermediate = convective_flux_intermediate->getPointer(2, ei);
                        double* S_intermediate = source_intermediate->getPointer(ei);
                        
                        Q_tmp[ei] = Q[ei][idx_cell] + beta*
                            (-(F_x_intermediate[idx_face_x] - F_x_intermediate[idx_face_x - 1])/dx[0] -
                            (F_y_intermediate[idx_face_y] - F_y_intermediate[idx_face_y - 1])/dx[1] -
                            (F_z_intermediate[idx_face_z] - F_z_intermediate[idx_face_z - 1])/dx[2] +
                            S_intermediate[idx_source]);
                    }
                    
                    double rho = 0;
                    double p   = 0;
                    
                    /*
                     * Check if the density or pressure are negative and if the mass or volume fractions
                     * are within the bound. If they are, switch to first order hyperbolic flux.
                     */
                    
                    switch (d_flow_model)
                    {
                        case SINGLE_SPECIES:
                        {
                            std::vector<const double*> m_ptr;
                            for (int di = 0; di < d_dim.getValue(); di++)
                            {
                                m_ptr.push_back(&Q_tmp[1 + di]);
                            }
                            
                            rho = Q_tmp[0];
                            const double& E = Q_tmp[1 + d_dim.getValue()];
                            
                            p = d_equation_of_state->getPressure(
                                &rho,
                                m_ptr,
                                &E);
                            
                            break;
                        }
                        case FOUR_EQN_SHYUE:
                        {
                            std::vector<const double*> m_ptr;
                            for (int di = 0; di < d_dim.getValue(); di++)
                            {
                                m_ptr.push_back(&Q_tmp[1 + di]);
                            }
                            
                            std::vector<const double*> Y_ptr;
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                Y_ptr.push_back(&Q_tmp[2 + d_dim.getValue() + si]);
                            }
                            
                            rho = Q_tmp[0];
                            const double& E = Q_tmp[1 + d_dim.getValue()];
                            
                            p = d_equation_of_state->getPressureWithMassFraction(
                                &rho,
                                m_ptr,
                                &E,
                                Y_ptr);
                            
                            break;
                        }
                        case FIVE_EQN_ALLAIRE:
                        {
                            std::vector<const double*> Z_rho_ptr;
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Z_rho_ptr.push_back(&Q_tmp[si]);
                            }
                            
                            std::vector<const double*> m_ptr;
                            for (int di = 0; di < d_dim.getValue(); di++)
                            {
                                m_ptr.push_back(&Q_tmp[d_num_species + di]);
                            }
                            
                            std::vector<const double*> Z_ptr;
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                Z_ptr.push_back(&Q_tmp[1 + d_num_species + d_dim.getValue() + si]);
                            }
                            
                            rho = d_equation_of_state->getTotalDensity(Z_rho_ptr);
                            const double& E = Q_tmp[d_num_species + d_dim.getValue()];
                            
                            p = d_equation_of_state->getPressureWithVolumeFraction(
                                &rho,
                                m_ptr,
                                &E,
                                Z_ptr);
                            
                            break;
                        }
                        default:
                        TBOX_ERROR(d_object_name
                                   << ": "
                                   << "Unknown d_flow_model."
                                   << std::endl);
                    }
                    
                    bool is_first_order = false;
                    
                    // Check if the density or pressure are negative.
                    
                    if (rho < 0 || p < 0)
                    {
                        is_first_order = true;
                    }
                    
                    // Check if the mass or volume fractions are within the bound.
                    
                    switch (d_flow_model)
                    {
                        case SINGLE_SPECIES:
                        {
                            break;
                        }
                        case FOUR_EQN_SHYUE:
                        {
                            double Y_last = 1.0;
                            
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                double Y = Q_tmp[2 + d_dim.getValue() + si];
                                
                                Y_last -= Y;
                                
                                if (Y < Y_bnd_lo || Y > Y_bnd_up)
                                {
                                    is_first_order = true;
                                }
                            }
                            
                            if (Y_last < Y_bnd_lo || Y_last > Y_bnd_up)
                            {
                                is_first_order = true;
                            }
                            
                            break;
                        }
                        case FIVE_EQN_ALLAIRE:
                        {
                            for (int si = 0; si < d_num_species; si++)
                            {
                                double Y = Q_tmp[si]/rho;
                                
                                if (Y < Y_bnd_lo || Y > Y_bnd_up)
                                {
                                    is_first_order = true;
                                    break;
                                }
                            }
                            
                            /*
                            double Z_last = 1.0;
                            
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                double Z = Q_tmp[1 + d_num_species + d_dim.getValue() + si];
                                
                                Z_last -= Z;
                                
                                if (Z < Z_bnd_lo || Z > Z_bnd_up)
                                {
                                    is_first_order = true;
                                }
                            }
                            
                            if (Z_last < Z_bnd_lo || Z_last > Z_bnd_up)
                            {
                                is_first_order = true;
                            }
                            */
                            
                            break;
                        }
                        default:
                                TBOX_ERROR(d_object_name
                                           << ": "
                                           << "Unknown d_flow_model."
                                           << std::endl);
                    }
                    
                    // Change the fluxes to first order fluxes if it is needed.
                    
                    if (is_first_order)
                    {
                        switch (d_flow_model)
                        {
                            case SINGLE_SPECIES:
                            {
                                const double& rho_cell = Q[0][idx_cell];
                            
                                const double& E_cell = Q[1 + d_dim.getValue()][idx_cell];
                                
                                const double u_cell = Q[1][idx_cell]/rho_cell;
                                
                                const double v_cell = Q[2][idx_cell]/rho_cell;
                                
                                const double w_cell = Q[3][idx_cell]/rho_cell;
                                
                                std::vector<const double*> m_ptr_cell;
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    m_ptr_cell.push_back(&Q[1 + di][idx_cell]);
                                }
                                
                                const double p_cell = d_equation_of_state->getPressure(
                                    &rho_cell,
                                    m_ptr_cell,
                                    &E_cell);
                                
                                const double c_cell = d_equation_of_state->getSoundSpeedWithPressure(
                                    &rho_cell,
                                    &p_cell);
                                
                                const double sp_x_cell = fabs(u_cell) + c_cell;
                                const double sp_y_cell = fabs(v_cell) + c_cell;
                                const double sp_z_cell = fabs(w_cell) + c_cell;
                                
                                // Correct the fluxes in the x direction.
                                {
                                    // Compute indices of left and right cells.
                                    const int idx_cell_L   = (i - 1 + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_cell_R   = (i + 1 + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const double& rho_L = Q[0][idx_cell_L];
                                    const double& rho_R = Q[0][idx_cell_R];
                                    
                                    const double& E_L = Q[1 + d_dim.getValue()][idx_cell_L];
                                    const double& E_R = Q[1 + d_dim.getValue()][idx_cell_R];
                                    
                                    const double u_L = Q[1][idx_cell_L]/rho_L;
                                    const double u_R = Q[1][idx_cell_R]/rho_R;
                                    
                                    const double v_L = Q[2][idx_cell_L]/rho_L;
                                    const double v_R = Q[2][idx_cell_R]/rho_R;
                                    
                                    const double w_L = Q[3][idx_cell_L]/rho_L;
                                    const double w_R = Q[3][idx_cell_R]/rho_R;
                                    
                                    std::vector<const double*> m_ptr_L;
                                    for (int di = 0; di < d_dim.getValue(); di++)
                                    {
                                        m_ptr_L.push_back(&Q[1 + di][idx_cell_L]);
                                    }
                                    
                                    std::vector<const double*> m_ptr_R;
                                    for (int di = 0; di < d_dim.getValue(); di++)
                                    {
                                        m_ptr_R.push_back(&Q[1 + di][idx_cell_R]);
                                    }
                                    
                                    const double p_L = d_equation_of_state->getPressure(
                                        &rho_L,
                                        m_ptr_L,
                                        &E_L);
                                    
                                    const double p_R = d_equation_of_state->getPressure(
                                        &rho_R,
                                        m_ptr_R,
                                        &E_R);
                                    
                                    const double c_L = d_equation_of_state->getSoundSpeedWithPressure(
                                        &rho_L,
                                        &p_L);
                                    
                                    const double c_R = d_equation_of_state->getSoundSpeedWithPressure(
                                        &rho_R,
                                        &p_R);
                                    
                                    const double sp_x_L = fabs(u_L) + c_L;
                                    const double sp_x_R = fabs(u_R) + c_R;
                                    
                                    std::vector<double> F_x_cell;
                                    std::vector<double> F_x_L;
                                    std::vector<double> F_x_R;
                                    
                                    F_x_cell.push_back(rho_cell*u_cell);
                                    F_x_cell.push_back(rho_cell*u_cell*u_cell + p_cell);
                                    F_x_cell.push_back(rho_cell*u_cell*v_cell);
                                    F_x_cell.push_back(rho_cell*u_cell*w_cell);
                                    F_x_cell.push_back(u_cell*(E_cell + p_cell));
                                    
                                    F_x_L.push_back(rho_L*u_L);
                                    F_x_L.push_back(rho_L*u_L*u_L + p_L);
                                    F_x_L.push_back(rho_L*u_L*v_L);
                                    F_x_L.push_back(rho_L*u_L*w_L);
                                    F_x_L.push_back(u_L*(E_L + p_L));
                                    
                                    F_x_R.push_back(rho_R*u_R);
                                    F_x_R.push_back(rho_R*u_R*u_R + p_R);
                                    F_x_R.push_back(rho_R*u_R*v_R);
                                    F_x_R.push_back(rho_R*u_R*w_R);
                                    F_x_R.push_back(u_R*(E_R + p_R));
                                    
                                    const double alpha_x_L = fmax(sp_x_L, sp_x_cell);
                                    const double alpha_x_R = fmax(sp_x_R, sp_x_cell);
                                    
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        double* F_x = convective_flux_intermediate->getPointer(0, ei);
                                        F_x[idx_face_x - 1] = 0.5*dt*(F_x_L[ei] + F_x_cell[ei] -
                                            alpha_x_L*(Q[ei][idx_cell] - Q[ei][idx_cell_L]));
                                        
                                        F_x[idx_face_x]     = 0.5*dt*(F_x_cell[ei] + F_x_R[ei] -
                                            alpha_x_R*(Q[ei][idx_cell_R] - Q[ei][idx_cell]));
                                    }
                                }
                                
                                // Correct the fluxes in the y direction.
                                {
                                    const int idx_cell_B   = (i + d_num_ghosts[0]) +
                                        (j - 1 + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_cell_T   = (i + d_num_ghosts[0]) +
                                        (j + 1 + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const double& rho_B = Q[0][idx_cell_B];
                                    const double& rho_T = Q[0][idx_cell_T];
                                    
                                    const double& E_B = Q[1 + d_dim.getValue()][idx_cell_B];
                                    const double& E_T = Q[1 + d_dim.getValue()][idx_cell_T];
                                    
                                    const double u_B = Q[1][idx_cell_B]/rho_B;
                                    const double u_T = Q[1][idx_cell_T]/rho_T;
                                    
                                    const double v_B = Q[2][idx_cell_B]/rho_B;
                                    const double v_T = Q[2][idx_cell_T]/rho_T;
                                    
                                    const double w_B = Q[3][idx_cell_B]/rho_B;
                                    const double w_T = Q[3][idx_cell_T]/rho_T;
                                    
                                    std::vector<const double*> m_ptr_B;
                                    for (int di = 0; di < d_dim.getValue(); di++)
                                    {
                                        m_ptr_B.push_back(&Q[1 + di][idx_cell_B]);
                                    }
                                    
                                    std::vector<const double*> m_ptr_T;
                                    for (int di = 0; di < d_dim.getValue(); di++)
                                    {
                                        m_ptr_T.push_back(&Q[1 + di][idx_cell_T]);
                                    }
                                    
                                    const double p_B = d_equation_of_state->getPressure(
                                        &rho_B,
                                        m_ptr_B,
                                        &E_B);
                                    
                                    const double p_T = d_equation_of_state->getPressure(
                                        &rho_T,
                                        m_ptr_T,
                                        &E_T);
                                    
                                    const double c_B = d_equation_of_state->getSoundSpeedWithPressure(
                                        &rho_B,
                                        &p_B);
                                    
                                    const double c_T = d_equation_of_state->getSoundSpeedWithPressure(
                                        &rho_T,
                                        &p_T);
                                    
                                    const double sp_y_B = fabs(v_B) + c_B;
                                    const double sp_y_T = fabs(v_T) + c_T;
                                    
                                    std::vector<double> F_y_cell;
                                    std::vector<double> F_y_B;
                                    std::vector<double> F_y_T;
                                    
                                    F_y_cell.push_back(rho_cell*v_cell);
                                    F_y_cell.push_back(rho_cell*v_cell*u_cell);
                                    F_y_cell.push_back(rho_cell*v_cell*v_cell + p_cell);
                                    F_y_cell.push_back(rho_cell*v_cell*w_cell);
                                    F_y_cell.push_back(v_cell*(E_cell + p_cell));
                                    
                                    F_y_B.push_back(rho_B*v_B);
                                    F_y_B.push_back(rho_B*v_B*u_B);
                                    F_y_B.push_back(rho_B*v_B*v_B + p_B);
                                    F_y_B.push_back(rho_B*v_B*w_B);
                                    F_y_B.push_back(v_B*(E_B + p_B));
                                    
                                    F_y_T.push_back(rho_T*v_T);
                                    F_y_T.push_back(rho_T*v_T*u_T);
                                    F_y_T.push_back(rho_T*v_T*v_T + p_T);
                                    F_y_T.push_back(rho_T*v_T*w_T);
                                    F_y_T.push_back(v_T*(E_T + p_T));
                                    
                                    const double alpha_y_B = fmax(sp_y_B, sp_y_cell);
                                    const double alpha_y_T = fmax(sp_y_T, sp_y_cell);
                                    
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        double* F_y = convective_flux_intermediate->getPointer(1, ei);
                                        F_y[idx_face_y - 1] = 0.5*dt*(F_y_B[ei] + F_y_cell[ei] -
                                            alpha_y_B*(Q[ei][idx_cell] - Q[ei][idx_cell_B]));
                                        
                                        F_y[idx_face_y]     = 0.5*dt*(F_y_cell[ei] + F_y_T[ei] -
                                            alpha_y_T*(Q[ei][idx_cell_T] - Q[ei][idx_cell]));
                                    }
                                }
                                
                                // Correct the fluxes in the z direction.
                                {
                                    const int idx_cell_B   = (i + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k - 1 + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_cell_F   = (i + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + 1 + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const double& rho_B = Q[0][idx_cell_B];
                                    const double& rho_F = Q[0][idx_cell_F];
                                    
                                    const double& E_B = Q[1 + d_dim.getValue()][idx_cell_B];
                                    const double& E_F = Q[1 + d_dim.getValue()][idx_cell_F];
                                    
                                    const double u_B = Q[1][idx_cell_B]/rho_B;
                                    const double u_F = Q[1][idx_cell_F]/rho_F;
                                    
                                    const double v_B = Q[2][idx_cell_B]/rho_B;
                                    const double v_F = Q[2][idx_cell_F]/rho_F;
                                    
                                    const double w_B = Q[3][idx_cell_B]/rho_B;
                                    const double w_F = Q[3][idx_cell_F]/rho_F;
                                    
                                    std::vector<const double*> m_ptr_B;
                                    for (int di = 0; di < d_dim.getValue(); di++)
                                    {
                                        m_ptr_B.push_back(&Q[1 + di][idx_cell_B]);
                                    }
                                    
                                    std::vector<const double*> m_ptr_F;
                                    for (int di = 0; di < d_dim.getValue(); di++)
                                    {
                                        m_ptr_F.push_back(&Q[1 + di][idx_cell_F]);
                                    }
                                    
                                    const double p_B = d_equation_of_state->getPressure(
                                        &rho_B,
                                        m_ptr_B,
                                        &E_B);
                                    
                                    const double p_F = d_equation_of_state->getPressure(
                                        &rho_F,
                                        m_ptr_F,
                                        &E_F);
                                    
                                    const double c_B = d_equation_of_state->getSoundSpeedWithPressure(
                                        &rho_B,
                                        &p_B);
                                    
                                    const double c_F = d_equation_of_state->getSoundSpeedWithPressure(
                                        &rho_F,
                                        &p_F);
                                    
                                    const double sp_z_B = fabs(w_B) + c_B;
                                    const double sp_z_F = fabs(w_F) + c_F;
                                    
                                    std::vector<double> F_z_cell;
                                    std::vector<double> F_z_B;
                                    std::vector<double> F_z_F;
                                    
                                    F_z_cell.push_back(rho_cell*w_cell);
                                    F_z_cell.push_back(rho_cell*w_cell*u_cell);
                                    F_z_cell.push_back(rho_cell*w_cell*v_cell);
                                    F_z_cell.push_back(rho_cell*w_cell*w_cell + p_cell);
                                    F_z_cell.push_back(w_cell*(E_cell + p_cell));
                                    
                                    F_z_B.push_back(rho_B*w_B);
                                    F_z_B.push_back(rho_B*w_B*u_B);
                                    F_z_B.push_back(rho_B*w_B*v_B);
                                    F_z_B.push_back(rho_B*w_B*w_B + p_B);
                                    F_z_B.push_back(w_B*(E_B + p_B));
                                    
                                    F_z_F.push_back(rho_F*w_F);
                                    F_z_F.push_back(rho_F*w_F*u_F);
                                    F_z_F.push_back(rho_F*w_F*v_F);
                                    F_z_F.push_back(rho_F*w_F*w_F + p_F);
                                    F_z_F.push_back(w_F*(E_F + p_F));
                                    
                                    const double alpha_z_B = fmax(sp_z_B, sp_z_cell);
                                    const double alpha_z_F = fmax(sp_z_F, sp_z_cell);
                                    
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        double* F_z = convective_flux_intermediate->getPointer(2, ei);
                                        F_z[idx_face_z - 1] = 0.5*dt*(F_z_B[ei] + F_z_cell[ei] -
                                            alpha_z_B*(Q[ei][idx_cell] - Q[ei][idx_cell_B]));
                                        
                                        F_z[idx_face_z]     = 0.5*dt*(F_z_cell[ei] + F_z_F[ei] -
                                            alpha_z_F*(Q[ei][idx_cell_F] - Q[ei][idx_cell]));
                                    }
                                }
                                
                                break;
                            }
                            case FOUR_EQN_SHYUE:
                            {
                                const double& rho_cell = Q[0][idx_cell];
                                
                                const double& E_cell = Q[1 + d_dim.getValue()][idx_cell];
                                
                                const double u_cell = Q[1][idx_cell]/rho_cell;
                                
                                const double v_cell = Q[2][idx_cell]/rho_cell;
                                
                                const double w_cell = Q[3][idx_cell]/rho_cell;
                                
                                std::vector<const double*> m_ptr_cell;
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    m_ptr_cell.push_back(&Q[1 + di][idx_cell]);
                                }
                                
                                std::vector<const double*> Y_ptr_cell;
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_ptr_cell.push_back(&Q[2 + d_dim.getValue() + si][idx_cell]);
                                }
                                
                                const double p_cell = d_equation_of_state->getPressureWithMassFraction(
                                    &rho_cell,
                                    m_ptr_cell,
                                    &E_cell,
                                    Y_ptr_cell);
                                
                                const double c_cell = d_equation_of_state->getSoundSpeedWithMassFractionAndPressure(
                                    &rho_cell,
                                    Y_ptr_cell,
                                    &p_cell);
                                
                                const double sp_x_cell = fabs(u_cell) + c_cell;
                                const double sp_y_cell = fabs(v_cell) + c_cell;
                                const double sp_z_cell = fabs(w_cell) + c_cell;
                                
                                // Correct the fluxes in the x direction.
                                {
                                    // Compute indices of left and right cells.
                                    const int idx_cell_L   = (i - 1 + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_cell_R   = (i + 1 + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const double& rho_L = Q[0][idx_cell_L];
                                    const double& rho_R = Q[0][idx_cell_R];
                                    
                                    const double& E_L = Q[1 + d_dim.getValue()][idx_cell_L];
                                    const double& E_R = Q[1 + d_dim.getValue()][idx_cell_R];
                                    
                                    const double u_L = Q[1][idx_cell_L]/rho_L;
                                    const double u_R = Q[1][idx_cell_R]/rho_R;
                                    
                                    const double v_L = Q[2][idx_cell_L]/rho_L;
                                    const double v_R = Q[2][idx_cell_R]/rho_R;
                                    
                                    const double w_L = Q[3][idx_cell_L]/rho_L;
                                    const double w_R = Q[3][idx_cell_R]/rho_R;
                                    
                                    std::vector<const double*> m_ptr_L;
                                    for (int di = 0; di < d_dim.getValue(); di++)
                                    {
                                        m_ptr_L.push_back(&Q[1 + di][idx_cell_L]);
                                    }
                                    
                                    std::vector<const double*> m_ptr_R;
                                    for (int di = 0; di < d_dim.getValue(); di++)
                                    {
                                        m_ptr_R.push_back(&Q[1 + di][idx_cell_R]);
                                    }
                                    
                                    std::vector<const double*> Y_ptr_L;
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Y_ptr_L.push_back(&Q[2 + d_dim.getValue() + si][idx_cell_L]);
                                    }
                                    
                                    std::vector<const double*> Y_ptr_R;
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Y_ptr_R.push_back(&Q[2 + d_dim.getValue() + si][idx_cell_R]);
                                    }
                                    
                                    const double p_L = d_equation_of_state->getPressureWithMassFraction(
                                        &rho_L,
                                        m_ptr_L,
                                        &E_L,
                                        Y_ptr_L);
                                    
                                    const double p_R = d_equation_of_state->getPressureWithMassFraction(
                                        &rho_R,
                                        m_ptr_R,
                                        &E_R,
                                        Y_ptr_R);
                                    
                                    const double c_L = d_equation_of_state->getSoundSpeedWithMassFractionAndPressure(
                                        &rho_L,
                                        Y_ptr_L,
                                        &p_L);
                                    
                                    const double c_R = d_equation_of_state->getSoundSpeedWithMassFractionAndPressure(
                                        &rho_R,
                                        Y_ptr_R,
                                        &p_R);
                                    
                                    const double sp_x_L = fabs(u_L) + c_L;
                                    const double sp_x_R = fabs(u_R) + c_R;
                                    
                                    std::vector<double> F_x_cell;
                                    std::vector<double> F_x_L;
                                    std::vector<double> F_x_R;
                                    
                                    F_x_cell.push_back(rho_cell*u_cell);
                                    F_x_cell.push_back(rho_cell*u_cell*u_cell + p_cell);
                                    F_x_cell.push_back(rho_cell*u_cell*v_cell);
                                    F_x_cell.push_back(rho_cell*u_cell*w_cell);
                                    F_x_cell.push_back(u_cell*(E_cell + p_cell));
                                    for (int si = 0; si < d_num_species - 1; si++)
                                    {
                                        F_x_cell.push_back(u_cell*(*Y_ptr_cell[si]));
                                    }
                                    
                                    F_x_L.push_back(rho_L*u_L);
                                    F_x_L.push_back(rho_L*u_L*u_L + p_L);
                                    F_x_L.push_back(rho_L*u_L*v_L);
                                    F_x_L.push_back(rho_L*u_L*w_L);
                                    F_x_L.push_back(u_L*(E_L + p_L));
                                    for (int si = 0; si < d_num_species - 1; si++)
                                    {
                                        F_x_L.push_back(u_L*(*Y_ptr_L[si]));
                                    }
                                    
                                    F_x_R.push_back(rho_R*u_R);
                                    F_x_R.push_back(rho_R*u_R*u_R + p_R);
                                    F_x_R.push_back(rho_R*u_R*v_R);
                                    F_x_R.push_back(rho_R*u_R*w_R);
                                    F_x_R.push_back(u_R*(E_R + p_R));
                                    for (int si = 0; si < d_num_species - 1; si++)
                                    {
                                        F_x_R.push_back(u_R*(*Y_ptr_R[si]));
                                    }
                                    
                                    const double alpha_x_L = fmax(sp_x_L, sp_x_cell);
                                    const double alpha_x_R = fmax(sp_x_R, sp_x_cell);
                                    
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        double* F_x = convective_flux_intermediate->getPointer(0, ei);
                                        F_x[idx_face_x - 1] = 0.5*dt*(F_x_L[ei] + F_x_cell[ei] -
                                            alpha_x_L*(Q[ei][idx_cell] - Q[ei][idx_cell_L]));
                                        
                                        F_x[idx_face_x]     = 0.5*dt*(F_x_cell[ei] + F_x_R[ei] -
                                            alpha_x_R*(Q[ei][idx_cell_R] - Q[ei][idx_cell]));
                                    }
                                }
                                
                                // Correct the fluxes in the y direction.
                                {
                                    const int idx_cell_B   = (i + d_num_ghosts[0]) +
                                        (j - 1 + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_cell_T   = (i + d_num_ghosts[0]) +
                                        (j + 1 + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const double& rho_B = Q[0][idx_cell_B];
                                    const double& rho_T = Q[0][idx_cell_T];
                                    
                                    const double& E_B = Q[1 + d_dim.getValue()][idx_cell_B];
                                    const double& E_T = Q[1 + d_dim.getValue()][idx_cell_T];
                                    
                                    const double u_B = Q[1][idx_cell_B]/rho_B;
                                    const double u_T = Q[1][idx_cell_T]/rho_T;
                                    
                                    const double v_B = Q[2][idx_cell_B]/rho_B;
                                    const double v_T = Q[2][idx_cell_T]/rho_T;
                                    
                                    const double w_B = Q[3][idx_cell_B]/rho_B;
                                    const double w_T = Q[3][idx_cell_T]/rho_T;
                                    
                                    std::vector<const double*> m_ptr_B;
                                    for (int di = 0; di < d_dim.getValue(); di++)
                                    {
                                        m_ptr_B.push_back(&Q[1 + di][idx_cell_B]);
                                    }
                                    
                                    std::vector<const double*> m_ptr_T;
                                    for (int di = 0; di < d_dim.getValue(); di++)
                                    {
                                        m_ptr_T.push_back(&Q[1 + di][idx_cell_T]);
                                    }
                                    
                                    std::vector<const double*> Y_ptr_B;
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Y_ptr_B.push_back(&Q[2 + d_dim.getValue() + si][idx_cell_B]);
                                    }
                                    
                                    std::vector<const double*> Y_ptr_T;
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Y_ptr_T.push_back(&Q[2 + d_dim.getValue() + si][idx_cell_T]);
                                    }
                                    
                                    const double p_B = d_equation_of_state->getPressureWithMassFraction(
                                        &rho_B,
                                        m_ptr_B,
                                        &E_B,
                                        Y_ptr_B);
                                    
                                    const double p_T = d_equation_of_state->getPressureWithMassFraction(
                                        &rho_T,
                                        m_ptr_T,
                                        &E_T,
                                        Y_ptr_T);
                                    
                                    const double c_B = d_equation_of_state->getSoundSpeedWithMassFractionAndPressure(
                                        &rho_B,
                                        Y_ptr_B,
                                        &p_B);
                                    
                                    const double c_T = d_equation_of_state->getSoundSpeedWithMassFractionAndPressure(
                                        &rho_T,
                                        Y_ptr_T,
                                        &p_T);
                                    
                                    const double sp_y_B = fabs(v_B) + c_B;
                                    const double sp_y_T = fabs(v_T) + c_T;
                                    
                                    std::vector<double> F_y_cell;
                                    std::vector<double> F_y_B;
                                    std::vector<double> F_y_T;
                                    
                                    F_y_cell.push_back(rho_cell*v_cell);
                                    F_y_cell.push_back(rho_cell*v_cell*u_cell);
                                    F_y_cell.push_back(rho_cell*v_cell*v_cell + p_cell);
                                    F_y_cell.push_back(rho_cell*v_cell*w_cell);
                                    F_y_cell.push_back(v_cell*(E_cell + p_cell));
                                    for (int si = 0; si < d_num_species - 1; si++)
                                    {
                                        F_y_cell.push_back(v_cell*(*Y_ptr_cell[si]));
                                    }
                                    
                                    F_y_B.push_back(rho_B*v_B);
                                    F_y_B.push_back(rho_B*v_B*u_B);
                                    F_y_B.push_back(rho_B*v_B*v_B + p_B);
                                    F_y_B.push_back(rho_B*v_B*w_B);
                                    F_y_B.push_back(v_B*(E_B + p_B));
                                    for (int si = 0; si < d_num_species - 1; si++)
                                    {
                                        F_y_B.push_back(v_B*(*Y_ptr_B[si]));
                                    }
                                    
                                    F_y_T.push_back(rho_T*v_T);
                                    F_y_T.push_back(rho_T*v_T*u_T);
                                    F_y_T.push_back(rho_T*v_T*v_T + p_T);
                                    F_y_T.push_back(rho_T*v_T*w_T);
                                    F_y_T.push_back(v_T*(E_T + p_T));
                                    for (int si = 0; si < d_num_species - 1; si++)
                                    {
                                        F_y_T.push_back(v_T*(*Y_ptr_T[si]));
                                    }
                                    
                                    const double alpha_y_B = fmax(sp_y_B, sp_y_cell);
                                    const double alpha_y_T = fmax(sp_y_T, sp_y_cell);
                                    
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        double* F_y = convective_flux_intermediate->getPointer(1, ei);
                                        F_y[idx_face_y - 1] = 0.5*dt*(F_y_B[ei] + F_y_cell[ei] -
                                            alpha_y_B*(Q[ei][idx_cell] - Q[ei][idx_cell_B]));
                                        
                                        F_y[idx_face_y]     = 0.5*dt*(F_y_cell[ei] + F_y_T[ei] -
                                            alpha_y_T*(Q[ei][idx_cell_T] - Q[ei][idx_cell]));
                                    }
                                }
                                
                                // Correct the fluxes in the z direction.
                                {
                                    const int idx_cell_B   = (i + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k - 1 + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_cell_F   = (i + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + 1 + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const double& rho_B = Q[0][idx_cell_B];
                                    const double& rho_F = Q[0][idx_cell_F];
                                    
                                    const double& E_B = Q[1 + d_dim.getValue()][idx_cell_B];
                                    const double& E_F = Q[1 + d_dim.getValue()][idx_cell_F];
                                    
                                    const double u_B = Q[1][idx_cell_B]/rho_B;
                                    const double u_F = Q[1][idx_cell_F]/rho_F;
                                    
                                    const double v_B = Q[2][idx_cell_B]/rho_B;
                                    const double v_F = Q[2][idx_cell_F]/rho_F;
                                    
                                    const double w_B = Q[3][idx_cell_B]/rho_B;
                                    const double w_F = Q[3][idx_cell_F]/rho_F;
                                    
                                    std::vector<const double*> m_ptr_B;
                                    for (int di = 0; di < d_dim.getValue(); di++)
                                    {
                                        m_ptr_B.push_back(&Q[1 + di][idx_cell_B]);
                                    }
                                    
                                    std::vector<const double*> m_ptr_F;
                                    for (int di = 0; di < d_dim.getValue(); di++)
                                    {
                                        m_ptr_F.push_back(&Q[1 + di][idx_cell_F]);
                                    }
                                    
                                    std::vector<const double*> Y_ptr_B;
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Y_ptr_B.push_back(&Q[2 + d_dim.getValue() + si][idx_cell_B]);
                                    }
                                    
                                    std::vector<const double*> Y_ptr_F;
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Y_ptr_F.push_back(&Q[2 + d_dim.getValue() + si][idx_cell_F]);
                                    }
                                    
                                    const double p_B = d_equation_of_state->getPressureWithMassFraction(
                                        &rho_B,
                                        m_ptr_B,
                                        &E_B,
                                        Y_ptr_B);
                                    
                                    const double p_F = d_equation_of_state->getPressureWithMassFraction(
                                        &rho_F,
                                        m_ptr_F,
                                        &E_F,
                                        Y_ptr_F);
                                    
                                    const double c_B = d_equation_of_state->getSoundSpeedWithMassFractionAndPressure(
                                        &rho_B,
                                        Y_ptr_B,
                                        &p_B);
                                    
                                    const double c_F = d_equation_of_state->getSoundSpeedWithMassFractionAndPressure(
                                        &rho_F,
                                        Y_ptr_F,
                                        &p_F);
                                    
                                    const double sp_z_B = fabs(w_B) + c_B;
                                    const double sp_z_F = fabs(w_F) + c_F;
                                    
                                    std::vector<double> F_z_cell;
                                    std::vector<double> F_z_B;
                                    std::vector<double> F_z_F;
                                    
                                    F_z_cell.push_back(rho_cell*w_cell);
                                    F_z_cell.push_back(rho_cell*w_cell*u_cell);
                                    F_z_cell.push_back(rho_cell*w_cell*v_cell);
                                    F_z_cell.push_back(rho_cell*w_cell*w_cell + p_cell);
                                    F_z_cell.push_back(w_cell*(E_cell + p_cell));
                                    for (int si = 0; si < d_num_species - 1; si++)
                                    {
                                        F_z_cell.push_back(w_cell*(*Y_ptr_cell[si]));
                                    }
                                    
                                    F_z_B.push_back(rho_B*w_B);
                                    F_z_B.push_back(rho_B*w_B*u_B);
                                    F_z_B.push_back(rho_B*w_B*v_B);
                                    F_z_B.push_back(rho_B*w_B*w_B + p_B);
                                    F_z_B.push_back(w_B*(E_B + p_B));
                                    for (int si = 0; si < d_num_species - 1; si++)
                                    {
                                        F_z_B.push_back(w_B*(*Y_ptr_B[si]));
                                    }
                                    
                                    F_z_F.push_back(rho_F*w_F);
                                    F_z_F.push_back(rho_F*w_F*u_F);
                                    F_z_F.push_back(rho_F*w_F*v_F);
                                    F_z_F.push_back(rho_F*w_F*w_F + p_F);
                                    F_z_F.push_back(w_F*(E_F + p_F));
                                    for (int si = 0; si < d_num_species - 1; si++)
                                    {
                                        F_z_F.push_back(w_F*(*Y_ptr_F[si]));
                                    }
                                    
                                    const double alpha_z_B = fmax(sp_z_B, sp_z_cell);
                                    const double alpha_z_F = fmax(sp_z_F, sp_z_cell);
                                    
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        double* F_z = convective_flux_intermediate->getPointer(2, ei);
                                        F_z[idx_face_z - 1] = 0.5*dt*(F_z_B[ei] + F_z_cell[ei] -
                                            alpha_z_B*(Q[ei][idx_cell] - Q[ei][idx_cell_B]));
                                        
                                        F_z[idx_face_z]     = 0.5*dt*(F_z_cell[ei] + F_z_F[ei] -
                                            alpha_z_F*(Q[ei][idx_cell_F] - Q[ei][idx_cell]));
                                    }
                                }
                                
                                break;
                            }
                            case FIVE_EQN_ALLAIRE:
                            {
                                std::vector<const double*> Z_rho_ptr_cell;
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Z_rho_ptr_cell.push_back(&Q[si][idx_cell]);
                                }
                                
                                const double rho_cell = d_equation_of_state->getTotalDensity(Z_rho_ptr_cell);
                                
                                const double& E_cell = Q[d_num_species + d_dim.getValue()][idx_cell];
                                
                                const double u_cell = Q[d_num_species][idx_cell]/rho_cell;
                                
                                const double v_cell = Q[d_num_species + 1][idx_cell]/rho_cell;
                                
                                const double w_cell = Q[d_num_species + 2][idx_cell]/rho_cell;
                                
                                std::vector<const double*> m_ptr_cell;
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    m_ptr_cell.push_back(&Q[d_num_species + di][idx_cell]);
                                }
                                
                                std::vector<const double*> Z_ptr_cell;
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Z_ptr_cell.push_back(&Q[1 + d_num_species + d_dim.getValue() + si][idx_cell]);
                                }
                                
                                const double p_cell = d_equation_of_state->getPressureWithVolumeFraction(
                                    &rho_cell,
                                    m_ptr_cell,
                                    &E_cell,
                                    Z_ptr_cell);
                                
                                const double c_cell = d_equation_of_state->getSoundSpeedWithVolumeFractionAndPressure(
                                    &rho_cell,
                                    Z_ptr_cell,
                                    &p_cell);
                                
                                const double sp_x_cell = fabs(u_cell) + c_cell;
                                const double sp_y_cell = fabs(v_cell) + c_cell;
                                const double sp_z_cell = fabs(w_cell) + c_cell;
                                
                                // Correct the fluxes in the x direction.
                                {
                                    // Compute indices of left and right cells.
                                    const int idx_cell_L   = (i - 1 + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_cell_R   = (i + 1 + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    std::vector<const double*> Z_rho_ptr_L;
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Z_rho_ptr_L.push_back(&Q[si][idx_cell_L]);
                                    }
                                    
                                    std::vector<const double*> Z_rho_ptr_R;
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Z_rho_ptr_R.push_back(&Q[si][idx_cell_R]);
                                    }
                                    
                                    const double& rho_L = d_equation_of_state->getTotalDensity(Z_rho_ptr_L);
                                    const double& rho_R = d_equation_of_state->getTotalDensity(Z_rho_ptr_R);
                                    
                                    const double& E_L = Q[d_num_species + d_dim.getValue()][idx_cell_L];
                                    const double& E_R = Q[d_num_species + d_dim.getValue()][idx_cell_R];
                                    
                                    const double u_L = Q[d_num_species][idx_cell_L]/rho_L;
                                    const double u_R = Q[d_num_species][idx_cell_R]/rho_R;
                                    
                                    const double v_L = Q[d_num_species + 1][idx_cell_L]/rho_L;
                                    const double v_R = Q[d_num_species + 1][idx_cell_R]/rho_R;
                                    
                                    const double w_L = Q[d_num_species + 2][idx_cell_L]/rho_L;
                                    const double w_R = Q[d_num_species + 2][idx_cell_R]/rho_R;
                                    
                                    std::vector<const double*> m_ptr_L;
                                    for (int di = 0; di < d_dim.getValue(); di++)
                                    {
                                        m_ptr_L.push_back(&Q[d_num_species + di][idx_cell_L]);
                                    }
                                    
                                    std::vector<const double*> m_ptr_R;
                                    for (int di = 0; di < d_dim.getValue(); di++)
                                    {
                                        m_ptr_R.push_back(&Q[d_num_species + di][idx_cell_R]);
                                    }
                                    
                                    std::vector<const double*> Z_ptr_L;
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Z_ptr_L.push_back(&Q[1 + d_num_species + d_dim.getValue() + si][idx_cell_L]);
                                    }
                                    
                                    std::vector<const double*> Z_ptr_R;
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Z_ptr_R.push_back(&Q[1 + d_num_species + d_dim.getValue() + si][idx_cell_R]);
                                    }
                                    
                                    const double p_L = d_equation_of_state->getPressureWithVolumeFraction(
                                        &rho_L,
                                        m_ptr_L,
                                        &E_L,
                                        Z_ptr_L);
                                    
                                    const double p_R = d_equation_of_state->getPressureWithVolumeFraction(
                                        &rho_R,
                                        m_ptr_R,
                                        &E_R,
                                        Z_ptr_R);
                                    
                                    const double c_L = d_equation_of_state->getSoundSpeedWithVolumeFractionAndPressure(
                                        &rho_L,
                                        Z_ptr_L,
                                        &p_L);
                                    
                                    const double c_R = d_equation_of_state->getSoundSpeedWithVolumeFractionAndPressure(
                                        &rho_R,
                                        Z_ptr_R,
                                        &p_R);
                                    
                                    const double sp_x_L = fabs(u_L) + c_L;
                                    const double sp_x_R = fabs(u_R) + c_R;
                                    
                                    std::vector<double> F_x_cell;
                                    std::vector<double> F_x_L;
                                    std::vector<double> F_x_R;
                                    
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        F_x_cell.push_back(u_cell*(*Z_rho_ptr_cell[si]));
                                    }
                                    F_x_cell.push_back(rho_cell*u_cell*u_cell + p_cell);
                                    F_x_cell.push_back(rho_cell*u_cell*v_cell);
                                    F_x_cell.push_back(rho_cell*u_cell*w_cell);
                                    F_x_cell.push_back(u_cell*(E_cell + p_cell));
                                    for (int si = 0; si < d_num_species - 1; si++)
                                    {
                                        F_x_cell.push_back(u_cell*(*Z_ptr_cell[si]));
                                    }
                                    
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        F_x_L.push_back(u_L*(*Z_rho_ptr_L[si]));
                                    }
                                    F_x_L.push_back(rho_L*u_L*u_L + p_L);
                                    F_x_L.push_back(rho_L*u_L*v_L);
                                    F_x_L.push_back(rho_L*u_L*w_L);
                                    F_x_L.push_back(u_L*(E_L + p_L));
                                    for (int si = 0; si < d_num_species - 1; si++)
                                    {
                                        F_x_L.push_back(u_L*(*Z_ptr_L[si]));
                                    }
                                    
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        F_x_R.push_back(u_R*(*Z_rho_ptr_R[si]));
                                    }
                                    F_x_R.push_back(rho_R*u_R*u_R + p_R);
                                    F_x_R.push_back(rho_R*u_R*v_R);
                                    F_x_R.push_back(rho_R*u_R*w_R);
                                    F_x_R.push_back(u_R*(E_R + p_R));
                                    for (int si = 0; si < d_num_species - 1; si++)
                                    {
                                        F_x_R.push_back(u_R*(*Z_ptr_R[si]));
                                    }
                                    
                                    const double alpha_x_L = fmax(sp_x_L, sp_x_cell);
                                    const double alpha_x_R = fmax(sp_x_R, sp_x_cell);
                                    
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        double* F_x = convective_flux_intermediate->getPointer(0, ei);
                                        F_x[idx_face_x - 1] = 0.5*dt*(F_x_L[ei] + F_x_cell[ei] -
                                            alpha_x_L*(Q[ei][idx_cell] - Q[ei][idx_cell_L]));
                                        
                                        F_x[idx_face_x]     = 0.5*dt*(F_x_cell[ei] + F_x_R[ei] -
                                            alpha_x_R*(Q[ei][idx_cell_R] - Q[ei][idx_cell]));
                                    }
                                }
                                
                                // Correct the fluxes in the y direction.
                                {
                                    const int idx_cell_B   = (i + d_num_ghosts[0]) +
                                        (j - 1 + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_cell_T   = (i + d_num_ghosts[0]) +
                                        (j + 1 + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    std::vector<const double*> Z_rho_ptr_B;
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Z_rho_ptr_B.push_back(&Q[si][idx_cell_B]);
                                    }
                                    
                                    std::vector<const double*> Z_rho_ptr_T;
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Z_rho_ptr_T.push_back(&Q[si][idx_cell_T]);
                                    }
                                    
                                    const double& rho_B = d_equation_of_state->getTotalDensity(Z_rho_ptr_B);
                                    const double& rho_T = d_equation_of_state->getTotalDensity(Z_rho_ptr_T);
                                    
                                    const double& E_B = Q[d_num_species + d_dim.getValue()][idx_cell_B];
                                    const double& E_T = Q[d_num_species + d_dim.getValue()][idx_cell_T];
                                    
                                    const double u_B = Q[d_num_species][idx_cell_B]/rho_B;
                                    const double u_T = Q[d_num_species][idx_cell_T]/rho_T;
                                    
                                    const double v_B = Q[d_num_species + 1][idx_cell_B]/rho_B;
                                    const double v_T = Q[d_num_species + 1][idx_cell_T]/rho_T;
                                    
                                    const double w_B = Q[d_num_species + 2][idx_cell_B]/rho_B;
                                    const double w_T = Q[d_num_species + 2][idx_cell_T]/rho_T;
                                    
                                    std::vector<const double*> m_ptr_B;
                                    for (int di = 0; di < d_dim.getValue(); di++)
                                    {
                                        m_ptr_B.push_back(&Q[d_num_species + di][idx_cell_B]);
                                    }
                                    
                                    std::vector<const double*> m_ptr_T;
                                    for (int di = 0; di < d_dim.getValue(); di++)
                                    {
                                        m_ptr_T.push_back(&Q[d_num_species + di][idx_cell_T]);
                                    }
                                    
                                    std::vector<const double*> Z_ptr_B;
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Z_ptr_B.push_back(&Q[1 + d_num_species + d_dim.getValue() + si][idx_cell_B]);
                                    }
                                    
                                    std::vector<const double*> Z_ptr_T;
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Z_ptr_T.push_back(&Q[1 + d_num_species + d_dim.getValue() + si][idx_cell_T]);
                                    }
                                    
                                    const double p_B = d_equation_of_state->getPressureWithVolumeFraction(
                                        &rho_B,
                                        m_ptr_B,
                                        &E_B,
                                        Z_ptr_B);
                                    
                                    const double p_T = d_equation_of_state->getPressureWithVolumeFraction(
                                        &rho_T,
                                        m_ptr_T,
                                        &E_T,
                                        Z_ptr_T);
                                    
                                    const double c_B = d_equation_of_state->getSoundSpeedWithVolumeFractionAndPressure(
                                        &rho_B,
                                        Z_ptr_B,
                                        &p_B);
                                    
                                    const double c_T = d_equation_of_state->getSoundSpeedWithVolumeFractionAndPressure(
                                        &rho_T,
                                        Z_ptr_T,
                                        &p_T);
                                    
                                    const double sp_y_B = fabs(v_B) + c_B;
                                    const double sp_y_T = fabs(v_T) + c_T;
                                    
                                    std::vector<double> F_y_cell;
                                    std::vector<double> F_y_B;
                                    std::vector<double> F_y_T;
                                    
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        F_y_cell.push_back(v_cell*(*Z_rho_ptr_cell[si]));
                                    }
                                    F_y_cell.push_back(rho_cell*v_cell*u_cell);
                                    F_y_cell.push_back(rho_cell*v_cell*v_cell + p_cell);
                                    F_y_cell.push_back(rho_cell*v_cell*w_cell);
                                    F_y_cell.push_back(v_cell*(E_cell + p_cell));
                                    for (int si = 0; si < d_num_species - 1; si++)
                                    {
                                        F_y_cell.push_back(v_cell*(*Z_ptr_cell[si]));
                                    }
                                    
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        F_y_B.push_back(v_B*(*Z_rho_ptr_B[si]));
                                    }
                                    F_y_B.push_back(rho_B*v_B*u_B);
                                    F_y_B.push_back(rho_B*v_B*v_B + p_B);
                                    F_y_B.push_back(rho_B*v_B*w_B);
                                    F_y_B.push_back(v_B*(E_B + p_B));
                                    for (int si = 0; si < d_num_species - 1; si++)
                                    {
                                        F_y_B.push_back(v_B*(*Z_ptr_B[si]));
                                    }
                                    
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        F_y_T.push_back(v_T*(*Z_rho_ptr_T[si]));
                                    }
                                    F_y_T.push_back(rho_T*v_T*u_T);
                                    F_y_T.push_back(rho_T*v_T*v_T + p_T);
                                    F_y_T.push_back(rho_T*v_T*w_T);
                                    F_y_T.push_back(v_T*(E_T + p_T));
                                    for (int si = 0; si < d_num_species - 1; si++)
                                    {
                                        F_y_T.push_back(v_T*(*Z_ptr_T[si]));
                                    }
                                    
                                    const double alpha_y_B = fmax(sp_y_B, sp_y_cell);
                                    const double alpha_y_T = fmax(sp_y_T, sp_y_cell);
                                    
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        double* F_y = convective_flux_intermediate->getPointer(1, ei);
                                        F_y[idx_face_y - 1] = 0.5*dt*(F_y_B[ei] + F_y_cell[ei] -
                                            alpha_y_B*(Q[ei][idx_cell] - Q[ei][idx_cell_B]));
                                        
                                        F_y[idx_face_y]     = 0.5*dt*(F_y_cell[ei] + F_y_T[ei] -
                                            alpha_y_T*(Q[ei][idx_cell_T] - Q[ei][idx_cell]));
                                    }
                                }
                                
                                // Correct the fluxes in the z direction.
                                {
                                    const int idx_cell_B   = (i + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k - 1 + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_cell_F   = (i + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + 1 + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    std::vector<const double*> Z_rho_ptr_B;
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Z_rho_ptr_B.push_back(&Q[si][idx_cell_B]);
                                    }
                                    
                                    std::vector<const double*> Z_rho_ptr_F;
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Z_rho_ptr_F.push_back(&Q[si][idx_cell_F]);
                                    }
                                    
                                    const double& rho_B = d_equation_of_state->getTotalDensity(Z_rho_ptr_B);
                                    const double& rho_F = d_equation_of_state->getTotalDensity(Z_rho_ptr_F);
                                    
                                    const double& E_B = Q[d_num_species + d_dim.getValue()][idx_cell_B];
                                    const double& E_F = Q[d_num_species + d_dim.getValue()][idx_cell_F];
                                    
                                    const double u_B = Q[d_num_species][idx_cell_B]/rho_B;
                                    const double u_F = Q[d_num_species][idx_cell_F]/rho_F;
                                    
                                    const double v_B = Q[d_num_species + 1][idx_cell_B]/rho_B;
                                    const double v_F = Q[d_num_species + 1][idx_cell_F]/rho_F;
                                    
                                    const double w_B = Q[d_num_species + 2][idx_cell_B]/rho_B;
                                    const double w_F = Q[d_num_species + 2][idx_cell_F]/rho_F;
                                    
                                    std::vector<const double*> m_ptr_B;
                                    for (int di = 0; di < d_dim.getValue(); di++)
                                    {
                                        m_ptr_B.push_back(&Q[d_num_species + di][idx_cell_B]);
                                    }
                                    
                                    std::vector<const double*> m_ptr_F;
                                    for (int di = 0; di < d_dim.getValue(); di++)
                                    {
                                        m_ptr_F.push_back(&Q[d_num_species + di][idx_cell_F]);
                                    }
                                    
                                    std::vector<const double*> Z_ptr_B;
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Z_ptr_B.push_back(&Q[1 + d_num_species + d_dim.getValue() + si][idx_cell_B]);
                                    }
                                    
                                    std::vector<const double*> Z_ptr_F;
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Z_ptr_F.push_back(&Q[1 + d_num_species + d_dim.getValue() + si][idx_cell_F]);
                                    }
                                    
                                    const double p_B = d_equation_of_state->getPressureWithVolumeFraction(
                                        &rho_B,
                                        m_ptr_B,
                                        &E_B,
                                        Z_ptr_B);
                                    
                                    const double p_F = d_equation_of_state->getPressureWithVolumeFraction(
                                        &rho_F,
                                        m_ptr_F,
                                        &E_F,
                                        Z_ptr_F);
                                    
                                    const double c_B = d_equation_of_state->getSoundSpeedWithVolumeFractionAndPressure(
                                        &rho_B,
                                        Z_ptr_B,
                                        &p_B);
                                    
                                    const double c_F = d_equation_of_state->getSoundSpeedWithVolumeFractionAndPressure(
                                        &rho_F,
                                        Z_ptr_F,
                                        &p_F);
                                    
                                    const double sp_z_B = fabs(w_B) + c_B;
                                    const double sp_z_F = fabs(w_F) + c_F;
                                    
                                    std::vector<double> F_z_cell;
                                    std::vector<double> F_z_B;
                                    std::vector<double> F_z_F;
                                    
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        F_z_cell.push_back(w_cell*(*Z_rho_ptr_cell[si]));
                                    }
                                    F_z_cell.push_back(rho_cell*w_cell*u_cell);
                                    F_z_cell.push_back(rho_cell*w_cell*v_cell);
                                    F_z_cell.push_back(rho_cell*w_cell*w_cell + p_cell);
                                    F_z_cell.push_back(w_cell*(E_cell + p_cell));
                                    for (int si = 0; si < d_num_species - 1; si++)
                                    {
                                        F_z_cell.push_back(w_cell*(*Z_ptr_cell[si]));
                                    }
                                    
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        F_z_B.push_back(w_B*(*Z_rho_ptr_B[si]));
                                    }
                                    F_z_B.push_back(rho_B*w_B*u_B);
                                    F_z_B.push_back(rho_B*w_B*v_B);
                                    F_z_B.push_back(rho_B*w_B*w_B + p_B);
                                    F_z_B.push_back(w_B*(E_B + p_B));
                                    for (int si = 0; si < d_num_species - 1; si++)
                                    {
                                        F_z_B.push_back(w_B*(*Z_ptr_B[si]));
                                    }
                                    
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        F_z_F.push_back(w_F*(*Z_rho_ptr_F[si]));
                                    }
                                    F_z_F.push_back(rho_F*w_F*u_F);
                                    F_z_F.push_back(rho_F*w_F*v_F);
                                    F_z_F.push_back(rho_F*w_F*w_F + p_F);
                                    F_z_F.push_back(w_F*(E_F + p_F));
                                    for (int si = 0; si < d_num_species - 1; si++)
                                    {
                                        F_z_F.push_back(w_F*(*Z_ptr_F[si]));
                                    }
                                    
                                    const double alpha_z_B = fmax(sp_z_B, sp_z_cell);
                                    const double alpha_z_F = fmax(sp_z_F, sp_z_cell);
                                    
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        double* F_z = convective_flux_intermediate->getPointer(2, ei);
                                        F_z[idx_face_z - 1] = 0.5*dt*(F_z_B[ei] + F_z_cell[ei] -
                                            alpha_z_B*(Q[ei][idx_cell] - Q[ei][idx_cell_B]));
                                        
                                        F_z[idx_face_z]     = 0.5*dt*(F_z_cell[ei] + F_z_F[ei] -
                                            alpha_z_F*(Q[ei][idx_cell_F] - Q[ei][idx_cell]));
                                    }
                                }
                                
                                break;
                            }
                            default:
                            TBOX_ERROR(d_object_name
                                       << ": "
                                       << "Unknown d_flow_model."
                                       << std::endl);
                        }
                    }
                }
            }
        }
    }
}
