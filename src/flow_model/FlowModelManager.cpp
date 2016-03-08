#include "flow_model/FlowModelManager.hpp"

FlowModelManager::FlowModelManager(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
            const std::string& flow_model_str,
            int& num_species):
                d_object_name(object_name),
                d_dim(dim),
                d_grid_geometry(grid_geometry),
                d_num_ghosts(hier::IntVector::getZero(d_dim)),
                d_num_species(num_species),
                d_equation_of_state_initialized(false)
{
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
    
    if (d_num_species > 1 && d_flow_model == SINGLE_SPECIES)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Number of species = "
            << d_num_species
            << " shouldn't use single-species model."
            << std::endl); 
    }
    
    /*
     * Initialize the conservative variables based on the flow model.
     */
    switch (d_flow_model)
    {
        case SINGLE_SPECIES:
        {
            d_density = boost::shared_ptr<pdat::CellVariable<double> > (
                new pdat::CellVariable<double>(d_dim, "density", 1));
            
            break;
        }
        case FOUR_EQN_SHYUE:
        {
            d_density = boost::shared_ptr<pdat::CellVariable<double> > (
                new pdat::CellVariable<double>(d_dim, "density", 1));
            
            break;
        }
        case FIVE_EQN_ALLAIRE:
        {
            d_partial_density = boost::shared_ptr<pdat::CellVariable<double> > (
                new pdat::CellVariable<double>(d_dim, "partial density", d_num_species));
            
            break;
        }
        default:
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "d_flow_model '"
                << d_flow_model
                << "' not yet implemented."
                << std::endl);
        }
    }
    
    d_momentum = boost::shared_ptr<pdat::CellVariable<double> > (
        new pdat::CellVariable<double>(d_dim, "momentum", d_dim.getValue()));
    
    d_total_energy = boost::shared_ptr<pdat::CellVariable<double> > (
        new pdat::CellVariable<double>(d_dim, "total energy", 1));
    
    switch (d_flow_model)
    {
        case SINGLE_SPECIES:
        {
            break;
        }
        case FOUR_EQN_SHYUE:
        {
            d_mass_fraction = boost::shared_ptr<pdat::CellVariable<double> > (
                new pdat::CellVariable<double>(d_dim, "mass fraction", d_num_species));
            
            break;
        }
        case FIVE_EQN_ALLAIRE:
        {
            d_volume_fraction = boost::shared_ptr<pdat::CellVariable<double> > (
                new pdat::CellVariable<double>(d_dim, "volume fraction", d_num_species));
            
            break;
        }
        default:
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "d_flow_model '"
                << d_flow_model
                << "' not yet implemented."
                << std::endl);
        }
    }
}


/*
 * Initialize d_equation_of_state.
 */
void
FlowModelManager::initializeEquationOfState(
    const std::string& equation_of_state_string,
    const boost::shared_ptr<tbox::Database>& equation_of_state_db,
    boost::shared_ptr<EquationOfState>& equation_of_state)
{
    if (equation_of_state_string == "IDEAL_GAS")
    {
        switch (d_flow_model)
        {
            case SINGLE_SPECIES:
            {
                equation_of_state.reset(new EquationOfStateIdealGas(
                    "ideal gas",
                    d_dim,
                    d_num_species,
                    equation_of_state_db,
                    NO_ASSUMPTION));
                
                break;
            }
            case FOUR_EQN_SHYUE:
            {
                equation_of_state.reset(new EquationOfStateIdealGas(
                    "ideal gas",
                    d_dim,
                    d_num_species,
                    equation_of_state_db,
                    ISOTHERMAL));
                
                break;
            }
            case FIVE_EQN_ALLAIRE:
            {
                equation_of_state.reset(new EquationOfStateIdealGas(
                    "ideal gas",
                    d_dim,
                    d_num_species,
                    equation_of_state_db,
                    ISOBARIC));
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "d_flow_model '"
                    << d_flow_model
                    << "' not yet implemented."
                    << std::endl);
            }
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Unknown equation_of_state string = '"
            << equation_of_state_string
            << "' found in data for Equation_of_state."
            << std::endl);  
    }
    
    d_equation_of_state = equation_of_state;
    d_equation_of_state_initialized = true;
}


/*
 * Initialize d_conv_flux_reconstructor.
 * The function also computes the number of ghost cells required.
 */
void
FlowModelManager::initializeConvectiveFluxReconstructor(
    const std::string& shock_capturing_scheme_str,
    const boost::shared_ptr<tbox::Database>& shock_capturing_scheme_db,
    boost::shared_ptr<ConvectiveFluxReconstructor>& conv_flux_reconstructor,
    boost::shared_ptr<pdat::FaceVariable<double> >& convective_flux,
    boost::shared_ptr<pdat::CellVariable<double> >& source)
{
    if (shock_capturing_scheme_str == "LLF")
    {
        conv_flux_reconstructor.reset(new ConvectiveFluxReconstructorLLF(
            "LLF",
            d_dim,
            d_grid_geometry,
            d_flow_model,
            d_num_eqn,
            d_num_species,
            d_equation_of_state,
            shock_capturing_scheme_db));
    }
    else if (shock_capturing_scheme_str == "FIRST_ORDER_HLLC")
    {
        conv_flux_reconstructor.reset(new ConvectiveFluxReconstructorFirstOrderHLLC(
            "first order HLLC",
            d_dim,
            d_grid_geometry,
            d_flow_model,
            d_num_eqn,
            d_num_species,
            d_equation_of_state,
            shock_capturing_scheme_db));
    }
    else if (shock_capturing_scheme_str == "WENO_JS5_LLF")
    {
        conv_flux_reconstructor.reset(new ConvectiveFluxReconstructorWENO_JS5_LLF(
            "WENO-JS5-LLF",
            d_dim,
            d_grid_geometry,
            d_flow_model,
            d_num_eqn,
            d_num_species,
            d_equation_of_state,
            shock_capturing_scheme_db));
    }
    else if (shock_capturing_scheme_str == "WENO_CU6_M2_LLF")
    {
        conv_flux_reconstructor.reset(new ConvectiveFluxReconstructorWENO_CU6_M2_LLF(
            "WENO-CU6-M2-LLF",
            d_dim,
            d_grid_geometry,
            d_flow_model,
            d_num_eqn,
            d_num_species,
            d_equation_of_state,
            shock_capturing_scheme_db));
    }
    else if (shock_capturing_scheme_str == "WCNS_JS5_HLLC_HLL")
    {
        conv_flux_reconstructor.reset(new ConvectiveFluxReconstructorWCNS_JS5_HLLC_HLL(
            "WCNS-JS5-HLLC-HLL",
            d_dim,
            d_grid_geometry,
            d_flow_model,
            d_num_eqn,
            d_num_species,
            d_equation_of_state,
            shock_capturing_scheme_db));
    }
    else if (shock_capturing_scheme_str == "WCNS_CU6_M2_HLLC_HLL")
    {
        conv_flux_reconstructor.reset(new ConvectiveFluxReconstructorWCNS_CU6_M2_HLLC_HLL(
            "WCNS-CU6-M2-HLLC-HLL",
            d_dim,
            d_grid_geometry,
            d_flow_model,
            d_num_eqn,
            d_num_species,
            d_equation_of_state,
            shock_capturing_scheme_db));
    }
    else if (shock_capturing_scheme_str == "WCNS_HW6_HLLC_HLL")
    {
        conv_flux_reconstructor.reset(new ConvectiveFluxReconstructorWCNS_HW6_HLLC_HLL(
            "WCNS-HW6-HLLC-HLL",
            d_dim,
            d_grid_geometry,
            d_flow_model,
            d_num_eqn,
            d_num_species,
            d_equation_of_state,
            shock_capturing_scheme_db));
    }
    else if (shock_capturing_scheme_str == "WCNS_HW6_LD_HLLC_HLL")
    {
        conv_flux_reconstructor.reset(new ConvectiveFluxReconstructorWCNS_HW6_LD_HLLC_HLL(
            "WCNS-HW6-LD-HLLC-HLL",
            d_dim,
            d_grid_geometry,
            d_flow_model,
            d_num_eqn,
            d_num_species,
            d_equation_of_state,
            shock_capturing_scheme_db));
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
    
    d_conv_flux_reconstructor = conv_flux_reconstructor;
    
    /*
     * Initialize the number of ghost cells and boost::shared_ptr of the variables
     * in d_conv_flux_reconstructor.
     */
    setVariablesForConvectiveFluxReconstructor(
        convective_flux,
        source);
}


/*
 * Initialize d_initial_conditions.
 */
void
FlowModelManager::initializeInitialConditions(
    const std::string& project_name,
    boost::shared_ptr<InitialConditions>& initial_conditions)
{
    if (!d_equation_of_state_initialized)
    {
        TBOX_ERROR(d_object_name
            << ": initializeInitialConditions()\n"
            << "d_equation_of_state is not initialized yet."
            << std::endl);
    }
    
    initial_conditions.reset(new InitialConditions(
        "initial conditions",
        project_name,
        d_dim,
        d_grid_geometry,
        d_flow_model,
        d_num_species,
        d_equation_of_state));
    
    d_initial_conditions = initial_conditions;
    
    /*
     * Initialize boost::shared_ptr of the variables
     * in d_initial_conditions.
     */
    setVariablesForInitialConditions();
}


/*
 * Initialize d_Euler_boundary_conditions.
 */
void
FlowModelManager::initializeEulerBoundaryConditions(
    const std::string& project_name,
    const boost::shared_ptr<tbox::Database>& Euler_boundary_conditions_db,
    const bool& Euler_boundary_conditions_db_is_from_restart,
    boost::shared_ptr<EulerBoundaryConditions>& Euler_boundary_conditions)
{
    if (!d_equation_of_state_initialized)
    {
        TBOX_ERROR(d_object_name
            << ": initializeEulerBoundaryConditions()\n"
            << "d_equation_of_state is not initialized yet."
            << std::endl);
    }
    
    Euler_boundary_conditions.reset(new EulerBoundaryConditions(
        "Euler boundary conditions",
        project_name,
        d_dim,
        d_grid_geometry,
        d_flow_model,
        d_num_species,
        d_equation_of_state,
        Euler_boundary_conditions_db,
        Euler_boundary_conditions_db_is_from_restart));
    
    d_Euler_boundary_conditions = Euler_boundary_conditions;
    
    /*
     * Initialize boost::shared_ptr of the variables
     * in d_Euler_boundary_conditions.
     */
    setVariablesForEulerBoundaryConditions();
}


/*
 * Initialize d_gradient_tagger.
 */
void
FlowModelManager::initializeGradientTagger(
    const boost::shared_ptr<tbox::Database>& gradient_tagger_db,
    boost::shared_ptr<GradientTagger>& gradient_tagger)
{
    if (!d_equation_of_state_initialized)
    {
        TBOX_ERROR(d_object_name
            << ": initializeGradientTagger()\n"
            << "d_equation_of_state is not initialized yet."
            << std::endl);
    }
    
    gradient_tagger.reset(new GradientTagger(
        "gradient tagger",
        d_dim,
        d_grid_geometry,
        d_flow_model,
        d_num_species,
        d_equation_of_state,
        gradient_tagger_db));
    
    d_gradient_tagger = gradient_tagger;
    
    /*
     * Initialize the number of ghost cells and boost::shared_ptr of the variables
     * in d_gradient_tagger.
     */
    setVariablesForGradientTagger();
}


/*
 * Initialize d_multiresolution_tagger.
 */
void
FlowModelManager::initializeMultiresolutionTagger(
   const boost::shared_ptr<tbox::Database>& multiresolution_tagger_db,
   boost::shared_ptr<MultiresolutionTagger>& multiresolution_tagger)
{
    if (!d_equation_of_state_initialized)
    {
        TBOX_ERROR(d_object_name
            << ": initializeMultiresolutionTagger()\n"
            << "d_equation_of_state is not initialized yet."
            << std::endl);
    }
    
    if (multiresolution_tagger_db != nullptr)
    {
        multiresolution_tagger.reset(new MultiresolutionTagger(
            "multiresolution tagger",
            d_dim,
            d_grid_geometry,
            d_flow_model,
            d_num_species,
            d_equation_of_state,
            multiresolution_tagger_db));
        
        d_multiresolution_tagger = multiresolution_tagger;
        
        /*
         * Initialize the number of ghost cells and boost::shared_ptr of the variables
         * in d_multiresolution_tagger.
         */
        setVariablesForMultiresolutionTagger();
    }
}


/*
 * Compute and set the number of ghost cells needed.
 */
void
FlowModelManager::setNumberOfGhostCells()
{
    /*
     * Compute the number of ghost cells required.
     */
    d_num_ghosts = d_conv_flux_reconstructor->
        getConvectiveFluxNumberOfGhostCells();
    
    if (d_gradient_tagger != nullptr)
    {
        d_num_ghosts = hier::IntVector::max(
            d_num_ghosts,
            d_gradient_tagger->getGradientTaggerNumberOfGhostCells());
    }
    
    if (d_multiresolution_tagger != nullptr)
    {
        d_num_ghosts = hier::IntVector::max(
            d_num_ghosts,
            d_multiresolution_tagger->getMultiresolutionTaggerNumberOfGhostCells());
    }
    
    /*
     * Set the number of ghost cells required.
     */
    d_conv_flux_reconstructor->setNumberOfGhostCells(d_num_ghosts);
    
    d_Euler_boundary_conditions->setNumberOfGhostCells(d_num_ghosts);
    
    if (d_gradient_tagger != nullptr)
    {
        d_gradient_tagger->setNumberOfGhostCells(d_num_ghosts);
    }
    
    if (d_multiresolution_tagger != nullptr)
    {
        d_multiresolution_tagger->setNumberOfGhostCells(d_num_ghosts);
    }
}


/*
 * Register the conservative variables.
 */
void
FlowModelManager::registerConservativeVariables(
    RungeKuttaLevelIntegrator* integrator)
{
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
                << "d_flow_model '"
                << d_flow_model
                << "' not yet implemented."
                << std::endl);
        }
    }
}


/*
 * Register the temporary variables used in refinement taggers.
 */
void
FlowModelManager::registerRefinementTaggerVariables(
   RungeKuttaLevelIntegrator* integrator)
{
    if (d_gradient_tagger != nullptr)
    {
        d_gradient_tagger->registerGradientTaggerVariables(integrator);
    }
    
    if (d_multiresolution_tagger != nullptr)
    {
        d_multiresolution_tagger->registerMultiresolutionTaggerVariables(integrator);
    }
}


/*
 * Set the plotting context.
 */
void
FlowModelManager::setPlotContext(
    const boost::shared_ptr<hier::VariableContext>& plot_context)
{
    d_plot_context = plot_context;
}


/*
 * Register the plotting quantities.
 */
#ifdef HAVE_HDF5
void
FlowModelManager::registerPlotQuantities(
    RungeKuttaLevelIntegrator* integrator,
    const boost::shared_ptr<appu::VisItDataWriter>& visit_writer)
{
    hier::VariableDatabase* vardb = hier::VariableDatabase::getDatabase();
    
    switch (d_flow_model)
    {
        case SINGLE_SPECIES:
        {
            visit_writer->registerPlotQuantity(
                "density",
                "SCALAR",
                vardb->mapVariableAndContextToIndex(
                   d_density,
                   d_plot_context));
            
            visit_writer->registerPlotQuantity(
                "momentum",
                "VECTOR",
                vardb->mapVariableAndContextToIndex(
                   d_momentum,
                   d_plot_context));
            
            visit_writer->registerPlotQuantity(
                "total energy",
                "SCALAR",
                vardb->mapVariableAndContextToIndex(
                   d_total_energy,
                   d_plot_context));
            
            visit_writer->registerDerivedPlotQuantity("pressure",
                "SCALAR",
                this);
            
            visit_writer->registerDerivedPlotQuantity("sound speed",
                "SCALAR",
                this);
            
            visit_writer->registerDerivedPlotQuantity("velocity",
                "VECTOR",
                this);
            
            break;
        }
        case FOUR_EQN_SHYUE:
        {
            visit_writer->registerPlotQuantity(
                "density",
                "SCALAR",
                vardb->mapVariableAndContextToIndex(
                   d_density,
                   d_plot_context));
            
            visit_writer->registerPlotQuantity(
                "momentum",
                "VECTOR",
                vardb->mapVariableAndContextToIndex(
                   d_momentum,
                   d_plot_context));
            
            visit_writer->registerPlotQuantity(
                "total energy",
                "SCALAR",
                vardb->mapVariableAndContextToIndex(
                   d_total_energy,
                   d_plot_context));
            
            for (int si = 0; si < d_num_species; si++)
            {
                std::string mass_fraction_name =
                    "mass fraction " + tbox::Utilities::intToString(si);
                    
                visit_writer->registerPlotQuantity(
                    mass_fraction_name,
                    "SCALAR",
                    vardb->mapVariableAndContextToIndex(
                       d_mass_fraction,
                       d_plot_context),
                    si);
            }
            
            visit_writer->registerDerivedPlotQuantity("pressure",
                "SCALAR",
                this);
            
            visit_writer->registerDerivedPlotQuantity("sound speed",
                "SCALAR",
                this);
            
            visit_writer->registerDerivedPlotQuantity("velocity",
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
                
                visit_writer->registerPlotQuantity(
                    partial_density_name,
                    "SCALAR",
                    vardb->mapVariableAndContextToIndex(
                        d_partial_density,
                        d_plot_context),
                    si);
            }
            
            visit_writer->registerPlotQuantity(
                "momentum",
                "VECTOR",
                vardb->mapVariableAndContextToIndex(
                   d_momentum,
                   d_plot_context));
            
            visit_writer->registerPlotQuantity(
                "total energy",
                "SCALAR",
                vardb->mapVariableAndContextToIndex(
                   d_total_energy,
                   d_plot_context));
            
            for (int si = 0; si < d_num_species; si++)
            {
                std::string volume_fraction_name =
                    "volume fraction " + tbox::Utilities::intToString(si);
                    
                visit_writer->registerPlotQuantity(
                    volume_fraction_name,
                    "SCALAR",
                    vardb->mapVariableAndContextToIndex(
                       d_volume_fraction,
                       d_plot_context),
                    si);
            }
            
            visit_writer->registerDerivedPlotQuantity("pressure",
                "SCALAR",
                this);
            
            visit_writer->registerDerivedPlotQuantity("sound speed",
                "SCALAR",
                this);
            
            visit_writer->registerDerivedPlotQuantity("velocity",
                "VECTOR",
                this);
            
            visit_writer->registerDerivedPlotQuantity("density",
                "SCALAR",
                this);
            
            for (int si = 0; si < d_num_species; si++)
            {
                std::string mass_fraction_name =
                    "mass fraction " + tbox::Utilities::intToString(si);
                    
                visit_writer->registerDerivedPlotQuantity(mass_fraction_name,
                                                            "SCALAR",
                                                            this);
            }
            
            break;
        }
        default:
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "d_flow_model '"
                << d_flow_model
                << "' not yet implemented."
                << std::endl);
        }
    }
    
    if (d_multiresolution_tagger != nullptr)
    {
        d_multiresolution_tagger->registerPlotQuantities(
            integrator,
            visit_writer,
            d_plot_context);
    }
    
    if (d_gradient_tagger != nullptr)
    {
        d_gradient_tagger->registerPlotQuantities(
            integrator,
            visit_writer,
            d_plot_context);
    }
}
#endif


/*
 * Compute the stable time increment for patch using a CFL
 * condition and return the computed dt.
 */
double
FlowModelManager::computeStableDtOnPatch(
    hier::Patch& patch,
    const bool initial_time,
    const double dt_time,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
    NULL_USE(initial_time);
    NULL_USE(dt_time);
    
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
                    patch.getPatchData(d_density, data_context)));
            
            boost::shared_ptr<pdat::CellData<double> > momentum(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_momentum, data_context)));
            
            boost::shared_ptr<pdat::CellData<double> > total_energy(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_total_energy, data_context)));
    
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
                    patch.getPatchData(d_density, data_context)));
            
            boost::shared_ptr<pdat::CellData<double> > momentum(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_momentum, data_context)));
            
            boost::shared_ptr<pdat::CellData<double> > total_energy(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_total_energy, data_context)));
            
            boost::shared_ptr<pdat::CellData<double> > mass_fraction(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_mass_fraction, data_context)));
            
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
                    patch.getPatchData(d_partial_density, data_context)));
            
            boost::shared_ptr<pdat::CellData<double> > momentum(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_momentum, data_context)));
            
            boost::shared_ptr<pdat::CellData<double> > total_energy(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_total_energy, data_context)));
            
            boost::shared_ptr<pdat::CellData<double> > volume_fraction(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_volume_fraction, data_context)));
            
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
                << "d_flow_model '"
                << d_flow_model
                << "' not yet implemented."
                << std::endl);
        }
    }
    
    stable_dt = 1.0/stable_spectral_radius;
    
    return stable_dt;
}


/*
 * Create a vector of pointers to conservative variables with the
 * given data context.
 */
std::vector<double*>
FlowModelManager::createConservativeVariableVector(
    hier::Patch& patch,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const bool fill_zero)
{
    std::vector<double*> Q;
    
    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(patch_geom);
#endif
    
    // Get the box that covers the interior of patch.
    const hier::Box interior_box = patch.getBox();
    
    switch (d_flow_model)
    {
        case SINGLE_SPECIES:
        {
            boost::shared_ptr<pdat::CellData<double> > density(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_density, data_context)));
            
            boost::shared_ptr<pdat::CellData<double> > momentum(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_momentum, data_context)));
            
            boost::shared_ptr<pdat::CellData<double> > total_energy(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_total_energy, data_context)));
            
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(density);
            TBOX_ASSERT(momentum);
            TBOX_ASSERT(total_energy);
            
            TBOX_ASSERT(density->getGhostCellWidth() == d_num_ghosts);
            TBOX_ASSERT(momentum->getGhostCellWidth() == d_num_ghosts);
            TBOX_ASSERT(total_energy->getGhostCellWidth() == d_num_ghosts);
#endif
            
            if (fill_zero)
            {
                // Initialize all time-dependent data within the interior box with zero values
                density->fillAll(0.0, interior_box);
                momentum->fillAll(0.0, interior_box);
                total_energy->fillAll(0.0, interior_box);
            }
            
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
                    patch.getPatchData(d_density, data_context)));
            
            boost::shared_ptr<pdat::CellData<double> > momentum(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_momentum, data_context)));
            
            boost::shared_ptr<pdat::CellData<double> > total_energy(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_total_energy, data_context)));
            
            boost::shared_ptr<pdat::CellData<double> > mass_fraction(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_mass_fraction, data_context)));
            
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
            
            if (fill_zero)
            {
                // Initialize all time-dependent data within the interior box with zero values
                density->fillAll(0.0, interior_box);
                momentum->fillAll(0.0, interior_box);
                total_energy->fillAll(0.0, interior_box);
                mass_fraction->fillAll(0.0, interior_box);
            }
            
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
                    patch.getPatchData(d_partial_density, data_context)));
            
            boost::shared_ptr<pdat::CellData<double> > momentum(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_momentum, data_context)));
            
            boost::shared_ptr<pdat::CellData<double> > total_energy(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_total_energy, data_context)));
            
            boost::shared_ptr<pdat::CellData<double> > volume_fraction(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch.getPatchData(d_volume_fraction, data_context)));
                        
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
            if (fill_zero)
            {
                // Initialize all time-dependent data within the interior box with zero values
                partial_density->fillAll(0.0, interior_box);
                momentum->fillAll(0.0, interior_box);
                total_energy->fillAll(0.0, interior_box);
                volume_fraction->fillAll(0.0, interior_box);
            }
            
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
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "d_flow_model '"
                << d_flow_model
                << "' not yet implemented."
                << std::endl);
        }
    }
    
    return Q;
}


/*
 * Update the mass fraction/volume fraction of the last species
 * in the conservative variable vector.
 */
void
FlowModelManager::updateConservativeVariableVector(
    hier::Patch& patch,
    std::vector<double*>& conservative_variable_vector)
{
    std::vector<double*> Q = conservative_variable_vector;
    
    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(patch_geom);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box dummy_box = patch.getBox();
    const hier::Box interior_box = dummy_box;
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // ghost cells.
    dummy_box.grow(d_num_ghosts);
    const hier::Box ghost_box = dummy_box;
    const hier::IntVector ghostcell_dims = ghost_box.numberCells();
    
    if (d_dim == tbox::Dimension(1))
    {
        switch (d_flow_model)
        {
            case SINGLE_SPECIES:
            {
                break;
            }
            case FOUR_EQN_SHYUE:
            {
                for (int i = 0; i < interior_dims[0]; i++)
                {
                    // Compute index of time-dependent data.
                    int idx_cell   = i + d_num_ghosts[0];
                    
                    Q[d_num_eqn][idx_cell] = 1.0;
                    
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q[d_num_eqn][idx_cell] -= Q[d_num_eqn - 1 - si][idx_cell];
                    }
                }
                
                break;
            }
            case FIVE_EQN_ALLAIRE:
            {
                for (int i = 0; i < interior_dims[0]; i++)
                {
                    // Compute index of time-dependent data.
                    int idx_cell   = i + d_num_ghosts[0];
                    
                    Q[d_num_eqn][idx_cell] = 1.0;
                    
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q[d_num_eqn][idx_cell] -= Q[d_num_eqn - 1 - si][idx_cell];
                    }
                }
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "d_flow_model '"
                    << d_flow_model
                    << "' not yet implemented."
                    << std::endl);
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        switch (d_flow_model)
        {
            case SINGLE_SPECIES:
            {
                break;
            }
            case FOUR_EQN_SHYUE:
            {
                for (int j = 0; j < interior_dims[1]; j++)
                {
                    for (int i = 0; i < interior_dims[0]; i++)
                    {
                        // Compute index of time-dependent data.
                        int idx_cell   = (i + d_num_ghosts[0]) +
                            (j + d_num_ghosts[1])*ghostcell_dims[0];
                        
                        Q[d_num_eqn][idx_cell] = 1.0;
                        
                        for (int si = 0; si < d_num_species - 1; si++)
                        {
                            Q[d_num_eqn][idx_cell] -= Q[d_num_eqn - 1 - si][idx_cell];
                        }
                    }
                }
                
                break;
            }
            case FIVE_EQN_ALLAIRE:
            {
                for (int j = 0; j < interior_dims[1]; j++)
                {
                    for (int i = 0; i < interior_dims[0]; i++)
                    {
                        // Compute index of time-dependent data.
                        int idx_cell   = (i + d_num_ghosts[0]) +
                            (j + d_num_ghosts[1])*ghostcell_dims[0];
                        
                        Q[d_num_eqn][idx_cell] = 1.0;
                        
                        for (int si = 0; si < d_num_species - 1; si++)
                        {
                            Q[d_num_eqn][idx_cell] -= Q[d_num_eqn - 1 - si][idx_cell];
                        }
                    }
                }
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "d_flow_model '"
                    << d_flow_model
                    << "' not yet implemented."
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
                break;
            }
            case FOUR_EQN_SHYUE:
            {
                for (int k = 0; k < interior_dims[2]; k++)
                {
                    for (int j = 0; j < interior_dims[1]; j++)
                    {
                        for (int i = 0; i < interior_dims[0]; i++)
                        {
                            // Compute index of time-dependent data.
                            int idx_cell   = (i + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            Q[d_num_eqn][idx_cell] = 1.0;
                            
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                Q[d_num_eqn][idx_cell] -= Q[d_num_eqn - 1 - si][idx_cell];
                            }
                        }
                    }
                }
                
                break;
            }
            case FIVE_EQN_ALLAIRE:
            {
                for (int k = 0; k < interior_dims[2]; k++)
                {
                    for (int j = 0; j < interior_dims[1]; j++)
                    {
                        for (int i = 0; i < interior_dims[0]; i++)
                        {
                            // Compute index of time-dependent data.
                            int idx_cell   = (i + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                            
                            Q[d_num_eqn][idx_cell] = 1.0;
                            
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                Q[d_num_eqn][idx_cell] -= Q[d_num_eqn - 1 - si][idx_cell];
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
                    << "d_flow_model '"
                    << d_flow_model
                    << "' not yet implemented."
                    << std::endl);
            }
        }
    }
}


/*
 * Put the d_flow_model into the restart database.
 */
void
FlowModelManager::putToRestart(
    const boost::shared_ptr<tbox::Database>& restart_db) const
{
    switch(d_flow_model)
    {
        case SINGLE_SPECIES:
        {
            restart_db->putString("d_flow_model", "SINGLE_SPECIES");
            break;
        }
        case FOUR_EQN_SHYUE:
        {
            restart_db->putString("d_flow_model", "FOUR_EQN_SHYUE");
            break;
        }
        case FIVE_EQN_ALLAIRE:
        {
            restart_db->putString("d_flow_model", "FIVE_EQN_ALLAIRE");
            break;
        }
        default:
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "d_flow_model '"
                << d_flow_model
                << "' not yet implemented."
                << std::endl);
        }
    }
}


/*
 * Compute derived plot quantities registered with the VisIt data
 * writers from data that is maintained on each patch in the
 * hierarchy.
 */
bool
FlowModelManager::packDerivedDataIntoDoubleBuffer(
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
                    << "d_flow_model '"
                    << d_flow_model
                    << "' not yet implemented."
                    << std::endl);
            }
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "variable_name '"
            << variable_name
            << "' is unknown."
            << std::endl);
    }
    
    return data_on_patch;
}


/*
 * Print all characteristics of d_flow_model.
 */
void
FlowModelManager::printClassData(std::ostream& os) const
{
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
                << "d_flow_model '"
                << d_flow_model
                << "' not yet implemented."
                << std::endl);
        }
    }
    
    os << "d_num_eqn = " << d_num_eqn << std::endl;
}


/*
 * Print all data statistics for the time-dependent variables.
 */
void
FlowModelManager::printDataStatistics(
    std::ostream& os,
    const boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy) const
{
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    math::HierarchyCellDataOpsReal<double> cell_double_operator(patch_hierarchy, 0, 0);
    
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
                << "d_flow_model '"
                << d_flow_model
                << "' not yet implemented."
                << std::endl);
        }
    }
}


/*
 * Initialize the number of ghost cells and boost::shared_ptr of the variables
 * in d_conv_flux_reconstructor.
 */
void
FlowModelManager::setVariablesForConvectiveFluxReconstructor(
    boost::shared_ptr<pdat::FaceVariable<double> >& convective_flux,
    boost::shared_ptr<pdat::CellVariable<double> >& source)
{
    switch (d_flow_model)
    {
        case SINGLE_SPECIES:
        {
            d_conv_flux_reconstructor->setVariablesForSingleSpecies(
                d_density,
                d_momentum,
                d_total_energy,
                convective_flux,
                source);
            
            break;
        }
        case FOUR_EQN_SHYUE:
        {
            d_conv_flux_reconstructor->setVariablesForFourEqnShyue(
                d_density,
                d_momentum,
                d_total_energy,
                d_mass_fraction,
                convective_flux,
                source);
            
            break;
        }
        case FIVE_EQN_ALLAIRE:
        {
            d_conv_flux_reconstructor->setVariablesForFiveEqnAllaire(
                d_partial_density,
                d_momentum,
                d_total_energy,
                d_volume_fraction,
                convective_flux,
                source);
            
            break;
        }
        default:
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "d_flow_model '"
                << d_flow_model
                << "' not yet implemented."
                << std::endl);
        }
    }
}


/*
 * Initialize boost::shared_ptr of the variables
 * in d_initial_conditions.
 */
void
FlowModelManager::setVariablesForInitialConditions()
{
    switch (d_flow_model)
    {
        case SINGLE_SPECIES:
        {
            d_initial_conditions->setVariablesForSingleSpecies(
                d_density,
                d_momentum,
                d_total_energy);
            
            break;
        }
        case FOUR_EQN_SHYUE:
        {
            d_initial_conditions->setVariablesForFourEqnShyue(
                d_density,
                d_momentum,
                d_total_energy,
                d_mass_fraction);
            
            break;
        }
        case FIVE_EQN_ALLAIRE:
        {
            d_initial_conditions->setVariablesForFiveEqnAllaire(
                d_partial_density,
                d_momentum,
                d_total_energy,
                d_volume_fraction);
            
            break;
        }
        default:
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "d_flow_model '"
                << d_flow_model
                << "' not yet implemented."
                << std::endl);
        }
    }
}


/*
 * Initialize boost::shared_ptr of the variables
 * in d_Euler_boundary_conditions.
 */
void
FlowModelManager::setVariablesForEulerBoundaryConditions()
{
    switch (d_flow_model)
    {
        case SINGLE_SPECIES:
        {
            d_Euler_boundary_conditions->setVariablesForSingleSpecies(
                d_density,
                d_momentum,
                d_total_energy);
            
            break;
        }
        case FOUR_EQN_SHYUE:
        {
            d_Euler_boundary_conditions->setVariablesForFourEqnShyue(
                d_density,
                d_momentum,
                d_total_energy,
                d_mass_fraction);
            
            break;
        }
        case FIVE_EQN_ALLAIRE:
        {
            d_Euler_boundary_conditions->setVariablesForFiveEqnAllaire(
                d_partial_density,
                d_momentum,
                d_total_energy,
                d_volume_fraction);
            
            break;
        }
        default:
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "d_flow_model '"
                << d_flow_model
                << "' not yet implemented."
                << std::endl);
        }
    }
}


/*
 * Initialize the number of ghost cells and boost::shared_ptr of the variables
 * in d_gradient_tagger.
 */
void
FlowModelManager::setVariablesForGradientTagger()
{
    if (d_gradient_tagger != nullptr)
    {
        switch (d_flow_model)
        {
            case SINGLE_SPECIES:
            {
                d_gradient_tagger->setVariablesForSingleSpecies(
                    d_density,
                    d_momentum,
                    d_total_energy);
                
                break;
            }
            case FOUR_EQN_SHYUE:
            {
                d_gradient_tagger->setVariablesForFourEqnShyue(
                    d_density,
                    d_momentum,
                    d_total_energy,
                    d_mass_fraction);
                
                break;
            }
            case FIVE_EQN_ALLAIRE:
            {
                d_gradient_tagger->setVariablesForFiveEqnAllaire(
                    d_partial_density,
                    d_momentum,
                    d_total_energy,
                    d_volume_fraction);
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "d_flow_model '"
                    << d_flow_model
                    << "' not yet implemented."
                    << std::endl);
            }
        }
    }
}


/*
 * Initialize the number of ghost cells and boost::shared_ptr of the variables
 * in d_multiresolution_tagger.
 */
void
FlowModelManager::setVariablesForMultiresolutionTagger()
{
    if (d_multiresolution_tagger != nullptr)
    {
        switch (d_flow_model)
        {
            case SINGLE_SPECIES:
            {
                d_multiresolution_tagger->setVariablesForSingleSpecies(
                    d_density,
                    d_momentum,
                    d_total_energy);
                
                break;
            }
            case FOUR_EQN_SHYUE:
            {
                d_multiresolution_tagger->setVariablesForFourEqnShyue(
                    d_density,
                    d_momentum,
                    d_total_energy,
                    d_mass_fraction);
                
                break;
            }
            case FIVE_EQN_ALLAIRE:
            {
                d_multiresolution_tagger->setVariablesForFiveEqnAllaire(
                    d_partial_density,
                    d_momentum,
                    d_total_energy,
                    d_volume_fraction);
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "d_flow_model '"
                    << d_flow_model
                    << "' not yet implemented."
                    << std::endl);
            }
        }
    }
}
