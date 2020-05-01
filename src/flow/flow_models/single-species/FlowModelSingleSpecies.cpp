#include "flow/flow_models/single-species/FlowModelSingleSpecies.hpp"

#include "flow/flow_models/single-species/FlowModelBasicUtilitiesSingleSpecies.hpp"
#include "flow/flow_models/single-species/FlowModelBoundaryUtilitiesSingleSpecies.hpp"
#include "flow/flow_models/single-species/FlowModelDiffusiveFluxUtilitiesSingleSpecies.hpp"
#include "flow/flow_models/single-species/FlowModelRiemannSolverSingleSpecies.hpp"
#include "flow/flow_models/single-species/FlowModelStatisticsUtilitiesSingleSpecies.hpp"

boost::shared_ptr<pdat::CellVariable<double> > FlowModelSingleSpecies::s_variable_density;
boost::shared_ptr<pdat::CellVariable<double> > FlowModelSingleSpecies::s_variable_momentum;
boost::shared_ptr<pdat::CellVariable<double> > FlowModelSingleSpecies::s_variable_total_energy;

FlowModelSingleSpecies::FlowModelSingleSpecies(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_species,
    const boost::shared_ptr<tbox::Database>& flow_model_db):
        FlowModel(
            object_name,
            dim,
            grid_geometry,
            num_species,
            2 + dim.getValue(),
            flow_model_db),
        d_num_subghosts_velocity(-hier::IntVector::getOne(d_dim)),
        d_num_subghosts_internal_energy(-hier::IntVector::getOne(d_dim)),
        d_num_subghosts_pressure(-hier::IntVector::getOne(d_dim)),
        d_num_subghosts_sound_speed(-hier::IntVector::getOne(d_dim)),
        d_num_subghosts_temperature(-hier::IntVector::getOne(d_dim)),
        d_num_subghosts_convective_flux_x(-hier::IntVector::getOne(d_dim)),
        d_num_subghosts_convective_flux_y(-hier::IntVector::getOne(d_dim)),
        d_num_subghosts_convective_flux_z(-hier::IntVector::getOne(d_dim)),
        d_num_subghosts_max_wave_speed_x(-hier::IntVector::getOne(d_dim)),
        d_num_subghosts_max_wave_speed_y(-hier::IntVector::getOne(d_dim)),
        d_num_subghosts_max_wave_speed_z(-hier::IntVector::getOne(d_dim)),
        d_num_subghosts_max_diffusivity(-hier::IntVector::getOne(d_dim)),
        d_subghost_box_velocity(hier::Box::getEmptyBox(d_dim)),
        d_subghost_box_internal_energy(hier::Box::getEmptyBox(d_dim)),
        d_subghost_box_pressure(hier::Box::getEmptyBox(d_dim)),
        d_subghost_box_sound_speed(hier::Box::getEmptyBox(d_dim)),
        d_subghost_box_temperature(hier::Box::getEmptyBox(d_dim)),
        d_subghost_box_convective_flux_x(hier::Box::getEmptyBox(d_dim)),
        d_subghost_box_convective_flux_y(hier::Box::getEmptyBox(d_dim)),
        d_subghost_box_convective_flux_z(hier::Box::getEmptyBox(d_dim)),
        d_subghost_box_max_wave_speed_x(hier::Box::getEmptyBox(d_dim)),
        d_subghost_box_max_wave_speed_y(hier::Box::getEmptyBox(d_dim)),
        d_subghost_box_max_wave_speed_z(hier::Box::getEmptyBox(d_dim)),
        d_subghost_box_max_diffusivity(hier::Box::getEmptyBox(d_dim)),
        d_subghostcell_dims_velocity(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_internal_energy(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_pressure(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_sound_speed(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_temperature(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_convective_flux_x(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_convective_flux_y(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_convective_flux_z(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_max_wave_speed_x(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_max_wave_speed_y(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_max_wave_speed_z(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_max_diffusivity(hier::IntVector::getZero(d_dim)),
        d_cell_data_velocity_computed(false),
        d_cell_data_internal_energy_computed(false),
        d_cell_data_pressure_computed(false),
        d_cell_data_sound_speed_computed(false),
        d_cell_data_temperature_computed(false),
        d_cell_data_convective_flux_x_computed(false),
        d_cell_data_convective_flux_y_computed(false),
        d_cell_data_convective_flux_z_computed(false),
        d_cell_data_max_wave_speed_x_computed(false),
        d_cell_data_max_wave_speed_y_computed(false),
        d_cell_data_max_wave_speed_z_computed(false),
        d_cell_data_max_diffusivity_computed(false)
{
    d_eqn_form.reserve(d_num_eqn);
    
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        d_eqn_form.push_back(EQN_FORM::CONSERVATIVE);
    }
    
    /*
     * Initialize the conservative variables.
     */
    
    s_variable_density = boost::shared_ptr<pdat::CellVariable<double> > (
        new pdat::CellVariable<double>(d_dim, "density", 1));
    
    s_variable_momentum = boost::shared_ptr<pdat::CellVariable<double> > (
        new pdat::CellVariable<double>(d_dim, "momentum", d_dim.getValue()));
    
    s_variable_total_energy = boost::shared_ptr<pdat::CellVariable<double> > (
        new pdat::CellVariable<double>(d_dim, "total energy", 1));
    
    /*
     * Initialize d_equation_of_state_mixing_rules_manager and get the equation of state
     * mixing rules object.
     */
    
    if (flow_model_db->keyExists("equation_of_state"))
    {
        d_equation_of_state_str = flow_model_db->getString("equation_of_state");
    }
    else if (flow_model_db->keyExists("d_equation_of_state_str"))
    {
        d_equation_of_state_str = flow_model_db->getString("d_equation_of_state_str");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'equation_of_state'/'d_equation_of_state_str' found in data for flow model"
            << std::endl);
    }
    
    boost::shared_ptr<tbox::Database> equation_of_state_mixing_rules_db;
    
    if (flow_model_db->keyExists("Equation_of_state_mixing_rules"))
    {
        equation_of_state_mixing_rules_db = flow_model_db->getDatabase("Equation_of_state_mixing_rules");
    }
    else if (flow_model_db->keyExists("d_equation_of_state_mixing_rules_db"))
    {
        equation_of_state_mixing_rules_db = flow_model_db->getDatabase("d_equation_of_state_mixing_rules_db");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'Equation_of_state_mixing_rules'/'d_equation_of_state_mixing_rules_db' found in data for flow model"
            << std::endl);
    }
    
    d_equation_of_state_mixing_rules_manager.reset(new EquationOfStateMixingRulesManager(
        "d_equation_of_state_mixing_rules_manager",
        d_dim,
        d_num_species,
        MIXING_CLOSURE_MODEL::NO_MODEL,
        equation_of_state_mixing_rules_db,
        d_equation_of_state_str));
    
    d_equation_of_state_mixing_rules =
        d_equation_of_state_mixing_rules_manager->getEquationOfStateMixingRules();
    
    std::vector<double*> thermo_properties_ptr;
    
    const int num_thermo_properties = d_equation_of_state_mixing_rules->
        getNumberOfSpeciesThermodynamicProperties();
    
    thermo_properties_ptr.reserve(num_thermo_properties);
    d_thermo_properties.resize(num_thermo_properties);
    
    for (int ti = 0; ti < num_thermo_properties; ti++)
    {
        thermo_properties_ptr.push_back(&d_thermo_properties[ti]);
    }
    
    d_equation_of_state_mixing_rules->getSpeciesThermodynamicProperties(
        thermo_properties_ptr,
        0);
    
    /*
     * Initialize d_equation_of_shear_viscosity_mixing_rules_manager and get the equation of
     * shear viscosity mixing rules object.
     */
    
    if ((flow_model_db->keyExists("equation_of_shear_viscosity")) ||
        (flow_model_db->keyExists("d_equation_of_shear_viscosity_str")))
    {
        if (flow_model_db->keyExists("equation_of_shear_viscosity"))
        {
            d_equation_of_shear_viscosity_str =
                flow_model_db->getString("equation_of_shear_viscosity");
        }
        else if (flow_model_db->keyExists("d_equation_of_shear_viscosity_str"))
        {
            d_equation_of_shear_viscosity_str =
                flow_model_db->getString("d_equation_of_shear_viscosity_str");
        }
        
        boost::shared_ptr<tbox::Database> equation_of_shear_viscosity_mixing_rules_db;
        
        if (flow_model_db->keyExists("Equation_of_shear_viscosity_mixing_rules"))
        {
            equation_of_shear_viscosity_mixing_rules_db =
                flow_model_db->getDatabase("Equation_of_shear_viscosity_mixing_rules");
        }
        else if (flow_model_db->keyExists("d_equation_of_shear_viscosity_mixing_rules_db"))
        {
            equation_of_shear_viscosity_mixing_rules_db =
                flow_model_db->getDatabase("d_equation_of_shear_viscosity_mixing_rules_db");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "No key 'Equation_of_shear_viscosity_mixing_rules'/"
                << "'d_equation_of_shear_viscosity_mixing_rules_db' found in data for flow model"
                << std::endl);
        }
        
        d_equation_of_shear_viscosity_mixing_rules_manager.reset(
            new EquationOfShearViscosityMixingRulesManager(
                "d_equation_of_shear_viscosity_mixing_rules_manager",
                d_dim,
                d_num_species,
                MIXING_CLOSURE_MODEL::NO_MODEL,
                equation_of_shear_viscosity_mixing_rules_db,
                d_equation_of_shear_viscosity_str));
        
        d_equation_of_shear_viscosity_mixing_rules =
            d_equation_of_shear_viscosity_mixing_rules_manager->
                getEquationOfShearViscosityMixingRules();
        
        std::vector<double*> molecular_properties_ptr;
        
        const int num_molecular_properties = d_equation_of_shear_viscosity_mixing_rules->
            getNumberOfSpeciesMolecularProperties();
        
        molecular_properties_ptr.reserve(num_molecular_properties);
        d_molecular_properties_shear_viscosity.resize(num_molecular_properties);
        
        for (int ti = 0; ti < num_molecular_properties; ti++)
        {
            molecular_properties_ptr.push_back(&d_molecular_properties_shear_viscosity[ti]);
        }
        
        d_equation_of_shear_viscosity_mixing_rules->getSpeciesMolecularProperties(
            molecular_properties_ptr,
            0);
    }
    
    /*
     * Initialize d_equation_of_bulk_viscosity_mixing_rules_manager and get the equation of
     * bulk viscosity mixing rules object.
     */
    
    if ((flow_model_db->keyExists("equation_of_bulk_viscosity")) ||
        (flow_model_db->keyExists("d_equation_of_bulk_viscosity_str")))
    {
        if (flow_model_db->keyExists("equation_of_bulk_viscosity"))
        {
            d_equation_of_bulk_viscosity_str =
                flow_model_db->getString("equation_of_bulk_viscosity");
        }
        else if (flow_model_db->keyExists("d_equation_of_bulk_viscosity_str"))
        {
            d_equation_of_bulk_viscosity_str =
                flow_model_db->getString("d_equation_of_bulk_viscosity_str");
        }
        
        boost::shared_ptr<tbox::Database> equation_of_bulk_viscosity_mixing_rules_db;
        
        if (flow_model_db->keyExists("Equation_of_bulk_viscosity_mixing_rules"))
        {
            equation_of_bulk_viscosity_mixing_rules_db =
                flow_model_db->getDatabase("Equation_of_bulk_viscosity_mixing_rules");
        }
        else if (flow_model_db->keyExists("d_equation_of_bulk_viscosity_mixing_rules_db"))
        {
            equation_of_bulk_viscosity_mixing_rules_db =
                flow_model_db->getDatabase("d_equation_of_bulk_viscosity_mixing_rules_db");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "No key 'Equation_of_bulk_viscosity_mixing_rules'/"
                << "'d_equation_of_bulk_viscosity_mixing_rules_db' found in data for flow model"
                << std::endl);
        }
        
        d_equation_of_bulk_viscosity_mixing_rules_manager.reset(
            new EquationOfBulkViscosityMixingRulesManager(
                "d_equation_of_bulk_viscosity_mixing_rules_manager",
                d_dim,
                d_num_species,
                MIXING_CLOSURE_MODEL::NO_MODEL,
                equation_of_bulk_viscosity_mixing_rules_db,
                d_equation_of_bulk_viscosity_str));
        
        d_equation_of_bulk_viscosity_mixing_rules =
            d_equation_of_bulk_viscosity_mixing_rules_manager->
                getEquationOfBulkViscosityMixingRules();
        
        std::vector<double*> molecular_properties_ptr;
        
        const int num_molecular_properties = d_equation_of_bulk_viscosity_mixing_rules->
            getNumberOfSpeciesMolecularProperties();
        
        molecular_properties_ptr.reserve(num_molecular_properties);
        d_molecular_properties_bulk_viscosity.resize(num_molecular_properties);
        
        for (int ti = 0; ti < num_molecular_properties; ti++)
        {
            molecular_properties_ptr.push_back(&d_molecular_properties_bulk_viscosity[ti]);
        }
        
        d_equation_of_bulk_viscosity_mixing_rules->getSpeciesMolecularProperties(
            molecular_properties_ptr,
            0);
    }
    
    /*
     * Initialize d_equation_of_thermal_conductivity_mixing_rules_manager and get the equation of
     * thermal conductivity mixing rules object.
     */
    
    if ((flow_model_db->keyExists("equation_of_thermal_conductivity")) ||
        (flow_model_db->keyExists("d_equation_of_thermal_conductivity_str")))
    {
        if (flow_model_db->keyExists("equation_of_thermal_conductivity"))
        {
            d_equation_of_thermal_conductivity_str =
                flow_model_db->getString("equation_of_thermal_conductivity");
        }
        else if (flow_model_db->keyExists("d_equation_of_thermal_conductivity_str"))
        {
            d_equation_of_thermal_conductivity_str =
                flow_model_db->getString("d_equation_of_thermal_conductivity_str");
        }
        
        boost::shared_ptr<tbox::Database> equation_of_thermal_conductivity_mixing_rules_db;
        
        if (flow_model_db->keyExists("Equation_of_thermal_conductivity_mixing_rules"))
        {
            equation_of_thermal_conductivity_mixing_rules_db =
                flow_model_db->getDatabase("Equation_of_thermal_conductivity_mixing_rules");
        }
        else if (flow_model_db->keyExists("d_equation_of_thermal_conductivity_mixing_rules_db"))
        {
            equation_of_thermal_conductivity_mixing_rules_db =
                flow_model_db->getDatabase("d_equation_of_thermal_conductivity_mixing_rules_db");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "No key 'Equation_of_thermal_conductivity_mixing_rules'/"
                << "'d_equation_of_thermal_conductivity_mixing_rules_db' found in data for flow model"
                << std::endl);
        }
        
        d_equation_of_thermal_conductivity_mixing_rules_manager.reset(
            new EquationOfThermalConductivityMixingRulesManager(
                "d_equation_of_thermal_conductivity_mixing_rules_manager",
                d_dim,
                d_num_species,
                MIXING_CLOSURE_MODEL::NO_MODEL,
                equation_of_thermal_conductivity_mixing_rules_db,
                d_equation_of_thermal_conductivity_str));
        
        d_equation_of_thermal_conductivity_mixing_rules =
            d_equation_of_thermal_conductivity_mixing_rules_manager->
                getEquationOfThermalConductivityMixingRules();
        
        std::vector<double*> molecular_properties_ptr;
        
        const int num_molecular_properties = d_equation_of_thermal_conductivity_mixing_rules->
            getNumberOfSpeciesMolecularProperties();
        
        molecular_properties_ptr.reserve(num_molecular_properties);
        d_molecular_properties_thermal_conductivity.resize(num_molecular_properties);
        
        for (int ti = 0; ti < num_molecular_properties; ti++)
        {
            molecular_properties_ptr.push_back(&d_molecular_properties_thermal_conductivity[ti]);
        }
        
        d_equation_of_thermal_conductivity_mixing_rules->getSpeciesMolecularProperties(
            molecular_properties_ptr,
            0);
    }
    
    /*
     * Initialize Riemann solver object.
     */
    d_flow_model_riemann_solver.reset(new FlowModelRiemannSolverSingleSpecies(
        "d_flow_model_riemann_solver",
        d_dim,
        d_grid_geometry,
        d_num_species));
    
    /*
     * Initialize basic utilities object.
     */
    d_flow_model_basic_utilities.reset(new FlowModelBasicUtilitiesSingleSpecies(
        "d_flow_model_basic_utilities",
        d_dim,
        d_grid_geometry,
        d_num_species,
        d_equation_of_state_mixing_rules));
    
    /*
     * Initialize diffusive flux utilities object.
     */
    d_flow_model_diffusive_flux_utilities.reset(new FlowModelDiffusiveFluxUtilitiesSingleSpecies(
        "d_flow_model_diffusive_flux_utilities",
        d_dim,
        d_grid_geometry,
        d_num_species,
        d_equation_of_shear_viscosity_mixing_rules,
        d_equation_of_bulk_viscosity_mixing_rules,
        d_equation_of_thermal_conductivity_mixing_rules));
    
    /*
     * Initialize statistics utilities object.
     */
    d_flow_model_statistics_utilities.reset(new FlowModelStatisticsUtilitiesSingleSpecies(
        "d_flow_model_statistics_utilities",
        d_dim,
        d_grid_geometry,
        d_num_species,
        flow_model_db,
        d_equation_of_state_mixing_rules,
        d_equation_of_shear_viscosity_mixing_rules,
        d_equation_of_bulk_viscosity_mixing_rules,
        d_equation_of_thermal_conductivity_mixing_rules));
    
    /*
     * Initialize boundary utilities object.
     */
    d_flow_model_boundary_utilities.reset(
        new FlowModelBoundaryUtilitiesSingleSpecies(
            "d_flow_model_boundary_utilities",
            d_dim,
            d_num_species,
            d_num_eqn,
            d_equation_of_state_mixing_rules));
}


/*
 * Print all characteristics of the flow model class.
 */
void
FlowModelSingleSpecies::printClassData(std::ostream& os) const
{
    os << "\nPrint FlowModelSingleSpecies object..."
       << std::endl;
    
    os << std::endl;
    
    os << "FlowModelSingleSpecies: this = "
       << (FlowModelSingleSpecies *)this
       << std::endl;
    
    os << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
    d_equation_of_state_mixing_rules_manager->printClassData(os);
}


/*
 * Put the characteristics of the flow model class into the restart database.
 */
void
FlowModelSingleSpecies::putToRestart(
    const boost::shared_ptr<tbox::Database>& restart_db) const
{
    /*
     * Put the properties of d_equation_of_state_mixing_rules into the restart database.
     */
    
    restart_db->putString("d_equation_of_state_str", d_equation_of_state_str);
    
    boost::shared_ptr<tbox::Database> restart_equation_of_state_mixing_rules_db =
        restart_db->putDatabase("d_equation_of_state_mixing_rules_db");
    d_equation_of_state_mixing_rules->putToRestart(restart_equation_of_state_mixing_rules_db);
    
    /*
     * Put the properties of d_equation_of_shear_viscosity_mixing_rules into the restart database.
     */
    
    if (d_equation_of_shear_viscosity_mixing_rules)
    {
        restart_db->putString("d_equation_of_shear_viscosity_str", d_equation_of_shear_viscosity_str);
        
        boost::shared_ptr<tbox::Database> restart_equation_of_shear_viscosity_mixing_rules_db =
            restart_db->putDatabase("d_equation_of_shear_viscosity_mixing_rules_db");
        d_equation_of_shear_viscosity_mixing_rules->
            putToRestart(restart_equation_of_shear_viscosity_mixing_rules_db);
    }
    
    /*
     * Put the properties of d_equation_of_bulk_viscosity_mixing_rules into the restart database.
     */
    
    if (d_equation_of_bulk_viscosity_mixing_rules)
    {
        restart_db->putString("d_equation_of_bulk_viscosity_str", d_equation_of_bulk_viscosity_str);
        
        boost::shared_ptr<tbox::Database> restart_equation_of_bulk_viscosity_mixing_rules_db =
            restart_db->putDatabase("d_equation_of_bulk_viscosity_mixing_rules_db");
        d_equation_of_bulk_viscosity_mixing_rules->
            putToRestart(restart_equation_of_bulk_viscosity_mixing_rules_db);
    }
    
    /*
     * Put the properties of d_equation_of_thermal_conductivity_mixing_rules into the restart database.
     */
    
    if (d_equation_of_thermal_conductivity_mixing_rules)
    {
        restart_db->putString("d_equation_of_thermal_conductivity_str", d_equation_of_thermal_conductivity_str);
        
        boost::shared_ptr<tbox::Database> restart_equation_of_thermal_conductivity_mixing_rules_db =
            restart_db->putDatabase("d_equation_of_thermal_conductivity_mixing_rules_db");
        d_equation_of_thermal_conductivity_mixing_rules->
            putToRestart(restart_equation_of_thermal_conductivity_mixing_rules_db);
    }
    
    /*
     * Put the properties of d_flow_model_statistics_utilities into the restart database.
     */
    d_flow_model_statistics_utilities->putToRestart(restart_db);
}


/*
 * Register the conservative variables.
 */
void
FlowModelSingleSpecies::registerConservativeVariables(
    RungeKuttaLevelIntegrator* integrator,
    const hier::IntVector& num_ghosts,
    const hier::IntVector& num_ghosts_intermediate)
{
    integrator->registerVariable(
        s_variable_density,
        num_ghosts,
        num_ghosts_intermediate,
        RungeKuttaLevelIntegrator::TIME_DEP,
        d_grid_geometry,
        "CONSERVATIVE_COARSEN",
        "CONSERVATIVE_LINEAR_REFINE");
    
    integrator->registerVariable(
        s_variable_momentum,
        num_ghosts,
        num_ghosts_intermediate,
        RungeKuttaLevelIntegrator::TIME_DEP,
        d_grid_geometry,
        "CONSERVATIVE_COARSEN",
        "CONSERVATIVE_LINEAR_REFINE");
    
    integrator->registerVariable(
        s_variable_total_energy,
        num_ghosts,
        num_ghosts_intermediate,
        RungeKuttaLevelIntegrator::TIME_DEP,
        d_grid_geometry,
        "CONSERVATIVE_COARSEN",
        "CONSERVATIVE_LINEAR_REFINE");
}


/*
 * Get the names of conservative variables.
 */
std::vector<std::string>
FlowModelSingleSpecies::getNamesOfConservativeVariables(bool have_underscores)
{
    std::vector<std::string> names;
    names.reserve(3);
    
    if (have_underscores)
    {
        names.push_back("density");
        names.push_back("momentum");
        names.push_back("total_energy");
    }
    else
    {
        names.push_back("density");
        names.push_back("momentum");
        names.push_back("total energy");
    }
    
    return names;
}


/*
 * Get the names of primitive variables.
 */
std::vector<std::string>
FlowModelSingleSpecies::getNamesOfPrimitiveVariables(bool have_underscores)
{
    NULL_USE(have_underscores);
    
    std::vector<std::string> names;
    names.reserve(3);
    
    names.push_back("density");
    names.push_back("velocity");
    names.push_back("pressure");
    
    return names;
}


/*
 * Get the variable types of conservative variables.
 */
std::vector<std::string>
FlowModelSingleSpecies::getVariableTypesOfConservativeVariables()
{
    std::vector<std::string> types;
    types.reserve(3);
    
    types.push_back("SCALAR");
    types.push_back("VECTOR");
    types.push_back("SCALAR");
    
    return types;
}


/*
 * Get the variable types of primitive variables.
 */
std::vector<std::string>
FlowModelSingleSpecies::getVariableTypesOfPrimitiveVariables()
{
    std::vector<std::string> types;
    types.reserve(3);
    
    types.push_back("SCALAR");
    types.push_back("VECTOR");
    types.push_back("SCALAR");
    
    return types;
}


/*
 * Get the conservative variables.
 */
std::vector<boost::shared_ptr<pdat::CellVariable<double> > >
FlowModelSingleSpecies::getConservativeVariables()
{
    std::vector<boost::shared_ptr<pdat::CellVariable<double> > > conservative_variables;
    conservative_variables.reserve(3);
    
    conservative_variables.push_back(s_variable_density);
    conservative_variables.push_back(s_variable_momentum);
    conservative_variables.push_back(s_variable_total_energy);
    
    return conservative_variables;
}


/*
 * Register a patch with a data context.
 */
void
FlowModelSingleSpecies::registerPatchWithDataContext(
    const hier::Patch& patch,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
    // Check whether the patch is already unregistered.
    if (d_patch)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::registerPatchWithDataContext()\n"
            << "The patch is not yet unregistered."
            << std::endl);
    }
    
    d_patch = &patch;
    
    setDataContext(data_context);
    
    /*
     * Set the number of ghost cells of conservative variables.
     */
    
    boost::shared_ptr<pdat::CellData<double> > data_density(
        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
            d_patch->getPatchData(s_variable_density, getDataContext())));
    
    d_num_ghosts = data_density->getGhostCellWidth();
    
    /*
     * Set the interior and ghost boxes with their dimensions for the conservative cell variables.
     */
    
    d_interior_box = d_patch->getBox();
    d_interior_dims = d_interior_box.numberCells();
    
    d_ghost_box = d_interior_box;
    d_ghost_box.grow(d_num_ghosts);
    d_ghostcell_dims = d_ghost_box.numberCells();
}


/*
 * Register different derived variables in the registered patch. The derived variables to be registered
 * are given as entires in a map of the variable name to the number of sub-ghost cells required.
 * If the variable to be registered is one of the conservative variable, the corresponding entry
 * in the map is ignored.
 */
void
FlowModelSingleSpecies::registerDerivedVariables(
    const std::unordered_map<std::string, hier::IntVector>& num_subghosts_of_data)

{
    // Check whether a patch is already registered.
    if (!d_patch)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::registerDerivedVariables()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    // Check whether all or part of derived cell data is alredy computed.
    if (d_derived_cell_data_computed)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::registerDerivedVariables()\n"
            << "Derived cell data is already computed."
            << std::endl);
    }
    
    for (std::unordered_map<std::string, hier::IntVector>::const_iterator it = num_subghosts_of_data.begin();
         it != num_subghosts_of_data.end();
         it++)
    {
        if ((it->second < hier::IntVector::getZero(d_dim)) ||
            (it->second > d_num_ghosts))
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::registerDerivedVariables()\n"
                << "The number of sub-ghost cells of variables '"
                << it->first
                << "' is not between zero and d_num_ghosts."
                << std::endl);
        }
    }
    
    if (num_subghosts_of_data.find("VELOCITY") != num_subghosts_of_data.end())
    {
        setNumberOfSubGhosts(
            num_subghosts_of_data.find("VELOCITY")->second,
            "VELOCITY",
            "VELOCITY");
    }
    
    if (num_subghosts_of_data.find("INTERNAL_ENERGY") != num_subghosts_of_data.end())
    {
        setNumberOfSubGhosts(
            num_subghosts_of_data.find("INTERNAL_ENERGY")->second,
            "INTERNAL_ENERGY",
            "INTERNAL_ENERGY");
    }
    
    if (num_subghosts_of_data.find("PRESSURE") != num_subghosts_of_data.end())
    {
        setNumberOfSubGhosts(
            num_subghosts_of_data.find("PRESSURE")->second,
            "PRESSURE",
            "PRESSURE");
    }
    
    if (num_subghosts_of_data.find("SOUND_SPEED") != num_subghosts_of_data.end())
    {
        setNumberOfSubGhosts(
            num_subghosts_of_data.find("SOUND_SPEED")->second,
            "SOUND_SPEED",
            "SOUND_SPEED");
    }
    
    if (num_subghosts_of_data.find("TEMPERATURE") != num_subghosts_of_data.end())
    {
        setNumberOfSubGhosts(
            num_subghosts_of_data.find("TEMPERATURE")->second,
            "TEMPERATURE",
            "TEMPERATURE");
    }
    
    if (num_subghosts_of_data.find("CONVECTIVE_FLUX_X") != num_subghosts_of_data.end())
    {
        setNumberOfSubGhosts(
            num_subghosts_of_data.find("CONVECTIVE_FLUX_X")->second,
            "CONVECTIVE_FLUX_X",
            "CONVECTIVE_FLUX_X");
    }
    
    if (num_subghosts_of_data.find("CONVECTIVE_FLUX_Y") != num_subghosts_of_data.end())
    {
        setNumberOfSubGhosts(
            num_subghosts_of_data.find("CONVECTIVE_FLUX_Y")->second,
            "CONVECTIVE_FLUX_Y",
            "CONVECTIVE_FLUX_Y");
    }
    
    if (num_subghosts_of_data.find("CONVECTIVE_FLUX_Z") != num_subghosts_of_data.end())
    {
        setNumberOfSubGhosts(
            num_subghosts_of_data.find("CONVECTIVE_FLUX_Z")->second,
            "CONVECTIVE_FLUX_Z",
            "CONVECTIVE_FLUX_Z");
    }
    
    if (num_subghosts_of_data.find("PRIMITIVE_VARIABLES") != num_subghosts_of_data.end())
    {
        setNumberOfSubGhosts(
            num_subghosts_of_data.find("PRIMITIVE_VARIABLES")->second,
            "PRIMITIVE_VARIABLES",
            "PRIMITIVE_VARIABLES");
    }
    
    if (num_subghosts_of_data.find("MAX_WAVE_SPEED_X") != num_subghosts_of_data.end())
    {
        setNumberOfSubGhosts(
            num_subghosts_of_data.find("MAX_WAVE_SPEED_X")->second,
            "MAX_WAVE_SPEED_X",
            "MAX_WAVE_SPEED_X");
    }
    
    if (num_subghosts_of_data.find("MAX_WAVE_SPEED_Y") != num_subghosts_of_data.end())
    {
        setNumberOfSubGhosts(
            num_subghosts_of_data.find("MAX_WAVE_SPEED_Y")->second,
            "MAX_WAVE_SPEED_Y",
            "MAX_WAVE_SPEED_Y");
    }
    
    if (num_subghosts_of_data.find("MAX_WAVE_SPEED_Z") != num_subghosts_of_data.end())
    {
        setNumberOfSubGhosts(
            num_subghosts_of_data.find("MAX_WAVE_SPEED_Z")->second,
            "MAX_WAVE_SPEED_Z",
            "MAX_WAVE_SPEED_Z");
    }
    
    if (num_subghosts_of_data.find("MAX_DIFFUSIVITY") != num_subghosts_of_data.end())
    {
        setNumberOfSubGhosts(
            num_subghosts_of_data.find("MAX_DIFFUSIVITY")->second,
            "MAX_DIFFUSIVITY",
            "MAX_DIFFUSIVITY");
    }
}


/*
 * Unregister the registered patch. The registered data context and all global derived
 * cell data in the patch are cleared.
 */
void
FlowModelSingleSpecies::unregisterPatch()
{
    // Check whether a patch is already registered.
    if (!d_patch)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::unregisterPatch()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    d_num_ghosts                      = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_velocity          = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_internal_energy   = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_pressure          = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_sound_speed       = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_temperature       = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_convective_flux_x = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_convective_flux_y = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_convective_flux_z = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_max_wave_speed_x  = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_max_wave_speed_y  = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_max_wave_speed_z  = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_max_diffusivity   = -hier::IntVector::getOne(d_dim);
    
    d_interior_box                   = hier::Box::getEmptyBox(d_dim);
    d_ghost_box                      = hier::Box::getEmptyBox(d_dim);
    d_subdomain_box                  = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_velocity          = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_internal_energy   = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_pressure          = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_sound_speed       = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_temperature       = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_convective_flux_x = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_convective_flux_y = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_convective_flux_z = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_max_wave_speed_x  = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_max_wave_speed_y  = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_max_wave_speed_z  = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_max_diffusivity   = hier::Box::getEmptyBox(d_dim);
    
    d_interior_dims                       = hier::IntVector::getZero(d_dim);
    d_ghostcell_dims                      = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_velocity          = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_internal_energy   = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_pressure          = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_sound_speed       = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_temperature       = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_convective_flux_x = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_convective_flux_y = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_convective_flux_z = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_max_wave_speed_x  = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_max_wave_speed_y  = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_max_wave_speed_z  = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_max_diffusivity   = hier::IntVector::getZero(d_dim);
    
    d_data_velocity.reset();
    d_data_internal_energy.reset();
    d_data_pressure.reset();
    d_data_sound_speed.reset();
    d_data_temperature.reset();
    d_data_convective_flux_x.reset();
    d_data_convective_flux_y.reset();
    d_data_convective_flux_z.reset();
    d_data_max_wave_speed_x.reset();
    d_data_max_wave_speed_y.reset();
    d_data_max_wave_speed_z.reset();
    d_data_max_diffusivity.reset();
    
    d_cell_data_velocity_computed          = false;
    d_cell_data_internal_energy_computed   = false;
    d_cell_data_pressure_computed          = false;
    d_cell_data_sound_speed_computed       = false;
    d_cell_data_temperature_computed       = false;
    d_cell_data_convective_flux_x_computed = false;
    d_cell_data_convective_flux_y_computed = false;
    d_cell_data_convective_flux_z_computed = false;
    d_cell_data_max_wave_speed_x_computed  = false;
    d_cell_data_max_wave_speed_y_computed  = false;
    d_cell_data_max_wave_speed_z_computed  = false;
    d_cell_data_max_diffusivity_computed   = false;
    
    d_flow_model_diffusive_flux_utilities->clearCellData();
    
    d_derived_cell_data_computed = false;
    
    d_patch = nullptr;
    clearDataContext();
}


/*
 * Allocate memory for cell data of different registered derived variables.
 */
void
FlowModelSingleSpecies::allocateMemoryForDerivedCellData()
{
    if (d_num_subghosts_velocity > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_velocity_computed)
        {
            if (!d_data_velocity)
            {
                // Create the cell data of velocity.
                d_data_velocity.reset(
                    new pdat::CellData<double>(d_interior_box, d_dim.getValue(), d_num_subghosts_velocity));
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::allocateMemoryForDerivedCellData()\n"
                << "Cell data of 'VELOCITY' is aleady computed."
                << std::endl);
        }
    }
    
    if (d_num_subghosts_internal_energy > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_internal_energy_computed)
        {
            if (!d_data_internal_energy)
            {
                // Create the cell data of internal energy.
                d_data_internal_energy.reset(
                    new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_internal_energy));
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::allocateMemoryForDerivedCellData()\n"
                << "Cell data of 'INTERNAL_ENERGY' is aleady computed."
                << std::endl);
        }
    }
    
    if (d_num_subghosts_pressure > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_pressure_computed)
        {
            if (!d_data_pressure)
            {
                // Create the cell data of pressure.
                d_data_pressure.reset(
                    new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_pressure));
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::allocateMemoryForDerivedCellData()\n"
                << "Cell data of 'PRESSURE' is aleady computed."
                << std::endl);
        }
    }
    
    if (d_num_subghosts_sound_speed > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_sound_speed_computed)
        {
            if (!d_data_sound_speed)
            {
                // Create the cell data of sound speed.
                d_data_sound_speed.reset(
                    new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_sound_speed));
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::allocateMemoryForDerivedCellData()\n"
                << "Cell data of 'SOUND_SPEED' is aleady computed."
                << std::endl);
        }
    }
    
    if (d_num_subghosts_temperature > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_temperature_computed)
        {
            if (!d_data_temperature)
            {
                // Create the cell data of temperature.
                d_data_temperature.reset(
                    new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_temperature));
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::allocateMemoryForDerivedCellData()\n"
                << "Cell data of 'TEMPERATURE' is aleady computed."
                << std::endl);
        }
    }
    
    if (d_num_subghosts_convective_flux_x > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_convective_flux_x_computed)
        {
            if (!d_data_convective_flux_x)
            {
                // Create the cell data of convective flux in the x-direction.
                d_data_convective_flux_x.reset(
                    new pdat::CellData<double>(d_interior_box, d_num_eqn, d_num_subghosts_convective_flux_x));
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::allocateMemoryForDerivedCellData()\n"
                << "Cell data of 'CONVECTIVE_FLUX_X' is aleady computed."
                << std::endl);
        }
    }
    
    if (d_num_subghosts_convective_flux_y > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_convective_flux_y_computed)
        {
            if (!d_data_convective_flux_y)
            {
                // Create the cell data of convective flux in the y-direction.
                d_data_convective_flux_y.reset(
                    new pdat::CellData<double>(d_interior_box, d_num_eqn, d_num_subghosts_convective_flux_y));
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::allocateMemoryForDerivedCellData()\n"
                << "Cell data of 'CONVECTIVE_FLUX_Y' is aleady computed."
                << std::endl);
        }
    }
    
    if (d_num_subghosts_convective_flux_z > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_convective_flux_z_computed)
        {
            if (!d_data_convective_flux_z)
            {
                // Create the cell data of convective flux in the z-direction.
                d_data_convective_flux_z.reset(
                    new pdat::CellData<double>(d_interior_box, d_num_eqn, d_num_subghosts_convective_flux_z));
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::allocateMemoryForDerivedCellData()\n"
                << "Cell data of 'CONVECTIVE_FLUX_Z' is aleady computed."
                << std::endl);
        }
    }
    
    if (d_num_subghosts_max_wave_speed_x > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_max_wave_speed_x_computed)
        {
            if (!d_data_max_wave_speed_x)
            {
                // Create the cell data of maximum wave speed in the x-direction.
                d_data_max_wave_speed_x.reset(
                    new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_max_wave_speed_x));
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::allocateMemoryForDerivedCellData()\n"
                << "Cell data of 'MAX_WAVE_SPEED_X' is aleady computed."
                << std::endl);
        }
    }
    
    if (d_num_subghosts_max_wave_speed_y > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_max_wave_speed_y_computed)
        {
            if (!d_data_max_wave_speed_y)
            {
                // Create the cell data of maximum wave speed in the y-direction.
                d_data_max_wave_speed_y.reset(
                    new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_max_wave_speed_y));
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::allocateMemoryForDerivedCellData()\n"
                << "Cell data of 'MAX_WAVE_SPEED_Y' is aleady computed."
                << std::endl);
        }
    }
    
    if (d_num_subghosts_max_wave_speed_z > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_max_wave_speed_z_computed)
        {
            if (!d_data_max_wave_speed_z)
            {
                // Create the cell data of maximum wave speed in the z-direction.
                d_data_max_wave_speed_z.reset(
                    new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_max_wave_speed_z));
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::allocateMemoryForDerivedCellData()\n"
                << "Cell data of 'MAX_WAVE_SPEED_Z' is aleady computed."
                << std::endl);
        }
    }
    
    if (d_num_subghosts_max_diffusivity > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_max_diffusivity_computed)
        {
            if (!d_data_max_diffusivity)
            {
                // Create the cell data of maximum diffusivity.
                d_data_max_diffusivity.reset(
                    new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_max_diffusivity));
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::allocateMemoryForDerivedCellData()\n"
                << "Cell data of 'MAX_DIFFUSIVITY' is aleady computed."
                << std::endl);
        }
    }
}


/*
 * Compute the cell data of different registered derived variables with the registered data context.
 */
void
FlowModelSingleSpecies::computeDerivedCellData()
{
    // Check whether a patch is already registered.
    if (!d_patch)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::computeDerivedCellData()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    /*
     * Set the boxes and their dimensions for the derived cell variables.
     */
    if (!d_derived_cell_data_computed)
    {
        setDerivedCellVariableGhostBoxes();
    }
    
    // Compute the velocity cell data.
    if (d_num_subghosts_velocity > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_velocity_computed)
        {
            computeCellDataOfVelocity(
                d_subdomain_box);
        }
    }
    
    // Compute the internal energy cell data.
    if (d_num_subghosts_internal_energy > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_internal_energy_computed)
        {
            computeCellDataOfInternalEnergyWithVelocity(
                d_subdomain_box);
        }
    }
    
    // Compute the pressure cell data.
    if (d_num_subghosts_pressure > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_pressure_computed)
        {
            computeCellDataOfPressureWithInternalEnergy(
                d_subdomain_box);
        }
    }
    
    // Compute the sound speed cell data.
    if (d_num_subghosts_sound_speed > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_sound_speed_computed)
        {
            computeCellDataOfSoundSpeedWithPressure(
                d_subdomain_box);
        }
    }
    
    // Compute the temperature cell data.
    if (d_num_subghosts_temperature > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_temperature_computed)
        {
            computeCellDataOfTemperatureWithPressure(
                d_subdomain_box);
        }
    }
    
    // Compute the x-direction convective flux cell data.
    if (d_num_subghosts_convective_flux_x > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_convective_flux_x_computed)
        {
            computeCellDataOfConvectiveFluxWithVelocityAndPressure(
                DIRECTION::X_DIRECTION,
                d_subdomain_box);
        }
    }
    
    // Compute the y-direction convective flux cell data.
    if (d_num_subghosts_convective_flux_y > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_convective_flux_y_computed)
        {
            computeCellDataOfConvectiveFluxWithVelocityAndPressure(
                DIRECTION::Y_DIRECTION,
                d_subdomain_box);
        }
    }
    
    // Compute the z-direction convective flux cell data.
    if (d_num_subghosts_convective_flux_z > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_convective_flux_z_computed)
        {
            computeCellDataOfConvectiveFluxWithVelocityAndPressure(
                DIRECTION::Z_DIRECTION,
                d_subdomain_box);
        }
    }
    
    // Compute the x-direction maximum wave speed cell data.
    if (d_num_subghosts_max_wave_speed_x > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_max_wave_speed_x_computed)
        {
            computeCellDataOfMaxWaveSpeedWithVelocityAndSoundSpeed(
                DIRECTION::X_DIRECTION,
                d_subdomain_box);
        }
    }
    
    // Compute the y-direction maximum wave speed cell data.
    if (d_num_subghosts_max_wave_speed_y > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_max_wave_speed_y_computed)
        {
            computeCellDataOfMaxWaveSpeedWithVelocityAndSoundSpeed(
                DIRECTION::Y_DIRECTION,
                d_subdomain_box);
        }
    }
    
    // Compute the z-direction maximum wave speed cell data.
    if (d_num_subghosts_max_wave_speed_z > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_max_wave_speed_z_computed)
        {
            computeCellDataOfMaxWaveSpeedWithVelocityAndSoundSpeed(
                DIRECTION::Z_DIRECTION,
                d_subdomain_box);
        }
    }
    
    // Compute the maximum diffusivity cell data.
    if (d_num_subghosts_max_diffusivity > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_max_diffusivity_computed)
        {
            computeCellDataOfMaxDiffusivityWithPressureAndTemperature(
                d_subdomain_box);
        }
    }
    
    d_derived_cell_data_computed = true;
}


/*
 * Get the cell data of one cell variable in the registered patch.
 */
boost::shared_ptr<pdat::CellData<double> >
FlowModelSingleSpecies::getCellData(const std::string& variable_key)
{
    // Check whether a patch is already registered.
    if (!d_patch)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::getCellData()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    boost::shared_ptr<pdat::CellData<double> > cell_data;
    
    if (variable_key == "DENSITY")
    {
        cell_data = getCellDataOfDensity();
    }
    else if (variable_key == "MOMENTUM")
    {
        cell_data = getCellDataOfMomentum();
    }
    else if (variable_key == "TOTAL_ENERGY")
    {
        cell_data = getCellDataOfTotalEnergy();
    }
    else if (variable_key == "VELOCITY")
    {
        if (!d_cell_data_velocity_computed)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::getCellData()\n"
                << "Cell data of 'VELOCITY' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_velocity;
    }
    else if (variable_key == "INTERNAL_ENERGY")
    {
        if (!d_cell_data_internal_energy_computed)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::getCellData()\n"
                << "Cell data of 'INTERNAL_ENERGY' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_internal_energy;
    }
    else if (variable_key == "PRESSURE")
    {
        if (!d_cell_data_pressure_computed)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::getCellData()\n"
                << "Cell data of 'PRESSURE' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_pressure;
    }
    else if (variable_key == "SOUND_SPEED")
    {
        if (!d_cell_data_sound_speed_computed)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::getCellData()\n"
                << "Cell data of 'SOUND_SPEED' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_sound_speed;
    }
    else if (variable_key == "TEMPERATURE")
    {
        if (!d_cell_data_temperature_computed)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::getCellData()\n"
                << "Cell data of 'TEMPERATURE' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_temperature;
    }
    else if (variable_key == "CONVECTIVE_FLUX_X")
    {
        if (!d_cell_data_convective_flux_x_computed)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::getCellData()\n"
                << "Cell data of 'CONVECTIVE_FLUX_X' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_convective_flux_x;
    }
    else if (variable_key == "CONVECTIVE_FLUX_Y")
    {
        if (!d_cell_data_convective_flux_y_computed)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::getCellData()\n"
                << "Cell data of 'CONVECTIVE_FLUX_Y' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_convective_flux_y;
    }
    else if (variable_key == "CONVECTIVE_FLUX_Z")
    {
        if (!d_cell_data_convective_flux_z_computed)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::getCellData()\n"
                << "Cell data of 'CONVECTIVE_FLUX_Z' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_convective_flux_z;
    }
    else if (variable_key == "MAX_WAVE_SPEED_X")
    {
        if (!d_cell_data_max_wave_speed_x_computed)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::getCellData()\n"
                << "Cell data of 'MAX_WAVE_SPEED_X' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_max_wave_speed_x;
    }
    else if (variable_key == "MAX_WAVE_SPEED_Y")
    {
        if (!d_cell_data_max_wave_speed_y_computed)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::getCellData()\n"
                << "Cell data of 'MAX_WAVE_SPEED_Y' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_max_wave_speed_y;
    }
    else if (variable_key == "MAX_WAVE_SPEED_Z")
    {
        if (!d_cell_data_max_wave_speed_z_computed)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::getCellData()\n"
                << "Cell data of 'MAX_WAVE_SPEED_Z' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_max_wave_speed_z;
    }
    else if (variable_key == "MAX_DIFFUSIVITY")
    {
        if (!d_cell_data_max_diffusivity_computed)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::getCellData()\n"
                << "Cell data of 'MAX_DIFFUSIVITY' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_max_diffusivity;
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::getCellData()\n"
            << "Unknown cell data with variable_key = '" << variable_key
            << "' requested."
            << std::endl);
    }
    
    return cell_data;
}


/*
 * Get the cell data of different cell variables in the registered patch.
 */
std::vector<boost::shared_ptr<pdat::CellData<double> > >
FlowModelSingleSpecies::getCellData(
    const std::vector<std::string>& variable_keys)
{
    std::vector<boost::shared_ptr<pdat::CellData<double> > > cell_data(
        static_cast<int>(variable_keys.size()));
    
    for (int vi = 0; static_cast<int>(variable_keys.size()); vi++)
    {
        cell_data[vi] = getCellData(variable_keys[vi]);
    }
    
    return cell_data;
}


/*
 * Get the cell data of species cell variables in the registered patch.
 */
std::vector<boost::shared_ptr<pdat::CellData<double> > >
FlowModelSingleSpecies::getSpeciesCellData(
    const std::string& variable_key)
{
    TBOX_ERROR(d_object_name
        << ": Method FlowModelSingleSpecies::getSpeciesCellData() is not implemented yet!"
        << std::endl);
    
    std::vector<boost::shared_ptr<pdat::CellData<double> > > tmp;
    return tmp;
}


/*
 * Fill the cell data of conservative variables in the interior box with value zero.
 */
void
FlowModelSingleSpecies::fillCellDataOfConservativeVariablesWithZero()
{
    // Check whether a patch is already registered.
    if (!d_patch)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::fillCellDataOfConservativeVariablesWithZero()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    boost::shared_ptr<pdat::CellData<double> > data_density = getCellDataOfDensity();
    boost::shared_ptr<pdat::CellData<double> > data_momentum = getCellDataOfMomentum();
    boost::shared_ptr<pdat::CellData<double> > data_total_energy = getCellDataOfTotalEnergy();
    
    data_density->fillAll(double(0), d_interior_box);
    data_momentum->fillAll(double(0), d_interior_box);
    data_total_energy->fillAll(double(0), d_interior_box);
}


/*
 * Update the cell data of conservative variables in the interior box after time advancement.
 */
void
FlowModelSingleSpecies::updateCellDataOfConservativeVariables()
{
    // Check whether a patch is already registered.
    if (!d_patch)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::updateCellDataOfConservativeVariables()\n"
            << "No patch is registered yet."
            << std::endl);
    }
}


/*
 * Get the cell data of the conservative variables in the registered patch.
 */
std::vector<boost::shared_ptr<pdat::CellData<double> > >
FlowModelSingleSpecies::getCellDataOfConservativeVariables()
{
    // Check whether a patch is already registered.
    if (!d_patch)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::getCellDataOfConservativeVariables()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    std::vector<boost::shared_ptr<pdat::CellData<double> > > cell_data;
    cell_data.reserve(3);
    
    cell_data.push_back(getCellDataOfDensity());
    cell_data.push_back(getCellDataOfMomentum());
    cell_data.push_back(getCellDataOfTotalEnergy());
    
    return cell_data;
}


/*
 * Get the cell data of the primitive variables in the registered patch.
 */
std::vector<boost::shared_ptr<pdat::CellData<double> > >
FlowModelSingleSpecies::getCellDataOfPrimitiveVariables()
{
    // Check whether a patch is already registered.
    if (!d_patch)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::getCellDataOfPrimitiveVariables()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    std::vector<boost::shared_ptr<pdat::CellData<double> > > cell_data;
    cell_data.reserve(3);
    
    cell_data.push_back(getCellDataOfDensity());
    if (!d_data_velocity)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::getCellDataOfPrimitiveVariables()\n"
            << "Cell data of 'VELOCITY' is not registered/computed yet."
            << std::endl);
    }
    cell_data.push_back(d_data_velocity);
    if (!d_data_pressure)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::getCellDataOfPrimitiveVariables()\n"
            << "Cell data of 'PRESSURE' is not registered/computed yet."
            << std::endl);
    }
    cell_data.push_back(d_data_pressure);
    
    return cell_data;
}


/*
 * Compute derived plot quantities registered with the VisIt data writers from data that
 * is maintained on each patch in the hierarchy.
 */
bool
FlowModelSingleSpecies::packDerivedDataIntoDoubleBuffer(
    double* buffer,
    const hier::Patch& patch,
    const hier::Box& region,
    const std::string& variable_name,
    int depth_id,
    double simulation_time) const
{
    if (!d_plot_context)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "packDerivedDataIntoDoubleBuffer()\n"
            << "The plotting context is not set yet."
            << std::endl);
    }
    
    NULL_USE(simulation_time);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT((region * patch.getBox()).isSpatiallyEqual(region));
#endif
    
    bool data_on_patch = false;
    
    // Get the dimensions of the region.
    const hier::IntVector region_dims = region.numberCells();
    
    if (variable_name == "pressure")
    {
        boost::shared_ptr<pdat::CellData<double> > data_density(
            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                patch.getPatchData(s_variable_density, d_plot_context)));
        
        boost::shared_ptr<pdat::CellData<double> > data_momentum(
            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                patch.getPatchData(s_variable_momentum, d_plot_context)));
        
        boost::shared_ptr<pdat::CellData<double> > data_total_energy(
            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                patch.getPatchData(s_variable_total_energy, d_plot_context)));
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density);
        TBOX_ASSERT(data_momentum);
        TBOX_ASSERT(data_total_energy);
        TBOX_ASSERT(data_density->getGhostBox().isSpatiallyEqual(patch.getBox()));
        TBOX_ASSERT(data_momentum->getGhostBox().isSpatiallyEqual(patch.getBox()));
        TBOX_ASSERT(data_total_energy->getGhostBox().isSpatiallyEqual(patch.getBox()));
#endif
        
        // Get the dimensions of box that covers the data.
        const hier::Box data_box = data_density->getGhostBox();
        const hier::IntVector data_dims = data_box.numberCells();
        
        // Get the pointers to the conservative variables.
        const double* const rho   = data_density->getPointer();
        const double* const rho_u = data_momentum->getPointer(0);
        const double* const rho_v = d_dim > tbox::Dimension(1) ? data_momentum->getPointer(1) : NULL;
        const double* const rho_w = d_dim > tbox::Dimension(2) ? data_momentum->getPointer(2) : NULL;
        const double* const E     = data_total_energy->getPointer(0);
        
        size_t offset_data = data_box.offset(region.lower());
        
        // Get the thermodynamic properties of the species.
        std::vector<const double*> thermo_properties_ptr;
        thermo_properties_ptr.reserve(static_cast<int> (d_thermo_properties.size()));
        for (int ti = 0; ti < static_cast<int> (d_thermo_properties.size()); ti++)
        {
            thermo_properties_ptr.push_back(&d_thermo_properties[ti]);
        }
        
        if (d_dim == tbox::Dimension(1))
        {
            for (int i = 0; i < region_dims[0]; i++)
            {
                // Compute the linear indices.
                size_t idx_data = offset_data + i;
                size_t idx_region = i;
                
                const double epsilon = (E[idx_data] -
                    double(1)/double(2)*rho_u[idx_data]*rho_u[idx_data]/rho[idx_data])/rho[idx_data];
                
                buffer[idx_region] = d_equation_of_state_mixing_rules->getEquationOfState()->
                    getPressure(
                        &rho[idx_data],
                        &epsilon,
                        thermo_properties_ptr);
            }
        }
        else if (d_dim == tbox::Dimension(2))
        {
            for (int j = 0; j < region_dims[1]; j++)
            {
                for (int i = 0; i < region_dims[0]; i++)
                {
                    // Compute the linear indices.
                    size_t idx_data = offset_data + i +
                        j*data_dims[0];
                    
                    size_t idx_region = i +
                        j*region_dims[0];
                    
                    const double epsilon = (E[idx_data] - double(1)/double(2)*(rho_u[idx_data]*rho_u[idx_data] +
                        rho_v[idx_data]*rho_v[idx_data])/rho[idx_data])/rho[idx_data];
                    
                    buffer[idx_region] = d_equation_of_state_mixing_rules->getEquationOfState()->
                        getPressure(
                            &rho[idx_data],
                            &epsilon,
                            thermo_properties_ptr);
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
                        // Compute the linear indices.
                        size_t idx_data = offset_data + i +
                            j*data_dims[0] +
                            k*data_dims[0]*data_dims[1];
                        
                        size_t idx_region = i +
                            j*region_dims[0] +
                            k*region_dims[0]*region_dims[1];
                        
                        const double epsilon = (E[idx_data] - double(1)/double(2)*(rho_u[idx_data]*rho_u[idx_data] +
                            rho_v[idx_data]*rho_v[idx_data] + rho_w[idx_data]*rho_w[idx_data])
                                /rho[idx_data])/rho[idx_data];
                        
                        buffer[idx_region] = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getPressure(
                                &rho[idx_data],
                                &epsilon,
                                thermo_properties_ptr);
                    }
                }
            }
        }
        
        data_on_patch = true;
    }
    else if (variable_name == "sound speed")
    {
        boost::shared_ptr<pdat::CellData<double> > data_density(
            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                patch.getPatchData(s_variable_density, d_plot_context)));
        
        boost::shared_ptr<pdat::CellData<double> > data_momentum(
            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                patch.getPatchData(s_variable_momentum, d_plot_context)));
        
        boost::shared_ptr<pdat::CellData<double> > data_total_energy(
            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                patch.getPatchData(s_variable_total_energy, d_plot_context)));
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density);
        TBOX_ASSERT(data_momentum);
        TBOX_ASSERT(data_total_energy);
        TBOX_ASSERT(data_density->getGhostBox().isSpatiallyEqual(patch.getBox()));
        TBOX_ASSERT(data_momentum->getGhostBox().isSpatiallyEqual(patch.getBox()));
        TBOX_ASSERT(data_total_energy->getGhostBox().isSpatiallyEqual(patch.getBox()));
#endif
        
        // Get the dimensions of box that covers the data.
        const hier::Box data_box = data_density->getGhostBox();
        const hier::IntVector data_dims = data_box.numberCells();
        
        // Get the pointers to the conservative variables.
        const double* const rho   = data_density->getPointer();
        const double* const rho_u = data_momentum->getPointer(0);
        const double* const rho_v = d_dim > tbox::Dimension(1) ? data_momentum->getPointer(1) : NULL;
        const double* const rho_w = d_dim > tbox::Dimension(2) ? data_momentum->getPointer(2) : NULL;
        const double* const E     = data_total_energy->getPointer(0);
        
        size_t offset_data = data_box.offset(region.lower());
        
        // Get the thermodynamic properties of the species.
        std::vector<const double*> thermo_properties_ptr;
        thermo_properties_ptr.reserve(static_cast<int> (d_thermo_properties.size()));
        for (int ti = 0; ti < static_cast<int> (d_thermo_properties.size()); ti++)
        {
            thermo_properties_ptr.push_back(&d_thermo_properties[ti]);
        }
        
        if (d_dim == tbox::Dimension(1))
        {
            for (int i = 0; i < region_dims[0]; i++)
            {
                // Compute the linear indices.
                size_t idx_data = offset_data + i;
                size_t idx_region = i;
                
                const double epsilon = (E[idx_data] - double(1)/double(2)*rho_u[idx_data]*rho_u[idx_data]/
                    rho[idx_data])/rho[idx_data];
                
                const double p = d_equation_of_state_mixing_rules->getEquationOfState()->
                    getPressure(
                        &rho[idx_data],
                        &epsilon,
                        thermo_properties_ptr);
                
                buffer[idx_region] = d_equation_of_state_mixing_rules->getEquationOfState()->
                    getSoundSpeed(
                        &rho[idx_data],
                        &p,
                        thermo_properties_ptr);
            }
        }
        else if (d_dim == tbox::Dimension(2))
        {
            for (int j = 0; j < region_dims[1]; j++)
            {
                for (int i = 0; i < region_dims[0]; i++)
                {
                    // Compute the linear indices.
                    size_t idx_data = offset_data + i +
                        j*data_dims[0];
                    
                    size_t idx_region = i +
                        j*region_dims[0];
                    
                    const double epsilon = (E[idx_data] - double(1)/double(2)*(rho_u[idx_data]*rho_u[idx_data] +
                        rho_v[idx_data]*rho_v[idx_data])/rho[idx_data])/rho[idx_data];
                    
                    const double p = d_equation_of_state_mixing_rules->getEquationOfState()->
                        getPressure(
                            &rho[idx_data],
                            &epsilon,
                            thermo_properties_ptr);
                    
                    buffer[idx_region] = d_equation_of_state_mixing_rules->getEquationOfState()->
                        getSoundSpeed(
                            &rho[idx_data],
                            &p,
                            thermo_properties_ptr);
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
                        // Compute the linear indices.
                        size_t idx_data = offset_data + i +
                            j*data_dims[0] +
                            k*data_dims[0]*data_dims[1];
                        
                        size_t idx_region = i +
                            j*region_dims[0] +
                            k*region_dims[0]*region_dims[1];
                        
                        const double epsilon = (E[idx_data] - double(1)/double(2)*(rho_u[idx_data]*rho_u[idx_data] +
                            rho_v[idx_data]*rho_v[idx_data] + rho_w[idx_data]*rho_w[idx_data])
                                /rho[idx_data])/rho[idx_data];
                        
                        const double p = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getPressure(
                                &rho[idx_data],
                                &epsilon,
                                thermo_properties_ptr);
                        
                        buffer[idx_region] = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getSoundSpeed(
                                &rho[idx_data],
                                &p,
                                thermo_properties_ptr);
                    }
                }
            }
        }
        
        data_on_patch = true;
    }
    else if (variable_name == "velocity")
    {
        boost::shared_ptr<pdat::CellData<double> > data_density(
            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                patch.getPatchData(s_variable_density, d_plot_context)));
        
        boost::shared_ptr<pdat::CellData<double> > data_momentum(
            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                patch.getPatchData(s_variable_momentum, d_plot_context)));
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_density);
        TBOX_ASSERT(data_momentum);
        TBOX_ASSERT(data_density->getGhostBox().isSpatiallyEqual(patch.getBox()));
        TBOX_ASSERT(data_momentum->getGhostBox().isSpatiallyEqual(patch.getBox()));
#endif
        
        // Get the dimensions of box that covers the data.
        const hier::Box data_box = data_density->getGhostBox();
        const hier::IntVector data_dims = data_box.numberCells();
        
        // Get the pointers to the conservative variables.
        const double* const rho   = data_density->getPointer();
        const double* const m     = data_momentum->getPointer(depth_id);
        
        size_t offset_data = data_box.offset(region.lower());
        
        if (d_dim == tbox::Dimension(1))
        {
            for (int i = 0; i < region_dims[0]; i++)
            {
                // Compute the linear indices.
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
                    // Compute the linear indices.
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
                        // Compute the linear indices.
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
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "packDerivedDataIntoDoubleBuffer()\n"
            << "Unknown/unsupported variable '"
            << variable_name
            << "' required to compute."
            << std::endl);
    }
    
    return data_on_patch;
}


/*
 * Register the plotting quantities.
 */
#ifdef HAVE_HDF5
void
FlowModelSingleSpecies::registerPlotQuantities(
    const boost::shared_ptr<ExtendedVisItDataWriter>& visit_writer)
{
    if (!d_plot_context)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "registerPlotQuantities()\n"
            << "The plotting context is not set yet."
            << std::endl);
    }
    
    hier::VariableDatabase* vardb = hier::VariableDatabase::getDatabase();
    
    visit_writer->registerPlotQuantity(
        "density",
        "SCALAR",
        vardb->mapVariableAndContextToIndex(
           s_variable_density,
           d_plot_context));
    
    /*
    visit_writer->registerPlotQuantity(
        "momentum",
        "VECTOR",
        vardb->mapVariableAndContextToIndex(
           s_variable_momentum,
           d_plot_context));
    
    visit_writer->registerPlotQuantity(
        "total energy",
        "SCALAR",
        vardb->mapVariableAndContextToIndex(
           s_variable_total_energy,
           d_plot_context));
    */
    
    visit_writer->registerDerivedPlotQuantity("pressure",
        "SCALAR",
        this);
    
    visit_writer->registerDerivedPlotQuantity("sound speed",
        "SCALAR",
        this);
    
    visit_writer->registerDerivedPlotQuantity("velocity",
        "VECTOR",
        this);
}
#endif


/*
 * Set the number of sub-ghost cells of a variable.
 * This function can be called recursively if the variables are computed recursively.
 */
void
FlowModelSingleSpecies::setNumberOfSubGhosts(
    const hier::IntVector& num_subghosts,
    const std::string& variable_name,
    const std::string& parent_variable_name)
{
    NULL_USE(parent_variable_name);
    
    if (variable_name == "VELOCITY")
    {
        if (d_num_subghosts_velocity > -hier::IntVector::getOne(d_dim))
        {
            if (num_subghosts > d_num_subghosts_velocity)
            {
                /*
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::setNumberOfSubGhosts()\n"
                    << "Number of ghosts of '"
                    << parent_variable_name
                    << "' exceeds"
                    << " number of ghosts of '"
                    << variable_name
                    << "'."
                    << std::endl);
                */
                
                d_num_subghosts_velocity = num_subghosts;
            }
        }
        else
        {
            d_num_subghosts_velocity = num_subghosts;
        }
    }
    else if (variable_name == "INTERNAL_ENERGY")
    {
        if (d_num_subghosts_internal_energy > -hier::IntVector::getOne(d_dim))
        {
            if (num_subghosts > d_num_subghosts_internal_energy)
            {
                /*
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::setNumberOfSubGhosts()\n"
                    << "Number of ghosts of '"
                    << parent_variable_name
                    << "' exceeds"
                    << " number of ghosts of '"
                    << variable_name
                    << "'."
                    << std::endl);
                */
                
                d_num_subghosts_internal_energy = num_subghosts;
            }
        }
        else
        {
            d_num_subghosts_internal_energy = num_subghosts;
        }
        
        setNumberOfSubGhosts(num_subghosts, "VELOCITY", parent_variable_name);
    }
    else if (variable_name == "PRESSURE")
    {
        if (d_num_subghosts_pressure > -hier::IntVector::getOne(d_dim))
        {
            if (num_subghosts > d_num_subghosts_pressure)
            {
                /*
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::setNumberOfSubGhosts()\n"
                    << "Number of ghosts of '"
                    << parent_variable_name
                    << "' exceeds"
                    << " number of ghosts of '"
                    << variable_name
                    << "'."
                    << std::endl);
                */
                
                d_num_subghosts_pressure = num_subghosts;
            }
        }
        else
        {
            d_num_subghosts_pressure = num_subghosts;
        }
        
        setNumberOfSubGhosts(num_subghosts, "INTERNAL_ENERGY", parent_variable_name);
    }
    else if (variable_name == "SOUND_SPEED")
    {
        if (d_num_subghosts_sound_speed > -hier::IntVector::getOne(d_dim))
        {
            if (num_subghosts > d_num_subghosts_sound_speed)
            {
                /*
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::setNumberOfSubGhosts()\n"
                    << "Number of ghosts of '"
                    << parent_variable_name
                    << "' exceeds"
                    << " number of ghosts of '"
                    << variable_name
                    << "'."
                    << std::endl);
                */
                
                d_num_subghosts_sound_speed = num_subghosts;
            }
        }
        else
        {
            d_num_subghosts_sound_speed = num_subghosts;
        }
        
        setNumberOfSubGhosts(num_subghosts, "PRESSURE", parent_variable_name);
    }
    else if (variable_name == "TEMPERATURE")
    {
        if (d_num_subghosts_temperature > -hier::IntVector::getOne(d_dim))
        {
            if (num_subghosts > d_num_subghosts_temperature)
            {
                /*
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::setNumberOfSubGhosts()\n"
                    << "Number of ghosts of '"
                    << parent_variable_name
                    << "' exceeds"
                    << " number of ghosts of '"
                    << variable_name
                    << "'."
                    << std::endl);
                */
                
                d_num_subghosts_temperature = num_subghosts;
            }
        }
        else
        {
            d_num_subghosts_temperature = num_subghosts;
        }
        
        setNumberOfSubGhosts(num_subghosts, "PRESSURE", parent_variable_name);
    }
    else if (variable_name == "CONVECTIVE_FLUX_X")
    {
        if (d_num_subghosts_convective_flux_x > -hier::IntVector::getOne(d_dim))
        {
            if (num_subghosts > d_num_subghosts_convective_flux_x)
            {
                /*
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::setNumberOfSubGhosts()\n"
                    << "Number of ghosts of '"
                    << parent_variable_name
                    << "' exceeds"
                    << " number of ghosts of '"
                    << variable_name
                    << "'."
                    << std::endl);
                */
                
                d_num_subghosts_convective_flux_x = num_subghosts;
            }
        }
        else
        {
            d_num_subghosts_convective_flux_x = num_subghosts;
        }
        
        setNumberOfSubGhosts(num_subghosts, "VELOCITY", parent_variable_name);
        setNumberOfSubGhosts(num_subghosts, "PRESSURE", parent_variable_name);
    }
    else if (variable_name == "CONVECTIVE_FLUX_Y")
    {
        if (d_num_subghosts_convective_flux_y > -hier::IntVector::getOne(d_dim))
        {
            if (num_subghosts > d_num_subghosts_convective_flux_y)
            {
                /*
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::setNumberOfSubGhosts()\n"
                    << "Number of ghosts of '"
                    << parent_variable_name
                    << "' exceeds"
                    << " number of ghosts of '"
                    << variable_name
                    << "'."
                    << std::endl);
                */
                
                d_num_subghosts_convective_flux_y = num_subghosts;
            }
        }
        else
        {
            d_num_subghosts_convective_flux_y = num_subghosts;
        }
        
        setNumberOfSubGhosts(num_subghosts, "VELOCITY", parent_variable_name);
        setNumberOfSubGhosts(num_subghosts, "PRESSURE", parent_variable_name);
    }
    else if (variable_name == "CONVECTIVE_FLUX_Z")
    {
        if (d_num_subghosts_convective_flux_z > -hier::IntVector::getOne(d_dim))
        {
            if (num_subghosts > d_num_subghosts_convective_flux_z)
            {
                /*
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::setNumberOfSubGhosts()\n"
                    << "Number of ghosts of '"
                    << parent_variable_name
                    << "' exceeds"
                    << " number of ghosts of '"
                    << variable_name
                    << "'."
                    << std::endl);
                */
                
                d_num_subghosts_convective_flux_z = num_subghosts;
            }
        }
        else
        {
            d_num_subghosts_convective_flux_z = num_subghosts;
        }
        
        setNumberOfSubGhosts(num_subghosts, "VELOCITY", parent_variable_name);
        setNumberOfSubGhosts(num_subghosts, "PRESSURE", parent_variable_name);
    }
    else if (variable_name == "PRIMITIVE_VARIABLES")
    {
        setNumberOfSubGhosts(num_subghosts, "VELOCITY", parent_variable_name);
        setNumberOfSubGhosts(num_subghosts, "PRESSURE", parent_variable_name);
    }
    else if (variable_name == "MAX_WAVE_SPEED_X")
    {
        if (d_num_subghosts_max_wave_speed_x > -hier::IntVector::getOne(d_dim))
        {
            if (num_subghosts > d_num_subghosts_max_wave_speed_x)
            {
                /*
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::setNumberOfSubGhosts()\n"
                    << "Number of ghosts of '"
                    << parent_variable_name
                    << "' exceeds"
                    << " number of ghosts of '"
                    << variable_name
                    << "'."
                    << std::endl);
                */
                
                d_num_subghosts_max_wave_speed_x = num_subghosts;
            }
        }
        else
        {
            d_num_subghosts_max_wave_speed_x = num_subghosts;
        }
        
        setNumberOfSubGhosts(num_subghosts, "VELOCITY", parent_variable_name);
        setNumberOfSubGhosts(num_subghosts, "SOUND_SPEED", parent_variable_name);
    }
    else if (variable_name == "MAX_WAVE_SPEED_Y")
    {
        if (d_num_subghosts_max_wave_speed_y > -hier::IntVector::getOne(d_dim))
        {
            if (num_subghosts > d_num_subghosts_max_wave_speed_y)
            {
                /*
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::setNumberOfSubGhosts()\n"
                    << "Number of ghosts of '"
                    << parent_variable_name
                    << "' exceeds"
                    << " number of ghosts of '"
                    << variable_name
                    << "'."
                    << std::endl);
                */
                
                d_num_subghosts_max_wave_speed_y = num_subghosts;
            }
        }
        else
        {
            d_num_subghosts_max_wave_speed_y = num_subghosts;
        }
        
        setNumberOfSubGhosts(num_subghosts, "VELOCITY", parent_variable_name);
        setNumberOfSubGhosts(num_subghosts, "SOUND_SPEED", parent_variable_name);
    }
    else if (variable_name == "MAX_WAVE_SPEED_Z")
    {
        if (d_num_subghosts_max_wave_speed_z > -hier::IntVector::getOne(d_dim))
        {
            if (num_subghosts > d_num_subghosts_max_wave_speed_z)
            {
                /*
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::setNumberOfSubGhosts()\n"
                    << "Number of ghosts of '"
                    << parent_variable_name
                    << "' exceeds"
                    << " number of ghosts of '"
                    << variable_name
                    << "'."
                    << std::endl);
                */
                
                d_num_subghosts_max_wave_speed_z = num_subghosts;
            }
        }
        else
        {
            d_num_subghosts_max_wave_speed_z = num_subghosts;
        }
        
        setNumberOfSubGhosts(num_subghosts, "VELOCITY", parent_variable_name);
        setNumberOfSubGhosts(num_subghosts, "SOUND_SPEED", parent_variable_name);
    }
    else if (variable_name == "MAX_DIFFUSIVITY")
    {
        if (d_num_subghosts_max_diffusivity > -hier::IntVector::getOne(d_dim))
        {
            if (num_subghosts > d_num_subghosts_max_diffusivity)
            {
                /*
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::setNumberOfSubGhosts()\n"
                    << "Number of ghosts of '"
                    << parent_variable_name
                    << "' exceeds"
                    << " number of ghosts of '"
                    << variable_name
                    << "'."
                    << std::endl);
                */
                
                d_num_subghosts_max_diffusivity = num_subghosts;
            }
        }
        else
        {
            d_num_subghosts_max_diffusivity = num_subghosts;
        }
        
        setNumberOfSubGhosts(num_subghosts, "PRESSURE", parent_variable_name);
        setNumberOfSubGhosts(num_subghosts, "TEMPERATURE", parent_variable_name);
    }
}


/*
 * Set the ghost boxes of derived cell variables.
 */
void
FlowModelSingleSpecies::setDerivedCellVariableGhostBoxes()
{
    if (d_num_subghosts_velocity > -hier::IntVector::getOne(d_dim))
    {
        d_subghost_box_velocity = d_interior_box;
        d_subghost_box_velocity.grow(d_num_subghosts_velocity);
        d_subghostcell_dims_velocity = d_subghost_box_velocity.numberCells();
    }
    
    if (d_num_subghosts_internal_energy > -hier::IntVector::getOne(d_dim))
    {
        d_subghost_box_internal_energy = d_interior_box;
        d_subghost_box_internal_energy.grow(d_num_subghosts_internal_energy);
        d_subghostcell_dims_internal_energy = d_subghost_box_internal_energy.numberCells();
    }
    
    if (d_num_subghosts_pressure > -hier::IntVector::getOne(d_dim))
    {
        d_subghost_box_pressure = d_interior_box;
        d_subghost_box_pressure.grow(d_num_subghosts_pressure);
        d_subghostcell_dims_pressure = d_subghost_box_pressure.numberCells();
    }
    
    if (d_num_subghosts_sound_speed > -hier::IntVector::getOne(d_dim))
    {
        d_subghost_box_sound_speed = d_interior_box;
        d_subghost_box_sound_speed.grow(d_num_subghosts_sound_speed);
        d_subghostcell_dims_sound_speed = d_subghost_box_sound_speed.numberCells();
    }
    
    if (d_num_subghosts_temperature > -hier::IntVector::getOne(d_dim))
    {
        d_subghost_box_temperature = d_interior_box;
        d_subghost_box_temperature.grow(d_num_subghosts_temperature);
        d_subghostcell_dims_temperature = d_subghost_box_temperature.numberCells();
    }
    
    if (d_num_subghosts_convective_flux_x > -hier::IntVector::getOne(d_dim))
    {
        d_subghost_box_convective_flux_x = d_interior_box;
        d_subghost_box_convective_flux_x.grow(d_num_subghosts_convective_flux_x);
        d_subghostcell_dims_convective_flux_x = d_subghost_box_convective_flux_x.numberCells();
    }
    
    if (d_num_subghosts_convective_flux_y > -hier::IntVector::getOne(d_dim))
    {
        d_subghost_box_convective_flux_y = d_interior_box;
        d_subghost_box_convective_flux_y.grow(d_num_subghosts_convective_flux_y);
        d_subghostcell_dims_convective_flux_y = d_subghost_box_convective_flux_y.numberCells();
    }
    
    if (d_num_subghosts_convective_flux_z > -hier::IntVector::getOne(d_dim))
    {
        d_subghost_box_convective_flux_z = d_interior_box;
        d_subghost_box_convective_flux_z.grow(d_num_subghosts_convective_flux_z);
        d_subghostcell_dims_convective_flux_z = d_subghost_box_convective_flux_z.numberCells();
    }
    
    if (d_num_subghosts_max_wave_speed_x > -hier::IntVector::getOne(d_dim))
    {
        d_subghost_box_max_wave_speed_x = d_interior_box;
        d_subghost_box_max_wave_speed_x.grow(d_num_subghosts_max_wave_speed_x);
        d_subghostcell_dims_max_wave_speed_x = d_subghost_box_max_wave_speed_x.numberCells();
    }
    
    if (d_num_subghosts_max_wave_speed_y > -hier::IntVector::getOne(d_dim))
    {
        d_subghost_box_max_wave_speed_y = d_interior_box;
        d_subghost_box_max_wave_speed_y.grow(d_num_subghosts_max_wave_speed_y);
        d_subghostcell_dims_max_wave_speed_y = d_subghost_box_max_wave_speed_y.numberCells();
    }
    
    if (d_num_subghosts_max_wave_speed_z > -hier::IntVector::getOne(d_dim))
    {
        d_subghost_box_max_wave_speed_z = d_interior_box;
        d_subghost_box_max_wave_speed_z.grow(d_num_subghosts_max_wave_speed_z);
        d_subghostcell_dims_max_wave_speed_z = d_subghost_box_max_wave_speed_z.numberCells();
    }
    
    if (d_num_subghosts_max_diffusivity > -hier::IntVector::getOne(d_dim))
    {
        d_subghost_box_max_diffusivity = d_interior_box;
        d_subghost_box_max_diffusivity.grow(d_num_subghosts_max_diffusivity);
        d_subghostcell_dims_max_diffusivity = d_subghost_box_max_diffusivity.numberCells();
    }
}


/*
 * Get the cell data of density in the registered patch.
 */
boost::shared_ptr<pdat::CellData<double> >
FlowModelSingleSpecies::getCellDataOfDensity()
{
    // Get the cell data of the registered variable density.
    boost::shared_ptr<pdat::CellData<double> > data_density(
        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
            d_patch->getPatchData(s_variable_density, getDataContext())));
    
    return data_density;
}


/*
 * Get the cell data of momentum in the registered patch.
 */
boost::shared_ptr<pdat::CellData<double> >
FlowModelSingleSpecies::getCellDataOfMomentum()
{
    // Get the cell data of the registered variable momentum.
    boost::shared_ptr<pdat::CellData<double> > data_momentum(
        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
            d_patch->getPatchData(s_variable_momentum, getDataContext())));
    
    return data_momentum;
}


/*
 * Get the cell data of total energy in the registered patch.
 */
boost::shared_ptr<pdat::CellData<double> >
FlowModelSingleSpecies::getCellDataOfTotalEnergy()
{
    // Get the cell data of the registered variable total energy.
    boost::shared_ptr<pdat::CellData<double> > data_total_energy(
        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
            d_patch->getPatchData(s_variable_total_energy, getDataContext())));
    
    return data_total_energy;
}


/*
 * Compute the cell data of velocity in the registered patch.
 */
void
FlowModelSingleSpecies::computeCellDataOfVelocity(
    const hier::Box& domain)
{
    if (d_num_subghosts_velocity > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_velocity_computed)
        {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(d_data_velocity);
#endif
            /*
             * Get the local lower index and number of cells in each direction of the domain.
             */
            
            hier::IntVector domain_lo(d_dim);
            hier::IntVector domain_dims(d_dim);
            
            if (domain.empty())
            {
                domain_lo = -d_num_subghosts_velocity;
                domain_dims = d_subghostcell_dims_velocity;
            }
            else
            {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(d_subghost_box_velocity.contains(domain));
#endif
                
                domain_lo = domain.lower() - d_interior_box.lower();
                domain_dims = domain.numberCells();
            }
            
            // Get the cell data of the variables density and momentum.
            boost::shared_ptr<pdat::CellData<double> > data_density =
                getCellDataOfDensity();
            
            boost::shared_ptr<pdat::CellData<double> > data_momentum =
                getCellDataOfMomentum();
            
            // Get the pointer to the cell data of density.
            double* rho = data_density->getPointer(0);
            
            if (d_dim == tbox::Dimension(1))
            {
                /*
                 * Get the local lower index, numbers of cells in each dimension and numbers of ghost cells.
                 */
                
                const int domain_lo_0 = domain_lo[0];
                const int domain_dim_0 = domain_dims[0];
                
                const int num_ghosts_0 = d_num_ghosts[0];
                const int num_subghosts_0_velocity = d_num_subghosts_velocity[0];
                
                // Get the pointer to the cell data of velocity.
                double* u = d_data_velocity->getPointer(0);
                
                // Get the pointer to the cell data of momentum.
                double* rho_u = data_momentum->getPointer(0);
                
                // Compute the velocity field.
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx = i + num_ghosts_0;
                    const int idx_velocity = i + num_subghosts_0_velocity;
                    
                    u[idx_velocity] = rho_u[idx]/rho[idx];
                }
            }
            else if (d_dim == tbox::Dimension(2))
            {
                /*
                 * Get the local lower indices, numbers of cells in each dimension and numbers of ghost cells.
                 */
                
                const int domain_lo_0 = domain_lo[0];
                const int domain_lo_1 = domain_lo[1];
                const int domain_dim_0 = domain_dims[0];
                const int domain_dim_1 = domain_dims[1];
                
                const int num_ghosts_0 = d_num_ghosts[0];
                const int num_ghosts_1 = d_num_ghosts[1];
                const int ghostcell_dim_0 = d_ghostcell_dims[0];
                
                const int num_subghosts_0_velocity = d_num_subghosts_velocity[0];
                const int num_subghosts_1_velocity = d_num_subghosts_velocity[1];
                const int subghostcell_dim_0_velocity = d_subghostcell_dims_velocity[0];
                
                // Get the pointers to the cell data of velocity.
                double* u = d_data_velocity->getPointer(0);
                double* v = d_data_velocity->getPointer(1);
                
                // Get the pointers to the cell data of momentum.
                double* rho_u = data_momentum->getPointer(0);
                double* rho_v = data_momentum->getPointer(1);
                
                // Compute the velocity field.
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx = (i + num_ghosts_0) +
                            (j + num_ghosts_1)*ghostcell_dim_0;
                        
                        const int idx_velocity = (i + num_subghosts_0_velocity) +
                            (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                        
                        u[idx_velocity] = rho_u[idx]/rho[idx];
                        v[idx_velocity] = rho_v[idx]/rho[idx];
                    }
                }
            }
            else if (d_dim == tbox::Dimension(3))
            {
                /*
                 * Get the local lower indices, numbers of cells in each dimension and numbers of ghost cells.
                 */
                
                const int domain_lo_0 = domain_lo[0];
                const int domain_lo_1 = domain_lo[1];
                const int domain_lo_2 = domain_lo[2];
                const int domain_dim_0 = domain_dims[0];
                const int domain_dim_1 = domain_dims[1];
                const int domain_dim_2 = domain_dims[2];
                
                const int num_ghosts_0 = d_num_ghosts[0];
                const int num_ghosts_1 = d_num_ghosts[1];
                const int num_ghosts_2 = d_num_ghosts[2];
                const int ghostcell_dim_0 = d_ghostcell_dims[0];
                const int ghostcell_dim_1 = d_ghostcell_dims[1];
                
                const int num_subghosts_0_velocity = d_num_subghosts_velocity[0];
                const int num_subghosts_1_velocity = d_num_subghosts_velocity[1];
                const int num_subghosts_2_velocity = d_num_subghosts_velocity[2];
                const int subghostcell_dim_0_velocity = d_subghostcell_dims_velocity[0];
                const int subghostcell_dim_1_velocity = d_subghostcell_dims_velocity[1];
                
                // Get the pointers to the cell data of velocity.
                double* u = d_data_velocity->getPointer(0);
                double* v = d_data_velocity->getPointer(1);
                double* w = d_data_velocity->getPointer(2);
                
                // Get the pointers to the cell data of momentum.
                double* rho_u = data_momentum->getPointer(0);
                double* rho_v = data_momentum->getPointer(1);
                double* rho_w = data_momentum->getPointer(2);
                
                // Compute the velocity field.
                for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                {
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx = (i + num_ghosts_0) +
                                (j + num_ghosts_1)*ghostcell_dim_0 +
                                (k + num_ghosts_2)*ghostcell_dim_0*ghostcell_dim_1;
                            
                            const int idx_velocity = (i + num_subghosts_0_velocity) +
                                (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                    subghostcell_dim_1_velocity;
                            
                            u[idx_velocity] = rho_u[idx]/rho[idx];
                            v[idx_velocity] = rho_v[idx]/rho[idx];
                            w[idx_velocity] = rho_w[idx]/rho[idx];
                        }
                    }
                }
            }
            
            d_cell_data_velocity_computed = true;
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::computeCellDataOfVelocity()\n"
            << "Cell data of 'VELOCITY' is not yet registered."
            << std::endl);
    }
}


/*
 * Compute the cell data of internal energy with velocity in the registered patch.
 */
void
FlowModelSingleSpecies::computeCellDataOfInternalEnergyWithVelocity(
    const hier::Box& domain)
{
    if (d_num_subghosts_internal_energy > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_internal_energy_computed)
        {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(d_data_internal_energy);
#endif
            
            /*
             * Get the local lower index and number of cells in each direction of the domain.
             */
            
            hier::IntVector domain_lo(d_dim);
            hier::IntVector domain_dims(d_dim);
            
            if (domain.empty())
            {
                domain_lo = -d_num_subghosts_internal_energy;
                domain_dims = d_subghostcell_dims_internal_energy;
            }
            else
            {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(d_subghost_box_internal_energy.contains(domain));
#endif
                
                domain_lo = domain.lower() - d_interior_box.lower();
                domain_dims = domain.numberCells();
            }
            
            // Get the cell data of the variables density, momentum and total energy.
            boost::shared_ptr<pdat::CellData<double> > data_density =
                getCellDataOfDensity();
            
            boost::shared_ptr<pdat::CellData<double> > data_total_energy =
                getCellDataOfTotalEnergy();
            
            if (!d_cell_data_velocity_computed)
            {
                computeCellDataOfVelocity(domain);
            }
            
            // Get the pointers to the cell data of internal energy, density and total energy.
            double* epsilon = d_data_internal_energy->getPointer(0);
            double* rho = data_density->getPointer(0);
            double* E   = data_total_energy->getPointer(0);
            
            // Get the thermodynamic properties of the species.
            std::vector<const double*> thermo_properties_ptr;
            thermo_properties_ptr.reserve(static_cast<int> (d_thermo_properties.size()));
            for (int ti = 0; ti < static_cast<int> (d_thermo_properties.size()); ti++)
            {
                thermo_properties_ptr.push_back(&d_thermo_properties[ti]);
            }
            
            if (d_dim == tbox::Dimension(1))
            {
                /*
                 * Get the local lower index, numbers of cells in each dimension and numbers of ghost cells.
                 */
                
                const int domain_lo_0 = domain_lo[0];
                const int domain_dim_0 = domain_dims[0];
                
                const int num_ghosts_0 = d_num_ghosts[0];
                const int num_subghosts_0_velocity = d_num_subghosts_velocity[0];
                const int num_subghosts_0_internal_energy = d_num_subghosts_internal_energy[0];
                
                // Get the pointer to cell data of velocity.
                double* u = d_data_velocity->getPointer(0);
                
                // Compute the internal energy field.
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx = i + num_ghosts_0;
                    const int idx_velocity = i + num_subghosts_0_velocity;
                    const int idx_internal_energy = i + num_subghosts_0_internal_energy;
                    
                    epsilon[idx_internal_energy] = E[idx]/rho[idx] -
                        double(1)/double(2)*u[idx_velocity]*u[idx_velocity];
                }
            }
            else if (d_dim == tbox::Dimension(2))
            {
                /*
                 * Get the local lower indices, numbers of cells in each dimension and numbers of ghost cells.
                 */
                
                const int domain_lo_0 = domain_lo[0];
                const int domain_lo_1 = domain_lo[1];
                const int domain_dim_0 = domain_dims[0];
                const int domain_dim_1 = domain_dims[1];
                
                const int num_ghosts_0 = d_num_ghosts[0];
                const int num_ghosts_1 = d_num_ghosts[1];
                const int ghostcell_dim_0 = d_ghostcell_dims[0];
                
                const int num_subghosts_0_velocity = d_num_subghosts_velocity[0];
                const int num_subghosts_1_velocity = d_num_subghosts_velocity[1];
                const int subghostcell_dim_0_velocity = d_subghostcell_dims_velocity[0];
                
                const int num_subghosts_0_internal_energy = d_num_subghosts_internal_energy[0];
                const int num_subghosts_1_internal_energy = d_num_subghosts_internal_energy[1];
                const int subghostcell_dim_0_internal_energy = d_subghostcell_dims_internal_energy[0];
                
                // Get the pointers to the cell data of velocity.
                double* u = d_data_velocity->getPointer(0);
                double* v = d_data_velocity->getPointer(1);
                
                // Compute the internal energy field.
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx = (i + num_ghosts_0) +
                            (j + num_ghosts_1)*ghostcell_dim_0;
                        
                        const int idx_velocity = (i + num_subghosts_0_velocity) +
                            (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                        
                        const int idx_internal_energy = (i + num_subghosts_0_internal_energy) +
                            (j + num_subghosts_1_internal_energy)*subghostcell_dim_0_internal_energy;
                        
                        epsilon[idx_internal_energy] = E[idx]/rho[idx] -
                            double(1)/double(2)*(u[idx_velocity]*u[idx_velocity] + v[idx_velocity]*v[idx_velocity]);
                    }
                }
            }
            else if (d_dim == tbox::Dimension(3))
            {
                /*
                 * Get the local lower indices, numbers of cells in each dimension and numbers of ghost cells.
                 */
                
                const int domain_lo_0 = domain_lo[0];
                const int domain_lo_1 = domain_lo[1];
                const int domain_lo_2 = domain_lo[2];
                const int domain_dim_0 = domain_dims[0];
                const int domain_dim_1 = domain_dims[1];
                const int domain_dim_2 = domain_dims[2];
                
                const int num_ghosts_0 = d_num_ghosts[0];
                const int num_ghosts_1 = d_num_ghosts[1];
                const int num_ghosts_2 = d_num_ghosts[2];
                const int ghostcell_dim_0 = d_ghostcell_dims[0];
                const int ghostcell_dim_1 = d_ghostcell_dims[1];
                
                const int num_subghosts_0_velocity = d_num_subghosts_velocity[0];
                const int num_subghosts_1_velocity = d_num_subghosts_velocity[1];
                const int num_subghosts_2_velocity = d_num_subghosts_velocity[2];
                const int subghostcell_dim_0_velocity = d_subghostcell_dims_velocity[0];
                const int subghostcell_dim_1_velocity = d_subghostcell_dims_velocity[1];
                
                const int num_subghosts_0_internal_energy = d_num_subghosts_internal_energy[0];
                const int num_subghosts_1_internal_energy = d_num_subghosts_internal_energy[1];
                const int num_subghosts_2_internal_energy = d_num_subghosts_internal_energy[2];
                const int subghostcell_dim_0_internal_energy = d_subghostcell_dims_internal_energy[0];
                const int subghostcell_dim_1_internal_energy = d_subghostcell_dims_internal_energy[1];
                
                // Get the pointers to the cell data of velocity.
                double* u = d_data_velocity->getPointer(0);
                double* v = d_data_velocity->getPointer(1);
                double* w = d_data_velocity->getPointer(2);
                
                // Compute the internal energy field.
                for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                {
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx = (i + num_ghosts_0) +
                                (j + num_ghosts_1)*ghostcell_dim_0 +
                                (k + num_ghosts_2)*ghostcell_dim_0*ghostcell_dim_1;
                            
                            const int idx_velocity = (i + num_subghosts_0_velocity) +
                                (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                    subghostcell_dim_1_velocity;
                            
                            const int idx_internal_energy = (i + num_subghosts_0_internal_energy) +
                                (j + num_subghosts_1_internal_energy)*subghostcell_dim_0_internal_energy +
                                (k + num_subghosts_2_internal_energy)*subghostcell_dim_0_internal_energy*
                                    subghostcell_dim_1_internal_energy;
                            
                            epsilon[idx_internal_energy] = E[idx]/rho[idx] -
                                double(1)/double(2)*(u[idx_velocity]*u[idx_velocity] + v[idx_velocity]*v[idx_velocity] +
                                    w[idx_velocity]*w[idx_velocity]);
                        }
                    }
                }
            }
            
            d_cell_data_internal_energy_computed = true;
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::computeCellDataOfInternalEnergyWithVelocity()\n"
            << "Cell data of 'INTERNAL_ENERGY' is not yet registered."
            << std::endl);
    }
}


/*
 * Compute the cell data of pressure with internal energy in the registered patch.
 */
void
FlowModelSingleSpecies::computeCellDataOfPressureWithInternalEnergy(
    const hier::Box& domain)
{
    if (d_num_subghosts_pressure > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_pressure_computed)
        {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(d_data_pressure);
#endif
            
            // Get the cell data of the variables density, momentum and total energy.
            boost::shared_ptr<pdat::CellData<double> > data_density =
                getCellDataOfDensity();
            
            if (!d_cell_data_internal_energy_computed)
            {
                computeCellDataOfInternalEnergyWithVelocity(domain);
            }
            
            // Get the thermodynamic properties of the species.
            std::vector<const double*> thermo_properties_ptr;
            thermo_properties_ptr.reserve(static_cast<int> (d_thermo_properties.size()));
            for (int ti = 0; ti < static_cast<int> (d_thermo_properties.size()); ti++)
            {
                thermo_properties_ptr.push_back(&d_thermo_properties[ti]);
            }
            
            // Compute the pressure field.
            d_equation_of_state_mixing_rules->getEquationOfState()->computePressure(
                d_data_pressure,
                data_density,
                d_data_internal_energy,
                thermo_properties_ptr,
                domain);
            
            d_cell_data_pressure_computed = true;
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::computeCellDataOfPressureWithInternalEnergy()\n"
            << "Cell data of 'PRESSURE' is not yet registered."
            << std::endl);
    }
}


/*
 * Compute the cell data of sound speed with pressure in the registered patch.
 */
void
FlowModelSingleSpecies::computeCellDataOfSoundSpeedWithPressure(
    const hier::Box& domain)
{
    if (d_num_subghosts_sound_speed > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_sound_speed_computed)
        {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(d_data_sound_speed);
#endif
            
            // Get the cell data of the variable density and pressure.
            boost::shared_ptr<pdat::CellData<double> > data_density =
                getCellDataOfDensity();
            
            if (!d_cell_data_pressure_computed)
            {
                computeCellDataOfPressureWithInternalEnergy(domain);
            }
            
            // Get the thermodynamic properties of the species.
            std::vector<const double*> thermo_properties_ptr;
            thermo_properties_ptr.reserve(static_cast<int> (d_thermo_properties.size()));
            for (int ti = 0; ti < static_cast<int> (d_thermo_properties.size()); ti++)
            {
                thermo_properties_ptr.push_back(&d_thermo_properties[ti]);
            }
            
            // Compute the sound speed field.
            d_equation_of_state_mixing_rules->getEquationOfState()->computeSoundSpeed(
                d_data_sound_speed,
                data_density,
                d_data_pressure,
                thermo_properties_ptr,
                domain);
            
            d_cell_data_sound_speed_computed = true;
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::computeCellDataOfSoundSpeedWithPressure()\n"
            << "Cell data of 'SOUND_SPEED' is not yet registered."
            << std::endl);
    }
}


/*
 * Compute the cell data of temperature with pressure in the registered patch.
 */
void
FlowModelSingleSpecies::computeCellDataOfTemperatureWithPressure(
    const hier::Box& domain)
{
    if (d_num_subghosts_temperature > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_temperature_computed)
        {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(d_data_temperature);
#endif
            
            // Get the cell data of the variable density and pressure.
            boost::shared_ptr<pdat::CellData<double> > data_density =
                getCellDataOfDensity();
            
            if (!d_cell_data_pressure_computed)
            {
                computeCellDataOfPressureWithInternalEnergy(domain);
            }
            
            // Get the thermodynamic properties of the species.
            std::vector<const double*> thermo_properties_ptr;
            thermo_properties_ptr.reserve(static_cast<int> (d_thermo_properties.size()));
            for (int ti = 0; ti < static_cast<int> (d_thermo_properties.size()); ti++)
            {
                thermo_properties_ptr.push_back(&d_thermo_properties[ti]);
            }
            
            // Compute the temperature field.
            d_equation_of_state_mixing_rules->getEquationOfState()->computeTemperature(
                d_data_temperature,
                data_density,
                d_data_pressure,
                thermo_properties_ptr,
                domain);
            
            d_cell_data_temperature_computed = true;
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::computeCellDataOfTemperatureWithPressure()\n"
            << "Cell data of 'TEMPERATURE' is not yet registered."
            << std::endl);
    }
}


/*
 * Compute the cell data of convective flux with velocity and pressure in the registered patch.
 */
void
FlowModelSingleSpecies::computeCellDataOfConvectiveFluxWithVelocityAndPressure(
    const DIRECTION::TYPE& direction,
    const hier::Box& domain)
{
    if (direction == DIRECTION::X_DIRECTION)
    {
        if (d_num_subghosts_convective_flux_x > -hier::IntVector::getOne(d_dim))
        {
            if (!d_cell_data_convective_flux_x_computed)
            {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(d_data_convective_flux_x);
#endif
                
                /*
                 * Get the local lower index and number of cells in each direction of the domain.
                 */
                
                hier::IntVector domain_lo(d_dim);
                hier::IntVector domain_dims(d_dim);
                
                if (domain.empty())
                {
                    domain_lo = -d_num_subghosts_convective_flux_x;
                    domain_dims = d_subghostcell_dims_convective_flux_x;
                }
                else
                {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
                    TBOX_ASSERT(d_subghost_box_convective_flux_x.contains(domain));
#endif
                    
                    domain_lo = domain.lower() - d_interior_box.lower();
                    domain_dims = domain.numberCells();
                }
                
                // Get the pointers to the components of the convective flux in the x-direction.
                std::vector<double*> F_x;
                F_x.reserve(d_num_eqn);
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    F_x.push_back(d_data_convective_flux_x->getPointer(ei));
                }
                
                boost::shared_ptr<pdat::CellData<double> > data_momentum =
                    getCellDataOfMomentum();
                
                boost::shared_ptr<pdat::CellData<double> > data_total_energy =
                    getCellDataOfTotalEnergy();
                
                if (!d_cell_data_velocity_computed)
                {
                    computeCellDataOfVelocity(domain);
                }
                
                if (!d_cell_data_pressure_computed)
                {
                    computeCellDataOfPressureWithInternalEnergy(domain);
                }
                
                // Get the pointers to the cell data of total energy and pressure.
                double* E = data_total_energy->getPointer(0);
                double* p = d_data_pressure->getPointer(0);
                
                if (d_dim == tbox::Dimension(1))
                {
                    /*
                     * Get the local lower index, numbers of cells in each dimension and numbers of ghost cells.
                     */
                    
                    const int domain_lo_0 = domain_lo[0];
                    const int domain_dim_0 = domain_dims[0];
                    
                    const int num_ghosts_0 = d_num_ghosts[0];
                    const int num_subghosts_0_pressure = d_num_subghosts_pressure[0];
                    const int num_subghosts_0_velocity = d_num_subghosts_velocity[0];
                    const int num_subghosts_0_convective_flux_x = d_num_subghosts_convective_flux_x[0];
                    
                    // Get the pointer to the cell data of momentum.
                    double* rho_u = data_momentum->getPointer(0);
                    
                    // Get the pointer to the cell data of velocity.
                    double* u = d_data_velocity->getPointer(0);
                    
                    // Compute the convective flux in the x-direction.
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx = i + num_ghosts_0;
                        const int idx_pressure = i + num_subghosts_0_pressure;
                        const int idx_velocity = i + num_subghosts_0_velocity;
                        const int idx_convective_flux_x = i + num_subghosts_0_convective_flux_x;
                        
                        F_x[0][idx_convective_flux_x] = rho_u[idx];
                        F_x[1][idx_convective_flux_x] = u[idx_velocity]*rho_u[idx] + p[idx_pressure];
                        F_x[2][idx_convective_flux_x] = u[idx_velocity]*(E[idx] + p[idx_pressure]);
                    }
                }
                else if (d_dim == tbox::Dimension(2))
                {
                    /*
                     * Get the local lower indices, numbers of cells in each dimension and numbers of ghost cells.
                     */
                    
                    const int domain_lo_0 = domain_lo[0];
                    const int domain_lo_1 = domain_lo[1];
                    const int domain_dim_0 = domain_dims[0];
                    const int domain_dim_1 = domain_dims[1];
                    
                    const int num_ghosts_0 = d_num_ghosts[0];
                    const int num_ghosts_1 = d_num_ghosts[1];
                    const int ghostcell_dim_0 = d_ghostcell_dims[0];
                    
                    const int num_subghosts_0_pressure = d_num_subghosts_pressure[0];
                    const int num_subghosts_1_pressure = d_num_subghosts_pressure[1];
                    const int subghostcell_dim_0_pressure = d_subghostcell_dims_pressure[0];
                    
                    const int num_subghosts_0_velocity = d_num_subghosts_velocity[0];
                    const int num_subghosts_1_velocity = d_num_subghosts_velocity[1];
                    const int subghostcell_dim_0_velocity = d_subghostcell_dims_velocity[0];
                    
                    const int num_subghosts_0_convective_flux_x = d_num_subghosts_convective_flux_x[0];
                    const int num_subghosts_1_convective_flux_x = d_num_subghosts_convective_flux_x[1];
                    const int subghostcell_dim_0_convective_flux_x = d_subghostcell_dims_convective_flux_x[0];
                    
                    // Get the pointers to the cell data of momentum.
                    double* rho_u = data_momentum->getPointer(0);
                    double* rho_v = data_momentum->getPointer(1);
                    
                    // Get the pointer to the cell data of velocity.
                    double* u = d_data_velocity->getPointer(0);
                    
                    // Compute the convective flux in the x-direction.
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            const int idx = (i + num_ghosts_0) +
                                (j + num_ghosts_1)*ghostcell_dim_0;
                            
                            const int idx_pressure = (i + num_subghosts_0_pressure) +
                                (j + num_subghosts_1_pressure)*subghostcell_dim_0_pressure;
                            
                            const int idx_velocity = (i + num_subghosts_0_velocity) +
                                (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                            
                            const int idx_convective_flux_x = (i + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                            
                            F_x[0][idx_convective_flux_x] = rho_u[idx];
                            F_x[1][idx_convective_flux_x] = u[idx_velocity]*rho_u[idx] + p[idx_pressure];
                            F_x[2][idx_convective_flux_x] = u[idx_velocity]*rho_v[idx];
                            F_x[3][idx_convective_flux_x] = u[idx_velocity]*(E[idx] + p[idx_pressure]);
                        }
                    }
                }
                else if (d_dim == tbox::Dimension(3))
                {
                    /*
                     * Get the local lower indices, numbers of cells in each dimension and numbers of ghost cells.
                     */
                    
                    const int domain_lo_0 = domain_lo[0];
                    const int domain_lo_1 = domain_lo[1];
                    const int domain_lo_2 = domain_lo[2];
                    const int domain_dim_0 = domain_dims[0];
                    const int domain_dim_1 = domain_dims[1];
                    const int domain_dim_2 = domain_dims[2];
                    
                    const int num_ghosts_0 = d_num_ghosts[0];
                    const int num_ghosts_1 = d_num_ghosts[1];
                    const int num_ghosts_2 = d_num_ghosts[2];
                    const int ghostcell_dim_0 = d_ghostcell_dims[0];
                    const int ghostcell_dim_1 = d_ghostcell_dims[1];
                    
                    const int num_subghosts_0_pressure = d_num_subghosts_pressure[0];
                    const int num_subghosts_1_pressure = d_num_subghosts_pressure[1];
                    const int num_subghosts_2_pressure = d_num_subghosts_pressure[2];
                    const int subghostcell_dim_0_pressure = d_subghostcell_dims_pressure[0];
                    const int subghostcell_dim_1_pressure = d_subghostcell_dims_pressure[1];
                    
                    const int num_subghosts_0_velocity = d_num_subghosts_velocity[0];
                    const int num_subghosts_1_velocity = d_num_subghosts_velocity[1];
                    const int num_subghosts_2_velocity = d_num_subghosts_velocity[2];
                    const int subghostcell_dim_0_velocity = d_subghostcell_dims_velocity[0];
                    const int subghostcell_dim_1_velocity = d_subghostcell_dims_velocity[1];
                    
                    const int num_subghosts_0_convective_flux_x = d_num_subghosts_convective_flux_x[0];
                    const int num_subghosts_1_convective_flux_x = d_num_subghosts_convective_flux_x[1];
                    const int num_subghosts_2_convective_flux_x = d_num_subghosts_convective_flux_x[2];
                    const int subghostcell_dim_0_convective_flux_x = d_subghostcell_dims_convective_flux_x[0];
                    const int subghostcell_dim_1_convective_flux_x = d_subghostcell_dims_convective_flux_x[1];
                    
                    // Get the pointers to the cell data of momentum.
                    double* rho_u = data_momentum->getPointer(0);
                    double* rho_v = data_momentum->getPointer(1);
                    double* rho_w = data_momentum->getPointer(2);
                    
                    // Get the pointer to the cell data of velocity.
                    double* u = d_data_velocity->getPointer(0);
                    
                    // Compute the convective flux in the x-direction.
                    for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                    {
                        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                            {
                                // Compute the linear indices.
                                const int idx = (i + num_ghosts_0) +
                                    (j + num_ghosts_1)*ghostcell_dim_0 +
                                    (k + num_ghosts_2)*ghostcell_dim_0*ghostcell_dim_1;
                                
                                const int idx_pressure = (i + num_subghosts_0_pressure) +
                                    (j + num_subghosts_1_pressure)*subghostcell_dim_0_pressure +
                                    (k + num_subghosts_2_pressure)*subghostcell_dim_0_pressure*
                                        subghostcell_dim_1_pressure;
                                
                                const int idx_velocity = (i + num_subghosts_0_velocity) +
                                    (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                    (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                        subghostcell_dim_1_velocity;
                                
                                const int idx_convective_flux_x = (i + num_subghosts_0_convective_flux_x) +
                                    (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                    (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                        subghostcell_dim_1_convective_flux_x;
                                
                                F_x[0][idx_convective_flux_x] = rho_u[idx];
                                F_x[1][idx_convective_flux_x] = u[idx_velocity]*rho_u[idx]+ p[idx_pressure];
                                F_x[2][idx_convective_flux_x] = u[idx_velocity]*rho_v[idx];
                                F_x[3][idx_convective_flux_x] = u[idx_velocity]*rho_w[idx];
                                F_x[4][idx_convective_flux_x] = u[idx_velocity]*(E[idx] + p[idx_pressure]);
                            }
                        }
                    }
                }
                
                d_cell_data_convective_flux_x_computed = true;
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::computeCellDataOfConvectiveFluxWithVelocityAndPressure()\n"
                << "Cell data of 'CONVECTIVE_FLUX_X' is not yet registered."
                << std::endl);
        }
    }
    else if (direction == DIRECTION::Y_DIRECTION)
    {
        if (d_num_subghosts_convective_flux_y > -hier::IntVector::getOne(d_dim))
        {
            if (!d_cell_data_convective_flux_y_computed)
            {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(d_data_convective_flux_y);
#endif
                /*
                 * Get the local lower index and number of cells in each direction of the domain.
                 */
                
                hier::IntVector domain_lo(d_dim);
                hier::IntVector domain_dims(d_dim);
                
                if (domain.empty())
                {
                    domain_lo = -d_num_subghosts_convective_flux_y;
                    domain_dims = d_subghostcell_dims_convective_flux_y;
                }
                else
                {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
                    TBOX_ASSERT(d_subghost_box_convective_flux_y.contains(domain));
#endif
                    
                    domain_lo = domain.lower() - d_interior_box.lower();
                    domain_dims = domain.numberCells();
                }
                
                // Get the pointers to the components of the convective flux in the y-direction.
                std::vector<double*> F_y;
                F_y.reserve(d_num_eqn);
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    F_y.push_back(d_data_convective_flux_y->getPointer(ei));
                }
                
                boost::shared_ptr<pdat::CellData<double> > data_momentum =
                    getCellDataOfMomentum();
                
                boost::shared_ptr<pdat::CellData<double> > data_total_energy =
                    getCellDataOfTotalEnergy();
                
                if (!d_cell_data_velocity_computed)
                {
                    computeCellDataOfVelocity(domain);
                }
                
                if (!d_cell_data_pressure_computed)
                {
                    computeCellDataOfPressureWithInternalEnergy(domain);
                }
                
                // Get the pointers to the cell data of total energy and pressure.
                double* E = data_total_energy->getPointer(0);
                double* p = d_data_pressure->getPointer(0);
                
                if (d_dim == tbox::Dimension(1))
                {
                    TBOX_ERROR(d_object_name
                        << ": FlowModelSingleSpecies::computeCellDataOfConvectiveFluxWithVelocityAndPressure()\n"
                        << "'CONVECTIVE_FLUX_Y' cannot be obtained for problem with dimension less than two."
                        << std::endl);
                }
                else if (d_dim == tbox::Dimension(2))
                {
                    /*
                     * Get the local lower indices, numbers of cells in each dimension and numbers of ghost cells.
                     */
                    
                    const int domain_lo_0 = domain_lo[0];
                    const int domain_lo_1 = domain_lo[1];
                    const int domain_dim_0 = domain_dims[0];
                    const int domain_dim_1 = domain_dims[1];
                    
                    const int num_ghosts_0 = d_num_ghosts[0];
                    const int num_ghosts_1 = d_num_ghosts[1];
                    const int ghostcell_dim_0 = d_ghostcell_dims[0];
                    
                    const int num_subghosts_0_pressure = d_num_subghosts_pressure[0];
                    const int num_subghosts_1_pressure = d_num_subghosts_pressure[1];
                    const int subghostcell_dim_0_pressure = d_subghostcell_dims_pressure[0];
                    
                    const int num_subghosts_0_velocity = d_num_subghosts_velocity[0];
                    const int num_subghosts_1_velocity = d_num_subghosts_velocity[1];
                    const int subghostcell_dim_0_velocity = d_subghostcell_dims_velocity[0];
                    
                    const int num_subghosts_0_convective_flux_y = d_num_subghosts_convective_flux_y[0];
                    const int num_subghosts_1_convective_flux_y = d_num_subghosts_convective_flux_y[1];
                    const int subghostcell_dim_0_convective_flux_y = d_subghostcell_dims_convective_flux_y[0];
                    
                    // Get the pointers to the cell data of momentum.
                    double* rho_u = data_momentum->getPointer(0);
                    double* rho_v = data_momentum->getPointer(1);
                    
                    // Get the pointer to the cell data of velocity.
                    double* v = d_data_velocity->getPointer(1);
                    
                    // Compute the convective flux in the y-direction.
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            const int idx = (i + num_ghosts_0) +
                                (j + num_ghosts_1)*ghostcell_dim_0;
                            
                            const int idx_pressure = (i + num_subghosts_0_pressure) +
                                (j + num_subghosts_1_pressure)*subghostcell_dim_0_pressure;
                            
                            const int idx_velocity = (i + num_subghosts_0_velocity) +
                                (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                            
                            const int idx_convective_flux_y = (i + num_subghosts_0_convective_flux_y) +
                                (j + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                            
                            F_y[0][idx_convective_flux_y] = rho_v[idx];
                            F_y[1][idx_convective_flux_y] = v[idx_velocity]*rho_u[idx];
                            F_y[2][idx_convective_flux_y] = v[idx_velocity]*rho_v[idx] + p[idx_pressure];
                            F_y[3][idx_convective_flux_y] = v[idx_velocity]*(E[idx] + p[idx_pressure]);
                        }
                    }
                }
                else if (d_dim == tbox::Dimension(3))
                {
                    /*
                     * Get the local lower indices, numbers of cells in each dimension and numbers of ghost cells.
                     */
                    
                    const int domain_lo_0 = domain_lo[0];
                    const int domain_lo_1 = domain_lo[1];
                    const int domain_lo_2 = domain_lo[2];
                    const int domain_dim_0 = domain_dims[0];
                    const int domain_dim_1 = domain_dims[1];
                    const int domain_dim_2 = domain_dims[2];
                    
                    const int num_ghosts_0 = d_num_ghosts[0];
                    const int num_ghosts_1 = d_num_ghosts[1];
                    const int num_ghosts_2 = d_num_ghosts[2];
                    const int ghostcell_dim_0 = d_ghostcell_dims[0];
                    const int ghostcell_dim_1 = d_ghostcell_dims[1];
                    
                    const int num_subghosts_0_pressure = d_num_subghosts_pressure[0];
                    const int num_subghosts_1_pressure = d_num_subghosts_pressure[1];
                    const int num_subghosts_2_pressure = d_num_subghosts_pressure[2];
                    const int subghostcell_dim_0_pressure = d_subghostcell_dims_pressure[0];
                    const int subghostcell_dim_1_pressure = d_subghostcell_dims_pressure[1];
                    
                    const int num_subghosts_0_velocity = d_num_subghosts_velocity[0];
                    const int num_subghosts_1_velocity = d_num_subghosts_velocity[1];
                    const int num_subghosts_2_velocity = d_num_subghosts_velocity[2];
                    const int subghostcell_dim_0_velocity = d_subghostcell_dims_velocity[0];
                    const int subghostcell_dim_1_velocity = d_subghostcell_dims_velocity[1];
                    
                    const int num_subghosts_0_convective_flux_y = d_num_subghosts_convective_flux_y[0];
                    const int num_subghosts_1_convective_flux_y = d_num_subghosts_convective_flux_y[1];
                    const int num_subghosts_2_convective_flux_y = d_num_subghosts_convective_flux_y[2];
                    const int subghostcell_dim_0_convective_flux_y = d_subghostcell_dims_convective_flux_y[0];
                    const int subghostcell_dim_1_convective_flux_y = d_subghostcell_dims_convective_flux_y[1];
                    
                    // Get the pointers to the cell data of momentum.
                    double* rho_u = data_momentum->getPointer(0);
                    double* rho_v = data_momentum->getPointer(1);
                    double* rho_w = data_momentum->getPointer(2);
                    
                    // Get the pointer to the cell data of velocity.
                    double* v = d_data_velocity->getPointer(1);
                    
                    // Compute the convective flux in the y-direction.
                    for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                    {
                        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                            {
                                // Compute the linear indices.
                                const int idx = (i + num_ghosts_0) +
                                    (j + num_ghosts_1)*ghostcell_dim_0 +
                                    (k + num_ghosts_2)*ghostcell_dim_0*ghostcell_dim_1;
                                
                                const int idx_pressure = (i + num_subghosts_0_pressure) +
                                    (j + num_subghosts_1_pressure)*subghostcell_dim_0_pressure +
                                    (k + num_subghosts_2_pressure)*subghostcell_dim_0_pressure*
                                        subghostcell_dim_1_pressure;
                                
                                const int idx_velocity = (i + num_subghosts_0_velocity) +
                                    (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                    (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                        subghostcell_dim_1_velocity;
                                
                                const int idx_convective_flux_y = (i + num_subghosts_0_convective_flux_y) +
                                    (j + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                    (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                        subghostcell_dim_1_convective_flux_y;
                                
                                F_y[0][idx_convective_flux_y] = rho_v[idx];
                                F_y[1][idx_convective_flux_y] = v[idx_velocity]*rho_u[idx];
                                F_y[2][idx_convective_flux_y] = v[idx_velocity]*rho_v[idx] + p[idx_pressure];
                                F_y[3][idx_convective_flux_y] = v[idx_velocity]*rho_w[idx];
                                F_y[4][idx_convective_flux_y] = v[idx_velocity]*(E[idx] + p[idx_pressure]);
                            }
                        }
                    }
                }
                
                d_cell_data_convective_flux_y_computed = true;
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::computeCellDataOfConvectiveFluxWithVelocityAndPressure()\n"
                << "Cell data of 'CONVECTIVE_FLUX_Y' is not yet registered."
                << std::endl);
        }
    }
    else if (direction == DIRECTION::Z_DIRECTION)
    {
        if (d_num_subghosts_convective_flux_z > -hier::IntVector::getOne(d_dim))
        {
            if (!d_cell_data_convective_flux_z_computed)
            {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(d_data_convective_flux_z);
#endif
                
                /*
                 * Get the local lower index and number of cells in each direction of the domain.
                 */
                
                hier::IntVector domain_lo(d_dim);
                hier::IntVector domain_dims(d_dim);
                
                if (domain.empty())
                {
                    domain_lo = -d_num_subghosts_convective_flux_z;
                    domain_dims = d_subghostcell_dims_convective_flux_z;
                }
                else
                {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
                    TBOX_ASSERT(d_subghost_box_convective_flux_z.contains(domain));
#endif
                    
                    domain_lo = domain.lower() - d_interior_box.lower();
                    domain_dims = domain.numberCells();
                }
                
                // Get the pointers to the components of the convective flux in the z-direction.
                std::vector<double*> F_z;
                F_z.reserve(d_num_eqn);
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    F_z.push_back(d_data_convective_flux_z->getPointer(ei));
                }
                
                boost::shared_ptr<pdat::CellData<double> > data_momentum =
                    getCellDataOfMomentum();
                
                boost::shared_ptr<pdat::CellData<double> > data_total_energy =
                    getCellDataOfTotalEnergy();
                
                if (!d_cell_data_velocity_computed)
                {
                    computeCellDataOfVelocity(domain);
                }
                
                if (!d_cell_data_pressure_computed)
                {
                    computeCellDataOfPressureWithInternalEnergy(domain);
                }
                
                // Get the pointers to the cell data of total energy and pressure.
                double* E = data_total_energy->getPointer(0);
                double* p = d_data_pressure->getPointer(0);
                
                if (d_dim == tbox::Dimension(1) || d_dim == tbox::Dimension(2))
                {
                    TBOX_ERROR(d_object_name
                        << ": FlowModelSingleSpecies::computeCellDataOfConvectiveFluxWithVelocityAndPressure()\n"
                        << "'CONVECTIVE_FLUX_Z' cannot be obtained for problem with dimension less than three."
                        << std::endl);
                }
                else if (d_dim == tbox::Dimension(3))
                {
                    /*
                     * Get the local lower indices, numbers of cells in each dimension and numbers of ghost cells.
                     */
                    
                    const int domain_lo_0 = domain_lo[0];
                    const int domain_lo_1 = domain_lo[1];
                    const int domain_lo_2 = domain_lo[2];
                    const int domain_dim_0 = domain_dims[0];
                    const int domain_dim_1 = domain_dims[1];
                    const int domain_dim_2 = domain_dims[2];
                    
                    const int num_ghosts_0 = d_num_ghosts[0];
                    const int num_ghosts_1 = d_num_ghosts[1];
                    const int num_ghosts_2 = d_num_ghosts[2];
                    const int ghostcell_dim_0 = d_ghostcell_dims[0];
                    const int ghostcell_dim_1 = d_ghostcell_dims[1];
                    
                    const int num_subghosts_0_pressure = d_num_subghosts_pressure[0];
                    const int num_subghosts_1_pressure = d_num_subghosts_pressure[1];
                    const int num_subghosts_2_pressure = d_num_subghosts_pressure[2];
                    const int subghostcell_dim_0_pressure = d_subghostcell_dims_pressure[0];
                    const int subghostcell_dim_1_pressure = d_subghostcell_dims_pressure[1];
                    
                    const int num_subghosts_0_velocity = d_num_subghosts_velocity[0];
                    const int num_subghosts_1_velocity = d_num_subghosts_velocity[1];
                    const int num_subghosts_2_velocity = d_num_subghosts_velocity[2];
                    const int subghostcell_dim_0_velocity = d_subghostcell_dims_velocity[0];
                    const int subghostcell_dim_1_velocity = d_subghostcell_dims_velocity[1];
                    
                    const int num_subghosts_0_convective_flux_z = d_num_subghosts_convective_flux_z[0];
                    const int num_subghosts_1_convective_flux_z = d_num_subghosts_convective_flux_z[1];
                    const int num_subghosts_2_convective_flux_z = d_num_subghosts_convective_flux_z[2];
                    const int subghostcell_dim_0_convective_flux_z = d_subghostcell_dims_convective_flux_z[0];
                    const int subghostcell_dim_1_convective_flux_z = d_subghostcell_dims_convective_flux_z[1];
                    
                    // Get the pointers to the cell data of momentum.
                    double* rho_u = data_momentum->getPointer(0);
                    double* rho_v = data_momentum->getPointer(1);
                    double* rho_w = data_momentum->getPointer(2);
                    
                    // Get the pointer to the cell data of velocity.
                    double* w = d_data_velocity->getPointer(2);
                    
                    // Compute the convective flux in the z-direction.
                    for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                    {
                        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                            {
                                // Compute the linear indices.
                                const int idx = (i + num_ghosts_0) +
                                    (j + num_ghosts_1)*ghostcell_dim_0 +
                                    (k + num_ghosts_2)*ghostcell_dim_0*ghostcell_dim_1;
                                
                                const int idx_pressure = (i + num_subghosts_0_pressure) +
                                    (j + num_subghosts_1_pressure)*subghostcell_dim_0_pressure +
                                    (k + num_subghosts_2_pressure)*subghostcell_dim_0_pressure*
                                        subghostcell_dim_1_pressure;
                                
                                const int idx_velocity = (i + num_subghosts_0_velocity) +
                                    (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                    (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                        subghostcell_dim_1_velocity;
                                
                                const int idx_convective_flux_z = (i + num_subghosts_0_convective_flux_z) +
                                    (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                    (k + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                        subghostcell_dim_1_convective_flux_z;
                                
                                F_z[0][idx_convective_flux_z] = rho_w[idx];
                                F_z[1][idx_convective_flux_z] = w[idx_velocity]*rho_u[idx];
                                F_z[2][idx_convective_flux_z] = w[idx_velocity]*rho_v[idx];
                                F_z[3][idx_convective_flux_z] = w[idx_velocity]*rho_w[idx] + p[idx_pressure];
                                F_z[4][idx_convective_flux_z] = w[idx_velocity]*(E[idx] + p[idx_pressure]);
                            }
                        }
                    }
                }
                
                d_cell_data_convective_flux_z_computed = true;
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::computeCellDataOfConvectiveFluxWithVelocityAndPressure()\n"
                << "Cell data of 'CONVECTIVE_FLUX_Z' is not yet registered."
                << std::endl);
        }
    }
}


/*
 * Compute the cell data of maximum wave speed with velocity and sound speed in the registered patch.
 */
void
FlowModelSingleSpecies::computeCellDataOfMaxWaveSpeedWithVelocityAndSoundSpeed(
    const DIRECTION::TYPE& direction,
    const hier::Box& domain)
{
    if (direction == DIRECTION::X_DIRECTION)
    {
        if (d_num_subghosts_max_wave_speed_x > -hier::IntVector::getOne(d_dim))
        {
            if (!d_cell_data_max_wave_speed_x_computed)
            {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(d_data_max_wave_speed_x);
#endif
                
                /*
                 * Get the local lower index and number of cells in each direction of the domain.
                 */
                
                hier::IntVector domain_lo(d_dim);
                hier::IntVector domain_dims(d_dim);
                
                if (domain.empty())
                {
                    domain_lo = -d_num_subghosts_max_wave_speed_x;
                    domain_dims = d_subghostcell_dims_max_wave_speed_x;
                }
                else
                {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
                    TBOX_ASSERT(d_subghost_box_max_wave_speed_x.contains(domain));
#endif
                    
                    domain_lo = domain.lower() - d_interior_box.lower();
                    domain_dims = domain.numberCells();
                }
                
                if (!d_cell_data_velocity_computed)
                {
                    computeCellDataOfVelocity(domain);
                }
                
                if (!d_cell_data_sound_speed_computed)
                {
                    computeCellDataOfSoundSpeedWithPressure(domain);
                }
                
                // Get the pointers to the cell data of maximum wave speed and velocity in x-direction, and sound speed.
                double* lambda_max_x = d_data_max_wave_speed_x->getPointer(0);
                double* u            = d_data_velocity->getPointer(0);
                double* c            = d_data_sound_speed->getPointer(0);
                
                if (d_dim == tbox::Dimension(1))
                {
                    /*
                     * Get the local lower index, numbers of cells in each dimension and numbers of ghost cells.
                     */
                    
                    const int domain_lo_0 = domain_lo[0];
                    const int domain_dim_0 = domain_dims[0];
                    
                    const int num_subghosts_0_sound_speed = d_num_subghosts_sound_speed[0];
                    const int num_subghosts_0_velocity = d_num_subghosts_velocity[0];
                    const int num_subghosts_0_max_wave_speed_x = d_num_subghosts_max_wave_speed_x[0];
                    
                    // Compute the maximum wave speed in the x-direction.
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_sound_speed = i + num_subghosts_0_sound_speed;
                        const int idx_velocity = i + num_subghosts_0_velocity;
                        const int idx_max_wave_speed_x = i + num_subghosts_0_max_wave_speed_x;
                        
                        lambda_max_x[idx_max_wave_speed_x] = fabs(u[idx_velocity]) + c[idx_sound_speed];
                    }
                }
                else if (d_dim == tbox::Dimension(2))
                {
                    /*
                     * Get the local lower indices, numbers of cells in each dimension and numbers of ghost cells.
                     */
                    
                    const int domain_lo_0 = domain_lo[0];
                    const int domain_lo_1 = domain_lo[1];
                    const int domain_dim_0 = domain_dims[0];
                    const int domain_dim_1 = domain_dims[1];
                    
                    const int num_subghosts_0_sound_speed = d_num_subghosts_sound_speed[0];
                    const int num_subghosts_1_sound_speed = d_num_subghosts_sound_speed[1];
                    const int subghostcell_dim_0_sound_speed = d_subghostcell_dims_sound_speed[0];
                    
                    const int num_subghosts_0_velocity = d_num_subghosts_velocity[0];
                    const int num_subghosts_1_velocity = d_num_subghosts_velocity[1];
                    const int subghostcell_dim_0_velocity = d_subghostcell_dims_velocity[0];
                    
                    const int num_subghosts_0_max_wave_speed_x = d_num_subghosts_max_wave_speed_x[0];
                    const int num_subghosts_1_max_wave_speed_x = d_num_subghosts_max_wave_speed_x[1];
                    const int subghostcell_dim_0_max_wave_speed_x = d_subghostcell_dims_max_wave_speed_x[0];
                    
                    // Compute the maximum wave speed in the x-direction.
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_sound_speed = (i + num_subghosts_0_sound_speed) +
                                (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed;
                            
                            const int idx_velocity = (i + num_subghosts_0_velocity) +
                                (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                            
                            const int idx_max_wave_speed_x = (i + num_subghosts_0_max_wave_speed_x) +
                                (j + num_subghosts_1_max_wave_speed_x)*subghostcell_dim_0_max_wave_speed_x;
                            
                            lambda_max_x[idx_max_wave_speed_x] = fabs(u[idx_velocity]) + c[idx_sound_speed];
                        }
                    }
                }
                else if (d_dim == tbox::Dimension(3))
                {
                    /*
                     * Get the local lower indices, numbers of cells in each dimension and numbers of ghost cells.
                     */
                    
                    const int domain_lo_0 = domain_lo[0];
                    const int domain_lo_1 = domain_lo[1];
                    const int domain_lo_2 = domain_lo[2];
                    const int domain_dim_0 = domain_dims[0];
                    const int domain_dim_1 = domain_dims[1];
                    const int domain_dim_2 = domain_dims[2];
                    
                    const int num_subghosts_0_sound_speed = d_num_subghosts_sound_speed[0];
                    const int num_subghosts_1_sound_speed = d_num_subghosts_sound_speed[1];
                    const int num_subghosts_2_sound_speed = d_num_subghosts_sound_speed[2];
                    const int subghostcell_dim_0_sound_speed = d_subghostcell_dims_sound_speed[0];
                    const int subghostcell_dim_1_sound_speed = d_subghostcell_dims_sound_speed[1];
                    
                    const int num_subghosts_0_velocity = d_num_subghosts_velocity[0];
                    const int num_subghosts_1_velocity = d_num_subghosts_velocity[1];
                    const int num_subghosts_2_velocity = d_num_subghosts_velocity[2];
                    const int subghostcell_dim_0_velocity = d_subghostcell_dims_velocity[0];
                    const int subghostcell_dim_1_velocity = d_subghostcell_dims_velocity[1];
                    
                    const int num_subghosts_0_max_wave_speed_x = d_num_subghosts_max_wave_speed_x[0];
                    const int num_subghosts_1_max_wave_speed_x = d_num_subghosts_max_wave_speed_x[1];
                    const int num_subghosts_2_max_wave_speed_x = d_num_subghosts_max_wave_speed_x[2];
                    const int subghostcell_dim_0_max_wave_speed_x = d_subghostcell_dims_max_wave_speed_x[0];
                    const int subghostcell_dim_1_max_wave_speed_x = d_subghostcell_dims_max_wave_speed_x[1];
                    
                    // Compute the maximum wave speed in the x-direction.
                    for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                    {
                        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                            {
                                // Compute the linear indices.
                                const int idx_sound_speed = (i + num_subghosts_0_sound_speed) +
                                    (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed +
                                    (k + num_subghosts_2_sound_speed)*subghostcell_dim_0_sound_speed*
                                        subghostcell_dim_1_sound_speed;
                                
                                const int idx_velocity = (i + num_subghosts_0_velocity) +
                                    (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                    (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                        subghostcell_dim_1_velocity;
                                
                                const int idx_max_wave_speed_x = (i + num_subghosts_0_max_wave_speed_x) +
                                    (j + num_subghosts_1_max_wave_speed_x)*subghostcell_dim_0_max_wave_speed_x +
                                    (k + num_subghosts_2_max_wave_speed_x)*subghostcell_dim_0_max_wave_speed_x*
                                        subghostcell_dim_1_max_wave_speed_x;
                                
                                lambda_max_x[idx_max_wave_speed_x] = fabs(u[idx_velocity]) + c[idx_sound_speed];
                            }
                        }
                    }
                }
                
                d_cell_data_max_wave_speed_x_computed = true;
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::computeCellDataOfMaxWaveSpeedWithVelocityAndSoundSpeed()\n"
                << "Cell data of 'MAX_WAVE_SPEED_X' is not yet registered."
                << std::endl);
        }
    }
    else if (direction == DIRECTION::Y_DIRECTION)
    {
        if (d_num_subghosts_max_wave_speed_y > -hier::IntVector::getOne(d_dim))
        {
            if (!d_cell_data_max_wave_speed_y_computed)
            {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(d_data_max_wave_speed_y);
#endif
                
                /*
                 * Get the local lower index and number of cells in each direction of the domain.
                 */
                
                hier::IntVector domain_lo(d_dim);
                hier::IntVector domain_dims(d_dim);
                
                if (domain.empty())
                {
                    domain_lo = -d_num_subghosts_max_wave_speed_y;
                    domain_dims = d_subghostcell_dims_max_wave_speed_y;
                }
                else
                {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
                    TBOX_ASSERT(d_subghost_box_max_wave_speed_y.contains(domain));
#endif
                    
                    domain_lo = domain.lower() - d_interior_box.lower();
                    domain_dims = domain.numberCells();
                }
                
                if (!d_cell_data_velocity_computed)
                {
                    computeCellDataOfVelocity(domain);
                }
                
                if (!d_cell_data_sound_speed_computed)
                {
                    computeCellDataOfSoundSpeedWithPressure(domain);
                }
                
                // Get the pointers to the cell data of maximum wave speed and velocity in y-direction, and sound speed.
                double* lambda_max_y = d_data_max_wave_speed_y->getPointer(0);
                double* v            = d_data_velocity->getPointer(1);
                double* c            = d_data_sound_speed->getPointer(0);
                
                if (d_dim == tbox::Dimension(1))
                {
                    TBOX_ERROR(d_object_name
                        << ": FlowModelSingleSpecies::computeCellDataOfMaxWaveSpeedWithVelocityAndSoundSpeed()\n"
                        << "'MAX_WAVE_SPEED_Y' cannot be obtained for problem with dimension less than two."
                        << std::endl);
                }
                else if (d_dim == tbox::Dimension(2))
                {
                    /*
                     * Get the local lower indices, numbers of cells in each dimension and numbers of ghost cells.
                     */
                    
                    const int domain_lo_0 = domain_lo[0];
                    const int domain_lo_1 = domain_lo[1];
                    const int domain_dim_0 = domain_dims[0];
                    const int domain_dim_1 = domain_dims[1];
                    
                    const int num_subghosts_0_sound_speed = d_num_subghosts_sound_speed[0];
                    const int num_subghosts_1_sound_speed = d_num_subghosts_sound_speed[1];
                    const int subghostcell_dim_0_sound_speed = d_subghostcell_dims_sound_speed[0];
                    
                    const int num_subghosts_0_velocity = d_num_subghosts_velocity[0];
                    const int num_subghosts_1_velocity = d_num_subghosts_velocity[1];
                    const int subghostcell_dim_0_velocity = d_subghostcell_dims_velocity[0];
                    
                    const int num_subghosts_0_max_wave_speed_y = d_num_subghosts_max_wave_speed_y[0];
                    const int num_subghosts_1_max_wave_speed_y = d_num_subghosts_max_wave_speed_y[1];
                    const int subghostcell_dim_0_max_wave_speed_y = d_subghostcell_dims_max_wave_speed_y[0];
                    
                    // Compute the maximum wave speed in the y-direction.
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_sound_speed = (i + num_subghosts_0_sound_speed) +
                                (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed;
                            
                            const int idx_velocity = (i + num_subghosts_0_velocity) +
                                (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                            
                            const int idx_max_wave_speed_y = (i + num_subghosts_0_max_wave_speed_y) +
                                (j + num_subghosts_1_max_wave_speed_y)*subghostcell_dim_0_max_wave_speed_y;
                            
                            lambda_max_y[idx_max_wave_speed_y] = fabs(v[idx_velocity]) + c[idx_sound_speed];
                        }
                    }
                }
                else if (d_dim == tbox::Dimension(3))
                {
                    /*
                     * Get the local lower indices, numbers of cells in each dimension and numbers of ghost cells.
                     */
                    
                    const int domain_lo_0 = domain_lo[0];
                    const int domain_lo_1 = domain_lo[1];
                    const int domain_lo_2 = domain_lo[2];
                    const int domain_dim_0 = domain_dims[0];
                    const int domain_dim_1 = domain_dims[1];
                    const int domain_dim_2 = domain_dims[2];
                    
                    const int num_subghosts_0_sound_speed = d_num_subghosts_sound_speed[0];
                    const int num_subghosts_1_sound_speed = d_num_subghosts_sound_speed[1];
                    const int num_subghosts_2_sound_speed = d_num_subghosts_sound_speed[2];
                    const int subghostcell_dim_0_sound_speed = d_subghostcell_dims_sound_speed[0];
                    const int subghostcell_dim_1_sound_speed = d_subghostcell_dims_sound_speed[1];
                    
                    const int num_subghosts_0_velocity = d_num_subghosts_velocity[0];
                    const int num_subghosts_1_velocity = d_num_subghosts_velocity[1];
                    const int num_subghosts_2_velocity = d_num_subghosts_velocity[2];
                    const int subghostcell_dim_0_velocity = d_subghostcell_dims_velocity[0];
                    const int subghostcell_dim_1_velocity = d_subghostcell_dims_velocity[1];
                    
                    const int num_subghosts_0_max_wave_speed_y = d_num_subghosts_max_wave_speed_y[0];
                    const int num_subghosts_1_max_wave_speed_y = d_num_subghosts_max_wave_speed_y[1];
                    const int num_subghosts_2_max_wave_speed_y = d_num_subghosts_max_wave_speed_y[2];
                    const int subghostcell_dim_0_max_wave_speed_y = d_subghostcell_dims_max_wave_speed_y[0];
                    const int subghostcell_dim_1_max_wave_speed_y = d_subghostcell_dims_max_wave_speed_y[1];
                    
                    // Compute the maximum wave speed in the y-direction.
                    for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                    {
                        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                            {
                                // Compute the linear indices.
                                const int idx_sound_speed = (i + num_subghosts_0_sound_speed) +
                                    (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed +
                                    (k + num_subghosts_2_sound_speed)*subghostcell_dim_0_sound_speed*
                                        subghostcell_dim_1_sound_speed;
                                
                                const int idx_velocity = (i + num_subghosts_0_velocity) +
                                    (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                    (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                        subghostcell_dim_1_velocity;
                                
                                const int idx_max_wave_speed_y = (i + num_subghosts_0_max_wave_speed_y) +
                                    (j + num_subghosts_1_max_wave_speed_y)*subghostcell_dim_0_max_wave_speed_y +
                                    (k + num_subghosts_2_max_wave_speed_y)*subghostcell_dim_0_max_wave_speed_y*
                                        subghostcell_dim_1_max_wave_speed_y;
                                
                                lambda_max_y[idx_max_wave_speed_y] = fabs(v[idx_velocity]) + c[idx_sound_speed];
                            }
                        }
                    }
                }
                
                d_cell_data_max_wave_speed_y_computed = true;
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::computeCellDataOfMaxWaveSpeedWithVelocityAndSoundSpeed()\n"
                << "Cell data of 'MAX_WAVE_SPEED_Y' is not yet registered."
                << std::endl);
        }
    }
    else if (direction == DIRECTION::Z_DIRECTION)
    {
        if (d_num_subghosts_max_wave_speed_z > -hier::IntVector::getOne(d_dim))
        {
            if (!d_cell_data_max_wave_speed_z_computed)
            {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(d_data_max_wave_speed_z);
#endif
                
                /*
                 * Get the local lower index and number of cells in each direction of the domain.
                 */
                
                hier::IntVector domain_lo(d_dim);
                hier::IntVector domain_dims(d_dim);
                
                if (domain.empty())
                {
                    domain_lo = -d_num_subghosts_max_wave_speed_z;
                    domain_dims = d_subghostcell_dims_max_wave_speed_z;
                }
                else
                {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
                    TBOX_ASSERT(d_subghost_box_max_wave_speed_z.contains(domain));
#endif
                    
                    domain_lo = domain.lower() - d_interior_box.lower();
                    domain_dims = domain.numberCells();
                }
                
                if (!d_cell_data_velocity_computed)
                {
                    computeCellDataOfVelocity(domain);
                }
                
                if (!d_cell_data_sound_speed_computed)
                {
                    computeCellDataOfSoundSpeedWithPressure(domain);
                }
                
                // Get the pointers to the cell data of maximum wave speed and velocity in z-direction, and sound speed.
                double* lambda_max_z = d_data_max_wave_speed_z->getPointer(0);
                double* w            = d_data_velocity->getPointer(2);
                double* c            = d_data_sound_speed->getPointer(0);
                
                if (d_dim == tbox::Dimension(1) || d_dim == tbox::Dimension(2))
                {
                    TBOX_ERROR(d_object_name
                        << ": FlowModelSingleSpecies::computeCellDataOfMaxWaveSpeedWithVelocityAndSoundSpeed()\n"
                        << "'MAX_WAVE_SPEED_Z' cannot be obtained for problem with dimension less than three."
                        << std::endl);
                }
                else if (d_dim == tbox::Dimension(3))
                {
                    /*
                     * Get the local lower indices, numbers of cells in each dimension and numbers of ghost cells.
                     */
                    
                    const int domain_lo_0 = domain_lo[0];
                    const int domain_lo_1 = domain_lo[1];
                    const int domain_lo_2 = domain_lo[2];
                    const int domain_dim_0 = domain_dims[0];
                    const int domain_dim_1 = domain_dims[1];
                    const int domain_dim_2 = domain_dims[2];
                    
                    const int num_subghosts_0_sound_speed = d_num_subghosts_sound_speed[0];
                    const int num_subghosts_1_sound_speed = d_num_subghosts_sound_speed[1];
                    const int num_subghosts_2_sound_speed = d_num_subghosts_sound_speed[2];
                    const int subghostcell_dim_0_sound_speed = d_subghostcell_dims_sound_speed[0];
                    const int subghostcell_dim_1_sound_speed = d_subghostcell_dims_sound_speed[1];
                    
                    const int num_subghosts_0_velocity = d_num_subghosts_velocity[0];
                    const int num_subghosts_1_velocity = d_num_subghosts_velocity[1];
                    const int num_subghosts_2_velocity = d_num_subghosts_velocity[2];
                    const int subghostcell_dim_0_velocity = d_subghostcell_dims_velocity[0];
                    const int subghostcell_dim_1_velocity = d_subghostcell_dims_velocity[1];
                    
                    const int num_subghosts_0_max_wave_speed_z = d_num_subghosts_max_wave_speed_z[0];
                    const int num_subghosts_1_max_wave_speed_z = d_num_subghosts_max_wave_speed_z[1];
                    const int num_subghosts_2_max_wave_speed_z = d_num_subghosts_max_wave_speed_z[2];
                    const int subghostcell_dim_0_max_wave_speed_z = d_subghostcell_dims_max_wave_speed_z[0];
                    const int subghostcell_dim_1_max_wave_speed_z = d_subghostcell_dims_max_wave_speed_z[1];
                    
                    // Compute the maximum wave speed in the z-direction.
                    for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                    {
                        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                            {
                                // Compute the linear indices.
                                const int idx_sound_speed = (i + num_subghosts_0_sound_speed) +
                                    (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed +
                                    (k + num_subghosts_2_sound_speed)*subghostcell_dim_0_sound_speed*
                                        subghostcell_dim_1_sound_speed;
                                
                                const int idx_velocity = (i + num_subghosts_0_velocity) +
                                    (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                    (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                        subghostcell_dim_1_velocity;
                                
                                const int idx_max_wave_speed_z = (i + num_subghosts_0_max_wave_speed_z) +
                                    (j + num_subghosts_1_max_wave_speed_z)*subghostcell_dim_0_max_wave_speed_z +
                                    (k + num_subghosts_2_max_wave_speed_z)*subghostcell_dim_0_max_wave_speed_z*
                                        subghostcell_dim_1_max_wave_speed_z;
                                
                                lambda_max_z[idx_max_wave_speed_z] = fabs(w[idx_velocity]) + c[idx_sound_speed];
                            }
                        }
                    }
                }
                
                d_cell_data_max_wave_speed_z_computed = true;
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::computeCellDataOfMaxWaveSpeedWithVelocityAndSoundSpeed()\n"
                << "Cell data of 'MAX_WAVE_SPEED_Z' is not yet registered."
                << std::endl);
        }
    }
}


/*
 * Compute the cell data of maximum diffusivity with pressure and temperature in the registered patch.
 */
void
FlowModelSingleSpecies::computeCellDataOfMaxDiffusivityWithPressureAndTemperature(
    const hier::Box& domain)
{
    if (!d_equation_of_shear_viscosity_mixing_rules ||
        !d_equation_of_bulk_viscosity_mixing_rules ||
        !d_equation_of_thermal_conductivity_mixing_rules)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "computeCellDataOfMaxDiffusivityWithPressureAndTemperature()\n"
            << "Either mixing rule of shear viscosity, bulk viscosity or"
            << " thermal conductivity is not initialized."
            << std::endl);
    }
    
    if (d_num_subghosts_max_diffusivity > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_max_diffusivity_computed)
        {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(d_data_max_diffusivity);
#endif
            
            /*
             * Get the local lower index and number of cells in each direction of the domain.
             */
            
            hier::IntVector domain_lo(d_dim);
            hier::IntVector domain_dims(d_dim);
            
            if (domain.empty())
            {
                domain_lo = -d_num_subghosts_max_diffusivity;
                domain_dims = d_subghostcell_dims_max_diffusivity;
            }
            else
            {
    #ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(d_subghost_box_max_diffusivity.contains(domain));
    #endif
                
                domain_lo = domain.lower() - d_interior_box.lower();
                domain_dims = domain.numberCells();
            }
            
            // Get the cell data of the variable density and pressure.
            boost::shared_ptr<pdat::CellData<double> > data_density =
                getCellDataOfDensity();
            
            if (!d_cell_data_pressure_computed)
            {
                computeCellDataOfPressureWithInternalEnergy(domain);
            }
            
            if (!d_cell_data_temperature_computed)
            {
                computeCellDataOfTemperatureWithPressure(domain);
            }
            
            /*
             * Create temporary cell data of isobaric specific heat capacity, shear viscosity, bulk
             * viscosity and thermal conductivity.
             */
            
            boost::shared_ptr<pdat::CellData<double> > data_isobaric_specific_heat_capacity(
                new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_max_diffusivity));
            
            boost::shared_ptr<pdat::CellData<double> > data_shear_viscosity(
                new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_max_diffusivity));
            
            boost::shared_ptr<pdat::CellData<double> > data_bulk_viscosity(
                new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_max_diffusivity));
            
            boost::shared_ptr<pdat::CellData<double> > data_thermal_conductivity(
                new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_max_diffusivity));
            
            // Get the pointers to the cell data of maximum diffusivity, density, isobaric specific heat
            // capacity, shear viscosity, bulk viscosity and thermal conductivity.
            double* D_max = d_data_max_diffusivity->getPointer(0);
            double* rho   = data_density->getPointer(0);
            double* c_p   = data_isobaric_specific_heat_capacity->getPointer(0);
            double* mu    = data_shear_viscosity->getPointer(0);
            double* mu_v  = data_bulk_viscosity->getPointer(0);
            double* kappa = data_thermal_conductivity->getPointer(0);
            
            // Get the thermodynamic properties of the species.
            std::vector<const double*> thermo_properties_ptr;
            thermo_properties_ptr.reserve(static_cast<int> (d_thermo_properties.size()));
            for (int ti = 0; ti < static_cast<int> (d_thermo_properties.size()); ti++)
            {
                thermo_properties_ptr.push_back(&d_thermo_properties[ti]);
            }
            
            // Get the molecular properties of the species for shear viscosity.
            std::vector<const double*> molecular_properties_shear_viscosity_ptr;
            molecular_properties_shear_viscosity_ptr.reserve(
                static_cast<int> (d_molecular_properties_shear_viscosity.size()));
            for (int ti = 0; ti < static_cast<int> (d_molecular_properties_shear_viscosity.size()); ti++)
            {
                molecular_properties_shear_viscosity_ptr.push_back(
                    &d_molecular_properties_shear_viscosity[ti]);
            }
            
            // Get the molecular properties of the species for bulk viscosity.
            std::vector<const double*> molecular_properties_bulk_viscosity_ptr;
            molecular_properties_bulk_viscosity_ptr.reserve(
                static_cast<int> (d_molecular_properties_bulk_viscosity.size()));
            for (int ti = 0; ti < static_cast<int> (d_molecular_properties_bulk_viscosity.size()); ti++)
            {
                molecular_properties_bulk_viscosity_ptr.push_back(
                    &d_molecular_properties_bulk_viscosity[ti]);
            }
            
            // Get the molecular properties of the species for thermal conductivity.
            std::vector<const double*> molecular_properties_thermal_conductivity_ptr;
            molecular_properties_thermal_conductivity_ptr.reserve(
                static_cast<int> (d_molecular_properties_thermal_conductivity.size()));
            for (int ti = 0; ti < static_cast<int> (d_molecular_properties_thermal_conductivity.size()); ti++)
            {
                molecular_properties_thermal_conductivity_ptr.push_back(
                    &d_molecular_properties_thermal_conductivity[ti]);
            }
            
            // Compute the isobaric specific heat capacity field.
            d_equation_of_state_mixing_rules->getEquationOfState()->
                computeIsobaricSpecificHeatCapacity(
                    data_isobaric_specific_heat_capacity,
                    data_density,
                    d_data_pressure,
                    thermo_properties_ptr,
                    domain);
            
            // Compute the shear viscosity field.
            d_equation_of_shear_viscosity_mixing_rules->getEquationOfShearViscosity()->
                computeShearViscosity(
                    data_shear_viscosity,
                    d_data_pressure,
                    d_data_temperature,
                    molecular_properties_shear_viscosity_ptr,
                    domain);
            
            // Compute the bulk viscosity field.
            d_equation_of_bulk_viscosity_mixing_rules->getEquationOfBulkViscosity()->
                computeBulkViscosity(
                    data_bulk_viscosity,
                    d_data_pressure,
                    d_data_temperature,
                    molecular_properties_bulk_viscosity_ptr,
                    domain);
            
            // Compute the thermal conductivity field.
            d_equation_of_thermal_conductivity_mixing_rules->getEquationOfThermalConductivity()->
                computeThermalConductivity(
                    data_thermal_conductivity,
                    d_data_pressure,
                    d_data_temperature,
                    molecular_properties_thermal_conductivity_ptr,
                    domain);
            
            if (d_dim == tbox::Dimension(1))
            {
                /*
                 * Get the local lower index, numbers of cells in each dimension and numbers of ghost cells.
                 */
                
                const int domain_lo_0 = domain_lo[0];
                const int domain_dim_0 = domain_dims[0];
                
                const int num_ghosts_0 = d_num_ghosts[0];
                const int num_subghosts_0_max_diffusivity = d_num_subghosts_max_diffusivity[0];
                
    #ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
    #endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx = i + num_ghosts_0;
                    const int idx_max_diffusivity = i + num_subghosts_0_max_diffusivity;
                    
                    D_max[idx_max_diffusivity] = fmax(mu[idx_max_diffusivity]/rho[idx],
                        mu_v[idx_max_diffusivity]/rho[idx]);
                    
                    D_max[idx_max_diffusivity] = fmax(D_max[idx_max_diffusivity],
                        kappa[idx_max_diffusivity]/(rho[idx]*c_p[idx_max_diffusivity]));
                }
            }
            else if (d_dim == tbox::Dimension(2))
            {
                /*
                 * Get the local lower indices, numbers of cells in each dimension and numbers of ghost cells.
                 */
                
                const int domain_lo_0 = domain_lo[0];
                const int domain_lo_1 = domain_lo[1];
                const int domain_dim_0 = domain_dims[0];
                const int domain_dim_1 = domain_dims[1];
                
                const int num_ghosts_0 = d_num_ghosts[0];
                const int num_ghosts_1 = d_num_ghosts[1];
                const int ghostcell_dim_0 = d_ghostcell_dims[0];
                
                const int num_subghosts_0_max_diffusivity = d_num_subghosts_max_diffusivity[0];
                const int num_subghosts_1_max_diffusivity = d_num_subghosts_max_diffusivity[1];
                const int subghostcell_dim_0_max_diffusivity = d_subghostcell_dims_max_diffusivity[0];
                
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
    #ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
    #endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx = (i + num_ghosts_0) +
                            (j + num_ghosts_1)*ghostcell_dim_0;
                        
                        const int idx_max_diffusivity = (i + num_subghosts_0_max_diffusivity) +
                            (j + num_subghosts_1_max_diffusivity)*subghostcell_dim_0_max_diffusivity;
                        
                        D_max[idx_max_diffusivity] = fmax(mu[idx_max_diffusivity]/rho[idx],
                            mu_v[idx_max_diffusivity]/rho[idx]);
                        
                        D_max[idx_max_diffusivity] = fmax(D_max[idx_max_diffusivity],
                            kappa[idx_max_diffusivity]/(rho[idx]*c_p[idx_max_diffusivity]));
                    }
                }
            }
            else if (d_dim == tbox::Dimension(3))
            {
                /*
                 * Get the local lower indices, numbers of cells in each dimension and numbers of ghost cells.
                 */
                
                const int domain_lo_0 = domain_lo[0];
                const int domain_lo_1 = domain_lo[1];
                const int domain_lo_2 = domain_lo[2];
                const int domain_dim_0 = domain_dims[0];
                const int domain_dim_1 = domain_dims[1];
                const int domain_dim_2 = domain_dims[2];
                
                const int num_ghosts_0 = d_num_ghosts[0];
                const int num_ghosts_1 = d_num_ghosts[1];
                const int num_ghosts_2 = d_num_ghosts[2];
                const int ghostcell_dim_0 = d_ghostcell_dims[0];
                const int ghostcell_dim_1 = d_ghostcell_dims[1];
                
                const int num_subghosts_0_max_diffusivity = d_num_subghosts_max_diffusivity[0];
                const int num_subghosts_1_max_diffusivity = d_num_subghosts_max_diffusivity[1];
                const int num_subghosts_2_max_diffusivity = d_num_subghosts_max_diffusivity[2];
                const int subghostcell_dim_0_max_diffusivity = d_subghostcell_dims_max_diffusivity[0];
                const int subghostcell_dim_1_max_diffusivity = d_subghostcell_dims_max_diffusivity[1];
                
                for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                {
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
    #ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
    #endif
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx = (i + num_ghosts_0) +
                                (j + num_ghosts_1)*ghostcell_dim_0 +
                                (k + num_ghosts_2)*ghostcell_dim_0*ghostcell_dim_1;
                            
                            const int idx_max_diffusivity = (i + num_subghosts_0_max_diffusivity) +
                                (j + num_subghosts_1_max_diffusivity)*subghostcell_dim_0_max_diffusivity +
                                (k + num_subghosts_2_max_diffusivity)*subghostcell_dim_0_max_diffusivity*
                                    subghostcell_dim_1_max_diffusivity;
                            
                            D_max[idx_max_diffusivity] = fmax(mu[idx_max_diffusivity]/rho[idx],
                                mu_v[idx_max_diffusivity]/rho[idx]);
                            
                            D_max[idx_max_diffusivity] = fmax(D_max[idx_max_diffusivity],
                                kappa[idx_max_diffusivity]/(rho[idx]*c_p[idx_max_diffusivity]));
                        }
                    }
                }
            }
            
            d_cell_data_max_diffusivity_computed = true;
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "computeCellDataOfMaxDiffusivityWithPressureAndTemperature()\n"
            << "Cell data of 'MAX_DIFFUSIVITY' is not yet registered."
            << std::endl);
    }
}
