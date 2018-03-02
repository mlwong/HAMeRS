#include "flow/flow_models/four-eqn_conservative/FlowModelFourEqnConservative.hpp"

#include "flow/flow_models/four-eqn_conservative/FlowModelBoundaryUtilitiesFourEqnConservative.hpp"
#include "flow/flow_models/four-eqn_conservative/FlowModelRiemannSolverFourEqnConservative.hpp"
#include "flow/flow_models/four-eqn_conservative/FlowModelStatisticsUtilitiesFourEqnConservative.hpp"

boost::shared_ptr<pdat::CellVariable<double> > FlowModelFourEqnConservative::s_variable_partial_densities;
boost::shared_ptr<pdat::CellVariable<double> > FlowModelFourEqnConservative::s_variable_momentum;
boost::shared_ptr<pdat::CellVariable<double> > FlowModelFourEqnConservative::s_variable_total_energy;

FlowModelFourEqnConservative::FlowModelFourEqnConservative(
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
            num_species + dim.getValue() + 1,
            flow_model_db),
        d_num_subghosts_density(-hier::IntVector::getOne(d_dim)),
        d_num_subghosts_mass_fractions(-hier::IntVector::getOne(d_dim)),
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
        d_num_subghosts_diffusivities(-hier::IntVector::getOne(d_dim)),
        d_subghost_box_density(hier::Box::getEmptyBox(dim)),
        d_subghost_box_mass_fractions(hier::Box::getEmptyBox(dim)),
        d_subghost_box_velocity(hier::Box::getEmptyBox(dim)),
        d_subghost_box_internal_energy(hier::Box::getEmptyBox(dim)),
        d_subghost_box_pressure(hier::Box::getEmptyBox(dim)),
        d_subghost_box_sound_speed(hier::Box::getEmptyBox(dim)),
        d_subghost_box_temperature(hier::Box::getEmptyBox(dim)),
        d_subghost_box_convective_flux_x(hier::Box::getEmptyBox(dim)),
        d_subghost_box_convective_flux_y(hier::Box::getEmptyBox(dim)),
        d_subghost_box_convective_flux_z(hier::Box::getEmptyBox(dim)),
        d_subghost_box_max_wave_speed_x(hier::Box::getEmptyBox(dim)),
        d_subghost_box_max_wave_speed_y(hier::Box::getEmptyBox(dim)),
        d_subghost_box_max_wave_speed_z(hier::Box::getEmptyBox(dim)),
        d_subghost_box_max_diffusivity(hier::Box::getEmptyBox(dim)),
        d_subghost_box_diffusivities(hier::Box::getEmptyBox(dim)),
        d_subghostcell_dims_density(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_mass_fractions(hier::IntVector::getZero(d_dim)),
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
        d_subghostcell_dims_diffusivities(hier::IntVector::getZero(d_dim))
{
    d_eqn_form.reserve(d_num_eqn);
    
    // Set the equation forms for partial densities.
    for (int si = 0; si < d_num_species; si++)
    {
        d_eqn_form.push_back(EQN_FORM::CONSERVATIVE);
    }
    
    // Set the equation forms for momentum.
    for (int di = 0; di < d_dim.getValue(); di++)
    {
        d_eqn_form.push_back(EQN_FORM::CONSERVATIVE);
    }
    
    // Set the equation form for total energy.
    d_eqn_form.push_back(EQN_FORM::CONSERVATIVE);
    
    // Set the bounds for the variables.
    d_Y_bound_lo = double(-0.001);
    d_Y_bound_up = double(1.001);
    
    /*
     * Initialize the conservative variables.
     */
    
    s_variable_partial_densities = boost::shared_ptr<pdat::CellVariable<double> > (
        new pdat::CellVariable<double>(d_dim, "partial densities", d_num_species));
    
    s_variable_momentum = boost::shared_ptr<pdat::CellVariable<double> > (
        new pdat::CellVariable<double>(d_dim, "momentum", d_dim.getValue()));
    
    s_variable_total_energy = boost::shared_ptr<pdat::CellVariable<double> > (
        new pdat::CellVariable<double>(d_dim, "total energy", 1));
    
    /*
     * Initialize d_equation_of_state_mixing_rules_manager and get the equation of state mixing rules object.
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
        MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC,
        equation_of_state_mixing_rules_db,
        d_equation_of_state_str));
    
    d_equation_of_state_mixing_rules =
        d_equation_of_state_mixing_rules_manager->getEquationOfStateMixingRules();
    
    /*
     * Initialize d_equation_of_mass_diffusivity_mixing_rules and get the equation of
     * mass diffusivity mixing rules object.
     */
    
    if ((flow_model_db->keyExists("equation_of_mass_diffusivity")) ||
        (flow_model_db->keyExists("d_equation_of_mass_diffusivity_str")))
    {
        if (flow_model_db->keyExists("equation_of_mass_diffusivity"))
        {
            d_equation_of_mass_diffusivity_str =
                flow_model_db->getString("equation_of_mass_diffusivity");
        }
        else if (flow_model_db->keyExists("d_equation_of_mass_diffusivity_str"))
        {
            d_equation_of_mass_diffusivity_str =
                flow_model_db->getString("d_equation_of_mass_diffusivity_str");
        }
        
        boost::shared_ptr<tbox::Database> equation_of_mass_diffusivity_mixing_rules_db;
        
        if (flow_model_db->keyExists("Equation_of_mass_diffusivity_mixing_rules"))
        {
            equation_of_mass_diffusivity_mixing_rules_db =
                flow_model_db->getDatabase("Equation_of_mass_diffusivity_mixing_rules");
        }
        else if (flow_model_db->keyExists("d_equation_of_mass_diffusivity_mixing_rules_db"))
        {
            equation_of_mass_diffusivity_mixing_rules_db =
                flow_model_db->getDatabase("d_equation_of_mass_diffusivity_mixing_rules_db");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "No key 'Equation_of_mass_diffusivity_mixing_rules'/"
                << "'d_equation_of_mass_diffusivity_mixing_rules_db' found in data for flow model"
                << std::endl);
        }
        
        d_equation_of_mass_diffusivity_mixing_rules_manager.reset(
            new EquationOfMassDiffusivityMixingRulesManager(
                "d_equation_of_mass_diffusivity_mixing_rules_manager",
                d_dim,
                d_num_species,
                MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC,
                equation_of_mass_diffusivity_mixing_rules_db,
                d_equation_of_mass_diffusivity_str));
        
        d_equation_of_mass_diffusivity_mixing_rules =
            d_equation_of_mass_diffusivity_mixing_rules_manager->
                getEquationOfMassDiffusivityMixingRules();
    }
    
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
                MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC,
                equation_of_shear_viscosity_mixing_rules_db,
                d_equation_of_shear_viscosity_str));
        
        d_equation_of_shear_viscosity_mixing_rules =
            d_equation_of_shear_viscosity_mixing_rules_manager->
                getEquationOfShearViscosityMixingRules();
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
                MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC,
                equation_of_bulk_viscosity_mixing_rules_db,
                d_equation_of_bulk_viscosity_str));
        
        d_equation_of_bulk_viscosity_mixing_rules =
            d_equation_of_bulk_viscosity_mixing_rules_manager->
                getEquationOfBulkViscosityMixingRules();
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
                MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC,
                equation_of_thermal_conductivity_mixing_rules_db,
                d_equation_of_thermal_conductivity_str));
        
        d_equation_of_thermal_conductivity_mixing_rules =
            d_equation_of_thermal_conductivity_mixing_rules_manager->
                getEquationOfThermalConductivityMixingRules();
    }
    
    /*
     * Initialize Riemann solver object.
     */
    d_flow_model_riemann_solver.reset(new FlowModelRiemannSolverFourEqnConservative(
        "d_flow_model_riemann_solver",
        d_dim,
        d_grid_geometry,
        d_num_species));
    
    /*
     * Initialize statistics utilities object.
     */
    d_flow_model_statistics_utilities.reset(new FlowModelStatisticsUtilitiesFourEqnConservative(
        "d_flow_model_statistics_utilities",
        d_dim,
        d_grid_geometry,
        d_num_species,
        flow_model_db,
        d_equation_of_state_mixing_rules,
        d_equation_of_mass_diffusivity_mixing_rules,
        d_equation_of_shear_viscosity_mixing_rules,
        d_equation_of_bulk_viscosity_mixing_rules,
        d_equation_of_thermal_conductivity_mixing_rules));
    
    /*
     * Initialize the Riemann solvers.
     */
    d_Riemann_solver_HLLC = boost::shared_ptr<RiemannSolverFourEqnConservativeHLLC> (
        new RiemannSolverFourEqnConservativeHLLC(
            d_object_name,
            d_dim,
            d_num_eqn,
            d_num_species,
            d_equation_of_state_mixing_rules));
    
    d_Riemann_solver_HLLC_HLL = boost::shared_ptr<RiemannSolverFourEqnConservativeHLLC_HLL> (
        new RiemannSolverFourEqnConservativeHLLC_HLL(
            d_object_name,
            d_dim,
            d_num_eqn,
            d_num_species,
            d_equation_of_state_mixing_rules));
    
    /*
     * Initialize boundary utilities object.
     */
    d_flow_model_boundary_utilities.reset(
        new FlowModelBoundaryUtilitiesFourEqnConservative(
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
FlowModelFourEqnConservative::printClassData(std::ostream& os) const
{
    os << "\nPrint FlowModelFourEqnConservative object..."
       << std::endl;
    
    os << std::endl;
    
    os << "FlowModelFourEqnConservative: this = "
       << (FlowModelFourEqnConservative *)this
       << std::endl;
    
    os << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
    d_equation_of_state_mixing_rules_manager->printClassData(os);
}


/*
 * Put the characteristics of the flow model class into the restart database.
 */
void
FlowModelFourEqnConservative::putToRestart(
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
     * Put the properties of d_equation_of_mass_diffusivity_mixing_rules into the restart database.
     */
    
    if (d_equation_of_mass_diffusivity_mixing_rules)
    {
        restart_db->putString("d_equation_of_mass_diffusivity_str", d_equation_of_mass_diffusivity_str);
        
        boost::shared_ptr<tbox::Database> restart_equation_of_mass_diffusivity_mixing_rules_db =
            restart_db->putDatabase("d_equation_of_mass_diffusivity_mixing_rules_db");
        d_equation_of_mass_diffusivity_mixing_rules->
            putToRestart(restart_equation_of_mass_diffusivity_mixing_rules_db);
    }
    
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
FlowModelFourEqnConservative::registerConservativeVariables(
    RungeKuttaLevelIntegrator* integrator,
    const hier::IntVector& num_ghosts,
    const hier::IntVector& num_ghosts_intermediate)
{
    integrator->registerVariable(
        s_variable_partial_densities,
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
FlowModelFourEqnConservative::getNamesOfConservativeVariables(bool have_underscores)
{
    std::vector<std::string> names;
    names.reserve(3);
    
    if (have_underscores)
    {
        names.push_back("partial_densities");
        names.push_back("momentum");
        names.push_back("total_energy");
    }
    else
    {
        names.push_back("partial densities");
        names.push_back("momentum");
        names.push_back("total energy");
    }
    
    return names;
}


/*
 * Get the names of primitive variables.
 */
std::vector<std::string> FlowModelFourEqnConservative::getNamesOfPrimitiveVariables(bool have_underscores)
{
    std::vector<std::string> names;
    names.reserve(3);
    
    if (have_underscores)
    {
        names.push_back("partial_densities");
        names.push_back("velocity");
        names.push_back("pressure");
    }
    else
    {
        names.push_back("partial densities");
        names.push_back("velocity");
        names.push_back("pressure");
    }
    
    return names;
}


/*
 * Get the variable types of conservative variables.
 */
std::vector<std::string>
FlowModelFourEqnConservative::getVariableTypesOfConservativeVariables()
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
FlowModelFourEqnConservative::getVariableTypesOfPrimitiveVariables()
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
FlowModelFourEqnConservative::getConservativeVariables()
{
    std::vector<boost::shared_ptr<pdat::CellVariable<double> > > conservative_variables;
    conservative_variables.reserve(3);
    
    conservative_variables.push_back(s_variable_partial_densities);
    conservative_variables.push_back(s_variable_momentum);
    conservative_variables.push_back(s_variable_total_energy);
    
    return conservative_variables;
}


/*
 * Register a patch with a data context.
 */
void
FlowModelFourEqnConservative::registerPatchWithDataContext(
    const hier::Patch& patch,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
    // Check whether the patch is already unregistered.
    if (d_patch)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::registerPatchWithDataContext()\n"
            << "The patch is not yet unregistered."
            << std::endl);
    }
    
    d_patch = &patch;
    
    setDataContext(data_context);
    
    /*
     * Set the number of ghost cells of conservative variables.
     */
    
    boost::shared_ptr<pdat::CellData<double> > data_partial_densities(
        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
            d_patch->getPatchData(s_variable_partial_densities, getDataContext())));
    
    d_num_ghosts = data_partial_densities->getGhostCellWidth();
    
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
FlowModelFourEqnConservative::registerDerivedCellVariable(
    const std::unordered_map<std::string, hier::IntVector>& num_subghosts_of_data)
{
    // Check whether a patch is already registered.
    if (!d_patch)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::registerDerivedCellVariable()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    // Check whether all or part of derived cell data is already computed.
    if (d_global_derived_cell_data_computed)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::registerDerivedCellVariable()\n"
            << "Global derived cell data is already computed."
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
                << ": FlowModelFourEqnConservative::registerDerivedCellVariable()\n"
                << "The number of sub-ghost cells of variables '"
                << it->first
                << "' is not between zero and d_num_ghosts."
                << std::endl);
        }
    }
    
    if (num_subghosts_of_data.find("DENSITY") != num_subghosts_of_data.end())
    {
        setNumberOfSubGhosts(
            num_subghosts_of_data.find("DENSITY")->second,
            "DENSITY",
            "DENSITY");
    }
    
    if (num_subghosts_of_data.find("MASS_FRACTIONS") != num_subghosts_of_data.end())
    {
        setNumberOfSubGhosts(
            num_subghosts_of_data.find("MASS_FRACTIONS")->second,
            "MASS_FRACTIONS",
            "MASS_FRACTIONS");
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
    
    if (num_subghosts_of_data.find("DILATATION") != num_subghosts_of_data.end())
    {
        setNumberOfSubGhosts(
            num_subghosts_of_data.find("DILATATION")->second,
            "DILATATION",
            "DILATATION");
    }
    
    if (num_subghosts_of_data.find("VORTICITY") != num_subghosts_of_data.end())
    {
        setNumberOfSubGhosts(
            num_subghosts_of_data.find("VORTICITY")->second,
            "VORTICITY",
            "VORTICITY");
    }
    
    if (num_subghosts_of_data.find("ENSTROPHY") != num_subghosts_of_data.end())
    {
        setNumberOfSubGhosts(
            num_subghosts_of_data.find("ENSTROPHY")->second,
            "ENSTROPHY",
            "ENSTROPHY");
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
 * Register the required derived variables for transformation between conservative
 * variables and characteristic variables.
 */
void
FlowModelFourEqnConservative::registerDerivedVariablesForCharacteristicProjectionOfConservativeVariables(
    const hier::IntVector& num_subghosts,
    const AVERAGING::TYPE& averaging)
{
    // Check whether a patch is already registered.
    if (!d_patch)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::"
            << "registerDerivedVariablesForCharacteristicProjectionOfConservativeVariables()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    // Check whether all or part of derived cell data is already computed.
    if (d_global_derived_cell_data_computed)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::"
            << "registerDerivedVariablesForCharacteristicProjectionOfConservativeVariables()\n"
            << "Global derived cell data is already computed."
            << std::endl);
    }
    
    NULL_USE(num_subghosts);
    
    d_proj_var_conservative_averaging = averaging;
}


/*
 * Register the required derived variables for transformation between primitive variables
 * and characteristic variables.
 */
void
FlowModelFourEqnConservative::registerDerivedVariablesForCharacteristicProjectionOfPrimitiveVariables(
    const hier::IntVector& num_subghosts,
    const AVERAGING::TYPE& averaging)
{
    // Check whether a patch is already registered.
    if (!d_patch)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::"
            << "registerDerivedVariablesForCharacteristicProjectionOfPrimitiveVariables()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    // Check whether all or part of derived cell data is already computed.
    if (d_global_derived_cell_data_computed)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::"
            << "registerDerivedVariablesForCharacteristicProjectionOfPrimitiveVariables()\n"
            << "Global derived cell data is already computed."
            << std::endl);
    }
    
    d_proj_var_primitive_averaging = averaging;
    
    setNumberOfSubGhosts(
        num_subghosts,
        "SOUND_SPEED",
        "PROJECTION_MATRICES");
}


/*
 * Register the required variables for the computation of diffusive flux in the
 * registered patch.
 */
void
FlowModelFourEqnConservative::registerDiffusiveFlux(const hier::IntVector& num_subghosts)
{
    // Check whether a patch is already registered.
    if (!d_patch)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::"
            << "registerDiffusiveFlux()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    // Check whether all or part of derived cell data is already computed.
    if (d_global_derived_cell_data_computed)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::"
            << "registerDiffusiveFlux()\n"
            << "Global derived cell data is already computed."
            << std::endl);
    }
    
    setNumberOfSubGhosts(
        num_subghosts,
        "DENSITY",
        "DIFFUSIVE_FLUX");
    
    setNumberOfSubGhosts(
        num_subghosts,
        "MASS_FRACTIONS",
        "DIFFUSIVE_FLUX");
    
    setNumberOfSubGhosts(
        num_subghosts,
        "VELOCITY",
        "DIFFUSIVE_FLUX");
    
    setNumberOfSubGhosts(
        num_subghosts,
        "PRESSURE",
        "DIFFUSIVE_FLUX");
    
    setNumberOfSubGhosts(
        num_subghosts,
        "TEMPERATURE",
        "DIFFUSIVE_FLUX");
    
    d_num_subghosts_diffusivities = 
        hier::IntVector::min(d_num_subghosts_density, d_num_subghosts_mass_fractions);
    
    d_num_subghosts_diffusivities =
        hier::IntVector::min(d_num_subghosts_diffusivities, d_num_subghosts_velocity);
    
    d_num_subghosts_diffusivities =
        hier::IntVector::min(d_num_subghosts_diffusivities, d_num_subghosts_pressure);
    
    d_num_subghosts_diffusivities =
        hier::IntVector::min(d_num_subghosts_diffusivities, d_num_subghosts_temperature);
}


/*
 * Unregister the registered patch. The registered data context and all global derived
 * cell data in the patch are dumped.
 */
void
FlowModelFourEqnConservative::unregisterPatch()
{
    // Check whether a patch is already registered.
    if (!d_patch)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::unregisterPatch()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    d_patch = nullptr;
    
    d_num_ghosts                      = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_density           = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_mass_fractions    = -hier::IntVector::getOne(d_dim);
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
    d_num_subghosts_diffusivities     = -hier::IntVector::getOne(d_dim);
    
    d_interior_box                   = hier::Box::getEmptyBox(d_dim);
    d_ghost_box                      = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_density           = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_mass_fractions    = hier::Box::getEmptyBox(d_dim);
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
    d_subghost_box_diffusivities     = hier::Box::getEmptyBox(d_dim);
    
    d_interior_dims                       = hier::IntVector::getZero(d_dim);
    d_ghostcell_dims                      = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_density           = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_mass_fractions    = hier::IntVector::getZero(d_dim);
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
    d_subghostcell_dims_diffusivities     = hier::IntVector::getZero(d_dim);
    
    d_data_density.reset();
    d_data_mass_fractions.reset();
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
    d_data_diffusivities.reset();
    
    d_global_derived_cell_data_computed = false;
    
    clearDataContext();
}


/*
 * Compute global cell data of different registered derived variables with the registered data context.
 */
void
FlowModelFourEqnConservative::computeGlobalDerivedCellData(const hier::Box& domain)
{
    // Check whether a patch is already registered.
    if (!d_patch)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::computeGlobalDerivedCellData()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    /*
     * Set the boxes and their dimensions for the derived cell variables.
     */
    if (!d_global_derived_cell_data_computed)
    {
        setGhostBoxesAndDimensionsDerivedCellVariables();
    }
    
    // Compute the total density cell data.
    if (d_num_subghosts_density > -hier::IntVector::getOne(d_dim))
    {
        if (!d_data_density)
        {
            computeGlobalCellDataDensity(
                domain);
        }
    }
    
    // Compute the mass fraction cell data.
    if (d_num_subghosts_mass_fractions > -hier::IntVector::getOne(d_dim))
    {
        if (!d_data_mass_fractions)
        {
            computeGlobalCellDataMassFractionsWithDensity(
                domain);
        }
    }
    
    // Compute the velocity cell data.
    if (d_num_subghosts_velocity > -hier::IntVector::getOne(d_dim))
    {
        if (!d_data_velocity)
        {
            computeGlobalCellDataVelocityWithDensity(
                domain);
        }
    }
    
    // Compute the internal energy cell data.
    if (d_num_subghosts_internal_energy > -hier::IntVector::getOne(d_dim))
    {
        if (!d_data_internal_energy)
        {
            computeGlobalCellDataInternalEnergyWithDensityAndVelocity(
                domain);
        }
    }
    
    // Compute the pressure cell data.
    if (d_num_subghosts_pressure > -hier::IntVector::getOne(d_dim))
    {
        if (!d_data_pressure)
        {
            computeGlobalCellDataPressureWithDensityMassFractionsAndInternalEnergy(
                domain);
        }
    }
    
    // Compute the sound speed cell data.
    if (d_num_subghosts_sound_speed > -hier::IntVector::getOne(d_dim))
    {
        if (!d_data_sound_speed)
        {
            computeGlobalCellDataSoundSpeedWithDensityMassFractionsAndPressure(
                domain);
        }
    }
    
    // Compute the temperature cell data.
    if (d_num_subghosts_temperature > -hier::IntVector::getOne(d_dim))
    {
        if (!d_data_temperature)
        {
            computeGlobalCellDataTemperatureWithDensityMassFractionsAndPressure(
                domain);
        }
    }
    
    // Compute the x-direction convective flux cell data.
    if (d_num_subghosts_convective_flux_x > -hier::IntVector::getOne(d_dim))
    {
        if (!d_data_convective_flux_x)
        {
            computeGlobalCellDataConvectiveFluxWithVelocityAndPressure(
                DIRECTION::X_DIRECTION,
                domain);
        }
    }
    
    // Compute the y-direction convective flux cell data.
    if (d_num_subghosts_convective_flux_y > -hier::IntVector::getOne(d_dim))
    {
        if (!d_data_convective_flux_y)
        {
            computeGlobalCellDataConvectiveFluxWithVelocityAndPressure(
                DIRECTION::Y_DIRECTION,
                domain);
        }
    }
    
    // Compute the z-direction convective flux cell data.
    if (d_num_subghosts_convective_flux_z > -hier::IntVector::getOne(d_dim))
    {
        if (!d_data_convective_flux_z)
        {
            computeGlobalCellDataConvectiveFluxWithVelocityAndPressure(
                DIRECTION::Z_DIRECTION,
                domain);
        }
    }
    
    // Compute the x-direction maximum wave speed cell data.
    if (d_num_subghosts_max_wave_speed_x > -hier::IntVector::getOne(d_dim))
    {
        if (!d_data_max_wave_speed_x)
        {
            computeGlobalCellDataMaxWaveSpeedWithVelocityAndSoundSpeed(
                DIRECTION::X_DIRECTION,
                domain);
        }
    }
    
    // Compute the y-direction maximum wave speed cell data.
    if (d_num_subghosts_max_wave_speed_y > -hier::IntVector::getOne(d_dim))
    {
        if (!d_data_max_wave_speed_y)
        {
            computeGlobalCellDataMaxWaveSpeedWithVelocityAndSoundSpeed(
                DIRECTION::Y_DIRECTION,
                domain);
        }
    }
    
    // Compute the z-direction maximum wave speed cell data.
    if (d_num_subghosts_max_wave_speed_z > -hier::IntVector::getOne(d_dim))
    {
        if (!d_data_max_wave_speed_z)
        {
            computeGlobalCellDataMaxWaveSpeedWithVelocityAndSoundSpeed(
                DIRECTION::Z_DIRECTION,
                domain);
        }
    }
    
    // Compute the maximum diffusivity cell data.
    if (d_num_subghosts_max_diffusivity > -hier::IntVector::getOne(d_dim))
    {
        if (!d_data_max_diffusivity)
        {
            computeGlobalCellDataMaxDiffusivityWithDensityMassFractionsPressureAndTemperature(
                domain);
        }
    }
    
    d_global_derived_cell_data_computed = true;
}


/*
 * Get the global cell data of one cell variable in the registered patch.
 */
boost::shared_ptr<pdat::CellData<double> >
FlowModelFourEqnConservative::getGlobalCellData(const std::string& variable_key)
{
    // Check whether a patch is already registered.
    if (!d_patch)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::getGlobalCellData()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    boost::shared_ptr<pdat::CellData<double> > cell_data;
    
    if (variable_key == "PARTIAL_DENSITY")
    {
        cell_data = getGlobalCellDataPartialDensities();
    }
    else if (variable_key == "MOMENTUM")
    {
        cell_data = getGlobalCellDataMomentum();
    }
    else if (variable_key == "TOTAL_ENERGY")
    {
        cell_data = getGlobalCellDataTotalEnergy();
    }
    else if (variable_key == "DENSITY")
    {
        if (!d_data_density)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::getGlobalCellData()\n"
                << "Cell data of 'DENSITY' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_density;
    }
    else if (variable_key == "MASS_FRACTIONS")
    {
        if (!d_data_mass_fractions)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::getGlobalCellData()\n"
                << "Cell data of 'MASS_FRACTION' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_mass_fractions;
    }
    else if (variable_key == "VELOCITY")
    {
        if (!d_data_velocity)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::getGlobalCellData()\n"
                << "Cell data of 'VELOCITY' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_velocity;
    }
    else if (variable_key == "INTERNAL_ENERGY")
    {
        if (!d_data_internal_energy)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::getGlobalCellData()\n"
                << "Cell data of 'INTERNAL_ENERGY' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_internal_energy;
    }
    else if (variable_key == "PRESSURE")
    {
        if (!d_data_pressure)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::getGlobalCellData()\n"
                << "Cell data of 'PRESSURE' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_pressure;
    }
    else if (variable_key == "SOUND_SPEED")
    {
        if (!d_data_sound_speed)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::getGlobalCellData()\n"
                << "Cell data of 'SOUND_SPEED' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_sound_speed;
    }
    else if (variable_key == "TEMPERATURE")
    {
        if (!d_data_temperature)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::getGlobalCellData()\n"
                << "Cell data of 'TEMPERATURE' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_temperature;
    }
    else if (variable_key == "CONVECTIVE_FLUX_X")
    {
        if (!d_data_convective_flux_x)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::getGlobalCellData()\n"
                << "Cell data of 'CONVECTIVE_FLUX_X' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_convective_flux_x;
    }
    else if (variable_key == "CONVECTIVE_FLUX_Y")
    {
        if (!d_data_convective_flux_y)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::getGlobalCellData()\n"
                << "Cell data of 'CONVECTIVE_FLUX_Y' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_convective_flux_y;
    }
    else if (variable_key == "CONVECTIVE_FLUX_Z")
    {
        if (!d_data_convective_flux_z)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::getGlobalCellData()\n"
                << "Cell data of 'CONVECTIVE_FLUX_Z' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_convective_flux_z;
    }
    else if (variable_key == "MAX_WAVE_SPEED_X")
    {
        if (!d_data_max_wave_speed_x)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::getGlobalCellData()\n"
                << "Cell data of 'MAX_WAVE_SPEED_X' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_max_wave_speed_x;
    }
    else if (variable_key == "MAX_WAVE_SPEED_Y")
    {
        if (!d_data_max_wave_speed_y)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::getGlobalCellData()\n"
                << "Cell data of 'MAX_WAVE_SPEED_Y' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_max_wave_speed_y;
    }
    else if (variable_key == "MAX_WAVE_SPEED_Z")
    {
        if (!d_data_max_wave_speed_z)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::getGlobalCellData()\n"
                << "Cell data of 'MAX_WAVE_SPEED_Z' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_max_wave_speed_z;
    }
    else if (variable_key == "MAX_DIFFUSIVITY")
    {
        if (!d_data_max_diffusivity)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::getGlobalCellData()\n"
                << "Cell data of 'MAX_DIFFUSIVITY' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_max_diffusivity;
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::getGlobalCellData()\n"
            << "Unknown cell data with variable_key = '" << variable_key
            << "' requested."
            << std::endl);
    }
    
    return cell_data;
}


/*
 * Get the global cell data of different cell variables in the registered patch.
 */
std::vector<boost::shared_ptr<pdat::CellData<double> > >
FlowModelFourEqnConservative::getGlobalCellData(
    const std::vector<std::string>& variable_keys)
{
    std::vector<boost::shared_ptr<pdat::CellData<double> > > cell_data(
        static_cast<int>(variable_keys.size()));
    
    for (int vi = 0; static_cast<int>(variable_keys.size()); vi++)
    {
        cell_data[vi] = getGlobalCellData(variable_keys[vi]);
    }
    
    return cell_data;
}


/*
 * Fill the interior global cell data of conservative variables with zeros.
 */
void
FlowModelFourEqnConservative::fillZeroGlobalCellDataConservativeVariables()
{
    // Check whether a patch is already registered.
    if (!d_patch)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::fillZeroGlobalCellDataConservativeVariables()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    boost::shared_ptr<pdat::CellData<double> > data_partial_densities = getGlobalCellDataPartialDensities();
    boost::shared_ptr<pdat::CellData<double> > data_momentum = getGlobalCellDataMomentum();
    boost::shared_ptr<pdat::CellData<double> > data_total_energy = getGlobalCellDataTotalEnergy();
    
    data_partial_densities->fillAll(double(0), d_interior_box);
    data_momentum->fillAll(double(0), d_interior_box);
    data_total_energy->fillAll(double(0), d_interior_box);
}


/*
 * Update the interior global cell data of conservative variables.
 */
void
FlowModelFourEqnConservative::updateGlobalCellDataConservativeVariables()
{
    // Check whether a patch is already registered.
    if (!d_patch)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::updateGlobalCellDataConservativeVariables()\n"
            << "No patch is registered yet."
            << std::endl);
    }
}


/*
 * Get the global cell data of the conservative variables in the registered patch.
 */
std::vector<boost::shared_ptr<pdat::CellData<double> > >
FlowModelFourEqnConservative::getGlobalCellDataConservativeVariables()
{
    // Check whether a patch is already registered.
    if (!d_patch)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::getGlobalCellDataConservativeVariables()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    std::vector<boost::shared_ptr<pdat::CellData<double> > > global_cell_data;
    global_cell_data.reserve(3);
    
    global_cell_data.push_back(getGlobalCellDataPartialDensities());
    global_cell_data.push_back(getGlobalCellDataMomentum());
    global_cell_data.push_back(getGlobalCellDataTotalEnergy());
    
    return global_cell_data;
}


/*
 * Get the global cell data of the primitive variables in the registered patch.
 */
std::vector<boost::shared_ptr<pdat::CellData<double> > >
FlowModelFourEqnConservative::getGlobalCellDataPrimitiveVariables()
{
    // Check whether a patch is already registered.
    if (!d_patch)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::getGlobalCellDataPrimitiveVariables()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    std::vector<boost::shared_ptr<pdat::CellData<double> > > global_cell_data;
    global_cell_data.reserve(3);
    
    global_cell_data.push_back(getGlobalCellDataPartialDensities());
    if (!d_data_velocity)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::getGlobalCellDataPrimitiveVariables()\n"
            << "Cell data of 'VELOCITY' is not registered/computed yet."
            << std::endl);
    }
    global_cell_data.push_back(d_data_velocity);
    if (!d_data_pressure)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::getGlobalCellDataPrimitiveVariables()\n"
            << "Cell data of 'PRESSURE' is not registered/computed yet."
            << std::endl);
    }
    global_cell_data.push_back(d_data_pressure);
    
    return global_cell_data;
}


/*
 * Get the number of projection variables for transformation between conservative
 * variables and characteristic variables.
 */
int
FlowModelFourEqnConservative::getNumberOfProjectionVariablesForConservativeVariables() const
{
    TBOX_ERROR(d_object_name
        << ": FlowModelFourEqnConservative::"
        << "getNumberOfProjectionVariablesForConservativeVariables()\n"
        << "Method getNumberOfProjectionVariablesForConservativeVariables()"
        << " is not yet implemented."
        << std::endl);
    
    return 0;
}

/*
 * Get the number of projection variables for transformation between primitive variables
 * and characteristic variables.
 */
int
FlowModelFourEqnConservative::getNumberOfProjectionVariablesForPrimitiveVariables() const
{
    return d_num_species + 2;
}


/*
 * Compute global side data of the projection variables for transformation between
 * conservative variables and characteristic variables.
 */
void
FlowModelFourEqnConservative::computeGlobalSideDataProjectionVariablesForConservativeVariables(
    std::vector<boost::shared_ptr<pdat::SideData<double> > >& projection_variables)
{
    NULL_USE(projection_variables);
    
    TBOX_ERROR(d_object_name
        << ": FlowModelFourEqnConservative::"
        << "computeGlobalSideDataProjectionVariablesForConservativeVariables()\n"
        << "Method computeGlobalSideDataProjectionVariablesForConservativeVariables()"
        << " is not yet implemented."
        << std::endl);
}


/*
 * Compute global side data of the projection variables for transformation between
 * primitive variables and characteristic variables.
 */
void
FlowModelFourEqnConservative::computeGlobalSideDataProjectionVariablesForPrimitiveVariables(
    std::vector<boost::shared_ptr<pdat::SideData<double> > >& projection_variables)
{
    // Create empty box.
    const hier::Box empty_box(d_dim);
    
    /*
     * Get the number of ghost cells and ghost cell dimension of projection variables.
     */
    
    const hier::IntVector num_ghosts_projection_var = projection_variables[0]->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_projection_var =
        projection_variables[0]->getGhostBox().numberCells();
    
    /*
     * Check the size of variables.
     */
    
    if (static_cast<int>(projection_variables.size()) != d_num_species + 2)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::"
            << "computeGlobalSideDataProjectionVariablesForPrimitiveVariables()\n"
            << "There should be number of species projection plus two variables."
            << std::endl);
    }
    
    /*
     * Check potential failures.
     */
    
    for (int vi = 0; vi < static_cast<int>(projection_variables.size()); vi++)
    {
        const hier::IntVector interior_dims_projection_var =
            projection_variables[vi]->getBox().numberCells();
        if (interior_dims_projection_var != d_interior_dims)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::"
                << "computeGlobalSideDataProjectionVariablesForPrimitiveVariables()\n"
                << "The interior dimension of the projection variables does not match that of patch."
                << std::endl);
        }
    }
    
    for (int vi = 1; vi < static_cast<int>(projection_variables.size()); vi++)
    {
        if (num_ghosts_projection_var != projection_variables[vi]->getGhostCellWidth())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::"
                << "computeGlobalSideDataProjectionVariablesForPrimitiveVariables()\n"
                << "The projection variables don't have same ghost cell width."
                << std::endl);
        }
    }
    
    if (num_ghosts_projection_var > d_num_ghosts)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::"
            << "computeGlobalSideDataProjectionVariablesForPrimitiveVariables()\n"
            << "The projection variables have ghost cell width larger than that of density."
            << std::endl);
    }
    
    if (num_ghosts_projection_var > d_num_subghosts_sound_speed)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::"
            << "computeGlobalSideDataProjectionVariablesForPrimitiveVariables()\n"
            << "The projection variables have ghost cell width larger than that of sound speed."
            << std::endl);
    }
    
    // Get the cell data of the variable partial densities.
    boost::shared_ptr<pdat::CellData<double> > data_partial_densities =
        getGlobalCellDataPartialDensities();
    
    // Get the pointers to the cell data of partial densities, total density and sound speed.
    std::vector<double*> rho_Y;
    rho_Y.reserve(d_num_species);
    for (int si = 0; si < d_num_species; si++)
    {
        rho_Y.push_back(data_partial_densities->getPointer(si));
    }
    if (!d_data_density)
    {
        computeGlobalCellDataDensity(empty_box);
    }
    double* rho = d_data_density->getPointer(0);
    if (!d_data_sound_speed)
    {
        computeGlobalCellDataSoundSpeedWithDensityMassFractionsAndPressure(empty_box);
    }
    double* c = d_data_sound_speed->getPointer(0);
    
    /*
     * Declare pointers to different data.
     */
    
    std::vector<double*> rho_Y_average;
    rho_Y_average.resize(d_num_species);
    double* rho_average = nullptr;
    double* c_average = nullptr;
    
    if (d_dim == tbox::Dimension(1))
    {
        const int interior_dim_0 = d_interior_dims[0];
        
        const int num_ghosts_0 = d_num_ghosts[0];
        const int num_ghosts_0_projection_var = num_ghosts_projection_var[0];
        const int num_subghosts_0_density = d_num_subghosts_density[0];
        const int num_subghosts_0_sound_speed = d_num_subghosts_sound_speed[0];
        
        switch (d_proj_var_primitive_averaging)
        {
            case AVERAGING::SIMPLE:
            {
                /*
                 * Compute the projection variables in the x-direction.
                 */
                
                for (int si = 0; si < d_num_species; si++)
                {
                    rho_Y_average[si] = projection_variables[si]->getPointer(0);
                }
                rho_average = projection_variables[d_num_species]->getPointer(0);
                c_average = projection_variables[d_num_species + 1]->getPointer(0);
                
                for (int si = 0; si < d_num_species; si++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = -num_ghosts_0_projection_var;
                         i < interior_dim_0 + 1 + num_ghosts_0_projection_var;
                         i++)
                    {
                        // Compute the linear indices.
                        const int idx_face_x = i + num_ghosts_0_projection_var;
                        const int idx_L = i - 1 + num_ghosts_0;
                        const int idx_R = i + num_ghosts_0;
                        
                        rho_Y_average[si][idx_face_x] = double(1)/double(2)*(rho_Y[si][idx_L] + rho_Y[si][idx_R]);
                    }
                }
                
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = -num_ghosts_0_projection_var;
                     i < interior_dim_0 + 1 + num_ghosts_0_projection_var;
                     i++)
                {
                    // Compute the linear indices.
                    const int idx_face_x = i + num_ghosts_0_projection_var;
                    const int idx_density_L = i - 1 + num_subghosts_0_density;
                    const int idx_density_R = i + num_subghosts_0_density;
                    const int idx_sound_speed_L = i - 1 + num_subghosts_0_sound_speed;
                    const int idx_sound_speed_R = i + num_subghosts_0_sound_speed;
                    
                    rho_average[idx_face_x] = double(1)/double(2)*(rho[idx_density_L] + rho[idx_density_R]);
                    c_average[idx_face_x] = double(1)/double(2)*(c[idx_sound_speed_L] + c[idx_sound_speed_R]);
                }
                
                break;
            }
            case AVERAGING::ROE:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelFourEqnConservative::"
                    << "computeGlobalSideDataProjectionVariablesForPrimitiveVariables()\n"
                    << "Roe averaging is not yet implemented."
                    << std::endl);
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelFourEqnConservative::"
                    << "computeGlobalSideDataProjectionVariablesForPrimitiveVariables()\n"
                    << "Unknown d_proj_var_primitive_averaging given."
                    << std::endl);
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int interior_dim_0 = d_interior_dims[0];
        const int interior_dim_1 = d_interior_dims[1];
        
        const int num_ghosts_0 = d_num_ghosts[0];
        const int num_ghosts_1 = d_num_ghosts[1];
        const int ghostcell_dim_0 = d_ghostcell_dims[0];
        
        const int num_ghosts_0_projection_var = num_ghosts_projection_var[0];
        const int num_ghosts_1_projection_var = num_ghosts_projection_var[1];
        const int ghostcell_dim_0_projection_var = ghostcell_dims_projection_var[0];
        
        const int num_subghosts_0_density = d_num_subghosts_density[0];
        const int num_subghosts_1_density= d_num_subghosts_density[1];
        const int subghostcell_dim_0_density = d_subghostcell_dims_density[0];
        
        const int num_subghosts_0_sound_speed = d_num_subghosts_sound_speed[0];
        const int num_subghosts_1_sound_speed = d_num_subghosts_sound_speed[1];
        const int subghostcell_dim_0_sound_speed = d_subghostcell_dims_sound_speed[0];
        
        switch (d_proj_var_primitive_averaging)
        {
            case AVERAGING::SIMPLE:
            {
                /*
                 * Compute the projection variables in the x-direction.
                 */
                
                for (int si = 0; si < d_num_species; si++)
                {
                    rho_Y_average[si] = projection_variables[si]->getPointer(0);
                }
                rho_average = projection_variables[d_num_species]->getPointer(0);
                c_average = projection_variables[d_num_species + 1]->getPointer(0);
                
                for (int si = 0; si < d_num_species; si++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = -num_ghosts_0_projection_var;
                             i < interior_dim_0 + 1 + num_ghosts_0_projection_var;
                             i++)
                        {
                            // Compute the linear indices.
                            const int idx_face_x = (i + num_ghosts_0_projection_var) +
                                (j + num_ghosts_1_projection_var)*(ghostcell_dim_0_projection_var + 1);
                            
                            const int idx_L = (i - 1 + num_ghosts_0) +
                                (j + num_ghosts_1)*ghostcell_dim_0;
                            
                            const int idx_R = (i + num_ghosts_0) +
                                (j + num_ghosts_1)*ghostcell_dim_0;
                            
                            rho_Y_average[si][idx_face_x] = double(1)/double(2)*(rho_Y[si][idx_L] + rho_Y[si][idx_R]);
                        }
                    }
                }
                
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = -num_ghosts_0_projection_var;
                         i < interior_dim_0 + 1 + num_ghosts_0_projection_var;
                         i++)
                    {
                        // Compute the linear indices.
                        const int idx_face_x = (i + num_ghosts_0_projection_var) +
                            (j + num_ghosts_1_projection_var)*(ghostcell_dim_0_projection_var + 1);
                        
                        const int idx_density_L = (i - 1 + num_subghosts_0_density) +
                            (j + num_subghosts_1_density)*subghostcell_dim_0_density;
                        
                        const int idx_density_R = (i + num_subghosts_0_density) +
                            (j + num_subghosts_1_density)*subghostcell_dim_0_density;
                        
                        const int idx_sound_speed_L = (i - 1 + num_subghosts_0_sound_speed) +
                            (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed;
                        
                        const int idx_sound_speed_R = (i + num_subghosts_0_sound_speed) +
                            (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed;
                        
                        rho_average[idx_face_x] = double(1)/double(2)*(rho[idx_density_L] + rho[idx_density_R]);
                        c_average[idx_face_x] = double(1)/double(2)*(c[idx_sound_speed_L] + c[idx_sound_speed_R]);
                    }
                }
                
                /*
                 * Compute the projection variables in the y-direction.
                 */
                
                for (int si = 0; si < d_num_species; si++)
                {
                    rho_Y_average[si] = projection_variables[si]->getPointer(1);
                }
                rho_average = projection_variables[d_num_species]->getPointer(1);
                c_average = projection_variables[d_num_species + 1]->getPointer(1);
                
                for (int si = 0; si < d_num_species; si++)
                {
                    for (int j = -num_ghosts_1_projection_var;
                         j < interior_dim_1 + 1 + num_ghosts_1_projection_var;
                         j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_face_y = (i + num_ghosts_0_projection_var) +
                                (j + num_ghosts_1_projection_var)*ghostcell_dim_0_projection_var;
                            
                            const int idx_B = (i + num_ghosts_0) +
                                (j - 1 + num_ghosts_1)*ghostcell_dim_0;
                            
                            const int idx_T = (i + num_ghosts_0) +
                                (j + num_ghosts_1)*ghostcell_dim_0;
                            
                            rho_Y_average[si][idx_face_y] = double(1)/double(2)*(rho_Y[si][idx_B] + rho_Y[si][idx_T]);
                        }
                    }
                }
                
                for (int j = -num_ghosts_1_projection_var;
                     j < interior_dim_1 + 1 + num_ghosts_1_projection_var;
                     j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_face_y = (i + num_ghosts_0_projection_var) +
                            (j + num_ghosts_1_projection_var)*ghostcell_dim_0_projection_var;
                        
                        const int idx_density_B = (i + num_subghosts_0_density) +
                            (j - 1 + num_subghosts_1_density)*subghostcell_dim_0_density;
                        
                        const int idx_density_T = (i + num_subghosts_0_density) +
                            (j + num_subghosts_1_density)*subghostcell_dim_0_density;
                        
                        const int idx_sound_speed_B = (i + num_subghosts_0_sound_speed) +
                            (j - 1 + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed;
                        
                        const int idx_sound_speed_T = (i + num_subghosts_0_sound_speed) +
                            (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed;
                        
                        rho_average[idx_face_y] = double(1)/double(2)*(rho[idx_density_B] + rho[idx_density_T]);
                        c_average[idx_face_y] = double(1)/double(2)*(c[idx_sound_speed_B] + c[idx_sound_speed_T]);
                    }
                }
                
                break;
            }
            case AVERAGING::ROE:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelFourEqnConservative::"
                    << "computeGlobalSideDataProjectionVariablesForPrimitiveVariables()\n"
                    << "Roe averaging is not yet implemented."
                    << std::endl);
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelFourEqnConservative::"
                    << "computeGlobalSideDataProjectionVariablesForPrimitiveVariables()\n"
                    << "Unknown d_proj_var_primitive_averaging given."
                    << std::endl);
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int interior_dim_0 = d_interior_dims[0];
        const int interior_dim_1 = d_interior_dims[1];
        const int interior_dim_2 = d_interior_dims[2];
        
        const int num_ghosts_0 = d_num_ghosts[0];
        const int num_ghosts_1 = d_num_ghosts[1];
        const int num_ghosts_2 = d_num_ghosts[2];
        const int ghostcell_dim_0 = d_ghostcell_dims[0];
        const int ghostcell_dim_1 = d_ghostcell_dims[1];
        
        const int num_ghosts_0_projection_var = num_ghosts_projection_var[0];
        const int num_ghosts_1_projection_var = num_ghosts_projection_var[1];
        const int num_ghosts_2_projection_var = num_ghosts_projection_var[2];
        const int ghostcell_dim_0_projection_var = ghostcell_dims_projection_var[0];
        const int ghostcell_dim_1_projection_var = ghostcell_dims_projection_var[1];
        
        const int num_subghosts_0_density = d_num_subghosts_density[0];
        const int num_subghosts_1_density = d_num_subghosts_density[1];
        const int num_subghosts_2_density = d_num_subghosts_density[2];
        const int subghostcell_dim_0_density = d_subghostcell_dims_density[0];
        const int subghostcell_dim_1_density = d_subghostcell_dims_density[1];
        
        const int num_subghosts_0_sound_speed = d_num_subghosts_sound_speed[0];
        const int num_subghosts_1_sound_speed = d_num_subghosts_sound_speed[1];
        const int num_subghosts_2_sound_speed = d_num_subghosts_sound_speed[2];
        const int subghostcell_dim_0_sound_speed = d_subghostcell_dims_sound_speed[0];
        const int subghostcell_dim_1_sound_speed = d_subghostcell_dims_sound_speed[1];
        
        switch (d_proj_var_primitive_averaging)
        {
            case AVERAGING::SIMPLE:
            {
                /*
                 * Compute the projection variables in the x-direction.
                 */
                
                for (int si = 0; si < d_num_species; si++)
                {
                    rho_Y_average[si] = projection_variables[si]->getPointer(0);
                }
                rho_average = projection_variables[d_num_species]->getPointer(0);
                c_average = projection_variables[d_num_species + 1]->getPointer(0);
                
                for (int si = 0; si < d_num_species; si++)
                {
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = -num_ghosts_0_projection_var;
                                 i < interior_dim_0 + 1 + num_ghosts_0_projection_var;
                                 i++)
                            {
                                // Compute the linear indices.
                                const int idx_face_x = (i + num_ghosts_0_projection_var) +
                                    (j + num_ghosts_1_projection_var)*(ghostcell_dim_0_projection_var + 1) +
                                    (k + num_ghosts_2_projection_var)*(ghostcell_dim_0_projection_var + 1)*
                                        ghostcell_dim_1_projection_var;
                                
                                const int idx_L = (i - 1 + num_ghosts_0) +
                                    (j + num_ghosts_1)*ghostcell_dim_0 +
                                    (k + num_ghosts_2)*ghostcell_dim_0*
                                        ghostcell_dim_1;
                                
                                const int idx_R = (i + num_ghosts_0) +
                                    (j + num_ghosts_1)*ghostcell_dim_0 +
                                    (k + num_ghosts_2)*ghostcell_dim_0*
                                        ghostcell_dim_1;
                                
                                rho_Y_average[si][idx_face_x] = double(1)/double(2)*(rho_Y[si][idx_L] + rho_Y[si][idx_R]);
                            }
                        }
                    }
                }
                
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = -num_ghosts_0_projection_var;
                             i < interior_dim_0 + 1 + num_ghosts_0_projection_var;
                             i++)
                        {
                            // Compute the linear indices.
                            const int idx_face_x = (i + num_ghosts_0_projection_var) +
                                (j + num_ghosts_1_projection_var)*(ghostcell_dim_0_projection_var + 1) +
                                (k + num_ghosts_2_projection_var)*(ghostcell_dim_0_projection_var + 1)*
                                    ghostcell_dim_1_projection_var;
                            
                            const int idx_density_L = (i - 1 + num_subghosts_0_density) +
                                (j + num_subghosts_1_density)*subghostcell_dim_0_density +
                                (k + num_subghosts_2_density)*subghostcell_dim_0_density*
                                    subghostcell_dim_1_density;
                            
                            const int idx_density_R = (i + num_subghosts_0_density) +
                                (j + num_subghosts_1_density)*subghostcell_dim_0_density +
                                (k + num_subghosts_2_density)*subghostcell_dim_0_density*
                                    subghostcell_dim_1_density;
                            
                            const int idx_sound_speed_L = (i - 1 + num_subghosts_0_sound_speed) +
                                (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed +
                                (k + num_subghosts_2_sound_speed)*subghostcell_dim_0_sound_speed*
                                    subghostcell_dim_1_sound_speed;
                            
                            const int idx_sound_speed_R = (i + num_subghosts_0_sound_speed) +
                                (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed +
                                (k + num_subghosts_2_sound_speed)*subghostcell_dim_0_sound_speed*
                                    subghostcell_dim_1_sound_speed;
                            
                            rho_average[idx_face_x] = double(1)/double(2)*(rho[idx_density_L] + rho[idx_density_R]);
                            c_average[idx_face_x] = double(1)/double(2)*(c[idx_sound_speed_L] + c[idx_sound_speed_R]);
                        }
                    }
                }
                
                /*
                 * Compute the projection variables in the y-direction.
                 */
                
                for (int si = 0; si < d_num_species; si++)
                {
                    rho_Y_average[si] = projection_variables[si]->getPointer(1);
                }
                rho_average = projection_variables[d_num_species]->getPointer(1);
                c_average = projection_variables[d_num_species + 1]->getPointer(1);
                
                for (int si = 0; si < d_num_species; si++)
                {
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = -num_ghosts_1_projection_var;
                             j < interior_dim_1 + 1 + num_ghosts_1_projection_var;
                             j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                // Compute the linear indices.
                                const int idx_face_y = (i + num_ghosts_0_projection_var) +
                                    (j + num_ghosts_1_projection_var)*ghostcell_dim_0_projection_var +
                                    (k + num_ghosts_2_projection_var)*ghostcell_dim_0_projection_var*
                                        (ghostcell_dim_1_projection_var + 1);
                                
                                const int idx_B = (i + num_ghosts_0) +
                                    (j - 1 + num_ghosts_1)*ghostcell_dim_0 +
                                    (k + num_ghosts_2)*ghostcell_dim_0*
                                        ghostcell_dim_1;
                                
                                const int idx_T = (i + num_ghosts_0) +
                                    (j + num_ghosts_1)*ghostcell_dim_0 +
                                    (k + num_ghosts_2)*ghostcell_dim_0*
                                        ghostcell_dim_1;
                                
                                rho_Y_average[si][idx_face_y] = double(1)/double(2)*(rho_Y[si][idx_B] + rho_Y[si][idx_T]);
                            }
                        }
                    }
                }
                
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = -num_ghosts_1_projection_var;
                         j < interior_dim_1 + 1 + num_ghosts_1_projection_var;
                         j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_face_y = (i + num_ghosts_0_projection_var) +
                                (j + num_ghosts_1_projection_var)*ghostcell_dim_0_projection_var +
                                (k + num_ghosts_2_projection_var)*ghostcell_dim_0_projection_var*
                                    (ghostcell_dim_1_projection_var + 1);
                            
                            const int idx_density_B = (i + num_subghosts_0_density) +
                                (j - 1 + num_subghosts_1_density)*subghostcell_dim_0_density +
                                (k + num_subghosts_2_density)*subghostcell_dim_0_density*
                                    subghostcell_dim_1_density;
                            
                            const int idx_density_T = (i + num_subghosts_0_density) +
                                (j + num_subghosts_1_density)*subghostcell_dim_0_density +
                                (k + num_subghosts_2_density)*subghostcell_dim_0_density*
                                    subghostcell_dim_1_density;
                            
                            const int idx_sound_speed_B = (i + num_subghosts_0_sound_speed) +
                                (j - 1 + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed +
                                (k + num_subghosts_2_sound_speed)*subghostcell_dim_0_sound_speed*
                                    subghostcell_dim_1_sound_speed;
                            
                            const int idx_sound_speed_T = (i + num_subghosts_0_sound_speed) +
                                (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed +
                                (k + num_subghosts_2_sound_speed)*subghostcell_dim_0_sound_speed*
                                    subghostcell_dim_1_sound_speed;
                            
                            rho_average[idx_face_y] = double(1)/double(2)*(rho[idx_density_B] + rho[idx_density_T]);
                            c_average[idx_face_y] = double(1)/double(2)*(c[idx_sound_speed_B] + c[idx_sound_speed_T]);
                        }
                    }
                }
                
                /*
                 * Compute the projection variables in the z-direction.
                 */
                
                for (int si = 0; si < d_num_species; si++)
                {
                    rho_Y_average[si] = projection_variables[si]->getPointer(2);
                }
                rho_average = projection_variables[d_num_species]->getPointer(2);
                c_average = projection_variables[d_num_species + 1]->getPointer(2);
                
                for (int si = 0; si < d_num_species; si++)
                {
                    for (int k = -num_ghosts_2_projection_var;
                         k < interior_dim_2 + 1 + num_ghosts_2_projection_var;
                         k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                // Compute the linear indices.
                                const int idx_face_z = (i + num_ghosts_0_projection_var) +
                                    (j + num_ghosts_1_projection_var)*ghostcell_dim_0_projection_var +
                                    (k + num_ghosts_2_projection_var)*ghostcell_dim_0_projection_var*
                                        ghostcell_dim_1_projection_var;
                                
                                const int idx_B = (i + num_ghosts_0) +
                                    (j + num_ghosts_1)*ghostcell_dim_0 +
                                    (k - 1 + num_ghosts_2)*ghostcell_dim_0*
                                        ghostcell_dim_1;
                                
                                const int idx_F = (i + num_ghosts_0) +
                                    (j + num_ghosts_1)*ghostcell_dim_0 +
                                    (k + num_ghosts_2)*ghostcell_dim_0*
                                        ghostcell_dim_1;
                                
                                rho_Y_average[si][idx_face_z] = double(1)/double(2)*(rho_Y[si][idx_B] + rho_Y[si][idx_F]);
                            }
                        }
                    }
                }
                
                for (int k = -num_ghosts_2_projection_var;
                     k < interior_dim_2 + 1 + num_ghosts_2_projection_var;
                     k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_face_z = (i + num_ghosts_0_projection_var) +
                                (j + num_ghosts_1_projection_var)*ghostcell_dim_0_projection_var +
                                (k + num_ghosts_2_projection_var)*ghostcell_dim_0_projection_var*
                                    ghostcell_dim_1_projection_var;
                            
                            const int idx_density_B = (i + num_subghosts_0_density) +
                                (j + num_subghosts_1_density)*subghostcell_dim_0_density +
                                (k - 1 + num_subghosts_2_density)*subghostcell_dim_0_density*
                                    subghostcell_dim_1_density;
                            
                            const int idx_density_F = (i + num_subghosts_0_density) +
                                (j + num_subghosts_1_density)*subghostcell_dim_0_density +
                                (k + num_subghosts_2_density)*subghostcell_dim_0_density*
                                    subghostcell_dim_1_density;
                            
                            const int idx_sound_speed_B = (i + num_subghosts_0_sound_speed) +
                                (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed +
                                (k - 1 + num_subghosts_2_sound_speed)*subghostcell_dim_0_sound_speed*
                                    subghostcell_dim_1_sound_speed;
                            
                            const int idx_sound_speed_F = (i + num_subghosts_0_sound_speed) +
                                (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed +
                                (k + num_subghosts_2_sound_speed)*subghostcell_dim_0_sound_speed*
                                    subghostcell_dim_1_sound_speed;
                            
                            rho_average[idx_face_z] = double(1)/double(2)*(rho[idx_density_B] + rho[idx_density_F]);
                            c_average[idx_face_z] = double(1)/double(2)*(c[idx_sound_speed_B] + c[idx_sound_speed_F]);
                        }
                    }
                }
                
                break;
            }
            case AVERAGING::ROE:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelFourEqnConservative::"
                    << "computeGlobalSideDataProjectionVariablesForPrimitiveVariables()\n"
                    << "Roe averaging is not yet implemented."
                    << std::endl);
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelFourEqnConservative::"
                    << "computeGlobalSideDataProjectionVariablesForPrimitiveVariables()\n"
                    << "Unknown d_proj_var_primitive_averaging given."
                    << std::endl);
            }
        }
    }
}


/*
 * Compute global side data of characteristic variables from conservative variables.
 */
void
FlowModelFourEqnConservative::computeGlobalSideDataCharacteristicVariablesFromConservativeVariables(
    std::vector<boost::shared_ptr<pdat::SideData<double> > >& characteristic_variables,
    const std::vector<boost::shared_ptr<pdat::CellData<double> > >& conservative_variables,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& projection_variables,
    const int& idx_offset)
{
    NULL_USE(characteristic_variables);
    NULL_USE(conservative_variables);
    NULL_USE(projection_variables);
    NULL_USE(idx_offset);
    
    TBOX_ERROR(d_object_name
        << ": FlowModelFourEqnConservative::"
        << "computeGlobalSideDataCharacteristicVariablesFromConservativeVariables()\n"
        << "Method computeGlobalSideDataCharacteristicVariablesFromConservativeVariables()"
        << " is not yet implemented."
        << std::endl);
}


/*
 * Compute global side data of characteristic variables from primitive variables.
 */
void
FlowModelFourEqnConservative::computeGlobalSideDataCharacteristicVariablesFromPrimitiveVariables(
    std::vector<boost::shared_ptr<pdat::SideData<double> > >& characteristic_variables,
    const std::vector<boost::shared_ptr<pdat::CellData<double> > >& primitive_variables,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& projection_variables,
    const int& idx_offset)
{
    /*
     * Get the numbers of ghost cells of the variables.
     */
    
    const hier::IntVector num_ghosts_characteristic_var = characteristic_variables[0]->
        getGhostCellWidth();
    
    std::vector<hier::IntVector> num_ghosts_primitive_var;
    num_ghosts_primitive_var.reserve(static_cast<int>(primitive_variables.size()));
    for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
    {
        num_ghosts_primitive_var.push_back(primitive_variables[vi]->
            getGhostCellWidth());
    }
    
    const hier::IntVector num_ghosts_projection_var = projection_variables[0]->getGhostCellWidth();
    
    /*
     * Get the ghost cell dimensions of characteristic and primitive variables.
     */
    
    const hier::IntVector ghostcell_dims_characteristic_var = characteristic_variables[0]->
        getGhostBox().numberCells();
    
    std::vector<hier::IntVector> ghostcell_dims_primitive_var;
    ghostcell_dims_primitive_var.reserve(static_cast<int>(primitive_variables.size()));
    for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
    {
        ghostcell_dims_primitive_var.push_back(primitive_variables[vi]->
            getGhostBox().numberCells());
    }
    
    /*
     * Check the size of variables.
     */
    
    if (static_cast<int>(characteristic_variables.size()) != d_num_eqn)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::"
            << "computeGlobalSideDataCharacteristicVariablesFromPrimitiveVariables()\n"
            << "The number of characteristic variables are incorrect."
            << std::endl);
    }
    if (static_cast<int>(primitive_variables.size()) != 3)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::"
            << "computeGlobalSideDataCharacteristicVariablesFromPrimitiveVariables()\n"
            << "The number of primitive variables are incorrect."
            << std::endl);
    }
    if (primitive_variables[0]->getDepth() != d_num_species ||
        primitive_variables[1]->getDepth() != d_dim.getValue() ||
        primitive_variables[2]->getDepth() != 1)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::"
            << "computeGlobalSideDataCharacteristicVariablesFromPrimitiveVariables()\n"
            << "The depths of one or more primitive variables are incorrect."
            << std::endl);
    }
    if (static_cast<int>(projection_variables.size()) != d_num_species + 2)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::"
            << "computeGlobalSideDataCharacteristicVariablesFromPrimitiveVariables()\n"
            << "There should be number of species projection plus two variables."
            << std::endl);
    }
    
    /*
     * Check potential failures.
     */
    
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        const hier::IntVector interior_dims_characteristic_var =
            characteristic_variables[ei]->getBox().numberCells();
        
        if (interior_dims_characteristic_var != d_interior_dims)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::"
                << "computeGlobalSideDataCharacteristicVariablesFromPrimitiveVariables()\n"
                << "The interior dimension of the characteristic variables does not match that of patch."
                << std::endl);
        }
    }
    for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
    {
        const hier::IntVector interior_dims_primitive_var =
            primitive_variables[vi]->getBox().numberCells();
        
        if (interior_dims_primitive_var != d_interior_dims)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::"
                << "computeGlobalSideDataCharacteristicVariablesFromPrimitiveVariables()\n"
                << "The interior dimension of the primitive variables does not match that of patch."
                << std::endl);
        }
    }
    for (int vi = 0; vi < d_num_species + 2; vi++)
    {
        const hier::IntVector interior_dims_projection_var = projection_variables[vi]->getBox().numberCells();
        if (interior_dims_projection_var != d_interior_dims)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::"
                << "computeGlobalSideDataCharacteristicVariablesFromPrimitiveVariables()\n"
                << "The interior dimension of the projection variables does not match that of patch."
                << std::endl);
        }
    }
    
    for (int ei = 1; ei < d_num_eqn; ei++)
    {
        if (num_ghosts_characteristic_var != characteristic_variables[ei]->getGhostCellWidth())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::"
                << "computeGlobalSideDataCharacteristicVariablesFromPrimitiveVariables()\n"
                << "The characteristic variables don't have same ghost cell width."
                << std::endl);
        }
    }
    for (int vi = 1; vi < d_num_species + 2; vi++)
    {
        if (num_ghosts_projection_var != projection_variables[vi]->getGhostCellWidth())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::"
                << "computeGlobalSideDataCharacteristicVariablesFromPrimitiveVariables()\n"
                << "The projection variables don't have same ghost cell width."
                << std::endl);
        }
    }
    
    if (num_ghosts_projection_var != num_ghosts_characteristic_var)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::"
            << "computeGlobalSideDataCharacteristicVariablesFromPrimitiveVariables()\n"
            << "The ghost cell width of the projection variables does not match that of"
            << " characteristic variables."
            << std::endl);
    }
    
    for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
    {
        if (num_ghosts_primitive_var[vi] - num_ghosts_characteristic_var
                + (hier::IntVector::getOne(d_dim))*idx_offset < hier::IntVector::getZero(d_dim) ||
            num_ghosts_characteristic_var - num_ghosts_primitive_var[vi]
                + (hier::IntVector::getOne(d_dim))*(idx_offset + 1) > hier::IntVector::getZero(d_dim))
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::"
                << "computeGlobalSideDataCharacteristicVariablesFromPrimitiveVariables()\n"
                << "The offset index is too large or the number of ghost of characteristic variable"
                << " is too large."
                << std::endl);
        }
    }
    
    /*
     * Declare containers to store pointers to different data.
     */
    
    std::vector<double*> V;
    V.reserve(d_num_eqn);
    
    for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
    {
        int depth = primitive_variables[vi]->getDepth();
        
        for (int di = 0; di < depth; di++)
        {
            V.push_back(primitive_variables[vi]->getPointer(di));
        }
    }
    
    std::vector<double*> W;
    W.resize(d_num_eqn);
    
    std::vector<double*> rho_Y_average;
    rho_Y_average.resize(d_num_species);
    double* rho_average = nullptr;
    double* c_average = nullptr;
    
    if (d_dim == tbox::Dimension(1))
    {
        const int interior_dim_0 = d_interior_dims[0];
        
        const int num_ghosts_0_characteristic_var = num_ghosts_characteristic_var[0];
        const int num_ghosts_0_rho_Y = num_ghosts_primitive_var[0][0];
        const int num_ghosts_0_vel = num_ghosts_primitive_var[1][0];
        const int num_ghosts_0_p = num_ghosts_primitive_var[2][0];
        
        const int idx_offset_rho_Y = idx_offset;
        const int idx_offset_vel = idx_offset;
        const int idx_offset_p = idx_offset;
        
        /*
         * Compute the characteristic variables in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(0);
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y_average[si] = projection_variables[si]->getPointer(0);
        }
        rho_average = projection_variables[d_num_species]->getPointer(0);
        c_average = projection_variables[d_num_species + 1]->getPointer(0);
        
        for (int si = 0; si < d_num_species; si++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = -num_ghosts_0_characteristic_var;
                 i < interior_dim_0 + 1 + num_ghosts_0_characteristic_var;
                 i++)
            {
                // Compute the linear indices.
                const int idx_face = i + num_ghosts_0_characteristic_var;
                const int idx_rho_Y = i + idx_offset_rho_Y + num_ghosts_0_rho_Y;
                const int idx_p = i + idx_offset_p + num_ghosts_0_p;
                
                W[1 + si][idx_face] = V[si][idx_rho_Y] - rho_Y_average[si][idx_face]/
                    (rho_average[idx_face]*c_average[idx_face]*
                        c_average[idx_face])*V[d_num_species + 1][idx_p];
            }
        }
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = -num_ghosts_0_characteristic_var;
             i < interior_dim_0 + 1 + num_ghosts_0_characteristic_var;
             i++)
        {
            // Compute the linear indices.
            const int idx_face = i + num_ghosts_0_characteristic_var;
            const int idx_vel = i + idx_offset_vel + num_ghosts_0_vel;
            const int idx_p = i + idx_offset_p + num_ghosts_0_p;
            
            W[0][idx_face] = V[d_num_species][idx_vel] -
                double(1)/(rho_average[idx_face]*c_average[idx_face])*V[d_num_species + 1][idx_p];
            
            W[d_num_species + 1][idx_face] = V[d_num_species][idx_vel] +
                double(1)/(rho_average[idx_face]*c_average[idx_face])*V[d_num_species + 1][idx_p];
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int interior_dim_0 = d_interior_dims[0];
        const int interior_dim_1 = d_interior_dims[1];
        
        const int num_ghosts_0_characteristic_var = num_ghosts_characteristic_var[0];
        const int num_ghosts_1_characteristic_var = num_ghosts_characteristic_var[1];
        const int ghostcell_dim_0_characteristic_var = ghostcell_dims_characteristic_var[0];
        
        const int num_ghosts_0_rho_Y = num_ghosts_primitive_var[0][0];
        const int num_ghosts_1_rho_Y = num_ghosts_primitive_var[0][1];
        const int ghostcell_dim_0_rho_Y = ghostcell_dims_primitive_var[0][0];
        
        const int num_ghosts_0_vel = num_ghosts_primitive_var[1][0];
        const int num_ghosts_1_vel = num_ghosts_primitive_var[1][1];
        const int ghostcell_dim_0_vel = ghostcell_dims_primitive_var[1][0];
        
        const int num_ghosts_0_p = num_ghosts_primitive_var[2][0];
        const int num_ghosts_1_p = num_ghosts_primitive_var[2][1];
        const int ghostcell_dim_0_p = ghostcell_dims_primitive_var[2][0];
        
        const int idx_offset_rho_Y = idx_offset;
        const int idx_offset_vel = idx_offset;
        const int idx_offset_p = idx_offset;
        
        /*
         * Compute the characteristic variables in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(0);
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y_average[si] = projection_variables[si]->getPointer(0);
        }
        rho_average = projection_variables[d_num_species]->getPointer(0);
        c_average = projection_variables[d_num_species + 1]->getPointer(0);
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = -num_ghosts_0_characteristic_var;
                     i < interior_dim_0 + 1 + num_ghosts_0_characteristic_var;
                     i++)
                {
                    // Compute the linear indices.
                    const int idx_face = (i + num_ghosts_0_characteristic_var) +
                        (j + num_ghosts_1_characteristic_var)*(ghostcell_dim_0_characteristic_var + 1);
                    
                    const int idx_rho_Y = (i + idx_offset_rho_Y + num_ghosts_0_rho_Y) +
                        (j + num_ghosts_1_rho_Y)*ghostcell_dim_0_rho_Y;
                    
                    const int idx_p = (i + idx_offset_p + num_ghosts_0_p) +
                        (j + num_ghosts_1_p)*ghostcell_dim_0_p;
                    
                    W[1 + si][idx_face] = V[si][idx_rho_Y] - rho_Y_average[si][idx_face]/
                        (rho_average[idx_face]*c_average[idx_face]*c_average[idx_face])*
                            V[d_num_species + 2][idx_p];
                }
            }
        }
        
        for (int j = 0; j < interior_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = -num_ghosts_0_characteristic_var;
                 i < interior_dim_0 + 1 + num_ghosts_0_characteristic_var;
                 i++)
            {
                // Compute the linear indices.
                const int idx_face = (i + num_ghosts_0_characteristic_var) +
                    (j + num_ghosts_1_characteristic_var)*(ghostcell_dim_0_characteristic_var + 1);
                
                const int idx_vel = (i + idx_offset_vel + num_ghosts_0_vel) +
                    (j + num_ghosts_1_vel)*ghostcell_dim_0_vel;
                
                const int idx_p = (i + idx_offset_p + num_ghosts_0_p) +
                    (j + num_ghosts_1_p)*ghostcell_dim_0_p;
                
                W[0][idx_face] = V[d_num_species][idx_vel] -
                    double(1)/(rho_average[idx_face]*c_average[idx_face])*V[d_num_species + 2][idx_p];
                
                W[d_num_species + 1][idx_face] = V[d_num_species + 1][idx_vel];
                
                W[d_num_species + 2][idx_face] = V[d_num_species][idx_vel] +
                    double(1)/(rho_average[idx_face]*c_average[idx_face])*V[d_num_species + 2][idx_p];
            }
        }
        
        /*
         * Compute the characteristic variables in the y-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(1);
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y_average[si] = projection_variables[si]->getPointer(1);
        }
        rho_average = projection_variables[d_num_species]->getPointer(1);
        c_average = projection_variables[d_num_species + 1]->getPointer(1);
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = -num_ghosts_1_characteristic_var;
                 j < interior_dim_1 + 1 + num_ghosts_1_characteristic_var;
                 j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_face = (i + num_ghosts_0_characteristic_var) +
                        (j + num_ghosts_1_characteristic_var)*ghostcell_dim_0_characteristic_var;
                    
                    const int idx_rho_Y = (i + num_ghosts_0_rho_Y) +
                        (j + idx_offset_rho_Y + num_ghosts_1_rho_Y)*ghostcell_dim_0_rho_Y;
                    
                    const int idx_p = (i + num_ghosts_0_p) +
                        (j + idx_offset_p + num_ghosts_1_p)*ghostcell_dim_0_p;
                    
                    W[1 + si][idx_face] = V[si][idx_rho_Y] - rho_Y_average[si][idx_face]/
                        (rho_average[idx_face]*c_average[idx_face]*c_average[idx_face])
                            *V[d_num_species + 2][idx_p];
                }
            }
        }
        
        for (int j = -num_ghosts_1_characteristic_var;
             j < interior_dim_1 + 1 + num_ghosts_1_characteristic_var;
             j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_face = (i + num_ghosts_0_characteristic_var) +
                    (j + num_ghosts_1_characteristic_var)*ghostcell_dim_0_characteristic_var;
                
                const int idx_vel = (i + num_ghosts_0_vel) +
                    (j + idx_offset_vel + num_ghosts_1_vel)*ghostcell_dim_0_vel;
                
                const int idx_p = (i + num_ghosts_0_p) +
                    (j + idx_offset_p + num_ghosts_1_p)*ghostcell_dim_0_p;
                
                W[0][idx_face] = V[d_num_species + 1][idx_vel] -
                    double(1)/(rho_average[idx_face]*c_average[idx_face])*V[d_num_species + 2][idx_p];
                
                W[d_num_species + 1][idx_face] = V[d_num_species][idx_vel];
                
                W[d_num_species + 2][idx_face] = V[d_num_species + 1][idx_vel] +
                    double(1)/(rho_average[idx_face]*c_average[idx_face])*V[d_num_species + 2][idx_p];
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int interior_dim_0 = d_interior_dims[0];
        const int interior_dim_1 = d_interior_dims[1];
        const int interior_dim_2 = d_interior_dims[2];
        
        const int num_ghosts_0_characteristic_var = num_ghosts_characteristic_var[0];
        const int num_ghosts_1_characteristic_var = num_ghosts_characteristic_var[1];
        const int num_ghosts_2_characteristic_var = num_ghosts_characteristic_var[2];
        const int ghostcell_dim_0_characteristic_var = ghostcell_dims_characteristic_var[0];
        const int ghostcell_dim_1_characteristic_var = ghostcell_dims_characteristic_var[1];
        
        const int num_ghosts_0_rho_Y = num_ghosts_primitive_var[0][0];
        const int num_ghosts_1_rho_Y = num_ghosts_primitive_var[0][1];
        const int num_ghosts_2_rho_Y = num_ghosts_primitive_var[0][2];
        const int ghostcell_dim_0_rho_Y = ghostcell_dims_primitive_var[0][0];
        const int ghostcell_dim_1_rho_Y = ghostcell_dims_primitive_var[0][1];
        
        const int num_ghosts_0_vel = num_ghosts_primitive_var[1][0];
        const int num_ghosts_1_vel = num_ghosts_primitive_var[1][1];
        const int num_ghosts_2_vel = num_ghosts_primitive_var[1][2];
        const int ghostcell_dim_0_vel = ghostcell_dims_primitive_var[1][0];
        const int ghostcell_dim_1_vel = ghostcell_dims_primitive_var[1][1];
        
        const int num_ghosts_0_p = num_ghosts_primitive_var[2][0];
        const int num_ghosts_1_p = num_ghosts_primitive_var[2][1];
        const int num_ghosts_2_p = num_ghosts_primitive_var[2][2];
        const int ghostcell_dim_0_p = ghostcell_dims_primitive_var[2][0];
        const int ghostcell_dim_1_p = ghostcell_dims_primitive_var[2][1];
        
        const int idx_offset_rho_Y = idx_offset;
        const int idx_offset_vel = idx_offset;
        const int idx_offset_p = idx_offset;
        
        /*
         * Compute the characteristic variables in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(0);
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y_average[si] = projection_variables[si]->getPointer(0);
        }
        rho_average = projection_variables[d_num_species]->getPointer(0);
        c_average = projection_variables[d_num_species + 1]->getPointer(0);
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = -num_ghosts_0_characteristic_var;
                         i < interior_dim_0 + 1 + num_ghosts_0_characteristic_var;
                         i++)
                    {
                        // Compute the linear indices.
                        const int idx_face = (i + num_ghosts_0_characteristic_var) +
                            (j + num_ghosts_1_characteristic_var)*(ghostcell_dim_0_characteristic_var + 1) +
                            (k + num_ghosts_2_characteristic_var)*(ghostcell_dim_0_characteristic_var + 1)*
                                ghostcell_dim_1_characteristic_var;
                        
                        const int idx_rho_Y = (i + idx_offset_rho_Y + num_ghosts_0_rho_Y) +
                            (j + num_ghosts_1_rho_Y)*ghostcell_dim_0_rho_Y +
                            (k + num_ghosts_2_rho_Y)*ghostcell_dim_0_rho_Y*
                                ghostcell_dim_1_rho_Y;
                        
                        const int idx_p = (i + idx_offset_p + num_ghosts_0_p) +
                            (j + num_ghosts_1_p)*ghostcell_dim_0_p +
                            (k + num_ghosts_2_p)*ghostcell_dim_0_p*
                                ghostcell_dim_1_p;
                        
                        W[1 + si][idx_face] = V[si][idx_rho_Y] - rho_Y_average[si][idx_face]/
                            (rho_average[idx_face]*c_average[idx_face]*c_average[idx_face])*
                                V[d_num_species + 3][idx_p];
                    }
                }
            }
        }
        
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = -num_ghosts_0_characteristic_var;
                     i < interior_dim_0 + 1 + num_ghosts_0_characteristic_var;
                     i++)
                {
                    // Compute the linear indices.
                    const int idx_face = (i + num_ghosts_0_characteristic_var) +
                        (j + num_ghosts_1_characteristic_var)*(ghostcell_dim_0_characteristic_var + 1) +
                        (k + num_ghosts_2_characteristic_var)*(ghostcell_dim_0_characteristic_var + 1)*
                            ghostcell_dim_1_characteristic_var;
                    
                    const int idx_vel = (i + idx_offset_vel + num_ghosts_0_vel) +
                        (j + num_ghosts_1_vel)*ghostcell_dim_0_vel +
                        (k + num_ghosts_2_vel)*ghostcell_dim_0_vel*
                            ghostcell_dim_1_vel;
                    
                    const int idx_p = (i + idx_offset_p + num_ghosts_0_p) +
                        (j + num_ghosts_1_p)*ghostcell_dim_0_p +
                        (k + num_ghosts_2_p)*ghostcell_dim_0_p*
                            ghostcell_dim_1_p;
                    
                    W[0][idx_face] = V[d_num_species][idx_vel] -
                        double(1)/(rho_average[idx_face]*c_average[idx_face])*V[d_num_species + 3][idx_p];
                    
                    W[d_num_species + 1][idx_face] = V[d_num_species + 1][idx_vel];
                    
                    W[d_num_species + 2][idx_face] = V[d_num_species + 2][idx_vel];
                    
                    W[d_num_species + 3][idx_face] = V[d_num_species][idx_vel] +
                        double(1)/(rho_average[idx_face]*c_average[idx_face])*V[d_num_species + 3][idx_p];
                }
            }
        }
        
        /*
         * Compute the characteristic variables in the y-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(1);
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y_average[si] = projection_variables[si]->getPointer(1);
        }
        rho_average = projection_variables[d_num_species]->getPointer(1);
        c_average = projection_variables[d_num_species + 1]->getPointer(1);
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = -num_ghosts_1_characteristic_var;
                     j < interior_dim_1 + 1 + num_ghosts_1_characteristic_var;
                     j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_face = (i + num_ghosts_0_characteristic_var) +
                            (j + num_ghosts_1_characteristic_var)*ghostcell_dim_0_characteristic_var +
                            (k + num_ghosts_2_characteristic_var)*ghostcell_dim_0_characteristic_var*
                                (ghostcell_dim_1_characteristic_var + 1);
                        
                        const int idx_rho_Y = (i + num_ghosts_0_rho_Y) +
                            (j + idx_offset_rho_Y + num_ghosts_1_rho_Y)*ghostcell_dim_0_rho_Y +
                            (k + num_ghosts_2_rho_Y)*ghostcell_dim_0_rho_Y*
                                ghostcell_dim_1_rho_Y;
                        
                        const int idx_p = (i + num_ghosts_0_p) +
                            (j + idx_offset_p + num_ghosts_1_p)*ghostcell_dim_0_p +
                            (k + num_ghosts_2_p)*ghostcell_dim_0_p*
                                ghostcell_dim_1_p;
                        
                        W[1 + si][idx_face] = V[si][idx_rho_Y] - rho_Y_average[si][idx_face]/
                            (rho_average[idx_face]*c_average[idx_face]*c_average[idx_face])*
                                V[d_num_species + 3][idx_p];
                    }
                }
            }
        }
        
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = -num_ghosts_1_characteristic_var;
                 j < interior_dim_1 + 1 + num_ghosts_1_characteristic_var;
                 j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_face = (i + num_ghosts_0_characteristic_var) +
                        (j + num_ghosts_1_characteristic_var)*ghostcell_dim_0_characteristic_var +
                        (k + num_ghosts_2_characteristic_var)*ghostcell_dim_0_characteristic_var*
                            (ghostcell_dim_1_characteristic_var + 1);
                    
                    const int idx_vel = (i + num_ghosts_0_vel) +
                        (j + idx_offset_vel + num_ghosts_1_vel)*ghostcell_dim_0_vel +
                        (k + num_ghosts_2_vel)*ghostcell_dim_0_vel*
                            ghostcell_dim_1_vel;
                    
                    const int idx_p = (i + num_ghosts_0_p) +
                        (j + idx_offset_p + num_ghosts_1_p)*ghostcell_dim_0_p +
                        (k + num_ghosts_2_p)*ghostcell_dim_0_p*
                            ghostcell_dim_1_p;
                    
                    W[0][idx_face] = V[d_num_species + 1][idx_vel] -
                        double(1)/(rho_average[idx_face]*c_average[idx_face])*V[d_num_species + 3][idx_p];
                    
                    W[d_num_species + 1][idx_face] = V[d_num_species][idx_vel];
                    
                    W[d_num_species + 2][idx_face] = V[d_num_species + 2][idx_vel];
                    
                    W[d_num_species + 3][idx_face] = V[d_num_species + 1][idx_vel] +
                        double(1)/(rho_average[idx_face]*c_average[idx_face])*V[d_num_species + 3][idx_p];
                }
            }
        }
        
        /*
         * Compute the characteristic variables in the z-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(2);
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y_average[si] = projection_variables[si]->getPointer(2);
        }
        rho_average = projection_variables[d_num_species]->getPointer(2);
        c_average = projection_variables[d_num_species + 1]->getPointer(2);
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = -num_ghosts_2_characteristic_var;
                 k < interior_dim_2 + 1 + num_ghosts_2_characteristic_var;
                 k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_face = (i + num_ghosts_0_characteristic_var) +
                            (j + num_ghosts_1_characteristic_var)*ghostcell_dim_0_characteristic_var +
                            (k + num_ghosts_2_characteristic_var)*ghostcell_dim_0_characteristic_var*
                                ghostcell_dim_1_characteristic_var;
                        
                        const int idx_rho_Y = (i + num_ghosts_0_rho_Y) +
                            (j + num_ghosts_1_rho_Y)*ghostcell_dim_0_rho_Y +
                            (k + idx_offset_rho_Y + num_ghosts_2_rho_Y)*ghostcell_dim_0_rho_Y*
                                ghostcell_dim_1_rho_Y;
                        
                        const int idx_p = (i + num_ghosts_0_p) +
                            (j + num_ghosts_1_p)*ghostcell_dim_0_p +
                            (k + idx_offset_p + num_ghosts_2_p)*ghostcell_dim_0_p*
                                ghostcell_dim_1_p;
                        
                        W[1 + si][idx_face] = V[si][idx_rho_Y] - rho_Y_average[si][idx_face]/
                            (rho_average[idx_face]*c_average[idx_face]*c_average[idx_face])*
                                V[d_num_species + 3][idx_p];
                    }
                }
            }
        }
        
        for (int k = -num_ghosts_2_characteristic_var;
             k < interior_dim_2 + 1 + num_ghosts_2_characteristic_var;
             k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_face = (i + num_ghosts_0_characteristic_var) +
                        (j + num_ghosts_1_characteristic_var)*ghostcell_dim_0_characteristic_var +
                        (k + num_ghosts_2_characteristic_var)*ghostcell_dim_0_characteristic_var*
                            ghostcell_dim_1_characteristic_var;
                    
                    const int idx_vel = (i + num_ghosts_0_vel) +
                        (j + num_ghosts_1_vel)*ghostcell_dim_0_vel +
                        (k + idx_offset_vel + num_ghosts_2_vel)*ghostcell_dim_0_vel*
                            ghostcell_dim_1_vel;
                    
                    const int idx_p = (i + num_ghosts_0_p) +
                        (j + num_ghosts_1_p)*ghostcell_dim_0_p +
                        (k + idx_offset_p + num_ghosts_2_p)*ghostcell_dim_0_p*
                            ghostcell_dim_1_p;
                    
                    W[0][idx_face] = V[d_num_species + 2][idx_vel] -
                        double(1)/(rho_average[idx_face]*c_average[idx_face])*V[d_num_species + 3][idx_p];
                    
                    W[d_num_species + 1][idx_face] = V[d_num_species][idx_vel];
                    
                    W[d_num_species + 2][idx_face] = V[d_num_species + 1][idx_vel];
                    
                    W[d_num_species + 3][idx_face] = V[d_num_species + 2][idx_vel] +
                        double(1)/(rho_average[idx_face]*c_average[idx_face])*V[d_num_species + 3][idx_p];
                }
            }
        }
    }
}


/*
 * Compute global side data of conservative variables from characteristic variables.
 */
void
FlowModelFourEqnConservative::computeGlobalSideDataConservativeVariablesFromCharacteristicVariables(
    std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& characteristic_variables,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& projection_variables)
{
    NULL_USE(conservative_variables);
    NULL_USE(characteristic_variables);
    NULL_USE(projection_variables);
    
    TBOX_ERROR(d_object_name
        << ": FlowModelFourEqnConservative::"
        << "computeGlobalSideDataConservativeVariablesFromCharacteristicVariables()\n"
        << "Method computeGlobalSideDataConservativeVariablesFromCharacteristicVariables()"
        << " is not yet implemented."
        << std::endl);
}


/*
 * Compute global side data of primitive variables from characteristic variables.
 */
void
FlowModelFourEqnConservative::computeGlobalSideDataPrimitiveVariablesFromCharacteristicVariables(
    std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& characteristic_variables,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& projection_variables)
{
    /*
     * Get the numbers of ghost cells of the variables.
     */
    
    const hier::IntVector num_ghosts_primitive_var = primitive_variables[0]->
        getGhostCellWidth();
    
    const hier::IntVector num_ghosts_characteristic_var = characteristic_variables[0]->
        getGhostCellWidth();
    
    const hier::IntVector num_ghosts_projection_var = projection_variables[0]->getGhostCellWidth();
    
    /*
     * Get the ghost cell dimensions of characteristic variables.
     */
    
    const hier::IntVector ghostcell_dims_characteristic_var = characteristic_variables[0]->
        getGhostBox().numberCells();
    
    /*
     * Check the size of variables.
     */
    
    if (static_cast<int>(primitive_variables.size()) != d_num_eqn)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::"
            << "computeGlobalSideDataPrimitiveVariablesFromCharacteristicVariables()\n"
            << "The number of characteristic variables are incorrect."
            << std::endl);
    }
    if (static_cast<int>(characteristic_variables.size()) != d_num_eqn)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::"
            << "computeGlobalSideDataPrimitiveVariablesFromCharacteristicVariables()\n"
            << "The number of primitive variables are incorrect."
            << std::endl);
    }
    if (static_cast<int>(projection_variables.size()) != d_num_species + 2)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::"
            << "computeGlobalSideDataPrimitiveVariablesFromCharacteristicVariables()\n"
            << "There should be number of species projection plus two variables."
            << std::endl);
    }
    
    /*
     * Check potential failures.
     */
    
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        const hier::IntVector interior_dims_primitive_var =
            primitive_variables[ei]->getBox().numberCells();
        
        if (interior_dims_primitive_var != d_interior_dims)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::"
                << "computeGlobalSideDataPrimitiveVariablesFromCharacteristicVariables()\n"
                << "The interior dimension of the primitive variables does not match that of patch."
                << std::endl);
        }
    }
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        const hier::IntVector interior_dims_characteristic_var =
            characteristic_variables[ei]->getBox().numberCells();
        
        if (interior_dims_characteristic_var != d_interior_dims)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::"
                << "computeGlobalSideDataPrimitiveVariablesFromCharacteristicVariables()\n"
                << "The interior dimension of the characteristic variables does not match that of patch."
                << std::endl);
        }
    }
    for (int vi = 0; vi < d_num_species + 2; vi++)
    {
        const hier::IntVector interior_dims_projection_var = projection_variables[vi]->getBox().numberCells();
        if (interior_dims_projection_var != d_interior_dims)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::"
                << "computeGlobalSideDataPrimitiveVariablesFromCharacteristicVariables()\n"
                << "The interior dimension of the projection variables does not match that of patch."
                << std::endl);
        }
    }
    
    for (int ei = 1; ei < d_num_eqn; ei++)
    {
        if (num_ghosts_primitive_var != primitive_variables[ei]->getGhostCellWidth())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::"
                << "computeGlobalSideDataPrimitiveVariablesFromCharacteristicVariables()\n"
                << "The primitive variables don't have same ghost cell width."
                << std::endl);
        }
    }
    for (int ei = 1; ei < d_num_eqn; ei++)
    {
        if (num_ghosts_characteristic_var != characteristic_variables[ei]->getGhostCellWidth())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::"
                << "computeGlobalSideDataPrimitiveVariablesFromCharacteristicVariables()\n"
                << "The characteristic variables don't have same ghost cell width."
                << std::endl);
        }
    }
    for (int vi = 1; vi < d_num_species + 2; vi++)
    {
        if (num_ghosts_projection_var != projection_variables[vi]->getGhostCellWidth())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::"
                << "computeGlobalSideDataPrimitiveVariablesFromCharacteristicVariables()\n"
                << "The projection variables don't have same ghost cell width."
                << std::endl);
        }
    }
    
    if (num_ghosts_projection_var != num_ghosts_primitive_var)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::"
            << "computeGlobalSideDataPrimitiveVariablesFromCharacteristicVariables()\n"
            << "The ghost cell width of the projection variables does not match that of"
            << " primitive variables."
            << std::endl);
    }
    if (num_ghosts_projection_var != num_ghosts_characteristic_var)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::"
            << "computeGlobalSideDataPrimitiveVariablesFromCharacteristicVariables()\n"
            << "The ghost cell width of the projection variables does not match that of"
            << " characteristic variables."
            << std::endl);
    }
    
    /*
     * Declare containers to store pointers to different data.
     */
    
    std::vector<double*> V;
    std::vector<double*> W;
    V.resize(d_num_eqn);
    W.resize(d_num_eqn);
    
    std::vector<double*> rho_Y_average;
    rho_Y_average.resize(d_num_species);
    double* rho_average = nullptr;
    double* c_average = nullptr;
    
    if (d_dim == tbox::Dimension(1))
    {
        const int interior_dim_0 = d_interior_dims[0];
        
        const int num_ghosts_0_characteristic_var = num_ghosts_characteristic_var[0];
        
        /*
         * Compute the primitive variables in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V[ei] = primitive_variables[ei]->getPointer(0);
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(0);
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y_average[si] = projection_variables[si]->getPointer(0);
        }
        rho_average = projection_variables[d_num_species]->getPointer(0);
        c_average = projection_variables[d_num_species + 1]->getPointer(0);
        
        for (int si = 0; si < d_num_species; si++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = -num_ghosts_0_characteristic_var;
                 i < interior_dim_0 + 1 + num_ghosts_0_characteristic_var;
                 i++)
            {
                // Compute the linear index.
                const int idx_face = i + num_ghosts_0_characteristic_var;
                
                V[si][idx_face] = -double(1)/double(2)*rho_Y_average[si][idx_face]/c_average[idx_face]*W[0][idx_face] +
                    W[si + 1][idx_face] + double(1)/double(2)*rho_Y_average[si][idx_face]/c_average[idx_face]*
                        W[d_num_eqn - 1][idx_face];
            }
        }
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = -num_ghosts_0_characteristic_var;
             i < interior_dim_0 + 1 + num_ghosts_0_characteristic_var;
             i++)
        {
            // Compute the linear index.
            const int idx_face = i + num_ghosts_0_characteristic_var;
            
            V[d_num_species][idx_face] = double(1)/double(2)*W[0][idx_face] +
                double(1)/double(2)*W[d_num_eqn - 1][idx_face];
            
            V[d_num_species + 1][idx_face] = -double(1)/double(2)*rho_average[idx_face]*c_average[idx_face]*
                W[0][idx_face] + double(1)/double(2)*rho_average[idx_face]*c_average[idx_face]*
                    W[d_num_eqn - 1][idx_face];
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int interior_dim_0 = d_interior_dims[0];
        const int interior_dim_1 = d_interior_dims[1];
        
        const int num_ghosts_0_characteristic_var = num_ghosts_characteristic_var[0];
        const int num_ghosts_1_characteristic_var = num_ghosts_characteristic_var[1];
        const int ghostcell_dim_0_characteristic_var = ghostcell_dims_characteristic_var[0];
        
        /*
         * Compute the primitive variables in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V[ei] = primitive_variables[ei]->getPointer(0);
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(0);
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y_average[si] = projection_variables[si]->getPointer(0);
        }
        rho_average = projection_variables[d_num_species]->getPointer(0);
        c_average = projection_variables[d_num_species + 1]->getPointer(0);
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = -num_ghosts_0_characteristic_var;
                     i < interior_dim_0 + 1 + num_ghosts_0_characteristic_var;
                     i++)
                {
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_characteristic_var) +
                        (j + num_ghosts_1_characteristic_var)*(ghostcell_dim_0_characteristic_var + 1);
                    
                    V[si][idx_face] = -double(1)/double(2)*rho_Y_average[si][idx_face]/c_average[idx_face]*W[0][idx_face] +
                        W[si + 1][idx_face] + double(1)/double(2)*rho_Y_average[si][idx_face]/c_average[idx_face]*
                            W[d_num_eqn - 1][idx_face];
                }
            }
        }
        
        for (int j = 0; j < interior_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = -num_ghosts_0_characteristic_var;
                 i < interior_dim_0 + 1 + num_ghosts_0_characteristic_var;
                 i++)
            {
                // Compute the linear index.
                const int idx_face = (i + num_ghosts_0_characteristic_var) +
                    (j + num_ghosts_1_characteristic_var)*(ghostcell_dim_0_characteristic_var + 1);
                
                V[d_num_species][idx_face] = double(1)/double(2)*W[0][idx_face] +
                    double(1)/double(2)*W[d_num_eqn - 1][idx_face];
                
                V[d_num_species + 1][idx_face] = W[d_num_species + 1][idx_face];
                
                V[d_num_species + 2][idx_face] = -double(1)/double(2)*rho_average[idx_face]*c_average[idx_face]*
                    W[0][idx_face] + double(1)/double(2)*rho_average[idx_face]*c_average[idx_face]*
                        W[d_num_eqn - 1][idx_face];
            }
        }
        
        /*
         * Compute the primitive variables in the y-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V[ei] = primitive_variables[ei]->getPointer(1);
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(1);
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y_average[si] = projection_variables[si]->getPointer(1);
        }
        rho_average = projection_variables[d_num_species]->getPointer(1);
        c_average = projection_variables[d_num_species + 1]->getPointer(1);
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = -num_ghosts_1_characteristic_var;
                 j < interior_dim_1 + 1 + num_ghosts_1_characteristic_var;
                 j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_characteristic_var) +
                        (j + num_ghosts_1_characteristic_var)*ghostcell_dim_0_characteristic_var;
                    
                    V[si][idx_face] = -double(1)/double(2)*rho_Y_average[si][idx_face]/c_average[idx_face]*W[0][idx_face] +
                        W[si + 1][idx_face] + double(1)/double(2)*rho_Y_average[si][idx_face]/c_average[idx_face]*
                            W[d_num_eqn - 1][idx_face];
                }
            }
        }
        
        for (int j = -num_ghosts_1_characteristic_var;
             j < interior_dim_1 + 1 + num_ghosts_1_characteristic_var;
             j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear index.
                const int idx_face = (i + num_ghosts_0_characteristic_var) +
                    (j + num_ghosts_1_characteristic_var)*ghostcell_dim_0_characteristic_var;
                
                V[d_num_species][idx_face] = W[d_num_species + 1][idx_face];
                
                V[d_num_species + 1][idx_face] = double(1)/double(2)*W[0][idx_face] +
                    double(1)/double(2)*W[d_num_eqn - 1][idx_face];
                
                V[d_num_species + 2][idx_face] = -double(1)/double(2)*rho_average[idx_face]*c_average[idx_face]*
                    W[0][idx_face] + double(1)/double(2)*rho_average[idx_face]*c_average[idx_face]*
                        W[d_num_eqn - 1][idx_face];
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int interior_dim_0 = d_interior_dims[0];
        const int interior_dim_1 = d_interior_dims[1];
        const int interior_dim_2 = d_interior_dims[2];
        
        const int num_ghosts_0_characteristic_var = num_ghosts_characteristic_var[0];
        const int num_ghosts_1_characteristic_var = num_ghosts_characteristic_var[1];
        const int num_ghosts_2_characteristic_var = num_ghosts_characteristic_var[2];
        const int ghostcell_dim_0_characteristic_var = ghostcell_dims_characteristic_var[0];
        const int ghostcell_dim_1_characteristic_var = ghostcell_dims_characteristic_var[1];
        
        /*
         * Compute the primitive variables in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V[ei] = primitive_variables[ei]->getPointer(0);
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(0);
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y_average[si] = projection_variables[si]->getPointer(0);
        }
        rho_average = projection_variables[d_num_species]->getPointer(0);
        c_average = projection_variables[d_num_species + 1]->getPointer(0);
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = -num_ghosts_0_characteristic_var;
                         i < interior_dim_0 + 1 + num_ghosts_0_characteristic_var;
                         i++)
                    {
                        // Compute the linear index.
                        const int idx_face = (i + num_ghosts_0_characteristic_var) +
                            (j + num_ghosts_1_characteristic_var)*(ghostcell_dim_0_characteristic_var + 1) +
                            (k + num_ghosts_2_characteristic_var)*(ghostcell_dim_0_characteristic_var + 1)*
                                ghostcell_dim_1_characteristic_var;
                        
                        V[si][idx_face] = -double(1)/double(2)*rho_Y_average[si][idx_face]/c_average[idx_face]*W[0][idx_face] +
                            W[si + 1][idx_face] + double(1)/double(2)*rho_Y_average[si][idx_face]/c_average[idx_face]*
                                W[d_num_eqn - 1][idx_face];
                    }
                }
            }
        }
        
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = -num_ghosts_0_characteristic_var;
                     i < interior_dim_0 + 1 + num_ghosts_0_characteristic_var;
                     i++)
                {
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_characteristic_var) +
                        (j + num_ghosts_1_characteristic_var)*(ghostcell_dim_0_characteristic_var + 1) +
                        (k + num_ghosts_2_characteristic_var)*(ghostcell_dim_0_characteristic_var + 1)*
                            ghostcell_dim_1_characteristic_var;
                    
                    V[d_num_species][idx_face] = double(1)/double(2)*W[0][idx_face] +
                        double(1)/double(2)*W[d_num_eqn - 1][idx_face];
                    
                    V[d_num_species + 1][idx_face] = W[d_num_species + 1][idx_face];
                    
                    V[d_num_species + 2][idx_face] = W[d_num_species + 2][idx_face];
                    
                    V[d_num_species + 3][idx_face] = -double(1)/double(2)*rho_average[idx_face]*c_average[idx_face]*
                        W[0][idx_face] + double(1)/double(2)*rho_average[idx_face]*c_average[idx_face]*
                            W[d_num_eqn - 1][idx_face];
                }
            }
        }
        
        /*
         * Compute the primitive variables in the y-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V[ei] = primitive_variables[ei]->getPointer(1);
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(1);
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y_average[si] = projection_variables[si]->getPointer(1);
        }
        rho_average = projection_variables[d_num_species]->getPointer(1);
        c_average = projection_variables[d_num_species + 1]->getPointer(1);
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = -num_ghosts_1_characteristic_var;
                     j < interior_dim_1 + 1 + num_ghosts_1_characteristic_var;
                     j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx_face = (i + num_ghosts_0_characteristic_var) +
                            (j + num_ghosts_1_characteristic_var)*ghostcell_dim_0_characteristic_var +
                            (k + num_ghosts_2_characteristic_var)*ghostcell_dim_0_characteristic_var*
                                (ghostcell_dim_1_characteristic_var + 1);
                        
                        V[si][idx_face] = -double(1)/double(2)*rho_Y_average[si][idx_face]/c_average[idx_face]*W[0][idx_face] +
                            W[si + 1][idx_face] + double(1)/double(2)*rho_Y_average[si][idx_face]/c_average[idx_face]*
                                W[d_num_eqn - 1][idx_face];
                    }
                }
            }
        }
        
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = -num_ghosts_1_characteristic_var;
                 j < interior_dim_1 + 1 + num_ghosts_1_characteristic_var;
                 j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_characteristic_var) +
                        (j + num_ghosts_1_characteristic_var)*ghostcell_dim_0_characteristic_var +
                        (k + num_ghosts_2_characteristic_var)*ghostcell_dim_0_characteristic_var*
                            (ghostcell_dim_1_characteristic_var + 1);
                    
                    V[d_num_species][idx_face] = W[d_num_species + 1][idx_face];
                    
                    V[d_num_species + 1][idx_face] = double(1)/double(2)*W[0][idx_face] +
                        double(1)/double(2)*W[d_num_eqn - 1][idx_face];
                    
                    V[d_num_species + 2][idx_face] = W[d_num_species + 2][idx_face];
                    
                    V[d_num_species + 3][idx_face] = -double(1)/double(2)*rho_average[idx_face]*c_average[idx_face]*
                        W[0][idx_face] + double(1)/double(2)*rho_average[idx_face]*c_average[idx_face]*
                            W[d_num_eqn - 1][idx_face];
                }
            }
        }
        
        /*
         * Compute the primitive variables in the z-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V[ei] = primitive_variables[ei]->getPointer(2);
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(2);
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y_average[si] = projection_variables[si]->getPointer(2);
        }
        rho_average = projection_variables[d_num_species]->getPointer(2);
        c_average = projection_variables[d_num_species + 1]->getPointer(2);
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = -num_ghosts_2_characteristic_var;
                 k < interior_dim_2 + 1 + num_ghosts_2_characteristic_var;
                 k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx_face = (i + num_ghosts_0_characteristic_var) +
                            (j + num_ghosts_1_characteristic_var)*ghostcell_dim_0_characteristic_var +
                            (k + num_ghosts_2_characteristic_var)*ghostcell_dim_0_characteristic_var*
                                ghostcell_dim_1_characteristic_var;
                        
                        V[si][idx_face] = -double(1)/double(2)*rho_Y_average[si][idx_face]/c_average[idx_face]*W[0][idx_face] +
                            W[si + 1][idx_face] + double(1)/double(2)*rho_Y_average[si][idx_face]/c_average[idx_face]*
                                W[d_num_eqn - 1][idx_face];
                    }
                }
            }
        }
        
        for (int k = -num_ghosts_2_characteristic_var;
             k < interior_dim_2 + 1 + num_ghosts_2_characteristic_var;
             k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_characteristic_var) +
                        (j + num_ghosts_1_characteristic_var)*ghostcell_dim_0_characteristic_var +
                        (k + num_ghosts_2_characteristic_var)*ghostcell_dim_0_characteristic_var*
                            ghostcell_dim_1_characteristic_var;
                    
                    V[d_num_species][idx_face] = W[d_num_species + 1][idx_face];
                    
                    V[d_num_species + 1][idx_face] = W[d_num_species + 2][idx_face];
                    
                    V[d_num_species + 2][idx_face] = double(1)/double(2)*W[0][idx_face] +
                        double(1)/double(2)*W[d_num_eqn - 1][idx_face];
                    
                    V[d_num_species + 3][idx_face] = -double(1)/double(2)*rho_average[idx_face]*c_average[idx_face]*W[0][idx_face] +
                        double(1)/double(2)*rho_average[idx_face]*c_average[idx_face]*W[d_num_eqn - 1][idx_face];
                }
            }
        }
    }
}


/*
 * Compute the local intercell quantities with conservative variables on each side of the face
 * from Riemann solver at face.
 * flux_face: Convective flux at face.
 * velocity_face: Velocity at face.
 */
void
FlowModelFourEqnConservative::computeLocalFaceFluxAndVelocityFromRiemannSolverWithConservativeVariables(
    std::vector<boost::reference_wrapper<double> >& flux_face,
    std::vector<boost::reference_wrapper<double> >& velocity_face,
    const std::vector<boost::reference_wrapper<double> >& conservative_variables_minus,
    const std::vector<boost::reference_wrapper<double> >& conservative_variables_plus,
    const DIRECTION::TYPE& direction,
    const RIEMANN_SOLVER::TYPE& Riemann_solver)
{
    switch (Riemann_solver)
    {
        case RIEMANN_SOLVER::HLLC:
        {
            d_Riemann_solver_HLLC->computeIntercellFluxFromConservativeVariables(
                flux_face,
                conservative_variables_minus,
                conservative_variables_plus,
                direction);
            
            break;
        }
        case RIEMANN_SOLVER::HLLC_HLL:
        {
            d_Riemann_solver_HLLC_HLL->computeIntercellFluxFromConservativeVariables(
                flux_face,
                conservative_variables_minus,
                conservative_variables_plus,
                direction);
            
            break;
        }
        default:
        {
            TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::"
            << "computeLocalFaceFluxAndVelocityFromRiemannSolverWithConservativeVariables()\n"
            << "Unknown Riemann solver required."
            << std::endl);
        }
    }
}


/*
 * Compute the local intercell quantities with primitive variables on each side of the face
 * from Riemann solver at face.
 * flux_face: Convective flux at face.
 * velocity_face: Velocity at face.
 */
void
FlowModelFourEqnConservative::computeLocalFaceFluxAndVelocityFromRiemannSolverWithPrimitiveVariables(
    std::vector<boost::reference_wrapper<double> >& flux_face,
    std::vector<boost::reference_wrapper<double> >& velocity_face,
    const std::vector<boost::reference_wrapper<double> >& primitive_variables_minus,
    const std::vector<boost::reference_wrapper<double> >& primitive_variables_plus,
    const DIRECTION::TYPE& direction,
    const RIEMANN_SOLVER::TYPE& Riemann_solver)
{
    switch (Riemann_solver)
    {
        case RIEMANN_SOLVER::HLLC:
        {
            d_Riemann_solver_HLLC->computeIntercellFluxFromPrimitiveVariables(
                flux_face,
                primitive_variables_minus,
                primitive_variables_plus,
                direction);
            
            break;
        }
        case RIEMANN_SOLVER::HLLC_HLL:
        {
            d_Riemann_solver_HLLC_HLL->computeIntercellFluxFromPrimitiveVariables(
                flux_face,
                primitive_variables_minus,
                primitive_variables_plus,
                direction);
            
            break;
        }
        default:
        {
            TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::"
            << "computeLocalFaceFluxAndVelocityFromRiemannSolverWithPrimitiveVariables()\n"
            << "Unknown Riemann solver required."
            << std::endl);
        }
    }
}


/*
 * Check whether the given side conservative variables are within the bounds.
 */
void
FlowModelFourEqnConservative::checkGlobalSideDataConservativeVariablesBounded(
    boost::shared_ptr<pdat::SideData<int> >& bounded_flag,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables)
{
    /*
     * Get the numbers of ghost cells of the variables.
     */
    
    const hier::IntVector num_ghosts_flag = bounded_flag->getGhostCellWidth();
    
    const hier::IntVector num_ghosts_conservative_var = conservative_variables[0]->
        getGhostCellWidth();
    
    /*
     * Get the ghost cell dimensions of conservative variables.
     */
    
    const hier::IntVector ghostcell_dims_conservative_var = conservative_variables[0]->
        getGhostBox().numberCells();
    
    /*
     * Check the size of variables.
     */
    
    if (static_cast<int>(conservative_variables.size()) != d_num_eqn)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::"
            << "checkGlobalSideDataConservativeVariablesBounded()\n"
            << "The number of conservative variables are incorrect."
            << std::endl);
    }
    
    /*
     * Check potential failures.
     */
    
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        const hier::IntVector interior_dims_conservative_var =
            conservative_variables[ei]->getBox().numberCells();
        
        if (interior_dims_conservative_var != d_interior_dims)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::"
                << "checkGlobalSideDataConservativeVariablesBounded()\n"
                << "The interior dimension of the conservative variables does not match that of patch."
                << std::endl);
        }
    }
    const hier::IntVector interior_dims_flag = bounded_flag->getBox().numberCells();
    if (interior_dims_flag != d_interior_dims)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::"
            << "checkGlobalSideDataConservativeVariablesBounded()\n"
            << "The interior dimension of the flag does not match that of patch."
            << std::endl);
    }
    
    for (int ei = 1; ei < d_num_eqn; ei++)
    {
        if (num_ghosts_conservative_var != conservative_variables[ei]->getGhostCellWidth())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::"
                << "checkGlobalSideDataConservativeVariablesBounded()\n"
                << "The conservative variables don't have same ghost cell width."
                << std::endl);
        }
    }
    
    if (num_ghosts_flag != num_ghosts_conservative_var)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::"
            << "checkGlobalSideDataConservativeVariablesBounded()\n"
            << "The ghost cell width of the flag does not match that of conservative variables."
            << std::endl);
    }
    
    bounded_flag->fillAll(1);
    
    // Create the side data of density.
    boost::shared_ptr<pdat::SideData<double> > data_density(
        new pdat::SideData<double>(d_interior_box, 1, num_ghosts_conservative_var));
    
    data_density->fillAll(double(0));
    
    /*
     * Declare containers to store pointers to different data.
     */
    
    int* are_bounded = nullptr;
    
    std::vector<double*> Q;
    Q.resize(d_num_eqn);
    
    double* rho = nullptr;
    
    if (d_dim == tbox::Dimension(1))
    {
        const int interior_dim_0 = d_interior_dims[0];
        
        const int num_ghosts_0_conservative_var = num_ghosts_conservative_var[0];
        
        /*
         * Check if conservative variables in the x-direction are bounded.
         */
        
        are_bounded = bounded_flag->getPointer(0);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            Q[ei] = conservative_variables[ei]->getPointer(0);
        }
        
        rho = data_density->getPointer(0);
        
        // Compute density.
        for (int si = 0; si < d_num_species; si++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = -num_ghosts_0_conservative_var;
                 i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                 i++)
            {
                // Compute the linear index.
                const int idx_face = i + num_ghosts_0_conservative_var;
                
                rho[idx_face] += Q[si][idx_face];
            }
        }
        
        // Check if mass fractions are bounded.
        for (int si = 0; si < d_num_species; si++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = -num_ghosts_0_conservative_var;
                 i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                 i++)
            {
                // Compute the linear index.
                const int idx_face = i + num_ghosts_0_conservative_var;
                
                const double Y = Q[si][idx_face]/rho[idx_face];
                
                if (Y > d_Y_bound_lo && Y < d_Y_bound_up)
                {
                    are_bounded[idx_face] &= 1;
                }
                else
                {
                    are_bounded[idx_face] &= 0;
                }
            }
        }
        
        // Check if density and total energy are bounded.
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = -num_ghosts_0_conservative_var;
             i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
             i++)
        {
            // Compute the linear index.
            const int idx_face = i + num_ghosts_0_conservative_var;
            
            if (rho[idx_face] > double(0))
            {
                are_bounded[idx_face] &= 1;
            }
            else
            {
                are_bounded[idx_face] &= 0;
            }
            
            if (Q[d_num_species + d_dim.getValue()][idx_face] > double(0))
            {
                are_bounded[idx_face] &= 1;
            }
            else
            {
                are_bounded[idx_face] &= 0;
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int interior_dim_0 = d_interior_dims[0];
        const int interior_dim_1 = d_interior_dims[1];
        
        const int num_ghosts_0_conservative_var = num_ghosts_conservative_var[0];
        const int num_ghosts_1_conservative_var = num_ghosts_conservative_var[1];
        const int ghostcell_dim_0_conservative_var = ghostcell_dims_conservative_var[0];
        
        /*
         * Check if conservative variables in the x-direction are bounded.
         */
        
        are_bounded = bounded_flag->getPointer(0);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            Q[ei] = conservative_variables[ei]->getPointer(0);
        }
        
        rho = data_density->getPointer(0);
        
        // Compute density.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = -num_ghosts_0_conservative_var;
                     i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                     i++)
                {
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1);
                    
                    rho[idx_face] += Q[si][idx_face];
                }
            }
        }
        
        // Check if mass fractions are bounded.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = -num_ghosts_0_conservative_var;
                     i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                     i++)
                {
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1);
                    
                    const double Y = Q[si][idx_face]/rho[idx_face];
                    
                    if (Y > d_Y_bound_lo && Y < d_Y_bound_up)
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                }
            }
        }
        
        // Check if density and total energy are bounded.
        for (int j = 0; j < interior_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = -num_ghosts_0_conservative_var;
                 i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                 i++)
            {
                // Compute the linear index.
                const int idx_face = (i + num_ghosts_0_conservative_var) +
                    (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1);
                
                if (rho[idx_face] > double(0))
                {
                    are_bounded[idx_face] &= 1;
                }
                else
                {
                    are_bounded[idx_face] &= 0;
                }
                
                if (Q[d_num_species + d_dim.getValue()][idx_face] > double(0))
                {
                    are_bounded[idx_face] &= 1;
                }
                else
                {
                    are_bounded[idx_face] &= 0;
                }
            }
        }
        
        /*
         * Check if conservative variables in the y-direction are bounded.
         */
        
        are_bounded = bounded_flag->getPointer(1);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            Q[ei] = conservative_variables[ei]->getPointer(1);
        }
        
        rho = data_density->getPointer(1);
        
        // Compute density.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = -num_ghosts_1_conservative_var;
                 j < interior_dim_1 + 1 + num_ghosts_1_conservative_var;
                 j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var;
                    
                    rho[idx_face] += Q[si][idx_face];
                }
            }
        }
        
        // Check if mass fractions are bounded.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = -num_ghosts_1_conservative_var;
                 j < interior_dim_1 + 1 + num_ghosts_1_conservative_var;
                 j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var;
                    
                    const double Y = Q[si][idx_face]/rho[idx_face];
                    
                    if (Y > d_Y_bound_lo && Y < d_Y_bound_up)
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                }
            }
        }
        
        // Check if density and total energy are bounded.
        for (int j = -num_ghosts_1_conservative_var;
             j < interior_dim_1 + 1 + num_ghosts_1_conservative_var;
             j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear index.
                const int idx_face = (i + num_ghosts_0_conservative_var) +
                    (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var;
                
                if (rho[idx_face] > double(0))
                {
                    are_bounded[idx_face] &= 1;
                }
                else
                {
                    are_bounded[idx_face] &= 0;
                }
                
                if (Q[d_num_species + d_dim.getValue()][idx_face] > double(0))
                {
                    are_bounded[idx_face] &= 1;
                }
                else
                {
                    are_bounded[idx_face] &= 0;
                }
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int interior_dim_0 = d_interior_dims[0];
        const int interior_dim_1 = d_interior_dims[1];
        const int interior_dim_2 = d_interior_dims[2];
        
        const int num_ghosts_0_conservative_var = num_ghosts_conservative_var[0];
        const int num_ghosts_1_conservative_var = num_ghosts_conservative_var[1];
        const int num_ghosts_2_conservative_var = num_ghosts_conservative_var[2];
        const int ghostcell_dim_0_conservative_var = ghostcell_dims_conservative_var[0];
        const int ghostcell_dim_1_conservative_var = ghostcell_dims_conservative_var[1];
        
        /*
         * Check if conservative variables in the x-direction are bounded.
         */
        
        are_bounded = bounded_flag->getPointer(0);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            Q[ei] = conservative_variables[ei]->getPointer(0);
        }
        
        rho = data_density->getPointer(0);
        
        // Compute density.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = -num_ghosts_0_conservative_var;
                         i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                         i++)
                    {
                        // Compute the linear indices.
                        const int idx_face = (i + num_ghosts_0_conservative_var) +
                            (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1) +
                            (k + num_ghosts_2_conservative_var)*(ghostcell_dim_0_conservative_var + 1)*
                                ghostcell_dim_1_conservative_var;
                        
                        rho[idx_face] += Q[si][idx_face];
                    }
                }
            }
        }
        
        // Check if mass fractions are bounded.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = -num_ghosts_0_conservative_var;
                         i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                         i++)
                    {
                        // Compute the linear indices.
                        const int idx_face = (i + num_ghosts_0_conservative_var) +
                            (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1) +
                            (k + num_ghosts_2_conservative_var)*(ghostcell_dim_0_conservative_var + 1)*
                                ghostcell_dim_1_conservative_var;
                        
                        const double Y = Q[si][idx_face]/rho[idx_face];
                        
                        if (Y > d_Y_bound_lo && Y < d_Y_bound_up)
                        {
                            are_bounded[idx_face] &= 1;
                        }
                        else
                        {
                            are_bounded[idx_face] &= 0;
                        }
                    }
                }
            }
        }
        
        // Check if density and total energy are bounded.
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = -num_ghosts_0_conservative_var;
                     i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                     i++)
                {
                    // Compute the linear indices.
                    const int idx_face = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1) +
                        (k + num_ghosts_2_conservative_var)*(ghostcell_dim_0_conservative_var + 1)*
                            ghostcell_dim_1_conservative_var;
                    
                    if (rho[idx_face] > double(0))
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                    
                    if (Q[d_num_species + d_dim.getValue()][idx_face] > double(0))
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                }
            }
        }
        
        /*
         * Check if conservative variables in the y-direction are bounded.
         */
        
        are_bounded = bounded_flag->getPointer(1);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            Q[ei] = conservative_variables[ei]->getPointer(1);
        }
        
        rho = data_density->getPointer(1);
        
        // Compute density.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = -num_ghosts_1_conservative_var;
                     j < interior_dim_1 + 1 + num_ghosts_1_conservative_var;
                     j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx_face = (i + num_ghosts_0_conservative_var) +
                            (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                            (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                                (ghostcell_dim_1_conservative_var + 1);
                        
                        rho[idx_face] += Q[si][idx_face];
                    }
                }
            }
        }
        
        // Check if mass fractions are bounded.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = -num_ghosts_1_conservative_var;
                     j < interior_dim_1 + 1 + num_ghosts_1_conservative_var;
                     j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_face = (i + num_ghosts_0_conservative_var) +
                            (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                            (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                                (ghostcell_dim_1_conservative_var + 1);
                        
                        const double Y = Q[si][idx_face]/rho[idx_face];
                        
                        if (Y > d_Y_bound_lo && Y < d_Y_bound_up)
                        {
                            are_bounded[idx_face] &= 1;
                        }
                        else
                        {
                            are_bounded[idx_face] &= 0;
                        }
                    }
                }
            }
        }
        
        // Check if density and total energy are bounded.
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = -num_ghosts_1_conservative_var;
                 j < interior_dim_1 + 1 + num_ghosts_1_conservative_var;
                 j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                        (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                            (ghostcell_dim_1_conservative_var + 1);
                    
                    if (rho[idx_face] > double(0))
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                    
                    if (Q[d_num_species + d_dim.getValue()][idx_face] > double(0))
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                }
            }
        }
        
        /*
         * Check if conservative variables in the z-direction are bounded.
         */
        
        are_bounded = bounded_flag->getPointer(2);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            Q[ei] = conservative_variables[ei]->getPointer(2);
        }
        
        rho = data_density->getPointer(2);
        
        // Compute density.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = -num_ghosts_2_conservative_var;
                 k < interior_dim_2 + 1 + num_ghosts_2_conservative_var;
                 k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx_face = (i + num_ghosts_0_conservative_var) +
                            (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                            (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                                ghostcell_dim_1_conservative_var;
                        
                        rho[idx_face] += Q[si][idx_face];
                    }
                }
            }
        }
        
        // Check if mass fractions are bounded.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = -num_ghosts_2_conservative_var;
                 k < interior_dim_2 + 1 + num_ghosts_2_conservative_var;
                 k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_face = (i + num_ghosts_0_conservative_var) +
                            (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                            (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                                ghostcell_dim_1_conservative_var;
                        
                        const double Y = Q[si][idx_face]/rho[idx_face];
                        
                        if (Y > d_Y_bound_lo && Y < d_Y_bound_up)
                        {
                            are_bounded[idx_face] &= 1;
                        }
                        else
                        {
                            are_bounded[idx_face] &= 0;
                        }
                    }
                }
            }
        }
        
        // Check if density and total energy are bounded.
        for (int k = -num_ghosts_2_conservative_var;
             k < interior_dim_2 + 1 + num_ghosts_2_conservative_var;
             k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                        (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                            ghostcell_dim_1_conservative_var;
                    
                    if (rho[idx_face] > double(0))
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                    
                    if (Q[d_num_species + d_dim.getValue()][idx_face] > double(0))
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                }
            }
        }
    }
}


/*
 * Check whether the given side primitive variables are within the bounds.
 */
void
FlowModelFourEqnConservative::checkGlobalSideDataPrimitiveVariablesBounded(
    boost::shared_ptr<pdat::SideData<int> >& bounded_flag,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables)
{
    /*
     * Get the numbers of ghost cells of the variables.
     */
    
    const hier::IntVector num_ghosts_flag = bounded_flag->getGhostCellWidth();
    
    const hier::IntVector num_ghosts_primitive_var = primitive_variables[0]->
        getGhostCellWidth();
    
    /*
     * Get the ghost cell dimensions of primitive variables.
     */
    
    const hier::IntVector ghostcell_dims_primitive_var = primitive_variables[0]->
        getGhostBox().numberCells();
    
    /*
     * Check the size of variables.
     */
    
    if (static_cast<int>(primitive_variables.size()) != d_num_eqn)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::"
            << "checkGlobalSideDataPrimitiveVariablesBounded()\n"
            << "The number of primitive variables are incorrect."
            << std::endl);
    }
    
    /*
     * Check potential failures.
     */
    
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        const hier::IntVector interior_dims_primitive_var =
            primitive_variables[ei]->getBox().numberCells();
        
        if (interior_dims_primitive_var != d_interior_dims)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::"
                << "checkGlobalSideDataPrimitiveVariablesBounded()\n"
                << "The interior dimension of the primitive variables does not match that of patch."
                << std::endl);
        }
    }
    const hier::IntVector interior_dims_flag = bounded_flag->getBox().numberCells();
    if (interior_dims_flag != d_interior_dims)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::"
            << "checkGlobalSideDataPrimitiveVariablesBounded()\n"
            << "The interior dimension of the flag does not match that of patch."
            << std::endl);
    }
    
    for (int ei = 1; ei < d_num_eqn; ei++)
    {
        if (num_ghosts_primitive_var != primitive_variables[ei]->getGhostCellWidth())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::"
                << "checkGlobalSideDataPrimitiveVariablesBounded()\n"
                << "The primitive variables don't have same ghost cell width."
                << std::endl);
        }
    }
    
    if (num_ghosts_flag != num_ghosts_primitive_var)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::"
            << "checkGlobalSideDataPrimitiveVariablesBounded()\n"
            << "The ghost cell width of the flag does not match that of primitive variables."
            << std::endl);
    }
    
    bounded_flag->fillAll(1);
    
    // Create the side data of density.
    boost::shared_ptr<pdat::SideData<double> > data_density(
        new pdat::SideData<double>(d_interior_box, 1, num_ghosts_primitive_var));
    
    data_density->fillAll(double(0));
    
    /*
     * Declare containers to store pointers to different data.
     */
    
    int* are_bounded = nullptr;
    
    std::vector<double*> V;
    V.resize(d_num_eqn);
    
    double* rho = nullptr;
    
    if (d_dim == tbox::Dimension(1))
    {
        const int interior_dim_0 = d_interior_dims[0];
        
        const int num_ghosts_0_primitive_var = num_ghosts_primitive_var[0];
        
        /*
         * Check if primitive variables in the x-direction are bounded.
         */
        
        are_bounded = bounded_flag->getPointer(0);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V[ei] = primitive_variables[ei]->getPointer(0);
        }
        
        rho = data_density->getPointer(0);
        
        // Compute density.
        for (int si = 0; si < d_num_species; si++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = -num_ghosts_0_primitive_var;
                 i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                 i++)
            {
                // Compute the linear index.
                const int idx_face = i + num_ghosts_0_primitive_var;
                
                rho[idx_face] += V[si][idx_face];
            }
        }
        
        // Check if mass fractions are bounded.
        for (int si = 0; si < d_num_species; si++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = -num_ghosts_0_primitive_var;
                 i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                 i++)
            {
                // Compute the linear index.
                const int idx_face = i + num_ghosts_0_primitive_var;
                
                const double Y = V[si][idx_face]/rho[idx_face];
                
                if (Y > d_Y_bound_lo && Y < d_Y_bound_up)
                {
                    are_bounded[idx_face] &= 1;
                }
                else
                {
                    are_bounded[idx_face] &= 0;
                }
            }
        }
        
        // Check if density and pressure are bounded.
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = -num_ghosts_0_primitive_var;
             i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
             i++)
        {
            // Compute the linear index.
            const int idx_face = i + num_ghosts_0_primitive_var;
            
            if (rho[idx_face] > double(0))
            {
                are_bounded[idx_face] &= 1;
            }
            else
            {
                are_bounded[idx_face] &= 0;
            }
            
            if (V[d_num_species + d_dim.getValue()][idx_face] > double(0))
            {
                are_bounded[idx_face] &= 1;
            }
            else
            {
                are_bounded[idx_face] &= 0;
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int interior_dim_0 = d_interior_dims[0];
        const int interior_dim_1 = d_interior_dims[1];
        
        const int num_ghosts_0_primitive_var = num_ghosts_primitive_var[0];
        const int num_ghosts_1_primitive_var = num_ghosts_primitive_var[1];
        const int ghostcell_dim_0_primitive_var = ghostcell_dims_primitive_var[0];
        
        /*
         * Check if primitive variables in the x-direction are bounded.
         */
        
        are_bounded = bounded_flag->getPointer(0);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V[ei] = primitive_variables[ei]->getPointer(0);
        }
        
        rho = data_density->getPointer(0);
        
        // Compute density.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = -num_ghosts_0_primitive_var;
                     i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                     i++)
                {
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1);
                    
                    rho[idx_face] += V[si][idx_face];
                }
            }
        }
        
        // Check if mass fractions are bounded.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = -num_ghosts_0_primitive_var;
                     i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                     i++)
                {
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1);
                    
                    const double Y = V[si][idx_face]/rho[idx_face];
                    
                    if (Y > d_Y_bound_lo && Y < d_Y_bound_up)
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                }
            }
        }
        
        // Check if density and pressure are bounded.
        for (int j = 0; j < interior_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = -num_ghosts_0_primitive_var;
                 i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                 i++)
            {
                // Compute the linear index.
                const int idx_face = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1);
                
                if (rho[idx_face] > double(0))
                {
                    are_bounded[idx_face] &= 1;
                }
                else
                {
                    are_bounded[idx_face] &= 0;
                }
                
                if (V[d_num_species + d_dim.getValue()][idx_face] > double(0))
                {
                    are_bounded[idx_face] &= 1;
                }
                else
                {
                    are_bounded[idx_face] &= 0;
                }
            }
        }
        
        /*
         * Check if primitive variables in the y-direction are bounded.
         */
        
        are_bounded = bounded_flag->getPointer(1);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V[ei] = primitive_variables[ei]->getPointer(1);
        }
        
        rho = data_density->getPointer(1);
        
        // Compute density.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = -num_ghosts_1_primitive_var;
                 j < interior_dim_1 + 1 + num_ghosts_1_primitive_var;
                 j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var;
                    
                    rho[idx_face] += V[si][idx_face];
                }
            }
        }
        
        // Check if mass fractions are bounded.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = -num_ghosts_1_primitive_var;
                 j < interior_dim_1 + 1 + num_ghosts_1_primitive_var;
                 j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var;
                    
                    const double Y = V[si][idx_face]/rho[idx_face];
                    
                    if (Y > d_Y_bound_lo && Y < d_Y_bound_up)
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                }
            }
        }
        
        // Check if density and pressure are bounded.
        for (int j = -num_ghosts_1_primitive_var;
             j < interior_dim_1 + 1 + num_ghosts_1_primitive_var;
             j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear index.
                const int idx_face = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var;
                
                if (rho[idx_face] > double(0))
                {
                    are_bounded[idx_face] &= 1;
                }
                else
                {
                    are_bounded[idx_face] &= 0;
                }
                
                if (V[d_num_species + d_dim.getValue()][idx_face] > double(0))
                {
                    are_bounded[idx_face] &= 1;
                }
                else
                {
                    are_bounded[idx_face] &= 0;
                }
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int interior_dim_0 = d_interior_dims[0];
        const int interior_dim_1 = d_interior_dims[1];
        const int interior_dim_2 = d_interior_dims[2];
        
        const int num_ghosts_0_primitive_var = num_ghosts_primitive_var[0];
        const int num_ghosts_1_primitive_var = num_ghosts_primitive_var[1];
        const int num_ghosts_2_primitive_var = num_ghosts_primitive_var[2];
        const int ghostcell_dim_0_primitive_var = ghostcell_dims_primitive_var[0];
        const int ghostcell_dim_1_primitive_var = ghostcell_dims_primitive_var[1];
        
        /*
         * Check if primitive variables in the x-direction are bounded.
         */
        
        are_bounded = bounded_flag->getPointer(0);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V[ei] = primitive_variables[ei]->getPointer(0);
        }
        
        rho = data_density->getPointer(0);
        
        // Compute density.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = -num_ghosts_0_primitive_var;
                         i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                         i++)
                    {
                        // Compute the linear index.
                        const int idx_face = (i + num_ghosts_0_primitive_var) +
                            (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1) +
                            (k + num_ghosts_2_primitive_var)*(ghostcell_dim_0_primitive_var + 1)*
                                ghostcell_dim_1_primitive_var;
                        
                        rho[idx_face] += V[si][idx_face];
                    }
                }
            }
        }
        
        // Check if mass fractions are bounded.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = -num_ghosts_0_primitive_var;
                         i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                         i++)
                    {
                        // Compute the linear index.
                        const int idx_face = (i + num_ghosts_0_primitive_var) +
                            (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1) +
                            (k + num_ghosts_2_primitive_var)*(ghostcell_dim_0_primitive_var + 1)*
                                ghostcell_dim_1_primitive_var;
                        
                        const double Y = V[si][idx_face]/rho[idx_face];
                        
                        if (Y > d_Y_bound_lo && Y < d_Y_bound_up)
                        {
                            are_bounded[idx_face] &= 1;
                        }
                        else
                        {
                            are_bounded[idx_face] &= 0;
                        }
                    }
                }
            }
        }
        
        // Check if density and pressure are bounded.
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = -num_ghosts_0_primitive_var;
                     i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                     i++)
                {
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1) +
                        (k + num_ghosts_2_primitive_var)*(ghostcell_dim_0_primitive_var + 1)*
                            ghostcell_dim_1_primitive_var;
                    
                    if (rho[idx_face] > double(0))
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                    
                    if (V[d_num_species + d_dim.getValue()][idx_face] > double(0))
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                }
            }
        }
        
        /*
         * Check if primitive variables in the y-direction are bounded.
         */
        
        are_bounded = bounded_flag->getPointer(1);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V[ei] = primitive_variables[ei]->getPointer(1);
        }
        
        rho = data_density->getPointer(1);
        
        // Compute density.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = -num_ghosts_1_primitive_var;
                     j < interior_dim_1 + 1 + num_ghosts_1_primitive_var;
                     j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx_face = (i + num_ghosts_0_primitive_var) +
                            (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                            (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                                (ghostcell_dim_1_primitive_var + 1);
                        
                        rho[idx_face] += V[si][idx_face];
                    }
                }
            }
        }
        
        // Check if mass fractions are bounded.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = -num_ghosts_1_primitive_var;
                     j < interior_dim_1 + 1 + num_ghosts_1_primitive_var;
                     j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx_face = (i + num_ghosts_0_primitive_var) +
                            (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                            (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                                (ghostcell_dim_1_primitive_var + 1);
                        
                        const double Y = V[si][idx_face]/rho[idx_face];
                        
                        if (Y > d_Y_bound_lo && Y < d_Y_bound_up)
                        {
                            are_bounded[idx_face] &= 1;
                        }
                        else
                        {
                            are_bounded[idx_face] &= 0;
                        }
                    }
                }
            }
        }
        
        // Check if density and pressure are bounded.
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = -num_ghosts_1_primitive_var;
                 j < interior_dim_1 + 1 + num_ghosts_1_primitive_var;
                 j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                        (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                            (ghostcell_dim_1_primitive_var + 1);
                    
                    if (rho[idx_face] > double(0))
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                    
                    if (V[d_num_species + d_dim.getValue()][idx_face] > double(0))
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                }
            }
        }
        
        /*
         * Check if primitive variables in the z-direction are bounded.
         */
        
        are_bounded = bounded_flag->getPointer(2);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V[ei] = primitive_variables[ei]->getPointer(2);
        }
        
        rho = data_density->getPointer(2);
        
        // Compute density.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = -num_ghosts_2_primitive_var;
                 k < interior_dim_2 + 1 + num_ghosts_2_primitive_var;
                 k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx_face = (i + num_ghosts_0_primitive_var) +
                            (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                            (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                                ghostcell_dim_1_primitive_var;
                        
                        rho[idx_face] += V[si][idx_face];
                    }
                }
            }
        }
        
        // Check if mass fractions are bounded.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = -num_ghosts_2_primitive_var;
                 k < interior_dim_2 + 1 + num_ghosts_2_primitive_var;
                 k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx_face = (i + num_ghosts_0_primitive_var) +
                            (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                            (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                                ghostcell_dim_1_primitive_var;
                        
                        const double Y = V[si][idx_face]/rho[idx_face];
                        
                        if (Y > d_Y_bound_lo && Y < d_Y_bound_up)
                        {
                            are_bounded[idx_face] &= 1;
                        }
                        else
                        {
                            are_bounded[idx_face] &= 0;
                        }
                    }
                }
            }
        }
        
        // Check if density and pressure are bounded.
        for (int k = -num_ghosts_2_primitive_var;
             k < interior_dim_2 + 1 + num_ghosts_2_primitive_var;
             k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                        (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                            ghostcell_dim_1_primitive_var;
                    
                    if (rho[idx_face] > double(0))
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                    
                    if (V[d_num_species + d_dim.getValue()][idx_face] > double(0))
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                }
            }
        }
    }
}


/*
 * Convert vector of pointers of conservative cell data to vectors of pointers of primitive cell data.
 */
void
FlowModelFourEqnConservative::convertLocalCellDataPointersConservativeVariablesToPrimitiveVariables(
    const std::vector<const double*>& conservative_variables,
    const std::vector<double*>& primitive_variables)
{
    const std::vector<const double*>& Q = conservative_variables;
    const std::vector<double*>&       V = primitive_variables;
    
    // Compute the mixture density.
    std::vector<const double*> rho_Y_ptr;
    rho_Y_ptr.reserve(d_num_species);
    for (int si = 0; si < d_num_species; si++)
    {
        rho_Y_ptr.push_back(Q[si]);
    }
    const double rho = d_equation_of_state_mixing_rules->getMixtureDensity(
        rho_Y_ptr);
    
    /*
     * Compute the internal energy.
     */
    double epsilon = double(0);
    if (d_dim == tbox::Dimension(1))
    {
        epsilon = ((*Q[d_num_species + d_dim.getValue()]) -
            double(1)/double(2)*((*Q[d_num_species])*(*Q[d_num_species]))/rho)/rho;
    }
    else if (d_dim == tbox::Dimension(2))
    {
        epsilon = ((*Q[d_num_species + d_dim.getValue()]) -
            double(1)/double(2)*((*Q[d_num_species])*(*Q[d_num_species]) +
            (*Q[d_num_species + 1])*(*Q[d_num_species + 1]))/rho)/rho;
    }
    else if (d_dim == tbox::Dimension(3))
    {
        epsilon = ((*Q[d_num_species + d_dim.getValue()]) -
            double(1)/double(2)*((*Q[d_num_species])*(*Q[d_num_species]) +
            (*Q[d_num_species + 1])*(*Q[d_num_species + 1]) +
            (*Q[d_num_species + 2])*(*Q[d_num_species + 2]))/rho)/rho;
    }
    
    /*
     * Compute the mass fractions.
     */
    double Y[d_num_species];
    for (int si = 0; si < d_num_species; si++)
    {
        Y[si] = (*Q[si])/rho;
    }
    
    /*
     * Get the pointers to the mass fractions.
     */
    std::vector<const double*> Y_ptr;
    Y_ptr.reserve(d_num_species);
    for (int si = 0; si < d_num_species; si++)
    {
        Y_ptr.push_back(&Y[si]);
    }
    
    // Compute the pressure.
    const double p = d_equation_of_state_mixing_rules->getPressure(
        &rho,
        &epsilon,
        Y_ptr);
    
    // Convert the conservative variables to primitive variables.
    for (int si = 0; si < d_num_species; si++)
    {
        *V[si] = *Q[si];
    }
    for (int di = 0; di < d_dim.getValue(); di++)
    {
        *V[d_num_species + di] = (*Q[d_num_species + di])/rho;
    }
    *V[d_num_species + d_dim.getValue()] = p;
}


/*
 * Convert vector of pointers of primitive cell data to vectors of pointers of conservative cell data.
 */
void
FlowModelFourEqnConservative::convertLocalCellDataPointersPrimitiveVariablesToConservativeVariables(
    const std::vector<const double*>& primitive_variables,
    const std::vector<double*>& conservative_variables)
{
    const std::vector<const double*>& V = primitive_variables;
    const std::vector<double*>&       Q = conservative_variables;
    
    // Compute the mixture density.
    std::vector<const double*> rho_Y_ptr;
    rho_Y_ptr.reserve(d_num_species);
    for (int si = 0; si < d_num_species; si++)
    {
        rho_Y_ptr.push_back(V[si]);
    }
    const double rho = d_equation_of_state_mixing_rules->getMixtureDensity(
        rho_Y_ptr);
    
    /*
     * Compute the mass fractions.
     */
    double Y[d_num_species];
    for (int si = 0; si < d_num_species; si++)
    {
        Y[si] = (*V[si])/rho;
    }
    
    /*
     * Get the pointers to the mass fractions.
     */
    std::vector<const double*> Y_ptr;
    Y_ptr.reserve(d_num_species);
    for (int si = 0; si < d_num_species; si++)
    {
        Y_ptr.push_back(&Y[si]);
    }
    
    // Compute the total energy.
    const double epsilon = d_equation_of_state_mixing_rules->getInternalEnergy(
        &rho,
        V[d_num_species + d_dim.getValue()],
        Y_ptr);
    
    const double E = epsilon + double(1)/double(2)*rho*((*Q[d_num_species])*(*Q[d_num_species]) +
        (*Q[d_num_species + 1])*(*Q[d_num_species + 1]) +
        (*Q[d_num_species + 2])*(*Q[d_num_species + 2]));
    
    // Convert the primitive variables to conservative variables.
    for (int si = 0; si < d_num_species; si++)
    {
        *Q[si] = *V[si];
    }
    for (int di = 0; di < d_dim.getValue(); di++)
    {
        *Q[d_num_species + di] = rho*(*V[d_num_species + di]);
    }
    *Q[d_num_species + d_dim.getValue()] = E;
}


/*
 * Get the variables for the derivatives in the diffusive fluxes.
 */
void
FlowModelFourEqnConservative::getDiffusiveFluxVariablesForDerivative(
    std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& derivative_var_data,
    std::vector<std::vector<int> >& derivative_var_component_idx,
    const DIRECTION::TYPE& flux_direction,
    const DIRECTION::TYPE& derivative_direction)
{
    // Create empty box.
    const hier::Box empty_box(d_dim);
    
    derivative_var_data.resize(d_num_eqn);
    derivative_var_component_idx.resize(d_num_eqn);
    
    if (!d_data_mass_fractions)
    {
        computeGlobalCellDataMassFractionsWithDensity(empty_box);
    }
    
    if (!d_data_velocity)
    {
        computeGlobalCellDataVelocityWithDensity(empty_box);
    }
    
    if (!d_data_temperature)
    {
        computeGlobalCellDataTemperatureWithDensityMassFractionsAndPressure(empty_box);
    }
    
    if (d_dim == tbox::Dimension(1))
    {
        switch (flux_direction)
        {
            case DIRECTION::X_DIRECTION:
            {
                switch (derivative_direction)
                {
                    case DIRECTION::X_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            derivative_var_data[si].resize(d_num_species + 1);
                            derivative_var_component_idx[si].resize(d_num_species + 1);
                            
                            derivative_var_data[si][0] = d_data_mass_fractions;
                            derivative_var_component_idx[si][0] = si;
                            
                            for (int sj = 0; sj < d_num_species; sj++)
                            {
                                derivative_var_data[si][1 + sj] = d_data_mass_fractions;
                                derivative_var_component_idx[si][1 + sj] = sj;
                            }
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[d_num_species].resize(1);
                        derivative_var_component_idx[d_num_species].resize(1);
                        
                        // Variable u.
                        derivative_var_data[d_num_species][0] = d_data_velocity;
                        derivative_var_component_idx[d_num_species][0] = 0;
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[d_num_species + 1].resize(
                            2 + d_num_species*(d_num_species + 1));
                        derivative_var_component_idx[d_num_species + 1].resize(
                            2 + d_num_species*(d_num_species + 1));
                        
                        // Variable u.
                        derivative_var_data[d_num_species + 1][0] = d_data_velocity;
                        derivative_var_component_idx[d_num_species + 1][0] = 0;
                        
                        // Variable T.
                        derivative_var_data[d_num_species + 1][1] = d_data_temperature;
                        derivative_var_component_idx[d_num_species + 1][1] = 0;
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            derivative_var_data[d_num_species + 1][2 + si*(d_num_species + 1)] =
                                d_data_mass_fractions;
                            derivative_var_component_idx[d_num_species + 1][2 + si*(d_num_species + 1)] =
                                si;
                            
                            for (int sj = 0; sj < d_num_species; sj++)
                            {
                                derivative_var_data[d_num_species + 1][3 + si*(d_num_species + 1) + sj] =
                                    d_data_mass_fractions;
                                derivative_var_component_idx[d_num_species + 1][3 + si*(d_num_species + 1) + sj] =
                                    sj;
                            }
                        }
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelFourEqnConservative::getDiffusiveFluxVariablesForDerivative()\n"
                            << "There are only x-direction for one-dimensional problem."
                            << std::endl);
                    }
                }
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelFourEqnConservative::getDiffusiveFluxVariablesForDerivative()\n"
                    << "There are only x-direction for one-dimensional problem."
                    << std::endl);
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        switch (flux_direction)
        {
            case DIRECTION::X_DIRECTION:
            {
                switch (derivative_direction)
                {
                    case DIRECTION::X_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            derivative_var_data[si].resize(d_num_species + 1);
                            derivative_var_component_idx[si].resize(d_num_species + 1);
                            
                            derivative_var_data[si][0] = d_data_mass_fractions;
                            derivative_var_component_idx[si][0] = si;
                            
                            for (int sj = 0; sj < d_num_species; sj++)
                            {
                                derivative_var_data[si][1 + sj] = d_data_mass_fractions;
                                derivative_var_component_idx[si][1 + sj] = sj;
                            }
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[d_num_species].resize(1);
                        derivative_var_component_idx[d_num_species].resize(1);
                        
                        // Variable u.
                        derivative_var_data[d_num_species][0] = d_data_velocity;
                        derivative_var_component_idx[d_num_species][0] = 0;
                        
                        derivative_var_data[d_num_species + 1].resize(1);
                        derivative_var_component_idx[d_num_species + 1].resize(1);
                        
                        // Variable v.
                        derivative_var_data[d_num_species + 1][0] = d_data_velocity;
                        derivative_var_component_idx[d_num_species + 1][0] = 1;
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[d_num_species + 2].resize(
                            3 + d_num_species*(d_num_species + 1));
                        derivative_var_component_idx[d_num_species + 2].resize(
                            3 + d_num_species*(d_num_species + 1));
                        
                        // Variable u.
                        derivative_var_data[d_num_species + 2][0] = d_data_velocity;
                        derivative_var_component_idx[d_num_species + 2][0] = 0;
                        
                        // Variable v.
                        derivative_var_data[d_num_species + 2][1] = d_data_velocity;
                        derivative_var_component_idx[d_num_species + 2][1] = 1;
                        
                        // Variable T.
                        derivative_var_data[d_num_species + 2][2] = d_data_temperature;
                        derivative_var_component_idx[d_num_species + 2][2] = 0;
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            derivative_var_data[d_num_species + 2][3 + si*(d_num_species + 1)] =
                                d_data_mass_fractions;
                            derivative_var_component_idx[d_num_species + 2][3 + si*(d_num_species + 1)] =
                                si;
                            
                            for (int sj = 0; sj < d_num_species; sj++)
                            {
                                derivative_var_data[d_num_species + 2][4 + si*(d_num_species + 1) + sj] =
                                    d_data_mass_fractions;
                                derivative_var_component_idx[d_num_species + 2][4 + si*(d_num_species + 1) + sj] =
                                    sj;
                            }
                        }
                        
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            derivative_var_data[si].resize(0);
                            derivative_var_component_idx[si].resize(0);
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[d_num_species].resize(1);
                        derivative_var_component_idx[d_num_species].resize(1);
                        
                        // Variable v.
                        derivative_var_data[d_num_species][0] = d_data_velocity;
                        derivative_var_component_idx[d_num_species][0] = 1;
                        
                        derivative_var_data[d_num_species + 1].resize(1);
                        derivative_var_component_idx[d_num_species + 1].resize(1);
                        
                        // Variable u.
                        derivative_var_data[d_num_species + 1][0] = d_data_velocity;
                        derivative_var_component_idx[d_num_species + 1][0] = 0;
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[d_num_species + 2].resize(2);
                        derivative_var_component_idx[d_num_species + 2].resize(2);
                        
                        // Variable u.
                        derivative_var_data[d_num_species + 2][0] = d_data_velocity;
                        derivative_var_component_idx[d_num_species + 2][0] = 0;
                        
                        // Variable v.
                        derivative_var_data[d_num_species + 2][1] = d_data_velocity;
                        derivative_var_component_idx[d_num_species + 2][1] = 1;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelFourEqnConservative::getDiffusiveFluxVariablesForDerivative()\n"
                            << "There are only x-direction and y-direction for two-dimensional problem."
                            << std::endl);
                    }
                }
                
                break;
            }
            case DIRECTION::Y_DIRECTION:
            {
                switch (derivative_direction)
                {
                    case DIRECTION::X_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            derivative_var_data[si].resize(0);
                            derivative_var_component_idx[si].resize(0);
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[d_num_species].resize(1);
                        derivative_var_component_idx[d_num_species].resize(1);
                        
                        // Variable v.
                        derivative_var_data[d_num_species][0] = d_data_velocity;
                        derivative_var_component_idx[d_num_species][0] = 1;
                        
                        derivative_var_data[d_num_species + 1].resize(1);
                        derivative_var_component_idx[d_num_species + 1].resize(1);
                        
                        // Variable u.
                        derivative_var_data[d_num_species + 1][0] = d_data_velocity;
                        derivative_var_component_idx[d_num_species + 1][0] = 0;
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[d_num_species + 2].resize(2);
                        derivative_var_component_idx[d_num_species + 2].resize(2);
                        
                        // Variable u.
                        derivative_var_data[d_num_species + 2][0] = d_data_velocity;
                        derivative_var_component_idx[d_num_species + 2][0] = 0;
                        
                        // Variable v.
                        derivative_var_data[d_num_species + 2][1] = d_data_velocity;
                        derivative_var_component_idx[d_num_species + 2][1] = 1;
                        
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            derivative_var_data[si].resize(d_num_species + 1);
                            derivative_var_component_idx[si].resize(d_num_species + 1);
                            
                            derivative_var_data[si][0] = d_data_mass_fractions;
                            derivative_var_component_idx[si][0] = si;
                            
                            for (int sj = 0; sj < d_num_species; sj++)
                            {
                                derivative_var_data[si][1 + sj] = d_data_mass_fractions;
                                derivative_var_component_idx[si][1 + sj] = sj;
                            }
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[d_num_species].resize(1);
                        derivative_var_component_idx[d_num_species].resize(1);
                        
                        // Variable u.
                        derivative_var_data[d_num_species][0] = d_data_velocity;
                        derivative_var_component_idx[d_num_species][0] = 0;
                        
                        derivative_var_data[d_num_species + 1].resize(1);
                        derivative_var_component_idx[d_num_species + 1].resize(1);
                        
                        // Variable v.
                        derivative_var_data[d_num_species + 1][0] = d_data_velocity;
                        derivative_var_component_idx[d_num_species + 1][0] = 1;
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[d_num_species + 2].resize(
                            3 + d_num_species*(d_num_species + 1));
                        derivative_var_component_idx[d_num_species + 2].resize(
                            3 + d_num_species*(d_num_species + 1));
                        
                        // Variable u.
                        derivative_var_data[d_num_species + 2][0] = d_data_velocity;
                        derivative_var_component_idx[d_num_species + 2][0] = 0;
                        
                        // Variable v.
                        derivative_var_data[d_num_species + 2][1] = d_data_velocity;
                        derivative_var_component_idx[d_num_species + 2][1] = 1;
                        
                        // Variable T.
                        derivative_var_data[d_num_species + 2][2] = d_data_temperature;
                        derivative_var_component_idx[d_num_species + 2][2] = 0;
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            derivative_var_data[d_num_species + 2][3 + si*(d_num_species + 1)] =
                                d_data_mass_fractions;
                            derivative_var_component_idx[d_num_species + 2][3 + si*(d_num_species + 1)] =
                                si;
                            
                            for (int sj = 0; sj < d_num_species; sj++)
                            {
                                derivative_var_data[d_num_species + 2][4 + si*(d_num_species + 1) + sj] =
                                    d_data_mass_fractions;
                                derivative_var_component_idx[d_num_species + 2][4 + si*(d_num_species + 1) + sj] =
                                    sj;
                            }
                        }
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelFourEqnConservative::getDiffusiveFluxVariablesForDerivative()\n"
                            << "There are only x-direction and y-direction for two-dimensional problem."
                            << std::endl);
                    }
                }
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelFourEqnConservative::getDiffusiveFluxVariablesForDerivative()\n"
                    << "There are only x-direction and y-direction for two-dimensional problem."
                    << std::endl);
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        switch (flux_direction)
        {
            case DIRECTION::X_DIRECTION:
            {
                switch (derivative_direction)
                {
                    case DIRECTION::X_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            derivative_var_data[si].resize(d_num_species + 1);
                            derivative_var_component_idx[si].resize(d_num_species + 1);
                            
                            derivative_var_data[si][0] = d_data_mass_fractions;
                            derivative_var_component_idx[si][0] = si;
                            
                            for (int sj = 0; sj < d_num_species; sj++)
                            {
                                derivative_var_data[si][1 + sj] = d_data_mass_fractions;
                                derivative_var_component_idx[si][1 + sj] = sj;
                            }
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[d_num_species].resize(1);
                        derivative_var_component_idx[d_num_species].resize(1);
                        
                        // Variable u.
                        derivative_var_data[d_num_species][0] = d_data_velocity;
                        derivative_var_component_idx[d_num_species][0] = 0;
                        
                        derivative_var_data[d_num_species + 1].resize(1);
                        derivative_var_component_idx[d_num_species + 1].resize(1);
                        
                        // Variable v.
                        derivative_var_data[d_num_species + 1][0] = d_data_velocity;
                        derivative_var_component_idx[d_num_species + 1][0] = 1;
                        
                        derivative_var_data[d_num_species + 2].resize(1);
                        derivative_var_component_idx[d_num_species + 2].resize(1);
                        
                        // Variable w.
                        derivative_var_data[d_num_species + 2][0] = d_data_velocity;
                        derivative_var_component_idx[d_num_species + 2][0] = 2;
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[d_num_species + 3].resize(
                            4 + d_num_species*(d_num_species + 1));
                        derivative_var_component_idx[d_num_species + 3].resize(
                            4 + d_num_species*(d_num_species + 1));
                        
                        // Variable u.
                        derivative_var_data[d_num_species + 3][0] = d_data_velocity;
                        derivative_var_component_idx[d_num_species + 3][0] = 0;
                        
                        // Variable v.
                        derivative_var_data[d_num_species + 3][1] = d_data_velocity;
                        derivative_var_component_idx[d_num_species + 3][1] = 1;
                        
                        // Variable w.
                        derivative_var_data[d_num_species + 3][2] = d_data_velocity;
                        derivative_var_component_idx[d_num_species + 3][2] = 2;
                        
                        // Variable T.
                        derivative_var_data[d_num_species + 3][3] = d_data_temperature;
                        derivative_var_component_idx[d_num_species + 3][3] = 0;
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            derivative_var_data[d_num_species + 3][4 + si*(d_num_species + 1)] =
                                d_data_mass_fractions;
                            derivative_var_component_idx[d_num_species + 3][4 + si*(d_num_species + 1)] =
                                si;
                            
                            for (int sj = 0; sj < d_num_species; sj++)
                            {
                                derivative_var_data[d_num_species + 3][5 + si*(d_num_species + 1) + sj] =
                                    d_data_mass_fractions;
                                derivative_var_component_idx[d_num_species + 3][5 + si*(d_num_species + 1) + sj] =
                                    sj;
                            }
                        }
                        
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            derivative_var_data[si].resize(0);
                            derivative_var_component_idx[si].resize(0);
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[d_num_species].resize(1);
                        derivative_var_component_idx[d_num_species].resize(1);
                        
                        // Variable v.
                        derivative_var_data[d_num_species][0] = d_data_velocity;
                        derivative_var_component_idx[d_num_species][0] = 1;
                        
                        derivative_var_data[d_num_species + 1].resize(1);
                        derivative_var_component_idx[d_num_species + 1].resize(1);
                        
                        // Variable u.
                        derivative_var_data[d_num_species + 1][0] = d_data_velocity;
                        derivative_var_component_idx[d_num_species + 1][0] = 0;
                        
                        derivative_var_data[d_num_species + 2].resize(0);
                        derivative_var_component_idx[d_num_species + 2].resize(0);
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[d_num_species + 3].resize(2);
                        derivative_var_component_idx[d_num_species + 3].resize(2);
                        
                        // Variable u.
                        derivative_var_data[d_num_species + 3][0] = d_data_velocity;
                        derivative_var_component_idx[d_num_species + 3][0] = 0;
                        
                        // Variable v.
                        derivative_var_data[d_num_species + 3][1] = d_data_velocity;
                        derivative_var_component_idx[d_num_species + 3][1] = 1;
                        
                        break;
                    }
                    case DIRECTION::Z_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            derivative_var_data[si].resize(0);
                            derivative_var_component_idx[si].resize(0);
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[d_num_species].resize(1);
                        derivative_var_component_idx[d_num_species].resize(1);
                        
                        // Variable w.
                        derivative_var_data[d_num_species][0] = d_data_velocity;
                        derivative_var_component_idx[d_num_species][0] = 2;
                        
                        derivative_var_data[d_num_species + 1].resize(0);
                        derivative_var_component_idx[d_num_species + 1].resize(0);
                        
                        derivative_var_data[d_num_species + 2].resize(1);
                        derivative_var_component_idx[d_num_species + 2].resize(1);
                        
                        // Variable u.
                        derivative_var_data[d_num_species + 2][0] = d_data_velocity;
                        derivative_var_component_idx[d_num_species + 2][0] = 0;
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[d_num_species + 3].resize(2);
                        derivative_var_component_idx[d_num_species + 3].resize(2);
                        
                        // Variable u.
                        derivative_var_data[d_num_species + 3][0] = d_data_velocity;
                        derivative_var_component_idx[d_num_species + 3][0] = 0;
                        
                        // Variable w.
                        derivative_var_data[d_num_species + 3][1] = d_data_velocity;
                        derivative_var_component_idx[d_num_species + 3][1] = 2;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelFourEqnConservative::getDiffusiveFluxVariablesForDerivative()\n"
                            << "There are only x-direction, y-direction and z-direction for three-dimensional problem."
                            << std::endl);
                    }
                }
                
                break;
            }
            case DIRECTION::Y_DIRECTION:
            {
                switch (derivative_direction)
                {
                    case DIRECTION::X_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            derivative_var_data[si].resize(0);
                            derivative_var_component_idx[si].resize(0);
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[d_num_species].resize(1);
                        derivative_var_component_idx[d_num_species].resize(1);
                        
                        // Variable v.
                        derivative_var_data[d_num_species][0] = d_data_velocity;
                        derivative_var_component_idx[d_num_species][0] = 1;
                        
                        derivative_var_data[d_num_species + 1].resize(1);
                        derivative_var_component_idx[d_num_species + 1].resize(1);
                        
                        // Variable u.
                        derivative_var_data[d_num_species + 1][0] = d_data_velocity;
                        derivative_var_component_idx[d_num_species + 1][0] = 0;
                        
                        derivative_var_data[d_num_species + 2].resize(0);
                        derivative_var_component_idx[d_num_species + 2].resize(0);
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[d_num_species + 3].resize(2);
                        derivative_var_component_idx[d_num_species + 3].resize(2);
                        
                        // Variable u.
                        derivative_var_data[d_num_species + 3][0] = d_data_velocity;
                        derivative_var_component_idx[d_num_species + 3][0] = 0;
                        
                        // Variable v.
                        derivative_var_data[d_num_species + 3][1] = d_data_velocity;
                        derivative_var_component_idx[d_num_species + 3][1] = 1;
                        
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            derivative_var_data[si].resize(d_num_species + 1);
                            derivative_var_component_idx[si].resize(d_num_species + 1);
                            
                            derivative_var_data[si][0] = d_data_mass_fractions;
                            derivative_var_component_idx[si][0] = si;
                            
                            for (int sj = 0; sj < d_num_species; sj++)
                            {
                                derivative_var_data[si][1 + sj] = d_data_mass_fractions;
                                derivative_var_component_idx[si][1 + sj] = sj;
                            }
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[d_num_species].resize(1);
                        derivative_var_component_idx[d_num_species].resize(1);
                        
                        // Variable u.
                        derivative_var_data[d_num_species][0] = d_data_velocity;
                        derivative_var_component_idx[d_num_species][0] = 0;
                        
                        derivative_var_data[d_num_species + 1].resize(1);
                        derivative_var_component_idx[d_num_species + 1].resize(1);
                        
                        // Variable v.
                        derivative_var_data[d_num_species + 1][0] = d_data_velocity;
                        derivative_var_component_idx[d_num_species + 1][0] = 1;
                        
                        derivative_var_data[d_num_species + 2].resize(1);
                        derivative_var_component_idx[d_num_species + 2].resize(1);
                        
                        // Variable w.
                        derivative_var_data[d_num_species + 2][0] = d_data_velocity;
                        derivative_var_component_idx[d_num_species + 2][0] = 2;
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[d_num_species + 3].resize(
                            4 + d_num_species*(d_num_species + 1));
                        derivative_var_component_idx[d_num_species + 3].resize(
                            4 + d_num_species*(d_num_species + 1));
                        
                        // Variable u.
                        derivative_var_data[d_num_species + 3][0] = d_data_velocity;
                        derivative_var_component_idx[d_num_species + 3][0] = 0;
                        
                        // Variable v.
                        derivative_var_data[d_num_species + 3][1] = d_data_velocity;
                        derivative_var_component_idx[d_num_species + 3][1] = 1;
                        
                        // Variable w.
                        derivative_var_data[d_num_species + 3][2] = d_data_velocity;
                        derivative_var_component_idx[d_num_species + 3][2] = 2;
                        
                        // Variable T.
                        derivative_var_data[d_num_species + 3][3] = d_data_temperature;
                        derivative_var_component_idx[d_num_species + 3][3] = 0;
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            derivative_var_data[d_num_species + 3][4 + si*(d_num_species + 1)] =
                                d_data_mass_fractions;
                            derivative_var_component_idx[d_num_species + 3][4 + si*(d_num_species + 1)] =
                                si;
                            
                            for (int sj = 0; sj < d_num_species; sj++)
                            {
                                derivative_var_data[d_num_species + 3][5 + si*(d_num_species + 1) + sj] =
                                    d_data_mass_fractions;
                                derivative_var_component_idx[d_num_species + 3][5 + si*(d_num_species + 1) + sj] =
                                    sj;
                            }
                        }
                        
                        break;
                    }
                    case DIRECTION::Z_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            derivative_var_data[si].resize(0);
                            derivative_var_component_idx[si].resize(0);
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[d_num_species].resize(0);
                        derivative_var_component_idx[d_num_species].resize(0);
                        
                        derivative_var_data[d_num_species + 1].resize(1);
                        derivative_var_component_idx[d_num_species + 1].resize(1);
                        
                        // Variable w.
                        derivative_var_data[d_num_species + 1][0] = d_data_velocity;
                        derivative_var_component_idx[d_num_species + 1][0] = 2;
                        
                        derivative_var_data[d_num_species + 2].resize(1);
                        derivative_var_component_idx[d_num_species + 2].resize(1);
                        
                        // Variable v.
                        derivative_var_data[d_num_species + 2][0] = d_data_velocity;
                        derivative_var_component_idx[d_num_species + 2][0] = 1;
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[d_num_species + 3].resize(2);
                        derivative_var_component_idx[d_num_species + 3].resize(2);
                        
                        // Variable v.
                        derivative_var_data[d_num_species + 3][0] = d_data_velocity;
                        derivative_var_component_idx[d_num_species + 3][0] = 1;
                        
                        // Variable w.
                        derivative_var_data[d_num_species + 3][1] = d_data_velocity;
                        derivative_var_component_idx[d_num_species + 3][1] = 2;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelFourEqnConservative::getDiffusiveFluxVariablesForDerivative()\n"
                            << "There are only x-direction, y-direction and z-direction for three-dimensional problem."
                            << std::endl);
                    }
                }
                
                break;
            }
            case DIRECTION::Z_DIRECTION:
            {
                switch (derivative_direction)
                {
                    case DIRECTION::X_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            derivative_var_data[si].resize(0);
                            derivative_var_component_idx[si].resize(0);
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[d_num_species].resize(1);
                        derivative_var_component_idx[d_num_species].resize(1);
                        
                        // Variable w.
                        derivative_var_data[d_num_species][0] = d_data_velocity;
                        derivative_var_component_idx[d_num_species][0] = 2;
                        
                        derivative_var_data[d_num_species + 1].resize(0);
                        derivative_var_component_idx[d_num_species + 1].resize(0);
                        
                        derivative_var_data[d_num_species + 2].resize(1);
                        derivative_var_component_idx[d_num_species + 2].resize(1);
                        
                        // Variable u.
                        derivative_var_data[d_num_species + 2][0] = d_data_velocity;
                        derivative_var_component_idx[d_num_species + 2][0] = 0;
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[d_num_species + 3].resize(2);
                        derivative_var_component_idx[d_num_species + 3].resize(2);
                        
                        // Variable u.
                        derivative_var_data[d_num_species + 3][0] = d_data_velocity;
                        derivative_var_component_idx[d_num_species + 3][0] = 0;
                        
                        // Variable w.
                        derivative_var_data[d_num_species + 3][1] = d_data_velocity;
                        derivative_var_component_idx[d_num_species + 3][1] = 2;
                        
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            derivative_var_data[si].resize(0);
                            derivative_var_component_idx[si].resize(0);
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[d_num_species].resize(0);
                        derivative_var_component_idx[d_num_species].resize(0);
                        
                        derivative_var_data[d_num_species + 1].resize(1);
                        derivative_var_component_idx[d_num_species + 1].resize(1);
                        
                        // Variable w.
                        derivative_var_data[d_num_species + 1][0] = d_data_velocity;
                        derivative_var_component_idx[d_num_species + 1][0] = 2;
                        
                        derivative_var_data[d_num_species + 2].resize(1);
                        derivative_var_component_idx[d_num_species + 2].resize(1);
                        
                        // Variable v.
                        derivative_var_data[d_num_species + 2][0] = d_data_velocity;
                        derivative_var_component_idx[d_num_species + 2][0] = 1;
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[d_num_species + 3].resize(2);
                        derivative_var_component_idx[d_num_species + 3].resize(2);
                        
                        // Variable v.
                        derivative_var_data[d_num_species + 3][0] = d_data_velocity;
                        derivative_var_component_idx[d_num_species + 3][0] = 1;
                        
                        // Variable w.
                        derivative_var_data[d_num_species + 3][1] = d_data_velocity;
                        derivative_var_component_idx[d_num_species + 3][1] = 2;
                        
                        break;
                    }
                    case DIRECTION::Z_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            derivative_var_data[si].resize(d_num_species + 1);
                            derivative_var_component_idx[si].resize(d_num_species + 1);
                            
                            derivative_var_data[si][0] = d_data_mass_fractions;
                            derivative_var_component_idx[si][0] = si;
                            
                            for (int sj = 0; sj < d_num_species; sj++)
                            {
                                derivative_var_data[si][1 + sj] = d_data_mass_fractions;
                                derivative_var_component_idx[si][1 + sj] = sj;
                            }
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[d_num_species].resize(1);
                        derivative_var_component_idx[d_num_species].resize(1);
                        
                        // Variable u.
                        derivative_var_data[d_num_species][0] = d_data_velocity;
                        derivative_var_component_idx[d_num_species][0] = 0;
                        
                        derivative_var_data[d_num_species + 1].resize(1);
                        derivative_var_component_idx[d_num_species + 1].resize(1);
                        
                        // Variable v.
                        derivative_var_data[d_num_species + 1][0] = d_data_velocity;
                        derivative_var_component_idx[d_num_species + 1][0] = 1;
                        
                        derivative_var_data[d_num_species + 2].resize(1);
                        derivative_var_component_idx[d_num_species + 2].resize(1);
                        
                        // Variable w.
                        derivative_var_data[d_num_species + 2][0] = d_data_velocity;
                        derivative_var_component_idx[d_num_species + 2][0] = 2;
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[d_num_species + 3].resize(
                            4 + d_num_species*(d_num_species + 1));
                        derivative_var_component_idx[d_num_species + 3].resize(
                            4 + d_num_species*(d_num_species + 1));
                        
                        // Variable u.
                        derivative_var_data[d_num_species + 3][0] = d_data_velocity;
                        derivative_var_component_idx[d_num_species + 3][0] = 0;
                        
                        // Variable v.
                        derivative_var_data[d_num_species + 3][1] = d_data_velocity;
                        derivative_var_component_idx[d_num_species + 3][1] = 1;
                        
                        // Variable w.
                        derivative_var_data[d_num_species + 3][2] = d_data_velocity;
                        derivative_var_component_idx[d_num_species + 3][2] = 2;
                        
                        // Variable T.
                        derivative_var_data[d_num_species + 3][3] = d_data_temperature;
                        derivative_var_component_idx[d_num_species + 3][3] = 0;
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            derivative_var_data[d_num_species + 3][4 + si*(d_num_species + 1)] =
                                d_data_mass_fractions;
                            derivative_var_component_idx[d_num_species + 3][4 + si*(d_num_species + 1)] =
                                si;
                            
                            for (int sj = 0; sj < d_num_species; sj++)
                            {
                                derivative_var_data[d_num_species + 3][5 + si*(d_num_species + 1) + sj] =
                                    d_data_mass_fractions;
                                derivative_var_component_idx[d_num_species + 3][5 + si*(d_num_species + 1) + sj] =
                                    sj;
                            }
                        }
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelFourEqnConservative::getDiffusiveFluxVariablesForDerivative()\n"
                            << "There are only x-direction, y-direction and z-direction for three-dimensional problem."
                            << std::endl);
                    }
                }
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelFourEqnConservative::getDiffusiveFluxVariablesForDerivative()\n"
                    << "There are only x-direction, y-direction and z-direction for three-dimensional problem."
                    << std::endl);
            }
        }
    }
    
    d_global_derived_cell_data_computed = true;
}


/*
 * Get the diffusivities in the diffusive flux.
 */
void
FlowModelFourEqnConservative::getDiffusiveFluxDiffusivities(
    std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& diffusivities_data,
    std::vector<std::vector<int> >& diffusivities_component_idx,
    const DIRECTION::TYPE& flux_direction,
    const DIRECTION::TYPE& derivative_direction)
{
    // Create empty box.
    const hier::Box empty_box(d_dim);
    
    if (!d_equation_of_mass_diffusivity_mixing_rules ||
        !d_equation_of_shear_viscosity_mixing_rules ||
        !d_equation_of_bulk_viscosity_mixing_rules ||
        !d_equation_of_thermal_conductivity_mixing_rules)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::getDiffusiveFluxDiffusivities()\n"
            << "Either mixing rule of mass diffusivity, shear viscosity, bulk viscosity or"
            << " thermal conductivity is not initialized."
            << std::endl);
    }
    
    diffusivities_data.resize(d_num_eqn);
    diffusivities_component_idx.resize(d_num_eqn);
        
    if (!d_data_diffusivities)
    {
        if (!d_data_density)
        {
            computeGlobalCellDataDensity(empty_box);
        }
        
        if (!d_data_mass_fractions)
        {
            computeGlobalCellDataMassFractionsWithDensity(empty_box);
        }
        
        if (!d_data_velocity)
        {
            computeGlobalCellDataVelocityWithDensity(empty_box);
        }
        
        if (!d_data_pressure)
        {
            computeGlobalCellDataPressureWithDensityMassFractionsAndInternalEnergy(empty_box);
        }
        
        if (!d_data_temperature)
        {
            computeGlobalCellDataTemperatureWithDensityMassFractionsAndPressure(empty_box);
        }
        
        /*
         * Create temporary cell data of mass diffusivities, shear viscosity, bulk viscosity and
         * thermal conductivity.
         */
        
        boost::shared_ptr<pdat::CellData<double> > data_mass_diffusivities(
            new pdat::CellData<double>(d_interior_box, d_num_species, d_num_subghosts_diffusivities));
        
        boost::shared_ptr<pdat::CellData<double> > data_shear_viscosity(
            new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_diffusivities));
        
        boost::shared_ptr<pdat::CellData<double> > data_bulk_viscosity(
            new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_diffusivities));
        
        boost::shared_ptr<pdat::CellData<double> > data_thermal_conductivity(
            new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_diffusivities));
        
        // Get the pointers to the cell data of density, mass fractions, mass diffusivities, shear
        // viscosity, bulk viscosity and thermal conductivity.
        double* rho = d_data_density->getPointer(0);
        std::vector<double*> Y;
        Y.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            Y.push_back(d_data_mass_fractions->getPointer(si));
        }
        std::vector<double*> D;
        D.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            D.push_back(data_mass_diffusivities->getPointer(si));
        }
        double* mu    = data_shear_viscosity->getPointer(0);
        double* mu_v  = data_bulk_viscosity->getPointer(0);
        double* kappa = data_thermal_conductivity->getPointer(0);
        
        // Compute the mass diffusivity fields.
        d_equation_of_mass_diffusivity_mixing_rules->computeMassDiffusivities(
            data_mass_diffusivities,
            d_data_pressure,
            d_data_temperature,
            d_data_mass_fractions,
            d_subghost_box_diffusivities);
        
        // Compute the shear viscosity field.
        d_equation_of_shear_viscosity_mixing_rules->computeShearViscosity(
            data_shear_viscosity,
            d_data_pressure,
            d_data_temperature,
            d_data_mass_fractions,
            d_subghost_box_diffusivities);
        
        // Compute the bulk viscosity field.
        d_equation_of_bulk_viscosity_mixing_rules->computeBulkViscosity(
            data_bulk_viscosity,
            d_data_pressure,
            d_data_temperature,
            d_data_mass_fractions,
            d_subghost_box_diffusivities);
        
        // Compute the thermal conductivity field.
        d_equation_of_thermal_conductivity_mixing_rules->computeThermalConductivity(
            data_thermal_conductivity,
            d_data_pressure,
            d_data_temperature,
            d_data_mass_fractions,
            d_subghost_box_diffusivities);
        
        /*
         * Compute enthalpy of each species.
         */
        
        std::vector<boost::shared_ptr<pdat::CellData<double> > > data_enthalpies;
        data_enthalpies.reserve(d_num_species);
        
        for (int si = 0; si < d_num_species; si++)
        {
            data_enthalpies.push_back(boost::make_shared<pdat::CellData<double> >(
                d_interior_box, 1, d_num_subghosts_diffusivities));
        }
        
        std::vector<double*> h;
        h.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            h.push_back(data_enthalpies[si]->getPointer(0));
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            data_enthalpies.push_back(boost::make_shared<pdat::CellData<double> >(
                d_interior_box, 1, d_num_subghosts_diffusivities));
        }        
        
        for (int si = 0; si < d_num_species; si++)
        {
            std::vector<double> species_thermo_properties;
            std::vector<double*> species_thermo_properties_ptr;
            std::vector<const double*> species_thermo_properties_const_ptr;
            
            const int num_thermo_properties = d_equation_of_state_mixing_rules->
                getNumberOfSpeciesThermodynamicProperties(si);
            
            species_thermo_properties.resize(num_thermo_properties);
            species_thermo_properties_ptr.reserve(num_thermo_properties);
            species_thermo_properties_const_ptr.reserve(num_thermo_properties);
            
            for (int ti = 0; ti < num_thermo_properties; ti++)
            {
                species_thermo_properties_ptr.push_back(&species_thermo_properties[ti]);
                species_thermo_properties_const_ptr.push_back(&species_thermo_properties[ti]);
            }
            
            d_equation_of_state_mixing_rules->getSpeciesThermodynamicProperties(
                species_thermo_properties_ptr,
                si);
            
            boost::shared_ptr<pdat::CellData<double> > data_density_species(
                new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_diffusivities));
            
            d_equation_of_state_mixing_rules->getEquationOfState(si)->
                computeDensity(
                    data_density_species,
                    d_data_pressure,
                    d_data_temperature,
                    species_thermo_properties_const_ptr,
                    d_subghost_box_diffusivities);
            
            d_equation_of_state_mixing_rules->getEquationOfState(si)->
                computeEnthalpy(
                    data_enthalpies[si],
                    data_density_species,
                    d_data_pressure,
                    species_thermo_properties_const_ptr,
                    d_subghost_box_diffusivities);
        }
        
        if (d_dim == tbox::Dimension(1))
        {
            d_data_diffusivities.reset(new pdat::CellData<double>(
                d_interior_box,
                2*d_num_species*(d_num_species + 1) + 3,
                d_num_subghosts_diffusivities));
            
            // Get the pointer to cell data of velocity and diffusivities.
            double* u = d_data_velocity->getPointer(0);
            
            std::vector<double*> D_ptr;
            D_ptr.reserve(2*d_num_species*(d_num_species + 1) + 3);
            
            for (int i = 0; i < 2*d_num_species*(d_num_species + 1) + 3; i++)
            {
                D_ptr.push_back(d_data_diffusivities->getPointer(i));
            }
            
            /*
             * Compute the diffusivities.
             */
            
            // Mass equations.
            for (int si = 0; si < d_num_species; si++)
            {
                const int component_idx = si*(d_num_species + 1);
                
                for (int i = -d_num_subghosts_diffusivities[0];
                     i < d_interior_dims[0] + d_num_subghosts_diffusivities[0];
                     i++)
                {
                    // Compute the linear indices.
                    const int idx_diffusivities = i + d_num_subghosts_diffusivities[0];
                    const int idx_density = i + d_num_subghosts_density[0];
                    
                    D_ptr[component_idx][idx_diffusivities] =
                        -rho[idx_density]*D[si][idx_diffusivities];
                }
            }
            
            for (int si = 0; si < d_num_species; si++)
            {
                for (int sj = 0; sj < d_num_species; sj++)
                {
                    const int component_idx = si*(d_num_species + 1) + sj + 1;
                    
                    for (int i = -d_num_subghosts_diffusivities[0];
                         i < d_interior_dims[0] + d_num_subghosts_diffusivities[0];
                         i++)
                    {
                        // Compute the linear indices.
                        const int idx_diffusivities = i + d_num_subghosts_diffusivities[0];
                        const int idx_density = i + d_num_subghosts_density[0];
                        const int idx_mass_fractions = i + d_num_subghosts_mass_fractions[0];
                        
                        D_ptr[component_idx][idx_diffusivities] =
                            rho[idx_density]*Y[si][idx_mass_fractions]*D[sj][idx_diffusivities];
                    }
                }
            }
            
            // Momentum and energy equations.
            for (int i = -d_num_subghosts_diffusivities[0];
                 i < d_interior_dims[0] + d_num_subghosts_diffusivities[0];
                 i++)
            {
                // Compute the linear indices.
                const int idx_diffusivities = i + d_num_subghosts_diffusivities[0];
                const int idx_velocity = i + d_num_subghosts_velocity[0];
                
                D_ptr[d_num_species*(d_num_species + 1)][idx_diffusivities] =
                    -(double(4)/double(3)*mu[idx_diffusivities] + mu_v[idx_diffusivities]);
                D_ptr[d_num_species*(d_num_species + 1) + 1][idx_diffusivities] =
                    -u[idx_velocity]*(double(4)/double(3)*mu[idx_diffusivities] + mu_v[idx_diffusivities]);
                D_ptr[d_num_species*(d_num_species + 1) + 2][idx_diffusivities] =
                    -kappa[idx_diffusivities];
            }
            
            // Energy equation.
            for (int si = 0; si < d_num_species; si++)
            {
                const int component_idx =
                    d_num_species*(d_num_species + 1) + 3 + si*(d_num_species + 1);
                
                for (int i = -d_num_subghosts_diffusivities[0];
                     i < d_interior_dims[0] + d_num_subghosts_diffusivities[0];
                     i++)
                {
                    // Compute the linear indices.
                    const int idx_diffusivities = i + d_num_subghosts_diffusivities[0];
                    const int idx_density = i + d_num_subghosts_density[0];
                    
                    D_ptr[component_idx][idx_diffusivities] =
                        -rho[idx_density]*D[si][idx_diffusivities]*
                            h[si][idx_diffusivities];
                }
            }
            
            for (int si = 0; si < d_num_species; si++)
            {
                for (int sj = 0; sj < d_num_species; sj++)
                {
                    const int component_idx =
                        d_num_species*(d_num_species + 1) + 4 + si*(d_num_species + 1) + sj;
                    
                    for (int i = -d_num_subghosts_diffusivities[0];
                         i < d_interior_dims[0] + d_num_subghosts_diffusivities[0];
                         i++)
                    {
                        // Compute the linear indices.
                        const int idx_diffusivities = i + d_num_subghosts_diffusivities[0];
                        const int idx_density = i + d_num_subghosts_density[0];
                        const int idx_mass_fractions = i + d_num_subghosts_mass_fractions[0];
                        
                        D_ptr[component_idx][idx_diffusivities] =
                            rho[idx_density]*Y[si][idx_mass_fractions]*D[sj][idx_diffusivities]*
                                h[si][idx_diffusivities];
                    }
                }
            }
        }
        else if (d_dim == tbox::Dimension(2))
        {
            d_data_diffusivities.reset(new pdat::CellData<double>(
                d_interior_box,
                2*d_num_species*(d_num_species + 1) + 10,
                d_num_subghosts_diffusivities));
            
            // Get the pointer to cell data of velocity and diffusivities.
            double* u = d_data_velocity->getPointer(0);
            double* v = d_data_velocity->getPointer(1);
            
            std::vector<double*> D_ptr;
            D_ptr.reserve(2*d_num_species*(d_num_species + 1) + 10);
            
            for (int i = 0; i < 2*d_num_species*(d_num_species + 1) + 10; i++)
            {
                D_ptr.push_back(d_data_diffusivities->getPointer(i));
            }
            
            /*
             * Compute the diffusivities.
             */
            
            // Mass equations.
            for (int si = 0; si < d_num_species; si++)
            {
                const int component_idx = si*(d_num_species + 1);
                
                for (int j = -d_num_subghosts_diffusivities[1];
                     j < d_interior_dims[1] + d_num_subghosts_diffusivities[1];
                     j++)
                {
                    for (int i = -d_num_subghosts_diffusivities[0];
                         i < d_interior_dims[0] + d_num_subghosts_diffusivities[0];
                         i++)
                    {
                        // Compute the linear indices.
                        const int idx_diffusivities = (i + d_num_subghosts_diffusivities[0]) +
                            (j + d_num_subghosts_diffusivities[1])*d_subghostcell_dims_diffusivities[0];
                        
                        const int idx_density = (i + d_num_subghosts_density[0]) +
                            (j + d_num_subghosts_density[1])*d_subghostcell_dims_density[0];
                        
                        D_ptr[component_idx][idx_diffusivities] =
                            -rho[idx_density]*D[si][idx_diffusivities];
                    }
                }
            }
            
            for (int si = 0; si < d_num_species; si++)
            {
                for (int sj = 0; sj < d_num_species; sj++)
                {
                    const int component_idx = si*(d_num_species + 1) + sj + 1;
                    
                    for (int j = -d_num_subghosts_diffusivities[1];
                         j < d_interior_dims[1] + d_num_subghosts_diffusivities[1];
                         j++)
                    {
                        for (int i = -d_num_subghosts_diffusivities[0];
                             i < d_interior_dims[0] + d_num_subghosts_diffusivities[0];
                             i++)
                        {
                            // Compute the linear indices.
                            const int idx_diffusivities = (i + d_num_subghosts_diffusivities[0]) +
                                (j + d_num_subghosts_diffusivities[1])*d_subghostcell_dims_diffusivities[0];
                            
                            const int idx_density = (i + d_num_subghosts_density[0]) +
                                (j + d_num_subghosts_density[1])*d_subghostcell_dims_density[0];
                            
                            const int idx_mass_fractions = (i + d_num_subghosts_mass_fractions[0]) +
                                (j + d_num_subghosts_mass_fractions[1])*d_subghostcell_dims_mass_fractions[0];
                            
                            D_ptr[component_idx][idx_diffusivities] =
                                rho[idx_density]*Y[si][idx_mass_fractions]*D[sj][idx_diffusivities];
                        }
                    }
                }
            }
            
            // Momentum and energy equations.
            for (int j = -d_num_subghosts_diffusivities[1];
                 j < d_interior_dims[1] + d_num_subghosts_diffusivities[1];
                 j++)
            {
                for (int i = -d_num_subghosts_diffusivities[0];
                     i < d_interior_dims[0] + d_num_subghosts_diffusivities[0];
                     i++)
                {
                    // Compute the linear indices.
                    const int idx_diffusivities = (i + d_num_subghosts_diffusivities[0]) +
                        (j + d_num_subghosts_diffusivities[1])*d_subghostcell_dims_diffusivities[0];
                    
                    const int idx_velocity = (i + d_num_subghosts_velocity[0]) +
                        (j + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0];
                    
                    D_ptr[d_num_species*(d_num_species + 1)][idx_diffusivities] =
                        -(double(4)/double(3)*mu[idx_diffusivities] + mu_v[idx_diffusivities]);
                    D_ptr[d_num_species*(d_num_species + 1) + 1][idx_diffusivities] =
                        double(2)/double(3)*mu[idx_diffusivities] - mu_v[idx_diffusivities];
                    D_ptr[d_num_species*(d_num_species + 1) + 2][idx_diffusivities] =
                        -mu[idx_diffusivities];
                    D_ptr[d_num_species*(d_num_species + 1) + 3][idx_diffusivities] =
                        -u[idx_velocity]*(double(4)/double(3)*mu[idx_diffusivities] + mu_v[idx_diffusivities]);
                    D_ptr[d_num_species*(d_num_species + 1) + 4][idx_diffusivities] =
                        -v[idx_velocity]*(double(4)/double(3)*mu[idx_diffusivities] + mu_v[idx_diffusivities]);
                    D_ptr[d_num_species*(d_num_species + 1) + 5][idx_diffusivities] =
                        u[idx_velocity]*(double(2)/double(3)*mu[idx_diffusivities] - mu_v[idx_diffusivities]);
                    D_ptr[d_num_species*(d_num_species + 1) + 6][idx_diffusivities] =
                        v[idx_velocity]*(double(2)/double(3)*mu[idx_diffusivities] - mu_v[idx_diffusivities]);
                    D_ptr[d_num_species*(d_num_species + 1) + 7][idx_diffusivities] =
                        -u[idx_velocity]*mu[idx_diffusivities];
                    D_ptr[d_num_species*(d_num_species + 1) + 8][idx_diffusivities] =
                        -v[idx_velocity]*mu[idx_diffusivities];
                    D_ptr[d_num_species*(d_num_species + 1) + 9][idx_diffusivities] =
                        -kappa[idx_diffusivities];
                }
            }
            
            // Energy equation.
            for (int si = 0; si < d_num_species; si++)
            {
                const int component_idx =
                    d_num_species*(d_num_species + 1) + 10 + si*(d_num_species + 1);
                
                for (int j = -d_num_subghosts_diffusivities[1];
                     j < d_interior_dims[1] + d_num_subghosts_diffusivities[1];
                     j++)
                {
                    for (int i = -d_num_subghosts_diffusivities[0];
                         i < d_interior_dims[0] + d_num_subghosts_diffusivities[0];
                         i++)
                    {
                        // Compute the linear indices.
                        const int idx_diffusivities = (i + d_num_subghosts_diffusivities[0]) +
                            (j + d_num_subghosts_diffusivities[1])*d_subghostcell_dims_diffusivities[0];
                        
                        const int idx_density = (i + d_num_subghosts_density[0]) +
                            (j + d_num_subghosts_density[1])*d_subghostcell_dims_density[0];
                        
                        D_ptr[component_idx][idx_diffusivities] =
                            -rho[idx_density]*D[si][idx_diffusivities]*
                                h[si][idx_diffusivities];
                    }
                }
            }
            
            for (int si = 0; si < d_num_species; si++)
            {
                for (int sj = 0; sj < d_num_species; sj++)
                {
                    const int component_idx =
                        d_num_species*(d_num_species + 1) + 11 + si*(d_num_species + 1) + sj;
                    
                    for (int j = -d_num_subghosts_diffusivities[1];
                         j < d_interior_dims[1] + d_num_subghosts_diffusivities[1];
                         j++)
                    {
                        for (int i = -d_num_subghosts_diffusivities[0];
                             i < d_interior_dims[0] + d_num_subghosts_diffusivities[0];
                             i++)
                        {
                            // Compute the linear indices.
                            const int idx_diffusivities = (i + d_num_subghosts_diffusivities[0]) +
                                (j + d_num_subghosts_diffusivities[1])*d_subghostcell_dims_diffusivities[0];
                            
                            const int idx_density = (i + d_num_subghosts_density[0]) +
                                (j + d_num_subghosts_density[1])*d_subghostcell_dims_density[0];
                            
                            const int idx_mass_fractions = (i + d_num_subghosts_mass_fractions[0]) +
                                (j + d_num_subghosts_mass_fractions[1])*d_subghostcell_dims_mass_fractions[0];
                            
                            D_ptr[component_idx][idx_diffusivities] =
                                rho[idx_density]*Y[si][idx_mass_fractions]*D[sj][idx_diffusivities]*
                                    h[si][idx_diffusivities];
                        }
                    }
                }
            }
        }
        else if (d_dim == tbox::Dimension(3))
        {
            d_data_diffusivities.reset(new pdat::CellData<double>(
                d_interior_box,
                2*d_num_species*(d_num_species + 1) + 13,
                d_num_subghosts_diffusivities));
            
            // Get the pointer to cell data of velocity and diffusivities.
            double* u = d_data_velocity->getPointer(0);
            double* v = d_data_velocity->getPointer(1);
            double* w = d_data_velocity->getPointer(2);
            
            std::vector<double*> D_ptr;
            D_ptr.reserve(2*d_num_species*(d_num_species + 1) + 13);
            
            for (int i = 0; i < 2*d_num_species*(d_num_species + 1) + 13; i++)
            {
                D_ptr.push_back(d_data_diffusivities->getPointer(i));
            }
            
            /*
             * Compute the diffusivities.
             */
            
            //Mass equations.
            for (int si = 0; si < d_num_species; si++)
            {
                const int component_idx = si*(d_num_species + 1);
                
                for (int k = -d_num_subghosts_diffusivities[2];
                     k < d_interior_dims[2] + d_num_subghosts_diffusivities[2];
                     k++)
                {
                    for (int j = -d_num_subghosts_diffusivities[1];
                         j < d_interior_dims[1] + d_num_subghosts_diffusivities[1];
                         j++)
                    {
                        for (int i = -d_num_subghosts_diffusivities[0];
                             i < d_interior_dims[0] + d_num_subghosts_diffusivities[0];
                             i++)
                        {
                            // Compute the linear indices.
                            const int idx_diffusivities = (i + d_num_subghosts_diffusivities[0]) +
                                (j + d_num_subghosts_diffusivities[1])*d_subghostcell_dims_diffusivities[0] +
                                (k + d_num_subghosts_diffusivities[2])*d_subghostcell_dims_diffusivities[0]*
                                    d_subghostcell_dims_diffusivities[1];
                            
                            const int idx_density = (i + d_num_subghosts_density[0]) +
                                (j + d_num_subghosts_density[1])*d_subghostcell_dims_density[0] +
                                (k + d_num_subghosts_density[2])*d_subghostcell_dims_density[0]*
                                    d_subghostcell_dims_density[1];
                            
                            D_ptr[component_idx][idx_diffusivities] =
                                -rho[idx_density]*D[si][idx_diffusivities];
                        }
                    }
                }
            }
            
            for (int si = 0; si < d_num_species; si++)
            {
                for (int sj = 0; sj < d_num_species; sj++)
                {
                    const int component_idx = si*(d_num_species + 1) + sj + 1;
                    
                    for (int k = -d_num_subghosts_diffusivities[2];
                         k < d_interior_dims[2] + d_num_subghosts_diffusivities[2];
                         k++)
                    {
                        for (int j = -d_num_subghosts_diffusivities[1];
                             j < d_interior_dims[1] + d_num_subghosts_diffusivities[1];
                             j++)
                        {
                            for (int i = -d_num_subghosts_diffusivities[0];
                                 i < d_interior_dims[0] + d_num_subghosts_diffusivities[0];
                                 i++)
                            {
                                // Compute the linear indices.
                                const int idx_diffusivities = (i + d_num_subghosts_diffusivities[0]) +
                                    (j + d_num_subghosts_diffusivities[1])*d_subghostcell_dims_diffusivities[0] +
                                    (k + d_num_subghosts_diffusivities[2])*d_subghostcell_dims_diffusivities[0]*
                                        d_subghostcell_dims_diffusivities[1];
                                
                                const int idx_density = (i + d_num_subghosts_density[0]) +
                                    (j + d_num_subghosts_density[1])*d_subghostcell_dims_density[0] +
                                    (k + d_num_subghosts_density[2])*d_subghostcell_dims_density[0]*
                                        d_subghostcell_dims_density[1];
                                
                                const int idx_mass_fractions = (i + d_num_subghosts_mass_fractions[0]) +
                                    (j + d_num_subghosts_mass_fractions[1])*d_subghostcell_dims_mass_fractions[0] +
                                    (k + d_num_subghosts_mass_fractions[2])*d_subghostcell_dims_mass_fractions[0]*
                                        d_subghostcell_dims_mass_fractions[1];
                                
                                D_ptr[component_idx][idx_diffusivities] =
                                    rho[idx_density]*Y[si][idx_mass_fractions]*D[sj][idx_diffusivities];
                            }
                        }
                    }
                }
            }
            
            // Momentum and energy equations.
            for (int k = -d_num_subghosts_diffusivities[2];
                 k < d_interior_dims[2] + d_num_subghosts_diffusivities[2];
                 k++)
            {
                for (int j = -d_num_subghosts_diffusivities[1];
                     j < d_interior_dims[1] + d_num_subghosts_diffusivities[1];
                     j++)
                {
                    for (int i = -d_num_subghosts_diffusivities[0];
                         i < d_interior_dims[0] + d_num_subghosts_diffusivities[0];
                         i++)
                    {
                        // Compute the linear indices.
                        const int idx_diffusivities = (i + d_num_subghosts_diffusivities[0]) +
                            (j + d_num_subghosts_diffusivities[1])*d_subghostcell_dims_diffusivities[0] +
                            (k + d_num_subghosts_diffusivities[2])*d_subghostcell_dims_diffusivities[0]*
                                d_subghostcell_dims_diffusivities[1];
                        
                        const int idx_velocity = (i + d_num_subghosts_velocity[0]) +
                            (j + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0] +
                            (k + d_num_subghosts_velocity[2])*d_subghostcell_dims_velocity[0]*
                                d_subghostcell_dims_velocity[1];
                        
                        D_ptr[d_num_species*(d_num_species + 1)][idx_diffusivities] =
                            -(double(4)/double(3)*mu[idx_diffusivities] + mu_v[idx_diffusivities]);
                        D_ptr[d_num_species*(d_num_species + 1) + 1][idx_diffusivities] =
                            double(2)/double(3)*mu[idx_diffusivities] - mu_v[idx_diffusivities];
                        D_ptr[d_num_species*(d_num_species + 1) + 2][idx_diffusivities] =
                            -mu[idx_diffusivities];
                        D_ptr[d_num_species*(d_num_species + 1) + 3][idx_diffusivities] =
                            -u[idx_velocity]*(double(4)/double(3)*mu[idx_diffusivities] + mu_v[idx_diffusivities]);
                        D_ptr[d_num_species*(d_num_species + 1) + 4][idx_diffusivities] =
                            -v[idx_velocity]*(double(4)/double(3)*mu[idx_diffusivities] + mu_v[idx_diffusivities]);
                        D_ptr[d_num_species*(d_num_species + 1) + 5][idx_diffusivities] =
                            -w[idx_velocity]*(double(4)/double(3)*mu[idx_diffusivities] + mu_v[idx_diffusivities]);
                        D_ptr[d_num_species*(d_num_species + 1) + 6][idx_diffusivities] =
                            u[idx_velocity]*(double(2)/double(3)*mu[idx_diffusivities] - mu_v[idx_diffusivities]);
                        D_ptr[d_num_species*(d_num_species + 1) + 7][idx_diffusivities] =
                            v[idx_velocity]*(double(2)/double(3)*mu[idx_diffusivities] - mu_v[idx_diffusivities]);
                        D_ptr[d_num_species*(d_num_species + 1) + 8][idx_diffusivities] =
                            w[idx_velocity]*(double(2)/double(3)*mu[idx_diffusivities] - mu_v[idx_diffusivities]);
                        D_ptr[d_num_species*(d_num_species + 1) + 9][idx_diffusivities] =
                            -u[idx_velocity]*mu[idx_diffusivities];
                        D_ptr[d_num_species*(d_num_species + 1) + 10][idx_diffusivities] =
                            -v[idx_velocity]*mu[idx_diffusivities];
                        D_ptr[d_num_species*(d_num_species + 1) + 11][idx_diffusivities] =
                            -w[idx_velocity]*mu[idx_diffusivities];
                        D_ptr[d_num_species*(d_num_species + 1) + 12][idx_diffusivities] =
                            -kappa[idx_diffusivities];
                    }
                }
            }
            
            // Energy equation.
            for (int si = 0; si < d_num_species; si++)
            {
                const int component_idx =
                    d_num_species*(d_num_species + 1) + 13 + si*(d_num_species + 1);
                
                for (int k = -d_num_subghosts_diffusivities[2];
                     k < d_interior_dims[2] + d_num_subghosts_diffusivities[2];
                     k++)
                {
                    for (int j = -d_num_subghosts_diffusivities[1];
                         j < d_interior_dims[1] + d_num_subghosts_diffusivities[1];
                         j++)
                    {
                        for (int i = -d_num_subghosts_diffusivities[0];
                             i < d_interior_dims[0] + d_num_subghosts_diffusivities[0];
                             i++)
                        {
                            // Compute the linear indices.
                            const int idx_diffusivities = (i + d_num_subghosts_diffusivities[0]) +
                                (j + d_num_subghosts_diffusivities[1])*d_subghostcell_dims_diffusivities[0] +
                                (k + d_num_subghosts_diffusivities[2])*d_subghostcell_dims_diffusivities[0]*
                                    d_subghostcell_dims_diffusivities[1];
                            
                            const int idx_density = (i + d_num_subghosts_density[0]) +
                                (j + d_num_subghosts_density[1])*d_subghostcell_dims_density[0] +
                                (k + d_num_subghosts_density[2])*d_subghostcell_dims_density[0]*
                                    d_subghostcell_dims_density[1];
                            
                            D_ptr[component_idx][idx_diffusivities] =
                                -rho[idx_density]*D[si][idx_diffusivities]*
                                    h[si][idx_diffusivities];
                        }
                    }
                }
            }
            
            for (int si = 0; si < d_num_species; si++)
            {
                for (int sj = 0; sj < d_num_species; sj++)
                {
                    const int component_idx =
                        d_num_species*(d_num_species + 1) + 14 + si*(d_num_species + 1) + sj;
                    
                    for (int k = -d_num_subghosts_diffusivities[2];
                         k < d_interior_dims[2] + d_num_subghosts_diffusivities[2];
                         k++)
                    {
                        for (int j = -d_num_subghosts_diffusivities[1];
                             j < d_interior_dims[1] + d_num_subghosts_diffusivities[1];
                             j++)
                        {
                            for (int i = -d_num_subghosts_diffusivities[0];
                                 i < d_interior_dims[0] + d_num_subghosts_diffusivities[0];
                                 i++)
                            {
                                // Compute the linear indices.
                                const int idx_diffusivities = (i + d_num_subghosts_diffusivities[0]) +
                                    (j + d_num_subghosts_diffusivities[1])*d_subghostcell_dims_diffusivities[0] +
                                    (k + d_num_subghosts_diffusivities[2])*d_subghostcell_dims_diffusivities[0]*
                                        d_subghostcell_dims_diffusivities[1];
                                
                                const int idx_density = (i + d_num_subghosts_density[0]) +
                                    (j + d_num_subghosts_density[1])*d_subghostcell_dims_density[0] +
                                    (k + d_num_subghosts_density[2])*d_subghostcell_dims_density[0]*
                                        d_subghostcell_dims_density[1];
                                
                                const int idx_mass_fractions = (i + d_num_subghosts_mass_fractions[0]) +
                                    (j + d_num_subghosts_mass_fractions[1])*d_subghostcell_dims_mass_fractions[0] +
                                    (k + d_num_subghosts_mass_fractions[2])*d_subghostcell_dims_mass_fractions[0]*
                                        d_subghostcell_dims_mass_fractions[1];
                                
                                D_ptr[component_idx][idx_diffusivities] =
                                    rho[idx_density]*Y[si][idx_mass_fractions]*D[sj][idx_diffusivities]*
                                        h[si][idx_diffusivities];
                            }
                        }
                    }
                }
            }
        }
        
        data_mass_diffusivities.reset();
        data_shear_viscosity.reset();
        data_bulk_viscosity.reset();
        data_thermal_conductivity.reset();
        data_enthalpies.clear();
    }
    
    if (d_dim == tbox::Dimension(1))
    {
        switch (flux_direction)
        {
            case DIRECTION::X_DIRECTION:
            {
                switch (derivative_direction)
                {
                    case DIRECTION::X_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            diffusivities_data[si].resize(d_num_species + 1);
                            diffusivities_component_idx[si].resize(d_num_species + 1);
                            
                            diffusivities_data[si][0] = d_data_diffusivities;
                            diffusivities_component_idx[si][0] = si*(d_num_species + 1);
                            
                            for (int sj = 0; sj < d_num_species; sj++)
                            {
                                diffusivities_data[si][1 + sj] = d_data_diffusivities;
                                diffusivities_component_idx[si][1 + sj] = si*(d_num_species + 1) + sj + 1;
                            }
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[d_num_species].resize(1);
                        diffusivities_component_idx[d_num_species].resize(1);
                        
                        // -(4/3*mu + mu_v).
                        diffusivities_data[d_num_species][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species][0] =
                            d_num_species*(d_num_species + 1);
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[d_num_species + 1].resize(
                            2 + d_num_species*(d_num_species + 1));
                        diffusivities_component_idx[d_num_species + 1].resize(
                            2 + d_num_species*(d_num_species + 1));
                        
                        // -u*(4/3*mu + mu_v).
                        diffusivities_data[d_num_species + 1][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 1][0] =
                            d_num_species*(d_num_species + 1) + 1;
                        
                        // -kappa.
                        diffusivities_data[d_num_species + 1][1] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 1][1] =
                            d_num_species*(d_num_species + 1) + 2;
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            diffusivities_data[d_num_species + 1][2 + si*(d_num_species + 1)] =
                                d_data_diffusivities;
                            diffusivities_component_idx[d_num_species + 1][2 + si*(d_num_species + 1)] =
                                d_num_species*(d_num_species + 1) + 3 + si*(d_num_species + 1);
                            
                            for (int sj = 0; sj < d_num_species; sj++)
                            {
                                diffusivities_data[d_num_species + 1][3 + si*(d_num_species + 1) + sj] =
                                    d_data_diffusivities;
                                diffusivities_component_idx[d_num_species + 1][3 + si*(d_num_species + 1) + sj] =
                                    d_num_species*(d_num_species + 1) + 4 + si*(d_num_species + 1) + sj;
                            }
                        }
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelFourEqnConservative::getDiffusiveFluxDiffusivities()\n"
                            << "There are only x-direction for one-dimensional problem."
                            << std::endl);
                    }
                }
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelFourEqnConservative::getDiffusiveFluxDiffusivities()\n"
                    << "There are only x-direction for one-dimensional problem."
                    << std::endl);
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        switch (flux_direction)
        {
            case DIRECTION::X_DIRECTION:
            {
                switch (derivative_direction)
                {
                    case DIRECTION::X_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            diffusivities_data[si].resize(d_num_species + 1);
                            diffusivities_component_idx[si].resize(d_num_species + 1);
                            
                            diffusivities_data[si][0] = d_data_diffusivities;
                            diffusivities_component_idx[si][0] = si*(d_num_species + 1);
                            
                            for (int sj = 0; sj < d_num_species; sj++)
                            {
                                diffusivities_data[si][1 + sj] = d_data_diffusivities;
                                diffusivities_component_idx[si][1 + sj] = si*(d_num_species + 1) + sj + 1;
                            }
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[d_num_species].resize(1);
                        diffusivities_component_idx[d_num_species].resize(1);
                        
                        // -(4/3*mu + mu_v).
                        diffusivities_data[d_num_species][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species][0] =
                            d_num_species*(d_num_species + 1);
                        
                        diffusivities_data[d_num_species + 1].resize(1);
                        diffusivities_component_idx[d_num_species + 1].resize(1);
                        
                        // -mu.
                        diffusivities_data[d_num_species + 1][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 1][0] =
                            d_num_species*(d_num_species + 1) + 2;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[d_num_species + 2].resize(
                            3 + d_num_species*(d_num_species + 1));
                        diffusivities_component_idx[d_num_species + 2].resize(
                            3 + d_num_species*(d_num_species + 1));
                        
                        // -u*(4/3*mu + mu_v).
                        diffusivities_data[d_num_species + 2][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 2][0] =
                            d_num_species*(d_num_species + 1) + 3;
                        
                        // -v*mu.
                        diffusivities_data[d_num_species + 2][1] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 2][1] =
                            d_num_species*(d_num_species + 1) + 8;
                        
                        // -kappa.
                        diffusivities_data[d_num_species + 2][2] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 2][2] =
                            d_num_species*(d_num_species + 1) + 9;
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            diffusivities_data[d_num_species + 2][3 + si*(d_num_species + 1)] =
                                d_data_diffusivities;
                            diffusivities_component_idx[d_num_species + 2][3 + si*(d_num_species + 1)] =
                                d_num_species*(d_num_species + 1) + 10 + si*(d_num_species + 1);
                            
                            for (int sj = 0; sj < d_num_species; sj++)
                            {
                                diffusivities_data[d_num_species + 2][4 + si*(d_num_species + 1) + sj] =
                                    d_data_diffusivities;
                                diffusivities_component_idx[d_num_species + 2][4 + si*(d_num_species + 1) + sj] =
                                    d_num_species*(d_num_species + 1) + 11 + si*(d_num_species + 1) + sj;
                            }
                        }
                        
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            diffusivities_data[si].resize(0);
                            diffusivities_component_idx[si].resize(0);
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[d_num_species].resize(1);
                        diffusivities_component_idx[d_num_species].resize(1);
                        
                        // 2/3*mu - mu_v.
                        diffusivities_data[d_num_species][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species][0] =
                            d_num_species*(d_num_species + 1) + 1;
                        
                        diffusivities_data[d_num_species + 1].resize(1);
                        diffusivities_component_idx[d_num_species + 1].resize(1);
                        
                        // -mu.
                        diffusivities_data[d_num_species + 1][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 1][0] =
                            d_num_species*(d_num_species + 1) + 2;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[d_num_species + 2].resize(2);
                        diffusivities_component_idx[d_num_species + 2].resize(2);
                        
                        // -v*mu.
                        diffusivities_data[d_num_species + 2][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 2][0] =
                            d_num_species*(d_num_species + 1) + 8;
                        
                        // u*(2/3*mu - mu_v).
                        diffusivities_data[d_num_species + 2][1] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 2][1] =
                            d_num_species*(d_num_species + 1) + 5;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelFourEqnConservative::getDiffusiveFluxDiffusivities()\n"
                            << "There are only x-direction and y-direction for two-dimensional problem."
                            << std::endl);
                    }
                }
                
                break;
            }
            case DIRECTION::Y_DIRECTION:
            {
                switch (derivative_direction)
                {
                    case DIRECTION::X_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            diffusivities_data[si].resize(0);
                            diffusivities_component_idx[si].resize(0);
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[d_num_species].resize(1);
                        diffusivities_component_idx[d_num_species].resize(1);
                        
                        // -mu.
                        diffusivities_data[d_num_species][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species][0] =
                            d_num_species*(d_num_species + 1) + 2;
                        
                        diffusivities_data[d_num_species + 1].resize(1);
                        diffusivities_component_idx[d_num_species + 1].resize(1);
                        
                        // 2/3*mu - mu_v.
                        diffusivities_data[d_num_species + 1][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 1][0] =
                            d_num_species*(d_num_species + 1) + 1;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[d_num_species + 2].resize(2);
                        diffusivities_component_idx[d_num_species + 2].resize(2);
                        
                        // v*(2/3*mu - mu_v).
                        diffusivities_data[d_num_species + 2][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 2][0] =
                            d_num_species*(d_num_species + 1) + 6;
                        
                        // -u*mu.
                        diffusivities_data[d_num_species + 2][1] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 2][1] =
                            d_num_species*(d_num_species + 1) + 7;
                        
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            diffusivities_data[si].resize(d_num_species + 1);
                            diffusivities_component_idx[si].resize(d_num_species + 1);
                            
                            diffusivities_data[si][0] = d_data_diffusivities;
                            diffusivities_component_idx[si][0] = si*(d_num_species + 1);
                            
                            for (int sj = 0; sj < d_num_species; sj++)
                            {
                                diffusivities_data[si][1 + sj] = d_data_diffusivities;
                                diffusivities_component_idx[si][1 + sj] = si*(d_num_species + 1) + sj + 1;
                            }
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[d_num_species].resize(1);
                        diffusivities_component_idx[d_num_species].resize(1);
                        
                        // -mu.
                        diffusivities_data[d_num_species][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species][0] =
                            d_num_species*(d_num_species + 1) + 2;
                        
                        diffusivities_data[d_num_species + 1].resize(1);
                        diffusivities_component_idx[d_num_species + 1].resize(1);
                        
                        // -(4/3*mu + mu_v).
                        diffusivities_data[d_num_species + 1][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 1][0] =
                            d_num_species*(d_num_species + 1);
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[d_num_species + 2].resize(
                            3 + d_num_species*(d_num_species + 1));
                        diffusivities_component_idx[d_num_species + 2].resize(
                            3 + d_num_species*(d_num_species + 1));
                        
                        // -u*mu.
                        diffusivities_data[d_num_species + 2][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 2][0] =
                            d_num_species*(d_num_species + 1) + 7;
                        
                        // -v*(4/3*mu + mu_v).
                        diffusivities_data[d_num_species + 2][1] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 2][1] =
                            d_num_species*(d_num_species + 1) + 4;
                        
                        // -kappa.
                        diffusivities_data[d_num_species + 2][2] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 2][2] =
                            d_num_species*(d_num_species + 1) + 9;
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            diffusivities_data[d_num_species + 2][3 + si*(d_num_species + 1)] =
                                d_data_diffusivities;
                            diffusivities_component_idx[d_num_species + 2][3 + si*(d_num_species + 1)] =
                                d_num_species*(d_num_species + 1) + 10 + si*(d_num_species + 1);
                            
                            for (int sj = 0; sj < d_num_species; sj++)
                            {
                                diffusivities_data[d_num_species + 2][4 + si*(d_num_species + 1) + sj] =
                                    d_data_diffusivities;
                                diffusivities_component_idx[d_num_species + 2][4 + si*(d_num_species + 1) + sj] =
                                    d_num_species*(d_num_species + 1) + 11 + si*(d_num_species + 1) + sj;
                            }
                        }
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelFourEqnConservative::getDiffusiveFluxDiffusivities()\n"
                            << "There are only x-direction and y-direction for two-dimensional problem."
                            << std::endl);
                    }
                }
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelFourEqnConservative::getDiffusiveFluxDiffusivities()\n"
                    << "There are only x-direction and y-direction for two-dimensional problem."
                    << std::endl);
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        switch (flux_direction)
        {
            case DIRECTION::X_DIRECTION:
            {
                switch (derivative_direction)
                {
                    case DIRECTION::X_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            diffusivities_data[si].resize(d_num_species + 1);
                            diffusivities_component_idx[si].resize(d_num_species + 1);
                            
                            diffusivities_data[si][0] = d_data_diffusivities;
                            diffusivities_component_idx[si][0] = si*(d_num_species + 1);
                            
                            for (int sj = 0; sj < d_num_species; sj++)
                            {
                                diffusivities_data[si][1 + sj] = d_data_diffusivities;
                                diffusivities_component_idx[si][1 + sj] = si*(d_num_species + 1) + sj + 1;
                            }
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[d_num_species].resize(1);
                        diffusivities_component_idx[d_num_species].resize(1);
                        
                        // -(4/3*mu + mu_v).
                        diffusivities_data[d_num_species][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species][0] =
                            d_num_species*(d_num_species + 1);
                        
                        diffusivities_data[d_num_species + 1].resize(1);
                        diffusivities_component_idx[d_num_species + 1].resize(1);
                        
                        // -mu.
                        diffusivities_data[d_num_species + 1][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 1][0] =
                            d_num_species*(d_num_species + 1) + 2;
                        
                        diffusivities_data[d_num_species + 2].resize(1);
                        diffusivities_component_idx[d_num_species + 2].resize(1);
                        
                        // -mu.
                        diffusivities_data[d_num_species + 2][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 2][0] =
                            d_num_species*(d_num_species + 1) + 2;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[d_num_species + 3].resize(
                            4 + d_num_species*(d_num_species + 1));
                        diffusivities_component_idx[d_num_species + 3].resize(
                            4 + d_num_species*(d_num_species + 1));
                        
                        // -u*(4/3*mu + mu_v).
                        diffusivities_data[d_num_species + 3][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][0] =
                            d_num_species*(d_num_species + 1) + 3;
                        
                        // -v*mu.
                        diffusivities_data[d_num_species + 3][1] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][1] =
                            d_num_species*(d_num_species + 1) + 10;
                        
                        // -w*mu.
                        diffusivities_data[d_num_species + 3][2] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][2] =
                            d_num_species*(d_num_species + 1) + 11;
                        
                        // -kappa.
                        diffusivities_data[d_num_species + 3][3] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][3] =
                            d_num_species*(d_num_species + 1) + 12;
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            diffusivities_data[d_num_species + 3][4 + si*(d_num_species + 1)] =
                                d_data_diffusivities;
                            diffusivities_component_idx[d_num_species + 3][4 + si*(d_num_species + 1)] =
                                d_num_species*(d_num_species + 1) + 13 + si*(d_num_species + 1);
                            
                            for (int sj = 0; sj < d_num_species; sj++)
                            {
                                diffusivities_data[d_num_species + 3][5 + si*(d_num_species + 1) + sj] =
                                    d_data_diffusivities;
                                diffusivities_component_idx[d_num_species + 3][5 + si*(d_num_species + 1) + sj] =
                                    d_num_species*(d_num_species + 1) + 14 + si*(d_num_species + 1) + sj;
                            }
                        }
                        
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            diffusivities_data[si].resize(0);
                            diffusivities_component_idx[si].resize(0);
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[d_num_species].resize(1);
                        diffusivities_component_idx[d_num_species].resize(1);
                        
                        // 2/3*mu - mu_v.
                        diffusivities_data[d_num_species][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species][0] =
                            d_num_species*(d_num_species + 1) + 1;
                        
                        diffusivities_data[d_num_species + 1].resize(1);
                        diffusivities_component_idx[d_num_species + 1].resize(1);
                        
                        // -mu.
                        diffusivities_data[d_num_species + 1][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 1][0] =
                            d_num_species*(d_num_species + 1) + 2;
                        
                        diffusivities_data[d_num_species + 2].resize(0);
                        diffusivities_component_idx[d_num_species + 2].resize(0);
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[d_num_species + 3].resize(2);
                        diffusivities_component_idx[d_num_species + 3].resize(2);
                        
                        // -v*mu.
                        diffusivities_data[d_num_species + 3][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][0] =
                            d_num_species*(d_num_species + 1) + 10;
                        
                        // u*(2/3*mu - mu_v).
                        diffusivities_data[d_num_species + 3][1] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][1] =
                            d_num_species*(d_num_species + 1) + 6;
                        
                        break;
                    }
                    case DIRECTION::Z_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            diffusivities_data[si].resize(0);
                            diffusivities_component_idx[si].resize(0);
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[d_num_species].resize(1);
                        diffusivities_component_idx[d_num_species].resize(1);
                        
                        // 2/3*mu - mu_v.
                        diffusivities_data[d_num_species][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species][0] =
                            d_num_species*(d_num_species + 1) + 1;
                        
                        diffusivities_data[d_num_species + 1].resize(0);
                        diffusivities_component_idx[d_num_species + 1].resize(0);
                        
                        diffusivities_data[d_num_species + 2].resize(1);
                        diffusivities_component_idx[d_num_species + 2].resize(1);
                        
                        // -mu.
                        diffusivities_data[d_num_species + 2][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 2][0] =
                            d_num_species*(d_num_species + 1) + 2;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[d_num_species + 3].resize(2);
                        diffusivities_component_idx[d_num_species + 3].resize(2);
                        
                        // -w*mu.
                        diffusivities_data[d_num_species + 3][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][0] =
                            d_num_species*(d_num_species + 1) + 11;
                        
                        // u*(2/3*mu - mu_v).
                        diffusivities_data[d_num_species + 3][1] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][1] =
                            d_num_species*(d_num_species + 1) + 6;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelFourEqnConservative::getDiffusiveFluxDiffusivities()\n"
                            << "There are only x-direction, y-direction and z-direction for three-dimensional problem."
                            << std::endl);
                    }
                }
                
                break;
            }
            case DIRECTION::Y_DIRECTION:
            {
                switch (derivative_direction)
                {
                    case DIRECTION::X_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            diffusivities_data[si].resize(0);
                            diffusivities_component_idx[si].resize(0);
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[d_num_species].resize(1);
                        diffusivities_component_idx[d_num_species].resize(1);
                        
                        // -mu.
                        diffusivities_data[d_num_species][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species][0] =
                            d_num_species*(d_num_species + 1) + 2;
                        
                        diffusivities_data[d_num_species + 1].resize(1);
                        diffusivities_component_idx[d_num_species + 1].resize(1);
                        
                        // 2/3*mu - mu_v.
                        diffusivities_data[d_num_species + 1][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 1][0] =
                            d_num_species*(d_num_species + 1) + 1;
                        
                        diffusivities_data[d_num_species + 2].resize(0);
                        diffusivities_component_idx[d_num_species + 2].resize(0);
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[d_num_species + 3].resize(2);
                        diffusivities_component_idx[d_num_species + 3].resize(2);
                        
                        // v*(2/3*mu - mu_v).
                        diffusivities_data[d_num_species + 3][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][0] =
                            d_num_species*(d_num_species + 1) + 7;
                        
                        // -u*mu.
                        diffusivities_data[d_num_species + 3][1] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][1] =
                            d_num_species*(d_num_species + 1) + 9;
                        
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            diffusivities_data[si].resize(d_num_species + 1);
                            diffusivities_component_idx[si].resize(d_num_species + 1);
                            
                            diffusivities_data[si][0] = d_data_diffusivities;
                            diffusivities_component_idx[si][0] = si*(d_num_species + 1);
                            
                            for (int sj = 0; sj < d_num_species; sj++)
                            {
                                diffusivities_data[si][1 + sj] = d_data_diffusivities;
                                diffusivities_component_idx[si][1 + sj] = si*(d_num_species + 1) + sj + 1;
                            }
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[d_num_species].resize(1);
                        diffusivities_component_idx[d_num_species].resize(1);
                        
                        // -mu.
                        diffusivities_data[d_num_species][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species][0] =
                            d_num_species*(d_num_species + 1) + 2;
                        
                        diffusivities_data[d_num_species + 1].resize(1);
                        diffusivities_component_idx[d_num_species + 1].resize(1);
                        
                        // -(4/3*mu + mu_v).
                        diffusivities_data[d_num_species + 1][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 1][0] =
                            d_num_species*(d_num_species + 1);
                        
                        diffusivities_data[d_num_species + 2].resize(1);
                        diffusivities_component_idx[d_num_species + 2].resize(1);
                        
                        // -mu.
                        diffusivities_data[d_num_species + 2][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 2][0] =
                            d_num_species*(d_num_species + 1) + 2;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[d_num_species + 3].resize(
                            4 + d_num_species*(d_num_species + 1));
                        diffusivities_component_idx[d_num_species + 3].resize(
                            4 + d_num_species*(d_num_species + 1));
                        
                        // -u*mu.
                        diffusivities_data[d_num_species + 3][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][0] =
                            d_num_species*(d_num_species + 1) + 9;
                        
                        // -v*(4/3*mu + mu_v).
                        diffusivities_data[d_num_species + 3][1] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][1] =
                            d_num_species*(d_num_species + 1) + 4;
                        
                        // -w*mu.
                        diffusivities_data[d_num_species + 3][2] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][2] =
                            d_num_species*(d_num_species + 1) + 11;
                        
                        // -kappa.
                        diffusivities_data[d_num_species + 3][3] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][3] =
                            d_num_species*(d_num_species + 1) + 12;
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            diffusivities_data[d_num_species + 3][4 + si*(d_num_species + 1)] =
                                d_data_diffusivities;
                            diffusivities_component_idx[d_num_species + 3][4 + si*(d_num_species + 1)] =
                                d_num_species*(d_num_species + 1) + 13 + si*(d_num_species + 1);
                            
                            for (int sj = 0; sj < d_num_species; sj++)
                            {
                                diffusivities_data[d_num_species + 3][5 + si*(d_num_species + 1) + sj] =
                                    d_data_diffusivities;
                                diffusivities_component_idx[d_num_species + 3][5 + si*(d_num_species + 1) + sj] =
                                    d_num_species*(d_num_species + 1) + 14 + si*(d_num_species + 1) + sj;
                            }
                        }
                        
                        break;
                    }
                    case DIRECTION::Z_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            diffusivities_data[si].resize(0);
                            diffusivities_component_idx[si].resize(0);
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[d_num_species].resize(0);
                        diffusivities_component_idx[d_num_species].resize(0);
                        
                        diffusivities_data[d_num_species + 1].resize(1);
                        diffusivities_component_idx[d_num_species + 1].resize(1);
                        
                        // 2/3*(mu - mu_v).
                        diffusivities_data[d_num_species + 1][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 1][0] =
                            d_num_species*(d_num_species + 1) + 1;
                        
                        diffusivities_data[d_num_species + 2].resize(1);
                        diffusivities_component_idx[d_num_species + 2].resize(1);
                        
                        // -mu.
                        diffusivities_data[d_num_species + 2][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 2][0] =
                            d_num_species*(d_num_species + 1) + 2;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[d_num_species + 3].resize(2);
                        diffusivities_component_idx[d_num_species + 3].resize(2);
                        
                        // -w*u.
                        diffusivities_data[d_num_species + 3][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][0] =
                            d_num_species*(d_num_species + 1) + 11;
                        
                        // v*(2/3*mu - mu_v).
                        diffusivities_data[d_num_species + 3][1] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][1] =
                            d_num_species*(d_num_species + 1) + 7;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelFourEqnConservative::getDiffusiveFluxDiffusivities()\n"
                            << "There are only x-direction, y-direction and z-direction for three-dimensional problem."
                            << std::endl);
                    }
                }
                
                break;
            }
            case DIRECTION::Z_DIRECTION:
            {
                switch (derivative_direction)
                {
                    case DIRECTION::X_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            diffusivities_data[si].resize(0);
                            diffusivities_component_idx[si].resize(0);
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[d_num_species].resize(1);
                        diffusivities_component_idx[d_num_species].resize(1);
                        
                        // -mu.
                        diffusivities_data[d_num_species][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species][0] =
                            d_num_species*(d_num_species + 1) + 2;
                        
                        diffusivities_data[d_num_species + 1].resize(0);
                        diffusivities_component_idx[d_num_species + 1].resize(0);
                        
                        diffusivities_data[d_num_species + 2].resize(1);
                        diffusivities_component_idx[d_num_species + 2].resize(1);
                        
                        // 2/3*mu - mu_v.
                        diffusivities_data[d_num_species + 2][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 2][0] =
                            d_num_species*(d_num_species + 1) + 1;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[d_num_species + 3].resize(2);
                        diffusivities_component_idx[d_num_species + 3].resize(2);
                        
                        // w*(2/3*mu - mu_v).
                        diffusivities_data[d_num_species + 3][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][0] =
                            d_num_species*(d_num_species + 1) + 8;
                        
                        // -u*mu.
                        diffusivities_data[d_num_species + 3][1] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][1] =
                            d_num_species*(d_num_species + 1) + 9;
                        
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            diffusivities_data[si].resize(0);
                            diffusivities_component_idx[si].resize(0);
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[d_num_species].resize(0);
                        diffusivities_component_idx[d_num_species].resize(0);
                        
                        diffusivities_data[d_num_species + 1].resize(1);
                        diffusivities_component_idx[d_num_species + 1].resize(1);
                        
                        // -mu.
                        diffusivities_data[d_num_species + 1][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 1][0] =
                            d_num_species*(d_num_species + 1) + 2;
                        
                        diffusivities_data[d_num_species + 2].resize(1);
                        diffusivities_component_idx[d_num_species + 2].resize(1);
                        
                        // 2/3*mu - mu_v.
                        diffusivities_data[d_num_species + 2][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 2][0] =
                            d_num_species*(d_num_species + 1) + 1;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[d_num_species + 3].resize(2);
                        diffusivities_component_idx[d_num_species + 3].resize(2);
                        
                        // w*(2/3*mu - mu_v).
                        diffusivities_data[d_num_species + 3][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][0] =
                            d_num_species*(d_num_species + 1) + 8;
                        
                        // -v*mu.
                        diffusivities_data[d_num_species + 3][1] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][1] =
                            d_num_species*(d_num_species + 1) + 10;
                        
                        break;
                    }
                    case DIRECTION::Z_DIRECTION:
                    {
                        /*
                         * Mass equations.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            diffusivities_data[si].resize(d_num_species + 1);
                            diffusivities_component_idx[si].resize(d_num_species + 1);
                            
                            diffusivities_data[si][0] = d_data_diffusivities;
                            diffusivities_component_idx[si][0] = si*(d_num_species + 1);
                            
                            for (int sj = 0; sj < d_num_species; sj++)
                            {
                                diffusivities_data[si][1 + sj] = d_data_diffusivities;
                                diffusivities_component_idx[si][1 + sj] = si*(d_num_species + 1) + sj + 1;
                            }
                        }
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[d_num_species].resize(1);
                        diffusivities_component_idx[d_num_species].resize(1);
                        
                        // -mu.
                        diffusivities_data[d_num_species][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species][0] =
                            d_num_species*(d_num_species + 1) + 2;
                        
                        diffusivities_data[d_num_species + 1].resize(1);
                        diffusivities_component_idx[d_num_species + 1].resize(1);
                        
                        // -mu.
                        diffusivities_data[d_num_species + 1][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 1][0] =
                            d_num_species*(d_num_species + 1) + 2;
                        
                        diffusivities_data[d_num_species + 2].resize(1);
                        diffusivities_component_idx[d_num_species + 2].resize(1);
                        
                        // -(4/3*mu + mu_v).
                        diffusivities_data[d_num_species + 2][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 2][0] =
                            d_num_species*(d_num_species + 1);
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[d_num_species + 3].resize(
                            4 + d_num_species*(d_num_species + 1));
                        diffusivities_component_idx[d_num_species + 3].resize(
                            4 + d_num_species*(d_num_species + 1));
                        
                        // -u*mu.
                        diffusivities_data[d_num_species + 3][0] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][0] =
                            d_num_species*(d_num_species + 1) + 9;
                        
                        // -v*mu.
                        diffusivities_data[d_num_species + 3][1] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][1] =
                            d_num_species*(d_num_species + 1) + 10;
                        
                        // -w*(4/3*mu + mu_v).
                        diffusivities_data[d_num_species + 3][2] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][2] =
                            d_num_species*(d_num_species + 1) + 5;
                        
                        // -kappa.
                        diffusivities_data[d_num_species + 3][3] = d_data_diffusivities;
                        diffusivities_component_idx[d_num_species + 3][3] =
                            d_num_species*(d_num_species + 1) + 12;
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            diffusivities_data[d_num_species + 3][4 + si*(d_num_species + 1)] =
                                d_data_diffusivities;
                            diffusivities_component_idx[d_num_species + 3][4 + si*(d_num_species + 1)] =
                                d_num_species*(d_num_species + 1) + 13 + si*(d_num_species + 1);
                            
                            for (int sj = 0; sj < d_num_species; sj++)
                            {
                                diffusivities_data[d_num_species + 3][5 + si*(d_num_species + 1) + sj] =
                                    d_data_diffusivities;
                                diffusivities_component_idx[d_num_species + 3][5 + si*(d_num_species + 1) + sj] =
                                    d_num_species*(d_num_species + 1) + 14 + si*(d_num_species + 1) + sj;
                            }
                        }
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelFourEqnConservative::getDiffusiveFluxDiffusivities()\n"
                            << "There are only x-direction, y-direction and z-direction for three-dimensional problem."
                            << std::endl);
                    }
                }
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelFourEqnConservative::getDiffusiveFluxDiffusivities()\n"
                    << "There are only x-direction, y-direction and z-direction for three-dimensional problem."
                    << std::endl);
            }
        }
    }
    
    d_global_derived_cell_data_computed = true;
}


/*
 * Compute derived plot quantities registered with the VisIt data writers from data that
 * is maintained on each patch in the hierarchy.
 */
bool
FlowModelFourEqnConservative::packDerivedDataIntoDoubleBuffer(
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
            << ": FlowModelFourEqnConservative::"
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
    
    if (variable_name == "density")
    {
        boost::shared_ptr<pdat::CellData<double> > data_partial_densities(
            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                patch.getPatchData(s_variable_partial_densities, d_plot_context)));
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_partial_densities);
        TBOX_ASSERT(data_partial_densities->getGhostBox().isSpatiallyEqual(patch.getBox()));
#endif
        
        // Get the dimensions of box that covers the data.
        const hier::Box data_box = data_partial_densities->getGhostBox();
        const hier::IntVector data_dims = data_box.numberCells();
        
        // Get the pointers to the conservative variables.
        std::vector<const double*> rho_Y;
        rho_Y.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y.push_back(data_partial_densities->getPointer(si));
        }
        
        size_t offset_data = data_box.offset(region.lower());
        
        if (d_dim == tbox::Dimension(1))
        {
            for (int i = 0; i < region_dims[0]; i++)
            {
                // Compute the linear indices.
                size_t idx_data = offset_data + i;
                size_t idx_region = i;
                
                std::vector<const double*> rho_Y_ptr;
                rho_Y_ptr.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    rho_Y_ptr.push_back(&rho_Y[si][idx_data]);
                }
                
                buffer[idx_region] = d_equation_of_state_mixing_rules->getMixtureDensity(rho_Y_ptr);
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
                    
                    std::vector<const double*> rho_Y_ptr;
                    rho_Y_ptr.reserve(d_num_species);
                    for (int si = 0; si < d_num_species; si++)
                    {
                        rho_Y_ptr.push_back(&rho_Y[si][idx_data]);
                    }
                    
                    buffer[idx_region] = d_equation_of_state_mixing_rules->getMixtureDensity(rho_Y_ptr);
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
                        
                        std::vector<const double*> rho_Y_ptr;
                        rho_Y_ptr.reserve(d_num_species);
                        for (int si = 0; si < d_num_species; si++)
                        {
                            rho_Y_ptr.push_back(&rho_Y[si][idx_data]);
                        }
                        
                        buffer[idx_region] = d_equation_of_state_mixing_rules->getMixtureDensity(rho_Y_ptr);
                    }
                }
            }
        }
        
        data_on_patch = true;
    }
    else if (variable_name == "pressure")
    {
        boost::shared_ptr<pdat::CellData<double> > data_partial_densities(
            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                patch.getPatchData(s_variable_partial_densities, d_plot_context)));
        
        boost::shared_ptr<pdat::CellData<double> > data_momentum(
            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                patch.getPatchData(s_variable_momentum, d_plot_context)));
        
        boost::shared_ptr<pdat::CellData<double> > data_total_energy(
            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                patch.getPatchData(s_variable_total_energy, d_plot_context)));
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_partial_densities);
        TBOX_ASSERT(data_momentum);
        TBOX_ASSERT(data_total_energy);
        TBOX_ASSERT(data_partial_densities->getGhostBox().isSpatiallyEqual(patch.getBox()));
        TBOX_ASSERT(data_momentum->getGhostBox().isSpatiallyEqual(patch.getBox()));
        TBOX_ASSERT(data_total_energy->getGhostBox().isSpatiallyEqual(patch.getBox()));
#endif
        
        // Get the dimensions of box that covers the data.
        const hier::Box data_box = data_partial_densities->getGhostBox();
        const hier::IntVector data_dims = data_box.numberCells();
        
        // Get the pointers to conservative variables
        std::vector<const double*> rho_Y;
        rho_Y.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y.push_back(data_partial_densities->getPointer(si));
        }
        const double* const rho_u = data_momentum->getPointer(0);
        const double* const rho_v = d_dim > tbox::Dimension(1) ? data_momentum->getPointer(1) : NULL;
        const double* const rho_w = d_dim > tbox::Dimension(2) ? data_momentum->getPointer(2) : NULL;
        const double* const E     = data_total_energy->getPointer(0);
        
        size_t offset_data = data_box.offset(region.lower());
        
        if (d_dim == tbox::Dimension(1))
        {
            for (int i = 0; i < region_dims[0]; i++)
            {
                // Compute the linear indices.
                size_t idx_data = offset_data + i;
                size_t idx_region = i;
                
                // Compute the mixture density.
                std::vector<const double*> rho_Y_ptr;
                rho_Y_ptr.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    rho_Y_ptr.push_back(&(rho_Y[si][idx_data]));
                }
                const double rho = d_equation_of_state_mixing_rules->getMixtureDensity(
                    rho_Y_ptr);
                
                /*
                 * Compute the internal energy.
                 */
                double epsilon = double(0);
                if (d_dim == tbox::Dimension(1))
                {
                    epsilon = (E[idx_data] -
                        double(1)/double(2)*rho_u[idx_data]*rho_u[idx_data]/rho)/rho;
                }
                else if (d_dim == tbox::Dimension(2))
                {
                    epsilon = (E[idx_data] -
                        double(1)/double(2)*(rho_u[idx_data]*rho_u[idx_data] +
                        rho_v[idx_data]*rho_v[idx_data])/rho)/rho;
                }
                else if (d_dim == tbox::Dimension(3))
                {
                    epsilon = (E[idx_data] -
                        double(1)/double(2)*(rho_u[idx_data]*rho_u[idx_data] +
                        rho_v[idx_data]*rho_v[idx_data] +
                        rho_w[idx_data]*rho_w[idx_data])/rho)/rho;
                }
                
                /*
                 * Compute the mass fractions.
                 */
                double Y[d_num_species];
                for (int si = 0; si < d_num_species; si++)
                {
                    Y[si] = rho_Y[si][idx_data]/rho;
                }
                
                /*
                 * Get the pointers to the mass fractions.
                 */
                std::vector<const double*> Y_ptr;
                Y_ptr.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    Y_ptr.push_back(&Y[si]);
                }
                
                buffer[idx_region] = d_equation_of_state_mixing_rules->
                    getPressure(
                        &rho,
                        &epsilon,
                        Y_ptr);
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
                    
                    // Compute the mixture density.
                    std::vector<const double*> rho_Y_ptr;
                    rho_Y_ptr.reserve(d_num_species);
                    for (int si = 0; si < d_num_species; si++)
                    {
                        rho_Y_ptr.push_back(&(rho_Y[si][idx_data]));
                    }
                    const double rho = d_equation_of_state_mixing_rules->getMixtureDensity(
                        rho_Y_ptr);
                    
                    /*
                     * Compute the internal energy.
                     */
                    double epsilon = double(0);
                    if (d_dim == tbox::Dimension(1))
                    {
                        epsilon = (E[idx_data] -
                            double(1)/double(2)*rho_u[idx_data]*rho_u[idx_data]/rho)/rho;
                    }
                    else if (d_dim == tbox::Dimension(2))
                    {
                        epsilon = (E[idx_data] -
                            double(1)/double(2)*(rho_u[idx_data]*rho_u[idx_data] +
                            rho_v[idx_data]*rho_v[idx_data])/rho)/rho;
                    }
                    else if (d_dim == tbox::Dimension(3))
                    {
                        epsilon = (E[idx_data] -
                            double(1)/double(2)*(rho_u[idx_data]*rho_u[idx_data] +
                            rho_v[idx_data]*rho_v[idx_data] +
                            rho_w[idx_data]*rho_w[idx_data])/rho)/rho;
                    }
                    
                    /*
                     * Compute the mass fractions.
                     */
                    double Y[d_num_species];
                    for (int si = 0; si < d_num_species; si++)
                    {
                        Y[si] = rho_Y[si][idx_data]/rho;
                    }
                    
                    /*
                     * Get the pointers to the mass fractions.
                     */
                    std::vector<const double*> Y_ptr;
                    Y_ptr.reserve(d_num_species);
                    for (int si = 0; si < d_num_species; si++)
                    {
                        Y_ptr.push_back(&Y[si]);
                    }
                    
                    buffer[idx_region] = d_equation_of_state_mixing_rules->
                        getPressure(
                            &rho,
                            &epsilon,
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
                        // Compute the linear indices.
                        size_t idx_data = offset_data + i +
                            j*data_dims[0] +
                            k*data_dims[0]*data_dims[1];
                        
                        size_t idx_region = i +
                            j*region_dims[0] +
                            k*region_dims[0]*region_dims[1];
                        
                        // Compute the mixture density.
                        std::vector<const double*> rho_Y_ptr;
                        rho_Y_ptr.reserve(d_num_species);
                        for (int si = 0; si < d_num_species; si++)
                        {
                            rho_Y_ptr.push_back(&(rho_Y[si][idx_data]));
                        }
                        const double rho = d_equation_of_state_mixing_rules->getMixtureDensity(
                            rho_Y_ptr);
                        
                        /*
                         * Compute the internal energy.
                         */
                        double epsilon = double(0);
                        if (d_dim == tbox::Dimension(1))
                        {
                            epsilon = (E[idx_data] -
                                double(1)/double(2)*rho_u[idx_data]*rho_u[idx_data]/rho)/rho;
                        }
                        else if (d_dim == tbox::Dimension(2))
                        {
                            epsilon = (E[idx_data] -
                                double(1)/double(2)*(rho_u[idx_data]*rho_u[idx_data] +
                                rho_v[idx_data]*rho_v[idx_data])/rho)/rho;
                        }
                        else if (d_dim == tbox::Dimension(3))
                        {
                            epsilon = (E[idx_data] -
                                double(1)/double(2)*(rho_u[idx_data]*rho_u[idx_data] +
                                rho_v[idx_data]*rho_v[idx_data] +
                                rho_w[idx_data]*rho_w[idx_data])/rho)/rho;
                        }
                        
                        /*
                         * Compute the mass fractions.
                         */
                        double Y[d_num_species];
                        for (int si = 0; si < d_num_species; si++)
                        {
                            Y[si] = rho_Y[si][idx_data]/rho;
                        }
                        
                        /*
                         * Get the pointers to the mass fractions.
                         */
                        std::vector<const double*> Y_ptr;
                        Y_ptr.reserve(d_num_species);
                        for (int si = 0; si < d_num_species; si++)
                        {
                            Y_ptr.push_back(&Y[si]);
                        }
                        
                        buffer[idx_region] = d_equation_of_state_mixing_rules->
                            getPressure(
                                &rho,
                                &epsilon,
                                Y_ptr);
                    }
                }
            }
        }
        
        data_on_patch = true;
    }
    else if (variable_name == "sound speed")
    {
        boost::shared_ptr<pdat::CellData<double> > data_partial_densities(
            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                patch.getPatchData(s_variable_partial_densities, d_plot_context)));
        
        boost::shared_ptr<pdat::CellData<double> > data_momentum(
            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                patch.getPatchData(s_variable_momentum, d_plot_context)));
        
        boost::shared_ptr<pdat::CellData<double> > data_total_energy(
            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                patch.getPatchData(s_variable_total_energy, d_plot_context)));
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_partial_densities);
        TBOX_ASSERT(data_momentum);
        TBOX_ASSERT(data_total_energy);
        TBOX_ASSERT(data_partial_densities->getGhostBox().isSpatiallyEqual(patch.getBox()));
        TBOX_ASSERT(data_momentum->getGhostBox().isSpatiallyEqual(patch.getBox()));
        TBOX_ASSERT(data_total_energy->getGhostBox().isSpatiallyEqual(patch.getBox()));
#endif
        
        // Get the dimensions of box that covers the data.
        const hier::Box data_box = data_partial_densities->getGhostBox();
        const hier::IntVector data_dims = data_box.numberCells();
        
        // Get the pointers to the conservative variables.
        std::vector<const double*> rho_Y;
        rho_Y.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y.push_back(data_partial_densities->getPointer(si));
        }
        const double* const rho_u = data_momentum->getPointer(0);
        const double* const rho_v = d_dim > tbox::Dimension(1) ? data_momentum->getPointer(1) : NULL;
        const double* const rho_w = d_dim > tbox::Dimension(2) ? data_momentum->getPointer(2) : NULL;
        const double* const E     = data_total_energy->getPointer(0);
        
        size_t offset_data = data_box.offset(region.lower());
        
        if (d_dim == tbox::Dimension(1))
        {
            for (int i = 0; i < region_dims[0]; i++)
            {
                // Compute the linear indices.
                size_t idx_data = offset_data + i;
                size_t idx_region = i;
                
                // Compute the mixture density.
                std::vector<const double*> rho_Y_ptr;
                rho_Y_ptr.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    rho_Y_ptr.push_back(&(rho_Y[si][idx_data]));
                }
                const double rho = d_equation_of_state_mixing_rules->getMixtureDensity(
                    rho_Y_ptr);
                
                /*
                 * Compute the internal energy.
                 */
                double epsilon = double(0);
                if (d_dim == tbox::Dimension(1))
                {
                    epsilon = (E[idx_data] -
                        double(1)/double(2)*rho_u[idx_data]*rho_u[idx_data]/rho)/rho;
                }
                else if (d_dim == tbox::Dimension(2))
                {
                    epsilon = (E[idx_data] -
                        double(1)/double(2)*(rho_u[idx_data]*rho_u[idx_data] +
                        rho_v[idx_data]*rho_v[idx_data])/rho)/rho;
                }
                else if (d_dim == tbox::Dimension(3))
                {
                    epsilon = (E[idx_data] -
                        double(1)/double(2)*(rho_u[idx_data]*rho_u[idx_data] +
                        rho_v[idx_data]*rho_v[idx_data] +
                        rho_w[idx_data]*rho_w[idx_data])/rho)/rho;
                }
                
                /*
                 * Compute the mass fractions.
                 */
                double Y[d_num_species];
                for (int si = 0; si < d_num_species; si++)
                {
                    Y[si] = rho_Y[si][idx_data]/rho;
                }
                
                /*
                 * Get the pointers to the mass fractions.
                 */
                std::vector<const double*> Y_ptr;
                Y_ptr.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    Y_ptr.push_back(&Y[si]);
                }
                
                const double p = d_equation_of_state_mixing_rules->
                    getPressure(
                        &rho,
                        &epsilon,
                        Y_ptr);
                
                buffer[idx_region] = d_equation_of_state_mixing_rules->
                    getSoundSpeed(
                        &rho,
                        &p,
                        Y_ptr);
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
                    
                    // Compute the mixture density.
                    std::vector<const double*> rho_Y_ptr;
                    rho_Y_ptr.reserve(d_num_species);
                    for (int si = 0; si < d_num_species; si++)
                    {
                        rho_Y_ptr.push_back(&(rho_Y[si][idx_data]));
                    }
                    const double rho = d_equation_of_state_mixing_rules->getMixtureDensity(
                        rho_Y_ptr);
                    
                    /*
                     * Compute the internal energy.
                     */
                    double epsilon = double(0);
                    if (d_dim == tbox::Dimension(1))
                    {
                        epsilon = (E[idx_data] -
                            double(1)/double(2)*rho_u[idx_data]*rho_u[idx_data]/rho)/rho;
                    }
                    else if (d_dim == tbox::Dimension(2))
                    {
                        epsilon = (E[idx_data] -
                            double(1)/double(2)*(rho_u[idx_data]*rho_u[idx_data] +
                            rho_v[idx_data]*rho_v[idx_data])/rho)/rho;
                    }
                    else if (d_dim == tbox::Dimension(3))
                    {
                        epsilon = (E[idx_data] -
                            double(1)/double(2)*(rho_u[idx_data]*rho_u[idx_data] +
                            rho_v[idx_data]*rho_v[idx_data] +
                            rho_w[idx_data]*rho_w[idx_data])/rho)/rho;
                    }
                    
                    /*
                     * Compute the mass fractions.
                     */
                    double Y[d_num_species];
                    for (int si = 0; si < d_num_species; si++)
                    {
                        Y[si] = rho_Y[si][idx_data]/rho;
                    }
                    
                    /*
                     * Get the pointers to the mass fractions.
                     */
                    std::vector<const double*> Y_ptr;
                    Y_ptr.reserve(d_num_species);
                    for (int si = 0; si < d_num_species; si++)
                    {
                        Y_ptr.push_back(&Y[si]);
                    }
                    
                    const double p = d_equation_of_state_mixing_rules->
                        getPressure(
                            &rho,
                            &epsilon,
                            Y_ptr);
                    
                    buffer[idx_region] = d_equation_of_state_mixing_rules->
                        getSoundSpeed(
                            &rho,
                            &p,
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
                        // Compute the linear indices.
                        size_t idx_data = offset_data + i +
                            j*data_dims[0] +
                            k*data_dims[0]*data_dims[1];
                        
                        size_t idx_region = i +
                            j*region_dims[0] +
                            k*region_dims[0]*region_dims[1];
                        
                        // Compute the mixture density.
                        std::vector<const double*> rho_Y_ptr;
                        rho_Y_ptr.reserve(d_num_species);
                        for (int si = 0; si < d_num_species; si++)
                        {
                            rho_Y_ptr.push_back(&(rho_Y[si][idx_data]));
                        }
                        const double rho = d_equation_of_state_mixing_rules->getMixtureDensity(
                            rho_Y_ptr);
                        
                        /*
                         * Compute the internal energy.
                         */
                        double epsilon = double(0);
                        if (d_dim == tbox::Dimension(1))
                        {
                            epsilon = (E[idx_data] -
                                double(1)/double(2)*rho_u[idx_data]*rho_u[idx_data]/rho)/rho;
                        }
                        else if (d_dim == tbox::Dimension(2))
                        {
                            epsilon = (E[idx_data] -
                                double(1)/double(2)*(rho_u[idx_data]*rho_u[idx_data] +
                                rho_v[idx_data]*rho_v[idx_data])/rho)/rho;
                        }
                        else if (d_dim == tbox::Dimension(3))
                        {
                            epsilon = (E[idx_data] -
                                double(1)/double(2)*(rho_u[idx_data]*rho_u[idx_data] +
                                rho_v[idx_data]*rho_v[idx_data] +
                                rho_w[idx_data]*rho_w[idx_data])/rho)/rho;
                        }
                        
                        /*
                         * Compute the mass fractions.
                         */
                        double Y[d_num_species];
                        for (int si = 0; si < d_num_species; si++)
                        {
                            Y[si] = rho_Y[si][idx_data]/rho;
                        }
                        
                        /*
                         * Get the pointers to the mass fractions.
                         */
                        std::vector<const double*> Y_ptr;
                        Y_ptr.reserve(d_num_species);
                        for (int si = 0; si < d_num_species; si++)
                        {
                            Y_ptr.push_back(&Y[si]);
                        }
                        
                        const double p = d_equation_of_state_mixing_rules->
                            getPressure(
                                &rho,
                                &epsilon,
                                Y_ptr);
                        
                        buffer[idx_region] = d_equation_of_state_mixing_rules->
                            getSoundSpeed(
                                &rho,
                                &p,
                                Y_ptr);
                    }
                }
            }
        }
        
        data_on_patch = true;
    }
    else if (variable_name == "velocity")
    {
        boost::shared_ptr<pdat::CellData<double> > data_partial_densities(
            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                patch.getPatchData(s_variable_partial_densities, d_plot_context)));
        
        boost::shared_ptr<pdat::CellData<double> > data_momentum(
            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                patch.getPatchData(s_variable_momentum, d_plot_context)));
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_partial_densities);
        TBOX_ASSERT(data_momentum);
        TBOX_ASSERT(data_partial_densities->getGhostBox().isSpatiallyEqual(patch.getBox()));
        TBOX_ASSERT(data_momentum->getGhostBox().isSpatiallyEqual(patch.getBox()));
#endif
        
        // Get the dimensions of box that covers the data.
        const hier::Box data_box = data_partial_densities->getGhostBox();
        const hier::IntVector data_dims = data_box.numberCells();
        
        // Get the pointers to the conservative variables.
        std::vector<const double*> rho_Y;
        rho_Y.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y.push_back(data_partial_densities->getPointer(si));
        }
        const double* const m = data_momentum->getPointer(depth_id);
        
        size_t offset_data = data_box.offset(region.lower());
        
        if (d_dim == tbox::Dimension(1))
        {
            for (int i = 0; i < region_dims[0]; i++)
            {
                // Compute the linear indices.
                size_t idx_data = offset_data + i;
                size_t idx_region = i;
                
                std::vector<const double*> rho_Y_ptr;
                rho_Y_ptr.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    rho_Y_ptr.push_back(&rho_Y[si][idx_data]);
                }
                
                double rho = d_equation_of_state_mixing_rules->getMixtureDensity(rho_Y_ptr);
                
                buffer[idx_region] = m[idx_data]/rho;
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
                    
                    std::vector<const double*> rho_Y_ptr;
                    rho_Y_ptr.reserve(d_num_species);
                    for (int si = 0; si < d_num_species; si++)
                    {
                        rho_Y_ptr.push_back(&rho_Y[si][idx_data]);
                    }
                    
                    double rho = d_equation_of_state_mixing_rules->getMixtureDensity(rho_Y_ptr);
                    
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
                        // Compute the linear indices.
                        size_t idx_data = offset_data + i +
                            j*data_dims[0] +
                            k*data_dims[0]*data_dims[1];
                        
                        size_t idx_region = i +
                            j*region_dims[0] +
                            k*region_dims[0]*region_dims[1];
                        
                        std::vector<const double*> rho_Y_ptr;
                        rho_Y_ptr.reserve(d_num_species);
                        for (int si = 0; si < d_num_species; si++)
                        {
                            rho_Y_ptr.push_back(&rho_Y[si][idx_data]);
                        }
                        
                        double rho = d_equation_of_state_mixing_rules->getMixtureDensity(rho_Y_ptr);
                        
                        buffer[idx_region] = m[idx_data]/rho;
                    }
                }
            }
        }
        
        data_on_patch = true;
    }
    else if (variable_name.find("mass fraction") != std::string::npos)
    {
        int species_idx = std::stoi(variable_name.substr(14));
        
        boost::shared_ptr<pdat::CellData<double> > data_partial_densities(
            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                patch.getPatchData(s_variable_partial_densities, d_plot_context)));
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_partial_densities);
        TBOX_ASSERT(data_partial_densities->getGhostBox().isSpatiallyEqual(patch.getBox()));
#endif
        
        // Get the dimensions of box that covers the data.
        const hier::Box data_box = data_partial_densities->getGhostBox();
        const hier::IntVector data_dims = data_box.numberCells();
        
        // Get the pointers to conservative variables.
        std::vector<const double*> rho_Y;
        rho_Y.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y.push_back(data_partial_densities->getPointer(si));
        }
        
        size_t offset_data = data_box.offset(region.lower());
        
        if (d_dim == tbox::Dimension(1))
        {
            for (int i = 0; i < region_dims[0]; i++)
            {
                // Compute the linear indices.
                size_t idx_data = offset_data + i;
                
                size_t idx_region = i;
                
                std::vector<const double*> rho_Y_ptr;
                rho_Y_ptr.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    rho_Y_ptr.push_back(&rho_Y[si][idx_data]);
                }
                
                double rho = d_equation_of_state_mixing_rules->getMixtureDensity(rho_Y_ptr);
                
                buffer[idx_region] = rho_Y[species_idx][idx_data]/rho;
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
                    
                    std::vector<const double*> rho_Y_ptr;
                    rho_Y_ptr.reserve(d_num_species);
                    for (int si = 0; si < d_num_species; si++)
                    {
                        rho_Y_ptr.push_back(&rho_Y[si][idx_data]);
                    }
                    
                    double rho = d_equation_of_state_mixing_rules->getMixtureDensity(rho_Y_ptr);
                    
                    buffer[idx_region] = rho_Y[species_idx][idx_data]/rho;
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
                        
                        std::vector<const double*> rho_Y_ptr;
                        rho_Y_ptr.reserve(d_num_species);
                        for (int si = 0; si < d_num_species; si++)
                        {
                            rho_Y_ptr.push_back(&rho_Y[si][idx_data]);
                        }
                        
                        double rho = d_equation_of_state_mixing_rules->getMixtureDensity(rho_Y_ptr);
                        
                        buffer[idx_region] = rho_Y[species_idx][idx_data]/rho;
                    }
                }
            }
        }
        
        data_on_patch = true;
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::"
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
FlowModelFourEqnConservative::registerPlotQuantities(
    const boost::shared_ptr<ExtendedVisItDataWriter>& visit_writer)
{
    if (!d_plot_context)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::"
            << "registerPlotQuantities()\n"
            << "The plotting context is not set yet."
            << std::endl);
    }
    
    /*
    hier::VariableDatabase* vardb = hier::VariableDatabase::getDatabase();
    
    for (int si = 0; si < d_num_species; si++)
    {
        std::string partial_densities_name =
            "partial density " + tbox::Utilities::intToString(si);
        
        visit_writer->registerPlotQuantity(
            partial_densities_name,
            "SCALAR",
            vardb->mapVariableAndContextToIndex(
                s_variable_partial_densities,
                d_plot_context),
            si);
    }
    
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
    
    visit_writer->registerDerivedPlotQuantity("density",
        "SCALAR",
        this);
    
    for (int si = 0; si < d_num_species - 1; si++)
    {
        std::string mass_fractions_name =
            "mass fraction " + tbox::Utilities::intToString(si);
            
        visit_writer->registerDerivedPlotQuantity(mass_fractions_name,
            "SCALAR",
            this);
    }
}
#endif


/*
 * Set the number of sub-ghost cells of a variable.
 * This function can be called recursively if the variables are computed recursively.
 */
void
FlowModelFourEqnConservative::setNumberOfSubGhosts(
    const hier::IntVector& num_subghosts,
    const std::string& variable_name,
    const std::string& parent_variable_name)
{
    NULL_USE(parent_variable_name);
    
    if (variable_name == "DENSITY")
    {
        if (d_num_subghosts_density > -hier::IntVector::getOne(d_dim))
        {
            if (num_subghosts > d_num_subghosts_density)
            {
                /*
                TBOX_ERROR(d_object_name
                    << ": FlowModelFourEqnConservative::setNumberOfSubGhosts()\n"
                    << "Number of ghosts of '"
                    << parent_variable_name
                    << "' exceeds"
                    << " number of ghosts of '"
                    << variable_name
                    << "'."
                    << std::endl);
                */
                
                d_num_subghosts_density = num_subghosts;
            }
        }
        else
        {
            d_num_subghosts_density = num_subghosts;
        }
    }
    else if (variable_name == "MASS_FRACTIONS")
    {
        if (d_num_subghosts_mass_fractions > -hier::IntVector::getOne(d_dim))
        {
            if (num_subghosts > d_num_subghosts_mass_fractions)
            {
                /*
                TBOX_ERROR(d_object_name
                    << ": FlowModelFourEqnConservative::setNumberOfSubGhosts()\n"
                    << "Number of ghosts of '"
                    << parent_variable_name
                    << "' exceeds"
                    << " number of ghosts of '"
                    << variable_name
                    << "'."
                    << std::endl);
                */
                
                d_num_subghosts_mass_fractions = num_subghosts;
            }
        }
        else
        {
            d_num_subghosts_mass_fractions = num_subghosts;
        }
        
        setNumberOfSubGhosts(num_subghosts, "DENSITY", parent_variable_name);
    }
    else if (variable_name == "VELOCITY")
    {
        if (d_num_subghosts_velocity > -hier::IntVector::getOne(d_dim))
        {
            if (num_subghosts > d_num_subghosts_velocity)
            {
                /*
                TBOX_ERROR(d_object_name
                    << ": FlowModelFourEqnConservative::setNumberOfSubGhosts()\n"
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
        
        setNumberOfSubGhosts(num_subghosts, "DENSITY", parent_variable_name);
    }
    else if (variable_name == "INTERNAL_ENERGY")
    {
        if (d_num_subghosts_internal_energy > -hier::IntVector::getOne(d_dim))
        {
            if (num_subghosts > d_num_subghosts_internal_energy)
            {
                /*
                TBOX_ERROR(d_object_name
                    << ": FlowModelFourEqnConservative::setNumberOfSubGhosts()\n"
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
        
        setNumberOfSubGhosts(num_subghosts, "DENSITY", parent_variable_name);
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
                    << ": FlowModelFourEqnConservative::setNumberOfSubGhosts()\n"
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
        
        setNumberOfSubGhosts(num_subghosts, "DENSITY", parent_variable_name);
        setNumberOfSubGhosts(num_subghosts, "MASS_FRACTIONS", parent_variable_name);
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
                    << ": FlowModelFourEqnConservative::setNumberOfSubGhosts()\n"
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
        
        setNumberOfSubGhosts(num_subghosts, "DENSITY", parent_variable_name);
        setNumberOfSubGhosts(num_subghosts, "MASS_FRACTIONS", parent_variable_name);
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
                    << ": FlowModelFourEqnConservative::setNumberOfSubGhosts()\n"
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
        
        setNumberOfSubGhosts(num_subghosts, "DENSITY", parent_variable_name);
        setNumberOfSubGhosts(num_subghosts, "MASS_FRACTIONS", parent_variable_name);
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
                    << ": FlowModelFourEqnConservative::setNumberOfSubGhosts()\n"
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
                    << ": FlowModelFourEqnConservative::setNumberOfSubGhosts()\n"
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
                    << ": FlowModelFourEqnConservative::setNumberOfSubGhosts()\n"
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
                    << ": FlowModelFourEqnConservative::setNumberOfSubGhosts()\n"
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
                    << ": FlowModelFourEqnConservative::setNumberOfSubGhosts()\n"
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
                    << ": FlowModelFourEqnConservative::setNumberOfSubGhosts()\n"
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
                    << ": FlowModelFourEqnConservative::setNumberOfSubGhosts()\n"
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
        
        setNumberOfSubGhosts(num_subghosts, "DENSITY", parent_variable_name);
        setNumberOfSubGhosts(num_subghosts, "MASS_FRACTIONS", parent_variable_name);
        setNumberOfSubGhosts(num_subghosts, "PRESSURE", parent_variable_name);
        setNumberOfSubGhosts(num_subghosts, "TEMPERATURE", parent_variable_name);
    }
}


/*
 * Set the ghost boxes and their dimensions of derived cell variables.
 */
void
FlowModelFourEqnConservative::setGhostBoxesAndDimensionsDerivedCellVariables()
{
    if (d_num_subghosts_density > -hier::IntVector::getOne(d_dim))
    {
        d_subghost_box_density = d_interior_box;
        d_subghost_box_density.grow(d_num_subghosts_density);
        d_subghostcell_dims_density = d_subghost_box_density.numberCells();
    }
    
    if (d_num_subghosts_mass_fractions > -hier::IntVector::getOne(d_dim))
    {
        d_subghost_box_mass_fractions = d_interior_box;
        d_subghost_box_mass_fractions.grow(d_num_subghosts_mass_fractions);
        d_subghostcell_dims_mass_fractions = d_subghost_box_mass_fractions.numberCells();
    }
    
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
    
    if (d_num_subghosts_diffusivities > -hier::IntVector::getOne(d_dim))
    {
        d_subghost_box_diffusivities = d_interior_box;
        d_subghost_box_diffusivities.grow(d_num_subghosts_diffusivities);
        d_subghostcell_dims_diffusivities = d_subghost_box_diffusivities.numberCells();
    }
}


/*
 * Get the global cell data of partial densities in the registered patch.
 */
boost::shared_ptr<pdat::CellData<double> >
FlowModelFourEqnConservative::getGlobalCellDataPartialDensities()
{
    // Get the cell data of the registered variable partial densities.
    boost::shared_ptr<pdat::CellData<double> > data_partial_densities(
        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
            d_patch->getPatchData(s_variable_partial_densities, getDataContext())));
    
    return data_partial_densities;
}


/*
 * Get the global cell data of momentum in the registered patch.
 */
boost::shared_ptr<pdat::CellData<double> >
FlowModelFourEqnConservative::getGlobalCellDataMomentum()
{
    // Get the cell data of the registered variable momentum.
    boost::shared_ptr<pdat::CellData<double> > data_momentum(
        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
            d_patch->getPatchData(s_variable_momentum, getDataContext())));
    
    return data_momentum;
}


/*
 * Get the global cell data of total energy in the registered patch.
 */
boost::shared_ptr<pdat::CellData<double> >
FlowModelFourEqnConservative::getGlobalCellDataTotalEnergy()
{
    // Get the cell data of the registered variable total energy.
    boost::shared_ptr<pdat::CellData<double> > data_total_energy(
        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
            d_patch->getPatchData(s_variable_total_energy, getDataContext())));
    
    return data_total_energy;
}


/*
 * Compute the global cell data of density in the registered patch.
 */
void
FlowModelFourEqnConservative::computeGlobalCellDataDensity(
    const hier::Box& domain)
{
    if (d_num_subghosts_density > -hier::IntVector::getOne(d_dim))
    {
        // Create the cell data of density.
        d_data_density.reset(
            new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_density));
        
        // Get the cell data of the variable partial densities.
        boost::shared_ptr<pdat::CellData<double> > data_partial_densities =
            getGlobalCellDataPartialDensities();
        
        // Compute the density field.
        d_equation_of_state_mixing_rules->computeMixtureDensity(
            d_data_density,
            data_partial_densities,
            domain);
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::computeGlobalCellDataDensity()\n"
            << "Cell data of 'DENSITY' is not yet registered."
            << std::endl);
    }
}


/*
 * Compute the global cell data of mass fractions with density in the registered patch.
 */
void
FlowModelFourEqnConservative::computeGlobalCellDataMassFractionsWithDensity(
    const hier::Box& domain)
{
    if (d_num_subghosts_mass_fractions > -hier::IntVector::getOne(d_dim))
    {
        // Create the cell data of mass fractions.
        d_data_mass_fractions.reset(
            new pdat::CellData<double>(d_interior_box, d_num_species, d_num_subghosts_mass_fractions));
        
        /*
         * Get the local lower indices and number of cells in each direction of the domain.
         */
        
        hier::IntVector domain_lo(d_dim);
        hier::IntVector domain_dims(d_dim);
        
        if (domain.empty())
        {
            domain_lo = -d_num_subghosts_mass_fractions;
            domain_dims = d_subghostcell_dims_mass_fractions;
        }
        else
        {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(d_subghost_box_mass_fractions.contains(domain));
#endif
            
            domain_lo = domain.lower() - d_interior_box.lower();
            domain_dims = domain.numberCells();
        }
        
        // Get the cell data of the variable partial densities.
        boost::shared_ptr<pdat::CellData<double> > data_partial_densities =
            getGlobalCellDataPartialDensities();
        
        if (!d_data_density)
        {
            computeGlobalCellDataDensity(domain);
        }
        
        // Get the pointers to the cell data of mass fractions, denisty and partial densities.
        std::vector<double*> Y;
        Y.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            Y.push_back(d_data_mass_fractions->getPointer(si));
        }
        double* rho = d_data_density->getPointer(0);
        std::vector<double*> rho_Y;
        rho_Y.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y.push_back(data_partial_densities->getPointer(si));
        }
        
        if (d_dim == tbox::Dimension(1))
        {
            /*
             * Get the local lower index, numbers of cells in each dimension and numbers of ghost cells.
             */
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_dim_0 = domain_dims[0];
            
            const int num_ghosts_0 = d_num_ghosts[0];
            const int num_subghosts_0_density = d_num_subghosts_density[0];
            const int num_subghosts_0_mass_fractions = d_num_subghosts_mass_fractions[0];
            
            // Compute the mass fraction field.
            for (int si = 0; si < d_num_species; si++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx = i + num_ghosts_0;
                    const int idx_density = i + num_subghosts_0_density;
                    const int idx_mass_fractions = i + num_subghosts_0_mass_fractions;
                    
                    Y[si][idx_mass_fractions] = rho_Y[si][idx]/rho[idx_density];
                }
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
            
            const int num_subghosts_0_density = d_num_subghosts_density[0];
            const int num_subghosts_1_density = d_num_subghosts_density[1];
            const int subghostcell_dim_0_density = d_subghostcell_dims_density[0];
            
            const int num_subghosts_0_mass_fractions = d_num_subghosts_mass_fractions[0];
            const int num_subghosts_1_mass_fractions = d_num_subghosts_mass_fractions[1];
            const int subghostcell_dim_0_mass_fractions = d_subghostcell_dims_mass_fractions[0];
            
            // Compute the mass fraction field.
            for (int si = 0; si < d_num_species; si++)
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
                            (j + num_ghosts_1)*ghostcell_dim_0;
                        
                        const int idx_density = (i + num_subghosts_0_density) +
                            (j + num_subghosts_1_density)*subghostcell_dim_0_density;
                        
                        const int idx_mass_fractions = (i + num_subghosts_0_mass_fractions) +
                            (j + num_subghosts_1_mass_fractions)*subghostcell_dim_0_mass_fractions;
                        
                        Y[si][idx_mass_fractions] = rho_Y[si][idx]/rho[idx_density];
                    }
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
            
            const int num_subghosts_0_density = d_num_subghosts_density[0];
            const int num_subghosts_1_density = d_num_subghosts_density[1];
            const int num_subghosts_2_density = d_num_subghosts_density[2];
            const int subghostcell_dim_0_density = d_subghostcell_dims_density[0];
            const int subghostcell_dim_1_density = d_subghostcell_dims_density[1];
            
            const int num_subghosts_0_mass_fractions = d_num_subghosts_mass_fractions[0];
            const int num_subghosts_1_mass_fractions = d_num_subghosts_mass_fractions[1];
            const int num_subghosts_2_mass_fractions = d_num_subghosts_mass_fractions[2];
            const int subghostcell_dim_0_mass_fractions = d_subghostcell_dims_mass_fractions[0];
            const int subghostcell_dim_1_mass_fractions = d_subghostcell_dims_mass_fractions[1];
            
            // Compute the mass fraction field.
            for (int si = 0; si < d_num_species; si++)
            {
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
                            
                            const int idx_density = (i + num_subghosts_0_density) +
                                (j + num_subghosts_1_density)*subghostcell_dim_0_density +
                                (k + num_subghosts_2_density)*subghostcell_dim_0_density*
                                    subghostcell_dim_1_density;
                            
                            const int idx_mass_fractions = (i + num_subghosts_0_mass_fractions) +
                                (j + num_subghosts_1_mass_fractions)*subghostcell_dim_0_mass_fractions +
                                (k + num_subghosts_2_mass_fractions)*subghostcell_dim_0_mass_fractions*
                                    subghostcell_dim_1_mass_fractions;
                            
                            Y[si][idx_mass_fractions] = rho_Y[si][idx]/rho[idx_density];
                        }
                    }
                }
            }
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::computeGlobalCellDataMassFractionsWithDensity()\n"
            << "Cell data of 'MASS_FRACTION' is not yet registered."
            << std::endl);
    }
}


/*
 * Compute the global cell data of velocity with density in the registered patch.
 */
void
FlowModelFourEqnConservative::computeGlobalCellDataVelocityWithDensity(
    const hier::Box& domain)
{
    if (d_num_subghosts_velocity > -hier::IntVector::getOne(d_dim))
    {
        // Create the cell data of velocity.
        d_data_velocity.reset(
            new pdat::CellData<double>(d_interior_box, d_dim.getValue(), d_num_subghosts_velocity));
        
        /*
         * Get the local lower indices and number of cells in each direction of the domain.
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
        
        // Get the cell data of the variable momentum.
        boost::shared_ptr<pdat::CellData<double> > data_momentum =
            getGlobalCellDataMomentum();
        
        if (!d_data_density)
        {
            computeGlobalCellDataDensity(domain);
        }
        
        // Get the pointer to the cell data of density.
        double* rho = d_data_density->getPointer(0);
        
        if (d_dim == tbox::Dimension(1))
        {
            /*
             * Get the local lower index, numbers of cells in each dimension and numbers of ghost cells.
             */
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_dim_0 = domain_dims[0];
            
            const int num_ghosts_0 = d_num_ghosts[0];
            const int num_subghosts_0_density = d_num_subghosts_density[0];
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
                const int idx_density = i + num_subghosts_0_density;
                const int idx_velocity = i + num_subghosts_0_velocity;
                
                u[idx_velocity] = rho_u[idx]/rho[idx_density];
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
            
            const int num_subghosts_0_density = d_num_subghosts_density[0];
            const int num_subghosts_1_density = d_num_subghosts_density[1];
            const int subghostcell_dim_0_density = d_subghostcell_dims_density[0];
            
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
                    
                    const int idx_density = (i + num_subghosts_0_density) +
                        (j + num_subghosts_1_density)*subghostcell_dim_0_density;
                    
                    const int idx_velocity = (i + num_subghosts_0_velocity) +
                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                    
                    u[idx_velocity] = rho_u[idx]/rho[idx_density];
                    v[idx_velocity] = rho_v[idx]/rho[idx_density];
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
            
            const int num_subghosts_0_density = d_num_subghosts_density[0];
            const int num_subghosts_1_density = d_num_subghosts_density[1];
            const int num_subghosts_2_density = d_num_subghosts_density[2];
            const int subghostcell_dim_0_density = d_subghostcell_dims_density[0];
            const int subghostcell_dim_1_density = d_subghostcell_dims_density[1];
            
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
                        
                        const int idx_density = (i + num_subghosts_0_density) +
                            (j + num_subghosts_1_density)*subghostcell_dim_0_density +
                            (k + num_subghosts_2_density)*subghostcell_dim_0_density*
                                subghostcell_dim_1_density;
                        
                        const int idx_velocity = (i + num_subghosts_0_velocity) +
                            (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                            (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                subghostcell_dim_1_velocity;
                        
                        u[idx_velocity] = rho_u[idx]/rho[idx_density];
                        v[idx_velocity] = rho_v[idx]/rho[idx_density];
                        w[idx_velocity] = rho_w[idx]/rho[idx_density];
                    }
                }
            }
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::computeGlobalCellDataVelocityWithDensity()\n"
            << "Cell data of 'VELOCITY' is not yet registered."
            << std::endl);
    }
}


/*
 * Compute the global cell data of internal energy with density and velocity in the registered
 * patch.
 */
void
FlowModelFourEqnConservative::computeGlobalCellDataInternalEnergyWithDensityAndVelocity(
    const hier::Box& domain)
{
    if (d_num_subghosts_internal_energy > -hier::IntVector::getOne(d_dim))
    {
        // Create the cell data of internal energy.
        d_data_internal_energy.reset(
            new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_internal_energy));
        
        /*
         * Get the local lower indices and number of cells in each direction of the domain.
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
        
        // Get the cell data of the variable total energy.
        boost::shared_ptr<pdat::CellData<double> > data_total_energy =
            getGlobalCellDataTotalEnergy();
        
        if (!d_data_density)
        {
            computeGlobalCellDataDensity(domain);
        }
        
        if (!d_data_velocity)
        {
            computeGlobalCellDataVelocityWithDensity(domain);
        }
        
        // Get the pointers to the cell data of internal energy, total energy and density.
        double* epsilon = d_data_internal_energy->getPointer(0);
        double* E = data_total_energy->getPointer(0);
        double* rho = d_data_density->getPointer(0);
        
        if (d_dim == tbox::Dimension(1))
        {
            /*
             * Get the local lower index, numbers of cells in each dimension and numbers of ghost cells.
             */
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_dim_0 = domain_dims[0];
            
            const int num_ghosts_0 = d_num_ghosts[0];
            const int num_subghosts_0_density = d_num_subghosts_density[0];
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
                const int idx_density = i + num_subghosts_0_density;
                const int idx_velocity = i + num_subghosts_0_velocity;
                const int idx_internal_energy = i + num_subghosts_0_internal_energy;
                
                epsilon[idx_internal_energy] = E[idx]/rho[idx_density] -
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
            
            const int num_subghosts_0_density = d_num_subghosts_density[0];
            const int num_subghosts_1_density = d_num_subghosts_density[1];
            const int subghostcell_dim_0_density = d_subghostcell_dims_density[0];
            
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
                    
                    const int idx_density = (i + num_subghosts_0_density) +
                        (j + num_subghosts_1_density)*subghostcell_dim_0_density;
                    
                    const int idx_velocity = (i + num_subghosts_0_velocity) +
                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                    
                    const int idx_internal_energy = (i + num_subghosts_0_internal_energy) +
                        (j + num_subghosts_1_internal_energy)*subghostcell_dim_0_internal_energy;
                    
                    epsilon[idx_internal_energy] = E[idx]/rho[idx_density] -
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
            
            const int num_subghosts_0_density = d_num_subghosts_density[0];
            const int num_subghosts_1_density = d_num_subghosts_density[1];
            const int num_subghosts_2_density = d_num_subghosts_density[2];
            const int subghostcell_dim_0_density = d_subghostcell_dims_density[0];
            const int subghostcell_dim_1_density = d_subghostcell_dims_density[1];
            
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
                        
                        const int idx_density = (i + num_subghosts_0_density) +
                            (j + num_subghosts_1_density)*subghostcell_dim_0_density +
                            (k + num_subghosts_2_density)*subghostcell_dim_0_density*
                                subghostcell_dim_1_density;
                        
                        const int idx_velocity = (i + num_subghosts_0_velocity) +
                            (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                            (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                subghostcell_dim_1_velocity;
                        
                        const int idx_internal_energy = (i + num_subghosts_0_internal_energy) +
                            (j + num_subghosts_1_internal_energy)*subghostcell_dim_0_internal_energy +
                            (k + num_subghosts_2_internal_energy)*subghostcell_dim_0_internal_energy*
                                subghostcell_dim_1_internal_energy;
                        
                        epsilon[idx_internal_energy] = E[idx]/rho[idx_density] -
                            double(1)/double(2)*(u[idx_velocity]*u[idx_velocity] + v[idx_velocity]*v[idx_velocity] +
                                w[idx_velocity]*w[idx_velocity]);
                    }
                }
            }
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::"
            << "computeGlobalCellDataInternalEnergyWithDensityAndVelocity()\n"
            << "Cell data of 'INTERNAL_ENERGY' is not yet registered."
            << std::endl);
    }
}


/*
 * Compute the global cell data of pressure with density, mass fractions and internal energy in
 * the registered patch.
 */
void
FlowModelFourEqnConservative::computeGlobalCellDataPressureWithDensityMassFractionsAndInternalEnergy(
    const hier::Box& domain)
{
    if (d_num_subghosts_pressure > -hier::IntVector::getOne(d_dim))
    {
        // Create the cell data of pressure.
        d_data_pressure.reset(
            new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_pressure));
        
        if (!d_data_density)
        {
            computeGlobalCellDataDensity(domain);
        }
        
        if (!d_data_mass_fractions)
        {
            computeGlobalCellDataMassFractionsWithDensity(domain);
        }
        
        if (!d_data_internal_energy)
        {
            computeGlobalCellDataInternalEnergyWithDensityAndVelocity(domain);
        }
        
        // Compute the pressure field.
        d_equation_of_state_mixing_rules->computePressure(
            d_data_pressure,
            d_data_density,
            d_data_internal_energy,
            d_data_mass_fractions,
            domain);
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::"
            << "computeGlobalCellDataPressureWithDensityMassFractionsAndInternalEnergy()\n"
            << "Cell data of 'PRESSURE' is not yet registered."
            << std::endl);
    }
}


/*
 * Compute the global cell data of sound speed with density, mass fractions and pressure in the
 * registered patch.
 */
void
FlowModelFourEqnConservative::computeGlobalCellDataSoundSpeedWithDensityMassFractionsAndPressure(
    const hier::Box& domain)
{
    if (d_num_subghosts_sound_speed > -hier::IntVector::getOne(d_dim))
    {
        // Create the cell data of sound speed.
        d_data_sound_speed.reset(
            new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_sound_speed));
        
        if (!d_data_density)
        {
            computeGlobalCellDataDensity(domain);
        }
        
        if (!d_data_mass_fractions)
        {
            computeGlobalCellDataMassFractionsWithDensity(domain);
        }
        
        if (!d_data_pressure)
        {
            computeGlobalCellDataPressureWithDensityMassFractionsAndInternalEnergy(domain);
        }
        
        // Compute the sound speed field.
        d_equation_of_state_mixing_rules->computeSoundSpeed(
            d_data_sound_speed,
            d_data_density,
            d_data_pressure,
            d_data_mass_fractions,
            domain);
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::"
            << "computeGlobalCellDataSoundSpeedWithDensityMassFractionsAndPressure()\n"
            << "Cell data of 'SOUND_SPEED' is not yet registered."
            << std::endl);
    }
}


/*
 * Compute the global cell data of temperature with density, mass fractions and pressure in the
 * registered patch.
 */
void
FlowModelFourEqnConservative::computeGlobalCellDataTemperatureWithDensityMassFractionsAndPressure(
    const hier::Box& domain)
{
    if (d_num_subghosts_temperature > -hier::IntVector::getOne(d_dim))
    {
        // Create the cell data of temperature.
        d_data_temperature.reset(
            new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_temperature));
        
        if (!d_data_density)
        {
            computeGlobalCellDataDensity(domain);
        }
        
        if (!d_data_mass_fractions)
        {
            computeGlobalCellDataMassFractionsWithDensity(domain);
        }
        
        if (!d_data_pressure)
        {
            computeGlobalCellDataPressureWithDensityMassFractionsAndInternalEnergy(domain);
        }
        
        // Compute the temperature field.
        d_equation_of_state_mixing_rules->computeTemperature(
            d_data_temperature,
            d_data_density,
            d_data_pressure,
            d_data_mass_fractions,
            domain);
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::"
            << "computeGlobalCellDataTemperatureWithDensityMassFractionsAndPressure()\n"
            << "Cell data of 'TEMPERATURE' is not yet registered."
            << std::endl);
    }
}


/*
 * Compute the global cell data of convective flux with velocity and pressure in the registered
 * patch.
 */
void
FlowModelFourEqnConservative::computeGlobalCellDataConvectiveFluxWithVelocityAndPressure(
    const DIRECTION::TYPE& direction,
    const hier::Box& domain)
{
    if (direction == DIRECTION::X_DIRECTION)
    {
        if (d_num_subghosts_convective_flux_x > -hier::IntVector::getOne(d_dim))
        {
            // Create the cell data of convective flux in the x-direction.
            d_data_convective_flux_x.reset(
                new pdat::CellData<double>(d_interior_box, d_num_eqn, d_num_subghosts_convective_flux_x));
            
            /*
             * Get the local lower indices and number of cells in each direction of the domain.
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
            
            boost::shared_ptr<pdat::CellData<double> > data_partial_densities =
                getGlobalCellDataPartialDensities();
            
            boost::shared_ptr<pdat::CellData<double> > data_momentum =
                getGlobalCellDataMomentum();
            
            boost::shared_ptr<pdat::CellData<double> > data_total_energy =
                getGlobalCellDataTotalEnergy();
            
            if (!d_data_velocity)
            {
                computeGlobalCellDataVelocityWithDensity(domain);
            }
            
            if (!d_data_pressure)
            {
                computeGlobalCellDataPressureWithDensityMassFractionsAndInternalEnergy(domain);
            }
            
            // Get the pointers to the cell data of partial densities, total energy and pressure.
            std::vector<double*> rho_Y;
            rho_Y.reserve(d_num_species);
            for (int si = 0; si < d_num_species; si++)
            {
                rho_Y.push_back(data_partial_densities->getPointer(si));
            }
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
                for (int si = 0; si < d_num_species; si++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx = i + num_ghosts_0;
                        const int idx_velocity = i + num_subghosts_0_velocity;
                        const int idx_convective_flux_x = i + num_subghosts_0_convective_flux_x;
                        
                        F_x[si][idx_convective_flux_x] = u[idx_velocity]*rho_Y[si][idx];
                    }
                }
                
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
                    
                    F_x[d_num_species][idx_convective_flux_x] = u[idx_velocity]*rho_u[idx] + p[idx_pressure];
                    F_x[d_num_species + 1][idx_convective_flux_x] = u[idx_velocity]*(E[idx] + p[idx_pressure]);
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
                for (int si = 0; si < d_num_species; si++)
                {
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            const int idx = (i + num_ghosts_0) +
                                (j + num_ghosts_1)*ghostcell_dim_0;
                            
                            const int idx_velocity = (i + num_subghosts_0_velocity) +
                                (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                            
                            const int idx_convective_flux_x = (i + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                            
                            F_x[si][idx_convective_flux_x] = u[idx_velocity]*rho_Y[si][idx];
                        }
                    }
                }
                
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
                        
                        F_x[d_num_species][idx_convective_flux_x] = u[idx_velocity]*rho_u[idx] + p[idx_pressure];
                        F_x[d_num_species + 1][idx_convective_flux_x] = u[idx_velocity]*rho_v[idx];
                        F_x[d_num_species + 2][idx_convective_flux_x] = u[idx_velocity]*(E[idx] + p[idx_pressure]);
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
                // Get the pointers to the cell data of momentum.
                double* rho_u = data_momentum->getPointer(0);
                double* rho_v = data_momentum->getPointer(1);
                double* rho_w = data_momentum->getPointer(2);
                
                // Get the pointer to the cell data of velocity.
                double* u = d_data_velocity->getPointer(0);
                
                // Compute the convective flux in the x-direction.
                for (int si = 0; si < d_num_species; si++)
                {
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
                                
                                const int idx_convective_flux_x = (i + num_subghosts_0_convective_flux_x) +
                                    (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                    (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                        subghostcell_dim_1_convective_flux_x;
                                
                                F_x[si][idx_convective_flux_x] = u[idx_velocity]*rho_Y[si][idx];
                            }
                        }
                    }
                }
                
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
                            
                            F_x[d_num_species][idx_convective_flux_x] = u[idx_velocity]*rho_u[idx] + p[idx_pressure];
                            F_x[d_num_species + 1][idx_convective_flux_x] = u[idx_velocity]*rho_v[idx];
                            F_x[d_num_species + 2][idx_convective_flux_x] = u[idx_velocity]*rho_w[idx];
                            F_x[d_num_species + 3][idx_convective_flux_x] = u[idx_velocity]*(E[idx] + p[idx_pressure]);
                        }
                    }
                }
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::"
                << "computeGlobalCellDataConvectiveFluxWithVelocityAndPressure()\n"
                << "Cell data of 'CONVECTIVE_FLUX_X' is not yet registered."
                << std::endl);
        }
    }
    else if (direction == DIRECTION::Y_DIRECTION)
    {
        if (d_num_subghosts_convective_flux_y > -hier::IntVector::getOne(d_dim))
        {
            // Create the cell data of convective flux in the y-direction.
            d_data_convective_flux_y.reset(
                new pdat::CellData<double>(d_interior_box, d_num_eqn, d_num_subghosts_convective_flux_y));
            
            /*
             * Get the local lower indices and number of cells in each direction of the domain.
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
            
            boost::shared_ptr<pdat::CellData<double> > data_partial_densities =
                getGlobalCellDataPartialDensities();
            
            boost::shared_ptr<pdat::CellData<double> > data_momentum =
                getGlobalCellDataMomentum();
            
            boost::shared_ptr<pdat::CellData<double> > data_total_energy =
                getGlobalCellDataTotalEnergy();
            
            if (!d_data_pressure)
            {
                computeGlobalCellDataPressureWithDensityMassFractionsAndInternalEnergy(domain);
            }
            
            if (!d_data_velocity)
            {
                computeGlobalCellDataVelocityWithDensity(domain);
            }
            
            // Get the pointers to the cell data of partial densities, total energy, volume fraction
            // and pressure.
            std::vector<double*> rho_Y;
            rho_Y.reserve(d_num_species);
            for (int si = 0; si < d_num_species; si++)
            {
                rho_Y.push_back(data_partial_densities->getPointer(si));
            }
            double* E = data_total_energy->getPointer(0);
            double* p = d_data_pressure->getPointer(0);
            
            if (d_dim == tbox::Dimension(1))
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelFourEqnConservative::"
                    << "computeGlobalCellDataConvectiveFluxWithVelocityAndPressure()\n"
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
                for (int si = 0; si < d_num_species; si++)
                {
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            const int idx = (i + num_ghosts_0) +
                                (j + num_ghosts_1)*ghostcell_dim_0;
                            
                            const int idx_velocity = (i + num_subghosts_0_velocity) +
                                (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                            
                            const int idx_convective_flux_y = (i + num_subghosts_0_convective_flux_y) +
                                (j + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                            
                            F_y[si][idx_convective_flux_y] = v[idx_velocity]*rho_Y[si][idx];
                        }
                    }
                }
                
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
                        
                        F_y[d_num_species][idx_convective_flux_y] = v[idx_velocity]*rho_u[idx];
                        F_y[d_num_species + 1][idx_convective_flux_y] = v[idx_velocity]*rho_v[idx] + p[idx_pressure];
                        F_y[d_num_species + 2][idx_convective_flux_y] = v[idx_velocity]*(E[idx] + p[idx_pressure]);
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
                for (int si = 0; si < d_num_species; si++)
                {
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
                                
                                const int idx_convective_flux_y = (i + num_subghosts_0_convective_flux_y) +
                                    (j + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                    (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                        subghostcell_dim_1_convective_flux_y;
                                
                                F_y[si][idx_convective_flux_y] = v[idx_velocity]*rho_Y[si][idx];
                            }
                        }
                    }
                }
                
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
                            
                            F_y[d_num_species][idx_convective_flux_y] = v[idx_velocity]*rho_u[idx];
                            F_y[d_num_species + 1][idx_convective_flux_y] = v[idx_velocity]*rho_v[idx] + p[idx_pressure];
                            F_y[d_num_species + 2][idx_convective_flux_y] = v[idx_velocity]*rho_w[idx];
                            F_y[d_num_species + 3][idx_convective_flux_y] = v[idx_velocity]*(E[idx] + p[idx_pressure]);
                        }
                    }
                }
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::"
                << "computeGlobalCellDataConvectiveFluxWithVelocityAndPressure()\n"
                << "Cell data of 'CONVECTIVE_FLUX_Y' is not yet registered."
                << std::endl);
        }
    }
    else if (direction == DIRECTION::Z_DIRECTION)
    {
        if (d_num_subghosts_convective_flux_z > -hier::IntVector::getOne(d_dim))
        {
            // Create the cell data of convective flux in the z-direction.
            d_data_convective_flux_z.reset(
                new pdat::CellData<double>(d_interior_box, d_num_eqn, d_num_subghosts_convective_flux_z));
            
            /*
             * Get the local lower indices and number of cells in each direction of the domain.
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
            
            boost::shared_ptr<pdat::CellData<double> > data_partial_densities =
                getGlobalCellDataPartialDensities();
            
            boost::shared_ptr<pdat::CellData<double> > data_momentum =
                getGlobalCellDataMomentum();
            
            boost::shared_ptr<pdat::CellData<double> > data_total_energy =
                getGlobalCellDataTotalEnergy();
            
            if (!d_data_pressure)
            {
                computeGlobalCellDataPressureWithDensityMassFractionsAndInternalEnergy(domain);
            }
            
            if (!d_data_velocity)
            {
                computeGlobalCellDataVelocityWithDensity(domain);
            }
            
            // Get the pointers to the cell data of partial densities, total energy, volume fraction
            // and pressure.
            std::vector<double*> rho_Y;
            rho_Y.reserve(d_num_species);
            for (int si = 0; si < d_num_species; si++)
            {
                rho_Y.push_back(data_partial_densities->getPointer(si));
            }
            double* E = data_total_energy->getPointer(0);
            double* p = d_data_pressure->getPointer(0);
            
            if (d_dim == tbox::Dimension(1) || d_dim == tbox::Dimension(2))
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelFourEqnConservative::"
                    << "computeGlobalCellDataConvectiveFluxWithVelocityAndPressure()\n"
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
                for (int si = 0; si < d_num_species; si++)
                {
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
                                
                                const int idx_convective_flux_z = (i + num_subghosts_0_convective_flux_z) +
                                    (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                    (k + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                        subghostcell_dim_1_convective_flux_z;
                                
                                F_z[si][idx_convective_flux_z] = w[idx_velocity]*rho_Y[si][idx];
                            }
                        }
                    }
                }
                
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
                            
                            F_z[d_num_species][idx_convective_flux_z] = w[idx_velocity]*rho_u[idx];
                            F_z[d_num_species + 1][idx_convective_flux_z] = w[idx_velocity]*rho_v[idx];
                            F_z[d_num_species + 2][idx_convective_flux_z] = w[idx_velocity]*rho_w[idx] + p[idx_pressure];
                            F_z[d_num_species + 3][idx_convective_flux_z] = w[idx_velocity]*(E[idx] + p[idx_pressure]);
                        }
                    }
                }
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::"
                << "computeGlobalCellDataConvectiveFluxWithVelocityAndPressure()\n"
                << "Cell data of 'CONVECTIVE_FLUX_Z' is not yet registered."
                << std::endl);
        }
    }
}


/*
 * Compute the global cell data of maximum wave speed with velocity and sound speed in the
 * registered patch.
 */
void
FlowModelFourEqnConservative::computeGlobalCellDataMaxWaveSpeedWithVelocityAndSoundSpeed(
    const DIRECTION::TYPE& direction,
    const hier::Box& domain)
{
    if (direction == DIRECTION::X_DIRECTION)
    {
        if (d_num_subghosts_max_wave_speed_x > -hier::IntVector::getOne(d_dim))
        {
            // Create the cell data of maximum wave speed in the x-direction.
            d_data_max_wave_speed_x.reset(
                new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_max_wave_speed_x));
            
            /*
             * Get the local lower indices and number of cells in each direction of the domain.
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
            
            if (!d_data_velocity)
            {
                computeGlobalCellDataVelocityWithDensity(domain);
            }
            
            if (!d_data_sound_speed)
            {
                computeGlobalCellDataSoundSpeedWithDensityMassFractionsAndPressure(domain);
            }
            
            // Get the pointers to the cell data of maximum wave speed and velocity in x-direction and sound speed.
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
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::"
                << "computeGlobalCellDataMaxWaveSpeedWithVelocityAndSoundSpeed()\n"
                << "Cell data of 'MAX_WAVE_SPEED_X' is not yet registered."
                << std::endl);
        }
    }
    else if (direction == DIRECTION::Y_DIRECTION)
    {
        if (d_num_subghosts_max_wave_speed_y > -hier::IntVector::getOne(d_dim))
        {
            // Create the cell data of maximum wave speed in the y-direction.
            d_data_max_wave_speed_y.reset(
                new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_max_wave_speed_y));
            
            /*
             * Get the local lower indices and number of cells in each direction of the domain.
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
            
            if (!d_data_sound_speed)
            {
                computeGlobalCellDataSoundSpeedWithDensityMassFractionsAndPressure(domain);
            }
            
            if (!d_data_velocity)
            {
                computeGlobalCellDataVelocityWithDensity(domain);
            }
            
            // Get the pointers to the cell data of maximum wave speed and velocity in y-direction, and sound speed.
            double* lambda_max_y = d_data_max_wave_speed_y->getPointer(0);
            double* v            = d_data_velocity->getPointer(1);
            double* c            = d_data_sound_speed->getPointer(0);
            
            if (d_dim == tbox::Dimension(1))
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelFourEqnConservative::"
                    << "computeGlobalCellDataMaxWaveSpeedWithVelocityAndSoundSpeed()\n"
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
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::"
                << "computeGlobalCellDataMaxWaveSpeedWithVelocityAndSoundSpeed()\n"
                << "Cell data of 'MAX_WAVE_SPEED_Y' is not yet registered."
                << std::endl);
        }
    }
    else if (direction == DIRECTION::Z_DIRECTION)
    {
        if (d_num_subghosts_max_wave_speed_z > -hier::IntVector::getOne(d_dim))
        {
            // Create the cell data of maximum wave speed in the z-direction.
            d_data_max_wave_speed_z.reset(
                new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_max_wave_speed_z));
            
            /*
             * Get the local lower indices and number of cells in each direction of the domain.
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
            
            if (!d_data_sound_speed)
            {
                computeGlobalCellDataSoundSpeedWithDensityMassFractionsAndPressure(domain);
            }
            
            if (!d_data_velocity)
            {
                computeGlobalCellDataVelocityWithDensity(domain);
            }
            
            // Get the pointers to the cell data of maximum wave speed and velocity in z-direction, and sound speed.
            double* lambda_max_z = d_data_max_wave_speed_z->getPointer(0);
            double* w            = d_data_velocity->getPointer(2);
            double* c            = d_data_sound_speed->getPointer(0);
            
            if (d_dim == tbox::Dimension(1) || d_dim == tbox::Dimension(2))
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelFourEqnConservative::"
                    << "computeGlobalCellDataMaxWaveSpeedWithVelocityAndSoundSpeed()\n"
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
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::"
                << "computeGlobalCellDataMaxWaveSpeedWithVelocityAndSoundSpeed()\n"
                << "Cell data of 'MAX_WAVE_SPEED_Z' is not yet registered."
                << std::endl);
        }
    }
}


/*
 * Compute the global cell data of maximum diffusivity with density, mass fractions, pressure
 * and temperature in the registered patch.
 */
void
FlowModelFourEqnConservative::computeGlobalCellDataMaxDiffusivityWithDensityMassFractionsPressureAndTemperature(
    const hier::Box& domain)
{
    if (!d_equation_of_mass_diffusivity_mixing_rules ||
        !d_equation_of_shear_viscosity_mixing_rules ||
        !d_equation_of_bulk_viscosity_mixing_rules ||
        !d_equation_of_thermal_conductivity_mixing_rules)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::"
            << "computeGlobalCellDataMaxDiffusivityWithDensityMassFractionsPressureAndTemperature()\n"
            << "Either mixing rule of mass diffusivity, shear viscosity, bulk viscosity or"
            << " thermal conductivity is not initialized."
            << std::endl);
    }
    
    if (d_num_subghosts_max_diffusivity > -hier::IntVector::getOne(d_dim))
    {
        // Create the cell data of maximum diffusivity.
        d_data_max_diffusivity.reset(
            new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_max_diffusivity));
        
        /*
         * Get the local lower indices and number of cells in each direction of the domain.
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
        
        if (!d_data_density)
        {
            computeGlobalCellDataDensity(domain);
        }
        
        if (!d_data_mass_fractions)
        {
            computeGlobalCellDataMassFractionsWithDensity(domain);
        }
        
        if (!d_data_pressure)
        {
            computeGlobalCellDataPressureWithDensityMassFractionsAndInternalEnergy(domain);
        }
        
        if (!d_data_temperature)
        {
            computeGlobalCellDataTemperatureWithDensityMassFractionsAndPressure(domain);
        }
        
        /*
         * Create temporary cell data of isobaric specific heat capacity, mass diffusivities, shear
         * viscosity, bulk viscosity and thermal conductivity.
         */
        
        boost::shared_ptr<pdat::CellData<double> > data_isobaric_specific_heat_capacity(
            new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_max_diffusivity));
        
        boost::shared_ptr<pdat::CellData<double> > data_mass_diffusivities(
            new pdat::CellData<double>(d_interior_box, d_num_species, d_num_subghosts_max_diffusivity));
        
        boost::shared_ptr<pdat::CellData<double> > data_shear_viscosity(
            new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_max_diffusivity));
        
        boost::shared_ptr<pdat::CellData<double> > data_bulk_viscosity(
            new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_max_diffusivity));
        
        boost::shared_ptr<pdat::CellData<double> > data_thermal_conductivity(
            new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_max_diffusivity));
        
        // Get the pointers to the cell data of maximum diffusivity, density, isobaric specific heat
        // capacity, mass diffusivities, shear viscosity, bulk viscosity and thermal conductivity.
        double* D_max = d_data_max_diffusivity->getPointer(0);
        double* rho   = d_data_density->getPointer(0);
        double* c_p   = data_isobaric_specific_heat_capacity->getPointer(0);
        std::vector<double*> D;
        D.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            D.push_back(data_mass_diffusivities->getPointer(si));
        }
        double* mu    = data_shear_viscosity->getPointer(0);
        double* mu_v  = data_bulk_viscosity->getPointer(0);
        double* kappa = data_thermal_conductivity->getPointer(0);
        
        // Compute the isobaric specific heat capacity field.
        d_equation_of_state_mixing_rules->computeIsobaricSpecificHeatCapacity(
            data_isobaric_specific_heat_capacity,
            d_data_density,
            d_data_pressure,
            d_data_mass_fractions,
            domain);
        
        // Compute the mass diffusivity fields.
        d_equation_of_mass_diffusivity_mixing_rules->computeMassDiffusivities(
            data_mass_diffusivities,
            d_data_pressure,
            d_data_temperature,
            d_data_mass_fractions,
            domain);
        
        // Compute the shear viscosity field.
        d_equation_of_shear_viscosity_mixing_rules->computeShearViscosity(
            data_shear_viscosity,
            d_data_pressure,
            d_data_temperature,
            d_data_mass_fractions,
            domain);
        
        // Compute the bulk viscosity field.
        d_equation_of_bulk_viscosity_mixing_rules->computeBulkViscosity(
            data_bulk_viscosity,
            d_data_pressure,
            d_data_temperature,
            d_data_mass_fractions,
            domain);
        
        // Compute the thermal conductivity field.
        d_equation_of_thermal_conductivity_mixing_rules->computeThermalConductivity(
            data_thermal_conductivity,
            d_data_pressure,
            d_data_temperature,
            d_data_mass_fractions,
            domain);
        
        if (d_dim == tbox::Dimension(1))
        {
            /*
             * Get the local lower index, numbers of cells in each dimension and numbers of ghost cells.
             */
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_dim_0 = domain_dims[0];
            
            const int num_subghosts_0_max_diffusivity = d_num_subghosts_max_diffusivity[0];
            const int num_subghosts_0_density = d_num_subghosts_density[0];
            
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_max_diffusivity = i + num_subghosts_0_max_diffusivity;
                const int idx_density = i + num_subghosts_0_density;
                
                D_max[idx_max_diffusivity] = fmax(mu[idx_max_diffusivity]/rho[idx_density],
                    mu_v[idx_max_diffusivity]/rho[idx_density]);
                
                D_max[idx_max_diffusivity] = fmax(D_max[idx_max_diffusivity],
                    kappa[idx_max_diffusivity]/(rho[idx_density]*c_p[idx_max_diffusivity]));
            }
            
            for (int si = 0; si < d_num_species; si++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_max_diffusivity = i + num_subghosts_0_max_diffusivity;
                    
                    D_max[idx_max_diffusivity] = fmax(D_max[idx_max_diffusivity],
                        D[si][idx_max_diffusivity]);
                }
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
            
            const int num_subghosts_0_max_diffusivity = d_num_subghosts_max_diffusivity[0];
            const int num_subghosts_1_max_diffusivity = d_num_subghosts_max_diffusivity[1];
            const int subghostcell_dim_0_max_diffusivity = d_subghostcell_dims_max_diffusivity[0];
            
            const int num_subghosts_0_density = d_num_subghosts_density[0];
            const int num_subghosts_1_density = d_num_subghosts_density[1];
            const int subghostcell_dim_0_density = d_subghostcell_dims_density[0];
            
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_max_diffusivity = (i + num_subghosts_0_max_diffusivity) +
                        (j + num_subghosts_1_max_diffusivity)*subghostcell_dim_0_max_diffusivity;
                    
                    const int idx_density = (i + num_subghosts_0_density) +
                        (j + num_subghosts_1_density)*subghostcell_dim_0_density;
                    
                    D_max[idx_max_diffusivity] = fmax(mu[idx_max_diffusivity]/rho[idx_density],
                        mu_v[idx_max_diffusivity]/rho[idx_density]);
                    
                    D_max[idx_max_diffusivity] = fmax(D_max[idx_max_diffusivity],
                        kappa[idx_max_diffusivity]/(rho[idx_density]*c_p[idx_max_diffusivity]));
                }
            }
            
            for (int si = 0; si < d_num_species; si++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx_max_diffusivity = (i + num_subghosts_0_max_diffusivity) +
                            (j + num_subghosts_1_max_diffusivity)*subghostcell_dim_0_max_diffusivity;
                        
                        D_max[idx_max_diffusivity] = fmax(D_max[idx_max_diffusivity],
                            D[si][idx_max_diffusivity]);
                    }
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
            
            const int num_subghosts_0_max_diffusivity = d_num_subghosts_max_diffusivity[0];
            const int num_subghosts_1_max_diffusivity = d_num_subghosts_max_diffusivity[1];
            const int num_subghosts_2_max_diffusivity = d_num_subghosts_max_diffusivity[2];
            const int subghostcell_dim_0_max_diffusivity = d_subghostcell_dims_max_diffusivity[0];
            const int subghostcell_dim_1_max_diffusivity = d_subghostcell_dims_max_diffusivity[1];
            
            const int num_subghosts_0_density = d_num_subghosts_density[0];
            const int num_subghosts_1_density = d_num_subghosts_density[1];
            const int num_subghosts_2_density = d_num_subghosts_density[2];
            const int subghostcell_dim_0_density = d_subghostcell_dims_density[0];
            const int subghostcell_dim_1_density = d_subghostcell_dims_density[1];
            
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
                        const int idx_max_diffusivity = (i + num_subghosts_0_max_diffusivity) +
                            (j + num_subghosts_1_max_diffusivity)*subghostcell_dim_0_max_diffusivity +
                            (k + num_subghosts_2_max_diffusivity)*subghostcell_dim_0_max_diffusivity*
                                subghostcell_dim_1_max_diffusivity;
                        
                        const int idx_density = (i + num_subghosts_0_density) +
                            (j + num_subghosts_1_density)*subghostcell_dim_0_density +
                            (k + num_subghosts_2_density)*subghostcell_dim_0_density*
                                subghostcell_dim_1_density;
                        
                        D_max[idx_max_diffusivity] = fmax(mu[idx_max_diffusivity]/rho[idx_density],
                            mu_v[idx_max_diffusivity]/rho[idx_density]);
                        
                        D_max[idx_max_diffusivity] = fmax(D_max[idx_max_diffusivity],
                            kappa[idx_max_diffusivity]/(rho[idx_density]*c_p[idx_max_diffusivity]));
                    }
                }
            }
            
            for (int si = 0; si < d_num_species; si++)
            {
                for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                {
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute the linear index.
                            const int idx_max_diffusivity = (i + num_subghosts_0_max_diffusivity) +
                                (j + num_subghosts_1_max_diffusivity)*subghostcell_dim_0_max_diffusivity +
                                (k + num_subghosts_2_max_diffusivity)*subghostcell_dim_0_max_diffusivity*
                                    subghostcell_dim_1_max_diffusivity;
                            
                            D_max[idx_max_diffusivity] = fmax(D_max[idx_max_diffusivity],
                                D[si][idx_max_diffusivity]);
                        }
                    }
                }            
            }
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::"
            << "computeGlobalCellDataMaxDiffusivityWithDensityMassFractionsPressureAndTemperature()\n"
            << "Cell data of 'MAX_DIFFUSIVITY' is not yet registered."
            << std::endl);
    }
}
