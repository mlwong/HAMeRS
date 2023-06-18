#include "flow/flow_models/four-eqn_conservative/FlowModelFourEqnConservative.hpp"

#include "flow/flow_models/four-eqn_conservative/FlowModelBasicUtilitiesFourEqnConservative.hpp"
#include "flow/flow_models/four-eqn_conservative/FlowModelBoundaryUtilitiesFourEqnConservative.hpp"
#include "flow/flow_models/four-eqn_conservative/FlowModelDiffusiveFluxUtilitiesFourEqnConservative.hpp"
#include "flow/flow_models/four-eqn_conservative/FlowModelImmersedBoundaryMethodFourEqnConservative.hpp"
#include "flow/flow_models/four-eqn_conservative/FlowModelMonitoringStatisticsUtilitiesFourEqnConservative.hpp"
#include "flow/flow_models/four-eqn_conservative/FlowModelRiemannSolverFourEqnConservative.hpp"
#include "flow/flow_models/four-eqn_conservative/FlowModelSourceUtilitiesFourEqnConservative.hpp"
#include "flow/flow_models/four-eqn_conservative/FlowModelStatisticsUtilitiesFourEqnConservative.hpp"

HAMERS_SHARED_PTR<pdat::CellVariable<Real> > FlowModelFourEqnConservative::s_variable_partial_densities;
HAMERS_SHARED_PTR<pdat::CellVariable<Real> > FlowModelFourEqnConservative::s_variable_momentum;
HAMERS_SHARED_PTR<pdat::CellVariable<Real> > FlowModelFourEqnConservative::s_variable_total_energy;

FlowModelFourEqnConservative::FlowModelFourEqnConservative(
    const std::string& object_name,
    const std::string& project_name,
    const tbox::Dimension& dim,
    const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_species,
    const HAMERS_SHARED_PTR<tbox::Database>& flow_model_db):
        FlowModel(
            object_name,
            project_name,
            dim,
            grid_geometry,
            num_species,
            num_species + dim.getValue() + 1,
            flow_model_db),
        d_num_subghosts_density(-hier::IntVector::getOne(d_dim)),
        d_num_subghosts_mass_fractions(-hier::IntVector::getOne(d_dim)),
        d_num_subghosts_mole_fractions(-hier::IntVector::getOne(d_dim)),
        d_num_subghosts_volume_fractions(-hier::IntVector::getOne(d_dim)),
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
        d_num_subghosts_species_densities(-hier::IntVector::getOne(d_dim)),
        d_num_subghosts_species_enthalpies(-hier::IntVector::getOne(d_dim)),
        d_subghost_box_density(hier::Box::getEmptyBox(d_dim)),
        d_subghost_box_mass_fractions(hier::Box::getEmptyBox(d_dim)),
        d_subghost_box_mole_fractions(hier::Box::getEmptyBox(d_dim)),
        d_subghost_box_volume_fractions(hier::Box::getEmptyBox(d_dim)),
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
        d_subghost_box_species_densities(hier::Box::getEmptyBox(d_dim)),
        d_subghost_box_species_enthalpies(hier::Box::getEmptyBox(d_dim)),
        d_subghostcell_dims_density(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_mass_fractions(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_mole_fractions(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_volume_fractions(hier::IntVector::getZero(d_dim)),
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
        d_subghostcell_dims_species_densities(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_species_enthalpies(hier::IntVector::getZero(d_dim)),
        d_cell_data_computed_density(false),
        d_cell_data_computed_mass_fractions(false),
        d_cell_data_computed_mole_fractions(false),
        d_cell_data_computed_volume_fractions(false),
        d_cell_data_computed_velocity(false),
        d_cell_data_computed_internal_energy(false),
        d_cell_data_computed_pressure(false),
        d_cell_data_computed_sound_speed(false),
        d_cell_data_computed_temperature(false),
        d_cell_data_computed_convective_flux_x(false),
        d_cell_data_computed_convective_flux_y(false),
        d_cell_data_computed_convective_flux_z(false),
        d_cell_data_computed_max_wave_speed_x(false),
        d_cell_data_computed_max_wave_speed_y(false),
        d_cell_data_computed_max_wave_speed_z(false),
        d_cell_data_computed_max_diffusivity(false),
        d_cell_data_computed_species_densities(false),
        d_cell_data_computed_species_enthalpies(false)
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
    
    /*
     * Initialize the conservative variables.
     */
    
    s_variable_partial_densities = HAMERS_SHARED_PTR<pdat::CellVariable<Real> > (
        new pdat::CellVariable<Real>(d_dim, "partial densities", d_num_species));
    
    s_variable_momentum = HAMERS_SHARED_PTR<pdat::CellVariable<Real> > (
        new pdat::CellVariable<Real>(d_dim, "momentum", d_dim.getValue()));
    
    s_variable_total_energy = HAMERS_SHARED_PTR<pdat::CellVariable<Real> > (
        new pdat::CellVariable<Real>(d_dim, "total energy", 1));
    
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
    
    HAMERS_SHARED_PTR<tbox::Database> equation_of_state_mixing_rules_db;
    
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
        
        HAMERS_SHARED_PTR<tbox::Database> equation_of_mass_diffusivity_mixing_rules_db;
        
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
        
        HAMERS_SHARED_PTR<tbox::Database> equation_of_shear_viscosity_mixing_rules_db;
        
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
        
        HAMERS_SHARED_PTR<tbox::Database> equation_of_bulk_viscosity_mixing_rules_db;
        
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
        
        HAMERS_SHARED_PTR<tbox::Database> equation_of_thermal_conductivity_mixing_rules_db;
        
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
     * Initialize basic utilities object.
     */
    d_flow_model_basic_utilities.reset(new FlowModelBasicUtilitiesFourEqnConservative(
        "d_flow_model_basic_utilities",
        d_dim,
        d_grid_geometry,
        d_num_species,
        d_equation_of_state_mixing_rules));
    
    /*
     * Initialize diffusive flux utilities object.
     */
    d_flow_model_diffusive_flux_utilities.reset(new FlowModelDiffusiveFluxUtilitiesFourEqnConservative(
        "d_flow_model_diffusive_flux_utilities",
        d_dim,
        d_grid_geometry,
        d_num_species,
        flow_model_db,
        d_equation_of_mass_diffusivity_mixing_rules,
        d_equation_of_shear_viscosity_mixing_rules,
        d_equation_of_bulk_viscosity_mixing_rules,
        d_equation_of_thermal_conductivity_mixing_rules));
    
    /*
     * Initialize source utilities object.
     */
    d_flow_model_source_utilities.reset(new FlowModelSourceUtilitiesFourEqnConservative(
        "d_flow_model_source_utilities",
        d_project_name,
        d_dim,
        d_grid_geometry,
        d_num_species,
        flow_model_db,
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
    
    /*
     * Initialize monitoring statistics utilities object.
     */
    d_flow_model_monitoring_statistics_utilities.reset(new FlowModelMonitoringStatisticsUtilitiesFourEqnConservative(
        "d_flow_model_monitoring_statistics_utilities",
        d_dim,
        d_grid_geometry,
        d_num_species,
        flow_model_db));
    
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
     * Initialize pointers to species cell data.
     */
    d_data_species_densities.resize(d_num_species, nullptr);
    d_data_species_enthalpies.resize(d_num_species, nullptr);
}


/*
 * Initialize the immersed boundary method object.
 */
void
FlowModelFourEqnConservative::initializeImmersedBoundaryMethod(
    const HAMERS_SHARED_PTR<ImmersedBoundaries>& immersed_boundaries,
    const HAMERS_SHARED_PTR<tbox::Database>& immersed_boundary_method_db)
{
    /*
     * Initialize immersed boundary method object.
     */
    d_flow_model_immersed_boundary_method.reset(
        new FlowModelImmersedBoundaryMethodFourEqnConservative(
            "d_flow_model_immersed_boundary_method",
            d_dim,
            d_grid_geometry,
            d_num_species,
            d_num_eqn,
            immersed_boundaries,
            immersed_boundary_method_db,
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
    const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const
{
    putToRestartBase(restart_db);
    
    /*
     * Put the properties of d_equation_of_state_mixing_rules into the restart database.
     */
    
    restart_db->putString("d_equation_of_state_str", d_equation_of_state_str);
    
    HAMERS_SHARED_PTR<tbox::Database> restart_equation_of_state_mixing_rules_db =
        restart_db->putDatabase("d_equation_of_state_mixing_rules_db");
    d_equation_of_state_mixing_rules->putToRestart(restart_equation_of_state_mixing_rules_db);
    
    /*
     * Put the properties of d_equation_of_mass_diffusivity_mixing_rules into the restart database.
     */
    
    if (d_equation_of_mass_diffusivity_mixing_rules)
    {
        restart_db->putString("d_equation_of_mass_diffusivity_str", d_equation_of_mass_diffusivity_str);
        
        HAMERS_SHARED_PTR<tbox::Database> restart_equation_of_mass_diffusivity_mixing_rules_db =
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
        
        HAMERS_SHARED_PTR<tbox::Database> restart_equation_of_shear_viscosity_mixing_rules_db =
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
        
        HAMERS_SHARED_PTR<tbox::Database> restart_equation_of_bulk_viscosity_mixing_rules_db =
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
        
        HAMERS_SHARED_PTR<tbox::Database> restart_equation_of_thermal_conductivity_mixing_rules_db =
            restart_db->putDatabase("d_equation_of_thermal_conductivity_mixing_rules_db");
        d_equation_of_thermal_conductivity_mixing_rules->
            putToRestart(restart_equation_of_thermal_conductivity_mixing_rules_db);
    }
    
    /*
     * Put the properties of d_flow_model_diffusive_flux_utilities into the restart database.
     */
    d_flow_model_diffusive_flux_utilities->putToRestart(restart_db);
    
    /*
     * Put the properties of d_flow_model_source_utilities into the restart database.
     */
    d_flow_model_source_utilities->putToRestart(restart_db);
    
    /*
     * Put the properties of d_flow_model_monitoring_statistics_utilities into the restart database.
     */
    d_flow_model_monitoring_statistics_utilities->putToRestart(restart_db);
    
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
std::vector<HAMERS_SHARED_PTR<pdat::CellVariable<Real> > >
FlowModelFourEqnConservative::getConservativeVariables()
{
    std::vector<HAMERS_SHARED_PTR<pdat::CellVariable<Real> > > conservative_variables;
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
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
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
    
    HAMERS_SHARED_PTR<pdat::CellData<Real> > data_partial_densities(
        HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
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
 * are given as entries in a map of the variable name to the number of sub-ghost cells required.
 * If the variable to be registered is one of the conservative variable, the corresponding entry
 * in the map is ignored.
 */
void
FlowModelFourEqnConservative::registerDerivedVariables(
    const std::unordered_map<std::string, hier::IntVector>& num_subghosts_of_data)
{
    // Check whether a patch is already registered.
    if (!d_patch)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::registerDerivedVariables()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    // Check whether all or part of derived cell data is already computed.
    if (d_derived_cell_data_computed)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::registerDerivedVariables()\n"
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
                << ": FlowModelFourEqnConservative::registerDerivedVariables()\n"
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
    
    if (num_subghosts_of_data.find("MOLE_FRACTIONS") != num_subghosts_of_data.end())
    {
        setNumberOfSubGhosts(
            num_subghosts_of_data.find("MOLE_FRACTIONS")->second,
            "MOLE_FRACTIONS",
            "MOLE_FRACTIONS");
    }
    
    if (num_subghosts_of_data.find("VOLUME_FRACTIONS") != num_subghosts_of_data.end())
    {
        setNumberOfSubGhosts(
            num_subghosts_of_data.find("VOLUME_FRACTIONS")->second,
            "VOLUME_FRACTIONS",
            "VOLUME_FRACTIONS");
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
    
    if (num_subghosts_of_data.find("SPECIES_DENSITIES") != num_subghosts_of_data.end())
    {
        setNumberOfSubGhosts(
            num_subghosts_of_data.find("SPECIES_DENSITIES")->second,
            "SPECIES_DENSITIES",
            "SPECIES_DENSITIES");
    }
    
    if (num_subghosts_of_data.find("SPECIES_ENTHALPIES") != num_subghosts_of_data.end())
    {
        setNumberOfSubGhosts(
            num_subghosts_of_data.find("SPECIES_ENTHALPIES")->second,
            "SPECIES_ENTHALPIES",
            "SPECIES_ENTHALPIES");
    }
}


/*
 * Unregister the registered patch. The registered data context and all global derived
 * cell data in the patch are cleared.
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
    
    d_num_ghosts                       = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_density            = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_mass_fractions     = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_mole_fractions     = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_volume_fractions   = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_velocity           = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_internal_energy    = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_pressure           = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_sound_speed        = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_temperature        = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_convective_flux_x  = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_convective_flux_y  = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_convective_flux_z  = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_max_wave_speed_x   = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_max_wave_speed_y   = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_max_wave_speed_z   = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_max_diffusivity    = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_species_densities  = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_species_enthalpies = -hier::IntVector::getOne(d_dim);
    
    d_interior_box                    = hier::Box::getEmptyBox(d_dim);
    d_ghost_box                       = hier::Box::getEmptyBox(d_dim);
    d_subdomain_box                   = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_density            = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_mass_fractions     = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_mole_fractions     = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_volume_fractions   = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_velocity           = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_internal_energy    = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_pressure           = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_sound_speed        = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_temperature        = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_convective_flux_x  = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_convective_flux_y  = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_convective_flux_z  = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_max_wave_speed_x   = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_max_wave_speed_y   = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_max_wave_speed_z   = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_max_diffusivity    = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_species_densities  = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_species_enthalpies = hier::Box::getEmptyBox(d_dim);
    
    d_interior_dims                        = hier::IntVector::getZero(d_dim);
    d_ghostcell_dims                       = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_density            = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_mass_fractions     = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_mole_fractions     = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_volume_fractions   = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_velocity           = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_internal_energy    = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_pressure           = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_sound_speed        = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_temperature        = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_convective_flux_x  = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_convective_flux_y  = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_convective_flux_z  = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_max_wave_speed_x   = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_max_wave_speed_y   = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_max_wave_speed_z   = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_max_diffusivity    = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_species_densities  = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_species_enthalpies = hier::IntVector::getZero(d_dim);
    
    d_data_density.reset();
    d_data_mass_fractions.reset();
    d_data_mole_fractions.reset();
    d_data_volume_fractions.reset();
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
    d_data_species_densities.assign(d_num_species, nullptr);
    d_data_species_enthalpies.assign(d_num_species, nullptr);
    
    d_cell_data_computed_density            = false;
    d_cell_data_computed_mass_fractions     = false;
    d_cell_data_computed_mole_fractions     = false;
    d_cell_data_computed_volume_fractions   = false;
    d_cell_data_computed_velocity           = false;
    d_cell_data_computed_internal_energy    = false;
    d_cell_data_computed_pressure           = false;
    d_cell_data_computed_sound_speed        = false;
    d_cell_data_computed_temperature        = false;
    d_cell_data_computed_convective_flux_x  = false;
    d_cell_data_computed_convective_flux_y  = false;
    d_cell_data_computed_convective_flux_z  = false;
    d_cell_data_computed_max_wave_speed_x   = false;
    d_cell_data_computed_max_wave_speed_y   = false;
    d_cell_data_computed_max_wave_speed_z   = false;
    d_cell_data_computed_max_diffusivity    = false;
    d_cell_data_computed_species_densities  = false;
    d_cell_data_computed_species_enthalpies = false;
    
    d_flow_model_diffusive_flux_utilities->clearCellAndSideData();
    d_flow_model_source_utilities->clearCellData();
    
    d_derived_cell_data_computed = false;
    
    d_patch = nullptr;
    clearDataContext();
}


/*
 * Allocate memory for cell data of different registered derived variables.
 */
void
FlowModelFourEqnConservative::allocateMemoryForDerivedCellData()
{
    if (d_num_subghosts_density > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_density)
        {
            if (!d_data_density)
            {
                // Create the cell data of density.
                d_data_density.reset(
                    new pdat::CellData<Real>(d_interior_box, 1, d_num_subghosts_density));
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::allocateMemoryForDerivedCellData()\n"
                << "Cell data of 'DENSITY' is aleady computed."
                << std::endl);
        }
    }
    
    if (d_num_subghosts_mass_fractions > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_mass_fractions)
        {
            if (!d_data_mass_fractions)
            {
                // Create the cell data of mass fractions.
                d_data_mass_fractions.reset(
                    new pdat::CellData<Real>(d_interior_box, d_num_species, d_num_subghosts_mass_fractions));
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::allocateMemoryForDerivedCellData()\n"
                << "Cell data of 'MASS_FRACTIONS' is aleady computed."
                << std::endl);
        }
    }
    
    if (d_num_subghosts_mole_fractions > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_mole_fractions)
        {
            if (!d_data_mole_fractions)
            {
                // Create the cell data of mole fractions.
                d_data_mole_fractions.reset(
                    new pdat::CellData<Real>(d_interior_box, d_num_species, d_num_subghosts_mole_fractions));
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::allocateMemoryForDerivedCellData()\n"
                << "Cell data of 'MOLE_FRACTIONS' is aleady computed."
                << std::endl);
        }
    }
    
    if (d_num_subghosts_volume_fractions > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_volume_fractions)
        {
            if (!d_data_volume_fractions)
            {
                // Create the cell data of volume fractions.
                d_data_volume_fractions.reset(
                    new pdat::CellData<Real>(d_interior_box, d_num_species, d_num_subghosts_volume_fractions));
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::allocateMemoryForDerivedCellData()\n"
                << "Cell data of 'VOLUME_FRACTIONS' is aleady computed."
                << std::endl);
        }
    }
    
    if (d_num_subghosts_velocity > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_velocity)
        {
            if (!d_data_velocity)
            {
                // Create the cell data of velocity.
                d_data_velocity.reset(
                    new pdat::CellData<Real>(d_interior_box, d_dim.getValue(), d_num_subghosts_velocity));
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::allocateMemoryForDerivedCellData()\n"
                << "Cell data of 'VELOCITY' is aleady computed."
                << std::endl);
        }
    }
    
    if (d_num_subghosts_internal_energy > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_internal_energy)
        {
            if (!d_data_internal_energy)
            {
                // Create the cell data of internal energy.
                d_data_internal_energy.reset(
                    new pdat::CellData<Real>(d_interior_box, 1, d_num_subghosts_internal_energy));
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::allocateMemoryForDerivedCellData()\n"
                << "Cell data of 'INTERNAL_ENERGY' is aleady computed."
                << std::endl);
        }
    }
    
    if (d_num_subghosts_pressure > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_pressure)
        {
            if (!d_data_pressure)
            {
                // Create the cell data of pressure.
                d_data_pressure.reset(
                    new pdat::CellData<Real>(d_interior_box, 1, d_num_subghosts_pressure));
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::allocateMemoryForDerivedCellData()\n"
                << "Cell data of 'PRESSURE' is aleady computed."
                << std::endl);
        }
    }
    
    if (d_num_subghosts_sound_speed > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_sound_speed)
        {
            if (!d_data_sound_speed)
            {
                // Create the cell data of sound speed.
                d_data_sound_speed.reset(
                    new pdat::CellData<Real>(d_interior_box, 1, d_num_subghosts_sound_speed));
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::allocateMemoryForDerivedCellData()\n"
                << "Cell data of 'SOUND_SPEED' is aleady computed."
                << std::endl);
        }
    }
    
    if (d_num_subghosts_temperature > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_temperature)
        {
            if (!d_data_temperature)
            {
                // Create the cell data of temperature.
                d_data_temperature.reset(
                    new pdat::CellData<Real>(d_interior_box, 1, d_num_subghosts_temperature));
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::allocateMemoryForDerivedCellData()\n"
                << "Cell data of 'TEMPERATURE' is aleady computed."
                << std::endl);
        }
    }
    
    if (d_num_subghosts_convective_flux_x > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_convective_flux_x)
        {
            if (!d_data_convective_flux_x)
            {
                // Create the cell data of convective flux in the x-direction.
                d_data_convective_flux_x.reset(
                    new pdat::CellData<Real>(d_interior_box, d_num_eqn, d_num_subghosts_convective_flux_x));
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::allocateMemoryForDerivedCellData()\n"
                << "Cell data of 'CONVECTIVE_FLUX_X' is aleady computed."
                << std::endl);
        }
    }
    
    if (d_num_subghosts_convective_flux_y > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_convective_flux_y)
        {
            if (!d_data_convective_flux_y)
            {
                // Create the cell data of convective flux in the y-direction.
                d_data_convective_flux_y.reset(
                    new pdat::CellData<Real>(d_interior_box, d_num_eqn, d_num_subghosts_convective_flux_y));
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::allocateMemoryForDerivedCellData()\n"
                << "Cell data of 'CONVECTIVE_FLUX_Y' is aleady computed."
                << std::endl);
        }
    }
    
    if (d_num_subghosts_convective_flux_z > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_convective_flux_z)
        {
            if (!d_data_convective_flux_z)
            {
                // Create the cell data of convective flux in the z-direction.
                d_data_convective_flux_z.reset(
                    new pdat::CellData<Real>(d_interior_box, d_num_eqn, d_num_subghosts_convective_flux_z));
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::allocateMemoryForDerivedCellData()\n"
                << "Cell data of 'CONVECTIVE_FLUX_Z' is aleady computed."
                << std::endl);
        }
    }
    
    if (d_num_subghosts_max_wave_speed_x > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_max_wave_speed_x)
        {
            if (!d_data_max_wave_speed_x)
            {
                // Create the cell data of maximum wave speed in the x-direction.
                d_data_max_wave_speed_x.reset(
                    new pdat::CellData<Real>(d_interior_box, 1, d_num_subghosts_max_wave_speed_x));
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::allocateMemoryForDerivedCellData()\n"
                << "Cell data of 'MAX_WAVE_SPEED_X' is aleady computed."
                << std::endl);
        }
    }
    
    if (d_num_subghosts_max_wave_speed_y > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_max_wave_speed_y)
        {
            if (!d_data_max_wave_speed_y)
            {
                // Create the cell data of maximum wave speed in the y-direction.
                d_data_max_wave_speed_y.reset(
                    new pdat::CellData<Real>(d_interior_box, 1, d_num_subghosts_max_wave_speed_y));
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::allocateMemoryForDerivedCellData()\n"
                << "Cell data of 'MAX_WAVE_SPEED_Y' is aleady computed."
                << std::endl);
        }
    }
    
    if (d_num_subghosts_max_wave_speed_z > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_max_wave_speed_z)
        {
            if (!d_data_max_wave_speed_z)
            {
                // Create the cell data of maximum wave speed in the z-direction.
                d_data_max_wave_speed_z.reset(
                    new pdat::CellData<Real>(d_interior_box, 1, d_num_subghosts_max_wave_speed_z));
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::allocateMemoryForDerivedCellData()\n"
                << "Cell data of 'MAX_WAVE_SPEED_Z' is aleady computed."
                << std::endl);
        }
    }
    
    if (d_num_subghosts_max_diffusivity > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_max_diffusivity)
        {
            if (!d_data_max_diffusivity)
            {
                // Create the cell data of maximum diffusivity.
                d_data_max_diffusivity.reset(
                    new pdat::CellData<Real>(d_interior_box, 1, d_num_subghosts_max_diffusivity));
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::allocateMemoryForDerivedCellData()\n"
                << "Cell data of 'MAX_DIFFUSIVITY' is aleady computed."
                << std::endl);
        }
    }
    
    if (d_num_subghosts_species_densities > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_species_densities)
        {
            // Create the cell data of species densities.
            for (int si = 0; si < d_num_species; si++)
            {
                if (!d_data_species_densities[si])
                {
                    d_data_species_densities[si].reset(
                        new pdat::CellData<Real>(d_interior_box, 1, d_num_subghosts_species_densities));
                }
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::allocateMemoryForDerivedCellData()\n"
                << "Cell data of 'SPECIES_DENSITIES' is aleady computed."
                << std::endl);
        }
    }
    
    if (d_num_subghosts_species_enthalpies > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_species_enthalpies)
        {
            // Create the cell data of species enthalpies.
            for (int si = 0; si < d_num_species; si++)
            {
                d_data_species_enthalpies[si].reset(
                    new pdat::CellData<Real>(d_interior_box, 1, d_num_subghosts_species_enthalpies));
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::allocateMemoryForDerivedCellData()\n"
                << "Cell data of 'SPECIES_ENTHALPIES' is aleady computed."
                << std::endl);
        }
    }
}


/*
 * Compute the cell data of different registered derived variables with the registered data context.
 */
void
FlowModelFourEqnConservative::computeDerivedCellData()
{
    // Check whether a patch is already registered.
    if (!d_patch)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::computeDerivedCellData()\n"
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
    
    // Compute the total density cell data.
    if (d_num_subghosts_density > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_density)
        {
            computeCellDataOfDensity(
                d_subdomain_box);
        }
    }
    
    // Compute the mass fraction cell data.
    if (d_num_subghosts_mass_fractions > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_mass_fractions)
        {
            computeCellDataOfMassFractionsWithDensity(
                d_subdomain_box);
        }
    }
    
    // Compute the mole fraction cell data.
    if (d_num_subghosts_mole_fractions > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_mole_fractions)
        {
            computeCellDataOfMoleFractionsWithMassFractions(
                d_subdomain_box);
        }
    }
    
    // Compute the volume fraction cell data.
    if (d_num_subghosts_volume_fractions > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_volume_fractions)
        {
            computeCellDataOfVolumeFractionsWithSpeciesDensities(
                d_subdomain_box);
        }
    }
    
    // Compute the velocity cell data.
    if (d_num_subghosts_velocity > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_velocity)
        {
            computeCellDataOfVelocityWithDensity(
                d_subdomain_box);
        }
    }
    
    // Compute the internal energy cell data.
    if (d_num_subghosts_internal_energy > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_internal_energy)
        {
            computeCellDataOfInternalEnergyWithDensityAndVelocity(
                d_subdomain_box);
        }
    }
    
    // Compute the pressure cell data.
    if (d_num_subghosts_pressure > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_pressure)
        {
            computeCellDataOfPressureWithDensityMassFractionsAndInternalEnergy(
                d_subdomain_box);
        }
    }
    
    // Compute the sound speed cell data.
    if (d_num_subghosts_sound_speed > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_sound_speed)
        {
            computeCellDataOfSoundSpeedWithDensityMassFractionsAndPressure(
                d_subdomain_box);
        }
    }
    
    // Compute the temperature cell data.
    if (d_num_subghosts_temperature > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_temperature)
        {
            computeCellDataOfTemperatureWithDensityMassFractionsAndPressure(
                d_subdomain_box);
        }
    }
    
    // Compute the x-direction convective flux cell data.
    if (d_num_subghosts_convective_flux_x > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_convective_flux_x)
        {
            computeCellDataOfConvectiveFluxWithVelocityAndPressure(
                DIRECTION::X_DIRECTION,
                d_subdomain_box);
        }
    }
    
    // Compute the y-direction convective flux cell data.
    if (d_num_subghosts_convective_flux_y > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_convective_flux_y)
        {
            computeCellDataOfConvectiveFluxWithVelocityAndPressure(
                DIRECTION::Y_DIRECTION,
                d_subdomain_box);
        }
    }
    
    // Compute the z-direction convective flux cell data.
    if (d_num_subghosts_convective_flux_z > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_convective_flux_z)
        {
            computeCellDataOfConvectiveFluxWithVelocityAndPressure(
                DIRECTION::Z_DIRECTION,
                d_subdomain_box);
        }
    }
    
    // Compute the x-direction maximum wave speed cell data.
    if (d_num_subghosts_max_wave_speed_x > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_max_wave_speed_x)
        {
            computeCellDataOfMaxWaveSpeedWithVelocityAndSoundSpeed(
                DIRECTION::X_DIRECTION,
                d_subdomain_box);
        }
    }
    
    // Compute the y-direction maximum wave speed cell data.
    if (d_num_subghosts_max_wave_speed_y > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_max_wave_speed_y)
        {
            computeCellDataOfMaxWaveSpeedWithVelocityAndSoundSpeed(
                DIRECTION::Y_DIRECTION,
                d_subdomain_box);
        }
    }
    
    // Compute the z-direction maximum wave speed cell data.
    if (d_num_subghosts_max_wave_speed_z > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_max_wave_speed_z)
        {
            computeCellDataOfMaxWaveSpeedWithVelocityAndSoundSpeed(
                DIRECTION::Z_DIRECTION,
                d_subdomain_box);
        }
    }
    
    // Compute the maximum diffusivity cell data.
    if (d_num_subghosts_max_diffusivity > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_max_diffusivity)
        {
            computeCellDataOfMaxDiffusivityWithDensityMassFractionsPressureAndTemperature(
                d_subdomain_box);
        }
    }
    
    // Compute the species densities cell data.
    if (d_num_subghosts_species_densities > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_species_densities)
        {
            computeCellDataOfSpeciesDensitiesWithPressureAndTemperature(
                d_subdomain_box);
        }
    }
    
    // Compute the species enthalpies cell data.
    if (d_num_subghosts_species_enthalpies > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_species_enthalpies)
        {
            computeCellDataOfSpeciesEnthalpiesWithSpeciesDensitiesAndPressure(
                d_subdomain_box);
        }
    }
    
    d_derived_cell_data_computed = true;
}


/*
 * Get the cell data of one cell variable in the registered patch.
 */
HAMERS_SHARED_PTR<pdat::CellData<Real> >
FlowModelFourEqnConservative::getCellData(const std::string& variable_key)
{
    // Check whether a patch is already registered.
    if (!d_patch)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::getCellData()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    HAMERS_SHARED_PTR<pdat::CellData<Real> > cell_data;
    
    if (variable_key == "PARTIAL_DENSITIES")
    {
        cell_data = getCellDataOfPartialDensities();
    }
    else if (variable_key == "MOMENTUM")
    {
        cell_data = getCellDataOfMomentum();
    }
    else if (variable_key == "TOTAL_ENERGY")
    {
        cell_data = getCellDataOfTotalEnergy();
    }
    else if (variable_key == "DENSITY")
    {
        if (!d_cell_data_computed_density)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::getCellData()\n"
                << "Cell data of 'DENSITY' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_density;
    }
    else if (variable_key == "MASS_FRACTIONS")
    {
        if (!d_cell_data_computed_mass_fractions)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::getCellData()\n"
                << "Cell data of 'MASS_FRACTIONS' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_mass_fractions;
    }
    else if (variable_key == "MOLE_FRACTIONS")
    {
        if (!d_cell_data_computed_mole_fractions)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::getCellData()\n"
                << "Cell data of 'MOLE_FRACTIONS' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_mole_fractions;
    }
    else if (variable_key == "VOLUME_FRACTIONS")
    {
        if (!d_cell_data_computed_volume_fractions)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::getCellData()\n"
                << "Cell data of 'VOLUME_FRACTIONS' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_volume_fractions;
    }
    else if (variable_key == "VELOCITY")
    {
        if (!d_cell_data_computed_velocity)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::getCellData()\n"
                << "Cell data of 'VELOCITY' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_velocity;
    }
    else if (variable_key == "INTERNAL_ENERGY")
    {
        if (!d_cell_data_computed_internal_energy)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::getCellData()\n"
                << "Cell data of 'INTERNAL_ENERGY' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_internal_energy;
    }
    else if (variable_key == "PRESSURE")
    {
        if (!d_cell_data_computed_pressure)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::getCellData()\n"
                << "Cell data of 'PRESSURE' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_pressure;
    }
    else if (variable_key == "SOUND_SPEED")
    {
        if (!d_cell_data_computed_sound_speed)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::getCellData()\n"
                << "Cell data of 'SOUND_SPEED' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_sound_speed;
    }
    else if (variable_key == "TEMPERATURE")
    {
        if (!d_cell_data_computed_temperature)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::getCellData()\n"
                << "Cell data of 'TEMPERATURE' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_temperature;
    }
    else if (variable_key == "CONVECTIVE_FLUX_X")
    {
        if (!d_cell_data_computed_convective_flux_x)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::getCellData()\n"
                << "Cell data of 'CONVECTIVE_FLUX_X' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_convective_flux_x;
    }
    else if (variable_key == "CONVECTIVE_FLUX_Y")
    {
        if (!d_cell_data_computed_convective_flux_y)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::getCellData()\n"
                << "Cell data of 'CONVECTIVE_FLUX_Y' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_convective_flux_y;
    }
    else if (variable_key == "CONVECTIVE_FLUX_Z")
    {
        if (!d_cell_data_computed_convective_flux_z)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::getCellData()\n"
                << "Cell data of 'CONVECTIVE_FLUX_Z' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_convective_flux_z;
    }
    else if (variable_key == "MAX_WAVE_SPEED_X")
    {
        if (!d_cell_data_computed_max_wave_speed_x)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::getCellData()\n"
                << "Cell data of 'MAX_WAVE_SPEED_X' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_max_wave_speed_x;
    }
    else if (variable_key == "MAX_WAVE_SPEED_Y")
    {
        if (!d_cell_data_computed_max_wave_speed_y)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::getCellData()\n"
                << "Cell data of 'MAX_WAVE_SPEED_Y' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_max_wave_speed_y;
    }
    else if (variable_key == "MAX_WAVE_SPEED_Z")
    {
        if (!d_cell_data_computed_max_wave_speed_z)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::getCellData()\n"
                << "Cell data of 'MAX_WAVE_SPEED_Z' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_max_wave_speed_z;
    }
    else if (variable_key == "MAX_DIFFUSIVITY")
    {
        if (!d_cell_data_computed_max_diffusivity)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::getCellData()\n"
                << "Cell data of 'MAX_DIFFUSIVITY' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_max_diffusivity;
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::getCellData()\n"
            << "Unknown cell data with variable_key = '" << variable_key
            << "' requested."
            << std::endl);
    }
    
    return cell_data;
}


/*
 * Get the cell data of different cell variables in the registered patch.
 */
std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >
FlowModelFourEqnConservative::getCellData(
    const std::vector<std::string>& variable_keys)
{
    std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > cell_data(
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
std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >
FlowModelFourEqnConservative::getSpeciesCellData(
    const std::string& variable_key)
{
    // Check whether the patch is already registered.
    if (!d_patch)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::getSpeciesCellData()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > species_cell_data;
    
    if (variable_key == "SPECIES_DENSITIES")
    {
        if (!d_cell_data_computed_species_densities)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::getSpeciesCellData()\n"
                << "Cell data of 'SPECIES_DENSITIES' is not registered/computed yet."
                << std::endl);
        }
        species_cell_data = d_data_species_densities;
    }
    else if (variable_key == "SPECIES_ENTHALPIES")
    {
        if (!d_cell_data_computed_species_enthalpies)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::getSpeciesCellData()\n"
                << "Cell data of 'SPECIES_ENTHALPIES' is not registered/computed yet."
                << std::endl);
        }
        species_cell_data = d_data_species_enthalpies;
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::getSpeciesCellData()\n"
            << "Unknown cell data with variable_key = '" << variable_key
            << "' requested."
            << std::endl);
    }
    
    return species_cell_data;
}


/*
 * Fill the cell data of conservative variables in the interior box with value zero.
 * Only fill the data when the mask has valid value if a mask cell data is given.
 */
void
FlowModelFourEqnConservative::fillCellDataOfConservativeVariablesWithZero(
    const HAMERS_SHARED_PTR<pdat::CellData<int> >& mask_cell_data,
    const int mask_valid_value)
{
    // Check whether a patch is already registered.
    if (!d_patch)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::fillCellDataOfConservativeVariablesWithZero()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    /*
     * Check whether the mask is used. Get the pointer to the cell data of the mask if the mask is used.
     * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
     */
    
    bool use_mask = false;
    int* mask = nullptr;
    hier::IntVector num_ghosts_mask(d_dim);
    hier::IntVector ghostcell_dims_mask(d_dim);
    
    if (mask_cell_data != nullptr)
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(mask_cell_data->getGhostCellWidth() >= d_num_ghosts);
#endif
        
        use_mask = true;
        mask = mask_cell_data->getPointer(0);
        
        num_ghosts_mask = mask_cell_data->getGhostCellWidth();
        ghostcell_dims_mask = mask_cell_data->getGhostBox().numberCells();
    }
    
    HAMERS_SHARED_PTR<pdat::CellData<Real> > data_partial_densities = getCellDataOfPartialDensities();
    HAMERS_SHARED_PTR<pdat::CellData<Real> > data_momentum = getCellDataOfMomentum();
    HAMERS_SHARED_PTR<pdat::CellData<Real> > data_total_energy = getCellDataOfTotalEnergy();
    
    // Get the pointers to the conservative variables.
    std::vector<Real*> rho_Y;
    rho_Y.reserve(d_num_species);
    for (int si = 0; si < d_num_species; si++)
    {
        rho_Y.push_back(data_partial_densities->getPointer(si));
    }
    Real* E = data_total_energy->getPointer(0);
    
    if (d_dim == tbox::Dimension(1))
    {
        // Get the pointer to the momentum component.
        Real* rho_u = data_momentum->getPointer(0);
        
        // Get the dimension and the number of ghost cells.
        
        const int interior_dim_0 = d_interior_dims[0];
        
        const int num_ghosts_0 = d_num_ghosts[0];
        
        if (use_mask)
        {
            const int num_ghosts_0_mask = num_ghosts_mask[0];
            
            for (int si = 0; si < d_num_species; si++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx = i + num_ghosts_0;
                    const int idx_mask = i + num_ghosts_0_mask;
                    
                    if (mask[idx_mask] == mask_valid_value)
                    {
                        rho_Y[si][idx] = Real(0);
                    }
                }
            }
            
            HAMERS_PRAGMA_SIMD
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx = i + num_ghosts_0;
                const int idx_mask = i + num_ghosts_0_mask;
                
                if (mask[idx_mask] == mask_valid_value)
                {
                    rho_u[idx] = Real(0);
                    E[idx]     = Real(0);
                }
            }
        }
        else
        {
            for (int si = 0; si < d_num_species; si++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx = i + num_ghosts_0;
                    
                    rho_Y[si][idx] = Real(0);
                }
            }
            
            HAMERS_PRAGMA_SIMD
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear index.
                const int idx = i + num_ghosts_0;
                
                rho_u[idx] = Real(0);
                E[idx]     = Real(0);
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        // Get the pointers to the momentum components.
        Real* rho_u = data_momentum->getPointer(0);
        Real* rho_v = data_momentum->getPointer(1);
        
        // Get the dimensions and the numbers of ghost cells.
        
        const int interior_dim_0 = d_interior_dims[0];
        const int interior_dim_1 = d_interior_dims[1];
        
        const int num_ghosts_0 = d_num_ghosts[0];
        const int num_ghosts_1 = d_num_ghosts[1];
        const int ghostcell_dim_0 = d_ghostcell_dims[0];
        
        if (use_mask)
        {
            const int num_ghosts_0_mask = num_ghosts_mask[0];
            const int num_ghosts_1_mask = num_ghosts_mask[1];
            const int ghostcell_dim_0_mask = ghostcell_dims_mask[0];
            
            for (int si = 0; si < d_num_species; si++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx  = (i + num_ghosts_0) +
                            (j + num_ghosts_1)*ghostcell_dim_0;
                        
                        const int idx_mask = (i + num_ghosts_0_mask) +
                            (j + num_ghosts_1_mask)*ghostcell_dim_0_mask;
                        
                        if (mask[idx_mask] == mask_valid_value)
                        {
                            rho_Y[si][idx] = Real(0);
                        }
                    }
                }
            }
            
            for (int j = 0; j < interior_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx  = (i + num_ghosts_0) +
                        (j + num_ghosts_1)*ghostcell_dim_0;
                    
                    const int idx_mask = (i + num_ghosts_0_mask) +
                        (j + num_ghosts_1_mask)*ghostcell_dim_0_mask;
                    
                    if (mask[idx_mask] == mask_valid_value)
                    {
                        rho_u[idx] = Real(0);
                        rho_v[idx] = Real(0);
                        E[idx]     = Real(0);
                    }
                }
            }
        }
        else
        {
            for (int si = 0; si < d_num_species; si++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx  = (i + num_ghosts_0) +
                            (j + num_ghosts_1)*ghostcell_dim_0;
                        
                        rho_Y[si][idx] = Real(0);
                    }
                }
            }
            
            for (int j = 0; j < interior_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx  = (i + num_ghosts_0) +
                        (j + num_ghosts_1)*ghostcell_dim_0;
                    
                    rho_u[idx] = Real(0);
                    rho_v[idx] = Real(0);
                    E[idx]     = Real(0);
                }
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        // Get the pointers to the momentum components.
        Real* rho_u = data_momentum->getPointer(0);
        Real* rho_v = data_momentum->getPointer(1);
        Real* rho_w = data_momentum->getPointer(2);
        
        // Get the dimensions and the numbers of ghost cells.
        
        const int interior_dim_0 = d_interior_dims[0];
        const int interior_dim_1 = d_interior_dims[1];
        const int interior_dim_2 = d_interior_dims[2];
        
        const int num_ghosts_0 = d_num_ghosts[0];
        const int num_ghosts_1 = d_num_ghosts[1];
        const int num_ghosts_2 = d_num_ghosts[2];
        const int ghostcell_dim_0 = d_ghostcell_dims[0];
        const int ghostcell_dim_1 = d_ghostcell_dims[1];
        
        if (use_mask)
        {
            const int num_ghosts_0_mask = num_ghosts_mask[0];
            const int num_ghosts_1_mask = num_ghosts_mask[1];
            const int num_ghosts_2_mask = num_ghosts_mask[2];
            const int ghostcell_dim_0_mask = ghostcell_dims_mask[0];
            const int ghostcell_dim_1_mask = ghostcell_dims_mask[1];
            
            for (int si = 0; si < d_num_species; si++)
            {
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx = (i + num_ghosts_0) +
                                (j + num_ghosts_1)*ghostcell_dim_0 +
                                (k + num_ghosts_2)*ghostcell_dim_0*ghostcell_dim_1;
                            
                            const int idx_mask = (i + num_ghosts_0_mask) +
                                (j + num_ghosts_1_mask)*ghostcell_dim_0_mask +
                                (k + num_ghosts_2_mask)*ghostcell_dim_0_mask*ghostcell_dim_1_mask;
                            
                            if (mask[idx_mask] == mask_valid_value)
                            {
                                rho_Y[si][idx] = Real(0);
                            }
                        }
                    }
                }
            }
            
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx = (i + num_ghosts_0) +
                            (j + num_ghosts_1)*ghostcell_dim_0 +
                            (k + num_ghosts_2)*ghostcell_dim_0*ghostcell_dim_1;
                        
                        const int idx_mask = (i + num_ghosts_0_mask) +
                            (j + num_ghosts_1_mask)*ghostcell_dim_0_mask +
                            (k + num_ghosts_2_mask)*ghostcell_dim_0_mask*ghostcell_dim_1_mask;
                        
                        if (mask[idx_mask] == mask_valid_value)
                        {
                            rho_u[idx] = Real(0);
                            rho_v[idx] = Real(0);
                            rho_w[idx] = Real(0);
                            E[idx]     = Real(0);
                        }
                    }
                }
            }
        }
        else
        {
            for (int si = 0; si < d_num_species; si++)
            {
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear index.
                            const int idx = (i + num_ghosts_0) +
                                (j + num_ghosts_1)*ghostcell_dim_0 +
                                (k + num_ghosts_2)*ghostcell_dim_0*ghostcell_dim_1;
                            
                            rho_Y[si][idx] = Real(0);
                        }
                    }
                }
            }
            
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx = (i + num_ghosts_0) +
                            (j + num_ghosts_1)*ghostcell_dim_0 +
                            (k + num_ghosts_2)*ghostcell_dim_0*ghostcell_dim_1;
                        
                        rho_u[idx] = Real(0);
                        rho_v[idx] = Real(0);
                        rho_w[idx] = Real(0);
                        E[idx]     = Real(0);
                    }
                }
            }
        }
    }
}


/*
 * Update the cell data of conservative variables in the interior box after time advancement.
 * Only update the data when the mask has valid value if a mask cell data is given.
 */
void
FlowModelFourEqnConservative::updateCellDataOfConservativeVariables(
    const HAMERS_SHARED_PTR<pdat::CellData<int> >& mask_cell_data,
    const int mask_valid_value)
{
    NULL_USE(mask_cell_data);
    NULL_USE(mask_valid_value);
    
    // Check whether a patch is already registered.
    if (!d_patch)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::updateCellDataOfConservativeVariables()\n"
            << "No patch is registered yet."
            << std::endl);
    }
}


/*
 * Get the cell data of the conservative variables in the registered patch.
 */
std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >
FlowModelFourEqnConservative::getCellDataOfConservativeVariables()
{
    // Check whether a patch is already registered.
    if (!d_patch)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::getCellDataOfConservativeVariables()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > cell_data;
    cell_data.reserve(3);
    
    cell_data.push_back(getCellDataOfPartialDensities());
    cell_data.push_back(getCellDataOfMomentum());
    cell_data.push_back(getCellDataOfTotalEnergy());
    
    return cell_data;
}


/*
 * Get the cell data of the primitive variables in the registered patch.
 */
std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >
FlowModelFourEqnConservative::getCellDataOfPrimitiveVariables()
{
    // Check whether a patch is already registered.
    if (!d_patch)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::getCellDataOfPrimitiveVariables()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > cell_data;
    cell_data.reserve(3);
    
    cell_data.push_back(getCellDataOfPartialDensities());
    if (!d_data_velocity)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::getCellDataOfPrimitiveVariables()\n"
            << "Cell data of 'VELOCITY' is not registered/computed yet."
            << std::endl);
    }
    cell_data.push_back(d_data_velocity);
    if (!d_data_pressure)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::getCellDataOfPrimitiveVariables()\n"
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
        HAMERS_SHARED_PTR<pdat::CellData<Real> > data_partial_densities(
            HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
                patch.getPatchData(s_variable_partial_densities, d_plot_context)));
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_partial_densities);
        TBOX_ASSERT(data_partial_densities->getGhostBox().isSpatiallyEqual(patch.getBox()));
#endif
        
        // Get the dimensions of box that covers the data.
        const hier::Box data_box = data_partial_densities->getGhostBox();
        const hier::IntVector data_dims = data_box.numberCells();
        
        // Get the pointers to the conservative variables.
        std::vector<const Real*> rho_Y;
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
                
                std::vector<const Real*> rho_Y_ptr;
                rho_Y_ptr.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    rho_Y_ptr.push_back(&rho_Y[si][idx_data]);
                }
                
                buffer[idx_region] = double(d_equation_of_state_mixing_rules->getMixtureDensity(rho_Y_ptr));
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
                    
                    std::vector<const Real*> rho_Y_ptr;
                    rho_Y_ptr.reserve(d_num_species);
                    for (int si = 0; si < d_num_species; si++)
                    {
                        rho_Y_ptr.push_back(&rho_Y[si][idx_data]);
                    }
                    
                    buffer[idx_region] = double(d_equation_of_state_mixing_rules->getMixtureDensity(rho_Y_ptr));
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
                        
                        std::vector<const Real*> rho_Y_ptr;
                        rho_Y_ptr.reserve(d_num_species);
                        for (int si = 0; si < d_num_species; si++)
                        {
                            rho_Y_ptr.push_back(&rho_Y[si][idx_data]);
                        }
                        
                        buffer[idx_region] = double(d_equation_of_state_mixing_rules->getMixtureDensity(rho_Y_ptr));
                    }
                }
            }
        }
        
        data_on_patch = true;
    }
    else if (variable_name == "pressure")
    {
        HAMERS_SHARED_PTR<pdat::CellData<Real> > data_partial_densities(
            HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
                patch.getPatchData(s_variable_partial_densities, d_plot_context)));
        
        HAMERS_SHARED_PTR<pdat::CellData<Real> > data_momentum(
            HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
                patch.getPatchData(s_variable_momentum, d_plot_context)));
        
        HAMERS_SHARED_PTR<pdat::CellData<Real> > data_total_energy(
            HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
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
        std::vector<const Real*> rho_Y;
        rho_Y.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y.push_back(data_partial_densities->getPointer(si));
        }
        const Real* const rho_u = data_momentum->getPointer(0);
        const Real* const rho_v = d_dim > tbox::Dimension(1) ? data_momentum->getPointer(1) : NULL;
        const Real* const rho_w = d_dim > tbox::Dimension(2) ? data_momentum->getPointer(2) : NULL;
        const Real* const E     = data_total_energy->getPointer(0);
        
        size_t offset_data = data_box.offset(region.lower());
        
        if (d_dim == tbox::Dimension(1))
        {
            for (int i = 0; i < region_dims[0]; i++)
            {
                // Compute the linear indices.
                size_t idx_data = offset_data + i;
                size_t idx_region = i;
                
                // Compute the mixture density.
                std::vector<const Real*> rho_Y_ptr;
                rho_Y_ptr.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    rho_Y_ptr.push_back(&(rho_Y[si][idx_data]));
                }
                const Real rho = d_equation_of_state_mixing_rules->getMixtureDensity(
                    rho_Y_ptr);
                
                /*
                 * Compute the internal energy.
                 */
                Real epsilon = Real(0);
                if (d_dim == tbox::Dimension(1))
                {
                    epsilon = (E[idx_data] -
                        Real(1)/Real(2)*rho_u[idx_data]*rho_u[idx_data]/rho)/rho;
                }
                else if (d_dim == tbox::Dimension(2))
                {
                    epsilon = (E[idx_data] -
                        Real(1)/Real(2)*(rho_u[idx_data]*rho_u[idx_data] +
                        rho_v[idx_data]*rho_v[idx_data])/rho)/rho;
                }
                else if (d_dim == tbox::Dimension(3))
                {
                    epsilon = (E[idx_data] -
                        Real(1)/Real(2)*(rho_u[idx_data]*rho_u[idx_data] +
                        rho_v[idx_data]*rho_v[idx_data] +
                        rho_w[idx_data]*rho_w[idx_data])/rho)/rho;
                }
                
                /*
                 * Compute the mass fractions.
                 */
                Real Y[d_num_species];
                for (int si = 0; si < d_num_species; si++)
                {
                    Y[si] = rho_Y[si][idx_data]/rho;
                }
                
                /*
                 * Get the pointers to the mass fractions.
                 */
                std::vector<const Real*> Y_ptr;
                Y_ptr.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    Y_ptr.push_back(&Y[si]);
                }
                
                buffer[idx_region] = double(d_equation_of_state_mixing_rules->
                    getPressure(
                        &rho,
                        &epsilon,
                        Y_ptr));
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
                    std::vector<const Real*> rho_Y_ptr;
                    rho_Y_ptr.reserve(d_num_species);
                    for (int si = 0; si < d_num_species; si++)
                    {
                        rho_Y_ptr.push_back(&(rho_Y[si][idx_data]));
                    }
                    const Real rho = d_equation_of_state_mixing_rules->getMixtureDensity(
                        rho_Y_ptr);
                    
                    /*
                     * Compute the internal energy.
                     */
                    Real epsilon = Real(0);
                    if (d_dim == tbox::Dimension(1))
                    {
                        epsilon = (E[idx_data] -
                            Real(1)/Real(2)*rho_u[idx_data]*rho_u[idx_data]/rho)/rho;
                    }
                    else if (d_dim == tbox::Dimension(2))
                    {
                        epsilon = (E[idx_data] -
                            Real(1)/Real(2)*(rho_u[idx_data]*rho_u[idx_data] +
                            rho_v[idx_data]*rho_v[idx_data])/rho)/rho;
                    }
                    else if (d_dim == tbox::Dimension(3))
                    {
                        epsilon = (E[idx_data] -
                            Real(1)/Real(2)*(rho_u[idx_data]*rho_u[idx_data] +
                            rho_v[idx_data]*rho_v[idx_data] +
                            rho_w[idx_data]*rho_w[idx_data])/rho)/rho;
                    }
                    
                    /*
                     * Compute the mass fractions.
                     */
                    Real Y[d_num_species];
                    for (int si = 0; si < d_num_species; si++)
                    {
                        Y[si] = rho_Y[si][idx_data]/rho;
                    }
                    
                    /*
                     * Get the pointers to the mass fractions.
                     */
                    std::vector<const Real*> Y_ptr;
                    Y_ptr.reserve(d_num_species);
                    for (int si = 0; si < d_num_species; si++)
                    {
                        Y_ptr.push_back(&Y[si]);
                    }
                    
                    buffer[idx_region] = double(d_equation_of_state_mixing_rules->
                        getPressure(
                            &rho,
                            &epsilon,
                            Y_ptr));
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
                        std::vector<const Real*> rho_Y_ptr;
                        rho_Y_ptr.reserve(d_num_species);
                        for (int si = 0; si < d_num_species; si++)
                        {
                            rho_Y_ptr.push_back(&(rho_Y[si][idx_data]));
                        }
                        const Real rho = d_equation_of_state_mixing_rules->getMixtureDensity(
                            rho_Y_ptr);
                        
                        /*
                         * Compute the internal energy.
                         */
                        Real epsilon = Real(0);
                        if (d_dim == tbox::Dimension(1))
                        {
                            epsilon = (E[idx_data] -
                                Real(1)/Real(2)*rho_u[idx_data]*rho_u[idx_data]/rho)/rho;
                        }
                        else if (d_dim == tbox::Dimension(2))
                        {
                            epsilon = (E[idx_data] -
                                Real(1)/Real(2)*(rho_u[idx_data]*rho_u[idx_data] +
                                rho_v[idx_data]*rho_v[idx_data])/rho)/rho;
                        }
                        else if (d_dim == tbox::Dimension(3))
                        {
                            epsilon = (E[idx_data] -
                                Real(1)/Real(2)*(rho_u[idx_data]*rho_u[idx_data] +
                                rho_v[idx_data]*rho_v[idx_data] +
                                rho_w[idx_data]*rho_w[idx_data])/rho)/rho;
                        }
                        
                        /*
                         * Compute the mass fractions.
                         */
                        Real Y[d_num_species];
                        for (int si = 0; si < d_num_species; si++)
                        {
                            Y[si] = rho_Y[si][idx_data]/rho;
                        }
                        
                        /*
                         * Get the pointers to the mass fractions.
                         */
                        std::vector<const Real*> Y_ptr;
                        Y_ptr.reserve(d_num_species);
                        for (int si = 0; si < d_num_species; si++)
                        {
                            Y_ptr.push_back(&Y[si]);
                        }
                        
                        buffer[idx_region] = double(d_equation_of_state_mixing_rules->
                            getPressure(
                                &rho,
                                &epsilon,
                                Y_ptr));
                    }
                }
            }
        }
        
        data_on_patch = true;
    }
    else if (variable_name == "sound speed")
    {
        HAMERS_SHARED_PTR<pdat::CellData<Real> > data_partial_densities(
            HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
                patch.getPatchData(s_variable_partial_densities, d_plot_context)));
        
        HAMERS_SHARED_PTR<pdat::CellData<Real> > data_momentum(
            HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
                patch.getPatchData(s_variable_momentum, d_plot_context)));
        
        HAMERS_SHARED_PTR<pdat::CellData<Real> > data_total_energy(
            HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
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
        std::vector<const Real*> rho_Y;
        rho_Y.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y.push_back(data_partial_densities->getPointer(si));
        }
        const Real* const rho_u = data_momentum->getPointer(0);
        const Real* const rho_v = d_dim > tbox::Dimension(1) ? data_momentum->getPointer(1) : NULL;
        const Real* const rho_w = d_dim > tbox::Dimension(2) ? data_momentum->getPointer(2) : NULL;
        const Real* const E     = data_total_energy->getPointer(0);
        
        size_t offset_data = data_box.offset(region.lower());
        
        if (d_dim == tbox::Dimension(1))
        {
            for (int i = 0; i < region_dims[0]; i++)
            {
                // Compute the linear indices.
                size_t idx_data = offset_data + i;
                size_t idx_region = i;
                
                // Compute the mixture density.
                std::vector<const Real*> rho_Y_ptr;
                rho_Y_ptr.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    rho_Y_ptr.push_back(&(rho_Y[si][idx_data]));
                }
                const Real rho = d_equation_of_state_mixing_rules->getMixtureDensity(
                    rho_Y_ptr);
                
                /*
                 * Compute the internal energy.
                 */
                Real epsilon = Real(0);
                if (d_dim == tbox::Dimension(1))
                {
                    epsilon = (E[idx_data] -
                        Real(1)/Real(2)*rho_u[idx_data]*rho_u[idx_data]/rho)/rho;
                }
                else if (d_dim == tbox::Dimension(2))
                {
                    epsilon = (E[idx_data] -
                        Real(1)/Real(2)*(rho_u[idx_data]*rho_u[idx_data] +
                        rho_v[idx_data]*rho_v[idx_data])/rho)/rho;
                }
                else if (d_dim == tbox::Dimension(3))
                {
                    epsilon = (E[idx_data] -
                        Real(1)/Real(2)*(rho_u[idx_data]*rho_u[idx_data] +
                        rho_v[idx_data]*rho_v[idx_data] +
                        rho_w[idx_data]*rho_w[idx_data])/rho)/rho;
                }
                
                /*
                 * Compute the mass fractions.
                 */
                Real Y[d_num_species];
                for (int si = 0; si < d_num_species; si++)
                {
                    Y[si] = rho_Y[si][idx_data]/rho;
                }
                
                /*
                 * Get the pointers to the mass fractions.
                 */
                std::vector<const Real*> Y_ptr;
                Y_ptr.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    Y_ptr.push_back(&Y[si]);
                }
                
                const Real p = d_equation_of_state_mixing_rules->
                    getPressure(
                        &rho,
                        &epsilon,
                        Y_ptr);
                
                /*
                 * Compute the Gruneisen parameter.
                 */
                const Real Gamma = d_equation_of_state_mixing_rules->
                    getGruneisenParameter(
                        &rho,
                        &p,
                        Y_ptr);
                
                /*
                 * Compute the partial pressure partial partial densities.
                 */
                std::vector<Real> Psi = d_equation_of_state_mixing_rules->
                    getPressureDerivativeWithPartialDensities(
                        &rho,
                        &p,
                        Y_ptr);
                
                /*
                 * Compute the sound speed.
                 */
                buffer[idx_region] = double(Gamma*p/rho);
                
                for (int si = 0; si < d_num_species; si++)
                {
                    buffer[idx_region] += double(Y[si]*Psi[si]);
                }
                
                buffer[idx_region] = double(std::sqrt(buffer[idx_region]));
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
                    std::vector<const Real*> rho_Y_ptr;
                    rho_Y_ptr.reserve(d_num_species);
                    for (int si = 0; si < d_num_species; si++)
                    {
                        rho_Y_ptr.push_back(&(rho_Y[si][idx_data]));
                    }
                    const Real rho = d_equation_of_state_mixing_rules->getMixtureDensity(
                        rho_Y_ptr);
                    
                    /*
                     * Compute the internal energy.
                     */
                    Real epsilon = Real(0);
                    if (d_dim == tbox::Dimension(1))
                    {
                        epsilon = (E[idx_data] -
                            Real(1)/Real(2)*rho_u[idx_data]*rho_u[idx_data]/rho)/rho;
                    }
                    else if (d_dim == tbox::Dimension(2))
                    {
                        epsilon = (E[idx_data] -
                            Real(1)/Real(2)*(rho_u[idx_data]*rho_u[idx_data] +
                            rho_v[idx_data]*rho_v[idx_data])/rho)/rho;
                    }
                    else if (d_dim == tbox::Dimension(3))
                    {
                        epsilon = (E[idx_data] -
                            Real(1)/Real(2)*(rho_u[idx_data]*rho_u[idx_data] +
                            rho_v[idx_data]*rho_v[idx_data] +
                            rho_w[idx_data]*rho_w[idx_data])/rho)/rho;
                    }
                    
                    /*
                     * Compute the mass fractions.
                     */
                    Real Y[d_num_species];
                    for (int si = 0; si < d_num_species; si++)
                    {
                        Y[si] = rho_Y[si][idx_data]/rho;
                    }
                    
                    /*
                     * Get the pointers to the mass fractions.
                     */
                    std::vector<const Real*> Y_ptr;
                    Y_ptr.reserve(d_num_species);
                    for (int si = 0; si < d_num_species; si++)
                    {
                        Y_ptr.push_back(&Y[si]);
                    }
                    
                    const Real p = d_equation_of_state_mixing_rules->
                        getPressure(
                            &rho,
                            &epsilon,
                            Y_ptr);
                    
                    /*
                     * Compute the Gruneisen parameter.
                     */
                    const Real Gamma = d_equation_of_state_mixing_rules->
                        getGruneisenParameter(
                            &rho,
                            &p,
                            Y_ptr);
                    
                    /*
                     * Compute the partial pressure partial partial densities.
                     */
                    std::vector<Real> Psi = d_equation_of_state_mixing_rules->
                        getPressureDerivativeWithPartialDensities(
                            &rho,
                            &p,
                            Y_ptr);
                    
                    /*
                     * Compute the sound speed.
                     */
                    buffer[idx_region] = double(Gamma*p/rho);
                    
                    for (int si = 0; si < d_num_species; si++)
                    {
                        buffer[idx_region] += double(Y[si]*Psi[si]);
                    }
                    
                    buffer[idx_region] = double(std::sqrt(buffer[idx_region]));
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
                        std::vector<const Real*> rho_Y_ptr;
                        rho_Y_ptr.reserve(d_num_species);
                        for (int si = 0; si < d_num_species; si++)
                        {
                            rho_Y_ptr.push_back(&(rho_Y[si][idx_data]));
                        }
                        const Real rho = d_equation_of_state_mixing_rules->getMixtureDensity(
                            rho_Y_ptr);
                        
                        /*
                         * Compute the internal energy.
                         */
                        Real epsilon = Real(0);
                        if (d_dim == tbox::Dimension(1))
                        {
                            epsilon = (E[idx_data] -
                                Real(1)/Real(2)*rho_u[idx_data]*rho_u[idx_data]/rho)/rho;
                        }
                        else if (d_dim == tbox::Dimension(2))
                        {
                            epsilon = (E[idx_data] -
                                Real(1)/Real(2)*(rho_u[idx_data]*rho_u[idx_data] +
                                rho_v[idx_data]*rho_v[idx_data])/rho)/rho;
                        }
                        else if (d_dim == tbox::Dimension(3))
                        {
                            epsilon = (E[idx_data] -
                                Real(1)/Real(2)*(rho_u[idx_data]*rho_u[idx_data] +
                                rho_v[idx_data]*rho_v[idx_data] +
                                rho_w[idx_data]*rho_w[idx_data])/rho)/rho;
                        }
                        
                        /*
                         * Compute the mass fractions.
                         */
                        Real Y[d_num_species];
                        for (int si = 0; si < d_num_species; si++)
                        {
                            Y[si] = rho_Y[si][idx_data]/rho;
                        }
                        
                        /*
                         * Get the pointers to the mass fractions.
                         */
                        std::vector<const Real*> Y_ptr;
                        Y_ptr.reserve(d_num_species);
                        for (int si = 0; si < d_num_species; si++)
                        {
                            Y_ptr.push_back(&Y[si]);
                        }
                        
                        const Real p = d_equation_of_state_mixing_rules->
                            getPressure(
                                &rho,
                                &epsilon,
                                Y_ptr);
                        
                        /*
                         * Compute the Gruneisen parameter.
                         */
                        const Real Gamma = d_equation_of_state_mixing_rules->
                            getGruneisenParameter(
                                &rho,
                                &p,
                                Y_ptr);
                        
                        /*
                         * Compute the partial pressure partial partial densities.
                         */
                        std::vector<Real> Psi = d_equation_of_state_mixing_rules->
                            getPressureDerivativeWithPartialDensities(
                                &rho,
                                &p,
                                Y_ptr);
                        
                        /*
                         * Compute the sound speed.
                         */
                        buffer[idx_region] = double(Gamma*p/rho);
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            buffer[idx_region] += double(Y[si]*Psi[si]);
                        }
                        
                        buffer[idx_region] = double(std::sqrt(buffer[idx_region]));
                    }
                }
            }
        }
        
        data_on_patch = true;
    }
    else if (variable_name == "velocity")
    {
        HAMERS_SHARED_PTR<pdat::CellData<Real> > data_partial_densities(
            HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
                patch.getPatchData(s_variable_partial_densities, d_plot_context)));
        
        HAMERS_SHARED_PTR<pdat::CellData<Real> > data_momentum(
            HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
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
        std::vector<const Real*> rho_Y;
        rho_Y.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y.push_back(data_partial_densities->getPointer(si));
        }
        const Real* const m = data_momentum->getPointer(depth_id);
        
        size_t offset_data = data_box.offset(region.lower());
        
        if (d_dim == tbox::Dimension(1))
        {
            for (int i = 0; i < region_dims[0]; i++)
            {
                // Compute the linear indices.
                size_t idx_data = offset_data + i;
                size_t idx_region = i;
                
                std::vector<const Real*> rho_Y_ptr;
                rho_Y_ptr.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    rho_Y_ptr.push_back(&rho_Y[si][idx_data]);
                }
                
                Real rho = d_equation_of_state_mixing_rules->getMixtureDensity(rho_Y_ptr);
                
                buffer[idx_region] = double(m[idx_data]/rho);
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
                    
                    std::vector<const Real*> rho_Y_ptr;
                    rho_Y_ptr.reserve(d_num_species);
                    for (int si = 0; si < d_num_species; si++)
                    {
                        rho_Y_ptr.push_back(&rho_Y[si][idx_data]);
                    }
                    
                    Real rho = d_equation_of_state_mixing_rules->getMixtureDensity(rho_Y_ptr);
                    
                    buffer[idx_region] = double(m[idx_data]/rho);
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
                        
                        std::vector<const Real*> rho_Y_ptr;
                        rho_Y_ptr.reserve(d_num_species);
                        for (int si = 0; si < d_num_species; si++)
                        {
                            rho_Y_ptr.push_back(&rho_Y[si][idx_data]);
                        }
                        
                        Real rho = d_equation_of_state_mixing_rules->getMixtureDensity(rho_Y_ptr);
                        
                        buffer[idx_region] = double(m[idx_data]/rho);
                    }
                }
            }
        }
        
        data_on_patch = true;
    }
    else if (variable_name.find("mass fraction") != std::string::npos)
    {
        int species_idx = std::stoi(variable_name.substr(14));
        
        HAMERS_SHARED_PTR<pdat::CellData<Real> > data_partial_densities(
            HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
                patch.getPatchData(s_variable_partial_densities, d_plot_context)));
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_partial_densities);
        TBOX_ASSERT(data_partial_densities->getGhostBox().isSpatiallyEqual(patch.getBox()));
#endif
        
        // Get the dimensions of box that covers the data.
        const hier::Box data_box = data_partial_densities->getGhostBox();
        const hier::IntVector data_dims = data_box.numberCells();
        
        // Get the pointers to conservative variables.
        std::vector<const Real*> rho_Y;
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
                
                std::vector<const Real*> rho_Y_ptr;
                rho_Y_ptr.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    rho_Y_ptr.push_back(&rho_Y[si][idx_data]);
                }
                
                Real rho = d_equation_of_state_mixing_rules->getMixtureDensity(rho_Y_ptr);
                
                buffer[idx_region] = double(rho_Y[species_idx][idx_data]/rho);
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
                    
                    std::vector<const Real*> rho_Y_ptr;
                    rho_Y_ptr.reserve(d_num_species);
                    for (int si = 0; si < d_num_species; si++)
                    {
                        rho_Y_ptr.push_back(&rho_Y[si][idx_data]);
                    }
                    
                    Real rho = d_equation_of_state_mixing_rules->getMixtureDensity(rho_Y_ptr);
                    
                    buffer[idx_region] = double(rho_Y[species_idx][idx_data]/rho);
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
                        
                        std::vector<const Real*> rho_Y_ptr;
                        rho_Y_ptr.reserve(d_num_species);
                        for (int si = 0; si < d_num_species; si++)
                        {
                            rho_Y_ptr.push_back(&rho_Y[si][idx_data]);
                        }
                        
                        Real rho = d_equation_of_state_mixing_rules->getMixtureDensity(rho_Y_ptr);
                        
                        buffer[idx_region] = double(rho_Y[species_idx][idx_data]/rho);
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
    const HAMERS_SHARED_PTR<ExtendedVisItDataWriter>& visit_writer)
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
    else if (variable_name == "MOLE_FRACTIONS")
    {
        if (d_num_subghosts_mole_fractions > -hier::IntVector::getOne(d_dim))
        {
            if (num_subghosts > d_num_subghosts_mole_fractions)
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
                
                d_num_subghosts_mole_fractions = num_subghosts;
            }
        }
        else
        {
            d_num_subghosts_mole_fractions = num_subghosts;
        }
        
        setNumberOfSubGhosts(num_subghosts, "MASS_FRACTIONS", parent_variable_name);
    }
    else if (variable_name == "VOLUME_FRACTIONS")
    {
        if (d_num_subghosts_volume_fractions > -hier::IntVector::getOne(d_dim))
        {
            if (num_subghosts > d_num_subghosts_volume_fractions)
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
                
                d_num_subghosts_volume_fractions = num_subghosts;
            }
        }
        else
        {
            d_num_subghosts_volume_fractions = num_subghosts;
        }
        
        setNumberOfSubGhosts(num_subghosts, "SPECIES_DENSITIES", parent_variable_name);
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
    else if (variable_name == "SPECIES_DENSITIES")
    {
        if (d_num_subghosts_species_densities > -hier::IntVector::getOne(d_dim))
        {
            if (num_subghosts > d_num_subghosts_species_densities)
            {
                /*
                TBOX_ERROR(d_object_name
                    << ": FlowModelFiveEqnAllaire::setNumberOfSubGhosts()\n"
                    << "Number of ghosts of '"
                    << parent_variable_name
                    << "' exceeds"
                    << " number of ghosts of '"
                    << variable_name
                    << "'."
                    << std::endl);
                */
                
                d_num_subghosts_species_densities = num_subghosts;
            }
        }
        else
        {
            d_num_subghosts_species_densities = num_subghosts;
        }
        
        setNumberOfSubGhosts(num_subghosts, "PRESSURE", parent_variable_name);
        setNumberOfSubGhosts(num_subghosts, "TEMPERATURE", parent_variable_name);
    }
    else if (variable_name == "SPECIES_ENTHALPIES")
    {
        if (d_num_subghosts_species_enthalpies > -hier::IntVector::getOne(d_dim))
        {
            if (num_subghosts > d_num_subghosts_species_enthalpies)
            {
                /*
                TBOX_ERROR(d_object_name
                    << ": FlowModelFiveEqnAllaire::setNumberOfSubGhosts()\n"
                    << "Number of ghosts of '"
                    << parent_variable_name
                    << "' exceeds"
                    << " number of ghosts of '"
                    << variable_name
                    << "'."
                    << std::endl);
                */
                
                d_num_subghosts_species_enthalpies = num_subghosts;
            }
        }
        else
        {
            d_num_subghosts_species_enthalpies = num_subghosts;
        }
        
        setNumberOfSubGhosts(num_subghosts, "SPECIES_DENSITIES", parent_variable_name);
        setNumberOfSubGhosts(num_subghosts, "PRESSURE", parent_variable_name);
    }
}


/*
 * Set the ghost boxes of derived cell variables.
 */
void
FlowModelFourEqnConservative::setDerivedCellVariableGhostBoxes()
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
    
    if (d_num_subghosts_mole_fractions > -hier::IntVector::getOne(d_dim))
    {
        d_subghost_box_mole_fractions = d_interior_box;
        d_subghost_box_mole_fractions.grow(d_num_subghosts_mole_fractions);
        d_subghostcell_dims_mole_fractions = d_subghost_box_mole_fractions.numberCells();
    }
    
    if (d_num_subghosts_volume_fractions > -hier::IntVector::getOne(d_dim))
    {
        d_subghost_box_volume_fractions = d_interior_box;
        d_subghost_box_volume_fractions.grow(d_num_subghosts_volume_fractions);
        d_subghostcell_dims_volume_fractions = d_subghost_box_volume_fractions.numberCells();
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
    
    if (d_num_subghosts_species_densities > -hier::IntVector::getOne(d_dim))
    {
        d_subghost_box_species_densities = d_interior_box;
        d_subghost_box_species_densities.grow(d_num_subghosts_species_densities);
        d_subghostcell_dims_species_densities = d_subghost_box_species_densities.numberCells();
    }
    
    if (d_num_subghosts_species_enthalpies > -hier::IntVector::getOne(d_dim))
    {
        d_subghost_box_species_enthalpies = d_interior_box;
        d_subghost_box_species_enthalpies.grow(d_num_subghosts_species_enthalpies);
        d_subghostcell_dims_species_enthalpies = d_subghost_box_species_enthalpies.numberCells();
    }
}


/*
 * Get the cell data of partial densities in the registered patch.
 */
HAMERS_SHARED_PTR<pdat::CellData<Real> >
FlowModelFourEqnConservative::getCellDataOfPartialDensities()
{
    // Get the cell data of the registered variable partial densities.
    HAMERS_SHARED_PTR<pdat::CellData<Real> > data_partial_densities(
        HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
            d_patch->getPatchData(s_variable_partial_densities, getDataContext())));
    
    return data_partial_densities;
}


/*
 * Get the cell data of momentum in the registered patch.
 */
HAMERS_SHARED_PTR<pdat::CellData<Real> >
FlowModelFourEqnConservative::getCellDataOfMomentum()
{
    // Get the cell data of the registered variable momentum.
    HAMERS_SHARED_PTR<pdat::CellData<Real> > data_momentum(
        HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
            d_patch->getPatchData(s_variable_momentum, getDataContext())));
    
    return data_momentum;
}


/*
 * Get the cell data of total energy in the registered patch.
 */
HAMERS_SHARED_PTR<pdat::CellData<Real> >
FlowModelFourEqnConservative::getCellDataOfTotalEnergy()
{
    // Get the cell data of the registered variable total energy.
    HAMERS_SHARED_PTR<pdat::CellData<Real> > data_total_energy(
        HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
            d_patch->getPatchData(s_variable_total_energy, getDataContext())));
    
    return data_total_energy;
}


/*
 * Compute the cell data of density in the registered patch.
 */
void
FlowModelFourEqnConservative::computeCellDataOfDensity(
    const hier::Box& domain)
{
    if (d_num_subghosts_density > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_density)
        {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(d_data_density);
#endif
            // Get the cell data of the variable partial densities.
            HAMERS_SHARED_PTR<pdat::CellData<Real> > data_partial_densities =
                getCellDataOfPartialDensities();
            
            // Compute the density field.
            d_equation_of_state_mixing_rules->computeMixtureDensity(
                d_data_density,
                data_partial_densities,
                domain);
            
            d_cell_data_computed_density = true;
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::computeCellDataOfDensity()\n"
            << "Cell data of 'DENSITY' is not yet registered."
            << std::endl);
    }
}


/*
 * Compute the cell data of mass fractions with density in the registered patch.
 */
void
FlowModelFourEqnConservative::computeCellDataOfMassFractionsWithDensity(
    const hier::Box& domain)
{
    if (d_num_subghosts_mass_fractions > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_mass_fractions)
        {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(d_data_mass_fractions);
#endif
            
            /*
             * Get the local lower index and number of cells in each direction of the domain.
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
            HAMERS_SHARED_PTR<pdat::CellData<Real> > data_partial_densities =
                getCellDataOfPartialDensities();
            
            if (!d_cell_data_computed_density)
            {
                computeCellDataOfDensity(domain);
            }
            
            // Get the pointers to the cell data of mass fractions, denisty and partial densities.
            std::vector<Real*> Y;
            Y.reserve(d_num_species);
            for (int si = 0; si < d_num_species; si++)
            {
                Y.push_back(d_data_mass_fractions->getPointer(si));
            }
            Real* rho = d_data_density->getPointer(0);
            std::vector<Real*> rho_Y;
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
                    HAMERS_PRAGMA_SIMD
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
                        HAMERS_PRAGMA_SIMD
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
                            HAMERS_PRAGMA_SIMD
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
            
            d_cell_data_computed_mass_fractions = true;
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::computeCellDataOfMassFractionsWithDensity()\n"
            << "Cell data of 'MASS_FRACTIONS' is not yet registered."
            << std::endl);
    }
}


/*
 * Compute the cell data of mole fractions with mass fractions in the registered patch.
 */
void
FlowModelFourEqnConservative::computeCellDataOfMoleFractionsWithMassFractions(
    const hier::Box& domain)
{
    if (d_num_subghosts_mole_fractions > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_mole_fractions)
        {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(d_data_mole_fractions);
#endif
            
            /*
             * Get the local lower index and number of cells in each direction of the domain.
             */
            
            hier::IntVector domain_lo(d_dim);
            hier::IntVector domain_dims(d_dim);
            
            if (domain.empty())
            {
                domain_lo = -d_num_subghosts_mole_fractions;
                domain_dims = d_subghostcell_dims_mole_fractions;
            }
            else
            {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(d_subghost_box_mole_fractions.contains(domain));
#endif
                
                domain_lo = domain.lower() - d_interior_box.lower();
                domain_dims = domain.numberCells();
            }
            
            if (!d_cell_data_computed_mass_fractions)
            {
                computeCellDataOfMassFractionsWithDensity(domain);
            }
            
            // Get the pointers to the cell data of mole fractions and mass fractions.
            std::vector<Real*> X;
            X.reserve(d_num_species);
            for (int si = 0; si < d_num_species; si++)
            {
                X.push_back(d_data_mole_fractions->getPointer(si));
            }
            std::vector<Real*> Y;
            Y.reserve(d_num_species);
            for (int si = 0; si < d_num_species; si++)
            {
                Y.push_back(d_data_mass_fractions->getPointer(si));
            }
            
            // Compute the molecular weight of the mixture.
            
            HAMERS_SHARED_PTR<pdat::CellData<Real> > data_molecular_weight_mixture(
                new pdat::CellData<Real>(d_interior_box, 1, d_num_subghosts_mole_fractions));
            
            d_equation_of_state_mixing_rules->computeMixtureMolecularWeight(
                data_molecular_weight_mixture,
                d_data_mass_fractions,
                domain);
            
            Real* M = data_molecular_weight_mixture->getPointer(0);
            
            // Compute the mole fraction of each species.
            
            if (d_dim == tbox::Dimension(1))
            {
                /*
                 * Get the local lower index, numbers of cells in each dimension and numbers of ghost cells.
                 */
                
                const int domain_lo_0 = domain_lo[0];
                const int domain_dim_0 = domain_dims[0];
                
                const int num_subghosts_0_mass_fractions = d_num_subghosts_mass_fractions[0];
                const int num_subghosts_0_mole_fractions = d_num_subghosts_mole_fractions[0];
                
                // Compute the mole fraction field.
                for (int si = 0; si < d_num_species; si++)
                {
                    const Real M_i = d_equation_of_state_mixing_rules->
                        getSpeciesMolecularWeight(si);
                    
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_mass_fractions = i + num_subghosts_0_mass_fractions;
                        const int idx_mole_fractions = i + num_subghosts_0_mole_fractions;
                        
                        X[si][idx_mole_fractions] = Y[si][idx_mass_fractions]*M[idx_mole_fractions]/M_i;
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
                
                const int num_subghosts_0_mass_fractions = d_num_subghosts_mass_fractions[0];
                const int num_subghosts_1_mass_fractions = d_num_subghosts_mass_fractions[1];
                const int subghostcell_dim_0_mass_fractions = d_subghostcell_dims_mass_fractions[0];
                
                const int num_subghosts_0_mole_fractions = d_num_subghosts_mole_fractions[0];
                const int num_subghosts_1_mole_fractions = d_num_subghosts_mole_fractions[1];
                const int subghostcell_dim_0_mole_fractions = d_subghostcell_dims_mole_fractions[0];
                
                // Compute the mole fraction field.
                for (int si = 0; si < d_num_species; si++)
                {
                    const Real M_i = d_equation_of_state_mixing_rules->
                        getSpeciesMolecularWeight(si);
                    
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_mass_fractions = (i + num_subghosts_0_mass_fractions) +
                                (j + num_subghosts_1_mass_fractions)*subghostcell_dim_0_mass_fractions;
                            
                            const int idx_mole_fractions = (i + num_subghosts_0_mole_fractions) +
                                (j + num_subghosts_1_mole_fractions)*subghostcell_dim_0_mole_fractions;
                            
                            X[si][idx_mole_fractions] = Y[si][idx_mass_fractions]*M[idx_mole_fractions]/M_i;
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
                
                const int num_subghosts_0_mass_fractions = d_num_subghosts_mass_fractions[0];
                const int num_subghosts_1_mass_fractions = d_num_subghosts_mass_fractions[1];
                const int num_subghosts_2_mass_fractions = d_num_subghosts_mass_fractions[2];
                const int subghostcell_dim_0_mass_fractions = d_subghostcell_dims_mass_fractions[0];
                const int subghostcell_dim_1_mass_fractions = d_subghostcell_dims_mass_fractions[1];
                
                const int num_subghosts_0_mole_fractions = d_num_subghosts_mole_fractions[0];
                const int num_subghosts_1_mole_fractions = d_num_subghosts_mole_fractions[1];
                const int num_subghosts_2_mole_fractions = d_num_subghosts_mole_fractions[2];
                const int subghostcell_dim_0_mole_fractions = d_subghostcell_dims_mole_fractions[0];
                const int subghostcell_dim_1_mole_fractions = d_subghostcell_dims_mole_fractions[1];
                
                // Compute the mole fraction field.
                for (int si = 0; si < d_num_species; si++)
                {
                    const Real M_i = d_equation_of_state_mixing_rules->
                        getSpeciesMolecularWeight(si);
                    
                    for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                    {
                        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                        {
                            HAMERS_PRAGMA_SIMD
                            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                            {
                                // Compute the linear indices.
                                const int idx_mass_fractions = (i + num_subghosts_0_mass_fractions) +
                                    (j + num_subghosts_1_mass_fractions)*subghostcell_dim_0_mass_fractions +
                                    (k + num_subghosts_2_mass_fractions)*subghostcell_dim_0_mass_fractions*
                                        subghostcell_dim_1_mass_fractions;
                                
                                const int idx_mole_fractions = (i + num_subghosts_0_mole_fractions) +
                                    (j + num_subghosts_1_mole_fractions)*subghostcell_dim_0_mole_fractions +
                                    (k + num_subghosts_2_mole_fractions)*subghostcell_dim_0_mole_fractions*
                                        subghostcell_dim_1_mole_fractions;
                                
                                X[si][idx_mole_fractions] = Y[si][idx_mass_fractions]*M[idx_mole_fractions]/M_i;
                            }
                        }
                    }
                }
            }
            
            d_cell_data_computed_mole_fractions = true;
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::computeCellDataOfMoleFractionsWithMassFractions()\n"
            << "Cell data of 'MOLE_FRACTIONS' is not yet registered."
            << std::endl);
    }
}


/*
 * Compute the cell data of volume fractions with species densities in the registered patch.
 */
void
FlowModelFourEqnConservative::computeCellDataOfVolumeFractionsWithSpeciesDensities(
    const hier::Box& domain)
{
    if (d_num_subghosts_volume_fractions > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_volume_fractions)
        {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(d_data_volume_fractions);
#endif
            
            /*
             * Get the local lower index and number of cells in each direction of the domain.
             */
            
            hier::IntVector domain_lo(d_dim);
            hier::IntVector domain_dims(d_dim);
            
            if (domain.empty())
            {
                domain_lo = -d_num_subghosts_volume_fractions;
                domain_dims = d_subghostcell_dims_volume_fractions;
            }
            else
            {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(d_subghost_box_volume_fractions.contains(domain));
#endif
                
                domain_lo = domain.lower() - d_interior_box.lower();
                domain_dims = domain.numberCells();
            }
            
            // Get the cell data of the variable partial densities.
            HAMERS_SHARED_PTR<pdat::CellData<Real> > data_partial_densities =
                getCellDataOfPartialDensities();
            
            if (!d_cell_data_computed_species_densities)
            {
                computeCellDataOfSpeciesDensitiesWithPressureAndTemperature(domain);
            }
            
            // Get the pointers to the cell data of volume fractions, partial densities and species densities.
            std::vector<Real*> Z;
            Z.reserve(d_num_species);
            for (int si = 0; si < d_num_species; si++)
            {
                Z.push_back(d_data_volume_fractions->getPointer(si));
            }
            std::vector<Real*> rho_Y;
            rho_Y.reserve(d_num_species);
            for (int si = 0; si < d_num_species; si++)
            {
                rho_Y.push_back(data_partial_densities->getPointer(si));
            }
            std::vector<Real*> rho_i;
            rho_i.reserve(d_num_species);
            for (int si = 0; si < d_num_species; si++)
            {
                rho_i.push_back(d_data_species_densities[si]->getPointer(0));
            }
            
            if (d_dim == tbox::Dimension(1))
            {
                /*
                 * Get the local lower index, numbers of cells in each dimension and numbers of ghost cells.
                 */
                
                const int domain_lo_0 = domain_lo[0];
                const int domain_dim_0 = domain_dims[0];
                
                const int num_ghosts_0 = d_num_ghosts[0];
                const int num_subghosts_0_species_densities = d_num_subghosts_species_densities[0];
                const int num_subghosts_0_volume_fractions = d_num_subghosts_volume_fractions[0];
                
                // Compute the volume fraction field.
                for (int si = 0; si < d_num_species; si++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx = i + num_ghosts_0;
                        const int idx_species_densities = i + num_subghosts_0_species_densities;
                        const int idx_volume_fractions = i + num_subghosts_0_volume_fractions;
                        
                        Z[si][idx_volume_fractions] = rho_Y[si][idx]/rho_i[si][idx_species_densities];
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
                
                const int num_subghosts_0_species_densities = d_num_subghosts_species_densities[0];
                const int num_subghosts_1_species_densities = d_num_subghosts_species_densities[1];
                const int subghostcell_dim_0_species_densities = d_subghostcell_dims_species_densities[0];
                
                const int num_subghosts_0_volume_fractions = d_num_subghosts_volume_fractions[0];
                const int num_subghosts_1_volume_fractions = d_num_subghosts_volume_fractions[1];
                const int subghostcell_dim_0_volume_fractions = d_subghostcell_dims_volume_fractions[0];
                
                // Compute the volume fraction field.
                for (int si = 0; si < d_num_species; si++)
                {
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx = (i + num_ghosts_0) +
                                (j + num_ghosts_1)*ghostcell_dim_0;
                            
                            const int idx_species_densities = (i + num_subghosts_0_species_densities) +
                                (j + num_subghosts_1_species_densities)*subghostcell_dim_0_species_densities;
                            
                            const int idx_volume_fractions = (i + num_subghosts_0_volume_fractions) +
                                (j + num_subghosts_1_volume_fractions)*subghostcell_dim_0_volume_fractions;
                            
                            Z[si][idx_volume_fractions] = rho_Y[si][idx]/rho_i[si][idx_species_densities];
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
                
                const int num_subghosts_0_species_densities = d_num_subghosts_species_densities[0];
                const int num_subghosts_1_species_densities = d_num_subghosts_species_densities[1];
                const int num_subghosts_2_species_densities = d_num_subghosts_species_densities[2];
                const int subghostcell_dim_0_species_densities = d_subghostcell_dims_species_densities[0];
                const int subghostcell_dim_1_species_densities = d_subghostcell_dims_species_densities[1];
                
                const int num_subghosts_0_volume_fractions = d_num_subghosts_volume_fractions[0];
                const int num_subghosts_1_volume_fractions = d_num_subghosts_volume_fractions[1];
                const int num_subghosts_2_volume_fractions = d_num_subghosts_volume_fractions[2];
                const int subghostcell_dim_0_volume_fractions = d_subghostcell_dims_volume_fractions[0];
                const int subghostcell_dim_1_volume_fractions = d_subghostcell_dims_volume_fractions[1];
                
                // Compute the volume fraction field.
                for (int si = 0; si < d_num_species; si++)
                {
                    for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                    {
                        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                        {
                            HAMERS_PRAGMA_SIMD
                            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                            {
                                // Compute the linear indices.
                                const int idx = (i + num_ghosts_0) +
                                    (j + num_ghosts_1)*ghostcell_dim_0 +
                                    (k + num_ghosts_2)*ghostcell_dim_0*ghostcell_dim_1;
                                
                                const int idx_species_densities = (i + num_subghosts_0_species_densities) +
                                    (j + num_subghosts_1_species_densities)*subghostcell_dim_0_species_densities +
                                    (k + num_subghosts_2_species_densities)*subghostcell_dim_0_species_densities*
                                        subghostcell_dim_1_species_densities;
                                
                                const int idx_volume_fractions = (i + num_subghosts_0_volume_fractions) +
                                    (j + num_subghosts_1_volume_fractions)*subghostcell_dim_0_volume_fractions +
                                    (k + num_subghosts_2_volume_fractions)*subghostcell_dim_0_volume_fractions*
                                        subghostcell_dim_1_volume_fractions;
                                
                                Z[si][idx_volume_fractions] = rho_Y[si][idx]/rho_i[si][idx_species_densities];
                            }
                        }
                    }
                }
            }
            
            d_cell_data_computed_volume_fractions = true;
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::computeCellDataOfVolumeFractionsWithSpeciesDensities()\n"
            << "Cell data of 'VOLUME_FRACTIONS' is not yet registered."
            << std::endl);
    }
}


/*
 * Compute the cell data of velocity with density in the registered patch.
 */
void
FlowModelFourEqnConservative::computeCellDataOfVelocityWithDensity(
    const hier::Box& domain)
{
    if (d_num_subghosts_velocity > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_velocity)
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
            
            // Get the cell data of the variable momentum.
            HAMERS_SHARED_PTR<pdat::CellData<Real> > data_momentum =
                getCellDataOfMomentum();
            
            if (!d_cell_data_computed_density)
            {
                computeCellDataOfDensity(domain);
            }
            
            // Get the pointer to the cell data of density.
            Real* rho = d_data_density->getPointer(0);
            
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
                Real* u = d_data_velocity->getPointer(0);
                
                // Get the pointer to the cell data of momentum.
                Real* rho_u = data_momentum->getPointer(0);
                
                // Compute the velocity field.
                HAMERS_PRAGMA_SIMD
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
                Real* u = d_data_velocity->getPointer(0);
                Real* v = d_data_velocity->getPointer(1);
                
                // Get the pointers to the cell data of momentum.
                Real* rho_u = data_momentum->getPointer(0);
                Real* rho_v = data_momentum->getPointer(1);
                
                // Compute the velocity field.
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
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
                Real* u = d_data_velocity->getPointer(0);
                Real* v = d_data_velocity->getPointer(1);
                Real* w = d_data_velocity->getPointer(2);
                
                // Get the pointers to the cell data of momentum.
                Real* rho_u = data_momentum->getPointer(0);
                Real* rho_v = data_momentum->getPointer(1);
                Real* rho_w = data_momentum->getPointer(2);
                
                // Compute the velocity field.
                for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                {
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
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
            
            d_cell_data_computed_velocity = true;
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::computeCellDataOfVelocityWithDensity()\n"
            << "Cell data of 'VELOCITY' is not yet registered."
            << std::endl);
    }
}


/*
 * Compute the cell data of internal energy with density and velocity in the registered
 * patch.
 */
void
FlowModelFourEqnConservative::computeCellDataOfInternalEnergyWithDensityAndVelocity(
    const hier::Box& domain)
{
    if (d_num_subghosts_internal_energy > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_internal_energy)
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
            
            // Get the cell data of the variable total energy.
            HAMERS_SHARED_PTR<pdat::CellData<Real> > data_total_energy =
                getCellDataOfTotalEnergy();
            
            if (!d_cell_data_computed_density)
            {
                computeCellDataOfDensity(domain);
            }
            
            if (!d_cell_data_computed_velocity)
            {
                computeCellDataOfVelocityWithDensity(domain);
            }
            
            // Get the pointers to the cell data of internal energy, total energy and density.
            Real* epsilon = d_data_internal_energy->getPointer(0);
            Real* E = data_total_energy->getPointer(0);
            Real* rho = d_data_density->getPointer(0);
            
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
                Real* u = d_data_velocity->getPointer(0);
                
                // Compute the internal energy field.
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx = i + num_ghosts_0;
                    const int idx_density = i + num_subghosts_0_density;
                    const int idx_velocity = i + num_subghosts_0_velocity;
                    const int idx_internal_energy = i + num_subghosts_0_internal_energy;
                    
                    epsilon[idx_internal_energy] = E[idx]/rho[idx_density] -
                        Real(1)/Real(2)*u[idx_velocity]*u[idx_velocity];
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
                Real* u = d_data_velocity->getPointer(0);
                Real* v = d_data_velocity->getPointer(1);
                
                // Compute the internal energy field.
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
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
                            Real(1)/Real(2)*(u[idx_velocity]*u[idx_velocity] + v[idx_velocity]*v[idx_velocity]);
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
                Real* u = d_data_velocity->getPointer(0);
                Real* v = d_data_velocity->getPointer(1);
                Real* w = d_data_velocity->getPointer(2);
                
                // Compute the internal energy field.
                for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                {
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
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
                                Real(1)/Real(2)*(u[idx_velocity]*u[idx_velocity] + v[idx_velocity]*v[idx_velocity] +
                                    w[idx_velocity]*w[idx_velocity]);
                        }
                    }
                }
            }
            
            d_cell_data_computed_internal_energy = true;
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::"
            << "computeCellDataOfInternalEnergyWithDensityAndVelocity()\n"
            << "Cell data of 'INTERNAL_ENERGY' is not yet registered."
            << std::endl);
    }
}


/*
 * Compute the cell data of pressure with density, mass fractions and internal energy in
 * the registered patch.
 */
void
FlowModelFourEqnConservative::computeCellDataOfPressureWithDensityMassFractionsAndInternalEnergy(
    const hier::Box& domain)
{
    if (d_num_subghosts_pressure > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_pressure)
        {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(d_data_pressure);
#endif
            
            if (!d_cell_data_computed_density)
            {
                computeCellDataOfDensity(domain);
            }
            
            if (!d_cell_data_computed_mass_fractions)
            {
                computeCellDataOfMassFractionsWithDensity(domain);
            }
            
            if (!d_cell_data_computed_internal_energy)
            {
                computeCellDataOfInternalEnergyWithDensityAndVelocity(domain);
            }
            
            // Compute the pressure field.
            d_equation_of_state_mixing_rules->computePressure(
                d_data_pressure,
                d_data_density,
                d_data_internal_energy,
                d_data_mass_fractions,
                domain);
            
            d_cell_data_computed_pressure = true;
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::"
            << "computeCellDataOfPressureWithDensityMassFractionsAndInternalEnergy()\n"
            << "Cell data of 'PRESSURE' is not yet registered."
            << std::endl);
    }
}


/*
 * Compute the cell data of sound speed with density, mass fractions and pressure in the
 * registered patch.
 */
void
FlowModelFourEqnConservative::computeCellDataOfSoundSpeedWithDensityMassFractionsAndPressure(
    const hier::Box& domain)
{
    if (d_num_subghosts_sound_speed > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_sound_speed)
        {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(d_data_sound_speed);
#endif
            
            /*
             * Get the local lower index and number of cells in each direction of the domain.
             */
            
            hier::IntVector domain_lo(d_dim);
            hier::IntVector domain_dims(d_dim);
            
            if (domain.empty())
            {
                domain_lo = -d_num_subghosts_sound_speed;
                domain_dims = d_subghostcell_dims_sound_speed;
            }
            else
            {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(d_subghost_box_sound_speed.contains(domain));
#endif
                
                domain_lo = domain.lower() - d_interior_box.lower();
                domain_dims = domain.numberCells();
            }
            
            if (!d_cell_data_computed_density)
            {
                computeCellDataOfDensity(domain);
            }
            
            if (!d_cell_data_computed_mass_fractions)
            {
                computeCellDataOfMassFractionsWithDensity(domain);
            }
            
            if (!d_cell_data_computed_pressure)
            {
                computeCellDataOfPressureWithDensityMassFractionsAndInternalEnergy(domain);
            }
            
            // Compute the partial derivatives.
            HAMERS_SHARED_PTR<pdat::CellData<Real> > data_gruneisen_parameter(
                new pdat::CellData<Real>(d_interior_box, 1, d_num_subghosts_sound_speed));
            
            HAMERS_SHARED_PTR<pdat::CellData<Real> > data_partial_pressure_partial_partial_densities(
                new pdat::CellData<Real>(d_interior_box, d_num_species, d_num_subghosts_sound_speed));
            
            d_equation_of_state_mixing_rules->computeGruneisenParameter(
                data_gruneisen_parameter,
                d_data_density,
                d_data_pressure,
                d_data_mass_fractions,
                domain);
            
            d_equation_of_state_mixing_rules->computePressureDerivativeWithPartialDensities(
                data_partial_pressure_partial_partial_densities,
                d_data_density,
                d_data_pressure,
                d_data_mass_fractions,
                domain);
            
            // Get the pointers to the cell data of sound speed, density, mass fractions, pressure,
            // Gruneisen parameter and partial pressure partial partial densities.
            Real* c     = d_data_sound_speed->getPointer(0);
            Real* rho   = d_data_density->getPointer(0);
            Real* p     = d_data_pressure->getPointer(0);
            Real* Gamma = data_gruneisen_parameter->getPointer(0);
            std::vector<Real*> Y;
            std::vector<Real*> Psi;
            Y.reserve(d_num_species);
            Psi.reserve(d_num_species);
            for (int si = 0; si < d_num_species; si++)
            {
                Y.push_back(d_data_mass_fractions->getPointer(si));
                Psi.push_back(data_partial_pressure_partial_partial_densities->getPointer(si));
            }
            
            // Compute the sound speed field.
            
            if (d_dim == tbox::Dimension(1))
            {
                /*
                 * Get the local lower index, numbers of cells in each dimension and numbers of ghost cells.
                 */
                
                const int domain_lo_0 = domain_lo[0];
                const int domain_dim_0 = domain_dims[0];
                
                const int num_subghosts_0_density = d_num_subghosts_density[0];
                const int num_subghosts_0_mass_fractions = d_num_subghosts_mass_fractions[0];
                const int num_subghosts_0_pressure = d_num_subghosts_pressure[0];
                const int num_subghosts_0_sound_speed = d_num_subghosts_sound_speed[0];
                
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_density = i + num_subghosts_0_density;
                    const int idx_pressure = i + num_subghosts_0_pressure;
                    const int idx_sound_speed = i + num_subghosts_0_sound_speed;
                    
                    c[idx_sound_speed] = Gamma[idx_sound_speed]*p[idx_pressure]/rho[idx_density];
                }
                
                for (int si = 0; si < d_num_species; si++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_mass_fractions = i + num_subghosts_0_mass_fractions;
                        const int idx_sound_speed = i + num_subghosts_0_sound_speed;
                        
                        c[idx_sound_speed] += Y[si][idx_mass_fractions]*Psi[si][idx_sound_speed];
                    }
                }
                
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_sound_speed = i + num_subghosts_0_sound_speed;
                    
                    c[idx_sound_speed] = std::sqrt(c[idx_sound_speed]);
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
                
                const int num_subghosts_0_density = d_num_subghosts_density[0];
                const int num_subghosts_1_density = d_num_subghosts_density[1];
                const int subghostcell_dim_0_density = d_subghostcell_dims_density[0];
                
                const int num_subghosts_0_mass_fractions = d_num_subghosts_mass_fractions[0];
                const int num_subghosts_1_mass_fractions = d_num_subghosts_mass_fractions[1];
                const int subghostcell_dim_0_mass_fractions = d_subghostcell_dims_mass_fractions[0];
                
                const int num_subghosts_0_pressure = d_num_subghosts_pressure[0];
                const int num_subghosts_1_pressure = d_num_subghosts_pressure[1];
                const int subghostcell_dim_0_pressure = d_subghostcell_dims_pressure[0];
                
                const int num_subghosts_0_sound_speed = d_num_subghosts_sound_speed[0];
                const int num_subghosts_1_sound_speed = d_num_subghosts_sound_speed[1];
                const int subghostcell_dim_0_sound_speed = d_subghostcell_dims_sound_speed[0];
                
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_density = (i + num_subghosts_0_density) +
                            (j + num_subghosts_1_density)*subghostcell_dim_0_density;
                        
                        const int idx_pressure = (i + num_subghosts_0_pressure) +
                            (j + num_subghosts_1_pressure)*subghostcell_dim_0_pressure;
                        
                        const int idx_sound_speed = (i + num_subghosts_0_sound_speed) +
                            (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed;
                        
                        c[idx_sound_speed] = Gamma[idx_sound_speed]*p[idx_pressure]/rho[idx_density];
                    }
                }
                
                for (int si = 0; si < d_num_species; si++)
                {
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_mass_fractions = (i + num_subghosts_0_mass_fractions) +
                                (j + num_subghosts_1_mass_fractions)*subghostcell_dim_0_mass_fractions;
                            
                            const int idx_sound_speed = (i + num_subghosts_0_sound_speed) +
                                (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed;
                            
                            c[idx_sound_speed] += Y[si][idx_mass_fractions]*Psi[si][idx_sound_speed];
                        }
                    }
                }
                
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx_sound_speed = (i + num_subghosts_0_sound_speed) +
                            (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed;
                        
                        c[idx_sound_speed] = std::sqrt(c[idx_sound_speed]);
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
                
                const int num_subghosts_0_pressure = d_num_subghosts_pressure[0];
                const int num_subghosts_1_pressure = d_num_subghosts_pressure[1];
                const int num_subghosts_2_pressure = d_num_subghosts_pressure[2];
                const int subghostcell_dim_0_pressure = d_subghostcell_dims_pressure[0];
                const int subghostcell_dim_1_pressure = d_subghostcell_dims_pressure[1];
                
                const int num_subghosts_0_sound_speed = d_num_subghosts_sound_speed[0];
                const int num_subghosts_1_sound_speed = d_num_subghosts_sound_speed[1];
                const int num_subghosts_2_sound_speed = d_num_subghosts_sound_speed[2];
                const int subghostcell_dim_0_sound_speed = d_subghostcell_dims_sound_speed[0];
                const int subghostcell_dim_1_sound_speed = d_subghostcell_dims_sound_speed[1];
                
                for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                {
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_density = (i + num_subghosts_0_density) +
                                (j + num_subghosts_1_density)*subghostcell_dim_0_density +
                                (k + num_subghosts_2_density)*subghostcell_dim_0_density*
                                    subghostcell_dim_1_density;
                            
                            const int idx_pressure = (i + num_subghosts_0_pressure) +
                                (j + num_subghosts_1_pressure)*subghostcell_dim_0_pressure +
                                (k + num_subghosts_2_pressure)*subghostcell_dim_0_pressure*
                                    subghostcell_dim_1_pressure;
                            
                            const int idx_sound_speed = (i + num_subghosts_0_sound_speed) +
                                (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed +
                                (k + num_subghosts_2_sound_speed)*subghostcell_dim_0_sound_speed*
                                    subghostcell_dim_1_sound_speed;
                            
                            c[idx_sound_speed] = Gamma[idx_sound_speed]*p[idx_pressure]/rho[idx_density];
                        }
                    }
                }
                
                for (int si = 0; si < d_num_species; si++)
                {
                    for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                    {
                        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                        {
                            HAMERS_PRAGMA_SIMD
                            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                            {
                                // Compute the linear indices.
                                const int idx_mass_fractions = (i + num_subghosts_0_mass_fractions) +
                                    (j + num_subghosts_1_mass_fractions)*subghostcell_dim_0_mass_fractions +
                                    (k + num_subghosts_2_mass_fractions)*subghostcell_dim_0_mass_fractions*
                                        subghostcell_dim_1_mass_fractions;
                                
                                const int idx_sound_speed = (i + num_subghosts_0_sound_speed) +
                                (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed +
                                (k + num_subghosts_2_sound_speed)*subghostcell_dim_0_sound_speed*
                                    subghostcell_dim_1_sound_speed;
                                
                                c[idx_sound_speed] += Y[si][idx_mass_fractions]*Psi[si][idx_sound_speed];
                            }
                        }
                    }
                }
                
                for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                {
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute the linear index.
                            const int idx_sound_speed = (i + num_subghosts_0_sound_speed) +
                                (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed +
                                (k + num_subghosts_2_sound_speed)*subghostcell_dim_0_sound_speed*
                                    subghostcell_dim_1_sound_speed;
                            
                            c[idx_sound_speed] = std::sqrt(c[idx_sound_speed]);
                        }
                    }
                }
            }
            
            d_cell_data_computed_sound_speed = true;
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::"
            << "computeCellDataOfSoundSpeedWithDensityMassFractionsAndPressure()\n"
            << "Cell data of 'SOUND_SPEED' is not yet registered."
            << std::endl);
    }
}


/*
 * Compute the cell data of temperature with density, mass fractions and pressure in the
 * registered patch.
 */
void
FlowModelFourEqnConservative::computeCellDataOfTemperatureWithDensityMassFractionsAndPressure(
    const hier::Box& domain)
{
    if (d_num_subghosts_temperature > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_temperature)
        {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(d_data_temperature);
#endif
            
            if (!d_cell_data_computed_density)
            {
                computeCellDataOfDensity(domain);
            }
            
            if (!d_cell_data_computed_mass_fractions)
            {
                computeCellDataOfMassFractionsWithDensity(domain);
            }
            
            if (!d_cell_data_computed_pressure)
            {
                computeCellDataOfPressureWithDensityMassFractionsAndInternalEnergy(domain);
            }
            
            // Compute the temperature field.
            d_equation_of_state_mixing_rules->computeTemperature(
                d_data_temperature,
                d_data_density,
                d_data_pressure,
                d_data_mass_fractions,
                domain);
            
            d_cell_data_computed_temperature = true;
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::"
            << "computeCellDataOfTemperatureWithDensityMassFractionsAndPressure()\n"
            << "Cell data of 'TEMPERATURE' is not yet registered."
            << std::endl);
    }
}


/*
 * Compute the cell data of convective flux with velocity and pressure in the registered
 * patch.
 */
void
FlowModelFourEqnConservative::computeCellDataOfConvectiveFluxWithVelocityAndPressure(
    const DIRECTION::TYPE& direction,
    const hier::Box& domain)
{
    if (direction == DIRECTION::X_DIRECTION)
    {
        if (d_num_subghosts_convective_flux_x > -hier::IntVector::getOne(d_dim))
        {
            if (!d_cell_data_computed_convective_flux_x)
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
                std::vector<Real*> F_x;
                F_x.reserve(d_num_eqn);
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    F_x.push_back(d_data_convective_flux_x->getPointer(ei));
                }
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_partial_densities =
                    getCellDataOfPartialDensities();
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_momentum =
                    getCellDataOfMomentum();
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_total_energy =
                    getCellDataOfTotalEnergy();
                
                if (!d_cell_data_computed_velocity)
                {
                    computeCellDataOfVelocityWithDensity(domain);
                }
                
                if (!d_cell_data_computed_pressure)
                {
                    computeCellDataOfPressureWithDensityMassFractionsAndInternalEnergy(domain);
                }
                
                // Get the pointers to the cell data of partial densities, total energy and pressure.
                std::vector<Real*> rho_Y;
                rho_Y.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    rho_Y.push_back(data_partial_densities->getPointer(si));
                }
                Real* E = data_total_energy->getPointer(0);
                Real* p = d_data_pressure->getPointer(0);
                
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
                    Real* rho_u = data_momentum->getPointer(0);
                    
                    // Get the pointer to the cell data of velocity.
                    Real* u = d_data_velocity->getPointer(0);
                    
                    // Compute the convective flux in the x-direction.
                    for (int si = 0; si < d_num_species; si++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx = i + num_ghosts_0;
                            const int idx_velocity = i + num_subghosts_0_velocity;
                            const int idx_convective_flux_x = i + num_subghosts_0_convective_flux_x;
                            
                            F_x[si][idx_convective_flux_x] = u[idx_velocity]*rho_Y[si][idx];
                        }
                    }
                    
                    HAMERS_PRAGMA_SIMD
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
                    Real* rho_u = data_momentum->getPointer(0);
                    Real* rho_v = data_momentum->getPointer(1);
                    
                    // Get the pointer to the cell data of velocity.
                    Real* u = d_data_velocity->getPointer(0);
                    
                    // Compute the convective flux in the x-direction.
                    for (int si = 0; si < d_num_species; si++)
                    {
                        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                        {
                            HAMERS_PRAGMA_SIMD
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
                        HAMERS_PRAGMA_SIMD
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
                    Real* rho_u = data_momentum->getPointer(0);
                    Real* rho_v = data_momentum->getPointer(1);
                    Real* rho_w = data_momentum->getPointer(2);
                    
                    // Get the pointer to the cell data of velocity.
                    Real* u = d_data_velocity->getPointer(0);
                    
                    // Compute the convective flux in the x-direction.
                    for (int si = 0; si < d_num_species; si++)
                    {
                        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                        {
                            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                            {
                                HAMERS_PRAGMA_SIMD
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
                            HAMERS_PRAGMA_SIMD
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
                
                d_cell_data_computed_convective_flux_x = true;
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::"
                << "computeCellDataOfConvectiveFluxWithVelocityAndPressure()\n"
                << "Cell data of 'CONVECTIVE_FLUX_X' is not yet registered."
                << std::endl);
        }
    }
    else if (direction == DIRECTION::Y_DIRECTION)
    {
        if (d_num_subghosts_convective_flux_y > -hier::IntVector::getOne(d_dim))
        {
            if (!d_cell_data_computed_convective_flux_y)
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
                std::vector<Real*> F_y;
                F_y.reserve(d_num_eqn);
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    F_y.push_back(d_data_convective_flux_y->getPointer(ei));
                }
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_partial_densities =
                    getCellDataOfPartialDensities();
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_momentum =
                    getCellDataOfMomentum();
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_total_energy =
                    getCellDataOfTotalEnergy();
                
                if (!d_cell_data_computed_pressure)
                {
                    computeCellDataOfPressureWithDensityMassFractionsAndInternalEnergy(domain);
                }
                
                if (!d_cell_data_computed_velocity)
                {
                    computeCellDataOfVelocityWithDensity(domain);
                }
                
                // Get the pointers to the cell data of partial densities, total energy, volume fraction
                // and pressure.
                std::vector<Real*> rho_Y;
                rho_Y.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    rho_Y.push_back(data_partial_densities->getPointer(si));
                }
                Real* E = data_total_energy->getPointer(0);
                Real* p = d_data_pressure->getPointer(0);
                
                if (d_dim == tbox::Dimension(1))
                {
                    TBOX_ERROR(d_object_name
                        << ": FlowModelFourEqnConservative::"
                        << "computeCellDataOfConvectiveFluxWithVelocityAndPressure()\n"
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
                    Real* rho_u = data_momentum->getPointer(0);
                    Real* rho_v = data_momentum->getPointer(1);
                    
                    // Get the pointer to the cell data of velocity.
                    Real* v = d_data_velocity->getPointer(1);
                    
                    // Compute the convective flux in the y-direction.
                    for (int si = 0; si < d_num_species; si++)
                    {
                        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                        {
                            HAMERS_PRAGMA_SIMD
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
                        HAMERS_PRAGMA_SIMD
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
                    Real* rho_u = data_momentum->getPointer(0);
                    Real* rho_v = data_momentum->getPointer(1);
                    Real* rho_w = data_momentum->getPointer(2);
                    
                    // Get the pointer to the cell data of velocity.
                    Real* v = d_data_velocity->getPointer(1);
                    
                    // Compute the convective flux in the y-direction.
                    for (int si = 0; si < d_num_species; si++)
                    {
                        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                        {
                            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                            {
                                HAMERS_PRAGMA_SIMD
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
                            HAMERS_PRAGMA_SIMD
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
                
                d_cell_data_computed_convective_flux_y = true;
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::"
                << "computeCellDataOfConvectiveFluxWithVelocityAndPressure()\n"
                << "Cell data of 'CONVECTIVE_FLUX_Y' is not yet registered."
                << std::endl);
        }
    }
    else if (direction == DIRECTION::Z_DIRECTION)
    {
        if (d_num_subghosts_convective_flux_z > -hier::IntVector::getOne(d_dim))
        {
            if (!d_cell_data_computed_convective_flux_z)
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
                std::vector<Real*> F_z;
                F_z.reserve(d_num_eqn);
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    F_z.push_back(d_data_convective_flux_z->getPointer(ei));
                }
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_partial_densities =
                    getCellDataOfPartialDensities();
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_momentum =
                    getCellDataOfMomentum();
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_total_energy =
                    getCellDataOfTotalEnergy();
                
                if (!d_cell_data_computed_pressure)
                {
                    computeCellDataOfPressureWithDensityMassFractionsAndInternalEnergy(domain);
                }
                
                if (!d_cell_data_computed_velocity)
                {
                    computeCellDataOfVelocityWithDensity(domain);
                }
                
                // Get the pointers to the cell data of partial densities, total energy, volume fraction
                // and pressure.
                std::vector<Real*> rho_Y;
                rho_Y.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    rho_Y.push_back(data_partial_densities->getPointer(si));
                }
                Real* E = data_total_energy->getPointer(0);
                Real* p = d_data_pressure->getPointer(0);
                
                if (d_dim == tbox::Dimension(1) || d_dim == tbox::Dimension(2))
                {
                    TBOX_ERROR(d_object_name
                        << ": FlowModelFourEqnConservative::"
                        << "computeCellDataOfConvectiveFluxWithVelocityAndPressure()\n"
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
                    Real* rho_u = data_momentum->getPointer(0);
                    Real* rho_v = data_momentum->getPointer(1);
                    Real* rho_w = data_momentum->getPointer(2);
                    
                    // Get the pointer to the cell data of velocity.
                    Real* w = d_data_velocity->getPointer(2);
                    
                    // Compute the convective flux in the z-direction.
                    for (int si = 0; si < d_num_species; si++)
                    {
                        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                        {
                            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                            {
                                HAMERS_PRAGMA_SIMD
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
                            HAMERS_PRAGMA_SIMD
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
                
                d_cell_data_computed_convective_flux_z = true;
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::"
                << "computeCellDataOfConvectiveFluxWithVelocityAndPressure()\n"
                << "Cell data of 'CONVECTIVE_FLUX_Z' is not yet registered."
                << std::endl);
        }
    }
}


/*
 * Compute the cell data of maximum wave speed with velocity and sound speed in the
 * registered patch.
 */
void
FlowModelFourEqnConservative::computeCellDataOfMaxWaveSpeedWithVelocityAndSoundSpeed(
    const DIRECTION::TYPE& direction,
    const hier::Box& domain)
{
    if (direction == DIRECTION::X_DIRECTION)
    {
        if (d_num_subghosts_max_wave_speed_x > -hier::IntVector::getOne(d_dim))
        {
            if (!d_cell_data_computed_max_wave_speed_x)
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
                
                if (!d_cell_data_computed_velocity)
                {
                    computeCellDataOfVelocityWithDensity(domain);
                }
                
                if (!d_cell_data_computed_sound_speed)
                {
                    computeCellDataOfSoundSpeedWithDensityMassFractionsAndPressure(domain);
                }
                
                // Get the pointers to the cell data of maximum wave speed and velocity in x-direction and sound speed.
                Real* lambda_max_x = d_data_max_wave_speed_x->getPointer(0);
                Real* u            = d_data_velocity->getPointer(0);
                Real* c            = d_data_sound_speed->getPointer(0);
                
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
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_sound_speed = i + num_subghosts_0_sound_speed;
                        const int idx_velocity = i + num_subghosts_0_velocity;
                        const int idx_max_wave_speed_x = i + num_subghosts_0_max_wave_speed_x;
                        
                        lambda_max_x[idx_max_wave_speed_x] = std::abs(u[idx_velocity]) + c[idx_sound_speed];
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
                        HAMERS_PRAGMA_SIMD
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_sound_speed = (i + num_subghosts_0_sound_speed) +
                                (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed;
                            
                            const int idx_velocity = (i + num_subghosts_0_velocity) +
                                (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                            
                            const int idx_max_wave_speed_x = (i + num_subghosts_0_max_wave_speed_x) +
                                (j + num_subghosts_1_max_wave_speed_x)*subghostcell_dim_0_max_wave_speed_x;
                            
                            lambda_max_x[idx_max_wave_speed_x] = std::abs(u[idx_velocity]) + c[idx_sound_speed];
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
                            HAMERS_PRAGMA_SIMD
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
                                
                                lambda_max_x[idx_max_wave_speed_x] = std::abs(u[idx_velocity]) + c[idx_sound_speed];
                            }
                        }
                    }
                }
                
                d_cell_data_computed_max_wave_speed_x = true;
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::"
                << "computeCellDataOfMaxWaveSpeedWithVelocityAndSoundSpeed()\n"
                << "Cell data of 'MAX_WAVE_SPEED_X' is not yet registered."
                << std::endl);
        }
    }
    else if (direction == DIRECTION::Y_DIRECTION)
    {
        if (d_num_subghosts_max_wave_speed_y > -hier::IntVector::getOne(d_dim))
        {
            if (!d_cell_data_computed_max_wave_speed_y)
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
                
                if (!d_cell_data_computed_sound_speed)
                {
                    computeCellDataOfSoundSpeedWithDensityMassFractionsAndPressure(domain);
                }
                
                if (!d_cell_data_computed_velocity)
                {
                    computeCellDataOfVelocityWithDensity(domain);
                }
                
                // Get the pointers to the cell data of maximum wave speed and velocity in y-direction, and sound speed.
                Real* lambda_max_y = d_data_max_wave_speed_y->getPointer(0);
                Real* v            = d_data_velocity->getPointer(1);
                Real* c            = d_data_sound_speed->getPointer(0);
                
                if (d_dim == tbox::Dimension(1))
                {
                    TBOX_ERROR(d_object_name
                        << ": FlowModelFourEqnConservative::"
                        << "computeCellDataOfMaxWaveSpeedWithVelocityAndSoundSpeed()\n"
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
                        HAMERS_PRAGMA_SIMD
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_sound_speed = (i + num_subghosts_0_sound_speed) +
                                (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed;
                            
                            const int idx_velocity = (i + num_subghosts_0_velocity) +
                                (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                            
                            const int idx_max_wave_speed_y = (i + num_subghosts_0_max_wave_speed_y) +
                                (j + num_subghosts_1_max_wave_speed_y)*subghostcell_dim_0_max_wave_speed_y;
                            
                            lambda_max_y[idx_max_wave_speed_y] = std::abs(v[idx_velocity]) + c[idx_sound_speed];
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
                            HAMERS_PRAGMA_SIMD
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
                                
                                lambda_max_y[idx_max_wave_speed_y] = std::abs(v[idx_velocity]) + c[idx_sound_speed];
                            }
                        }
                    }
                }
                
                d_cell_data_computed_max_wave_speed_y = true;
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::"
                << "computeCellDataOfMaxWaveSpeedWithVelocityAndSoundSpeed()\n"
                << "Cell data of 'MAX_WAVE_SPEED_Y' is not yet registered."
                << std::endl);
        }
    }
    else if (direction == DIRECTION::Z_DIRECTION)
    {
        if (d_num_subghosts_max_wave_speed_z > -hier::IntVector::getOne(d_dim))
        {
            if (!d_cell_data_computed_max_wave_speed_z)
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
                
                if (!d_cell_data_computed_sound_speed)
                {
                    computeCellDataOfSoundSpeedWithDensityMassFractionsAndPressure(domain);
                }
                
                if (!d_cell_data_computed_velocity)
                {
                    computeCellDataOfVelocityWithDensity(domain);
                }
                
                // Get the pointers to the cell data of maximum wave speed and velocity in z-direction, and sound speed.
                Real* lambda_max_z = d_data_max_wave_speed_z->getPointer(0);
                Real* w            = d_data_velocity->getPointer(2);
                Real* c            = d_data_sound_speed->getPointer(0);
                
                if (d_dim == tbox::Dimension(1) || d_dim == tbox::Dimension(2))
                {
                    TBOX_ERROR(d_object_name
                        << ": FlowModelFourEqnConservative::"
                        << "computeCellDataOfMaxWaveSpeedWithVelocityAndSoundSpeed()\n"
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
                            HAMERS_PRAGMA_SIMD
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
                                
                                lambda_max_z[idx_max_wave_speed_z] = std::abs(w[idx_velocity]) + c[idx_sound_speed];
                            }
                        }
                    }
                }
                
                d_cell_data_computed_max_wave_speed_z = true;
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFourEqnConservative::"
                << "computeCellDataOfMaxWaveSpeedWithVelocityAndSoundSpeed()\n"
                << "Cell data of 'MAX_WAVE_SPEED_Z' is not yet registered."
                << std::endl);
        }
    }
}


/*
 * Compute the cell data of maximum diffusivity with density, mass fractions, pressure
 * and temperature in the registered patch.
 */
void
FlowModelFourEqnConservative::computeCellDataOfMaxDiffusivityWithDensityMassFractionsPressureAndTemperature(
    const hier::Box& domain)
{
    if (!d_equation_of_mass_diffusivity_mixing_rules ||
        !d_equation_of_shear_viscosity_mixing_rules ||
        !d_equation_of_bulk_viscosity_mixing_rules ||
        !d_equation_of_thermal_conductivity_mixing_rules)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::"
            << "computeCellDataOfMaxDiffusivityWithDensityMassFractionsPressureAndTemperature()\n"
            << "Either mixing rule of mass diffusivity, shear viscosity, bulk viscosity or"
            << " thermal conductivity is not initialized."
            << std::endl);
    }
    
    if (d_num_subghosts_max_diffusivity > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_max_diffusivity)
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
            
            if (!d_cell_data_computed_density)
            {
                computeCellDataOfDensity(domain);
            }
            
            if (!d_cell_data_computed_mass_fractions)
            {
                computeCellDataOfMassFractionsWithDensity(domain);
            }
            
            if (!d_cell_data_computed_pressure)
            {
                computeCellDataOfPressureWithDensityMassFractionsAndInternalEnergy(domain);
            }
            
            if (!d_cell_data_computed_temperature)
            {
                computeCellDataOfTemperatureWithDensityMassFractionsAndPressure(domain);
            }
            
            /*
             * Create temporary cell data of isobaric specific heat capacity, mass diffusivities, shear
             * viscosity, bulk viscosity and thermal conductivity.
             */
            
            HAMERS_SHARED_PTR<pdat::CellData<Real> > data_isobaric_specific_heat_capacity(
                new pdat::CellData<Real>(d_interior_box, 1, d_num_subghosts_max_diffusivity));
            
            HAMERS_SHARED_PTR<pdat::CellData<Real> > data_mass_diffusivities(
                new pdat::CellData<Real>(d_interior_box, d_num_species, d_num_subghosts_max_diffusivity));
            
            HAMERS_SHARED_PTR<pdat::CellData<Real> > data_shear_viscosity(
                new pdat::CellData<Real>(d_interior_box, 1, d_num_subghosts_max_diffusivity));
            
            HAMERS_SHARED_PTR<pdat::CellData<Real> > data_bulk_viscosity(
                new pdat::CellData<Real>(d_interior_box, 1, d_num_subghosts_max_diffusivity));
            
            HAMERS_SHARED_PTR<pdat::CellData<Real> > data_thermal_conductivity(
                new pdat::CellData<Real>(d_interior_box, 1, d_num_subghosts_max_diffusivity));
            
            // Get the pointers to the cell data of maximum diffusivity, density, isobaric specific heat
            // capacity, mass diffusivities, shear viscosity, bulk viscosity and thermal conductivity.
            Real* D_max = d_data_max_diffusivity->getPointer(0);
            Real* rho   = d_data_density->getPointer(0);
            Real* c_p   = data_isobaric_specific_heat_capacity->getPointer(0);
            std::vector<Real*> D;
            D.reserve(d_num_species);
            for (int si = 0; si < d_num_species; si++)
            {
                D.push_back(data_mass_diffusivities->getPointer(si));
            }
            Real* mu    = data_shear_viscosity->getPointer(0);
            Real* mu_v  = data_bulk_viscosity->getPointer(0);
            Real* kappa = data_thermal_conductivity->getPointer(0);
            
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
                
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_max_diffusivity = i + num_subghosts_0_max_diffusivity;
                    const int idx_density = i + num_subghosts_0_density;
                    
                    D_max[idx_max_diffusivity] = std::max(mu[idx_max_diffusivity]/rho[idx_density],
                        mu_v[idx_max_diffusivity]/rho[idx_density]);
                    
                    D_max[idx_max_diffusivity] = std::max(D_max[idx_max_diffusivity],
                        kappa[idx_max_diffusivity]/(rho[idx_density]*c_p[idx_max_diffusivity]));
                }
                
                for (int si = 0; si < d_num_species; si++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx_max_diffusivity = i + num_subghosts_0_max_diffusivity;
                        
                        D_max[idx_max_diffusivity] = std::max(D_max[idx_max_diffusivity],
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
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_max_diffusivity = (i + num_subghosts_0_max_diffusivity) +
                            (j + num_subghosts_1_max_diffusivity)*subghostcell_dim_0_max_diffusivity;
                        
                        const int idx_density = (i + num_subghosts_0_density) +
                            (j + num_subghosts_1_density)*subghostcell_dim_0_density;
                        
                        D_max[idx_max_diffusivity] = std::max(mu[idx_max_diffusivity]/rho[idx_density],
                            mu_v[idx_max_diffusivity]/rho[idx_density]);
                        
                        D_max[idx_max_diffusivity] = std::max(D_max[idx_max_diffusivity],
                            kappa[idx_max_diffusivity]/(rho[idx_density]*c_p[idx_max_diffusivity]));
                    }
                }
                
                for (int si = 0; si < d_num_species; si++)
                {
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute the linear index.
                            const int idx_max_diffusivity = (i + num_subghosts_0_max_diffusivity) +
                                (j + num_subghosts_1_max_diffusivity)*subghostcell_dim_0_max_diffusivity;
                            
                            D_max[idx_max_diffusivity] = std::max(D_max[idx_max_diffusivity],
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
                        HAMERS_PRAGMA_SIMD
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
                            
                            D_max[idx_max_diffusivity] = std::max(mu[idx_max_diffusivity]/rho[idx_density],
                                mu_v[idx_max_diffusivity]/rho[idx_density]);
                            
                            D_max[idx_max_diffusivity] = std::max(D_max[idx_max_diffusivity],
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
                            HAMERS_PRAGMA_SIMD
                            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                            {
                                // Compute the linear index.
                                const int idx_max_diffusivity = (i + num_subghosts_0_max_diffusivity) +
                                    (j + num_subghosts_1_max_diffusivity)*subghostcell_dim_0_max_diffusivity +
                                    (k + num_subghosts_2_max_diffusivity)*subghostcell_dim_0_max_diffusivity*
                                        subghostcell_dim_1_max_diffusivity;
                                
                                D_max[idx_max_diffusivity] = std::max(D_max[idx_max_diffusivity],
                                    D[si][idx_max_diffusivity]);
                            }
                        }
                    }
                }
            }
            
            d_cell_data_computed_max_diffusivity = true;
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::"
            << "computeCellDataOfMaxDiffusivityWithDensityMassFractionsPressureAndTemperature()\n"
            << "Cell data of 'MAX_DIFFUSIVITY' is not yet registered."
            << std::endl);
    }
}


/*
 * Compute the cell data of species densities with pressure and temperature in the registered patch.
 */
void
FlowModelFourEqnConservative::computeCellDataOfSpeciesDensitiesWithPressureAndTemperature(
    const hier::Box& domain)
{
    if (d_num_subghosts_species_densities > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_species_densities)
        {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
            for (int si = 0; si < d_num_species; si++)
            {
                TBOX_ASSERT(d_data_species_densities[si]);
            }
#endif
            
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
            if (!domain.empty())
            {
                TBOX_ASSERT(d_subghost_box_species_densities.contains(domain));
            }
#endif
            
            if (!d_cell_data_computed_pressure)
            {
                computeCellDataOfPressureWithDensityMassFractionsAndInternalEnergy(domain);
            }
            
            if (!d_cell_data_computed_temperature)
            {
                computeCellDataOfTemperatureWithDensityMassFractionsAndPressure(domain);
            }
            
            // Compute the density of each species.
            
            for (int si = 0; si < d_num_species; si++)
            {
                std::vector<Real> species_thermo_properties;
                std::vector<Real*> species_thermo_properties_ptr;
                std::vector<const Real*> species_thermo_properties_const_ptr;
                
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
                
                d_equation_of_state_mixing_rules->getEquationOfState(si)->
                    computeDensity(
                        d_data_species_densities[si],
                        d_data_pressure,
                        d_data_temperature,
                        species_thermo_properties_const_ptr,
                        domain);
            }
            
            d_cell_data_computed_species_densities = true;
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::"
            << "computeCellDataOfSpeciesDensitiesWithPressureAndTemperature()\n"
            << "Cell data of 'SPECIES_DENSITIES' is not yet registered."
            << std::endl);
    }
}


/*
 * Compute the cell data of species enthalpies with species densities and pressure in the registered patch.
 */
void
FlowModelFourEqnConservative::computeCellDataOfSpeciesEnthalpiesWithSpeciesDensitiesAndPressure(
    const hier::Box& domain)
{
    if (d_num_subghosts_species_enthalpies > -hier::IntVector::getOne(d_dim))
    {
        if (!d_cell_data_computed_species_enthalpies)
        {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
            for (int si = 0; si < d_num_species; si++)
            {
                TBOX_ASSERT(d_data_species_enthalpies[si]);
            }
#endif
            
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
            if (!domain.empty())
            {
                TBOX_ASSERT(d_subghost_box_species_enthalpies.contains(domain));
            }
#endif
            
            if (!d_cell_data_computed_species_densities)
            {
                computeCellDataOfSpeciesDensitiesWithPressureAndTemperature(domain);
            }
            
            if (!d_cell_data_computed_pressure)
            {
                computeCellDataOfPressureWithDensityMassFractionsAndInternalEnergy(domain);
            }
            
            // Compute the temperature of each species.
            
            for (int si = 0; si < d_num_species; si++)
            {
                std::vector<Real> species_thermo_properties;
                std::vector<Real*> species_thermo_properties_ptr;
                std::vector<const Real*> species_thermo_properties_const_ptr;
                
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
                
                d_equation_of_state_mixing_rules->getEquationOfState(si)->
                    computeEnthalpy(
                        d_data_species_enthalpies[si],
                        d_data_species_densities[si],
                        d_data_pressure,
                        species_thermo_properties_const_ptr,
                        domain);
            }
            
            d_cell_data_computed_species_enthalpies = true;
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFourEqnConservative::"
            << "computeCellDataOfSpeciesEnthalpiesWithSpeciesDensitiesAndPressure()\n"
            << "Cell data of 'SPECIES_ENTHALPIES' is not yet registered."
            << std::endl);
    }
}
