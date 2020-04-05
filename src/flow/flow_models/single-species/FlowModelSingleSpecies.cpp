#include "flow/flow_models/single-species/FlowModelSingleSpecies.hpp"

#include "flow/flow_models/single-species/FlowModelBoundaryUtilitiesSingleSpecies.hpp"
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
        d_num_subghosts_diffusivities(-hier::IntVector::getOne(d_dim)),
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
    if (d_global_derived_cell_data_computed)
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
FlowModelSingleSpecies::registerDerivedVariablesForCharacteristicProjectionOfConservativeVariables(
    const hier::IntVector& num_subghosts,
    const AVERAGING::TYPE& averaging_type)
{
    NULL_USE(num_subghosts);
    
    // Check whether a patch is already registered.
    if (!d_patch)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "registerDerivedVariablesForCharacteristicProjectionOfConservativeVariables()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    // Check whether all or part of derived cell data is already computed.
    if (d_global_derived_cell_data_computed)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "registerDerivedVariablesForCharacteristicProjectionOfConservativeVariables()\n"
            << "Derived cell data is already computed."
            << std::endl);
    }
    
    d_proj_var_conservative_averaging_type = averaging_type;
}


/*
 * Register the required derived variables for transformation between primitive variables
 * and characteristic variables.
 */
void
FlowModelSingleSpecies::registerDerivedVariablesForCharacteristicProjectionOfPrimitiveVariables(
    const hier::IntVector& num_subghosts,
    const AVERAGING::TYPE& averaging_type)
{
    // Check whether a patch is already registered.
    if (!d_patch)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "registerDerivedVariablesForCharacteristicProjectionOfPrimitiveVariables()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    // Check whether all or part of derived cell data is computed.
    if (d_global_derived_cell_data_computed)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "registerDerivedVariablesForCharacteristicProjectionOfPrimitiveVariables()\n"
            << "Derived cell data is already computed."
            << std::endl);
    }
    
    d_proj_var_primitive_averaging_type = averaging_type;
    
    setNumberOfSubGhosts(
        num_subghosts,
        "SOUND_SPEED",
        "PROJECTION_MATRICES");
}


/*
 * Register the required variables for the computation of diffusive fluxes in the registered patch.
 */
void
FlowModelSingleSpecies::registerDiffusiveFluxes(const hier::IntVector& num_subghosts)
{
    // Check whether a patch is already registered.
    if (!d_patch)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "registerDiffusiveFluxes()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    // Check whether all or part of derived cell data is already computed.
    if (d_global_derived_cell_data_computed)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "registerDiffusiveFluxes()\n"
            << "Derived cell data is already computed."
            << std::endl);
    }
    
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
        hier::IntVector::min(d_num_subghosts_velocity, d_num_subghosts_pressure);
    
    d_num_subghosts_diffusivities =
        hier::IntVector::min(d_num_subghosts_diffusivities, d_num_subghosts_temperature);
}


/*
 * Unregister the registered patch. The registered data context and all global derived
 * cell data in the patch are dumped.
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
    
    d_patch = nullptr;
    
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
    d_num_subghosts_diffusivities     = -hier::IntVector::getOne(d_dim);
    
    d_interior_box                   = hier::Box::getEmptyBox(d_dim);
    d_ghost_box                      = hier::Box::getEmptyBox(d_dim);
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
 * Compute the cell data of different registered derived variables with the registered data context.
 */
void
FlowModelSingleSpecies::computeDerivedCellData(const hier::Box& domain)
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
    if (!d_global_derived_cell_data_computed)
    {
        setDerivedCellVariableGhostBoxes();
    }
    
    // Compute the velocity cell data.
    if (d_num_subghosts_velocity > -hier::IntVector::getOne(d_dim))
    {
        if (!d_data_velocity)
        {
            computeCellDataOfVelocity(
                domain);
        }
    }
    
    // Compute the internal energy cell data.
    if (d_num_subghosts_internal_energy > -hier::IntVector::getOne(d_dim))
    {
        if (!d_data_internal_energy)
        {
            computeCellDataOfInternalEnergyWithVelocity(
                domain);
        }
    }
    
    // Compute the pressure cell data.
    if (d_num_subghosts_pressure > -hier::IntVector::getOne(d_dim))
    {
        if (!d_data_pressure)
        {
            computeCellDataOfPressureWithInternalEnergy(
                domain);
        }
    }
    
    // Compute the sound speed cell data.
    if (d_num_subghosts_sound_speed > -hier::IntVector::getOne(d_dim))
    {
        if (!d_data_sound_speed)
        {
            computeCellDataOfSoundSpeedWithPressure(
                domain);
        }
    }
    
    // Compute the temperature cell data.
    if (d_num_subghosts_temperature > -hier::IntVector::getOne(d_dim))
    {
        if (!d_data_temperature)
        {
            computeCellDataOfTemperatureWithPressure(
                domain);
        }
    }
    
    // Compute the x-direction convective flux cell data.
    if (d_num_subghosts_convective_flux_x > -hier::IntVector::getOne(d_dim))
    {
        if (!d_data_convective_flux_x)
        {
            computeCellDataOfConvectiveFluxWithVelocityAndPressure(
                DIRECTION::X_DIRECTION,
                domain);
        }
    }
    
    // Compute the y-direction convective flux cell data.
    if (d_num_subghosts_convective_flux_y > -hier::IntVector::getOne(d_dim))
    {
        if (!d_data_convective_flux_y)
        {
            computeCellDataOfConvectiveFluxWithVelocityAndPressure(
                DIRECTION::Y_DIRECTION,
                domain);
        }
    }
    
    // Compute the z-direction convective flux cell data.
    if (d_num_subghosts_convective_flux_z > -hier::IntVector::getOne(d_dim))
    {
        if (!d_data_convective_flux_z)
        {
            computeCellDataOfConvectiveFluxWithVelocityAndPressure(
                DIRECTION::Z_DIRECTION,
                domain);
        }
    }
    
    // Compute the x-direction maximum wave speed cell data.
    if (d_num_subghosts_max_wave_speed_x > -hier::IntVector::getOne(d_dim))
    {
        if (!d_data_max_wave_speed_x)
        {
            computeCellDataOfMaxWaveSpeedWithVelocityAndSoundSpeed(
                DIRECTION::X_DIRECTION,
                domain);
        }
    }
    
    // Compute the y-direction maximum wave speed cell data.
    if (d_num_subghosts_max_wave_speed_y > -hier::IntVector::getOne(d_dim))
    {
        if (!d_data_max_wave_speed_y)
        {
            computeCellDataOfMaxWaveSpeedWithVelocityAndSoundSpeed(
                DIRECTION::Y_DIRECTION,
                domain);
        }
    }
    
    // Compute the z-direction maximum wave speed cell data.
    if (d_num_subghosts_max_wave_speed_z > -hier::IntVector::getOne(d_dim))
    {
        if (!d_data_max_wave_speed_z)
        {
            computeCellDataOfMaxWaveSpeedWithVelocityAndSoundSpeed(
                DIRECTION::Z_DIRECTION,
                domain);
        }
    }
    
    // Compute the maximum diffusivity cell data.
    if (d_num_subghosts_max_diffusivity > -hier::IntVector::getOne(d_dim))
    {
        if (!d_data_max_diffusivity)
        {
            computeCellDataOfMaxDiffusivityWithPressureAndTemperature(
                domain);
        }
    }
    
    d_global_derived_cell_data_computed = true;
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
        if (!d_data_velocity)
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
        if (!d_data_internal_energy)
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
        if (!d_data_pressure)
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
        if (!d_data_sound_speed)
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
        if (!d_data_temperature)
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
        if (!d_data_convective_flux_x)
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
        if (!d_data_convective_flux_y)
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
        if (!d_data_convective_flux_z)
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
        if (!d_data_max_wave_speed_x)
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
        if (!d_data_max_wave_speed_y)
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
        if (!d_data_max_wave_speed_z)
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
        if (!d_data_max_diffusivity)
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
 * Get the number of projection variables for transformation between conservative
 * variables and characteristic variables.
 */
int
FlowModelSingleSpecies::getNumberOfProjectionVariablesForConservativeVariables() const
{
    return d_dim.getValue() + 5;
}

/*
 * Get the number of projection variables for transformation between primitive variables
 * and characteristic variables.
 */
int
FlowModelSingleSpecies::getNumberOfProjectionVariablesForPrimitiveVariables() const
{
    return 2;
}


/*
 * Compute the side data of the projection variables for transformation between conservative variables and
 * characteristic variables.
 */
void
FlowModelSingleSpecies::computeSideDataOfProjectionVariablesForConservativeVariables(
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
    
    if (static_cast<int>(projection_variables.size()) != d_dim.getValue() + 5)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "computeSideDataOfProjectionVariablesForConservativeVariables()\n"
            << "The number of projection variables are incorrect."
            << std::endl);
    }
    
    /*
     * Check potential failures.
     */
    
    for (int vi = 0; vi < d_dim.getValue() + 5; vi++)
    {
        const hier::IntVector interior_dims_projection_var =
            projection_variables[vi]->getBox().numberCells();
        if (interior_dims_projection_var != d_interior_dims)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::"
                << "computeSideDataOfProjectionVariablesForConservativeVariables()\n"
                << "The interior dimension of the projection variables does not match that of patch."
                << std::endl);
        }
    }
    
    for (int vi = 1; vi < d_dim.getValue() + 5; vi++)
    {
        if (num_ghosts_projection_var != projection_variables[vi]->getGhostCellWidth())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::"
                << "computeSideDataOfProjectionVariablesForConservativeVariables()\n"
                << "The projection variables don't have same ghost cell width."
                << std::endl);
        }
    }
    
    if (num_ghosts_projection_var > d_num_ghosts)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "computeSideDataOfProjectionVariablesForConservativeVariables()\n"
            << "The projection variables have ghost cell width larger than that of density."
            << std::endl);
    }
    
    // Get the cell data of the conservative variables.
    
    boost::shared_ptr<pdat::CellData<double> > data_density =
        getCellDataOfDensity();
    
    boost::shared_ptr<pdat::CellData<double> > data_momentum =
        getCellDataOfMomentum();
    
    boost::shared_ptr<pdat::CellData<double> > data_total_energy =
        getCellDataOfTotalEnergy();
    
    // Get the pointers to the cell data of density and total energy.
    
    double* rho = data_density->getPointer(0);
    double* E   = data_total_energy->getPointer(0);
    
    /*
     * Create temporary side data of averaged conservative variables.
     */
    
    boost::shared_ptr<pdat::SideData<double> > data_density_averaged(
        new pdat::SideData<double>(d_interior_box, 1, num_ghosts_projection_var));
    
    boost::shared_ptr<pdat::SideData<double> > data_momentum_averaged(
        new pdat::SideData<double>(d_interior_box, d_dim.getValue(), num_ghosts_projection_var));
    
    boost::shared_ptr<pdat::SideData<double> > data_total_energy_averaged(
        new pdat::SideData<double>(d_interior_box, 1, num_ghosts_projection_var));
    
    /*
     * Create other temporay side data.
     */
    
    boost::shared_ptr<pdat::SideData<double> > data_internal_energy_averaged(
        new pdat::SideData<double>(d_interior_box, 1, num_ghosts_projection_var));
    
    boost::shared_ptr<pdat::SideData<double> > data_pressure_averaged(
        new pdat::SideData<double>(d_interior_box, 1, num_ghosts_projection_var));
    
    /*
     * Declare pointers to the side data of averaged density and total energy.
     */
    
    double* rho_average     = nullptr;
    double* E_average       = nullptr;
    double* e_average       = nullptr;
    double* epsilon_average = nullptr;
    double* p_average       = nullptr;
    double* H_average       = nullptr;
    
    // Get the thermodynamic properties of the species.
    std::vector<const double*> thermo_properties_ptr;
    thermo_properties_ptr.reserve(static_cast<int> (d_thermo_properties.size()));
    for (int ti = 0; ti < static_cast<int> (d_thermo_properties.size()); ti++)
    {
        thermo_properties_ptr.push_back(&d_thermo_properties[ti]);
    }
    
    if (d_dim == tbox::Dimension(1))
    {
        const int interior_dim_0 = d_interior_dims[0];
        
        const int num_ghosts_0 = d_num_ghosts[0];
        const int num_ghosts_0_projection_var = num_ghosts_projection_var[0];
        
        // Get the pointer to the cell data of momentum.
        
        double* rho_u = data_momentum->getPointer(0);
        
        // Declare pointers to the side data of averaged momentum and velocity.
        
        double* rho_u_average = nullptr;
        double* u_average     = nullptr;
        
        switch (d_proj_var_conservative_averaging_type)
        {
            case AVERAGING::SIMPLE:
            {
                /*
                 * Compute the averaged conservative variables in the x-direction.
                 */
                
                rho_average   = data_density_averaged->getPointer(0);
                rho_u_average = data_momentum_averaged->getPointer(0, 0);
                E_average     = data_total_energy_averaged->getPointer(0);
                
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
                    
                    rho_average[idx_face_x]   = double(1)/double(2)*(rho[idx_L] + rho[idx_R]);
                    rho_u_average[idx_face_x] = double(1)/double(2)*(rho_u[idx_L] + rho_u[idx_R]);
                    E_average[idx_face_x]     = double(1)/double(2)*(E[idx_L] + E[idx_R]);
                }
                
                /*
                 * Compute the projection variables in the x-direction.
                 */
                
                u_average       = projection_variables[0]->getPointer(0);
                e_average       = projection_variables[1]->getPointer(0);
                epsilon_average = data_internal_energy_averaged->getPointer(0);
                p_average       = data_pressure_averaged->getPointer(0);
                H_average       = projection_variables[2]->getPointer(0);
                
                // Compute the velocity and internal energy.
                
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = -num_ghosts_0_projection_var;
                     i < interior_dim_0 + 1 + num_ghosts_0_projection_var;
                     i++)
                {
                    // Compute the linear indices.
                    const int idx_face_x = i + num_ghosts_0_projection_var;
                    
                    u_average[idx_face_x] = rho_u_average[idx_face_x]/rho_average[idx_face_x];
                    
                    e_average[idx_face_x] = E_average[idx_face_x]/rho_average[idx_face_x];
                    
                    epsilon_average[idx_face_x] = e_average[idx_face_x] -
                        double(1)/double(2)*u_average[idx_face_x]*u_average[idx_face_x];
                }
                
                // Compute the presure, sound speed, partial pressure partial density and Gruneisen parameter.
                
                d_equation_of_state_mixing_rules->getEquationOfState()->computePressure(
                    data_pressure_averaged,
                    data_density_averaged,
                    data_internal_energy_averaged,
                    thermo_properties_ptr,
                    0);
                
                d_equation_of_state_mixing_rules->getEquationOfState()->computeSoundSpeed(
                    projection_variables[3],
                    data_density_averaged,
                    data_pressure_averaged,
                    thermo_properties_ptr,
                    0);
                
                d_equation_of_state_mixing_rules->getEquationOfState()->computePressureDerivativeWithDensity(
                    projection_variables[4],
                    data_density_averaged,
                    data_pressure_averaged,
                    thermo_properties_ptr,
                    0);
                
                d_equation_of_state_mixing_rules->getEquationOfState()->computeGruneisenParameter(
                    projection_variables[5],
                    data_density_averaged,
                    data_pressure_averaged,
                    thermo_properties_ptr,
                    0);
                
                // Compute the total specific enthalpy.
                
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = -num_ghosts_0_projection_var;
                     i < interior_dim_0 + 1 + num_ghosts_0_projection_var;
                     i++)
                {
                    // Compute the linear indices.
                    const int idx_face_x = i + num_ghosts_0_projection_var;
                    
                    H_average[idx_face_x] = e_average[idx_face_x] +
                        p_average[idx_face_x]/rho_average[idx_face_x];
                }
                
                break;
            }
            case AVERAGING::ROE:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::"
                    << "computeSideDataOfProjectionVariablesForConservativeVariables()\n"
                    << "Roe averaging is not yet implemented."
                    << std::endl);
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::"
                    << "computeSideDataOfProjectionVariablesForConservativeVariables()\n"
                    << "Unknown d_proj_var_conservative_averaging_type given."
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
        
        // Get the pointers to the cell data of momentum.
        
        double* rho_u = data_momentum->getPointer(0);
        double* rho_v = data_momentum->getPointer(1);
        
        // Declare pointers to the side data of averaged momentum and velocity.
        
        double* rho_u_average = nullptr;
        double* rho_v_average = nullptr;
        double* u_average     = nullptr;
        double* v_average     = nullptr;
        
        switch (d_proj_var_conservative_averaging_type)
        {
            case AVERAGING::SIMPLE:
            {
                /*
                 * Compute the averaged conservative variables in the x-direction.
                 */
                
                rho_average   = data_density_averaged->getPointer(0);
                rho_u_average = data_momentum_averaged->getPointer(0, 0);
                rho_v_average = data_momentum_averaged->getPointer(0, 1);
                E_average     = data_total_energy_averaged->getPointer(0);
                
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
                        
                        rho_average[idx_face_x]   = double(1)/double(2)*(rho[idx_L] + rho[idx_R]);
                        rho_u_average[idx_face_x] = double(1)/double(2)*(rho_u[idx_L] + rho_u[idx_R]);
                        rho_v_average[idx_face_x] = double(1)/double(2)*(rho_v[idx_L] + rho_v[idx_R]);
                        E_average[idx_face_x]     = double(1)/double(2)*(E[idx_L] + E[idx_R]);
                    }
                }
                
                /*
                 * Compute the projection variables in the x-direction.
                 */
                
                u_average       = projection_variables[0]->getPointer(0);
                v_average       = projection_variables[1]->getPointer(0);
                e_average       = projection_variables[2]->getPointer(0);
                epsilon_average = data_internal_energy_averaged->getPointer(0);
                p_average       = data_pressure_averaged->getPointer(0);
                H_average       = projection_variables[3]->getPointer(0);
                
                // Compute the velocity and internal energy.
                
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
                        
                        u_average[idx_face_x] = rho_u_average[idx_face_x]/rho_average[idx_face_x];
                        v_average[idx_face_x] = rho_v_average[idx_face_x]/rho_average[idx_face_x];
                        
                        e_average[idx_face_x] = E_average[idx_face_x]/rho_average[idx_face_x];
                        
                        epsilon_average[idx_face_x] = e_average[idx_face_x] -
                            double(1)/double(2)*(u_average[idx_face_x]*u_average[idx_face_x] +
                                v_average[idx_face_x]*v_average[idx_face_x]);
                    }
                }
                
                // Compute the presure, sound speed, partial pressure partial density and Gruneisen parameter.
                
                d_equation_of_state_mixing_rules->getEquationOfState()->computePressure(
                    data_pressure_averaged,
                    data_density_averaged,
                    data_internal_energy_averaged,
                    thermo_properties_ptr,
                    0);
                
                d_equation_of_state_mixing_rules->getEquationOfState()->computeSoundSpeed(
                    projection_variables[4],
                    data_density_averaged,
                    data_pressure_averaged,
                    thermo_properties_ptr,
                    0);
                
                d_equation_of_state_mixing_rules->getEquationOfState()->computePressureDerivativeWithDensity(
                    projection_variables[5],
                    data_density_averaged,
                    data_pressure_averaged,
                    thermo_properties_ptr,
                    0);
                
                d_equation_of_state_mixing_rules->getEquationOfState()->computeGruneisenParameter(
                    projection_variables[6],
                    data_density_averaged,
                    data_pressure_averaged,
                    thermo_properties_ptr,
                    0);
                
                // Compute the total specific enthalpy.
                
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
                        
                        H_average[idx_face_x] = e_average[idx_face_x] +
                            p_average[idx_face_x]/rho_average[idx_face_x];
                    }
                }
                
                /*
                 * Compute the averaged conservative variables in the y-direction.
                 */
                
                rho_average   = data_density_averaged->getPointer(1);
                rho_u_average = data_momentum_averaged->getPointer(1, 0);
                rho_v_average = data_momentum_averaged->getPointer(1, 1);
                E_average     = data_total_energy_averaged->getPointer(1);
                
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
                        
                        rho_average[idx_face_y]   = double(1)/double(2)*(rho[idx_B] + rho[idx_T]);
                        rho_u_average[idx_face_y] = double(1)/double(2)*(rho_u[idx_B] + rho_u[idx_T]);
                        rho_v_average[idx_face_y] = double(1)/double(2)*(rho_v[idx_B] + rho_v[idx_T]);
                        E_average[idx_face_y]     = double(1)/double(2)*(E[idx_B] + E[idx_T]);
                    }
                }
                
                /*
                 * Compute the projection variables in the y-direction.
                 */
                
                u_average       = projection_variables[0]->getPointer(1);
                v_average       = projection_variables[1]->getPointer(1);
                e_average       = projection_variables[2]->getPointer(1);
                epsilon_average = data_internal_energy_averaged->getPointer(1);
                p_average       = data_pressure_averaged->getPointer(1);
                H_average       = projection_variables[3]->getPointer(1);
                
                // Compute the velocity and internal energy.
                
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
                        
                        u_average[idx_face_y] = rho_u_average[idx_face_y]/rho_average[idx_face_y];
                        v_average[idx_face_y] = rho_v_average[idx_face_y]/rho_average[idx_face_y];
                        
                        e_average[idx_face_y] = E_average[idx_face_y]/rho_average[idx_face_y];
                        
                        epsilon_average[idx_face_y] = e_average[idx_face_y] -
                            double(1)/double(2)*(u_average[idx_face_y]*u_average[idx_face_y] +
                                v_average[idx_face_y]*v_average[idx_face_y]);
                    }
                }
                
                // Compute the presure, sound speed, partial pressure partial density and Gruneisen parameter.
                
                d_equation_of_state_mixing_rules->getEquationOfState()->computePressure(
                    data_pressure_averaged,
                    data_density_averaged,
                    data_internal_energy_averaged,
                    thermo_properties_ptr,
                    1);
                
                d_equation_of_state_mixing_rules->getEquationOfState()->computeSoundSpeed(
                    projection_variables[4],
                    data_density_averaged,
                    data_pressure_averaged,
                    thermo_properties_ptr,
                    1);
                
                d_equation_of_state_mixing_rules->getEquationOfState()->computePressureDerivativeWithDensity(
                    projection_variables[5],
                    data_density_averaged,
                    data_pressure_averaged,
                    thermo_properties_ptr,
                    1);
                
                d_equation_of_state_mixing_rules->getEquationOfState()->computeGruneisenParameter(
                    projection_variables[6],
                    data_density_averaged,
                    data_pressure_averaged,
                    thermo_properties_ptr,
                    1);
                
                // Compute the total specific enthalpy.
                
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
                        
                        H_average[idx_face_y] = e_average[idx_face_y] +
                            p_average[idx_face_y]/rho_average[idx_face_y];
                    }
                }
                
                break;
            }
            case AVERAGING::ROE:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::"
                    << "computeSideDataOfProjectionVariablesForConservativeVariables()\n"
                    << "Roe averaging is not yet implemented."
                    << std::endl);
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::"
                    << "computeSideDataOfProjectionVariablesForConservativeVariables()\n"
                    << "Unknown d_proj_var_conservative_averaging_type given."
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
        
        // Get the pointers to the cell data of momentum.
        
        double* rho_u = data_momentum->getPointer(0);
        double* rho_v = data_momentum->getPointer(1);
        double* rho_w = data_momentum->getPointer(2);
        
        // Declare pointers to the side data of averaged momentum and velocity.
        
        double* rho_u_average = nullptr;
        double* rho_v_average = nullptr;
        double* rho_w_average = nullptr;
        double* u_average     = nullptr;
        double* v_average     = nullptr;
        double* w_average     = nullptr;
        
        switch (d_proj_var_conservative_averaging_type)
        {
            case AVERAGING::SIMPLE:
            {
                /*
                 * Compute the averaged conservative variables in the x-direction.
                 */
                
                rho_average   = data_density_averaged->getPointer(0);
                rho_u_average = data_momentum_averaged->getPointer(0, 0);
                rho_v_average = data_momentum_averaged->getPointer(0, 1);
                rho_w_average = data_momentum_averaged->getPointer(0, 2);
                E_average     = data_total_energy_averaged->getPointer(0);
                
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
                            
                            rho_average[idx_face_x]   = double(1)/double(2)*(rho[idx_L] + rho[idx_R]);
                            rho_u_average[idx_face_x] = double(1)/double(2)*(rho_u[idx_L] + rho_u[idx_R]);
                            rho_v_average[idx_face_x] = double(1)/double(2)*(rho_v[idx_L] + rho_v[idx_R]);
                            rho_w_average[idx_face_x] = double(1)/double(2)*(rho_w[idx_L] + rho_w[idx_R]);
                            E_average[idx_face_x]     = double(1)/double(2)*(E[idx_L] + E[idx_R]);
                        }
                    }
                }
                
                /*
                 * Compute the projection variables in the x-direction.
                 */
                
                u_average       = projection_variables[0]->getPointer(0);
                v_average       = projection_variables[1]->getPointer(0);
                w_average       = projection_variables[2]->getPointer(0);
                e_average       = projection_variables[3]->getPointer(0);
                epsilon_average = data_internal_energy_averaged->getPointer(0);
                p_average       = data_pressure_averaged->getPointer(0);
                H_average       = projection_variables[4]->getPointer(0);
                
                // Compute the velocity and internal energy.
                
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
                            
                            u_average[idx_face_x] = rho_u_average[idx_face_x]/rho_average[idx_face_x];
                            v_average[idx_face_x] = rho_v_average[idx_face_x]/rho_average[idx_face_x];
                            w_average[idx_face_x] = rho_w_average[idx_face_x]/rho_average[idx_face_x];
                            
                            e_average[idx_face_x] = E_average[idx_face_x]/rho_average[idx_face_x];
                            
                            epsilon_average[idx_face_x] = e_average[idx_face_x] -
                                double(1)/double(2)*(u_average[idx_face_x]*u_average[idx_face_x] +
                                    v_average[idx_face_x]*v_average[idx_face_x] +
                                    w_average[idx_face_x]*w_average[idx_face_x]);
                        }
                    }
                }
                
                // Compute the presure, sound speed, partial pressure partial density and Gruneisen parameter.
                
                d_equation_of_state_mixing_rules->getEquationOfState()->computePressure(
                    data_pressure_averaged,
                    data_density_averaged,
                    data_internal_energy_averaged,
                    thermo_properties_ptr,
                    0);
                
                d_equation_of_state_mixing_rules->getEquationOfState()->computeSoundSpeed(
                    projection_variables[5],
                    data_density_averaged,
                    data_pressure_averaged,
                    thermo_properties_ptr,
                    0);
                
                d_equation_of_state_mixing_rules->getEquationOfState()->computePressureDerivativeWithDensity(
                    projection_variables[6],
                    data_density_averaged,
                    data_pressure_averaged,
                    thermo_properties_ptr,
                    0);
                
                d_equation_of_state_mixing_rules->getEquationOfState()->computeGruneisenParameter(
                    projection_variables[7],
                    data_density_averaged,
                    data_pressure_averaged,
                    thermo_properties_ptr,
                    0);
                
                // Compute the total specific enthalpy.
                
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
                            
                            H_average[idx_face_x] = e_average[idx_face_x] +
                                p_average[idx_face_x]/rho_average[idx_face_x];
                        }
                    }
                }
                
                /*
                 * Compute the averaged conservative variables in the y-direction.
                 */
                
                rho_average   = data_density_averaged->getPointer(1);
                rho_u_average = data_momentum_averaged->getPointer(1, 0);
                rho_v_average = data_momentum_averaged->getPointer(1, 1);
                rho_w_average = data_momentum_averaged->getPointer(1, 2);
                E_average     = data_total_energy_averaged->getPointer(1);
                
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
                            
                            rho_average[idx_face_y]   = double(1)/double(2)*(rho[idx_B] + rho[idx_T]);
                            rho_u_average[idx_face_y] = double(1)/double(2)*(rho_u[idx_B] + rho_u[idx_T]);
                            rho_v_average[idx_face_y] = double(1)/double(2)*(rho_v[idx_B] + rho_v[idx_T]);
                            rho_w_average[idx_face_y] = double(1)/double(2)*(rho_w[idx_B] + rho_w[idx_T]);
                            E_average[idx_face_y]     = double(1)/double(2)*(E[idx_B] + E[idx_T]);
                        }
                    }
                }
                
                /*
                 * Compute the projection variables in the y-direction.
                 */
                
                u_average       = projection_variables[0]->getPointer(1);
                v_average       = projection_variables[1]->getPointer(1);
                w_average       = projection_variables[2]->getPointer(1);
                e_average       = projection_variables[3]->getPointer(1);
                epsilon_average = data_internal_energy_averaged->getPointer(1);
                p_average       = data_pressure_averaged->getPointer(1);
                H_average       = projection_variables[4]->getPointer(1);
                
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
                            
                            u_average[idx_face_y] = rho_u_average[idx_face_y]/rho_average[idx_face_y];
                            v_average[idx_face_y] = rho_v_average[idx_face_y]/rho_average[idx_face_y];
                            w_average[idx_face_y] = rho_w_average[idx_face_y]/rho_average[idx_face_y];
                            
                            e_average[idx_face_y] = E_average[idx_face_y]/rho_average[idx_face_y];
                            
                            epsilon_average[idx_face_y] = e_average[idx_face_y] -
                                double(1)/double(2)*(u_average[idx_face_y]*u_average[idx_face_y] +
                                    v_average[idx_face_y]*v_average[idx_face_y] +
                                    w_average[idx_face_y]*w_average[idx_face_y]);
                        }
                    }
                }
                
                // Compute the presure, sound speed, partial pressure partial density and Gruneisen parameter.
                
                d_equation_of_state_mixing_rules->getEquationOfState()->computePressure(
                    data_pressure_averaged,
                    data_density_averaged,
                    data_internal_energy_averaged,
                    thermo_properties_ptr,
                    1);
                
                d_equation_of_state_mixing_rules->getEquationOfState()->computeSoundSpeed(
                    projection_variables[5],
                    data_density_averaged,
                    data_pressure_averaged,
                    thermo_properties_ptr,
                    1);
                
                d_equation_of_state_mixing_rules->getEquationOfState()->computePressureDerivativeWithDensity(
                    projection_variables[6],
                    data_density_averaged,
                    data_pressure_averaged,
                    thermo_properties_ptr,
                    1);
                
                d_equation_of_state_mixing_rules->getEquationOfState()->computeGruneisenParameter(
                    projection_variables[7],
                    data_density_averaged,
                    data_pressure_averaged,
                    thermo_properties_ptr,
                    1);
                
                // Compute the total specific enthalpy.
                
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
                            
                            H_average[idx_face_y] = e_average[idx_face_y] +
                                p_average[idx_face_y]/rho_average[idx_face_y];
                        }
                    }
                }
                
                /*
                 * Compute the averaged conservative variables in the z-direction.
                 */
                
                rho_average   = data_density_averaged->getPointer(2);
                rho_u_average = data_momentum_averaged->getPointer(2, 0);
                rho_v_average = data_momentum_averaged->getPointer(2, 1);
                rho_w_average = data_momentum_averaged->getPointer(2, 2);
                E_average     = data_total_energy_averaged->getPointer(2);
                
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
                            
                            rho_average[idx_face_z]   = double(1)/double(2)*(rho[idx_B] + rho[idx_F]);
                            rho_u_average[idx_face_z] = double(1)/double(2)*(rho_u[idx_B] + rho_u[idx_F]);
                            rho_v_average[idx_face_z] = double(1)/double(2)*(rho_v[idx_B] + rho_v[idx_F]);
                            rho_w_average[idx_face_z] = double(1)/double(2)*(rho_w[idx_B] + rho_w[idx_F]);
                            E_average[idx_face_z]     = double(1)/double(2)*(E[idx_B] + E[idx_F]);
                        }
                    }
                }
                
                /*
                 * Compute the projection variables in the z-direction.
                 */
                
                u_average       = projection_variables[0]->getPointer(2);
                v_average       = projection_variables[1]->getPointer(2);
                w_average       = projection_variables[2]->getPointer(2);
                e_average       = projection_variables[3]->getPointer(2);
                epsilon_average = data_internal_energy_averaged->getPointer(2);
                p_average       = data_pressure_averaged->getPointer(2);
                H_average       = projection_variables[4]->getPointer(2);
                
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
                            
                            u_average[idx_face_z] = rho_u_average[idx_face_z]/rho_average[idx_face_z];
                            v_average[idx_face_z] = rho_v_average[idx_face_z]/rho_average[idx_face_z];
                            w_average[idx_face_z] = rho_w_average[idx_face_z]/rho_average[idx_face_z];
                            
                            e_average[idx_face_z] = E_average[idx_face_z]/rho_average[idx_face_z];
                            
                            epsilon_average[idx_face_z] = e_average[idx_face_z] -
                                double(1)/double(2)*(u_average[idx_face_z]*u_average[idx_face_z] +
                                    v_average[idx_face_z]*v_average[idx_face_z] +
                                    w_average[idx_face_z]*w_average[idx_face_z]);
                        }
                    }
                }
                
                // Compute the presure, sound speed, partial pressure partial density and Gruneisen parameter.
                
                d_equation_of_state_mixing_rules->getEquationOfState()->computePressure(
                    data_pressure_averaged,
                    data_density_averaged,
                    data_internal_energy_averaged,
                    thermo_properties_ptr,
                    2);
                
                d_equation_of_state_mixing_rules->getEquationOfState()->computeSoundSpeed(
                    projection_variables[5],
                    data_density_averaged,
                    data_pressure_averaged,
                    thermo_properties_ptr,
                    2);
                
                d_equation_of_state_mixing_rules->getEquationOfState()->computePressureDerivativeWithDensity(
                    projection_variables[6],
                    data_density_averaged,
                    data_pressure_averaged,
                    thermo_properties_ptr,
                    2);
                
                d_equation_of_state_mixing_rules->getEquationOfState()->computeGruneisenParameter(
                    projection_variables[7],
                    data_density_averaged,
                    data_pressure_averaged,
                    thermo_properties_ptr,
                    2);
                
                // Compute the total specific enthalpy.
                
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
                            
                            H_average[idx_face_z] = e_average[idx_face_z] +
                                p_average[idx_face_z]/rho_average[idx_face_z];
                        }
                    }
                }
                
                break;
            }
            case AVERAGING::ROE:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::"
                    << "computeSideDataOfProjectionVariablesForConservativeVariables()\n"
                    << "Roe averaging is not yet implemented."
                    << std::endl);
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::"
                    << "computeSideDataOfProjectionVariablesForConservativeVariables()\n"
                    << "Unknown d_proj_var_conservative_averaging_type given."
                    << std::endl);
            }
        }
    }
}


/*
 * Compute the side data of the projection variables for transformation between primitive variables and characteristic
 * variables.
 */
void
FlowModelSingleSpecies::computeSideDataOfProjectionVariablesForPrimitiveVariables(
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
    
    if (static_cast<int>(projection_variables.size()) != 2)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "computeSideDataOfProjectionVariablesForPrimitiveVariables()\n"
            << "There should be two projection variables."
            << std::endl);
    }
    
    /*
     * Check potential failures.
     */
    
    for (int vi = 0; vi < 2; vi++)
    {
        const hier::IntVector interior_dims_projection_var =
            projection_variables[vi]->getBox().numberCells();
        if (interior_dims_projection_var != d_interior_dims)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::"
                << "computeSideDataOfProjectionVariablesForPrimitiveVariables()\n"
                << "The interior dimension of the projection variables does not match that of patch."
                << std::endl);
        }
    }
    
    if (num_ghosts_projection_var != projection_variables[1]->getGhostCellWidth())
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "computeSideDataOfProjectionVariablesForPrimitiveVariables()\n"
            << "The projection variables don't have same ghost cell width."
            << std::endl);
    }
    
    if (num_ghosts_projection_var > d_num_ghosts)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "computeSideDataOfProjectionVariablesForPrimitiveVariables()\n"
            << "The projection variables have ghost cell width larger than that of density."
            << std::endl);
    }
    
    if (num_ghosts_projection_var > d_num_subghosts_sound_speed)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "computeSideDataOfProjectionVariablesForPrimitiveVariables()\n"
            << "The projection variables have ghost cell width larger than that of sound speed."
            << std::endl);
    }
    
    // Get the cell data of the variable density.
    boost::shared_ptr<pdat::CellData<double> > data_density =
        getCellDataOfDensity();
    
    // Get the pointers to the cell data of density and sound speed.
    double* rho = data_density->getPointer(0);
    if (!d_data_sound_speed)
    {
        computeCellDataOfSoundSpeedWithPressure(empty_box);
    }
    double* c = d_data_sound_speed->getPointer(0);
    
    /*
     * Declare pointers to different data.
     */
    
    double* rho_average = nullptr;
    double* c_average = nullptr;
    
    if (d_dim == tbox::Dimension(1))
    {
        const int interior_dim_0 = d_interior_dims[0];
        
        const int num_ghosts_0 = d_num_ghosts[0];
        const int num_ghosts_0_projection_var = num_ghosts_projection_var[0];
        const int num_subghosts_0_sound_speed = d_num_subghosts_sound_speed[0];
        
        switch (d_proj_var_primitive_averaging_type)
        {
            case AVERAGING::SIMPLE:
            {
                /*
                 * Compute the projection variables in the x-direction.
                 */
                
                rho_average = projection_variables[0]->getPointer(0);
                c_average = projection_variables[1]->getPointer(0);
                
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
                    const int idx_sound_speed_L = i - 1 + num_subghosts_0_sound_speed;
                    const int idx_sound_speed_R = i + num_subghosts_0_sound_speed;
                    
                    rho_average[idx_face_x] = double(1)/double(2)*(rho[idx_L] + rho[idx_R]);
                    c_average[idx_face_x] = double(1)/double(2)*(c[idx_sound_speed_L] + c[idx_sound_speed_R]);
                }
                
                break;
            }
            case AVERAGING::ROE:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::"
                    << "computeSideDataOfProjectionVariablesForPrimitiveVariables()\n"
                    << "Roe averaging is not yet implemented."
                    << std::endl);
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::"
                    << "computeSideDataOfProjectionVariablesForPrimitiveVariables()\n"
                    << "Unknown d_proj_var_primitive_averaging_type given."
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
        
        const int num_subghosts_0_sound_speed = d_num_subghosts_sound_speed[0];
        const int num_subghosts_1_sound_speed = d_num_subghosts_sound_speed[1];
        const int subghostcell_dim_0_sound_speed = d_subghostcell_dims_sound_speed[0];
        
        switch (d_proj_var_primitive_averaging_type)
        {
            case AVERAGING::SIMPLE:
            {
                /*
                 * Compute the projection variables in the x-direction.
                 */
                
                rho_average = projection_variables[0]->getPointer(0);
                c_average = projection_variables[1]->getPointer(0);
                
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
                        
                        const int idx_sound_speed_L = (i - 1 + num_subghosts_0_sound_speed) +
                            (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed;
                        
                        const int idx_sound_speed_R = (i + num_subghosts_0_sound_speed) +
                            (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed;
                        
                        rho_average[idx_face_x] = double(1)/double(2)*(rho[idx_L] + rho[idx_R]);
                        c_average[idx_face_x] = double(1)/double(2)*(c[idx_sound_speed_L] + c[idx_sound_speed_R]);
                    }
                }
                
                /*
                 * Compute the projection variables in the y-direction.
                 */
                rho_average = projection_variables[0]->getPointer(1);
                c_average = projection_variables[1]->getPointer(1);
                
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
                        
                        const int idx_sound_speed_B = (i + num_subghosts_0_sound_speed) +
                            (j - 1 + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed;
                        
                        const int idx_sound_speed_T = (i + num_subghosts_0_sound_speed) +
                            (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed;
                        
                        rho_average[idx_face_y] = double(1)/double(2)*(rho[idx_B] + rho[idx_T]);
                        c_average[idx_face_y] = double(1)/double(2)*(c[idx_sound_speed_B] + c[idx_sound_speed_T]);
                    }
                }
                
                break;
            }
            case AVERAGING::ROE:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::"
                    << "computeSideDataOfProjectionVariablesForPrimitiveVariables()\n"
                    << "Roe averaging is not yet implemented."
                    << std::endl);
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::"
                    << "computeSideDataOfProjectionVariablesForPrimitiveVariables()\n"
                    << "Unknown d_proj_var_primitive_averaging_type given."
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
        
        const int num_subghosts_0_sound_speed = d_num_subghosts_sound_speed[0];
        const int num_subghosts_1_sound_speed = d_num_subghosts_sound_speed[1];
        const int num_subghosts_2_sound_speed = d_num_subghosts_sound_speed[2];
        const int subghostcell_dim_0_sound_speed = d_subghostcell_dims_sound_speed[0];
        const int subghostcell_dim_1_sound_speed = d_subghostcell_dims_sound_speed[1];
        
        switch (d_proj_var_primitive_averaging_type)
        {
            case AVERAGING::SIMPLE:
            {
                /*
                 * Compute the projection variables in the x-direction.
                 */
                
                rho_average = projection_variables[0]->getPointer(0);
                c_average = projection_variables[1]->getPointer(0);
                
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
                            
                            const int idx_sound_speed_L = (i - 1 + num_subghosts_0_sound_speed) +
                                (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed +
                                (k + num_subghosts_2_sound_speed)*subghostcell_dim_0_sound_speed*
                                    subghostcell_dim_1_sound_speed;
                            
                            const int idx_sound_speed_R = (i + num_subghosts_0_sound_speed) +
                                (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed +
                                (k + num_subghosts_2_sound_speed)*subghostcell_dim_0_sound_speed*
                                    subghostcell_dim_1_sound_speed;
                            
                            rho_average[idx_face_x] = double(1)/double(2)*(rho[idx_L] + rho[idx_R]);
                            c_average[idx_face_x] = double(1)/double(2)*(c[idx_sound_speed_L] + c[idx_sound_speed_R]);
                        }
                    }
                }
                
                /*
                 * Compute the projection variables in the y-direction.
                 */
                
                rho_average = projection_variables[0]->getPointer(1);
                c_average = projection_variables[1]->getPointer(1);
                
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
                            
                            const int idx_sound_speed_B = (i + num_subghosts_0_sound_speed) +
                                (j - 1 + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed +
                                (k + num_subghosts_2_sound_speed)*subghostcell_dim_0_sound_speed*
                                    subghostcell_dim_1_sound_speed;
                            
                            const int idx_sound_speed_T = (i + num_subghosts_0_sound_speed) +
                                (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed +
                                (k + num_subghosts_2_sound_speed)*subghostcell_dim_0_sound_speed*
                                    subghostcell_dim_1_sound_speed;
                            
                            rho_average[idx_face_y] = double(1)/double(2)*(rho[idx_B] + rho[idx_T]);
                            c_average[idx_face_y] = double(1)/double(2)*(c[idx_sound_speed_B] + c[idx_sound_speed_T]);
                        }
                    }
                }
                
                /*
                 * Compute the projection variables in the z-direction.
                 */
                
                rho_average = projection_variables[0]->getPointer(2);
                c_average = projection_variables[1]->getPointer(2);
                
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
                            
                            const int idx_sound_speed_B = (i + num_subghosts_0_sound_speed) +
                                (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed +
                                (k - 1 + num_subghosts_2_sound_speed)*subghostcell_dim_0_sound_speed*
                                    subghostcell_dim_1_sound_speed;
                            
                            const int idx_sound_speed_F = (i + num_subghosts_0_sound_speed) +
                                (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed +
                                (k + num_subghosts_2_sound_speed)*subghostcell_dim_0_sound_speed*
                                    subghostcell_dim_1_sound_speed;
                            
                            rho_average[idx_face_z] = double(1)/double(2)*(rho[idx_B] + rho[idx_F]);
                            c_average[idx_face_z] = double(1)/double(2)*(c[idx_sound_speed_B] + c[idx_sound_speed_F]);
                        }
                    }
                }
                
                break;
            }
            case AVERAGING::ROE:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::"
                    << "computeSideDataOfProjectionVariablesForPrimitiveVariables()\n"
                    << "Roe averaging is not yet implemented."
                    << std::endl);
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::"
                    << "computeSideDataOfProjectionVariablesForPrimitiveVariables()\n"
                    << "Unknown d_proj_var_primitive_averaging_type given."
                    << std::endl);
            }
        }
    }
}


/*
 * Compute the side data of characteristic variables from conservative variables.
 */
void
FlowModelSingleSpecies::computeSideDataOfCharacteristicVariablesFromConservativeVariables(
    std::vector<boost::shared_ptr<pdat::SideData<double> > >& characteristic_variables,
    const std::vector<boost::shared_ptr<pdat::CellData<double> > >& conservative_variables,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& projection_variables,
    const int& idx_offset)
{
    /*
     * Get the numbers of ghost cells of the variables.
     */
    
    const hier::IntVector num_ghosts_characteristic_var = characteristic_variables[0]->
        getGhostCellWidth();
    
    std::vector<hier::IntVector> num_ghosts_conservative_var;
    num_ghosts_conservative_var.reserve(static_cast<int>(conservative_variables.size()));
    for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
    {
        num_ghosts_conservative_var.push_back(conservative_variables[vi]->
            getGhostCellWidth());
    }
    
    const hier::IntVector num_ghosts_projection_var = projection_variables[0]->getGhostCellWidth();
    
    /*
     * Get the ghost cell dimensions of characteristic and conservative variables.
     */
    
    const hier::IntVector ghostcell_dims_characteristic_var = characteristic_variables[0]->
        getGhostBox().numberCells();
    
    std::vector<hier::IntVector> ghostcell_dims_conservative_var;
    ghostcell_dims_conservative_var.reserve(static_cast<int>(conservative_variables.size()));
    for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
    {
        ghostcell_dims_conservative_var.push_back(conservative_variables[vi]->
            getGhostBox().numberCells());
    }
    
    /*
     * Check the size of variables.
     */
    
    if (static_cast<int>(characteristic_variables.size()) != d_num_eqn)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "computeSideDataOfCharacteristicVariablesFromConservativeVariables()\n"
            << "The number of characteristic variables are incorrect."
            << std::endl);
    }
    if (static_cast<int>(conservative_variables.size()) != 3)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "computeSideDataOfCharacteristicVariablesFromConservativeVariables()\n"
            << "The number of conservative variables are incorrect."
            << std::endl);
    }
    if (conservative_variables[0]->getDepth() != 1 ||
        conservative_variables[1]->getDepth() != d_dim.getValue() ||
        conservative_variables[2]->getDepth() != 1)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "computeSideDataOfCharacteristicVariablesFromConservativeVariables()\n"
            << "The depths of one or more conservative variables are incorrect."
            << std::endl);
    }
    if (static_cast<int>(projection_variables.size()) != d_dim.getValue() + 5)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "computeSideDataOfCharacteristicVariablesFromConservativeVariables()\n"
            << "The number of projection variables are incorrect."
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
                << ": FlowModelSingleSpecies::"
                << "computeSideDataOfCharacteristicVariablesFromConservativeVariables()\n"
                << "The interior dimension of the characteristic variables does not match that of patch."
                << std::endl);
        }
    }
    for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
    {
        const hier::IntVector interior_dims_conservative_var =
            conservative_variables[vi]->getBox().numberCells();
        
        if (interior_dims_conservative_var != d_interior_dims)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::"
                << "computeSideDataOfCharacteristicVariablesFromConservativeVariables()\n"
                << "The interior dimension of the conservative variables does not match that of patch."
                << std::endl);
        }
    }
    for (int vi = 0; vi < d_dim.getValue() + 5; vi++)
    {
        const hier::IntVector interior_dims_projection_var = projection_variables[vi]->getBox().numberCells();
        if (interior_dims_projection_var != d_interior_dims)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::"
                << "computeSideDataOfCharacteristicVariablesFromConservativeVariables()\n"
                << "The interior dimension of the projection variables does not match that of patch."
                << std::endl);
        }
    }
    
    for (int ei = 1; ei < d_num_eqn; ei++)
    {
        if (num_ghosts_characteristic_var != characteristic_variables[ei]->getGhostCellWidth())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::"
                << "computeSideDataOfCharacteristicVariablesFromConservativeVariables()\n"
                << "The characteristic variables don't have same ghost cell width."
                << std::endl);
        }
    }
    
    for (int ei = 1; ei < d_dim.getValue() + 5; ei++)
    {
        if (num_ghosts_projection_var != projection_variables[ei]->getGhostCellWidth())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::"
                << "computeSideDataOfCharacteristicVariablesFromConservativeVariables()\n"
                << "The projection variables don't have same ghost cell width."
                << std::endl);
        }
    }
    
    if (num_ghosts_projection_var != num_ghosts_characteristic_var)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "computeSideDataOfCharacteristicVariablesFromConservativeVariables()\n"
            << "The ghost cell width of the projection variables does not match that of"
            << " characteristic variables."
            << std::endl);
    }
    
    for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
    {
        if (num_ghosts_conservative_var[vi] - num_ghosts_characteristic_var
                + (hier::IntVector::getOne(d_dim))*idx_offset < hier::IntVector::getZero(d_dim) ||
            num_ghosts_characteristic_var - num_ghosts_conservative_var[vi]
                + (hier::IntVector::getOne(d_dim))*(idx_offset + 1) > hier::IntVector::getZero(d_dim))
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::"
                << "computeSideDataOfCharacteristicVariablesFromConservativeVariables()\n"
                << "The offset index is too large or the number of ghost of characteristic variable"
                << " is too large."
                << std::endl);
        }
    }
    
    /*
     * Declare containers to store pointers to different data.
     */
    
    std::vector<double*> Q;
    Q.reserve(d_num_eqn);
    
    for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
    {
        int depth = conservative_variables[vi]->getDepth();
        
        for (int di = 0; di < depth; di++)
        {
            Q.push_back(conservative_variables[vi]->getPointer(di));
        }
    }
    
    std::vector<double*> W;
    W.resize(d_num_eqn);
    
    double* e_average     = nullptr;
    double* c_average     = nullptr;
    double* Psi_average   = nullptr;
    double* Gamma_average = nullptr;
    
    if (d_dim == tbox::Dimension(1))
    {
        const int interior_dim_0 = d_interior_dims[0];
        
        const int num_ghosts_0_characteristic_var = num_ghosts_characteristic_var[0];
        const int num_ghosts_0_rho = num_ghosts_conservative_var[0][0];
        const int num_ghosts_0_mom = num_ghosts_conservative_var[1][0];
        const int num_ghosts_0_E = num_ghosts_conservative_var[2][0];
        
        const int idx_offset_rho = idx_offset;
        const int idx_offset_mom = idx_offset;
        const int idx_offset_E = idx_offset;
        
        // Declare pointer to side data of velocity.
        
        double* u_average = nullptr;
        
        /*
         * Compute the characteristic variables in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(0);
        }
        
        u_average     = projection_variables[0]->getPointer(0);
        e_average     = projection_variables[1]->getPointer(0);
        c_average     = projection_variables[3]->getPointer(0);
        Psi_average   = projection_variables[4]->getPointer(0);
        Gamma_average = projection_variables[5]->getPointer(0);
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = -num_ghosts_0_characteristic_var;
             i < interior_dim_0 + 1 + num_ghosts_0_characteristic_var;
             i++)
        {
            // Compute the linear indices.
            const int idx_face = i + num_ghosts_0_characteristic_var;
            const int idx_rho = i + idx_offset_rho + num_ghosts_0_rho;
            const int idx_mom = i + idx_offset_mom + num_ghosts_0_mom;
            const int idx_E = i + idx_offset_E + num_ghosts_0_E;
            
            W[0][idx_face] = ((Gamma_average[idx_face]*(u_average[idx_face]*u_average[idx_face] -
                e_average[idx_face]) + Psi_average[idx_face] + u_average[idx_face]*c_average[idx_face])*
                Q[0][idx_rho] -
                (Gamma_average[idx_face]*u_average[idx_face] + c_average[idx_face])*Q[1][idx_mom] +
                Gamma_average[idx_face]*Q[2][idx_E])/
                    (double(2)*c_average[idx_face]*c_average[idx_face]);
            
            W[1][idx_face] = (-(Gamma_average[idx_face]*(u_average[idx_face]*u_average[idx_face] -
                e_average[idx_face]) + Psi_average[idx_face] - c_average[idx_face]*c_average[idx_face])*
                Q[0][idx_rho] +
                Gamma_average[idx_face]*u_average[idx_face]*Q[1][idx_mom] - Gamma_average[idx_face]*Q[2][idx_E])/
                    (c_average[idx_face]*c_average[idx_face]);
            
            W[2][idx_face] = ((Gamma_average[idx_face]*(u_average[idx_face]*u_average[idx_face] -
                e_average[idx_face]) + Psi_average[idx_face] - u_average[idx_face]*c_average[idx_face])*
                Q[0][idx_rho] -
                (Gamma_average[idx_face]*u_average[idx_face] - c_average[idx_face])*Q[1][idx_mom] +
                Gamma_average[idx_face]*Q[2][idx_E])/
                    (double(2)*c_average[idx_face]*c_average[idx_face]);
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int interior_dim_0 = d_interior_dims[0];
        const int interior_dim_1 = d_interior_dims[1];
        
        const int num_ghosts_0_characteristic_var = num_ghosts_characteristic_var[0];
        const int num_ghosts_1_characteristic_var = num_ghosts_characteristic_var[1];
        const int ghostcell_dim_0_characteristic_var = ghostcell_dims_characteristic_var[0];
        
        const int num_ghosts_0_rho = num_ghosts_conservative_var[0][0];
        const int num_ghosts_1_rho = num_ghosts_conservative_var[0][1];
        const int ghostcell_dim_0_rho = ghostcell_dims_conservative_var[0][0];
        
        const int num_ghosts_0_mom = num_ghosts_conservative_var[1][0];
        const int num_ghosts_1_mom = num_ghosts_conservative_var[1][1];
        const int ghostcell_dim_0_mom = ghostcell_dims_conservative_var[1][0];
        
        const int num_ghosts_0_E = num_ghosts_conservative_var[2][0];
        const int num_ghosts_1_E = num_ghosts_conservative_var[2][1];
        const int ghostcell_dim_0_E = ghostcell_dims_conservative_var[2][0];
        
        const int idx_offset_rho = idx_offset;
        const int idx_offset_mom = idx_offset;
        const int idx_offset_E = idx_offset;
        
        // Declare pointers to side data of velocity.
        
        double* u_average = nullptr;
        double* v_average = nullptr;
        
        /*
         * Compute the characteristic variables in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(0);
        }
        
        u_average     = projection_variables[0]->getPointer(0);
        v_average     = projection_variables[1]->getPointer(0);
        e_average     = projection_variables[2]->getPointer(0);
        c_average     = projection_variables[4]->getPointer(0);
        Psi_average   = projection_variables[5]->getPointer(0);
        Gamma_average = projection_variables[6]->getPointer(0);
        
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
                
                const int idx_rho = (i + idx_offset_rho + num_ghosts_0_rho) +
                    (j + num_ghosts_1_rho)*ghostcell_dim_0_rho;
                
                const int idx_mom = (i + idx_offset_mom + num_ghosts_0_mom) +
                    (j + num_ghosts_1_mom)*ghostcell_dim_0_mom;
                
                const int idx_E = (i + idx_offset_E + num_ghosts_0_E) +
                    (j + num_ghosts_1_E)*ghostcell_dim_0_E;
                
                W[0][idx_face] = ((Gamma_average[idx_face]*(u_average[idx_face]*u_average[idx_face] +
                    v_average[idx_face]*v_average[idx_face] - e_average[idx_face]) + Psi_average[idx_face] +
                    u_average[idx_face]*c_average[idx_face])*Q[0][idx_rho] -
                    (Gamma_average[idx_face]*u_average[idx_face] + c_average[idx_face])*Q[1][idx_mom] -
                    Gamma_average[idx_face]*v_average[idx_face]*Q[2][idx_mom] +
                    Gamma_average[idx_face]*Q[3][idx_E])/
                        (double(2)*c_average[idx_face]*c_average[idx_face]);
                
                W[1][idx_face] = (-(Gamma_average[idx_face]*(u_average[idx_face]*u_average[idx_face] +
                    v_average[idx_face]*v_average[idx_face] - e_average[idx_face]) + Psi_average[idx_face] -
                    c_average[idx_face]*c_average[idx_face])*Q[0][idx_rho] +
                    Gamma_average[idx_face]*u_average[idx_face]*Q[1][idx_mom] +
                    Gamma_average[idx_face]*v_average[idx_face]*Q[2][idx_mom] -
                    Gamma_average[idx_face]*Q[3][idx_E])/
                        (c_average[idx_face]*c_average[idx_face]);
                
                W[2][idx_face] = -Q[0][idx_rho] + Q[2][idx_mom]/v_average[idx_face];
                
                W[3][idx_face] = ((Gamma_average[idx_face]*(u_average[idx_face]*u_average[idx_face] +
                    v_average[idx_face]*v_average[idx_face] - e_average[idx_face]) + Psi_average[idx_face] -
                    u_average[idx_face]*c_average[idx_face])*Q[0][idx_rho] -
                    (Gamma_average[idx_face]*u_average[idx_face] - c_average[idx_face])*Q[1][idx_mom] -
                    Gamma_average[idx_face]*v_average[idx_face]*Q[2][idx_mom] +
                    Gamma_average[idx_face]*Q[3][idx_E])/
                        (double(2)*c_average[idx_face]*c_average[idx_face]);
            }
        }
        
        /*
         * Compute the characteristic variables in the y-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(1);
        }
        
        u_average     = projection_variables[0]->getPointer(1);
        v_average     = projection_variables[1]->getPointer(1);
        e_average     = projection_variables[2]->getPointer(1);
        c_average     = projection_variables[4]->getPointer(1);
        Psi_average   = projection_variables[5]->getPointer(1);
        Gamma_average = projection_variables[6]->getPointer(1);
        
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
                
                const int idx_rho = (i + num_ghosts_0_rho) +
                    (j + idx_offset_rho + num_ghosts_1_rho)*ghostcell_dim_0_rho;
                
                const int idx_mom = (i + num_ghosts_0_mom) +
                    (j + idx_offset_mom + num_ghosts_1_mom)*ghostcell_dim_0_mom;
                
                const int idx_E = (i + num_ghosts_0_E) +
                    (j + idx_offset_E + num_ghosts_1_E)*ghostcell_dim_0_E;
                
                W[0][idx_face] = ((Gamma_average[idx_face]*(u_average[idx_face]*u_average[idx_face] +
                    v_average[idx_face]*v_average[idx_face] - e_average[idx_face]) + Psi_average[idx_face] +
                    v_average[idx_face]*c_average[idx_face])*Q[0][idx_rho] -
                    Gamma_average[idx_face]*u_average[idx_face]*Q[1][idx_mom] -
                    (Gamma_average[idx_face]*v_average[idx_face] + c_average[idx_face])*Q[2][idx_mom] +
                    Gamma_average[idx_face]*Q[3][idx_E])/
                        (double(2)*c_average[idx_face]*c_average[idx_face]);
                
                W[1][idx_face] = -Q[0][idx_rho] + Q[1][idx_mom]/u_average[idx_face];
                
                W[2][idx_face] = (-(Gamma_average[idx_face]*(u_average[idx_face]*u_average[idx_face] +
                    v_average[idx_face]*v_average[idx_face] - e_average[idx_face]) + Psi_average[idx_face] -
                    c_average[idx_face]*c_average[idx_face])*Q[0][idx_rho] +
                    Gamma_average[idx_face]*u_average[idx_face]*Q[1][idx_mom] +
                    Gamma_average[idx_face]*v_average[idx_face]*Q[2][idx_mom] -
                    Gamma_average[idx_face]*Q[3][idx_E])/
                        (c_average[idx_face]*c_average[idx_face]);
                
                W[3][idx_face] = ((Gamma_average[idx_face]*(u_average[idx_face]*u_average[idx_face] +
                    v_average[idx_face]*v_average[idx_face] - e_average[idx_face]) + Psi_average[idx_face] -
                    v_average[idx_face]*c_average[idx_face])*Q[0][idx_rho] -
                    Gamma_average[idx_face]*u_average[idx_face]*Q[1][idx_mom] -
                    (Gamma_average[idx_face]*v_average[idx_face] - c_average[idx_face])*Q[2][idx_mom] +
                    Gamma_average[idx_face]*Q[3][idx_E])/
                        (double(2)*c_average[idx_face]*c_average[idx_face]);
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
        
        const int num_ghosts_0_rho = num_ghosts_conservative_var[0][0];
        const int num_ghosts_1_rho = num_ghosts_conservative_var[0][1];
        const int num_ghosts_2_rho = num_ghosts_conservative_var[0][2];
        const int ghostcell_dim_0_rho = ghostcell_dims_conservative_var[0][0];
        const int ghostcell_dim_1_rho = ghostcell_dims_conservative_var[0][1];
        
        const int num_ghosts_0_mom = num_ghosts_conservative_var[1][0];
        const int num_ghosts_1_mom = num_ghosts_conservative_var[1][1];
        const int num_ghosts_2_mom = num_ghosts_conservative_var[1][2];
        const int ghostcell_dim_0_mom = ghostcell_dims_conservative_var[1][0];
        const int ghostcell_dim_1_mom = ghostcell_dims_conservative_var[1][1];
        
        const int num_ghosts_0_E = num_ghosts_conservative_var[2][0];
        const int num_ghosts_1_E = num_ghosts_conservative_var[2][1];
        const int num_ghosts_2_E = num_ghosts_conservative_var[2][2];
        const int ghostcell_dim_0_E = ghostcell_dims_conservative_var[2][0];
        const int ghostcell_dim_1_E = ghostcell_dims_conservative_var[2][1];
        
        const int idx_offset_rho = idx_offset;
        const int idx_offset_mom = idx_offset;
        const int idx_offset_E = idx_offset;
        
        // Declare pointers to side data of velocity.
        
        double* u_average = nullptr;
        double* v_average = nullptr;
        double* w_average = nullptr;
        
        /*
         * Compute the characteristic variables in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(0);
        }
        
        u_average     = projection_variables[0]->getPointer(0);
        v_average     = projection_variables[1]->getPointer(0);
        w_average     = projection_variables[2]->getPointer(0);
        e_average     = projection_variables[3]->getPointer(0);
        c_average     = projection_variables[5]->getPointer(0);
        Psi_average   = projection_variables[6]->getPointer(0);
        Gamma_average = projection_variables[7]->getPointer(0);
        
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
                    
                    const int idx_rho = (i + idx_offset_rho + num_ghosts_0_rho) +
                        (j + num_ghosts_1_rho)*ghostcell_dim_0_rho +
                        (k + num_ghosts_2_rho)*ghostcell_dim_0_rho*
                            ghostcell_dim_1_rho;
                    
                    const int idx_mom = (i + idx_offset_mom + num_ghosts_0_mom) +
                        (j + num_ghosts_1_mom)*ghostcell_dim_0_mom +
                        (k + num_ghosts_2_mom)*ghostcell_dim_0_mom*
                            ghostcell_dim_1_mom;
                    
                    const int idx_E = (i + idx_offset_E + num_ghosts_0_E) +
                        (j + num_ghosts_1_E)*ghostcell_dim_0_E +
                        (k + num_ghosts_2_E)*ghostcell_dim_0_E*
                            ghostcell_dim_1_E;
                    
                    W[0][idx_face] = ((Gamma_average[idx_face]*(u_average[idx_face]*u_average[idx_face] +
                        v_average[idx_face]*v_average[idx_face] + w_average[idx_face]*w_average[idx_face] -
                        e_average[idx_face]) + Psi_average[idx_face] + u_average[idx_face]*c_average[idx_face])*
                        Q[0][idx_rho] -
                        (Gamma_average[idx_face]*u_average[idx_face] + c_average[idx_face])*Q[1][idx_mom] -
                        Gamma_average[idx_face]*v_average[idx_face]*Q[2][idx_mom] -
                        Gamma_average[idx_face]*w_average[idx_face]*Q[3][idx_mom] +
                        Gamma_average[idx_face]*Q[4][idx_E])/
                            (double(2)*c_average[idx_face]*c_average[idx_face]);
                    
                    W[1][idx_face] = (-(Gamma_average[idx_face]*(u_average[idx_face]*u_average[idx_face] +
                        v_average[idx_face]*v_average[idx_face] + w_average[idx_face]*w_average[idx_face] -
                        e_average[idx_face]) + Psi_average[idx_face] - c_average[idx_face]*c_average[idx_face])*
                        Q[0][idx_rho] +
                        Gamma_average[idx_face]*u_average[idx_face]*Q[1][idx_mom] +
                        Gamma_average[idx_face]*v_average[idx_face]*Q[2][idx_mom] +
                        Gamma_average[idx_face]*w_average[idx_face]*Q[3][idx_mom] -
                        Gamma_average[idx_face]*Q[4][idx_E])/
                            (c_average[idx_face]*c_average[idx_face]);
                    
                    W[2][idx_face] = -Q[0][idx_rho] + Q[2][idx_mom]/v_average[idx_face];
                    
                    W[3][idx_face] = -Q[0][idx_rho] + Q[3][idx_mom]/w_average[idx_face];
                    
                    W[4][idx_face] = ((Gamma_average[idx_face]*(u_average[idx_face]*u_average[idx_face] +
                        v_average[idx_face]*v_average[idx_face] + w_average[idx_face]*w_average[idx_face] -
                        e_average[idx_face]) + Psi_average[idx_face] - u_average[idx_face]*c_average[idx_face])*
                        Q[0][idx_rho] -
                        (Gamma_average[idx_face]*u_average[idx_face] - c_average[idx_face])*Q[1][idx_mom] -
                        Gamma_average[idx_face]*v_average[idx_face]*Q[2][idx_mom] -
                        Gamma_average[idx_face]*w_average[idx_face]*Q[3][idx_mom] +
                        Gamma_average[idx_face]*Q[4][idx_E])/
                            (double(2)*c_average[idx_face]*c_average[idx_face]);
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
        
        u_average     = projection_variables[0]->getPointer(1);
        v_average     = projection_variables[1]->getPointer(1);
        w_average     = projection_variables[2]->getPointer(1);
        e_average     = projection_variables[3]->getPointer(1);
        c_average     = projection_variables[5]->getPointer(1);
        Psi_average   = projection_variables[6]->getPointer(1);
        Gamma_average = projection_variables[7]->getPointer(1);
        
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
                    
                    const int idx_rho = (i + num_ghosts_0_rho) +
                        (j + idx_offset_rho + num_ghosts_1_rho)*ghostcell_dim_0_rho +
                        (k + num_ghosts_2_rho)*ghostcell_dim_0_rho*
                            ghostcell_dim_1_rho;
                    
                    const int idx_mom = (i + num_ghosts_0_mom) +
                        (j + idx_offset_mom + num_ghosts_1_mom)*ghostcell_dim_0_mom +
                        (k + num_ghosts_2_mom)*ghostcell_dim_0_mom*
                            ghostcell_dim_1_mom;
                    
                    const int idx_E = (i + num_ghosts_0_E) +
                        (j + idx_offset_E + num_ghosts_1_E)*ghostcell_dim_0_E +
                        (k + num_ghosts_2_E)*ghostcell_dim_0_E*
                            ghostcell_dim_1_E;
                    
                    W[0][idx_face] = ((Gamma_average[idx_face]*(u_average[idx_face]*u_average[idx_face] +
                        v_average[idx_face]*v_average[idx_face] + w_average[idx_face]*w_average[idx_face] -
                        e_average[idx_face]) + Psi_average[idx_face] + w_average[idx_face]*c_average[idx_face])*
                        Q[0][idx_rho] -
                        Gamma_average[idx_face]*u_average[idx_face]*Q[1][idx_mom] -
                        (Gamma_average[idx_face]*v_average[idx_face] + c_average[idx_face])*Q[2][idx_mom] -
                        Gamma_average[idx_face]*w_average[idx_face]*Q[3][idx_mom] +
                        Gamma_average[idx_face]*Q[4][idx_E])/
                            (double(2)*c_average[idx_face]*c_average[idx_face]);
                    
                    W[1][idx_face] = -Q[0][idx_rho] + Q[1][idx_mom]/u_average[idx_face];
                    
                    W[2][idx_face] = (-(Gamma_average[idx_face]*(u_average[idx_face]*u_average[idx_face] +
                        v_average[idx_face]*v_average[idx_face] + w_average[idx_face]*w_average[idx_face] -
                        e_average[idx_face]) + Psi_average[idx_face] - c_average[idx_face]*c_average[idx_face])*
                        Q[0][idx_rho] +
                        Gamma_average[idx_face]*u_average[idx_face]*Q[1][idx_mom] +
                        Gamma_average[idx_face]*v_average[idx_face]*Q[2][idx_mom] +
                        Gamma_average[idx_face]*w_average[idx_face]*Q[3][idx_mom] -
                        Gamma_average[idx_face]*Q[4][idx_E])/
                            (c_average[idx_face]*c_average[idx_face]);
                    
                    W[3][idx_face] = -Q[0][idx_rho] + Q[3][idx_mom]/w_average[idx_face];
                    
                    W[4][idx_face] = ((Gamma_average[idx_face]*(u_average[idx_face]*u_average[idx_face] +
                        v_average[idx_face]*v_average[idx_face] + w_average[idx_face]*w_average[idx_face] -
                        e_average[idx_face]) + Psi_average[idx_face] - w_average[idx_face]*c_average[idx_face])*
                        Q[0][idx_rho] -
                        Gamma_average[idx_face]*u_average[idx_face]*Q[1][idx_mom] -
                        (Gamma_average[idx_face]*v_average[idx_face] - c_average[idx_face])*Q[2][idx_mom] -
                        Gamma_average[idx_face]*w_average[idx_face]*Q[3][idx_mom] +
                        Gamma_average[idx_face]*Q[4][idx_E])/
                            (double(2)*c_average[idx_face]*c_average[idx_face]);
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
        
        u_average     = projection_variables[0]->getPointer(2);
        v_average     = projection_variables[1]->getPointer(2);
        w_average     = projection_variables[2]->getPointer(2);
        e_average     = projection_variables[3]->getPointer(2);
        c_average     = projection_variables[5]->getPointer(2);
        Psi_average   = projection_variables[6]->getPointer(2);
        Gamma_average = projection_variables[7]->getPointer(2);
        
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
                    
                    const int idx_rho = (i + num_ghosts_0_rho) +
                        (j + num_ghosts_1_rho)*ghostcell_dim_0_rho +
                        (k + idx_offset_rho + num_ghosts_2_rho)*ghostcell_dim_0_rho*
                            ghostcell_dim_1_rho;
                    
                    const int idx_mom = (i + num_ghosts_0_mom) +
                        (j + num_ghosts_1_mom)*ghostcell_dim_0_mom +
                        (k + idx_offset_mom + num_ghosts_2_mom)*ghostcell_dim_0_mom*
                            ghostcell_dim_1_mom;
                    
                    const int idx_E = (i + num_ghosts_0_E) +
                        (j + num_ghosts_1_E)*ghostcell_dim_0_E +
                        (k + idx_offset_E + num_ghosts_2_E)*ghostcell_dim_0_E*
                            ghostcell_dim_1_E;
                    
                    W[0][idx_face] = ((Gamma_average[idx_face]*(u_average[idx_face]*u_average[idx_face] +
                        v_average[idx_face]*v_average[idx_face] + w_average[idx_face]*w_average[idx_face] -
                        e_average[idx_face]) + Psi_average[idx_face] + w_average[idx_face]*c_average[idx_face])*
                        Q[0][idx_rho] -
                        Gamma_average[idx_face]*u_average[idx_face]*Q[1][idx_mom] -
                        Gamma_average[idx_face]*v_average[idx_face]*Q[2][idx_mom] -
                        (Gamma_average[idx_face]*w_average[idx_face] + c_average[idx_face])*Q[3][idx_mom] +
                        Gamma_average[idx_face]*Q[4][idx_E])/
                            (double(2)*c_average[idx_face]*c_average[idx_face]);
                    
                    W[1][idx_face] = -Q[0][idx_rho] + Q[1][idx_mom]/u_average[idx_face];
                    
                    W[2][idx_face] = -Q[0][idx_rho] + Q[2][idx_mom]/v_average[idx_face];
                    
                    W[3][idx_face] = (-(Gamma_average[idx_face]*(u_average[idx_face]*u_average[idx_face] +
                        v_average[idx_face]*v_average[idx_face] + w_average[idx_face]*w_average[idx_face] -
                        e_average[idx_face]) + Psi_average[idx_face] - c_average[idx_face]*c_average[idx_face])*
                        Q[0][idx_rho] +
                        Gamma_average[idx_face]*u_average[idx_face]*Q[1][idx_mom] +
                        Gamma_average[idx_face]*v_average[idx_face]*Q[2][idx_mom] +
                        Gamma_average[idx_face]*w_average[idx_face]*Q[3][idx_mom] -
                        Gamma_average[idx_face]*Q[4][idx_E])/
                            (c_average[idx_face]*c_average[idx_face]);
                    
                    W[4][idx_face] = ((Gamma_average[idx_face]*(u_average[idx_face]*u_average[idx_face] +
                        v_average[idx_face]*v_average[idx_face] + w_average[idx_face]*w_average[idx_face] -
                        e_average[idx_face]) + Psi_average[idx_face] - w_average[idx_face]*c_average[idx_face])*
                        Q[0][idx_rho] -
                        Gamma_average[idx_face]*u_average[idx_face]*Q[1][idx_mom] -
                        Gamma_average[idx_face]*v_average[idx_face]*Q[2][idx_mom] -
                        (Gamma_average[idx_face]*w_average[idx_face] - c_average[idx_face])*Q[3][idx_mom] +
                        Gamma_average[idx_face]*Q[4][idx_E])/
                            (double(2)*c_average[idx_face]*c_average[idx_face]);
                }
            }
        }
    }
}


/*
 * Compute the side data of characteristic variables from primitive variables.
 */
void
FlowModelSingleSpecies::computeSideDataOfCharacteristicVariablesFromPrimitiveVariables(
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
            << ": FlowModelSingleSpecies::"
            << "computeSideDataOfCharacteristicVariablesFromPrimitiveVariables()\n"
            << "The number of characteristic variables are incorrect."
            << std::endl);
    }
    if (static_cast<int>(primitive_variables.size()) != 3)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "computeSideDataOfCharacteristicVariablesFromPrimitiveVariables()\n"
            << "The number of primitive variables are incorrect."
            << std::endl);
    }
    if (primitive_variables[0]->getDepth() != 1 ||
        primitive_variables[1]->getDepth() != d_dim.getValue() ||
        primitive_variables[2]->getDepth() != 1)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "computeSideDataOfCharacteristicVariablesFromPrimitiveVariables()\n"
            << "The depths of one or more primitive variables are incorrect."
            << std::endl);
    }
    if (static_cast<int>(projection_variables.size()) != 2)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "computeSideDataOfCharacteristicVariablesFromPrimitiveVariables()\n"
            << "There should be two projection variables."
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
                << ": FlowModelSingleSpecies::"
                << "computeSideDataOfCharacteristicVariablesFromPrimitiveVariables()\n"
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
                << ": FlowModelSingleSpecies::"
                << "computeSideDataOfCharacteristicVariablesFromPrimitiveVariables()\n"
                << "The interior dimension of the primitive variables does not match that of patch."
                << std::endl);
        }
    }
    for (int vi = 0; vi < 2; vi++)
    {
        const hier::IntVector interior_dims_projection_var = projection_variables[vi]->getBox().numberCells();
        if (interior_dims_projection_var != d_interior_dims)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::"
                << "computeSideDataOfCharacteristicVariablesFromPrimitiveVariables()\n"
                << "The interior dimension of the projection variables does not match that of patch."
                << std::endl);
        }
    }
    
    for (int ei = 1; ei < d_num_eqn; ei++)
    {
        if (num_ghosts_characteristic_var != characteristic_variables[ei]->getGhostCellWidth())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::"
                << "computeSideDataOfCharacteristicVariablesFromPrimitiveVariables()\n"
                << "The characteristic variables don't have same ghost cell width."
                << std::endl);
        }
    }
    if (num_ghosts_projection_var != projection_variables[1]->getGhostCellWidth())
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "computeSideDataOfCharacteristicVariablesFromPrimitiveVariables()\n"
            << "The projection variables don't have same ghost cell width."
            << std::endl);
    }
    
    if (num_ghosts_projection_var != num_ghosts_characteristic_var)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "computeSideDataOfCharacteristicVariablesFromPrimitiveVariables()\n"
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
                << ": FlowModelSingleSpecies::"
                << "computeSideDataOfCharacteristicVariablesFromPrimitiveVariables()\n"
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
    
    double* rho_average = nullptr;
    double* c_average = nullptr;
    
    if (d_dim == tbox::Dimension(1))
    {
        const int interior_dim_0 = d_interior_dims[0];
        
        const int num_ghosts_0_characteristic_var = num_ghosts_characteristic_var[0];
        const int num_ghosts_0_rho = num_ghosts_primitive_var[0][0];
        const int num_ghosts_0_vel = num_ghosts_primitive_var[1][0];
        const int num_ghosts_0_p = num_ghosts_primitive_var[2][0];
        
        const int idx_offset_rho = idx_offset;
        const int idx_offset_vel = idx_offset;
        const int idx_offset_p = idx_offset;
        
        /*
         * Compute the characteristic variables in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(0);
        }
        
        rho_average = projection_variables[0]->getPointer(0);
        c_average = projection_variables[1]->getPointer(0);
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = -num_ghosts_0_characteristic_var;
             i < interior_dim_0 + 1 + num_ghosts_0_characteristic_var;
             i++)
        {
            // Compute the linear indices.
            const int idx_face = i + num_ghosts_0_characteristic_var;
            const int idx_rho = i + idx_offset_rho + num_ghosts_0_rho;
            const int idx_vel = i + idx_offset_vel + num_ghosts_0_vel;
            const int idx_p = i + idx_offset_p + num_ghosts_0_p;
            
            W[0][idx_face] = -double(1)/double(2)*rho_average[idx_face]*c_average[idx_face]*V[1][idx_vel] +
                double(1)/double(2)*V[2][idx_p];
            W[1][idx_face] = V[0][idx_rho] - double(1)/(c_average[idx_face]*c_average[idx_face])*V[2][idx_p];
            W[2][idx_face] = double(1)/double(2)*rho_average[idx_face]*c_average[idx_face]*V[1][idx_vel] +
                double(1)/double(2)*V[2][idx_p];
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int interior_dim_0 = d_interior_dims[0];
        const int interior_dim_1 = d_interior_dims[1];
        
        const int num_ghosts_0_characteristic_var = num_ghosts_characteristic_var[0];
        const int num_ghosts_1_characteristic_var = num_ghosts_characteristic_var[1];
        const int ghostcell_dim_0_characteristic_var = ghostcell_dims_characteristic_var[0];
        
        const int num_ghosts_0_rho = num_ghosts_primitive_var[0][0];
        const int num_ghosts_1_rho = num_ghosts_primitive_var[0][1];
        const int ghostcell_dim_0_rho = ghostcell_dims_primitive_var[0][0];
        
        const int num_ghosts_0_vel = num_ghosts_primitive_var[1][0];
        const int num_ghosts_1_vel = num_ghosts_primitive_var[1][1];
        const int ghostcell_dim_0_vel = ghostcell_dims_primitive_var[1][0];
        
        const int num_ghosts_0_p = num_ghosts_primitive_var[2][0];
        const int num_ghosts_1_p = num_ghosts_primitive_var[2][1];
        const int ghostcell_dim_0_p = ghostcell_dims_primitive_var[2][0];
        
        const int idx_offset_rho = idx_offset;
        const int idx_offset_vel = idx_offset;
        const int idx_offset_p = idx_offset;
        
        /*
         * Compute the characteristic variables in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(0);
        }
        
        rho_average = projection_variables[0]->getPointer(0);
        c_average = projection_variables[1]->getPointer(0);
        
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
                
                const int idx_rho = (i + idx_offset_rho + num_ghosts_0_rho) +
                    (j + num_ghosts_1_rho)*ghostcell_dim_0_rho;
                
                const int idx_vel = (i + idx_offset_vel + num_ghosts_0_vel) +
                    (j + num_ghosts_1_vel)*ghostcell_dim_0_vel;
                
                const int idx_p = (i + idx_offset_p + num_ghosts_0_p) +
                    (j + num_ghosts_1_p)*ghostcell_dim_0_p;
                
                W[0][idx_face] = -double(1)/double(2)*rho_average[idx_face]*c_average[idx_face]*V[1][idx_vel] +
                    double(1)/double(2)*V[3][idx_p];
                W[1][idx_face] = V[0][idx_rho] - double(1)/(c_average[idx_face]*c_average[idx_face])*V[3][idx_p];
                W[2][idx_face] = V[2][idx_vel];
                W[3][idx_face] = double(1)/double(2)*rho_average[idx_face]*c_average[idx_face]*V[1][idx_vel] +
                    double(1)/double(2)*V[3][idx_p];
            }
        }
        
        /*
         * Compute the characteristic variables in the y-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(1);
        }
        
        rho_average = projection_variables[0]->getPointer(1);
        c_average = projection_variables[1]->getPointer(1);
        
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
                
                const int idx_rho = (i + num_ghosts_0_rho) +
                    (j + idx_offset_rho + num_ghosts_1_rho)*ghostcell_dim_0_rho;
                
                const int idx_vel = (i + num_ghosts_0_vel) +
                    (j + idx_offset_vel + num_ghosts_1_vel)*ghostcell_dim_0_vel;
                
                const int idx_p = (i + num_ghosts_0_p) +
                    (j + idx_offset_p + num_ghosts_1_p)*ghostcell_dim_0_p;
                
                W[0][idx_face] = -double(1)/double(2)*rho_average[idx_face]*c_average[idx_face]*V[2][idx_vel] +
                    double(1)/double(2)*V[3][idx_p];
                W[1][idx_face] = V[0][idx_rho] - double(1)/(c_average[idx_face]*c_average[idx_face])*V[3][idx_p];
                W[2][idx_face] = V[1][idx_vel];
                W[3][idx_face] = double(1)/double(2)*rho_average[idx_face]*c_average[idx_face]*V[2][idx_vel] +
                    double(1)/double(2)*V[3][idx_p];
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
        
        const int num_ghosts_0_rho = num_ghosts_primitive_var[0][0];
        const int num_ghosts_1_rho = num_ghosts_primitive_var[0][1];
        const int num_ghosts_2_rho = num_ghosts_primitive_var[0][2];
        const int ghostcell_dim_0_rho = ghostcell_dims_primitive_var[0][0];
        const int ghostcell_dim_1_rho = ghostcell_dims_primitive_var[0][1];
        
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
        
        const int idx_offset_rho = idx_offset;
        const int idx_offset_vel = idx_offset;
        const int idx_offset_p = idx_offset;
        
        /*
         * Compute the characteristic variables in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(0);
        }
        
        rho_average = projection_variables[0]->getPointer(0);
        c_average = projection_variables[1]->getPointer(0);
        
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
                    
                    const int idx_rho = (i + idx_offset_rho + num_ghosts_0_rho) +
                        (j + num_ghosts_1_rho)*ghostcell_dim_0_rho +
                        (k + num_ghosts_2_rho)*ghostcell_dim_0_rho*
                            ghostcell_dim_1_rho;
                    
                    const int idx_vel = (i + idx_offset_vel + num_ghosts_0_vel) +
                        (j + num_ghosts_1_vel)*ghostcell_dim_0_vel +
                        (k + num_ghosts_2_vel)*ghostcell_dim_0_vel*
                            ghostcell_dim_1_vel;
                    
                    const int idx_p = (i + idx_offset_p + num_ghosts_0_p) +
                        (j + num_ghosts_1_p)*ghostcell_dim_0_p +
                        (k + num_ghosts_2_p)*ghostcell_dim_0_p*
                            ghostcell_dim_1_p;
                    
                    W[0][idx_face] = -double(1)/double(2)*rho_average[idx_face]*c_average[idx_face]*V[1][idx_vel] +
                        double(1)/double(2)*V[4][idx_p];
                    W[1][idx_face] = V[0][idx_rho] - double(1)/(c_average[idx_face]*c_average[idx_face])*V[4][idx_p];
                    W[2][idx_face] = V[2][idx_vel];
                    W[3][idx_face] = V[3][idx_vel];
                    W[4][idx_face] = double(1)/double(2)*rho_average[idx_face]*c_average[idx_face]*V[1][idx_vel] +
                        double(1)/double(2)*V[4][idx_p];
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
        
        rho_average = projection_variables[0]->getPointer(1);
        c_average = projection_variables[1]->getPointer(1);
        
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
                    
                    const int idx_rho = (i + num_ghosts_0_rho) +
                        (j + idx_offset_rho + num_ghosts_1_rho)*ghostcell_dim_0_rho +
                        (k + num_ghosts_2_rho)*ghostcell_dim_0_rho*
                            ghostcell_dim_1_rho;
                    
                    const int idx_vel = (i + num_ghosts_0_vel) +
                        (j + idx_offset_vel + num_ghosts_1_vel)*ghostcell_dim_0_vel +
                        (k + num_ghosts_2_vel)*ghostcell_dim_0_vel*
                            ghostcell_dim_1_vel;
                    
                    const int idx_p = (i + num_ghosts_0_p) +
                        (j + idx_offset_p + num_ghosts_1_p)*ghostcell_dim_0_p +
                        (k + num_ghosts_2_p)*ghostcell_dim_0_p*
                            ghostcell_dim_1_p;
                    
                    W[0][idx_face] = -double(1)/double(2)*rho_average[idx_face]*c_average[idx_face]*V[2][idx_vel] +
                        double(1)/double(2)*V[4][idx_p];
                    W[1][idx_face] = V[0][idx_rho] - double(1)/(c_average[idx_face]*c_average[idx_face])*V[4][idx_p];
                    W[2][idx_face] = V[1][idx_vel];
                    W[3][idx_face] = V[3][idx_vel];
                    W[4][idx_face] = double(1)/double(2)*rho_average[idx_face]*c_average[idx_face]*V[2][idx_vel] +
                        double(1)/double(2)*V[4][idx_p];
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
        
        rho_average = projection_variables[0]->getPointer(2);
        c_average = projection_variables[1]->getPointer(2);
        
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
                    
                    const int idx_rho = (i + num_ghosts_0_rho) +
                        (j + num_ghosts_1_rho)*ghostcell_dim_0_rho +
                        (k + idx_offset_rho + num_ghosts_2_rho)*ghostcell_dim_0_rho*
                            ghostcell_dim_1_rho;
                    
                    const int idx_vel = (i + num_ghosts_0_vel) +
                        (j + num_ghosts_1_vel)*ghostcell_dim_0_vel +
                        (k + idx_offset_vel + num_ghosts_2_vel)*ghostcell_dim_0_vel*
                            ghostcell_dim_1_vel;
                    
                    const int idx_p = (i + num_ghosts_0_p) +
                        (j + num_ghosts_1_p)*ghostcell_dim_0_p +
                        (k + idx_offset_p + num_ghosts_2_p)*ghostcell_dim_0_p*
                            ghostcell_dim_1_p;
                    
                    W[0][idx_face] = -double(1)/double(2)*rho_average[idx_face]*c_average[idx_face]*V[3][idx_vel] +
                        double(1)/double(2)*V[4][idx_p];
                    W[1][idx_face] = V[0][idx_rho] - double(1)/(c_average[idx_face]*c_average[idx_face])*V[4][idx_p];
                    W[2][idx_face] = V[1][idx_vel];
                    W[3][idx_face] = V[2][idx_vel];
                    W[4][idx_face] = double(1)/double(2)*rho_average[idx_face]*c_average[idx_face]*V[3][idx_vel] +
                        double(1)/double(2)*V[4][idx_p];
                }
            }
        }
    }
}


/*
 * Compute the side data of conservative variables from characteristic variables.
 */
void
FlowModelSingleSpecies::computeSideDataOfConservativeVariablesFromCharacteristicVariables(
    std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& characteristic_variables,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& projection_variables)
{
    /*
     * Get the numbers of ghost cells of the variables.
     */
    
    const hier::IntVector num_ghosts_conservative_var = conservative_variables[0]->
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
    
    if (static_cast<int>(conservative_variables.size()) != d_num_eqn)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "computeSideDataOfConservativeVariablesFromCharacteristicVariables()\n"
            << "The number of characteristic variables are incorrect."
            << std::endl);
    }
    if (static_cast<int>(characteristic_variables.size()) != d_num_eqn)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "computeSideDataOfConservativeVariablesFromCharacteristicVariables()\n"
            << "The number of conservative variables are incorrect."
            << std::endl);
    }
    if (static_cast<int>(projection_variables.size()) != d_dim.getValue() + 5)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "computeSideDataOfConservativeVariablesFromCharacteristicVariables()\n"
            << "The number of projection variables are incorrect."
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
                << ": FlowModelSingleSpecies::"
                << "computeSideDataOfConservativeVariablesFromCharacteristicVariables()\n"
                << "The interior dimension of the conservative variables does not match that of patch."
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
                << ": FlowModelSingleSpecies::"
                << "computeSideDataOfConservativeVariablesFromCharacteristicVariables()\n"
                << "The interior dimension of the characteristic variables does not match that of patch."
                << std::endl);
        }
    }
    for (int vi = 0; vi < d_dim.getValue() + 5; vi++)
    {
        const hier::IntVector interior_dims_projection_var = projection_variables[vi]->getBox().numberCells();
        if (interior_dims_projection_var != d_interior_dims)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::"
                << "computeSideDataOfConservativeVariablesFromCharacteristicVariables()\n"
                << "The interior dimension of the projection variables does not match that of patch."
                << std::endl);
        }
    }
    
    for (int ei = 1; ei < d_num_eqn; ei++)
    {
        if (num_ghosts_conservative_var != conservative_variables[ei]->getGhostCellWidth())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::"
                << "computeSideDataOfConservativeVariablesFromCharacteristicVariables()\n"
                << "The conservative variables don't have same ghost cell width."
                << std::endl);
        }
    }
    for (int ei = 1; ei < d_num_eqn; ei++)
    {
        if (num_ghosts_characteristic_var != characteristic_variables[ei]->getGhostCellWidth())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::"
                << "computeSideDataOfConservativeVariablesFromCharacteristicVariables()\n"
                << "The characteristic variables don't have same ghost cell width."
                << std::endl);
        }
    }
    for (int ei = 1; ei < d_dim.getValue() + 5; ei++)
    {
        if (num_ghosts_projection_var != projection_variables[ei]->getGhostCellWidth())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::"
                << "computeSideDataOfConservativeVariablesFromCharacteristicVariables()\n"
                << "The projection variables don't have same ghost cell width."
                << std::endl);
        }
    }
    
    if (num_ghosts_projection_var != num_ghosts_conservative_var)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "computeSideDataOfConservativeVariablesFromCharacteristicVariables()\n"
            << "The ghost cell width of the projection variables does not match that of conservative"
            << " variables."
            << std::endl);
    }
    if (num_ghosts_projection_var != num_ghosts_characteristic_var)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "computeSideDataOfConservativeVariablesFromCharacteristicVariables()\n"
            << "The ghost cell width of the projection variables does not match that of characteristic"
            << " variables."
            << std::endl);
    }
    
    /*
     * Declare containers to store pointers to different data.
     */
    
    std::vector<double*> Q;
    std::vector<double*> W;
    Q.resize(d_num_eqn);
    W.resize(d_num_eqn);
    
    double* e_average     = nullptr;
    double* H_average     = nullptr;
    double* c_average     = nullptr;
    double* Psi_average   = nullptr;
    double* Gamma_average = nullptr;
    
    if (d_dim == tbox::Dimension(1))
    {
        const int interior_dim_0 = d_interior_dims[0];
        
        const int num_ghosts_0_characteristic_var = num_ghosts_characteristic_var[0];
        
        // Declare pointer to side data of velocity.
        
        double* u_average = nullptr;
        
        /*
         * Compute the conservative variables in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            Q[ei] = conservative_variables[ei]->getPointer(0);
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(0);
        }
        
        u_average     = projection_variables[0]->getPointer(0);
        e_average     = projection_variables[1]->getPointer(0);
        H_average     = projection_variables[2]->getPointer(0);
        c_average     = projection_variables[3]->getPointer(0);
        Psi_average   = projection_variables[4]->getPointer(0);
        Gamma_average = projection_variables[5]->getPointer(0);
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = -num_ghosts_0_characteristic_var;
             i < interior_dim_0 + 1 + num_ghosts_0_characteristic_var;
             i++)
        {
            // Compute the linear index.
            const int idx_face = i + num_ghosts_0_characteristic_var;
            
            Q[0][idx_face] = W[0][idx_face] + W[1][idx_face] + W[2][idx_face];
            
            Q[1][idx_face] = (u_average[idx_face] - c_average[idx_face])*W[0][idx_face] +
                u_average[idx_face]*W[1][idx_face] +
                (u_average[idx_face] + c_average[idx_face])*W[2][idx_face];
            
            Q[2][idx_face] = (H_average[idx_face] - u_average[idx_face]*c_average[idx_face])*
                W[0][idx_face] +
                (e_average[idx_face] - Psi_average[idx_face]/Gamma_average[idx_face])*W[1][idx_face] +
                (H_average[idx_face] + u_average[idx_face]*c_average[idx_face])*
                W[2][idx_face];
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int interior_dim_0 = d_interior_dims[0];
        const int interior_dim_1 = d_interior_dims[1];
        
        const int num_ghosts_0_characteristic_var = num_ghosts_characteristic_var[0];
        const int num_ghosts_1_characteristic_var = num_ghosts_characteristic_var[1];
        const int ghostcell_dim_0_characteristic_var = ghostcell_dims_characteristic_var[0];
        
        // Declare pointers to side data of velocity.
        
        double* u_average = nullptr;
        double* v_average = nullptr;
        
        /*
         * Compute the conservative variables in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            Q[ei] = conservative_variables[ei]->getPointer(0);
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(0);
        }
        
        u_average     = projection_variables[0]->getPointer(0);
        v_average     = projection_variables[1]->getPointer(0);
        e_average     = projection_variables[2]->getPointer(0);
        H_average     = projection_variables[3]->getPointer(0);
        c_average     = projection_variables[4]->getPointer(0);
        Psi_average   = projection_variables[5]->getPointer(0);
        Gamma_average = projection_variables[6]->getPointer(0);
        
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
                
                Q[0][idx_face] = W[0][idx_face] + W[1][idx_face] + W[3][idx_face];
                
                Q[1][idx_face] = (u_average[idx_face] - c_average[idx_face])*W[0][idx_face] +
                    u_average[idx_face]*W[1][idx_face] +
                    (u_average[idx_face] + c_average[idx_face])*W[3][idx_face];
                
                Q[2][idx_face] = v_average[idx_face]*(W[0][idx_face] + W[1][idx_face] + W[2][idx_face] +
                    W[3][idx_face]);
                
                Q[3][idx_face] = (H_average[idx_face] - u_average[idx_face]*c_average[idx_face])*
                    W[0][idx_face] +
                    (e_average[idx_face] - Psi_average[idx_face]/Gamma_average[idx_face])*W[1][idx_face] +
                    v_average[idx_face]*v_average[idx_face]*W[2][idx_face] +
                    (H_average[idx_face] + u_average[idx_face]*c_average[idx_face])*
                    W[3][idx_face];
            }
        }
        
        /*
         * Compute the conservative variables in the y-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            Q[ei] = conservative_variables[ei]->getPointer(1);
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(1);
        }
        
        u_average     = projection_variables[0]->getPointer(1);
        v_average     = projection_variables[1]->getPointer(1);
        e_average     = projection_variables[2]->getPointer(1);
        H_average     = projection_variables[3]->getPointer(1);
        c_average     = projection_variables[4]->getPointer(1);
        Psi_average   = projection_variables[5]->getPointer(1);
        Gamma_average = projection_variables[6]->getPointer(1);
        
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
                
                Q[0][idx_face] = W[0][idx_face] + W[2][idx_face] + W[3][idx_face];
                
                Q[1][idx_face] = u_average[idx_face]*(W[0][idx_face] + W[1][idx_face] + W[2][idx_face] +
                    W[3][idx_face]);
                
                Q[2][idx_face] = (v_average[idx_face] - c_average[idx_face])*W[0][idx_face] +
                    v_average[idx_face]*W[2][idx_face] +
                    (v_average[idx_face] + c_average[idx_face])*W[3][idx_face];
                
                Q[3][idx_face] = (H_average[idx_face] - v_average[idx_face]*c_average[idx_face])*
                    W[0][idx_face] +
                    u_average[idx_face]*u_average[idx_face]*W[1][idx_face] +
                    (e_average[idx_face] - Psi_average[idx_face]/Gamma_average[idx_face])*W[2][idx_face] +
                    (H_average[idx_face] + v_average[idx_face]*c_average[idx_face])*
                    W[3][idx_face];
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
        
        // Declare pointers to side data of velocity.
        
        double* u_average = nullptr;
        double* v_average = nullptr;
        double* w_average = nullptr;
        
        /*
         * Compute the conservative variables in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            Q[ei] = conservative_variables[ei]->getPointer(0);
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(0);
        }
        
        u_average     = projection_variables[0]->getPointer(0);
        v_average     = projection_variables[1]->getPointer(0);
        w_average     = projection_variables[2]->getPointer(0);
        e_average     = projection_variables[3]->getPointer(0);
        H_average     = projection_variables[4]->getPointer(0);
        c_average     = projection_variables[5]->getPointer(0);
        Psi_average   = projection_variables[6]->getPointer(0);
        Gamma_average = projection_variables[7]->getPointer(0);
        
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
                    
                    Q[0][idx_face] = W[0][idx_face] + W[1][idx_face] + W[4][idx_face];
                    
                    Q[1][idx_face] = (u_average[idx_face] - c_average[idx_face])*W[0][idx_face] +
                        u_average[idx_face]*W[1][idx_face] +
                        (u_average[idx_face] + c_average[idx_face])*W[4][idx_face];
                    
                    Q[2][idx_face] = v_average[idx_face]*(W[0][idx_face] + W[1][idx_face] + W[2][idx_face] +
                        W[4][idx_face]);
                    
                    Q[3][idx_face] = w_average[idx_face]*(W[0][idx_face] + W[1][idx_face] + W[3][idx_face] +
                        W[4][idx_face]);
                    
                    Q[4][idx_face] = (H_average[idx_face] - u_average[idx_face]*c_average[idx_face])*
                        W[0][idx_face] +
                        (e_average[idx_face] - Psi_average[idx_face]/Gamma_average[idx_face])*W[1][idx_face] +
                        v_average[idx_face]*v_average[idx_face]*W[2][idx_face] +
                        w_average[idx_face]*w_average[idx_face]*W[3][idx_face] +
                        (H_average[idx_face] + u_average[idx_face]*c_average[idx_face])*
                        W[4][idx_face];
                }
            }
        }
        
        /*
         * Compute the conservative variables in the y-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            Q[ei] = conservative_variables[ei]->getPointer(1);
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(1);
        }
        
        u_average     = projection_variables[0]->getPointer(1);
        v_average     = projection_variables[1]->getPointer(1);
        w_average     = projection_variables[2]->getPointer(1);
        e_average     = projection_variables[3]->getPointer(1);
        H_average     = projection_variables[4]->getPointer(1);
        c_average     = projection_variables[5]->getPointer(1);
        Psi_average   = projection_variables[6]->getPointer(1);
        Gamma_average = projection_variables[7]->getPointer(1);
        
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
                    
                    Q[0][idx_face] = W[0][idx_face] + W[2][idx_face] + W[4][idx_face];
                    
                    Q[1][idx_face] = u_average[idx_face]*(W[0][idx_face] + W[1][idx_face] + W[2][idx_face] +
                        W[4][idx_face]);
                    
                    Q[2][idx_face] = (v_average[idx_face] - c_average[idx_face])*W[0][idx_face] +
                        v_average[idx_face]*W[2][idx_face] +
                        (v_average[idx_face] + c_average[idx_face])*W[4][idx_face];
                    
                    Q[3][idx_face] = w_average[idx_face]*(W[0][idx_face] + W[2][idx_face] + W[3][idx_face] +
                        W[4][idx_face]);
                    
                    Q[4][idx_face] = (H_average[idx_face] - v_average[idx_face]*c_average[idx_face])*
                        W[0][idx_face] +
                        u_average[idx_face]*u_average[idx_face]*W[1][idx_face] +
                        (e_average[idx_face] - Psi_average[idx_face]/Gamma_average[idx_face])*W[2][idx_face] +
                        w_average[idx_face]*w_average[idx_face]*W[3][idx_face] +
                        (H_average[idx_face] + v_average[idx_face]*c_average[idx_face])*
                        W[4][idx_face];
                }
            }
        }
        
        /*
         * Compute the conservative variables in the z-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            Q[ei] = conservative_variables[ei]->getPointer(2);
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(2);
        }
        
        u_average     = projection_variables[0]->getPointer(2);
        v_average     = projection_variables[1]->getPointer(2);
        w_average     = projection_variables[2]->getPointer(2);
        e_average     = projection_variables[3]->getPointer(2);
        H_average     = projection_variables[4]->getPointer(2);
        c_average     = projection_variables[5]->getPointer(2);
        Psi_average   = projection_variables[6]->getPointer(2);
        Gamma_average = projection_variables[7]->getPointer(2);
        
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
                    
                    Q[0][idx_face] = W[0][idx_face] + W[3][idx_face] + W[4][idx_face];
                    
                    Q[1][idx_face] = u_average[idx_face]*(W[0][idx_face] + W[1][idx_face] + W[3][idx_face] +
                        W[4][idx_face]);
                    
                    Q[2][idx_face] = v_average[idx_face]*(W[0][idx_face] + W[2][idx_face] + W[3][idx_face] +
                        W[4][idx_face]);
                    
                    Q[3][idx_face] = (w_average[idx_face] - c_average[idx_face])*W[0][idx_face] +
                        w_average[idx_face]*W[3][idx_face] +
                        (w_average[idx_face] + c_average[idx_face])*W[4][idx_face];
                    
                    Q[4][idx_face] = (H_average[idx_face] - w_average[idx_face]*c_average[idx_face])*
                        W[0][idx_face] +
                        u_average[idx_face]*u_average[idx_face]*W[1][idx_face] +
                        v_average[idx_face]*v_average[idx_face]*W[2][idx_face] +
                        (e_average[idx_face] - Psi_average[idx_face]/Gamma_average[idx_face])*W[3][idx_face] +
                        (H_average[idx_face] + w_average[idx_face]*c_average[idx_face])*
                        W[4][idx_face];
                }
            }
        }
    }
}


/*
 * Compute the side data of primitive variables from characteristic variables.
 */
void
FlowModelSingleSpecies::computeSideDataOfPrimitiveVariablesFromCharacteristicVariables(
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
            << ": FlowModelSingleSpecies::"
            << "computeSideDataOfPrimitiveVariablesFromCharacteristicVariables()\n"
            << "The number of characteristic variables are incorrect."
            << std::endl);
    }
    if (static_cast<int>(characteristic_variables.size()) != d_num_eqn)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "computeSideDataOfPrimitiveVariablesFromCharacteristicVariables()\n"
            << "The number of primitive variables are incorrect."
            << std::endl);
    }
    if (static_cast<int>(projection_variables.size()) != 2)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "computeSideDataOfPrimitiveVariablesFromCharacteristicVariables()\n"
            << "There should be two projection variables."
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
                << ": FlowModelSingleSpecies::"
                << "computeSideDataOfPrimitiveVariablesFromCharacteristicVariables()\n"
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
                << ": FlowModelSingleSpecies::"
                << "computeSideDataOfPrimitiveVariablesFromCharacteristicVariables()\n"
                << "The interior dimension of the characteristic variables does not match that of patch."
                << std::endl);
        }
    }
    for (int vi = 0; vi < 2; vi++)
    {
        const hier::IntVector interior_dims_projection_var = projection_variables[vi]->getBox().numberCells();
        if (interior_dims_projection_var != d_interior_dims)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::"
                << "computeSideDataOfPrimitiveVariablesFromCharacteristicVariables()\n"
                << "The interior dimension of the projection variables does not match that of patch."
                << std::endl);
        }
    }
    
    for (int ei = 1; ei < d_num_eqn; ei++)
    {
        if (num_ghosts_primitive_var != primitive_variables[ei]->getGhostCellWidth())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::"
                << "computeSideDataOfPrimitiveVariablesFromCharacteristicVariables()\n"
                << "The primitive variables don't have same ghost cell width."
                << std::endl);
        }
    }
    for (int ei = 1; ei < d_num_eqn; ei++)
    {
        if (num_ghosts_characteristic_var != characteristic_variables[ei]->getGhostCellWidth())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::"
                << "computeSideDataOfPrimitiveVariablesFromCharacteristicVariables()\n"
                << "The characteristic variables don't have same ghost cell width."
                << std::endl);
        }
    }
    if (num_ghosts_projection_var != projection_variables[1]->getGhostCellWidth())
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "computeSideDataOfPrimitiveVariablesFromCharacteristicVariables()\n"
            << "The projection variables don't have same ghost cell width."
            << std::endl);
    }
    
    if (num_ghosts_projection_var != num_ghosts_primitive_var)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "computeSideDataOfPrimitiveVariablesFromCharacteristicVariables()\n"
            << "The ghost cell width of the projection variables does not match that of primitive"
            << " variables."
            << std::endl);
    }
    if (num_ghosts_projection_var != num_ghosts_characteristic_var)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "computeSideDataOfPrimitiveVariablesFromCharacteristicVariables()\n"
            << "The ghost cell width of the projection variables does not match that of characteristic"
            << " variables."
            << std::endl);
    }
    
    /*
     * Declare containers to store pointers to different data.
     */
    
    std::vector<double*> V;
    std::vector<double*> W;
    V.resize(d_num_eqn);
    W.resize(d_num_eqn);
    
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
        
        rho_average = projection_variables[0]->getPointer(0);
        c_average = projection_variables[1]->getPointer(0);
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = -num_ghosts_0_characteristic_var;
             i < interior_dim_0 + 1 + num_ghosts_0_characteristic_var;
             i++)
        {
            // Compute the linear index.
            const int idx_face = i + num_ghosts_0_characteristic_var;
            
            V[0][idx_face] = double(1)/(c_average[idx_face]*c_average[idx_face])*W[0][idx_face] +
                W[1][idx_face] + double(1)/(c_average[idx_face]*c_average[idx_face])*W[2][idx_face];
            V[1][idx_face] = -double(1)/(rho_average[idx_face]*c_average[idx_face])*W[0][idx_face] +
                double(1)/(rho_average[idx_face]*c_average[idx_face])*W[2][idx_face];
            V[2][idx_face] = W[0][idx_face] + W[2][idx_face];
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
        
        rho_average = projection_variables[0]->getPointer(0);
        c_average = projection_variables[1]->getPointer(0);
        
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
                
                V[0][idx_face] = double(1)/(c_average[idx_face]*c_average[idx_face])*W[0][idx_face] +
                    W[1][idx_face] + double(1)/(c_average[idx_face]*c_average[idx_face])*W[3][idx_face];
                V[1][idx_face] = -double(1)/(rho_average[idx_face]*c_average[idx_face])*W[0][idx_face] +
                    double(1)/(rho_average[idx_face]*c_average[idx_face])*W[3][idx_face];
                V[2][idx_face] = W[2][idx_face];
                V[3][idx_face] = W[0][idx_face] + W[3][idx_face];
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
        
        rho_average = projection_variables[0]->getPointer(1);
        c_average = projection_variables[1]->getPointer(1);
        
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
                
                V[0][idx_face] = double(1)/(c_average[idx_face]*c_average[idx_face])*W[0][idx_face] + W[1][idx_face] +
                    double(1)/(c_average[idx_face]*c_average[idx_face])*W[3][idx_face];
                V[1][idx_face] = W[2][idx_face];
                V[2][idx_face] = -double(1)/(rho_average[idx_face]*c_average[idx_face])*W[0][idx_face] +
                    double(1)/(rho_average[idx_face]*c_average[idx_face])*W[3][idx_face];
                V[3][idx_face] = W[0][idx_face] + W[3][idx_face];
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
        
        rho_average = projection_variables[0]->getPointer(0);
        c_average = projection_variables[1]->getPointer(0);
        
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
                    
                    V[0][idx_face] = double(1)/(c_average[idx_face]*c_average[idx_face])*W[0][idx_face] +
                        W[1][idx_face] + double(1)/(c_average[idx_face]*c_average[idx_face])*W[4][idx_face];
                    V[1][idx_face] = -double(1)/(rho_average[idx_face]*c_average[idx_face])*W[0][idx_face] +
                        double(1)/(rho_average[idx_face]*c_average[idx_face])*W[4][idx_face];
                    V[2][idx_face] = W[2][idx_face];
                    V[3][idx_face] = W[3][idx_face];
                    V[4][idx_face] = W[0][idx_face] + W[4][idx_face];
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
        
        rho_average = projection_variables[0]->getPointer(1);
        c_average = projection_variables[1]->getPointer(1);
        
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
                    
                    V[0][idx_face] = double(1)/(c_average[idx_face]*c_average[idx_face])*W[0][idx_face] +
                        W[1][idx_face] + double(1)/(c_average[idx_face]*c_average[idx_face])*W[4][idx_face];
                    V[1][idx_face] = W[2][idx_face];
                    V[2][idx_face] = -double(1)/(rho_average[idx_face]*c_average[idx_face])*W[0][idx_face] +
                        double(1)/(rho_average[idx_face]*c_average[idx_face])*W[4][idx_face];
                    V[3][idx_face] = W[3][idx_face];
                    V[4][idx_face] = W[0][idx_face] + W[4][idx_face];
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
        
        rho_average = projection_variables[0]->getPointer(2);
        c_average = projection_variables[1]->getPointer(2);
        
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
                    
                    V[0][idx_face] = double(1)/(c_average[idx_face]*c_average[idx_face])*W[0][idx_face] +
                        W[1][idx_face] + double(1)/(c_average[idx_face]*c_average[idx_face])*W[4][idx_face];
                    V[1][idx_face] = W[2][idx_face];
                    V[2][idx_face] = W[3][idx_face];
                    V[3][idx_face] = -double(1)/(rho_average[idx_face]*c_average[idx_face])*W[0][idx_face] +
                        double(1)/(rho_average[idx_face]*c_average[idx_face])*W[4][idx_face];
                    V[4][idx_face] = W[0][idx_face] + W[4][idx_face];
                }
            }
        }
    }
}


/*
 * Check whether the given cell conservative variables are within the bounds.
 */
void
FlowModelSingleSpecies::checkCellDataOfConservativeVariablesBounded(
    boost::shared_ptr<pdat::CellData<int> >& bounded_flag,
    const std::vector<boost::shared_ptr<pdat::CellData<double> > >& conservative_variables)
{
    // NEED IMPLEMENTATION!
}


/*
 * Check whether the given side conservative variables are within the bounds.
 */
void
FlowModelSingleSpecies::checkSideDataOfConservativeVariablesBounded(
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
            << ": FlowModelSingleSpecies::"
            << "checkSideDataOfConservativeVariablesBounded()\n"
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
                << ": FlowModelSingleSpecies::"
                << "checkSideDataOfConservativeVariablesBounded()\n"
                << "The interior dimension of the conservative variables does not match that of patch."
                << std::endl);
        }
    }
    const hier::IntVector interior_dims_flag = bounded_flag->getBox().numberCells();
    if (interior_dims_flag != d_interior_dims)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "checkSideDataOfConservativeVariablesBounded()\n"
            << "The interior dimension of the flag does not match that of patch."
            << std::endl);
    }
    
    for (int ei = 1; ei < d_num_eqn; ei++)
    {
        if (num_ghosts_conservative_var != conservative_variables[ei]->getGhostCellWidth())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::"
                << "checkSideDataOfConservativeVariablesBounded()\n"
                << "The conservative variables don't have same ghost cell width."
                << std::endl);
        }
    }
    
    if (num_ghosts_flag != num_ghosts_conservative_var)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "checkSideDataOfConservativeVariablesBounded()\n"
            << "The ghost cell width of the flag does not match that of conservative variables."
            << std::endl);
    }
    
    bounded_flag->fillAll(1);
    
    /*
     * Declare containers to store pointers to different data.
     */
    
    std::vector<double*> Q;
    Q.resize(d_num_eqn);
    
    int* are_bounded = nullptr;
    
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
            
            if (Q[0][idx_face] > double(0))
            {
                are_bounded[idx_face] &= 1;
            }
            else
            {
                are_bounded[idx_face] &= 0;
            }
            
            if (Q[d_num_eqn - 1][idx_face] > double(0))
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
                
                if (Q[0][idx_face] > double(0))
                {
                    are_bounded[idx_face] &= 1;
                }
                else
                {
                    are_bounded[idx_face] &= 0;
                }
                
                if (Q[d_num_eqn - 1][idx_face] > double(0))
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
                
                if (Q[0][idx_face] > double(0))
                {
                    are_bounded[idx_face] &= 1;
                }
                else
                {
                    are_bounded[idx_face] &= 0;
                }
                
                if (Q[d_num_eqn - 1][idx_face] > double(0))
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
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1) +
                        (k + num_ghosts_2_conservative_var)*(ghostcell_dim_0_conservative_var + 1)*
                            ghostcell_dim_1_conservative_var;
                    
                    if (Q[0][idx_face] > double(0))
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                    
                    if (Q[d_num_eqn - 1][idx_face] > double(0))
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
                    
                    if (Q[0][idx_face] > double(0))
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                    
                    if (Q[d_num_eqn - 1][idx_face] > double(0))
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
                    
                    if (Q[0][idx_face] > double(0))
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                    
                    if (Q[d_num_eqn - 1][idx_face] > double(0))
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
 * Check whether the given cell primitive variables are within the bounds.
 */
void
FlowModelSingleSpecies::checkCellDataOfPrimitiveVariablesBounded(
    boost::shared_ptr<pdat::CellData<int> >& bounded_flag,
    const std::vector<boost::shared_ptr<pdat::CellData<double> > >& primitive_variables)
{
    // NEED IMPLEMENTATION!
}


/*
 * Check whether the given side primitive variables are within the bounds.
 */
void
FlowModelSingleSpecies::checkSideDataOfPrimitiveVariablesBounded(
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
            << ": FlowModelSingleSpecies::"
            << "checkSideDataOfPrimitiveVariablesBounded()\n"
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
                << ": FlowModelSingleSpecies::"
                << "checkSideDataOfPrimitiveVariablesBounded()\n"
                << "The interior dimension of the primitive variables does not match that of patch."
                << std::endl);
        }
    }
    const hier::IntVector interior_dims_flag = bounded_flag->getBox().numberCells();
    if (interior_dims_flag != d_interior_dims)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "checkSideDataOfPrimitiveVariablesBounded()\n"
            << "The interior dimension of the flag does not match that of patch."
            << std::endl);
    }
    
    for (int ei = 1; ei < d_num_eqn; ei++)
    {
        if (num_ghosts_primitive_var != primitive_variables[ei]->getGhostCellWidth())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::"
                << "checkSideDataOfPrimitiveVariablesBounded()\n"
                << "The primitive variables don't have same ghost cell width."
                << std::endl);
        }
    }
    
    if (num_ghosts_flag != num_ghosts_primitive_var)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "checkSideDataOfPrimitiveVariablesBounded()\n"
            << "The ghost cell width of the flag does not match that of primitive variables."
            << std::endl);
    }
    
    bounded_flag->fillAll(1);
    
    /*
     * Declare containers to store pointers to different data.
     */
    
    std::vector<double*> V;
    V.resize(d_num_eqn);
    
    int* are_bounded = nullptr;
    
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
            
            if (V[0][idx_face] > double(0))
            {
                are_bounded[idx_face] &= 1;
            }
            else
            {
                are_bounded[idx_face] &= 0;
            }
            
            if (V[d_num_eqn - 1][idx_face] > double(0))
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
                
                if (V[0][idx_face] > double(0))
                {
                    are_bounded[idx_face] &= 1;
                }
                else
                {
                    are_bounded[idx_face] &= 0;
                }
                
                if (V[d_num_eqn - 1][idx_face] > double(0))
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
                
                if (V[0][idx_face] > double(0))
                {
                    are_bounded[idx_face] &= 1;
                }
                else
                {
                    are_bounded[idx_face] &= 0;
                }
                
                if (V[d_num_eqn - 1][idx_face] > double(0))
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
                    
                    if (V[0][idx_face] > double(0))
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                    
                    if (V[d_num_eqn - 1][idx_face] > double(0))
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
                    
                    if (V[0][idx_face] > double(0))
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                    
                    if (V[d_num_eqn - 1][idx_face] > double(0))
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
                    
                    if (V[0][idx_face] > double(0))
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                    
                    if (V[d_num_eqn - 1][idx_face] > double(0))
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
 * Convert conservative variables to primitive variables.
 */
void
FlowModelSingleSpecies::convertConservativeVariablesToPrimitiveVariables(
    std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables)
{
    /*
     * Get the numbers of ghost cells of the variables.
     */
    
    const hier::IntVector num_ghosts_primitive_var = primitive_variables[0]->
        getGhostCellWidth();
    
    const hier::IntVector num_ghosts_conservative_var = conservative_variables[0]->
        getGhostCellWidth();
    
    /*
     * Get the ghost cell dimensions of of the variables.
     */
    
    const hier::IntVector ghostcell_dims_primitive_var = primitive_variables[0]->
        getGhostBox().numberCells();
    
    const hier::IntVector ghostcell_dims_conservative_var = conservative_variables[0]->
        getGhostBox().numberCells();
    
    /*
     * Get the size of variables.
     */
    
    int num_eqn_primitive_var = 0;
    int num_eqn_conservative_var = 0;
    
    /*
     * Check the size of variables.
     */
    
    for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
    {
        num_eqn_primitive_var += primitive_variables[vi]->getDepth();
    }
    
    for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
    {
        num_eqn_conservative_var += conservative_variables[vi]->getDepth();
    }
    
    if (num_eqn_primitive_var != d_num_eqn)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "convertConservativeVariablesToPrimitiveVariables()\n"
            << "The number of primitive variables are incorrect."
            << std::endl);
    }
    
    if (num_eqn_conservative_var != d_num_eqn)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "convertConservativeVariablesToPrimitiveVariables()\n"
            << "The number of conservative variables are incorrect."
            << std::endl);
    }
    
    /*
     * Check potential failures.
     */
    
    for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
    {
        const hier::IntVector interior_dims_primitive_var =
            primitive_variables[vi]->getBox().numberCells();
        
        if (interior_dims_primitive_var != d_interior_dims)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::"
                << "convertConservativeVariablesToPrimitiveVariables()\n"
                << "The interior dimension of the primitive variables does not match that of patch."
                << std::endl);
        }
    }
    
    for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
    {
        const hier::IntVector interior_dims_conservative_var =
            conservative_variables[vi]->getBox().numberCells();
        
        if (interior_dims_conservative_var != d_interior_dims)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::"
                << "convertConservativeVariablesToPrimitiveVariables()\n"
                << "The interior dimension of the conservative variables does not match that of patch."
                << std::endl);
        }
    }
    
    for (int vi = 1; vi < static_cast<int>(primitive_variables.size()); vi++)
    {
        if (num_ghosts_primitive_var != primitive_variables[vi]->getGhostCellWidth())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::"
                << "convertConservativeVariablesToPrimitiveVariables()\n"
                << "The primitive variables don't have same ghost cell width."
                << std::endl);
        }
    }
    
    for (int vi = 1; vi < static_cast<int>(conservative_variables.size()); vi++)
    {
        if (num_ghosts_conservative_var != conservative_variables[vi]->getGhostCellWidth())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::"
                << "convertConservativeVariablesToPrimitiveVariables()\n"
                << "The conservative variables don't have same ghost cell width."
                << std::endl);
        }
    }
    
    if (num_ghosts_primitive_var > num_ghosts_conservative_var)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "convertConservativeVariablesToPrimitiveVariables()\n"
            << "The ghost cell width of primitive variables is larger than that of conservative variables."
            << std::endl);
    }
    
    /*
     * Declare the pointers to the primitive variables and conservative variables.
     */
    
    std::vector<double*> V;
    V.resize(num_eqn_primitive_var);
    
    std::vector<double*> Q;
    Q.resize(num_eqn_conservative_var);
    
    int count_eqn = 0;
    
    /*
     * Convert conservative variables to primitive variables.
     */
    
    // Create the temporary side data.
    boost::shared_ptr<pdat::SideData<double> > data_density(
        new pdat::SideData<double>(d_interior_box, 1, num_ghosts_conservative_var));
    
    boost::shared_ptr<pdat::SideData<double> > data_internal_energy(
        new pdat::SideData<double>(d_interior_box, 1, num_ghosts_conservative_var));
    
    boost::shared_ptr<pdat::SideData<double> > data_pressure(
        new pdat::SideData<double>(d_interior_box, 1, num_ghosts_conservative_var));
    
    double* rho     = nullptr;
    double* epsilon = nullptr;
    double* p       = nullptr;
    
    // Get the thermodynamic properties of the species.
    std::vector<const double*> thermo_properties_ptr;
    thermo_properties_ptr.reserve(static_cast<int> (d_thermo_properties.size()));
    for (int ti = 0; ti < static_cast<int> (d_thermo_properties.size()); ti++)
    {
        thermo_properties_ptr.push_back(&d_thermo_properties[ti]);
    }
    
    if (d_dim == tbox::Dimension(1))
    {
        const int interior_dim_0 = d_interior_dims[0];
        
        const int num_ghosts_0_primitive_var    = num_ghosts_primitive_var[0];
        const int num_ghosts_0_conservative_var = num_ghosts_conservative_var[0];
        
        /*
         * Convert conservative variables to primitive variables in the x-direction.
         */
        
        count_eqn = 0;
        for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
        {
            int depth = primitive_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                V[count_eqn] = primitive_variables[vi]->getPointer(0, di);
                count_eqn++;
            }
        }
        
        count_eqn = 0;
        for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
        {
            int depth = conservative_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                Q[count_eqn] = conservative_variables[vi]->getPointer(0, di);
                count_eqn++;
            }
        }
        
        rho     = data_density->getPointer(0, 0);
        epsilon = data_internal_energy->getPointer(0, 0);
        p       = data_pressure->getPointer(0, 0);
        
        // Get the density.
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = -num_ghosts_0_conservative_var;
                i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                i++)
        {
            // Compute the linear index.
            const int idx_conservative_var = i + num_ghosts_0_conservative_var;
            
            rho[idx_conservative_var] = Q[0][idx_conservative_var];
        }
        
        // Compute the internal energy.
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = -num_ghosts_0_conservative_var;
                i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                i++)
        {
            // Compute the linear index.
            const int idx_conservative_var = i + num_ghosts_0_conservative_var;
            
            epsilon[idx_conservative_var] = (Q[1 + d_dim.getValue()][idx_conservative_var] -
                double(1)/double(2)*(Q[1][idx_conservative_var]*Q[1][idx_conservative_var])/
                rho[idx_conservative_var])/rho[idx_conservative_var];
        }
        
        // Compute the pressure.
        d_equation_of_state_mixing_rules->getEquationOfState()->computePressure(
            data_pressure,
            data_density,
            data_internal_energy,
            thermo_properties_ptr,
            0);
        
        // Set the density.
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = -num_ghosts_0_primitive_var;
                i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                i++)
        {
            // Compute the linear indices.
            const int idx_primitive_var    = i + num_ghosts_0_primitive_var;
            const int idx_conservative_var = i + num_ghosts_0_conservative_var;
            
            V[0][idx_primitive_var] = Q[0][idx_conservative_var];
        }
        
        // Set the velocity.
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = -num_ghosts_0_primitive_var;
                i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                i++)
        {
            // Compute the linear indices.
            const int idx_primitive_var    = i + num_ghosts_0_primitive_var;
            const int idx_conservative_var = i + num_ghosts_0_conservative_var;
            
            V[1][idx_primitive_var] = Q[1][idx_conservative_var]/rho[idx_conservative_var];
        }
        
        // Set the pressure.
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = -num_ghosts_0_primitive_var;
                i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                i++)
        {
            // Compute the linear indices.
            const int idx_primitive_var    = i + num_ghosts_0_primitive_var;
            const int idx_conservative_var = i + num_ghosts_0_conservative_var;
            
            V[1 + d_dim.getValue()][idx_primitive_var] = p[idx_conservative_var];
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int interior_dim_0 = d_interior_dims[0];
        const int interior_dim_1 = d_interior_dims[1];
        
        const int num_ghosts_0_primitive_var = num_ghosts_primitive_var[0];
        const int num_ghosts_1_primitive_var = num_ghosts_primitive_var[1];
        const int ghostcell_dim_0_primitive_var = ghostcell_dims_primitive_var[0];
        
        const int num_ghosts_0_conservative_var = num_ghosts_conservative_var[0];
        const int num_ghosts_1_conservative_var = num_ghosts_conservative_var[1];
        const int ghostcell_dim_0_conservative_var = ghostcell_dims_conservative_var[0];
        
        /*
         * Convert conservative variables to primitive variables in the x-direction.
         */
        
        count_eqn = 0;
        for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
        {
            int depth = primitive_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                V[count_eqn] = primitive_variables[vi]->getPointer(0, di);
                count_eqn++;
            }
        }
        
        count_eqn = 0;
        for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
        {
            int depth = conservative_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                Q[count_eqn] = conservative_variables[vi]->getPointer(0, di);
                count_eqn++;
            }
        }
        
        rho     = data_density->getPointer(0, 0);
        epsilon = data_internal_energy->getPointer(0, 0);
        p       = data_pressure->getPointer(0, 0);
        
        // Get the density.
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
                const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                    (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1);
                
                rho[idx_conservative_var] = Q[0][idx_conservative_var];
            }
        }
        
        // Compute the internal energy.
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
                const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                    (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1);
                
                epsilon[idx_conservative_var] = (Q[1 + d_dim.getValue()][idx_conservative_var] -
                    double(1)/double(2)*(Q[1][idx_conservative_var]*Q[1][idx_conservative_var] +
                    Q[2][idx_conservative_var]*Q[2][idx_conservative_var])/
                    rho[idx_conservative_var])/rho[idx_conservative_var];
            }
        }
        
        // Compute the pressure.
        d_equation_of_state_mixing_rules->getEquationOfState()->computePressure(
            data_pressure,
            data_density,
            data_internal_energy,
            thermo_properties_ptr,
            0);
        
        // Set the density.
        for (int j = 0; j < interior_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = -num_ghosts_0_primitive_var;
                    i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                    i++)
            {
                // Compute the linear indices.
                const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1);
                
                const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                    (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1);
                
                V[0][idx_primitive_var] = Q[0][idx_conservative_var];
            }
        }
        
        // Set the velocity.
        for (int j = 0; j < interior_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = -num_ghosts_0_primitive_var;
                    i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                    i++)
            {
                // Compute the linear indices.
                const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1);
                
                const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                    (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1);
                
                V[1][idx_primitive_var] = Q[1][idx_conservative_var]/rho[idx_conservative_var];
                V[2][idx_primitive_var] = Q[2][idx_conservative_var]/rho[idx_conservative_var];
            }
        }
        
        // Set the pressure.
        for (int j = 0; j < interior_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = -num_ghosts_0_primitive_var;
                    i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                    i++)
            {
                // Compute the linear indices.
                const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1);
                
                const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                    (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1);
                
                V[1 + d_dim.getValue()][idx_primitive_var] = p[idx_conservative_var];
            }
        }
        
        /*
         * Convert conservative variables to primitive variables in the y-direction.
         */
        
        count_eqn = 0;
        for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
        {
            int depth = primitive_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                V[count_eqn] = primitive_variables[vi]->getPointer(1, di);
                count_eqn++;
            }
        }
        
        count_eqn = 0;
        for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
        {
            int depth = conservative_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                Q[count_eqn] = conservative_variables[vi]->getPointer(1, di);
                count_eqn++;
            }
        }
        
        rho     = data_density->getPointer(1, 0);
        epsilon = data_internal_energy->getPointer(1, 0);
        p       = data_pressure->getPointer(1, 0);
        
        // Get the density.
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
                const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                    (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var;
                
                rho[idx_conservative_var] = Q[0][idx_conservative_var];
            }
        }
        
        // Compute the internal energy.
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
                const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                    (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var;
                
                epsilon[idx_conservative_var] = (Q[1 + d_dim.getValue()][idx_conservative_var] -
                    double(1)/double(2)*(Q[1][idx_conservative_var]*Q[1][idx_conservative_var] +
                    Q[2][idx_conservative_var]*Q[2][idx_conservative_var])/
                    rho[idx_conservative_var])/rho[idx_conservative_var];
            }
        }
        
        // Compute the pressure.
        d_equation_of_state_mixing_rules->getEquationOfState()->computePressure(
            data_pressure,
            data_density,
            data_internal_energy,
            thermo_properties_ptr,
            1);
        
        // Set the density.
        for (int j = -num_ghosts_1_primitive_var;
            j < interior_dim_1 + 1 + num_ghosts_1_primitive_var;
            j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var;
                
                const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                    (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var;
                
                V[0][idx_primitive_var] = Q[0][idx_conservative_var];
            }
        }
        
        // Set the velocity.
        for (int j = -num_ghosts_1_primitive_var;
            j < interior_dim_1 + 1 + num_ghosts_1_primitive_var;
            j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var;
                
                const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                    (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var;
                
                V[1][idx_primitive_var] = Q[1][idx_conservative_var]/rho[idx_conservative_var];
                V[2][idx_primitive_var] = Q[2][idx_conservative_var]/rho[idx_conservative_var];
            }
        }
        
        // Set the pressure.
        for (int j = -num_ghosts_1_primitive_var;
            j < interior_dim_1 + 1 + num_ghosts_1_primitive_var;
            j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var;
                
                const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                    (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var;
                
                V[1 + d_dim.getValue()][idx_primitive_var] = p[idx_conservative_var];
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
        
        const int num_ghosts_0_conservative_var = num_ghosts_conservative_var[0];
        const int num_ghosts_1_conservative_var = num_ghosts_conservative_var[1];
        const int num_ghosts_2_conservative_var = num_ghosts_conservative_var[2];
        const int ghostcell_dim_0_conservative_var = ghostcell_dims_conservative_var[0];
        const int ghostcell_dim_1_conservative_var = ghostcell_dims_conservative_var[1];
        
        /*
         * Convert conservative variables to primitive variables in the x-direction.
         */
        
        count_eqn = 0;
        for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
        {
            int depth = primitive_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                V[count_eqn] = primitive_variables[vi]->getPointer(0, di);
                count_eqn++;
            }
        }
        
        count_eqn = 0;
        for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
        {
            int depth = conservative_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                Q[count_eqn] = conservative_variables[vi]->getPointer(0, di);
                count_eqn++;
            }
        }
        
        rho     = data_density->getPointer(0, 0);
        epsilon = data_internal_energy->getPointer(0, 0);
        p       = data_pressure->getPointer(0, 0);
        
        // Get the density.
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
                    // Compute the linear index.
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1) +
                        (k + num_ghosts_2_conservative_var)*(ghostcell_dim_0_conservative_var + 1)*
                            ghostcell_dim_1_conservative_var;
                    
                    rho[idx_conservative_var] = Q[0][idx_conservative_var];
                }
            }
        }
        
        // Compute the internal energy.
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
                    // Compute the linear index.
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1) +
                        (k + num_ghosts_2_conservative_var)*(ghostcell_dim_0_conservative_var + 1)*
                            ghostcell_dim_1_conservative_var;
                    
                    epsilon[idx_conservative_var] = (Q[1 + d_dim.getValue()][idx_conservative_var] -
                        double(1)/double(2)*(Q[1][idx_conservative_var]*Q[1][idx_conservative_var] +
                        Q[2][idx_conservative_var]*Q[2][idx_conservative_var] +
                        Q[3][idx_conservative_var]*Q[3][idx_conservative_var])/
                        rho[idx_conservative_var])/rho[idx_conservative_var];
                }
            }
        }
        
        // Compute the pressure.
        d_equation_of_state_mixing_rules->getEquationOfState()->computePressure(
            data_pressure,
            data_density,
            data_internal_energy,
            thermo_properties_ptr,
            0);
        
        // Set the density.
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
                    // Compute the linear indices.
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1) +
                        (k + num_ghosts_2_primitive_var)*(ghostcell_dim_0_primitive_var + 1)*
                            ghostcell_dim_1_primitive_var;
                    
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1) +
                        (k + num_ghosts_2_conservative_var)*(ghostcell_dim_0_conservative_var + 1)*
                            ghostcell_dim_1_conservative_var;
                    
                    V[0][idx_primitive_var] = Q[0][idx_conservative_var];
                }
            }
        }
        
        // Set the velocity.
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
                    // Compute the linear indices.
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1) +
                        (k + num_ghosts_2_primitive_var)*(ghostcell_dim_0_primitive_var + 1)*
                            ghostcell_dim_1_primitive_var;
                    
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1) +
                        (k + num_ghosts_2_conservative_var)*(ghostcell_dim_0_conservative_var + 1)*
                            ghostcell_dim_1_conservative_var;
                    
                    V[1][idx_primitive_var] = Q[1][idx_conservative_var]/rho[idx_conservative_var];
                    V[2][idx_primitive_var] = Q[2][idx_conservative_var]/rho[idx_conservative_var];
                    V[3][idx_primitive_var] = Q[3][idx_conservative_var]/rho[idx_conservative_var];
                }
            }
        }
        
        // Set the pressure.
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
                    // Compute the linear indices.
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1) +
                        (k + num_ghosts_2_primitive_var)*(ghostcell_dim_0_primitive_var + 1)*
                            ghostcell_dim_1_primitive_var;
                    
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1) +
                        (k + num_ghosts_2_conservative_var)*(ghostcell_dim_0_conservative_var + 1)*
                            ghostcell_dim_1_conservative_var;
                    
                    V[1 + d_dim.getValue()][idx_primitive_var] = p[idx_conservative_var];
                }
            }
        }
        
        /*
         * Convert conservative variables to primitive variables in the y-direction.
         */
        
        count_eqn = 0;
        for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
        {
            int depth = primitive_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                V[count_eqn] = primitive_variables[vi]->getPointer(1, di);
                count_eqn++;
            }
        }
        
        count_eqn = 0;
        for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
        {
            int depth = conservative_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                Q[count_eqn] = conservative_variables[vi]->getPointer(1, di);
                count_eqn++;
            }
        }
        
        rho     = data_density->getPointer(1, 0);
        epsilon = data_internal_energy->getPointer(1, 0);
        p       = data_pressure->getPointer(1, 0);
        
        // Get the density.
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
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                        (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                            (ghostcell_dim_1_conservative_var + 1);
                    
                    rho[idx_conservative_var] = Q[0][idx_conservative_var];
                }
            }
        }
        
        // Compute the internal energy.
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
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                        (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                            (ghostcell_dim_1_conservative_var + 1);
                    
                    epsilon[idx_conservative_var] = (Q[1 + d_dim.getValue()][idx_conservative_var] -
                        double(1)/double(2)*(Q[1][idx_conservative_var]*Q[1][idx_conservative_var] +
                        Q[2][idx_conservative_var]*Q[2][idx_conservative_var] +
                        Q[3][idx_conservative_var]*Q[3][idx_conservative_var])/
                        rho[idx_conservative_var])/rho[idx_conservative_var];
                }
            }
        }
        
        // Compute the pressure.
        d_equation_of_state_mixing_rules->getEquationOfState()->computePressure(
            data_pressure,
            data_density,
            data_internal_energy,
            thermo_properties_ptr,
            1);
        
        // Set the density.
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
                    // Compute the linear indices.
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                        (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                            (ghostcell_dim_1_primitive_var + 1);
                    
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                        (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                            (ghostcell_dim_1_conservative_var + 1);
                    
                    V[0][idx_primitive_var] = Q[0][idx_conservative_var];
                }
            }
        }
        
        // Set the velocity.
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
                    // Compute the linear indices.
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                        (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                            (ghostcell_dim_1_primitive_var + 1);
                    
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                        (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                            (ghostcell_dim_1_conservative_var + 1);
                    
                    V[1][idx_primitive_var] = Q[1][idx_conservative_var]/rho[idx_conservative_var];
                    V[2][idx_primitive_var] = Q[2][idx_conservative_var]/rho[idx_conservative_var];
                    V[3][idx_primitive_var] = Q[3][idx_conservative_var]/rho[idx_conservative_var];
                }
            }
        }
        
        // Set the pressure.
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
                    // Compute the linear indices.
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                        (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                            (ghostcell_dim_1_primitive_var + 1);
                    
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                        (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                            (ghostcell_dim_1_conservative_var + 1);
                    
                    V[1 + d_dim.getValue()][idx_primitive_var] = p[idx_conservative_var];
                }
            }
        }
        
        /*
         * Convert conservative variables to primitive variables in the z-direction.
         */
        
        count_eqn = 0;
        for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
        {
            int depth = primitive_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                V[count_eqn] = primitive_variables[vi]->getPointer(2, di);
                count_eqn++;
            }
        }
        
        count_eqn = 0;
        for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
        {
            int depth = conservative_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                Q[count_eqn] = conservative_variables[vi]->getPointer(2, di);
                count_eqn++;
            }
        }
        
        rho     = data_density->getPointer(2, 0);
        epsilon = data_internal_energy->getPointer(2, 0);
        p       = data_pressure->getPointer(2, 0);
        
        // Get the density.
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
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                        (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                            ghostcell_dim_1_conservative_var;
                    
                    rho[idx_conservative_var] = Q[0][idx_conservative_var];
                }
            }
        }
        
        // Compute the internal energy.
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
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                        (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                            ghostcell_dim_1_conservative_var;
                    
                    epsilon[idx_conservative_var] = (Q[1 + d_dim.getValue()][idx_conservative_var] -
                        double(1)/double(2)*(Q[1][idx_conservative_var]*Q[1][idx_conservative_var] +
                        Q[2][idx_conservative_var]*Q[2][idx_conservative_var] +
                        Q[3][idx_conservative_var]*Q[3][idx_conservative_var])/
                        rho[idx_conservative_var])/rho[idx_conservative_var];
                }
            }
        }
        
        // Compute the pressure.
        d_equation_of_state_mixing_rules->getEquationOfState()->computePressure(
            data_pressure,
            data_density,
            data_internal_energy,
            thermo_properties_ptr,
            2);
        
        // Set the density.
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
                    // Compute the linear indices.
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                        (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                            ghostcell_dim_1_primitive_var;
                    
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                        (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                            ghostcell_dim_1_conservative_var;
                    
                    V[0][idx_primitive_var] = Q[0][idx_conservative_var];
                }
            }
        }
        
        // Set the velocity.
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
                    // Compute the linear indices.
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                        (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                            ghostcell_dim_1_primitive_var;
                    
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                        (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                            ghostcell_dim_1_conservative_var;
                    
                    V[1][idx_primitive_var] = Q[1][idx_conservative_var]/rho[idx_conservative_var];
                    V[2][idx_primitive_var] = Q[2][idx_conservative_var]/rho[idx_conservative_var];
                    V[3][idx_primitive_var] = Q[3][idx_conservative_var]/rho[idx_conservative_var];
                }
            }
        }
        
        // Set the pressure.
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
                    // Compute the linear indices.
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                        (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                            ghostcell_dim_1_primitive_var;
                    
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                        (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                            ghostcell_dim_1_conservative_var;
                    
                    V[1 + d_dim.getValue()][idx_primitive_var] = p[idx_conservative_var];
                }
            }
        }
    }
}


/*
 * Convert primitive variables to conservative variables.
 */
void
FlowModelSingleSpecies::convertPrimitiveVariablesToConservativeVariables(
    std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables)
{
    /*
     * Get the numbers of ghost cells of the variables.
     */
    
    const hier::IntVector num_ghosts_conservative_var = conservative_variables[0]->
        getGhostCellWidth();
    
    const hier::IntVector num_ghosts_primitive_var = primitive_variables[0]->
        getGhostCellWidth();
    
    /*
     * Get the ghost cell dimensions of of the variables.
     */
    
    const hier::IntVector ghostcell_dims_conservative_var = conservative_variables[0]->
        getGhostBox().numberCells();
    
    const hier::IntVector ghostcell_dims_primitive_var = primitive_variables[0]->
        getGhostBox().numberCells();
    
    /*
     * Get the size of variables.
     */
    
    int num_eqn_conservative_var = 0;
    int num_eqn_primitive_var = 0;
    
    /*
     * Check the size of variables.
     */
    
    for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
    {
        num_eqn_conservative_var += conservative_variables[vi]->getDepth();
    }
    
    for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
    {
        num_eqn_primitive_var += primitive_variables[vi]->getDepth();
    }
    
    if (num_eqn_conservative_var != d_num_eqn)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "convertPrimitiveVariablesToConservativeVariables()\n"
            << "The number of conservative variables are incorrect."
            << std::endl);
    }
    
    if (num_eqn_primitive_var != d_num_eqn)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "convertPrimitiveVariablesToConservativeVariables()\n"
            << "The number of primitive variables are incorrect."
            << std::endl);
    }
    
    /*
     * Check potential failures.
     */
    
    for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
    {
        const hier::IntVector interior_dims_conservative_var =
            conservative_variables[vi]->getBox().numberCells();
        
        if (interior_dims_conservative_var != d_interior_dims)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::"
                << "convertPrimitiveVariablesToConservativeVariables()\n"
                << "The interior dimension of the conservative variables does not match that of patch."
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
                << ": FlowModelSingleSpecies::"
                << "convertPrimitiveVariablesToConservativeVariables()\n"
                << "The interior dimension of the primitive variables does not match that of patch."
                << std::endl);
        }
    }
    
    for (int vi = 1; vi < static_cast<int>(conservative_variables.size()); vi++)
    {
        if (num_ghosts_conservative_var != conservative_variables[vi]->getGhostCellWidth())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::"
                << "convertPrimitiveVariablesToConservativeVariables()\n"
                << "The conservative variables don't have same ghost cell width."
                << std::endl);
        }
    }
    
    for (int vi = 1; vi < static_cast<int>(primitive_variables.size()); vi++)
    {
        if (num_ghosts_primitive_var != primitive_variables[vi]->getGhostCellWidth())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::"
                << "convertPrimitiveVariablesToConservativeVariables()\n"
                << "The primitive variables don't have same ghost cell width."
                << std::endl);
        }
    }
    
    if (num_ghosts_conservative_var > num_ghosts_primitive_var)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "convertPrimitiveVariablesToConservativeVariables()\n"
            << "The ghost cell width of conservative variables is larger than that of primitive variables."
            << std::endl);
    }
    
    /*
     * Declare the pointers to the conservative variables and primitive variables.
     */
    
    std::vector<double*> Q;
    Q.resize(num_eqn_conservative_var);
    
    std::vector<double*> V;
    V.resize(num_eqn_primitive_var);
    
    int count_eqn = 0;
    
    /*
     * Convert primitive variables to conservative variables.
     */
    
    // Create the temporary side data.
    boost::shared_ptr<pdat::SideData<double> > data_density(
        new pdat::SideData<double>(d_interior_box, 1, num_ghosts_primitive_var));
    
    boost::shared_ptr<pdat::SideData<double> > data_pressure(
        new pdat::SideData<double>(d_interior_box, 1, num_ghosts_primitive_var));
    
    boost::shared_ptr<pdat::SideData<double> > data_internal_energy(
        new pdat::SideData<double>(d_interior_box, 1, num_ghosts_primitive_var));
    
    double* rho     = nullptr;
    double* p       = nullptr;
    double* epsilon = nullptr;
    
    // Get the thermodynamic properties of the species.
    std::vector<const double*> thermo_properties_ptr;
    thermo_properties_ptr.reserve(static_cast<int> (d_thermo_properties.size()));
    for (int ti = 0; ti < static_cast<int> (d_thermo_properties.size()); ti++)
    {
        thermo_properties_ptr.push_back(&d_thermo_properties[ti]);
    }
    
    if (d_dim == tbox::Dimension(1))
    {
        const int interior_dim_0 = d_interior_dims[0];
        
        const int num_ghosts_0_conservative_var = num_ghosts_conservative_var[0];
        const int num_ghosts_0_primitive_var    = num_ghosts_primitive_var[0];
        
        /*
         * Convert primitive variables to conservative variables in the x-direction.
         */
        
        count_eqn = 0;
        for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
        {
            int depth = conservative_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                Q[count_eqn] = conservative_variables[vi]->getPointer(0, di);
                count_eqn++;
            }
        }
        
        count_eqn = 0;
        for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
        {
            int depth = primitive_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                V[count_eqn] = primitive_variables[vi]->getPointer(0, di);
                count_eqn++;
            }
        }
        
        rho     = data_density->getPointer(0, 0);
        epsilon = data_internal_energy->getPointer(0, 0);
        p       = data_pressure->getPointer(0, 0);
        
        // Get the density.
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = -num_ghosts_0_primitive_var;
                i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                i++)
        {
            // Compute the linear index.
            const int idx_primitive_var = i + num_ghosts_0_primitive_var;
            
            rho[idx_primitive_var] = V[0][idx_primitive_var];
        }
        
        // Get the pressure.
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = -num_ghosts_0_primitive_var;
                i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                i++)
        {
            // Compute the linear index.
            const int idx_primitive_var = i + num_ghosts_0_primitive_var;
            
            p[idx_primitive_var] = V[1 + d_dim.getValue()][idx_primitive_var];
        }
        
        // Compute the specific internal energy.
        d_equation_of_state_mixing_rules->getEquationOfState()->computeInternalEnergy(
            data_internal_energy,
            data_density,
            data_pressure,
            thermo_properties_ptr,
            0);
        
        // Set the density.
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = -num_ghosts_0_conservative_var;
                i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                i++)
        {
            // Compute the linear indices.
            const int idx_conservative_var = i + num_ghosts_0_conservative_var;
            const int idx_primitive_var    = i + num_ghosts_0_primitive_var;
            
            Q[0][idx_conservative_var] = V[0][idx_primitive_var];
        }
        
        // Set the momentum.
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = -num_ghosts_0_conservative_var;
                i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                i++)
        {
            // Compute the linear indices.
            const int idx_conservative_var = i + num_ghosts_0_conservative_var;
            const int idx_primitive_var    = i + num_ghosts_0_primitive_var;
            
            Q[1][idx_conservative_var] = rho[idx_primitive_var]*V[1][idx_primitive_var];
        }
        
        // Set the total energy.
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = -num_ghosts_0_conservative_var;
                i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                i++)
        {
            // Compute the linear indices.
            const int idx_conservative_var = i + num_ghosts_0_conservative_var;
            const int idx_primitive_var    = i + num_ghosts_0_primitive_var;
            
            Q[2][idx_conservative_var] = rho[idx_primitive_var]*
                (epsilon[idx_primitive_var] + double(1)/double(2)*
                V[1][idx_primitive_var]*V[1][idx_primitive_var]);
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int interior_dim_0 = d_interior_dims[0];
        const int interior_dim_1 = d_interior_dims[1];
        
        const int num_ghosts_0_conservative_var = num_ghosts_conservative_var[0];
        const int num_ghosts_1_conservative_var = num_ghosts_conservative_var[1];
        const int ghostcell_dim_0_conservative_var = ghostcell_dims_conservative_var[0];
        
        const int num_ghosts_0_primitive_var = num_ghosts_primitive_var[0];
        const int num_ghosts_1_primitive_var = num_ghosts_primitive_var[1];
        const int ghostcell_dim_0_primitive_var = ghostcell_dims_primitive_var[0];
        
        /*
         * Convert primitive variables to conservative variables in the x-direction.
         */
        
        count_eqn = 0;
        for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
        {
            int depth = conservative_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                Q[count_eqn] = conservative_variables[vi]->getPointer(0, di);
                count_eqn++;
            }
        }
        
        count_eqn = 0;
        for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
        {
            int depth = primitive_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                V[count_eqn] = primitive_variables[vi]->getPointer(0, di);
                count_eqn++;
            }
        }
        
        rho     = data_density->getPointer(0, 0);
        epsilon = data_internal_energy->getPointer(0, 0);
        p       = data_pressure->getPointer(0, 0);
        
        // Get the density.
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
                const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1);
                
                rho[idx_primitive_var] = V[0][idx_primitive_var];
            }
        }
        
        // Get the pressure.
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
                const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1);
                
                p[idx_primitive_var] = V[1 + d_dim.getValue()][idx_primitive_var];
            }
        }
        
        // Compute the specific internal energy.
        d_equation_of_state_mixing_rules->getEquationOfState()->computeInternalEnergy(
            data_internal_energy,
            data_density,
            data_pressure,
            thermo_properties_ptr,
            0);
        
        // Set the density.
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
                const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                    (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1);
                
                const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1);
                
                Q[0][idx_conservative_var] = V[0][idx_primitive_var];
            }
        }
        
        // Set the momentum.
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
                const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                    (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1);
                
                const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1);
                
                Q[1][idx_conservative_var] = rho[idx_primitive_var]*V[1][idx_primitive_var];
                Q[2][idx_conservative_var] = rho[idx_primitive_var]*V[2][idx_primitive_var];
            }
        }
        
        // Set the total energy.
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
                const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                    (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1);
                
                const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1);
                
                Q[3][idx_conservative_var] = rho[idx_primitive_var]*
                    (epsilon[idx_primitive_var] + double(1)/double(2)*(
                    V[1][idx_primitive_var]*V[1][idx_primitive_var] +
                    V[2][idx_primitive_var]*V[2][idx_primitive_var]));
            }
        }
        
        /*
         * Convert primitive variables to conservative variables in the y-direction.
         */
        
        count_eqn = 0;
        for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
        {
            int depth = conservative_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                Q[count_eqn] = conservative_variables[vi]->getPointer(1, di);
                count_eqn++;
            }
        }
        
        count_eqn = 0;
        for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
        {
            int depth = primitive_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                V[count_eqn] = primitive_variables[vi]->getPointer(1, di);
                count_eqn++;
            }
        }
        
        rho     = data_density->getPointer(1, 0);
        epsilon = data_internal_energy->getPointer(1, 0);
        p       = data_pressure->getPointer(1, 0);
        
        // Get the density.
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
                const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var;
                
                rho[idx_primitive_var] = V[0][idx_primitive_var];
            }
        }
        
        // Get the pressure.
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
                const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var;
                
                p[idx_primitive_var] = V[1 + d_dim.getValue()][idx_primitive_var];
            }
        }
        
        // Compute the specific internal energy.
        d_equation_of_state_mixing_rules->getEquationOfState()->computeInternalEnergy(
            data_internal_energy,
            data_density,
            data_pressure,
            thermo_properties_ptr,
            1);
        
        // Set the density.
        for (int j = -num_ghosts_1_primitive_var;
            j < interior_dim_1 + 1 + num_ghosts_1_primitive_var;
            j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                    (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var;
                
                const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var;
                
                Q[0][idx_conservative_var] = V[0][idx_primitive_var];
            }
        }
        
        // Set the momentum.
        for (int j = -num_ghosts_1_primitive_var;
            j < interior_dim_1 + 1 + num_ghosts_1_primitive_var;
            j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                    (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var;
                
                const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var;
                
                Q[1][idx_conservative_var] = rho[idx_primitive_var]*V[1][idx_primitive_var];
                Q[2][idx_conservative_var] = rho[idx_primitive_var]*V[2][idx_primitive_var];
            }
        }
        
        // Set the total energy.
        for (int j = -num_ghosts_1_primitive_var;
            j < interior_dim_1 + 1 + num_ghosts_1_primitive_var;
            j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                    (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var;
                
                const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var;
                
                Q[3][idx_conservative_var] = rho[idx_primitive_var]*
                    (epsilon[idx_primitive_var] + double(1)/double(2)*(
                    V[1][idx_primitive_var]*V[1][idx_primitive_var] +
                    V[2][idx_primitive_var]*V[2][idx_primitive_var]));
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
        
        const int num_ghosts_0_primitive_var = num_ghosts_primitive_var[0];
        const int num_ghosts_1_primitive_var = num_ghosts_primitive_var[1];
        const int num_ghosts_2_primitive_var = num_ghosts_primitive_var[2];
        const int ghostcell_dim_0_primitive_var = ghostcell_dims_primitive_var[0];
        const int ghostcell_dim_1_primitive_var = ghostcell_dims_primitive_var[1];
        
        /*
         * Convert primitive variables to conservative variables in the x-direction.
         */
        
        count_eqn = 0;
        for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
        {
            int depth = conservative_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                Q[count_eqn] = conservative_variables[vi]->getPointer(0, di);
                count_eqn++;
            }
        }
        
        count_eqn = 0;
        for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
        {
            int depth = primitive_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                V[count_eqn] = primitive_variables[vi]->getPointer(0, di);
                count_eqn++;
            }
        }
        
        rho     = data_density->getPointer(0, 0);
        epsilon = data_internal_energy->getPointer(0, 0);
        p       = data_pressure->getPointer(0, 0);
        
        // Get the density.
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
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1) +
                        (k + num_ghosts_2_primitive_var)*(ghostcell_dim_0_primitive_var + 1)*
                            ghostcell_dim_1_primitive_var;
                    
                    rho[idx_primitive_var] = V[0][idx_primitive_var];
                }
            }
        }
        
        // Get the pressure.
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
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1) +
                        (k + num_ghosts_2_primitive_var)*(ghostcell_dim_0_primitive_var + 1)*
                            ghostcell_dim_1_primitive_var;
                    
                    p[idx_primitive_var] = V[1 + d_dim.getValue()][idx_primitive_var];
                }
            }
        }
        
        // Compute the specific internal energy.
        d_equation_of_state_mixing_rules->getEquationOfState()->computeInternalEnergy(
            data_internal_energy,
            data_density,
            data_pressure,
            thermo_properties_ptr,
            0);
        
        // Set the density.
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
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1) +
                        (k + num_ghosts_2_conservative_var)*(ghostcell_dim_0_conservative_var + 1)*
                            ghostcell_dim_1_conservative_var;
                    
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1) +
                        (k + num_ghosts_2_primitive_var)*(ghostcell_dim_0_primitive_var + 1)*
                            ghostcell_dim_1_primitive_var;
                    
                    Q[0][idx_conservative_var] = V[0][idx_primitive_var];
                }
            }
        }
        
        // Set the momentum.
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
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1) +
                        (k + num_ghosts_2_conservative_var)*(ghostcell_dim_0_conservative_var + 1)*
                            ghostcell_dim_1_conservative_var;
                    
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1) +
                        (k + num_ghosts_2_primitive_var)*(ghostcell_dim_0_primitive_var + 1)*
                            ghostcell_dim_1_primitive_var;
                    
                    Q[1][idx_conservative_var] = rho[idx_primitive_var]*V[1][idx_primitive_var];
                    Q[2][idx_conservative_var] = rho[idx_primitive_var]*V[2][idx_primitive_var];
                    Q[3][idx_conservative_var] = rho[idx_primitive_var]*V[3][idx_primitive_var];
                }
            }
        }
        
        // Set the total energy.
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
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1) +
                        (k + num_ghosts_2_conservative_var)*(ghostcell_dim_0_conservative_var + 1)*
                            ghostcell_dim_1_conservative_var;
                    
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1) +
                        (k + num_ghosts_2_primitive_var)*(ghostcell_dim_0_primitive_var + 1)*
                            ghostcell_dim_1_primitive_var;
                    
                    Q[4][idx_conservative_var] = rho[idx_primitive_var]*
                        (epsilon[idx_primitive_var] + double(1)/double(2)*(
                        V[1][idx_primitive_var]*V[1][idx_primitive_var] + 
                        V[2][idx_primitive_var]*V[2][idx_primitive_var] +
                        V[3][idx_primitive_var]*V[3][idx_primitive_var]));
                }
            }
        }
        
        /*
         * Convert primitive variables to conservative variables in the y-direction.
         */
        
        count_eqn = 0;
        for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
        {
            int depth = conservative_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                Q[count_eqn] = conservative_variables[vi]->getPointer(1, di);
                count_eqn++;
            }
        }
        
        count_eqn = 0;
        for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
        {
            int depth = primitive_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                V[count_eqn] = primitive_variables[vi]->getPointer(1, di);
                count_eqn++;
            }
        }
        
        rho     = data_density->getPointer(1, 0);
        epsilon = data_internal_energy->getPointer(1, 0);
        p       = data_pressure->getPointer(1, 0);
        
        // Get the density.
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
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                        (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                            (ghostcell_dim_1_primitive_var + 1);
                    
                    rho[idx_primitive_var] = V[0][idx_primitive_var];
                }
            }
        }
        
        // Get the pressure.
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
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                        (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                            (ghostcell_dim_1_primitive_var + 1);
                    
                    p[idx_primitive_var] = V[1 + d_dim.getValue()][idx_primitive_var];
                }
            }
        }
        
        // Compute the specific internal energy.
        d_equation_of_state_mixing_rules->getEquationOfState()->computeInternalEnergy(
            data_internal_energy,
            data_density,
            data_pressure,
            thermo_properties_ptr,
            1);
        
        // Set the density.
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
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                        (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                            (ghostcell_dim_1_conservative_var + 1);
                    
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                        (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                            (ghostcell_dim_1_primitive_var + 1);
                    
                    Q[0][idx_conservative_var] = V[0][idx_primitive_var];
                }
            }
        }
        
        // Set the momentum.
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
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                        (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                            (ghostcell_dim_1_conservative_var + 1);
                    
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                        (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                            (ghostcell_dim_1_primitive_var + 1);
                    
                    Q[1][idx_conservative_var] = rho[idx_primitive_var]*V[1][idx_primitive_var];
                    Q[2][idx_conservative_var] = rho[idx_primitive_var]*V[2][idx_primitive_var];
                    Q[3][idx_conservative_var] = rho[idx_primitive_var]*V[3][idx_primitive_var];
                }
            }
        }
        
        // Set the total energy.
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
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                        (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                            (ghostcell_dim_1_conservative_var + 1);
                    
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                        (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                            (ghostcell_dim_1_primitive_var + 1);
                    
                    Q[4][idx_conservative_var] = rho[idx_primitive_var]*
                        (epsilon[idx_primitive_var] + double(1)/double(2)*(
                        V[1][idx_primitive_var]*V[1][idx_primitive_var] + 
                        V[2][idx_primitive_var]*V[2][idx_primitive_var] +
                        V[3][idx_primitive_var]*V[3][idx_primitive_var]));
                }
            }
        }
        
        /*
         * Convert primitive variables to conservative variables in the z-direction.
         */
        
        count_eqn = 0;
        for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
        {
            int depth = conservative_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                Q[count_eqn] = conservative_variables[vi]->getPointer(2, di);
                count_eqn++;
            }
        }
        
        count_eqn = 0;
        for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
        {
            int depth = primitive_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                V[count_eqn] = primitive_variables[vi]->getPointer(2, di);
                count_eqn++;
            }
        }
        
        rho     = data_density->getPointer(2, 0);
        epsilon = data_internal_energy->getPointer(2, 0);
        p       = data_pressure->getPointer(2, 0);
        
        // Get the density.
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
                        const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                            (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                            (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                                ghostcell_dim_1_primitive_var;
                        
                        rho[idx_primitive_var] = V[0][idx_primitive_var];
                    }
                }
            }
        }
        
        // Get the pressure.
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
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                        (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                            ghostcell_dim_1_primitive_var;
                    
                    p[idx_primitive_var] = V[1 + d_dim.getValue()][idx_primitive_var];
                }
            }
        }
        
        // Compute the specific internal energy.
        d_equation_of_state_mixing_rules->getEquationOfState()->computeInternalEnergy(
            data_internal_energy,
            data_density,
            data_pressure,
            thermo_properties_ptr,
            2);
        
        // Set the density.
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
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                        (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                            ghostcell_dim_1_conservative_var;
                    
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                        (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                            ghostcell_dim_1_primitive_var;
                    
                    Q[0][idx_conservative_var] = V[0][idx_primitive_var];
                }
            }
        }
        
        // Set the momentum.
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
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                        (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                            ghostcell_dim_1_conservative_var;
                    
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                        (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                            ghostcell_dim_1_primitive_var;
                    
                    Q[1][idx_conservative_var] = rho[idx_primitive_var]*V[1][idx_primitive_var];
                    Q[2][idx_conservative_var] = rho[idx_primitive_var]*V[2][idx_primitive_var];
                    Q[3][idx_conservative_var] = rho[idx_primitive_var]*V[3][idx_primitive_var];
                }
            }
        }
        
        // Set the total energy.
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
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                        (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                            ghostcell_dim_1_conservative_var;
                    
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                        (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                            ghostcell_dim_1_primitive_var;
                    
                    Q[4][idx_conservative_var] = rho[idx_primitive_var]*
                        (epsilon[idx_primitive_var] + double(1)/double(2)*(
                        V[1][idx_primitive_var]*V[1][idx_primitive_var] + 
                        V[2][idx_primitive_var]*V[2][idx_primitive_var] +
                        V[3][idx_primitive_var]*V[3][idx_primitive_var]));
                }
            }
        }
    }
}


/*
 * Convert conservative variables to primitive variables.
 */
void
FlowModelSingleSpecies::convertConservativeVariablesToPrimitiveVariables(
    const std::vector<const double*>& conservative_variables,
    const std::vector<double*>& primitive_variables)
{
    const std::vector<const double*>& Q = conservative_variables;
    const std::vector<double*>&       V = primitive_variables;
    
    // Get the pointers to the momentum components.
    std::vector<const double*> m_ptr;
    m_ptr.reserve(d_dim.getValue());
    for (int di = 0; di < d_dim.getValue(); di++)
    {
        m_ptr.push_back(Q[1 + di]);
    }
    
    // Get the thermodynamic properties of the species.
    std::vector<const double*> thermo_properties_ptr;
    thermo_properties_ptr.reserve(static_cast<int> (d_thermo_properties.size()));
    for (int ti = 0; ti < static_cast<int> (d_thermo_properties.size()); ti++)
    {
        thermo_properties_ptr.push_back(&d_thermo_properties[ti]);
    }
    
    // Compute the specific internal energy.
    double epsilon = double(0);
    if (d_dim == tbox::Dimension(1))
    {
        epsilon = (*Q[2] - double(1)/double(2)*(*Q[1])*(*Q[1])/(*Q[0]))/(*Q[0]);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        epsilon = (*Q[3] - double(1)/double(2)*((*Q[1])*(*Q[1]) + (*Q[2])*(*Q[2]))/(*Q[0]))/(*Q[0]);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        epsilon = (*Q[4] - double(1)/double(2)*((*Q[1])*(*Q[1]) + (*Q[2])*(*Q[2]) + (*Q[3])*(*Q[3]))/(*Q[0]))/(*Q[0]);
    }
    
    // Compute the pressure.
    const double p = d_equation_of_state_mixing_rules->getEquationOfState()->
        getPressure(
            Q[0],
            &epsilon,
            thermo_properties_ptr);
    
    // Convert the conservative variables to primitive variables.
    *V[0] = *Q[0];
    for (int di = 0; di < d_dim.getValue(); di++)
    {
        *V[1 + di] = (*Q[1 + di])/(*Q[0]);
    }
    *V[1 + d_dim.getValue()] = p;
}


/*
 * Convert primitive variables to conservative variables.
 */
void
FlowModelSingleSpecies::convertPrimitiveVariablesToConservativeVariables(
    const std::vector<const double*>& primitive_variables,
    const std::vector<double*>& conservative_variables)
{
    const std::vector<const double*>& V = primitive_variables;
    const std::vector<double*>&       Q = conservative_variables;
    
    // Get the thermodynamic properties of the species.
    std::vector<const double*> thermo_properties_ptr;
    thermo_properties_ptr.reserve(static_cast<int> (d_thermo_properties.size()));
    for (int ti = 0; ti < static_cast<int> (d_thermo_properties.size()); ti++)
    {
        thermo_properties_ptr.push_back(&d_thermo_properties[ti]);
    }
    
    // Compute the total energy.
    const double epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
        getInternalEnergy(
            V[0],
            V[1 + d_dim.getValue()],
            thermo_properties_ptr);
    
    double E = double(0);
    if (d_dim == tbox::Dimension(1))
    {
        E = (*V[0])*(epsilon + double(1)/double(2)*(*V[1])*(*V[1]));
    }
    else if (d_dim == tbox::Dimension(2))
    {
        E = (*V[0])*(epsilon + double(1)/double(2)*((*V[1])*(*V[1]) + (*V[2])*(*V[2])));
    }
    else if (d_dim == tbox::Dimension(3))
    {
        E = (*V[0])*(epsilon + double(1)/double(2)*((*V[1])*(*V[1]) + (*V[2])*(*V[2]) + (*V[3])*(*V[3])));
    }
    
    // Convert the primitive variables to conservative variables.
    *Q[0] = *V[0];
    for (int di = 0; di < d_dim.getValue(); di++)
    {
        *Q[1 + di] = (*V[0])*(*V[1 + di]);
    }
    *Q[1 + d_dim.getValue()] = E;
}


/*
 * Get the variables for the derivatives in the diffusive fluxes.
 */
void
FlowModelSingleSpecies::getDiffusiveFluxVariablesForDerivative(
    std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& derivative_var_data,
    std::vector<std::vector<int> >& derivative_var_component_idx,
    const DIRECTION::TYPE& flux_direction,
    const DIRECTION::TYPE& derivative_direction)
{
    // Create empty box.
    const hier::Box empty_box(d_dim);
    
    derivative_var_data.resize(d_num_eqn);
    derivative_var_component_idx.resize(d_num_eqn);
    
    if (!d_data_velocity)
    {
        computeCellDataOfVelocity(empty_box);
    }
    
    if (!d_data_temperature)
    {
        computeCellDataOfTemperatureWithPressure(empty_box);
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
                         * Mass equation.
                         */
                        
                        derivative_var_data[0].resize(0);
                        derivative_var_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[1].resize(1);
                        derivative_var_component_idx[1].resize(1);
                        
                        // Variable u.
                        derivative_var_data[1][0] = d_data_velocity;
                        derivative_var_component_idx[1][0] = 0;
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[2].resize(2);
                        derivative_var_component_idx[2].resize(2);
                        
                        // Variable u.
                        derivative_var_data[2][0] = d_data_velocity;
                        derivative_var_component_idx[2][0] = 0;
                        
                        // Variable T.
                        derivative_var_data[2][1] = d_data_temperature;
                        derivative_var_component_idx[2][1] = 0;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelSingleSpecies::getDiffusiveFluxVariablesForDerivative()\n"
                            << "There are only x-direction for one-dimensional problem."
                            << std::endl);
                    }
                }
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::getDiffusiveFluxVariablesForDerivative()\n"
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
                         * Mass equation.
                         */
                        
                        derivative_var_data[0].resize(0);
                        derivative_var_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[1].resize(1);
                        derivative_var_component_idx[1].resize(1);
                        
                        // Variable u.
                        derivative_var_data[1][0] = d_data_velocity;
                        derivative_var_component_idx[1][0] = 0;
                        
                        derivative_var_data[2].resize(1);
                        derivative_var_component_idx[2].resize(1);
                        
                        // Variable v.
                        derivative_var_data[2][0] = d_data_velocity;
                        derivative_var_component_idx[2][0] = 1;
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[3].resize(3);
                        derivative_var_component_idx[3].resize(3);
                        
                        // Variable u.
                        derivative_var_data[3][0] = d_data_velocity;
                        derivative_var_component_idx[3][0] = 0;
                        
                        // Variable v.
                        derivative_var_data[3][1] = d_data_velocity;
                        derivative_var_component_idx[3][1] = 1;
                        
                        // Variable T.
                        derivative_var_data[3][2] = d_data_temperature;
                        derivative_var_component_idx[3][2] = 0;
                        
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        /*
                         * Mass equation.
                         */
                        
                        derivative_var_data[0].resize(0);
                        derivative_var_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[1].resize(1);
                        derivative_var_component_idx[1].resize(1);
                        
                        // Variable v.
                        derivative_var_data[1][0] = d_data_velocity;
                        derivative_var_component_idx[1][0] = 1;
                        
                        derivative_var_data[2].resize(1);
                        derivative_var_component_idx[2].resize(1);
                        
                        // Variable u.
                        derivative_var_data[2][0] = d_data_velocity;
                        derivative_var_component_idx[2][0] = 0;
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[3].resize(2);
                        derivative_var_component_idx[3].resize(2);
                        
                        // Variable u.
                        derivative_var_data[3][0] = d_data_velocity;
                        derivative_var_component_idx[3][0] = 0;
                        
                        // Variable v.
                        derivative_var_data[3][1] = d_data_velocity;
                        derivative_var_component_idx[3][1] = 1;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelSingleSpecies::getDiffusiveFluxVariablesForDerivative()\n"
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
                         * Mass equation.
                         */
                        
                        derivative_var_data[0].resize(0);
                        derivative_var_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[1].resize(1);
                        derivative_var_component_idx[1].resize(1);
                        
                        // Variable v.
                        derivative_var_data[1][0] = d_data_velocity;
                        derivative_var_component_idx[1][0] = 1;
                        
                        derivative_var_data[2].resize(1);
                        derivative_var_component_idx[2].resize(1);
                        
                        // Variable u.
                        derivative_var_data[2][0] = d_data_velocity;
                        derivative_var_component_idx[2][0] = 0;
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[3].resize(2);
                        derivative_var_component_idx[3].resize(2);
                        
                        // Variable u.
                        derivative_var_data[3][0] = d_data_velocity;
                        derivative_var_component_idx[3][0] = 0;
                        
                        // Variable v.
                        derivative_var_data[3][1] = d_data_velocity;
                        derivative_var_component_idx[3][1] = 1;
                        
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        /*
                         * Mass equation.
                         */
                        
                        derivative_var_data[0].resize(0);
                        derivative_var_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[1].resize(1);
                        derivative_var_component_idx[1].resize(1);
                        
                        // Variable u.
                        derivative_var_data[1][0] = d_data_velocity;
                        derivative_var_component_idx[1][0] = 0;
                        
                        derivative_var_data[2].resize(1);
                        derivative_var_component_idx[2].resize(1);
                        
                        // Variable v.
                        derivative_var_data[2][0] = d_data_velocity;
                        derivative_var_component_idx[2][0] = 1;
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[3].resize(3);
                        derivative_var_component_idx[3].resize(3);
                        
                        // Variable u.
                        derivative_var_data[3][0] = d_data_velocity;
                        derivative_var_component_idx[3][0] = 0;
                        
                        // Variable v.
                        derivative_var_data[3][1] = d_data_velocity;
                        derivative_var_component_idx[3][1] = 1;
                        
                        // Variable T.
                        derivative_var_data[3][2] = d_data_temperature;
                        derivative_var_component_idx[3][2] = 0;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelSingleSpecies::getDiffusiveFluxVariablesForDerivative()\n"
                            << "There are only x-direction and y-direction for two-dimensional problem."
                            << std::endl);
                    }
                }
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::getDiffusiveFluxVariablesForDerivative()\n"
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
                         * Mass equation.
                         */
                        
                        derivative_var_data[0].resize(0);
                        derivative_var_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[1].resize(1);
                        derivative_var_component_idx[1].resize(1);
                        
                        // Variable u.
                        derivative_var_data[1][0] = d_data_velocity;
                        derivative_var_component_idx[1][0] = 0;
                        
                        derivative_var_data[2].resize(1);
                        derivative_var_component_idx[2].resize(1);
                        
                        // Variable v.
                        derivative_var_data[2][0] = d_data_velocity;
                        derivative_var_component_idx[2][0] = 1;
                        
                        derivative_var_data[3].resize(1);
                        derivative_var_component_idx[3].resize(1);
                        
                        // Variable w.
                        derivative_var_data[3][0] = d_data_velocity;
                        derivative_var_component_idx[3][0] = 2;
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[4].resize(4);
                        derivative_var_component_idx[4].resize(4);
                        
                        // Variable u.
                        derivative_var_data[4][0] = d_data_velocity;
                        derivative_var_component_idx[4][0] = 0;
                        
                        // Variable v.
                        derivative_var_data[4][1] = d_data_velocity;
                        derivative_var_component_idx[4][1] = 1;
                        
                        // Variable w.
                        derivative_var_data[4][2] = d_data_velocity;
                        derivative_var_component_idx[4][2] = 2;
                        
                        // Variable T.
                        derivative_var_data[4][3] = d_data_temperature;
                        derivative_var_component_idx[4][3] = 0;
                        
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        /*
                         * Mass equation.
                         */
                        
                        derivative_var_data[0].resize(0);
                        derivative_var_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[1].resize(1);
                        derivative_var_component_idx[1].resize(1);
                        
                        // Variable v.
                        derivative_var_data[1][0] = d_data_velocity;
                        derivative_var_component_idx[1][0] = 1;
                        
                        derivative_var_data[2].resize(1);
                        derivative_var_component_idx[2].resize(1);
                        
                        // Variable u.
                        derivative_var_data[2][0] = d_data_velocity;
                        derivative_var_component_idx[2][0] = 0;
                        
                        derivative_var_data[3].resize(0);
                        derivative_var_component_idx[3].resize(0);
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[4].resize(2);
                        derivative_var_component_idx[4].resize(2);
                        
                        // Variable u.
                        derivative_var_data[4][0] = d_data_velocity;
                        derivative_var_component_idx[4][0] = 0;
                        
                        // Variable v.
                        derivative_var_data[4][1] = d_data_velocity;
                        derivative_var_component_idx[4][1] = 1;
                        
                        break;
                    }
                    case DIRECTION::Z_DIRECTION:
                    {
                        /*
                         * Mass equation.
                         */
                        
                        derivative_var_data[0].resize(0);
                        derivative_var_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[1].resize(1);
                        derivative_var_component_idx[1].resize(1);
                        
                        // Variable w.
                        derivative_var_data[1][0] = d_data_velocity;
                        derivative_var_component_idx[1][0] = 2;
                        
                        derivative_var_data[2].resize(0);
                        derivative_var_component_idx[2].resize(0);
                        
                        derivative_var_data[3].resize(1);
                        derivative_var_component_idx[3].resize(1);
                        
                        // Variable u.
                        derivative_var_data[3][0] = d_data_velocity;
                        derivative_var_component_idx[3][0] = 0;
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[4].resize(2);
                        derivative_var_component_idx[4].resize(2);
                        
                        // Variable u.
                        derivative_var_data[4][0] = d_data_velocity;
                        derivative_var_component_idx[4][0] = 0;
                        
                        // Variable w.
                        derivative_var_data[4][1] = d_data_velocity;
                        derivative_var_component_idx[4][1] = 2;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelSingleSpecies::getDiffusiveFluxVariablesForDerivative()\n"
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
                         * Mass equation.
                         */
                        
                        derivative_var_data[0].resize(0);
                        derivative_var_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[1].resize(1);
                        derivative_var_component_idx[1].resize(1);
                        
                        // Variable v.
                        derivative_var_data[1][0] = d_data_velocity;
                        derivative_var_component_idx[1][0] = 1;
                        
                        derivative_var_data[2].resize(1);
                        derivative_var_component_idx[2].resize(1);
                        
                        // Variable u.
                        derivative_var_data[2][0] = d_data_velocity;
                        derivative_var_component_idx[2][0] = 0;
                        
                        derivative_var_data[3].resize(0);
                        derivative_var_component_idx[3].resize(0);
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[4].resize(2);
                        derivative_var_component_idx[4].resize(2);
                        
                        // Variable u.
                        derivative_var_data[4][0] = d_data_velocity;
                        derivative_var_component_idx[4][0] = 0;
                        
                        // Variable v.
                        derivative_var_data[4][1] = d_data_velocity;
                        derivative_var_component_idx[4][1] = 1;
                        
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        /*
                         * Mass equation.
                         */
                        
                        derivative_var_data[0].resize(0);
                        derivative_var_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[1].resize(1);
                        derivative_var_component_idx[1].resize(1);
                        
                        // Variable u.
                        derivative_var_data[1][0] = d_data_velocity;
                        derivative_var_component_idx[1][0] = 0;
                        
                        derivative_var_data[2].resize(1);
                        derivative_var_component_idx[2].resize(1);
                        
                        // Variable v.
                        derivative_var_data[2][0] = d_data_velocity;
                        derivative_var_component_idx[2][0] = 1;
                        
                        derivative_var_data[3].resize(1);
                        derivative_var_component_idx[3].resize(1);
                        
                        // Variable w.
                        derivative_var_data[3][0] = d_data_velocity;
                        derivative_var_component_idx[3][0] = 2;
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[4].resize(4);
                        derivative_var_component_idx[4].resize(4);
                        
                        // Variable u.
                        derivative_var_data[4][0] = d_data_velocity;
                        derivative_var_component_idx[4][0] = 0;
                        
                        // Variable v.
                        derivative_var_data[4][1] = d_data_velocity;
                        derivative_var_component_idx[4][1] = 1;
                        
                        // Variable w.
                        derivative_var_data[4][2] = d_data_velocity;
                        derivative_var_component_idx[4][2] = 2;
                        
                        // Variable T.
                        derivative_var_data[4][3] = d_data_temperature;
                        derivative_var_component_idx[4][3] = 0;
                        
                        break;
                    }
                    case DIRECTION::Z_DIRECTION:
                    {
                        /*
                         * Mass equation.
                         */
                        
                        derivative_var_data[0].resize(0);
                        derivative_var_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[1].resize(0);
                        derivative_var_component_idx[1].resize(0);
                        
                        derivative_var_data[2].resize(1);
                        derivative_var_component_idx[2].resize(1);
                        
                        // Variable w.
                        derivative_var_data[2][0] = d_data_velocity;
                        derivative_var_component_idx[2][0] = 2;
                        
                        derivative_var_data[3].resize(1);
                        derivative_var_component_idx[3].resize(1);
                        
                        // Variable v.
                        derivative_var_data[3][0] = d_data_velocity;
                        derivative_var_component_idx[3][0] = 1;
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[4].resize(2);
                        derivative_var_component_idx[4].resize(2);
                        
                        // Variable v.
                        derivative_var_data[4][0] = d_data_velocity;
                        derivative_var_component_idx[4][0] = 1;
                        
                        // Variable w.
                        derivative_var_data[4][1] = d_data_velocity;
                        derivative_var_component_idx[4][1] = 2;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelSingleSpecies::getDiffusiveFluxVariablesForDerivative()\n"
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
                         * Mass equation.
                         */
                        
                        derivative_var_data[0].resize(0);
                        derivative_var_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[1].resize(1);
                        derivative_var_component_idx[1].resize(1);
                        
                        // Variable w.
                        derivative_var_data[1][0] = d_data_velocity;
                        derivative_var_component_idx[1][0] = 2;
                        
                        derivative_var_data[2].resize(0);
                        derivative_var_component_idx[2].resize(0);
                        
                        derivative_var_data[3].resize(1);
                        derivative_var_component_idx[3].resize(1);
                        
                        // Variable u.
                        derivative_var_data[3][0] = d_data_velocity;
                        derivative_var_component_idx[3][0] = 0;
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[4].resize(2);
                        derivative_var_component_idx[4].resize(2);
                        
                        // Variable u.
                        derivative_var_data[4][0] = d_data_velocity;
                        derivative_var_component_idx[4][0] = 0;
                        
                        // Variable w.
                        derivative_var_data[4][1] = d_data_velocity;
                        derivative_var_component_idx[4][1] = 2;
                        
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        /*
                         * Mass equation.
                         */
                        
                        derivative_var_data[0].resize(0);
                        derivative_var_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[1].resize(0);
                        derivative_var_component_idx[1].resize(0);
                        
                        derivative_var_data[2].resize(1);
                        derivative_var_component_idx[2].resize(1);
                        
                        // Variable w.
                        derivative_var_data[2][0] = d_data_velocity;
                        derivative_var_component_idx[2][0] = 2;
                        
                        derivative_var_data[3].resize(1);
                        derivative_var_component_idx[3].resize(1);
                        
                        // Variable v.
                        derivative_var_data[3][0] = d_data_velocity;
                        derivative_var_component_idx[3][0] = 1;
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[4].resize(2);
                        derivative_var_component_idx[4].resize(2);
                        
                        // Variable v.
                        derivative_var_data[4][0] = d_data_velocity;
                        derivative_var_component_idx[4][0] = 1;
                        
                        // Variable w.
                        derivative_var_data[4][1] = d_data_velocity;
                        derivative_var_component_idx[4][1] = 2;
                        
                        break;
                    }
                    case DIRECTION::Z_DIRECTION:
                    {
                        /*
                         * Mass equation.
                         */
                        
                        derivative_var_data[0].resize(0);
                        derivative_var_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        derivative_var_data[1].resize(1);
                        derivative_var_component_idx[1].resize(1);
                        
                        // Variable u.
                        derivative_var_data[1][0] = d_data_velocity;
                        derivative_var_component_idx[1][0] = 0;
                        
                        derivative_var_data[2].resize(1);
                        derivative_var_component_idx[2].resize(1);
                        
                        // Variable v.
                        derivative_var_data[2][0] = d_data_velocity;
                        derivative_var_component_idx[2][0] = 1;
                        
                        derivative_var_data[3].resize(1);
                        derivative_var_component_idx[3].resize(1);
                        
                        // Variable w.
                        derivative_var_data[3][0] = d_data_velocity;
                        derivative_var_component_idx[3][0] = 2;
                        
                        /*
                         * Energy equation.
                         */
                        
                        derivative_var_data[4].resize(4);
                        derivative_var_component_idx[4].resize(4);
                        
                        // Variable u.
                        derivative_var_data[4][0] = d_data_velocity;
                        derivative_var_component_idx[4][0] = 0;
                        
                        // Variable v.
                        derivative_var_data[4][1] = d_data_velocity;
                        derivative_var_component_idx[4][1] = 1;
                        
                        // Variable w.
                        derivative_var_data[4][2] = d_data_velocity;
                        derivative_var_component_idx[4][2] = 2;
                        
                        // Variable T.
                        derivative_var_data[4][3] = d_data_temperature;
                        derivative_var_component_idx[4][3] = 0;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelSingleSpecies::getDiffusiveFluxVariablesForDerivative()\n"
                            << "There are only x-direction, y-direction and z-direction for three-dimensional problem."
                            << std::endl);
                    }
                }
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::getDiffusiveFluxVariablesForDerivative()\n"
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
FlowModelSingleSpecies::getDiffusiveFluxDiffusivities(
    std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& diffusivities_data,
    std::vector<std::vector<int> >& diffusivities_component_idx,
    const DIRECTION::TYPE& flux_direction,
    const DIRECTION::TYPE& derivative_direction)
{
    // Create empty box.
    const hier::Box empty_box(d_dim);
    
    if (!d_equation_of_shear_viscosity_mixing_rules ||
        !d_equation_of_bulk_viscosity_mixing_rules ||
        !d_equation_of_thermal_conductivity_mixing_rules)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::getDiffusiveFluxDiffusivities()\n"
            << "Either mixing rule of shear viscosity, bulk viscosity or"
            << " thermal conductivity is not initialized."
            << std::endl);
    }
    
    diffusivities_data.resize(d_num_eqn);
    diffusivities_component_idx.resize(d_num_eqn);
    
    if (!d_data_diffusivities)
    {
        if (!d_data_velocity)
        {
            computeCellDataOfVelocity(empty_box);
        }
        
        if (!d_data_pressure)
        {
            computeCellDataOfPressureWithInternalEnergy(empty_box);
        }
        
        if (!d_data_temperature)
        {
            computeCellDataOfTemperatureWithPressure(empty_box);
        }
        
        /*
         * Create temporary cell data of shear viscosity, bulk viscosity and thermal conductivity.
         */
        
        boost::shared_ptr<pdat::CellData<double> > data_shear_viscosity(
            new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_diffusivities));
        
        boost::shared_ptr<pdat::CellData<double> > data_bulk_viscosity(
            new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_diffusivities));
        
        boost::shared_ptr<pdat::CellData<double> > data_thermal_conductivity(
            new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_diffusivities));
        
        // Get the pointers to the cell data of shear viscosity, bulk viscosity and thermal conductivity.
        double* mu    = data_shear_viscosity->getPointer(0);
        double* mu_v  = data_bulk_viscosity->getPointer(0);
        double* kappa = data_thermal_conductivity->getPointer(0);
        
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
        
        // Compute the shear viscosity field.
        d_equation_of_shear_viscosity_mixing_rules->getEquationOfShearViscosity()->
            computeShearViscosity(
                data_shear_viscosity,
                d_data_pressure,
                d_data_temperature,
                molecular_properties_shear_viscosity_ptr,
                empty_box);
        
        // Compute the bulk viscosity field.
        d_equation_of_bulk_viscosity_mixing_rules->getEquationOfBulkViscosity()->
            computeBulkViscosity(
                data_bulk_viscosity,
                d_data_pressure,
                d_data_temperature,
                molecular_properties_bulk_viscosity_ptr,
                empty_box);
        
        // Compute the thermal conductivity field.
        d_equation_of_thermal_conductivity_mixing_rules->getEquationOfThermalConductivity()->
            computeThermalConductivity(
                data_thermal_conductivity,
                d_data_pressure,
                d_data_temperature,
                molecular_properties_thermal_conductivity_ptr,
                empty_box);
        
        if (d_dim == tbox::Dimension(1))
        {
            d_data_diffusivities.reset(new pdat::CellData<double>(
                d_interior_box, 3, d_num_subghosts_diffusivities));
            
            // Get the pointer to cell data of velocity and diffusivities.
            double* u = d_data_velocity->getPointer(0);
            
            double* D_00 = d_data_diffusivities->getPointer(0);
            double* D_01 = d_data_diffusivities->getPointer(1);
            double* D_02 = d_data_diffusivities->getPointer(2);
            
            /*
             * Compute the diffusivities.
             */
            for (int i = -d_num_subghosts_diffusivities[0];
                 i < d_interior_dims[0] + d_num_subghosts_diffusivities[0];
                 i++)
            {
                // Compute the linear indices.
                const int idx_diffusivities = i + d_num_subghosts_diffusivities[0];
                const int idx_velocity = i + d_num_subghosts_velocity[0];
                
                D_00[idx_diffusivities] = -(double(4)/double(3)*mu[idx_diffusivities] + mu_v[idx_diffusivities]);
                D_01[idx_diffusivities] = -u[idx_velocity]*(double(4)/double(3)*mu[idx_diffusivities] +
                    mu_v[idx_diffusivities]);
                D_02[idx_diffusivities] = -kappa[idx_diffusivities];
            }
        }
        else if (d_dim == tbox::Dimension(2))
        {
            d_data_diffusivities.reset(new pdat::CellData<double>(
                d_interior_box, 10, d_num_subghosts_diffusivities));
            
            // Get the pointer to cell data of velocity and diffusivities.
            double* u = d_data_velocity->getPointer(0);
            double* v = d_data_velocity->getPointer(1);
            
            double* D_00 = d_data_diffusivities->getPointer(0);
            double* D_01 = d_data_diffusivities->getPointer(1);
            double* D_02 = d_data_diffusivities->getPointer(2);
            double* D_03 = d_data_diffusivities->getPointer(3);
            double* D_04 = d_data_diffusivities->getPointer(4);
            double* D_05 = d_data_diffusivities->getPointer(5);
            double* D_06 = d_data_diffusivities->getPointer(6);
            double* D_07 = d_data_diffusivities->getPointer(7);
            double* D_08 = d_data_diffusivities->getPointer(8);
            double* D_09 = d_data_diffusivities->getPointer(9);
            
            /*
             * Compute the diffusivities.
             */
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
                    
                    D_00[idx_diffusivities] = -(double(4)/double(3)*mu[idx_diffusivities] + mu_v[idx_diffusivities]);
                    D_01[idx_diffusivities] = double(2)/double(3)*mu[idx_diffusivities] - mu_v[idx_diffusivities];
                    D_02[idx_diffusivities] = -mu[idx_diffusivities];
                    D_03[idx_diffusivities] = -u[idx_velocity]*(double(4)/double(3)*mu[idx_diffusivities] +
                        mu_v[idx_diffusivities]);
                    D_04[idx_diffusivities] = -v[idx_velocity]*(double(4)/double(3)*mu[idx_diffusivities] +
                        mu_v[idx_diffusivities]);
                    D_05[idx_diffusivities] = u[idx_velocity]*(double(2)/double(3)*mu[idx_diffusivities] -
                        mu_v[idx_diffusivities]);
                    D_06[idx_diffusivities] = v[idx_velocity]*(double(2)/double(3)*mu[idx_diffusivities] -
                        mu_v[idx_diffusivities]);
                    D_07[idx_diffusivities] = -u[idx_velocity]*mu[idx_diffusivities];
                    D_08[idx_diffusivities] = -v[idx_velocity]*mu[idx_diffusivities];
                    D_09[idx_diffusivities] = -kappa[idx_diffusivities];
                }
            }
        }
        else if (d_dim == tbox::Dimension(3))
        {
            d_data_diffusivities.reset(new pdat::CellData<double>(
                d_interior_box, 13, d_num_subghosts_diffusivities));
            
            // Get the pointer to cell data of velocity and diffusivities.
            double* u = d_data_velocity->getPointer(0);
            double* v = d_data_velocity->getPointer(1);
            double* w = d_data_velocity->getPointer(2);
            
            double* D_00 = d_data_diffusivities->getPointer(0);
            double* D_01 = d_data_diffusivities->getPointer(1);
            double* D_02 = d_data_diffusivities->getPointer(2);
            double* D_03 = d_data_diffusivities->getPointer(3);
            double* D_04 = d_data_diffusivities->getPointer(4);
            double* D_05 = d_data_diffusivities->getPointer(5);
            double* D_06 = d_data_diffusivities->getPointer(6);
            double* D_07 = d_data_diffusivities->getPointer(7);
            double* D_08 = d_data_diffusivities->getPointer(8);
            double* D_09 = d_data_diffusivities->getPointer(9);
            double* D_10 = d_data_diffusivities->getPointer(10);
            double* D_11 = d_data_diffusivities->getPointer(11);
            double* D_12 = d_data_diffusivities->getPointer(12);
            
            /*
             * Compute the diffusivities.
             */
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
                        const int idx_diffusivities = (i + d_num_subghosts_diffusivities[0]) +
                            (j + d_num_subghosts_diffusivities[1])*d_subghostcell_dims_diffusivities[0] +
                            (k + d_num_subghosts_diffusivities[2])*d_subghostcell_dims_diffusivities[0]*
                                d_subghostcell_dims_diffusivities[1];
                        
                        const int idx_velocity = (i + d_num_subghosts_velocity[0]) +
                            (j + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0] +
                            (k + d_num_subghosts_velocity[2])*d_subghostcell_dims_velocity[0]*
                                d_subghostcell_dims_velocity[1];
                        
                        D_00[idx_diffusivities] = -(double(4)/double(3)*mu[idx_diffusivities] + mu_v[idx_diffusivities]);
                        D_01[idx_diffusivities] = double(2)/double(3)*mu[idx_diffusivities] - mu_v[idx_diffusivities];
                        D_02[idx_diffusivities] = -mu[idx_diffusivities];
                        D_03[idx_diffusivities] = -u[idx_velocity]*(double(4)/double(3)*mu[idx_diffusivities] +
                            mu_v[idx_diffusivities]);
                        D_04[idx_diffusivities] = -v[idx_velocity]*(double(4)/double(3)*mu[idx_diffusivities] +
                            mu_v[idx_diffusivities]);
                        D_05[idx_diffusivities] = -w[idx_velocity]*(double(4)/double(3)*mu[idx_diffusivities] +
                            mu_v[idx_diffusivities]);
                        D_06[idx_diffusivities] = u[idx_velocity]*(double(2)/double(3)*mu[idx_diffusivities] -
                            mu_v[idx_diffusivities]);
                        D_07[idx_diffusivities] = v[idx_velocity]*(double(2)/double(3)*mu[idx_diffusivities] -
                            mu_v[idx_diffusivities]);
                        D_08[idx_diffusivities] = w[idx_velocity]*(double(2)/double(3)*mu[idx_diffusivities] -
                            mu_v[idx_diffusivities]);
                        D_09[idx_diffusivities] = -u[idx_velocity]*mu[idx_diffusivities];
                        D_10[idx_diffusivities] = -v[idx_velocity]*mu[idx_diffusivities];
                        D_11[idx_diffusivities] = -w[idx_velocity]*mu[idx_diffusivities];
                        D_12[idx_diffusivities] = -kappa[idx_diffusivities];
                    }
                }
            }
        }
        
        data_shear_viscosity.reset();
        data_bulk_viscosity.reset();
        data_thermal_conductivity.reset();
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
                         * Mass equation.
                         */
                        
                        diffusivities_data[0].resize(0);
                        diffusivities_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[1].resize(1);
                        diffusivities_component_idx[1].resize(1);
                        
                        // -(4/3*mu + mu_v).
                        diffusivities_data[1][0] = d_data_diffusivities;
                        diffusivities_component_idx[1][0] = 0;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[2].resize(2);
                        diffusivities_component_idx[2].resize(2);
                        
                        // -u*(4/3*mu + mu_v).
                        diffusivities_data[2][0] = d_data_diffusivities;
                        diffusivities_component_idx[2][0] = 1;
                        
                        // -kappa.
                        diffusivities_data[2][1] = d_data_diffusivities;
                        diffusivities_component_idx[2][1] = 2;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelSingleSpecies::getDiffusiveFluxDiffusivities()\n"
                            << "There are only x-direction for one-dimensional problem."
                            << std::endl);
                    }
                }
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::getDiffusiveFluxDiffusivities()\n"
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
                         * Mass equation.
                         */
                        
                        diffusivities_data[0].resize(0);
                        diffusivities_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[1].resize(1);
                        diffusivities_component_idx[1].resize(1);
                        
                        // -(4/3*mu + mu_v).
                        diffusivities_data[1][0] = d_data_diffusivities;
                        diffusivities_component_idx[1][0] = 0;
                        
                        diffusivities_data[2].resize(1);
                        diffusivities_component_idx[2].resize(1);
                        
                        // -mu.
                        diffusivities_data[2][0] = d_data_diffusivities;
                        diffusivities_component_idx[2][0] = 2;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[3].resize(3);
                        diffusivities_component_idx[3].resize(3);
                        
                        // -u*(4/3*mu + mu_v).
                        diffusivities_data[3][0] = d_data_diffusivities;
                        diffusivities_component_idx[3][0] = 3;
                        
                        // -v*mu.
                        diffusivities_data[3][1] = d_data_diffusivities;
                        diffusivities_component_idx[3][1] = 8;
                        
                        // -kappa.
                        diffusivities_data[3][2] = d_data_diffusivities;
                        diffusivities_component_idx[3][2] = 9;
                        
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        /*
                         * Mass equation.
                         */
                        
                        diffusivities_data[0].resize(0);
                        diffusivities_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[1].resize(1);
                        diffusivities_component_idx[1].resize(1);
                        
                        // 2/3*mu - mu_v.
                        diffusivities_data[1][0] = d_data_diffusivities;
                        diffusivities_component_idx[1][0] = 1;
                        
                        diffusivities_data[2].resize(1);
                        diffusivities_component_idx[2].resize(1);
                        
                        // -mu.
                        diffusivities_data[2][0] = d_data_diffusivities;
                        diffusivities_component_idx[2][0] = 2;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[3].resize(2);
                        diffusivities_component_idx[3].resize(2);
                        
                        // -v*mu.
                        diffusivities_data[3][0] = d_data_diffusivities;
                        diffusivities_component_idx[3][0] = 8;
                        
                        // u*(2/3*mu - mu_v).
                        diffusivities_data[3][1] = d_data_diffusivities;
                        diffusivities_component_idx[3][1] = 5;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelSingleSpecies::getDiffusiveFluxDiffusivities()\n"
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
                         * Mass equation.
                         */
                        
                        diffusivities_data[0].resize(0);
                        diffusivities_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[1].resize(1);
                        diffusivities_component_idx[1].resize(1);
                        
                        // -mu.
                        diffusivities_data[1][0] = d_data_diffusivities;
                        diffusivities_component_idx[1][0] = 2;
                        
                        diffusivities_data[2].resize(1);
                        diffusivities_component_idx[2].resize(1);
                        
                        // 2/3*mu - mu_v.
                        diffusivities_data[2][0] = d_data_diffusivities;
                        diffusivities_component_idx[2][0] = 1;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[3].resize(2);
                        diffusivities_component_idx[3].resize(2);
                        
                        // v*(2/3*mu - mu_v).
                        diffusivities_data[3][0] = d_data_diffusivities;
                        diffusivities_component_idx[3][0] = 6;
                        
                        // -u*mu.
                        diffusivities_data[3][1] = d_data_diffusivities;
                        diffusivities_component_idx[3][1] = 7;
                        
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        /*
                         * Mass equation.
                         */
                        
                        diffusivities_data[0].resize(0);
                        diffusivities_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[1].resize(1);
                        diffusivities_component_idx[1].resize(1);
                        
                        // -mu.
                        diffusivities_data[1][0] = d_data_diffusivities;
                        diffusivities_component_idx[1][0] = 2;
                        
                        diffusivities_data[2].resize(1);
                        diffusivities_component_idx[2].resize(1);
                        
                        // -(4/3*mu + mu_v).
                        diffusivities_data[2][0] = d_data_diffusivities;
                        diffusivities_component_idx[2][0] = 0;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[3].resize(3);
                        diffusivities_component_idx[3].resize(3);
                        
                        // -u*mu.
                        diffusivities_data[3][0] = d_data_diffusivities;
                        diffusivities_component_idx[3][0] = 7;
                        
                        // -v*(4/3*mu + mu_v).
                        diffusivities_data[3][1] = d_data_diffusivities;
                        diffusivities_component_idx[3][1] = 4;
                        
                        // -kappa.
                        diffusivities_data[3][2] = d_data_diffusivities;
                        diffusivities_component_idx[3][2] = 9;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelSingleSpecies::getDiffusiveFluxDiffusivities()\n"
                            << "There are only x-direction and y-direction for two-dimensional problem."
                            << std::endl);
                    }
                }
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::getDiffusiveFluxDiffusivities()\n"
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
                         * Mass equation.
                         */
                        
                        diffusivities_data[0].resize(0);
                        diffusivities_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[1].resize(1);
                        diffusivities_component_idx[1].resize(1);
                        
                        // -(4/3*mu + mu_v).
                        diffusivities_data[1][0] = d_data_diffusivities;
                        diffusivities_component_idx[1][0] = 0;
                        
                        diffusivities_data[2].resize(1);
                        diffusivities_component_idx[2].resize(1);
                        
                        // -mu.
                        diffusivities_data[2][0] = d_data_diffusivities;
                        diffusivities_component_idx[2][0] = 2;
                        
                        diffusivities_data[3].resize(1);
                        diffusivities_component_idx[3].resize(1);
                        
                        // -mu.
                        diffusivities_data[3][0] = d_data_diffusivities;
                        diffusivities_component_idx[3][0] = 2;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[4].resize(4);
                        diffusivities_component_idx[4].resize(4);
                        
                        // -u*(4/3*mu + mu_v).
                        diffusivities_data[4][0] = d_data_diffusivities;
                        diffusivities_component_idx[4][0] = 3;
                        
                        // -v*mu.
                        diffusivities_data[4][1] = d_data_diffusivities;
                        diffusivities_component_idx[4][1] = 10;
                        
                        // -w*mu.
                        diffusivities_data[4][2] = d_data_diffusivities;
                        diffusivities_component_idx[4][2] = 11;
                        
                        // -kappa.
                        diffusivities_data[4][3] = d_data_diffusivities;
                        diffusivities_component_idx[4][3] = 12;
                        
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        /*
                         * Mass equation.
                         */
                        
                        diffusivities_data[0].resize(0);
                        diffusivities_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[1].resize(1);
                        diffusivities_component_idx[1].resize(1);
                        
                        // 2/3*mu - mu_v.
                        diffusivities_data[1][0] = d_data_diffusivities;
                        diffusivities_component_idx[1][0] = 1;
                        
                        diffusivities_data[2].resize(1);
                        diffusivities_component_idx[2].resize(1);
                        
                        // -mu.
                        diffusivities_data[2][0] = d_data_diffusivities;
                        diffusivities_component_idx[2][0] = 2;
                        
                        diffusivities_data[3].resize(0);
                        diffusivities_component_idx[3].resize(0);
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[4].resize(2);
                        diffusivities_component_idx[4].resize(2);
                        
                        // -v*mu.
                        diffusivities_data[4][0] = d_data_diffusivities;
                        diffusivities_component_idx[4][0] = 10;
                        
                        // u*(2/3*mu - mu_v).
                        diffusivities_data[4][1] = d_data_diffusivities;
                        diffusivities_component_idx[4][1] = 6;
                        
                        break;
                    }
                    case DIRECTION::Z_DIRECTION:
                    {
                        /*
                         * Mass equation.
                         */
                        
                        diffusivities_data[0].resize(0);
                        diffusivities_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[1].resize(1);
                        diffusivities_component_idx[1].resize(1);
                        
                        // 2/3*mu - mu_v.
                        diffusivities_data[1][0] = d_data_diffusivities;
                        diffusivities_component_idx[1][0] = 1;
                        
                        diffusivities_data[2].resize(0);
                        diffusivities_component_idx[2].resize(0);
                        
                        diffusivities_data[3].resize(1);
                        diffusivities_component_idx[3].resize(1);
                        
                        // -mu.
                        diffusivities_data[3][0] = d_data_diffusivities;
                        diffusivities_component_idx[3][0] = 2;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[4].resize(2);
                        diffusivities_component_idx[4].resize(2);
                        
                        // -w*mu.
                        diffusivities_data[4][0] = d_data_diffusivities;
                        diffusivities_component_idx[4][0] = 11;
                        
                        // u*(2/3*mu - mu_v).
                        diffusivities_data[4][1] = d_data_diffusivities;
                        diffusivities_component_idx[4][1] = 6;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelSingleSpecies::getDiffusiveFluxDiffusivities()\n"
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
                         * Mass equation.
                         */
                        
                        diffusivities_data[0].resize(0);
                        diffusivities_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[1].resize(1);
                        diffusivities_component_idx[1].resize(1);
                        
                        // -mu.
                        diffusivities_data[1][0] = d_data_diffusivities;
                        diffusivities_component_idx[1][0] = 2;
                        
                        diffusivities_data[2].resize(1);
                        diffusivities_component_idx[2].resize(1);
                        
                        // 2/3*mu - mu_v.
                        diffusivities_data[2][0] = d_data_diffusivities;
                        diffusivities_component_idx[2][0] = 1;
                        
                        diffusivities_data[3].resize(0);
                        diffusivities_component_idx[3].resize(0);
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[4].resize(2);
                        diffusivities_component_idx[4].resize(2);
                        
                        // v*(2/3*mu - mu_v).
                        diffusivities_data[4][0] = d_data_diffusivities;
                        diffusivities_component_idx[4][0] = 7;
                        
                        // -u*mu.
                        diffusivities_data[4][1] = d_data_diffusivities;
                        diffusivities_component_idx[4][1] = 9;
                        
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        /*
                         * Mass equation.
                         */
                        
                        diffusivities_data[0].resize(0);
                        diffusivities_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[1].resize(1);
                        diffusivities_component_idx[1].resize(1);
                        
                        // -mu.
                        diffusivities_data[1][0] = d_data_diffusivities;
                        diffusivities_component_idx[1][0] = 2;
                        
                        diffusivities_data[2].resize(1);
                        diffusivities_component_idx[2].resize(1);
                        
                        // -(4/3*mu + mu_v).
                        diffusivities_data[2][0] = d_data_diffusivities;
                        diffusivities_component_idx[2][0] = 0;
                        
                        diffusivities_data[3].resize(1);
                        diffusivities_component_idx[3].resize(1);
                        
                        // -mu.
                        diffusivities_data[3][0] = d_data_diffusivities;
                        diffusivities_component_idx[3][0] = 2;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[4].resize(4);
                        diffusivities_component_idx[4].resize(4);
                        
                        // -u*mu.
                        diffusivities_data[4][0] = d_data_diffusivities;
                        diffusivities_component_idx[4][0] = 9;
                        
                        // -v*(4/3*mu + mu_v).
                        diffusivities_data[4][1] = d_data_diffusivities;
                        diffusivities_component_idx[4][1] = 4;
                        
                        // -w*mu.
                        diffusivities_data[4][2] = d_data_diffusivities;
                        diffusivities_component_idx[4][2] = 11;
                        
                        // -kappa.
                        diffusivities_data[4][3] = d_data_diffusivities;
                        diffusivities_component_idx[4][3] = 12;
                        
                        break;
                    }
                    case DIRECTION::Z_DIRECTION:
                    {
                        /*
                         * Mass equation.
                         */
                        
                        diffusivities_data[0].resize(0);
                        diffusivities_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[1].resize(0);
                        diffusivities_component_idx[1].resize(0);
                        
                        diffusivities_data[2].resize(1);
                        diffusivities_component_idx[2].resize(1);
                        
                        // 2/3*(mu - mu_v).
                        diffusivities_data[2][0] = d_data_diffusivities;
                        diffusivities_component_idx[2][0] = 1;
                        
                        diffusivities_data[3].resize(1);
                        diffusivities_component_idx[3].resize(1);
                        
                        // -mu.
                        diffusivities_data[3][0] = d_data_diffusivities;
                        diffusivities_component_idx[3][0] = 2;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[4].resize(2);
                        diffusivities_component_idx[4].resize(2);
                        
                        // -w*u.
                        diffusivities_data[4][0] = d_data_diffusivities;
                        diffusivities_component_idx[4][0] = 11;
                        
                        // v*(2/3*mu - mu_v).
                        diffusivities_data[4][1] = d_data_diffusivities;
                        diffusivities_component_idx[4][1] = 7;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelSingleSpecies::getDiffusiveFluxDiffusivities()\n"
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
                         * Mass equation.
                         */
                        
                        diffusivities_data[0].resize(0);
                        diffusivities_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[1].resize(1);
                        diffusivities_component_idx[1].resize(1);
                        
                        // -mu.
                        diffusivities_data[1][0] = d_data_diffusivities;
                        diffusivities_component_idx[1][0] = 2;
                        
                        diffusivities_data[2].resize(0);
                        diffusivities_component_idx[2].resize(0);
                        
                        diffusivities_data[3].resize(1);
                        diffusivities_component_idx[3].resize(1);
                        
                        // 2/3*mu - mu_v.
                        diffusivities_data[3][0] = d_data_diffusivities;
                        diffusivities_component_idx[3][0] = 1;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[4].resize(2);
                        diffusivities_component_idx[4].resize(2);
                        
                        // w*(2/3*mu - mu_v).
                        diffusivities_data[4][0] = d_data_diffusivities;
                        diffusivities_component_idx[4][0] = 8;
                        
                        // -u*mu.
                        diffusivities_data[4][1] = d_data_diffusivities;
                        diffusivities_component_idx[4][1] = 9;
                        
                        break;
                    }
                    case DIRECTION::Y_DIRECTION:
                    {
                        /*
                         * Mass equation.
                         */
                        
                        diffusivities_data[0].resize(0);
                        diffusivities_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[1].resize(0);
                        diffusivities_component_idx[1].resize(0);
                        
                        diffusivities_data[2].resize(1);
                        diffusivities_component_idx[2].resize(1);
                        
                        // -mu.
                        diffusivities_data[2][0] = d_data_diffusivities;
                        diffusivities_component_idx[2][0] = 2;
                        
                        diffusivities_data[3].resize(1);
                        diffusivities_component_idx[3].resize(1);
                        
                        // 2/3*mu - mu_v.
                        diffusivities_data[3][0] = d_data_diffusivities;
                        diffusivities_component_idx[3][0] = 1;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[4].resize(2);
                        diffusivities_component_idx[4].resize(2);
                        
                        // w*(2/3*mu - mu_v).
                        diffusivities_data[4][0] = d_data_diffusivities;
                        diffusivities_component_idx[4][0] = 8;
                        
                        // -v*mu.
                        diffusivities_data[4][1] = d_data_diffusivities;
                        diffusivities_component_idx[4][1] = 10;
                        
                        break;
                    }
                    case DIRECTION::Z_DIRECTION:
                    {
                        /*
                         * Mass equation.
                         */
                        
                        diffusivities_data[0].resize(0);
                        diffusivities_component_idx[0].resize(0);
                        
                        /*
                         * Momentum equation.
                         */
                        
                        diffusivities_data[1].resize(1);
                        diffusivities_component_idx[1].resize(1);
                        
                        // -mu.
                        diffusivities_data[1][0] = d_data_diffusivities;
                        diffusivities_component_idx[1][0] = 2;
                        
                        diffusivities_data[2].resize(1);
                        diffusivities_component_idx[2].resize(1);
                        
                        // -mu.
                        diffusivities_data[2][0] = d_data_diffusivities;
                        diffusivities_component_idx[2][0] = 2;
                        
                        diffusivities_data[3].resize(1);
                        diffusivities_component_idx[3].resize(1);
                        
                        // -(4/3*mu + mu_v).
                        diffusivities_data[3][0] = d_data_diffusivities;
                        diffusivities_component_idx[3][0] = 0;
                        
                        /*
                         * Energy equation.
                         */
                        
                        diffusivities_data[4].resize(4);
                        diffusivities_component_idx[4].resize(4);
                        
                        // -u*mu.
                        diffusivities_data[4][0] = d_data_diffusivities;
                        diffusivities_component_idx[4][0] = 9;
                        
                        // -v*mu.
                        diffusivities_data[4][1] = d_data_diffusivities;
                        diffusivities_component_idx[4][1] = 10;
                        
                        // -w*(4/3*mu + mu_v).
                        diffusivities_data[4][2] = d_data_diffusivities;
                        diffusivities_component_idx[4][2] = 5;
                        
                        // -kappa.
                        diffusivities_data[4][3] = d_data_diffusivities;
                        diffusivities_component_idx[4][3] = 12;
                        
                        break;
                    }
                    default:
                    {
                        TBOX_ERROR(d_object_name
                            << ": FlowModelSingleSpecies::getDiffusiveFluxDiffusivities()\n"
                            << "There are only x-direction, y-direction and z-direction for three-dimensional problem."
                            << std::endl);
                    }
                }
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::getDiffusiveFluxDiffusivities()\n"
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
    
    if (d_num_subghosts_diffusivities > -hier::IntVector::getOne(d_dim))
    {
        d_subghost_box_diffusivities = d_interior_box;
        d_subghost_box_diffusivities.grow(d_num_subghosts_diffusivities);
        d_subghostcell_dims_diffusivities = d_subghost_box_diffusivities.numberCells();
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
        // Create the cell data of velocity.
        d_data_velocity.reset(
            new pdat::CellData<double>(d_interior_box, d_dim.getValue(), d_num_subghosts_velocity));
        
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
        // Create the cell data of internal energy.
        d_data_internal_energy.reset(
            new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_internal_energy));
        
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
        
        if (!d_data_velocity)
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
        // Create the cell data of pressure.
        d_data_pressure.reset(
            new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_pressure));
        
        // Get the cell data of the variables density, momentum and total energy.
        boost::shared_ptr<pdat::CellData<double> > data_density =
            getCellDataOfDensity();
        
        if (!d_data_internal_energy)
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
        // Create the cell data of sound speed.
        d_data_sound_speed.reset(
            new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_sound_speed));
        
        // Get the cell data of the variable density and pressure.
        boost::shared_ptr<pdat::CellData<double> > data_density =
            getCellDataOfDensity();
        
        if (!d_data_pressure)
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
        // Create the cell data of temperature.
        d_data_temperature.reset(
            new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_temperature));
        
        // Get the cell data of the variable density and pressure.
        boost::shared_ptr<pdat::CellData<double> > data_density =
            getCellDataOfDensity();
        
        if (!d_data_pressure)
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
            // Create the cell data of convective flux in the x-direction.
            d_data_convective_flux_x.reset(
                new pdat::CellData<double>(d_interior_box, d_num_eqn, d_num_subghosts_convective_flux_x));
            
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
            
            if (!d_data_velocity)
            {
                computeCellDataOfVelocity(domain);
            }
            
            if (!d_data_pressure)
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
            // Create the cell data of convective flux in the y-direction.
            d_data_convective_flux_y.reset(
                new pdat::CellData<double>(d_interior_box, d_num_eqn, d_num_subghosts_convective_flux_y));
            
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
            
            if (!d_data_velocity)
            {
                computeCellDataOfVelocity(domain);
            }
            
            if (!d_data_pressure)
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
            // Create the cell data of convective flux in the z-direction.
            d_data_convective_flux_z.reset(
                new pdat::CellData<double>(d_interior_box, d_num_eqn, d_num_subghosts_convective_flux_z));
            
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
            
            if (!d_data_velocity)
            {
                computeCellDataOfVelocity(domain);
            }
            
            if (!d_data_pressure)
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
            // Create the cell data of maximum wave speed in the x-direction.
            d_data_max_wave_speed_x.reset(
                new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_max_wave_speed_x));
            
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
            
            if (!d_data_velocity)
            {
                computeCellDataOfVelocity(domain);
            }
            
            if (!d_data_sound_speed)
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
            // Create the cell data of maximum wave speed in the y-direction.
            d_data_max_wave_speed_y.reset(
                new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_max_wave_speed_y));
            
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
            
            if (!d_data_velocity)
            {
                computeCellDataOfVelocity(domain);
            }
            
            if (!d_data_sound_speed)
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
            // Create the cell data of maximum wave speed in the z-direction.
            d_data_max_wave_speed_z.reset(
                new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_max_wave_speed_z));
            
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
            
            if (!d_data_velocity)
            {
                computeCellDataOfVelocity(domain);
            }
            
            if (!d_data_sound_speed)
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
        // Create the cell data of maximum diffusivity.
        d_data_max_diffusivity.reset(
            new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_max_diffusivity));
        
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
        
        if (!d_data_pressure)
        {
            computeCellDataOfPressureWithInternalEnergy(domain);
        }
        
        if (!d_data_temperature)
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
