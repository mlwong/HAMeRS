#include "flow/flow_models/single-species/FlowModelSingleSpecies.hpp"

#include "flow/flow_models/single-species/FlowModelBoundaryUtilitiesSingleSpecies.hpp"
#include "flow/flow_models/single-species/FlowModelStatisticsUtilitiesSingleSpecies.hpp"

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
        d_num_subghosts_dilatation(-hier::IntVector::getOne(d_dim)),
        d_num_subghosts_vorticity(-hier::IntVector::getOne(d_dim)),
        d_num_subghosts_enstrophy(-hier::IntVector::getOne(d_dim)),
        d_num_subghosts_convective_flux_x(-hier::IntVector::getOne(d_dim)),
        d_num_subghosts_convective_flux_y(-hier::IntVector::getOne(d_dim)),
        d_num_subghosts_convective_flux_z(-hier::IntVector::getOne(d_dim)),
        d_num_subghosts_max_wave_speed_x(-hier::IntVector::getOne(d_dim)),
        d_num_subghosts_max_wave_speed_y(-hier::IntVector::getOne(d_dim)),
        d_num_subghosts_max_wave_speed_z(-hier::IntVector::getOne(d_dim)),
        d_num_subghosts_diffusivities(-hier::IntVector::getOne(d_dim)),
        d_subghost_box_velocity(hier::Box::getEmptyBox(dim)),
        d_subghost_box_internal_energy(hier::Box::getEmptyBox(dim)),
        d_subghost_box_pressure(hier::Box::getEmptyBox(dim)),
        d_subghost_box_sound_speed(hier::Box::getEmptyBox(dim)),
        d_subghost_box_temperature(hier::Box::getEmptyBox(dim)),
        d_subghost_box_dilatation(hier::Box::getEmptyBox(dim)),
        d_subghost_box_vorticity(hier::Box::getEmptyBox(dim)),
        d_subghost_box_enstrophy(hier::Box::getEmptyBox(dim)),
        d_subghost_box_convective_flux_x(hier::Box::getEmptyBox(dim)),
        d_subghost_box_convective_flux_y(hier::Box::getEmptyBox(dim)),
        d_subghost_box_convective_flux_z(hier::Box::getEmptyBox(dim)),
        d_subghost_box_max_wave_speed_x(hier::Box::getEmptyBox(dim)),
        d_subghost_box_max_wave_speed_y(hier::Box::getEmptyBox(dim)),
        d_subghost_box_max_wave_speed_z(hier::Box::getEmptyBox(dim)),
        d_subghost_box_diffusivities(hier::Box::getEmptyBox(dim)),
        d_subghostcell_dims_velocity(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_internal_energy(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_pressure(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_sound_speed(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_temperature(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_dilatation(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_vorticity(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_enstrophy(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_convective_flux_x(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_convective_flux_y(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_convective_flux_z(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_max_wave_speed_x(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_max_wave_speed_y(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_max_wave_speed_z(hier::IntVector::getZero(d_dim)),
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
    
    d_variable_density = boost::shared_ptr<pdat::CellVariable<double> > (
        new pdat::CellVariable<double>(d_dim, "density", 1));
    
    d_variable_momentum = boost::shared_ptr<pdat::CellVariable<double> > (
        new pdat::CellVariable<double>(d_dim, "momentum", d_dim.getValue()));
    
    d_variable_total_energy = boost::shared_ptr<pdat::CellVariable<double> > (
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
     * Initialize statistics utilities object.
     */
    d_flow_model_statistics_utilities.reset(new FlowModelStatisticsUtilitiesSingleSpecies(
        "d_flow_model_statistics_utilities",
        d_dim,
        d_grid_geometry,
        d_num_species,
        flow_model_db));
    
    /*
     * Initialize the Riemann solvers.
     */
    d_Riemann_solver_HLLC = boost::shared_ptr<RiemannSolverSingleSpeciesHLLC> (
        new RiemannSolverSingleSpeciesHLLC(
            d_object_name,
            d_dim,
            d_num_eqn,
            d_num_species,
            d_equation_of_state_mixing_rules));
    
    d_Riemann_solver_HLLC_HLL = boost::shared_ptr<RiemannSolverSingleSpeciesHLLC_HLL> (
        new RiemannSolverSingleSpeciesHLLC_HLL(
            d_object_name,
            d_dim,
            d_num_eqn,
            d_num_species,
            d_equation_of_state_mixing_rules));
    
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
        d_variable_density,
        num_ghosts,
        num_ghosts_intermediate,
        RungeKuttaLevelIntegrator::TIME_DEP,
        d_grid_geometry,
        "CONSERVATIVE_COARSEN",
        "CONSERVATIVE_LINEAR_REFINE");
    
    integrator->registerVariable(
        d_variable_momentum,
        num_ghosts,
        num_ghosts_intermediate,
        RungeKuttaLevelIntegrator::TIME_DEP,
        d_grid_geometry,
        "CONSERVATIVE_COARSEN",
        "CONSERVATIVE_LINEAR_REFINE");
    
    integrator->registerVariable(
        d_variable_total_energy,
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
    
    conservative_variables.push_back(d_variable_density);
    conservative_variables.push_back(d_variable_momentum);
    conservative_variables.push_back(d_variable_total_energy);
    
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
            d_patch->getPatchData(d_variable_density, getDataContext())));
    
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
FlowModelSingleSpecies::registerDerivedCellVariable(
    const std::unordered_map<std::string, hier::IntVector>& num_subghosts_of_data)

{
    // Check whether a patch is already registered.
    if (!d_patch)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::registerDerivedCellVariable()\n"
            << "No patch is registered yet."
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
                << ": FlowModelSingleSpecies::registerDerivedCellVariable()\n"
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
}


/*
 * Register the required derived variables for transformation between conservative
 * variables and characteristic variables.
 */
void
FlowModelSingleSpecies::registerDerivedVariablesForCharacteristicProjectionOfConservativeVariables(
    const hier::IntVector& num_subghosts,
    const AVERAGING::TYPE& averaging)
{
    // Check whether a patch is already registered.
    if (!d_patch)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "registerDerivedVariablesForCharacteristicProjectionOfConservativeVariables()\n"
            << "No patch is registered yet."
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
FlowModelSingleSpecies::registerDerivedVariablesForCharacteristicProjectionOfPrimitiveVariables(
    const hier::IntVector& num_subghosts,
    const AVERAGING::TYPE& averaging)
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
FlowModelSingleSpecies::registerDiffusiveFlux(const hier::IntVector& num_subghosts)
{
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
    
    d_diffusive_flux_var_registered = true;
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
    d_num_subghosts_dilatation        = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_vorticity         = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_enstrophy         = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_convective_flux_x = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_convective_flux_y = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_convective_flux_z = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_max_wave_speed_x  = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_max_wave_speed_y  = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_max_wave_speed_z  = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_diffusivities     = -hier::IntVector::getOne(d_dim);
    
    d_interior_box                   = hier::Box::getEmptyBox(d_dim);
    d_ghost_box                      = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_velocity          = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_internal_energy   = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_pressure          = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_sound_speed       = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_temperature       = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_dilatation        = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_vorticity         = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_enstrophy         = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_convective_flux_x = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_convective_flux_y = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_convective_flux_z = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_max_wave_speed_x  = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_max_wave_speed_y  = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_max_wave_speed_z  = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_diffusivities     = hier::Box::getEmptyBox(d_dim);
    
    d_interior_dims                       = hier::IntVector::getZero(d_dim);
    d_ghostcell_dims                      = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_velocity          = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_internal_energy   = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_pressure          = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_sound_speed       = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_temperature       = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_dilatation        = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_vorticity         = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_enstrophy         = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_convective_flux_x = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_convective_flux_y = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_convective_flux_z = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_max_wave_speed_x  = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_max_wave_speed_y  = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_max_wave_speed_z  = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_diffusivities     = hier::IntVector::getZero(d_dim);
    
    d_data_velocity.reset();
    d_data_internal_energy.reset();
    d_data_pressure.reset();
    d_data_sound_speed.reset();
    d_data_temperature.reset();
    d_data_dilatation.reset();
    d_data_vorticity.reset();
    d_data_enstrophy.reset();
    d_data_convective_flux_x.reset();
    d_data_convective_flux_y.reset();
    d_data_convective_flux_z.reset();
    d_data_max_wave_speed_x.reset();
    d_data_max_wave_speed_y.reset();
    d_data_max_wave_speed_z.reset();
    d_data_diffusivities.reset();
    
    clearDataContext();
}


/*
 * Compute global cell data of different registered derived variables with the registered data context.
 */
void
FlowModelSingleSpecies::computeGlobalDerivedCellData()
{
    // Check whether a patch is already registered.
    if (!d_patch)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::computeGlobalDerivedCellData()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    /*
     * Set the boxes and their dimensions for the derived cell variables.
     */
    setGhostBoxesAndDimensionsDerivedCellVariables();
    
    // Compute the velocity cell data.
    if (d_num_subghosts_velocity > -hier::IntVector::getOne(d_dim))
    {
        if (!d_data_velocity)
        {
            computeGlobalCellDataVelocity();
        }
    }
    
    // Compute the internal energy cell data.
    if (d_num_subghosts_internal_energy > -hier::IntVector::getOne(d_dim))
    {
        if (!d_data_internal_energy)
        {
            computeGlobalCellDataInternalEnergyWithVelocity();
        }
    }
    
    // Compute the pressure cell data.
    if (d_num_subghosts_pressure > -hier::IntVector::getOne(d_dim))
    {
        if (!d_data_pressure)
        {
            computeGlobalCellDataPressureWithInternalEnergy();
        }
    }
    
    // Compute the sound speed cell data.
    if (d_num_subghosts_sound_speed > -hier::IntVector::getOne(d_dim))
    {
        if (!d_data_sound_speed)
        {
            computeGlobalCellDataSoundSpeedWithPressure();
        }
    }
    
    // Compute the temperature cell data.
    if (d_num_subghosts_temperature > -hier::IntVector::getOne(d_dim))
    {
        if (!d_data_temperature)
        {
            computeGlobalCellDataTemperatureWithPressure();
        }
    }
    
    // Compute the dilatation cell data.
    if (d_num_subghosts_dilatation > -hier::IntVector::getOne(d_dim))
    {
        if (!d_data_dilatation)
        {
            computeGlobalCellDataDilatationWithVelocity();
        }
    }
    
    // Compute the vorticity cell data.
    if (d_num_subghosts_vorticity > -hier::IntVector::getOne(d_dim))
    {
        if (!d_data_vorticity)
        {
            computeGlobalCellDataVorticityWithVelocity();
        }
    }
    
    // Compute the enstrophy cell data.
    if (d_num_subghosts_enstrophy > -hier::IntVector::getOne(d_dim))
    {
        if (!d_data_enstrophy)
        {
            computeGlobalCellDataEnstrophyWithVorticity();
        }
    }
    
    // Compute the x-direction convective flux cell data.
    if (d_num_subghosts_convective_flux_x > -hier::IntVector::getOne(d_dim))
    {
        if (!d_data_convective_flux_x)
        {
            computeGlobalCellDataConvectiveFluxWithVelocityAndPressure(DIRECTION::X_DIRECTION);
        }
    }
    
    // Compute the y-direction convective flux cell data.
    if (d_num_subghosts_convective_flux_y > -hier::IntVector::getOne(d_dim))
    {
        if (!d_data_convective_flux_y)
        {
            computeGlobalCellDataConvectiveFluxWithVelocityAndPressure(DIRECTION::Y_DIRECTION);
        }
    }
    
    // Compute the z-direction convective flux cell data.
    if (d_num_subghosts_convective_flux_z > -hier::IntVector::getOne(d_dim))
    {
        if (!d_data_convective_flux_z)
        {
            computeGlobalCellDataConvectiveFluxWithVelocityAndPressure(DIRECTION::Z_DIRECTION);
        }
    }
    
    // Compute the x-direction maximum wave speed cell data.
    if (d_num_subghosts_max_wave_speed_x > -hier::IntVector::getOne(d_dim))
    {
        if (!d_data_max_wave_speed_x)
        {
            computeGlobalCellDataMaxWaveSpeedWithVelocityAndSoundSpeed(DIRECTION::X_DIRECTION);
        }
    }
    
    // Compute the y-direction maximum wave speed cell data.
    if (d_num_subghosts_max_wave_speed_y > -hier::IntVector::getOne(d_dim))
    {
        if (!d_data_max_wave_speed_y)
        {
            computeGlobalCellDataMaxWaveSpeedWithVelocityAndSoundSpeed(DIRECTION::Y_DIRECTION);
        }
    }
    
    // Compute the z-direction maximum wave speed cell data.
    if (d_num_subghosts_max_wave_speed_z > -hier::IntVector::getOne(d_dim))
    {
        if (!d_data_max_wave_speed_z)
        {
            computeGlobalCellDataMaxWaveSpeedWithVelocityAndSoundSpeed(DIRECTION::Z_DIRECTION);
        }
    }
}


/*
 * Get the global cell data of one cell variable in the registered patch.
 */
boost::shared_ptr<pdat::CellData<double> >
FlowModelSingleSpecies::getGlobalCellData(const std::string& variable_key)
{
    // Check whether a patch is already registered.
    if (!d_patch)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::getGlobalCellData()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    boost::shared_ptr<pdat::CellData<double> > cell_data;
    
    if (variable_key == "DENSITY")
    {
        cell_data = getGlobalCellDataDensity();
    }
    else if (variable_key == "MOMENTUM")
    {
        cell_data = getGlobalCellDataMomentum();
    }
    else if (variable_key == "TOTAL_ENERGY")
    {
        cell_data = getGlobalCellDataTotalEnergy();
    }
    else if (variable_key == "VELOCITY")
    {
        if (!d_data_velocity)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::getGlobalCellData()\n"
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
                << ": FlowModelSingleSpecies::getGlobalCellData()\n"
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
                << ": FlowModelSingleSpecies::getGlobalCellData()\n"
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
                << ": FlowModelSingleSpecies::getGlobalCellData()\n"
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
                << ": FlowModelSingleSpecies::getGlobalCellData()\n"
                << "Cell data of 'TEMPERATURE' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_temperature;
    }
    else if (variable_key == "DILATATION")
    {
        if (!d_data_dilatation)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::getGlobalCellData()\n"
                << "Cell data of 'DILATATION' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_dilatation;
    }
    else if (variable_key == "VORTICITY")
    {
        if (!d_data_vorticity)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::getGlobalCellData()\n"
                << "Cell data of 'VORTICITY' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_vorticity;
    }
    else if (variable_key == "ENSTROPHY")
    {
        if (!d_data_enstrophy)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::getGlobalCellData()\n"
                << "Cell data of 'ENSTROPHY' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_enstrophy;
    }
    else if (variable_key == "CONVECTIVE_FLUX_X")
    {
        if (!d_data_convective_flux_x)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::getGlobalCellData()\n"
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
                << ": FlowModelSingleSpecies::getGlobalCellData()\n"
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
                << ": FlowModelSingleSpecies::getGlobalCellData()\n"
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
                << ": FlowModelSingleSpecies::getGlobalCellData()\n"
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
                << ": FlowModelSingleSpecies::getGlobalCellData()\n"
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
                << ": FlowModelSingleSpecies::getGlobalCellData()\n"
                << "Cell data of 'MAX_WAVE_SPEED_Z' is not registered/computed yet."
                << std::endl);
        }
        cell_data = d_data_max_wave_speed_z;
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::getGlobalCellData()\n"
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
FlowModelSingleSpecies::getGlobalCellData(
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
FlowModelSingleSpecies::fillZeroGlobalCellDataConservativeVariables()
{
    // Check whether a patch is already registered.
    if (!d_patch)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::fillZeroGlobalCellDataConservativeVariables()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    boost::shared_ptr<pdat::CellData<double> > data_density = getGlobalCellDataDensity();
    boost::shared_ptr<pdat::CellData<double> > data_momentum = getGlobalCellDataMomentum();
    boost::shared_ptr<pdat::CellData<double> > data_total_energy = getGlobalCellDataTotalEnergy();
    
    data_density->fillAll(0.0, d_interior_box);
    data_momentum->fillAll(0.0, d_interior_box);
    data_total_energy->fillAll(0.0, d_interior_box);
}


/*
 * Update the interior global cell data of conservative variables.
 */
void
FlowModelSingleSpecies::updateGlobalCellDataConservativeVariables()
{
    // Check whether a patch is already registered.
    if (!d_patch)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::updateGlobalCellDataConservativeVariables()\n"
            << "No patch is registered yet."
            << std::endl);
    }
}


/*
 * Get the global cell data of the conservative variables in the registered patch.
 */
std::vector<boost::shared_ptr<pdat::CellData<double> > >
FlowModelSingleSpecies::getGlobalCellDataConservativeVariables()
{
    // Check whether a patch is already registered.
    if (!d_patch)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::getGlobalCellDataConservativeVariables()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    std::vector<boost::shared_ptr<pdat::CellData<double> > > global_cell_data;
    global_cell_data.reserve(3);
    
    global_cell_data.push_back(getGlobalCellDataDensity());
    global_cell_data.push_back(getGlobalCellDataMomentum());
    global_cell_data.push_back(getGlobalCellDataTotalEnergy());
    
    return global_cell_data;
}


/*
 * Get the global cell data of the primitive variables in the registered patch.
 */
std::vector<boost::shared_ptr<pdat::CellData<double> > >
FlowModelSingleSpecies::getGlobalCellDataPrimitiveVariables()
{
    // Check whether a patch is already registered.
    if (!d_patch)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::getGlobalCellDataPrimitiveVariables()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    std::vector<boost::shared_ptr<pdat::CellData<double> > > global_cell_data;
    global_cell_data.reserve(3);
    
    global_cell_data.push_back(getGlobalCellDataDensity());
    if (!d_data_velocity)
    {
        computeGlobalCellDataVelocity();
    }
    global_cell_data.push_back(d_data_velocity);
    if (!d_data_pressure)
    {
        computeGlobalCellDataPressureWithInternalEnergy();
    }
    global_cell_data.push_back(d_data_pressure);
    
    return global_cell_data;
}


/*
 * Get the number of projection variables for transformation between conservative
 * variables and characteristic variables.
 */
int
FlowModelSingleSpecies::getNumberOfProjectionVariablesForConservativeVariables() const
{
    TBOX_ERROR(d_object_name
        << ": FlowModelSingleSpecies::getNumberOfProjectionVariablesForConservativeVariables()\n"
        << "Method getNumberOfProjectionVariablesForConservativeVariables() is not yet implemented."
        << std::endl);
    
    return 0;
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
 * Compute global side data of the projection variables for transformation between
 * conservative variables and characteristic variables.
 */
void
FlowModelSingleSpecies::computeGlobalSideDataProjectionVariablesForConservativeVariables(
    std::vector<boost::shared_ptr<pdat::SideData<double> > >& projection_variables)
{
    NULL_USE(projection_variables);
    
    TBOX_ERROR(d_object_name
        << ": FlowModelSingleSpecies::"
        << "computeGlobalSideDataProjectionVariablesForConservativeVariables()\n"
        << "Method computeGlobalSideDataProjectionVariablesForConservativeVariables() is not"
        << " yet implemented."
        << std::endl);
}


/*
 * Compute global side data of the projection variables for transformation between
 * primitive variables and characteristic variables.
 */
void
FlowModelSingleSpecies::computeGlobalSideDataProjectionVariablesForPrimitiveVariables(
    std::vector<boost::shared_ptr<pdat::SideData<double> > >& projection_variables)
{
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
            << "computeGlobalSideDataProjectionVariablesForPrimitiveVariables()\n"
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
                << "computeGlobalSideDataProjectionVariablesForPrimitiveVariables()\n"
                << "The interior dimension of the projection variables does not match that of patch."
                << std::endl);
        }
    }
    
    if (num_ghosts_projection_var != projection_variables[1]->getGhostCellWidth())
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "computeGlobalSideDataProjectionVariablesForPrimitiveVariables()\n"
            << "The projection variables don't have same ghost cell width."
            << std::endl);
    }
    
    if (num_ghosts_projection_var > d_num_ghosts)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "computeGlobalSideDataProjectionVariablesForPrimitiveVariables()\n"
            << "The projection variables have ghost cell width larger than that of density."
            << std::endl);
    }
    
    if (num_ghosts_projection_var > d_num_subghosts_sound_speed)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "computeGlobalSideDataProjectionVariablesForPrimitiveVariables()\n"
            << "The projection variables have ghost cell width larger than that of sound speed."
            << std::endl);
    }
    
    // Get the cell data of the variable density.
    boost::shared_ptr<pdat::CellData<double> > data_density =
        getGlobalCellDataDensity();
    
    // Get the pointers to the cell data of density and sound speed.
    double* rho = data_density->getPointer(0);
    if (!d_data_sound_speed)
    {
        computeGlobalCellDataSoundSpeedWithPressure();
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
        
        switch (d_proj_var_primitive_averaging)
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
                    
                    rho_average[idx_face_x] = 0.5*(rho[idx_L] + rho[idx_R]);
                    c_average[idx_face_x] = 0.5*(c[idx_sound_speed_L] + c[idx_sound_speed_R]);
                }
                
                break;
            }
            case AVERAGING::ROE:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::"
                    << "computeGlobalSideDataProjectionVariablesForPrimitiveVariables()\n"
                    << "Roe averaging is not yet implemented."
                    << std::endl);
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::"
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
                        
                        rho_average[idx_face_x] = 0.5*(rho[idx_L] + rho[idx_R]);
                        c_average[idx_face_x] = 0.5*(c[idx_sound_speed_L] + c[idx_sound_speed_R]);
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
                        
                        rho_average[idx_face_y] = 0.5*(rho[idx_B] + rho[idx_T]);
                        c_average[idx_face_y] = 0.5*(c[idx_sound_speed_B] + c[idx_sound_speed_T]);
                    }
                }
                
                break;
            }
            case AVERAGING::ROE:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::"
                    << "computeGlobalSideDataProjectionVariablesForPrimitiveVariables()\n"
                    << "Roe averaging is not yet implemented."
                    << std::endl);
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::"
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
                            
                            rho_average[idx_face_x] = 0.5*(rho[idx_L] + rho[idx_R]);
                            c_average[idx_face_x] = 0.5*(c[idx_sound_speed_L] + c[idx_sound_speed_R]);
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
                            
                            rho_average[idx_face_y] = 0.5*(rho[idx_B] + rho[idx_T]);
                            c_average[idx_face_y] = 0.5*(c[idx_sound_speed_B] + c[idx_sound_speed_T]);
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
                            
                            rho_average[idx_face_z] = 0.5*(rho[idx_B] + rho[idx_F]);
                            c_average[idx_face_z] = 0.5*(c[idx_sound_speed_B] + c[idx_sound_speed_F]);
                        }
                    }
                }
                
                break;
            }
            case AVERAGING::ROE:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::"
                    << "computeGlobalSideDataProjectionVariablesForPrimitiveVariables()\n"
                    << "Roe averaging is not yet implemented."
                    << std::endl);
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::"
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
FlowModelSingleSpecies::computeGlobalSideDataCharacteristicVariablesFromConservativeVariables(
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
        << ": FlowModelSingleSpecies::"
        << "computeGlobalSideDataCharacteristicVariablesFromConservativeVariables()\n"
        << "Method computeGlobalSideDataCharacteristicVariablesFromConservativeVariables()"
        << " is not yet implemented."
        << std::endl);
}


/*
 * Compute global side data of characteristic variables from primitive variables.
 */
void
FlowModelSingleSpecies::computeGlobalSideDataCharacteristicVariablesFromPrimitiveVariables(
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
            << "computeGlobalSideDataCharacteristicVariablesFromPrimitiveVariables()\n"
            << "The number of characteristic variables are incorrect."
            << std::endl);
    }
    if (static_cast<int>(primitive_variables.size()) != 3)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "computeGlobalSideDataCharacteristicVariablesFromPrimitiveVariables()\n"
            << "The number of primitive variables are incorrect."
            << std::endl);
    }
    if (primitive_variables[0]->getDepth() != 1 ||
        primitive_variables[1]->getDepth() != d_dim.getValue() ||
        primitive_variables[2]->getDepth() != 1)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "computeGlobalSideDataCharacteristicVariablesFromPrimitiveVariables()\n"
            << "The depths of one or more primitive variables are incorrect."
            << std::endl);
    }
    if (static_cast<int>(projection_variables.size()) != 2)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "computeGlobalSideDataCharacteristicVariablesFromPrimitiveVariables()\n"
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
                << ": FlowModelSingleSpecies::"
                << "computeGlobalSideDataCharacteristicVariablesFromPrimitiveVariables()\n"
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
                << ": FlowModelSingleSpecies::"
                << "computeGlobalSideDataCharacteristicVariablesFromPrimitiveVariables()\n"
                << "The characteristic variables don't have same ghost cell width."
                << std::endl);
        }
    }
    if (num_ghosts_projection_var != projection_variables[1]->getGhostCellWidth())
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "computeGlobalSideDataCharacteristicVariablesFromPrimitiveVariables()\n"
            << "The projection variables don't have same ghost cell width."
            << std::endl);
    }
    
    if (num_ghosts_projection_var != num_ghosts_characteristic_var)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
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
                << ": FlowModelSingleSpecies::"
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
            
            W[0][idx_face] = -0.5*rho_average[idx_face]*c_average[idx_face]*V[1][idx_vel] + 0.5*V[2][idx_p];
            W[1][idx_face] = V[0][idx_rho] - 1.0/(c_average[idx_face]*c_average[idx_face])*V[2][idx_p];
            W[2][idx_face] = 0.5*rho_average[idx_face]*c_average[idx_face]*V[1][idx_vel] + 0.5*V[2][idx_p];
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
                
                W[0][idx_face] = -0.5*rho_average[idx_face]*c_average[idx_face]*V[1][idx_vel] + 0.5*V[3][idx_p];
                W[1][idx_face] = V[0][idx_rho] - 1.0/(c_average[idx_face]*c_average[idx_face])*V[3][idx_p];
                W[2][idx_face] = V[2][idx_vel];
                W[3][idx_face] = 0.5*rho_average[idx_face]*c_average[idx_face]*V[1][idx_vel] + 0.5*V[3][idx_p];
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
                
                W[0][idx_face] = -0.5*rho_average[idx_face]*c_average[idx_face]*V[2][idx_vel] + 0.5*V[3][idx_p];
                W[1][idx_face] = V[0][idx_rho] - 1.0/(c_average[idx_face]*c_average[idx_face])*V[3][idx_p];
                W[2][idx_face] = V[1][idx_vel];
                W[3][idx_face] = 0.5*rho_average[idx_face]*c_average[idx_face]*V[2][idx_vel] + 0.5*V[3][idx_p];
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
                    
                    W[0][idx_face] = -0.5*rho_average[idx_face]*c_average[idx_face]*V[1][idx_vel] + 0.5*V[4][idx_p];
                    W[1][idx_face] = V[0][idx_rho] - 1.0/(c_average[idx_face]*c_average[idx_face])*V[4][idx_p];
                    W[2][idx_face] = V[2][idx_vel];
                    W[3][idx_face] = V[3][idx_vel];
                    W[4][idx_face] = 0.5*rho_average[idx_face]*c_average[idx_face]*V[1][idx_vel] + 0.5*V[4][idx_p];
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
                    
                    W[0][idx_face] = -0.5*rho_average[idx_face]*c_average[idx_face]*V[2][idx_vel] + 0.5*V[4][idx_p];
                    W[1][idx_face] = V[0][idx_rho] - 1.0/(c_average[idx_face]*c_average[idx_face])*V[4][idx_p];
                    W[2][idx_face] = V[1][idx_vel];
                    W[3][idx_face] = V[3][idx_vel];
                    W[4][idx_face] = 0.5*rho_average[idx_face]*c_average[idx_face]*V[2][idx_vel] + 0.5*V[4][idx_p];
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
                    
                    W[0][idx_face] = -0.5*rho_average[idx_face]*c_average[idx_face]*V[3][idx_vel] + 0.5*V[4][idx_p];
                    W[1][idx_face] = V[0][idx_rho] - 1.0/(c_average[idx_face]*c_average[idx_face])*V[4][idx_p];
                    W[2][idx_face] = V[1][idx_vel];
                    W[3][idx_face] = V[2][idx_vel];
                    W[4][idx_face] = 0.5*rho_average[idx_face]*c_average[idx_face]*V[3][idx_vel] + 0.5*V[4][idx_p];
                }
            }
        }
    }
}


/*
 * Compute global side data of conservative variables from characteristic variables.
 */
void
FlowModelSingleSpecies::computeGlobalSideDataConservativeVariablesFromCharacteristicVariables(
    std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& characteristic_variables,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& projection_variables)
{
    NULL_USE(conservative_variables);
    NULL_USE(characteristic_variables);
    NULL_USE(projection_variables);
    
    TBOX_ERROR(d_object_name
        << ": FlowModelSingleSpecies::"
        << "computeGlobalSideDataConservativeVariablesFromCharacteristicVariables()\n"
        << "Method computeGlobalSideDataConservativeVariablesFromCharacteristicVariables() is not"
        << " yet implemented."
        << std::endl);
}


/*
 * Compute global side data of primitive variables from characteristic variables.
 */
void
FlowModelSingleSpecies::computeGlobalSideDataPrimitiveVariablesFromCharacteristicVariables(
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
            << "computeGlobalSideDataPrimitiveVariablesFromCharacteristicVariables()\n"
            << "The number of characteristic variables are incorrect."
            << std::endl);
    }
    if (static_cast<int>(characteristic_variables.size()) != d_num_eqn)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "computeGlobalSideDataPrimitiveVariablesFromCharacteristicVariables()\n"
            << "The number of primitive variables are incorrect."
            << std::endl);
    }
    if (static_cast<int>(projection_variables.size()) != 2)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "computeGlobalSideDataPrimitiveVariablesFromCharacteristicVariables()\n"
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
                << ": FlowModelSingleSpecies::"
                << "computeGlobalSideDataPrimitiveVariablesFromCharacteristicVariables()\n"
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
                << ": FlowModelSingleSpecies::"
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
                << ": FlowModelSingleSpecies::"
                << "computeGlobalSideDataPrimitiveVariablesFromCharacteristicVariables()\n"
                << "The characteristic variables don't have same ghost cell width."
                << std::endl);
        }
    }
    if (num_ghosts_projection_var != projection_variables[1]->getGhostCellWidth())
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "computeGlobalSideDataPrimitiveVariablesFromCharacteristicVariables()\n"
            << "The projection variables don't have same ghost cell width."
            << std::endl);
    }
    
    if (num_ghosts_projection_var != num_ghosts_primitive_var)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "computeGlobalSideDataPrimitiveVariablesFromCharacteristicVariables()\n"
            << "The ghost cell width of the projection variables does not match that of primitive"
            << " variables."
            << std::endl);
    }
    if (num_ghosts_projection_var != num_ghosts_characteristic_var)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "computeGlobalSideDataPrimitiveVariablesFromCharacteristicVariables()\n"
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
            
            V[0][idx_face] = 1.0/(c_average[idx_face]*c_average[idx_face])*W[0][idx_face] +
                W[1][idx_face] + 1.0/(c_average[idx_face]*c_average[idx_face])*W[2][idx_face];
            V[1][idx_face] = -1.0/(rho_average[idx_face]*c_average[idx_face])*W[0][idx_face] +
                1.0/(rho_average[idx_face]*c_average[idx_face])*W[2][idx_face];
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
                
                V[0][idx_face] = 1.0/(c_average[idx_face]*c_average[idx_face])*W[0][idx_face] +
                    W[1][idx_face] + 1.0/(c_average[idx_face]*c_average[idx_face])*W[3][idx_face];
                V[1][idx_face] = -1.0/(rho_average[idx_face]*c_average[idx_face])*W[0][idx_face] +
                    1.0/(rho_average[idx_face]*c_average[idx_face])*W[3][idx_face];
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
                
                V[0][idx_face] = 1.0/(c_average[idx_face]*c_average[idx_face])*W[0][idx_face] + W[1][idx_face] +
                    1.0/(c_average[idx_face]*c_average[idx_face])*W[3][idx_face];
                V[1][idx_face] = W[2][idx_face];
                V[2][idx_face] = -1.0/(rho_average[idx_face]*c_average[idx_face])*W[0][idx_face] +
                    1.0/(rho_average[idx_face]*c_average[idx_face])*W[3][idx_face];
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
                    
                    V[0][idx_face] = 1.0/(c_average[idx_face]*c_average[idx_face])*W[0][idx_face] +
                        W[1][idx_face] + 1.0/(c_average[idx_face]*c_average[idx_face])*W[4][idx_face];
                    V[1][idx_face] = -1.0/(rho_average[idx_face]*c_average[idx_face])*W[0][idx_face] +
                        1.0/(rho_average[idx_face]*c_average[idx_face])*W[4][idx_face];
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
                    
                    V[0][idx_face] = 1.0/(c_average[idx_face]*c_average[idx_face])*W[0][idx_face] +
                        W[1][idx_face] + 1.0/(c_average[idx_face]*c_average[idx_face])*W[4][idx_face];
                    V[1][idx_face] = W[2][idx_face];
                    V[2][idx_face] = -1.0/(rho_average[idx_face]*c_average[idx_face])*W[0][idx_face] +
                        1.0/(rho_average[idx_face]*c_average[idx_face])*W[4][idx_face];
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
                    
                    V[0][idx_face] = 1.0/(c_average[idx_face]*c_average[idx_face])*W[0][idx_face] +
                        W[1][idx_face] + 1.0/(c_average[idx_face]*c_average[idx_face])*W[4][idx_face];
                    V[1][idx_face] = W[2][idx_face];
                    V[2][idx_face] = W[3][idx_face];
                    V[3][idx_face] = -1.0/(rho_average[idx_face]*c_average[idx_face])*W[0][idx_face] +
                        1.0/(rho_average[idx_face]*c_average[idx_face])*W[4][idx_face];
                    V[4][idx_face] = W[0][idx_face] + W[4][idx_face];
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
 * The FlowModelSingleSpecies class modifies nothing for velocity_face.
 */
void
FlowModelSingleSpecies::computeLocalFaceFluxAndVelocityFromRiemannSolverWithConservativeVariables(
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
            << ": FlowModelSingleSpecies::"
            << "computeLocalIntercellQuantitiesFromRiemannSolverWithConservativeVariables()\n"
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
 * The FlowModelSingleSpecies class modifies nothing for velocity_face.
 */
void
FlowModelSingleSpecies::computeLocalFaceFluxAndVelocityFromRiemannSolverWithPrimitiveVariables(
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
            << ": FlowModelSingleSpecies::"
            << "computeLocalIntercellQuantitiesFromRiemannSolverWithPrimitiveVariables()\n"
            << "Unknown Riemann solver required."
            << std::endl);
        }
    }
}


/*
 * Check whether the given side conservative variables are within the bounds.
 */
void
FlowModelSingleSpecies::checkGlobalSideDataConservativeVariablesBounded(
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
                << ": FlowModelSingleSpecies::"
                << "checkGlobalSideDataConservativeVariablesBounded()\n"
                << "The interior dimension of the conservative variables does not match that of patch."
                << std::endl);
        }
    }
    const hier::IntVector interior_dims_flag = bounded_flag->getBox().numberCells();
    if (interior_dims_flag != d_interior_dims)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "checkGlobalSideDataConservativeVariablesBounded()\n"
            << "The interior dimension of the flag does not match that of patch."
            << std::endl);
    }
    
    for (int ei = 1; ei < d_num_eqn; ei++)
    {
        if (num_ghosts_conservative_var != conservative_variables[ei]->getGhostCellWidth())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::"
                << "checkGlobalSideDataConservativeVariablesBounded()\n"
                << "The conservative variables don't have same ghost cell width."
                << std::endl);
        }
    }
    
    if (num_ghosts_flag != num_ghosts_conservative_var)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "checkGlobalSideDataConservativeVariablesBounded()\n"
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
            
            if (Q[0][idx_face] > 0.0)
            {
                are_bounded[idx_face] &= 1;
            }
            else
            {
                are_bounded[idx_face] &= 0;
            }
            
            if (Q[d_num_eqn - 1][idx_face] > 0.0)
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
                
                if (Q[0][idx_face] > 0.0)
                {
                    are_bounded[idx_face] &= 1;
                }
                else
                {
                    are_bounded[idx_face] &= 0;
                }
                
                if (Q[d_num_eqn - 1][idx_face] > 0.0)
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
                
                if (Q[0][idx_face] > 0.0)
                {
                    are_bounded[idx_face] &= 1;
                }
                else
                {
                    are_bounded[idx_face] &= 0;
                }
                
                if (Q[d_num_eqn - 1][idx_face] > 0.0)
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
                    
                    if (Q[0][idx_face] > 0.0)
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                    
                    if (Q[d_num_eqn - 1][idx_face] > 0.0)
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
                    
                    if (Q[0][idx_face] > 0.0)
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                    
                    if (Q[d_num_eqn - 1][idx_face] > 0.0)
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
                    
                    if (Q[0][idx_face] > 0.0)
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                    
                    if (Q[d_num_eqn - 1][idx_face] > 0.0)
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
FlowModelSingleSpecies::checkGlobalSideDataPrimitiveVariablesBounded(
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
                << ": FlowModelSingleSpecies::"
                << "checkGlobalSideDataPrimitiveVariablesBounded()\n"
                << "The interior dimension of the primitive variables does not match that of patch."
                << std::endl);
        }
    }
    const hier::IntVector interior_dims_flag = bounded_flag->getBox().numberCells();
    if (interior_dims_flag != d_interior_dims)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "checkGlobalSideDataPrimitiveVariablesBounded()\n"
            << "The interior dimension of the flag does not match that of patch."
            << std::endl);
    }
    
    for (int ei = 1; ei < d_num_eqn; ei++)
    {
        if (num_ghosts_primitive_var != primitive_variables[ei]->getGhostCellWidth())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::"
                << "checkGlobalSideDataPrimitiveVariablesBounded()\n"
                << "The primitive variables don't have same ghost cell width."
                << std::endl);
        }
    }
    
    if (num_ghosts_flag != num_ghosts_primitive_var)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::"
            << "checkGlobalSideDataPrimitiveVariablesBounded()\n"
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
            
            if (V[0][idx_face] > 0.0)
            {
                are_bounded[idx_face] &= 1;
            }
            else
            {
                are_bounded[idx_face] &= 0;
            }
            
            if (V[d_num_eqn - 1][idx_face] > 0.0)
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
                
                if (V[0][idx_face] > 0.0)
                {
                    are_bounded[idx_face] &= 1;
                }
                else
                {
                    are_bounded[idx_face] &= 0;
                }
                
                if (V[d_num_eqn - 1][idx_face] > 0.0)
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
                
                if (V[0][idx_face] > 0.0)
                {
                    are_bounded[idx_face] &= 1;
                }
                else
                {
                    are_bounded[idx_face] &= 0;
                }
                
                if (V[d_num_eqn - 1][idx_face] > 0.0)
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
                    
                    if (V[0][idx_face] > 0.0)
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                    
                    if (V[d_num_eqn - 1][idx_face] > 0.0)
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
                    
                    if (V[0][idx_face] > 0.0)
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                    
                    if (V[d_num_eqn - 1][idx_face] > 0.0)
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
                    
                    if (V[0][idx_face] > 0.0)
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                    
                    if (V[d_num_eqn - 1][idx_face] > 0.0)
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
FlowModelSingleSpecies::convertLocalCellDataPointersConservativeVariablesToPrimitiveVariables(
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
    
    double epsilon = 0.0;
    if (d_dim == tbox::Dimension(1))
    {
        epsilon = (*Q[2] - 0.5*(*Q[1])*(*Q[1])/(*Q[0]))/(*Q[0]);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        epsilon = (*Q[3] - 0.5*((*Q[1])*(*Q[1]) + (*Q[2])*(*Q[2]))/(*Q[0]))/(*Q[0]);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        epsilon = (*Q[4] - 0.5*((*Q[1])*(*Q[1]) + (*Q[2])*(*Q[2]) + (*Q[3])*(*Q[3]))/(*Q[0]))/(*Q[0]);
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
 * Convert vector of pointers of primitive cell data to vectors of pointers of conservative cell data.
 */
void
FlowModelSingleSpecies::convertLocalCellDataPointersPrimitiveVariablesToConservativeVariables(
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
    
    // Compute the specific internal energy.
    const double epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
        getInternalEnergy(
            V[0],
            V[1 + d_dim.getValue()],
            thermo_properties_ptr);
    
    double E = 0.0;
    if (d_dim == tbox::Dimension(1))
    {
        E = (*V[0])*(epsilon + 0.5*(*V[1])*(*V[1]));
    }
    else if (d_dim == tbox::Dimension(2))
    {
        E = (*V[0])*(epsilon + 0.5*((*V[1])*(*V[1]) + (*V[2])*(*V[2])));
    }
    else if (d_dim == tbox::Dimension(3))
    {
        E = (*V[0])*(epsilon + 0.5*((*V[1])*(*V[1]) + (*V[2])*(*V[2]) + (*V[3])*(*V[3])));
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
    derivative_var_data.resize(d_num_eqn);
    derivative_var_component_idx.resize(d_num_eqn);
    
    if (!d_data_velocity)
    {
        computeGlobalCellDataVelocity();
    }
    
    if (!d_data_temperature)
    {
        computeGlobalCellDataTemperatureWithPressure();
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
            computeGlobalCellDataVelocity();
        }
        
        if (!d_data_pressure)
        {
            computeGlobalCellDataPressureWithInternalEnergy();
        }
        
        if (!d_data_temperature)
        {
            computeGlobalCellDataTemperatureWithPressure();
        }
        
        // Get the pointers to the cell data of pressure and temperature.
        double* p = d_data_pressure->getPointer(0);
        double* T = d_data_temperature->getPointer(0);
        
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
        
        /*
         * Compute shear viscosity.
         */
        
        boost::shared_ptr<pdat::CellData<double> > data_shear_viscosity(
            new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_diffusivities));
        
        double* mu = data_shear_viscosity->getPointer(0);
        
        if (d_dim == tbox::Dimension(1))
        {
            for (int i = -d_num_subghosts_diffusivities[0];
                 i < d_interior_dims[0] + d_num_subghosts_diffusivities[0];
                 i++)
            {
                // Compute the linear indices.
                const int idx_diffusivities = i + d_num_subghosts_diffusivities[0];
                const int idx_pressure = i + d_num_subghosts_pressure[0];
                const int idx_temperature = i + d_num_subghosts_temperature[0];
                
                mu[idx_diffusivities] = d_equation_of_shear_viscosity_mixing_rules->
                    getEquationOfShearViscosity()->
                        getShearViscosity(
                            &p[idx_pressure],
                            &T[idx_temperature],
                            molecular_properties_shear_viscosity_ptr);
            }
        }
        else if (d_dim == tbox::Dimension(2))
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
                        (j + d_num_subghosts_diffusivities[1])*d_subghostcell_dims_diffusivities[0];
                    
                    const int idx_pressure = (i + d_num_subghosts_pressure[0]) +
                        (j + d_num_subghosts_pressure[1])*d_subghostcell_dims_pressure[0];
                    
                    const int idx_temperature = (i + d_num_subghosts_temperature[0]) +
                        (j + d_num_subghosts_temperature[1])*d_subghostcell_dims_temperature[0];
                    
                    mu[idx_diffusivities] = d_equation_of_shear_viscosity_mixing_rules->
                        getEquationOfShearViscosity()->
                            getShearViscosity(
                                &p[idx_pressure],
                                &T[idx_temperature],
                                molecular_properties_shear_viscosity_ptr);
                }
            }
        }
        else if (d_dim == tbox::Dimension(3))
        {
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
                        
                        const int idx_pressure = (i + d_num_subghosts_pressure[0]) +
                            (j + d_num_subghosts_pressure[1])*d_subghostcell_dims_pressure[0] +
                            (k + d_num_subghosts_pressure[2])*d_subghostcell_dims_pressure[0]*
                                d_subghostcell_dims_pressure[1];
                        
                        const int idx_temperature = (i + d_num_subghosts_temperature[0]) +
                            (j + d_num_subghosts_temperature[1])*d_subghostcell_dims_temperature[0] +
                            (k + d_num_subghosts_temperature[2])*d_subghostcell_dims_temperature[0]*
                                d_subghostcell_dims_temperature[1];
                        
                        mu[idx_diffusivities] = d_equation_of_shear_viscosity_mixing_rules->
                            getEquationOfShearViscosity()->
                                getShearViscosity(
                                    &p[idx_pressure],
                                    &T[idx_temperature],
                                    molecular_properties_shear_viscosity_ptr);
                    }
                }
            }
        }
        
        /*
         * Compute bulk viscosity.
         */
        
        boost::shared_ptr<pdat::CellData<double> > data_bulk_viscosity(
            new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_diffusivities));
        
        double* mu_v = data_bulk_viscosity->getPointer(0);
        
        if (d_dim == tbox::Dimension(1))
        {
            for (int i = -d_num_subghosts_diffusivities[0];
                 i < d_interior_dims[0] + d_num_subghosts_diffusivities[0];
                 i++)
            {
                // Compute the linear indices.
                const int idx_diffusivities = i + d_num_subghosts_diffusivities[0];
                const int idx_pressure = i + d_num_subghosts_pressure[0];
                const int idx_temperature = i + d_num_subghosts_temperature[0];
                
                mu_v[idx_diffusivities] = d_equation_of_bulk_viscosity_mixing_rules->
                    getEquationOfBulkViscosity()->
                        getBulkViscosity(
                            &p[idx_pressure],
                            &T[idx_temperature],
                            molecular_properties_bulk_viscosity_ptr);
            }
        }
        else if (d_dim == tbox::Dimension(2))
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
                        (j + d_num_subghosts_diffusivities[1])*d_subghostcell_dims_diffusivities[0];
                    
                    const int idx_pressure = (i + d_num_subghosts_pressure[0]) +
                        (j + d_num_subghosts_pressure[1])*d_subghostcell_dims_pressure[0];
                    
                    const int idx_temperature = (i + d_num_subghosts_temperature[0]) +
                        (j + d_num_subghosts_temperature[1])*d_subghostcell_dims_temperature[0];
                    
                    mu_v[idx_diffusivities] = d_equation_of_bulk_viscosity_mixing_rules->
                        getEquationOfBulkViscosity()->
                            getBulkViscosity(
                                &p[idx_pressure],
                                &T[idx_temperature],
                                molecular_properties_bulk_viscosity_ptr);
                }
            }
        }
        else if (d_dim == tbox::Dimension(3))
        {
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
                        
                        const int idx_pressure = (i + d_num_subghosts_pressure[0]) +
                            (j + d_num_subghosts_pressure[1])*d_subghostcell_dims_pressure[0] +
                            (k + d_num_subghosts_pressure[2])*d_subghostcell_dims_pressure[0]*
                                d_subghostcell_dims_pressure[1];
                        
                        const int idx_temperature = (i + d_num_subghosts_temperature[0]) +
                            (j + d_num_subghosts_temperature[1])*d_subghostcell_dims_temperature[0] +
                            (k + d_num_subghosts_temperature[2])*d_subghostcell_dims_temperature[0]*
                                d_subghostcell_dims_temperature[1];
                        
                        mu_v[idx_diffusivities] = d_equation_of_bulk_viscosity_mixing_rules->
                            getEquationOfBulkViscosity()->
                                getBulkViscosity(
                                    &p[idx_pressure],
                                    &T[idx_temperature],
                                    molecular_properties_bulk_viscosity_ptr);
                    }
                }
            }
        }
        
        /*
         * Compute thermal conductivity.
         */
        
        boost::shared_ptr<pdat::CellData<double> > data_thermal_conductivity(new pdat::CellData<double>(
            d_interior_box, 1, d_num_subghosts_diffusivities));
        
        double* kappa = data_thermal_conductivity->getPointer(0);
        
        if (d_dim == tbox::Dimension(1))
        {
            for (int i = -d_num_subghosts_diffusivities[0];
                 i < d_interior_dims[0] + d_num_subghosts_diffusivities[0];
                 i++)
            {
                // Compute the linear indices.
                const int idx_diffusivities = i + d_num_subghosts_diffusivities[0];
                const int idx_pressure = i + d_num_subghosts_pressure[0];
                const int idx_temperature = i + d_num_subghosts_temperature[0];
                
                kappa[idx_diffusivities] = d_equation_of_thermal_conductivity_mixing_rules->
                    getEquationOfThermalConductivity()->
                        getThermalConductivity(
                            &p[idx_pressure],
                            &T[idx_temperature],
                            molecular_properties_thermal_conductivity_ptr);
            }
        }
        else if (d_dim == tbox::Dimension(2))
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
                        (j + d_num_subghosts_diffusivities[1])*d_subghostcell_dims_diffusivities[0];
                    
                    const int idx_pressure = (i + d_num_subghosts_pressure[0]) +
                        (j + d_num_subghosts_pressure[1])*d_subghostcell_dims_pressure[0];
                    
                    const int idx_temperature = (i + d_num_subghosts_temperature[0]) +
                        (j + d_num_subghosts_temperature[1])*d_subghostcell_dims_temperature[0];
                    
                    kappa[idx_diffusivities] = d_equation_of_thermal_conductivity_mixing_rules->
                        getEquationOfThermalConductivity()->
                            getThermalConductivity(
                                &p[idx_pressure],
                                &T[idx_temperature],
                                molecular_properties_thermal_conductivity_ptr);
                }
            }
        }
        else if (d_dim == tbox::Dimension(3))
        {
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
                        
                        const int idx_pressure = (i + d_num_subghosts_pressure[0]) +
                            (j + d_num_subghosts_pressure[1])*d_subghostcell_dims_pressure[0] +
                            (k + d_num_subghosts_pressure[2])*d_subghostcell_dims_pressure[0]*
                                d_subghostcell_dims_pressure[1];
                        
                        const int idx_temperature = (i + d_num_subghosts_temperature[0]) +
                            (j + d_num_subghosts_temperature[1])*d_subghostcell_dims_temperature[0] +
                            (k + d_num_subghosts_temperature[2])*d_subghostcell_dims_temperature[0]*
                                d_subghostcell_dims_temperature[1];
                        
                        kappa[idx_diffusivities] = d_equation_of_thermal_conductivity_mixing_rules->
                            getEquationOfThermalConductivity()->
                                getThermalConductivity(
                                    &p[idx_pressure],
                                    &T[idx_temperature],
                                    molecular_properties_thermal_conductivity_ptr);
                    }
                }
            }
        }
        
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
                
                D_00[idx_diffusivities] = -(4.0/3.0*mu[idx_diffusivities] + mu_v[idx_diffusivities]);
                D_01[idx_diffusivities] = -u[idx_velocity]*(4.0/3.0*mu[idx_diffusivities] +
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
                    
                    D_00[idx_diffusivities] = -(4.0/3.0*mu[idx_diffusivities] + mu_v[idx_diffusivities]);
                    D_01[idx_diffusivities] = 2.0/3.0*mu[idx_diffusivities] - mu_v[idx_diffusivities];
                    D_02[idx_diffusivities] = -mu[idx_diffusivities];
                    D_03[idx_diffusivities] = -u[idx_velocity]*(4.0/3.0*mu[idx_diffusivities] +
                        mu_v[idx_diffusivities]);
                    D_04[idx_diffusivities] = -v[idx_velocity]*(4.0/3.0*mu[idx_diffusivities] +
                        mu_v[idx_diffusivities]);
                    D_05[idx_diffusivities] = u[idx_velocity]*(2.0/3.0*mu[idx_diffusivities] -
                        mu_v[idx_diffusivities]);
                    D_06[idx_diffusivities] = v[idx_velocity]*(2.0/3.0*mu[idx_diffusivities] -
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
                        
                        D_00[idx_diffusivities] = -(4.0/3.0*mu[idx_diffusivities] + mu_v[idx_diffusivities]);
                        D_01[idx_diffusivities] = 2.0/3.0*mu[idx_diffusivities] - mu_v[idx_diffusivities];
                        D_02[idx_diffusivities] = -mu[idx_diffusivities];
                        D_03[idx_diffusivities] = -u[idx_velocity]*(4.0/3.0*mu[idx_diffusivities] +
                            mu_v[idx_diffusivities]);
                        D_04[idx_diffusivities] = -v[idx_velocity]*(4.0/3.0*mu[idx_diffusivities] +
                            mu_v[idx_diffusivities]);
                        D_05[idx_diffusivities] = -w[idx_velocity]*(4.0/3.0*mu[idx_diffusivities] +
                            mu_v[idx_diffusivities]);
                        D_06[idx_diffusivities] = u[idx_velocity]*(2.0/3.0*mu[idx_diffusivities] -
                            mu_v[idx_diffusivities]);
                        D_07[idx_diffusivities] = v[idx_velocity]*(2.0/3.0*mu[idx_diffusivities] -
                            mu_v[idx_diffusivities]);
                        D_08[idx_diffusivities] = w[idx_velocity]*(2.0/3.0*mu[idx_diffusivities] -
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
                patch.getPatchData(d_variable_density, d_plot_context)));
        
        boost::shared_ptr<pdat::CellData<double> > data_momentum(
            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                patch.getPatchData(d_variable_momentum, d_plot_context)));
        
        boost::shared_ptr<pdat::CellData<double> > data_total_energy(
            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                patch.getPatchData(d_variable_total_energy, d_plot_context)));
        
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
                
                const double epsilon = (E[idx_data] - 0.5*rho_u[idx_data]*rho_u[idx_data]/rho[idx_data])/rho[idx_data];
                
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
                    
                    const double epsilon = (E[idx_data] - 0.5*(rho_u[idx_data]*rho_u[idx_data] +
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
                        
                        const double epsilon = (E[idx_data] - 0.5*(rho_u[idx_data]*rho_u[idx_data] +
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
                patch.getPatchData(d_variable_density, d_plot_context)));
        
        boost::shared_ptr<pdat::CellData<double> > data_momentum(
            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                patch.getPatchData(d_variable_momentum, d_plot_context)));
        
        boost::shared_ptr<pdat::CellData<double> > data_total_energy(
            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                patch.getPatchData(d_variable_total_energy, d_plot_context)));
        
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
                
                const double epsilon = (E[idx_data] - 0.5*rho_u[idx_data]*rho_u[idx_data]/
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
                    
                    const double epsilon = (E[idx_data] - 0.5*(rho_u[idx_data]*rho_u[idx_data] +
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
                        
                        const double epsilon = (E[idx_data] - 0.5*(rho_u[idx_data]*rho_u[idx_data] +
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
                patch.getPatchData(d_variable_density, d_plot_context)));
        
        boost::shared_ptr<pdat::CellData<double> > data_momentum(
            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                patch.getPatchData(d_variable_momentum, d_plot_context)));
        
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
    const boost::shared_ptr<appu::VisItDataWriter>& visit_writer)
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
           d_variable_density,
           d_plot_context));
    
    /*
    visit_writer->registerPlotQuantity(
        "momentum",
        "VECTOR",
        vardb->mapVariableAndContextToIndex(
           d_variable_momentum,
           d_plot_context));
    
    visit_writer->registerPlotQuantity(
        "total energy",
        "SCALAR",
        vardb->mapVariableAndContextToIndex(
           d_variable_total_energy,
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
    if (variable_name == "VELOCITY")
    {
        if (d_num_subghosts_velocity > -hier::IntVector::getOne(d_dim))
        {
            if (num_subghosts > d_num_subghosts_velocity)
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::setNumberOfSubGhosts()\n"
                    << "Number of ghosts of '"
                    << parent_variable_name
                    << "' exceeds"
                    << " number of ghosts of '"
                    << variable_name
                    << "'."
                    << std::endl);
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
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::setNumberOfSubGhosts()\n"
                    << "Number of ghosts of '"
                    << parent_variable_name
                    << "' exceeds"
                    << " number of ghosts of '"
                    << variable_name
                    << "'."
                    << std::endl);
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
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::setNumberOfSubGhosts()\n"
                    << "Number of ghosts of '"
                    << parent_variable_name
                    << "' exceeds"
                    << " number of ghosts of '"
                    << variable_name
                    << "'."
                    << std::endl);
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
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::setNumberOfSubGhosts()\n"
                    << "Number of ghosts of '"
                    << parent_variable_name
                    << "' exceeds"
                    << " number of ghosts of '"
                    << variable_name
                    << "'."
                    << std::endl);
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
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::setNumberOfSubGhosts()\n"
                    << "Number of ghosts of '"
                    << parent_variable_name
                    << "' exceeds"
                    << " number of ghosts of '"
                    << variable_name
                    << "'."
                    << std::endl);
            }
        }
        else
        {
            d_num_subghosts_temperature = num_subghosts;
        }
        
        setNumberOfSubGhosts(num_subghosts, "PRESSURE", parent_variable_name);
    }
    else if (variable_name == "DILATATION")
    {
        if (d_num_subghosts_dilatation > -hier::IntVector::getOne(d_dim))
        {
            if (num_subghosts > d_num_subghosts_dilatation)
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::setNumberOfSubGhosts()\n"
                    << "Number of ghosts of '"
                    << parent_variable_name
                    << "' exceeds"
                    << " number of ghosts of '"
                    << variable_name
                    << "'."
                    << std::endl);
            }
        }
        else
        {
            d_num_subghosts_dilatation = num_subghosts;
        }
        
        setNumberOfSubGhosts(num_subghosts, "VELOCITY", parent_variable_name);
    }
    else if (variable_name == "VORTICITY")
    {
        if (d_num_subghosts_vorticity > -hier::IntVector::getOne(d_dim))
        {
            if (num_subghosts > d_num_subghosts_vorticity)
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::setNumberOfSubGhosts()\n"
                    << "Number of ghosts of '"
                    << parent_variable_name
                    << "' exceeds"
                    << " number of ghosts of '"
                    << variable_name
                    << "'."
                    << std::endl);
            }
        }
        else
        {
            d_num_subghosts_vorticity = num_subghosts;
        }
        
        setNumberOfSubGhosts(num_subghosts, "VELOCITY", parent_variable_name);
    }
    else if (variable_name == "ENSTROPHY")
    {
        if (d_num_subghosts_enstrophy > -hier::IntVector::getOne(d_dim))
        {
            if (num_subghosts > d_num_subghosts_enstrophy)
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::setNumberOfSubGhosts()\n"
                    << "Number of ghosts of '"
                    << parent_variable_name
                    << "' exceeds"
                    << " number of ghosts of '"
                    << variable_name
                    << "'."
                    << std::endl);
            }
        }
        else
        {
            d_num_subghosts_enstrophy = num_subghosts;
        }
        
        setNumberOfSubGhosts(num_subghosts, "VORTICITY", parent_variable_name);
    }
    else if (variable_name == "CONVECTIVE_FLUX_X")
    {
        if (d_num_subghosts_convective_flux_x > -hier::IntVector::getOne(d_dim))
        {
            if (num_subghosts > d_num_subghosts_convective_flux_x)
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::setNumberOfSubGhosts()\n"
                    << "Number of ghosts of '"
                    << parent_variable_name
                    << "' exceeds"
                    << " number of ghosts of '"
                    << variable_name
                    << "'."
                    << std::endl);
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
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::setNumberOfSubGhosts()\n"
                    << "Number of ghosts of '"
                    << parent_variable_name
                    << "' exceeds"
                    << " number of ghosts of '"
                    << variable_name
                    << "'."
                    << std::endl);
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
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::setNumberOfSubGhosts()\n"
                    << "Number of ghosts of '"
                    << parent_variable_name
                    << "' exceeds"
                    << " number of ghosts of '"
                    << variable_name
                    << "'."
                    << std::endl);
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
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::setNumberOfSubGhosts()\n"
                    << "Number of ghosts of '"
                    << parent_variable_name
                    << "' exceeds"
                    << " number of ghosts of '"
                    << variable_name
                    << "'."
                    << std::endl);
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
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::setNumberOfSubGhosts()\n"
                    << "Number of ghosts of '"
                    << parent_variable_name
                    << "' exceeds"
                    << " number of ghosts of '"
                    << variable_name
                    << "'."
                    << std::endl);
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
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::setNumberOfSubGhosts()\n"
                    << "Number of ghosts of '"
                    << parent_variable_name
                    << "' exceeds"
                    << " number of ghosts of '"
                    << variable_name
                    << "'."
                    << std::endl);
            }
        }
        else
        {
            d_num_subghosts_max_wave_speed_z = num_subghosts;
        }
        
        setNumberOfSubGhosts(num_subghosts, "VELOCITY", parent_variable_name);
        setNumberOfSubGhosts(num_subghosts, "SOUND_SPEED", parent_variable_name);
    }
}


/*
 * Set the ghost boxes and their dimensions of derived cell variables.
 */
void
FlowModelSingleSpecies::setGhostBoxesAndDimensionsDerivedCellVariables()
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
    
    if (d_num_subghosts_dilatation > -hier::IntVector::getOne(d_dim))
    {
        d_subghost_box_dilatation = d_interior_box;
        d_subghost_box_dilatation.grow(d_num_subghosts_dilatation);
        d_subghostcell_dims_dilatation = d_subghost_box_dilatation.numberCells();
    }
    
    if (d_num_subghosts_vorticity > -hier::IntVector::getOne(d_dim))
    {
        d_subghost_box_vorticity = d_interior_box;
        d_subghost_box_vorticity.grow(d_num_subghosts_vorticity);
        d_subghostcell_dims_vorticity = d_subghost_box_vorticity.numberCells();
    }
    
    if (d_num_subghosts_enstrophy > -hier::IntVector::getOne(d_dim))
    {
        d_subghost_box_enstrophy = d_interior_box;
        d_subghost_box_enstrophy.grow(d_num_subghosts_enstrophy);
        d_subghostcell_dims_enstrophy = d_subghost_box_enstrophy.numberCells();
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
    
    if (d_num_subghosts_diffusivities > -hier::IntVector::getOne(d_dim))
    {
        d_subghost_box_diffusivities = d_interior_box;
        d_subghost_box_diffusivities.grow(d_num_subghosts_diffusivities);
        d_subghostcell_dims_diffusivities = d_subghost_box_diffusivities.numberCells();
    }
}


/*
 * Get the global cell data of density in the registered patch.
 */
boost::shared_ptr<pdat::CellData<double> >
FlowModelSingleSpecies::getGlobalCellDataDensity()
{
    // Get the cell data of the registered variable density.
    boost::shared_ptr<pdat::CellData<double> > data_density(
        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
            d_patch->getPatchData(d_variable_density, getDataContext())));
    
    return data_density;
}


/*
 * Get the global cell data of momentum in the registered patch.
 */
boost::shared_ptr<pdat::CellData<double> >
FlowModelSingleSpecies::getGlobalCellDataMomentum()
{
    // Get the cell data of the registered variable momentum.
    boost::shared_ptr<pdat::CellData<double> > data_momentum(
        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
            d_patch->getPatchData(d_variable_momentum, getDataContext())));
    
    return data_momentum;
}


/*
 * Get the global cell data of total energy in the registered patch.
 */
boost::shared_ptr<pdat::CellData<double> >
FlowModelSingleSpecies::getGlobalCellDataTotalEnergy()
{
    // Get the cell data of the registered variable total energy.
    boost::shared_ptr<pdat::CellData<double> > data_total_energy(
        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
            d_patch->getPatchData(d_variable_total_energy, getDataContext())));
    
    return data_total_energy;
}


/*
 * Compute the global cell data of velocity in the registered patch.
 */
void
FlowModelSingleSpecies::computeGlobalCellDataVelocity()
{
    if (d_num_subghosts_velocity > -hier::IntVector::getOne(d_dim))
    {
        // Create the cell data of velocity.
        d_data_velocity.reset(
            new pdat::CellData<double>(d_interior_box, d_dim.getValue(), d_num_subghosts_velocity));
        
        // Get the cell data of the variables density and momentum.
        boost::shared_ptr<pdat::CellData<double> > data_density =
            getGlobalCellDataDensity();
        
        boost::shared_ptr<pdat::CellData<double> > data_momentum =
            getGlobalCellDataMomentum();
        
        // Get the pointer to the cell data of density.
        double* rho = data_density->getPointer(0);
        
        if (d_dim == tbox::Dimension(1))
        {
            // Get the pointer to the cell data of velocity.
            double* u = d_data_velocity->getPointer(0);
            
            // Get the pointer to the cell data of momentum.
            double* rho_u = data_momentum->getPointer(0);
            
            // Compute the velocity field.
            for (int i = -d_num_subghosts_velocity[0];
                 i < d_interior_dims[0] + d_num_subghosts_velocity[0];
                 i++)
            {
                // Compute the linear indices.
                const int idx = i + d_num_ghosts[0];
                const int idx_velocity = i + d_num_subghosts_velocity[0];
                
                u[idx_velocity] = rho_u[idx]/rho[idx];
            }
        }
        else if (d_dim == tbox::Dimension(2))
        {
            // Get the pointers to the cell data of velocity.
            double* u = d_data_velocity->getPointer(0);
            double* v = d_data_velocity->getPointer(1);
            
            // Get the pointers to the cell data of momentum.
            double* rho_u = data_momentum->getPointer(0);
            double* rho_v = data_momentum->getPointer(1);
            
            // Compute the velocity field.
            for (int j = -d_num_subghosts_velocity[1];
                 j < d_interior_dims[1] + d_num_subghosts_velocity[1];
                 j++)
            {
                for (int i = -d_num_subghosts_velocity[0];
                     i < d_interior_dims[0] + d_num_subghosts_velocity[0];
                     i++)
                {
                    // Compute the linear indices.
                    const int idx = (i + d_num_ghosts[0]) +
                        (j + d_num_ghosts[1])*d_ghostcell_dims[0];
                    
                    const int idx_velocity = (i + d_num_subghosts_velocity[0]) +
                        (j + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0];
                    
                    u[idx_velocity] = rho_u[idx]/rho[idx];
                    v[idx_velocity] = rho_v[idx]/rho[idx];
                }
            }
        }
        else if (d_dim == tbox::Dimension(3))
        {
            // Get the pointers to the cell data of velocity.
            double* u = d_data_velocity->getPointer(0);
            double* v = d_data_velocity->getPointer(1);
            double* w = d_data_velocity->getPointer(2);
            
            // Get the pointers to the cell data of momentum.
            double* rho_u = data_momentum->getPointer(0);
            double* rho_v = data_momentum->getPointer(1);
            double* rho_w = data_momentum->getPointer(2);
            
            // Compute the velocity field.
            for (int k = -d_num_subghosts_velocity[2];
                 k < d_interior_dims[2] + d_num_subghosts_velocity[2];
                 k++)
            {
                for (int j = -d_num_subghosts_velocity[1];
                     j < d_interior_dims[1] + d_num_subghosts_velocity[1];
                     j++)
                {
                    for (int i = -d_num_subghosts_velocity[0];
                         i < d_interior_dims[0] + d_num_subghosts_velocity[0];
                         i++)
                    {
                        // Compute the linear indices.
                        const int idx = (i + d_num_ghosts[0]) +
                            (j + d_num_ghosts[1])*d_ghostcell_dims[0] +
                            (k + d_num_ghosts[2])*d_ghostcell_dims[0]*d_ghostcell_dims[1];
                        
                        const int idx_velocity = (i + d_num_subghosts_velocity[0]) +
                            (j + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0] +
                            (k + d_num_subghosts_velocity[2])*d_subghostcell_dims_velocity[0]*d_subghostcell_dims_velocity[1];
                        
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
            << ": FlowModelSingleSpecies::computeGlobalCellDataVelocity()\n"
            << "Cell data of 'VELOCITY' is not yet registered."
            << std::endl);
    }
}


/*
 * Compute the global cell data of internal energy with velocity in the registered patch.
 */
void
FlowModelSingleSpecies::computeGlobalCellDataInternalEnergyWithVelocity()
{
    if (d_num_subghosts_internal_energy > -hier::IntVector::getOne(d_dim))
    {
        // Create the cell data of internal energy.
        d_data_internal_energy.reset(
            new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_internal_energy));
        
        // Get the cell data of the variables density, momentum and total energy.
        boost::shared_ptr<pdat::CellData<double> > data_density =
            getGlobalCellDataDensity();
        
        boost::shared_ptr<pdat::CellData<double> > data_total_energy =
            getGlobalCellDataTotalEnergy();
        
        if (!d_data_velocity)
        {
            computeGlobalCellDataVelocity();
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
            // Get the pointer to cell data of velocity.
            double* u = d_data_velocity->getPointer(0);
            
            // Compute the internal energy field.
            for (int i = -d_num_subghosts_pressure[0];
                 i < d_interior_dims[0] + d_num_subghosts_pressure[0];
                 i++)
            {
                // Compute the linear indices.
                const int idx = i + d_num_ghosts[0];
                const int idx_velocity = i + d_num_subghosts_velocity[0];
                const int idx_internal_energy = i + d_num_subghosts_internal_energy[0];
                
                epsilon[idx_internal_energy] =
                    E[idx]/rho[idx] - 0.5*u[idx_velocity]*u[idx_velocity];
            }
        }
        else if (d_dim == tbox::Dimension(2))
        {
            // Get the pointers to the cell data of velocity.
            double* u = d_data_velocity->getPointer(0);
            double* v = d_data_velocity->getPointer(1);
            
            // Compute the internal energy field.
            for (int j = -d_num_subghosts_pressure[1];
                 j < d_interior_dims[1] + d_num_subghosts_pressure[1];
                 j++)
            {
                for (int i = -d_num_subghosts_pressure[0];
                     i < d_interior_dims[0] + d_num_subghosts_pressure[0];
                     i++)
                {
                    // Compute the linear indices.
                    const int idx = (i + d_num_ghosts[0]) +
                        (j + d_num_ghosts[1])*d_ghostcell_dims[0];
                    
                    const int idx_velocity = (i + d_num_subghosts_velocity[0]) +
                        (j + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0];
                    
                    const int idx_internal_energy = (i + d_num_subghosts_internal_energy[0]) +
                        (j + d_num_subghosts_internal_energy[1])*d_subghostcell_dims_internal_energy[0];
                    
                    epsilon[idx_internal_energy] =
                        E[idx]/rho[idx] - 0.5*(u[idx_velocity]*u[idx_velocity] + v[idx_velocity]*v[idx_velocity]);
                }
            }
        }
        else if (d_dim == tbox::Dimension(3))
        {
            // Get the pointers to the cell data of velocity.
            double* u = d_data_velocity->getPointer(0);
            double* v = d_data_velocity->getPointer(1);
            double* w = d_data_velocity->getPointer(2);
            
            // Compute the internal energy field.
            for (int k = -d_num_subghosts_pressure[2];
                 k < d_interior_dims[2] + d_num_subghosts_pressure[2];
                 k++)
            {
                for (int j = -d_num_subghosts_pressure[1];
                     j < d_interior_dims[1] + d_num_subghosts_pressure[1];
                     j++)
                {
                    for (int i = -d_num_subghosts_pressure[0];
                         i < d_interior_dims[0] + d_num_subghosts_pressure[0];
                         i++)
                    {
                        // Compute the linear indices.
                        const int idx = (i + d_num_ghosts[0]) +
                            (j + d_num_ghosts[1])*d_ghostcell_dims[0] +
                            (k + d_num_ghosts[2])*d_ghostcell_dims[0]*d_ghostcell_dims[1];
                        
                        const int idx_velocity = (i + d_num_subghosts_velocity[0]) +
                            (j + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0] +
                            (k + d_num_subghosts_velocity[2])*d_subghostcell_dims_velocity[0]*
                                d_subghostcell_dims_velocity[1];
                        
                        const int idx_internal_energy = (i + d_num_subghosts_internal_energy[0]) +
                            (j + d_num_subghosts_internal_energy[1])*d_subghostcell_dims_internal_energy[0] +
                            (k + d_num_subghosts_internal_energy[2])*d_subghostcell_dims_internal_energy[0]*
                                d_subghostcell_dims_internal_energy[1];
                        
                        epsilon[idx_internal_energy] =
                            E[idx]/rho[idx] - 0.5*(u[idx_velocity]*u[idx_velocity] + v[idx_velocity]*v[idx_velocity] +
                                w[idx_velocity]*w[idx_velocity]);
                    }
                }
            }
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::computeGlobalCellDataPressureWithInternalEnergy()\n"
            << "Cell data of 'INTERNAL_ENERGY' is not yet registered."
            << std::endl);
    }
}


/*
 * Compute the global cell data of pressure with internal energy in the registered patch.
 */
void
FlowModelSingleSpecies::computeGlobalCellDataPressureWithInternalEnergy()
{
    if (d_num_subghosts_pressure > -hier::IntVector::getOne(d_dim))
    {
        // Create the cell data of pressure.
        d_data_pressure.reset(
            new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_pressure));
        
        // Get the cell data of the variables density, momentum and total energy.
        boost::shared_ptr<pdat::CellData<double> > data_density =
            getGlobalCellDataDensity();
        
        if (!d_data_internal_energy)
        {
            computeGlobalCellDataInternalEnergyWithVelocity();
        }
        
        // Get the pointers to the cell data of pressure, density and internal energy.
        double* p = d_data_pressure->getPointer(0);
        double* rho = data_density->getPointer(0);
        double* epsilon = d_data_internal_energy->getPointer(0);
        
        // Get the thermodynamic properties of the species.
        std::vector<const double*> thermo_properties_ptr;
        thermo_properties_ptr.reserve(static_cast<int> (d_thermo_properties.size()));
        for (int ti = 0; ti < static_cast<int> (d_thermo_properties.size()); ti++)
        {
            thermo_properties_ptr.push_back(&d_thermo_properties[ti]);
        }
        
        if (d_dim == tbox::Dimension(1))
        {
            // Compute the pressure field.
            for (int i = -d_num_subghosts_pressure[0];
                 i < d_interior_dims[0] + d_num_subghosts_pressure[0];
                 i++)
            {
                // Compute the linear indices.
                const int idx = i + d_num_ghosts[0];
                const int idx_internal_energy = i + d_num_subghosts_internal_energy[0];
                const int idx_pressure = i + d_num_subghosts_pressure[0];
                
                p[idx_pressure] = d_equation_of_state_mixing_rules->getEquationOfState()->
                    getPressure(
                        &rho[idx],
                        &epsilon[idx_internal_energy],
                        thermo_properties_ptr);
            }
        }
        else if (d_dim == tbox::Dimension(2))
        {
            // Compute the pressure field.
            for (int j = -d_num_subghosts_pressure[1];
                 j < d_interior_dims[1] + d_num_subghosts_pressure[1];
                 j++)
            {
                for (int i = -d_num_subghosts_pressure[0];
                     i < d_interior_dims[0] + d_num_subghosts_pressure[0];
                     i++)
                {
                    // Compute the linear indices.
                    const int idx = (i + d_num_ghosts[0]) +
                        (j + d_num_ghosts[1])*d_ghostcell_dims[0];
                    
                    const int idx_internal_energy = (i + d_num_subghosts_internal_energy[0]) +
                        (j + d_num_subghosts_internal_energy[1])*d_subghostcell_dims_internal_energy[0];
                    
                    const int idx_pressure = (i + d_num_subghosts_pressure[0]) +
                        (j + d_num_subghosts_pressure[1])*d_subghostcell_dims_pressure[0];
                    
                    p[idx_pressure] = d_equation_of_state_mixing_rules->getEquationOfState()->
                        getPressure(
                            &rho[idx],
                            &epsilon[idx_internal_energy],
                            thermo_properties_ptr);
                }
            }
        }
        else if (d_dim == tbox::Dimension(3))
        {
            // Compute the pressure field.
            for (int k = -d_num_subghosts_pressure[2];
                 k < d_interior_dims[2] + d_num_subghosts_pressure[2];
                 k++)
            {
                for (int j = -d_num_subghosts_pressure[1];
                     j < d_interior_dims[1] + d_num_subghosts_pressure[1];
                     j++)
                {
                    for (int i = -d_num_subghosts_pressure[0];
                         i < d_interior_dims[0] + d_num_subghosts_pressure[0];
                         i++)
                    {
                        // Compute the linear indices.
                        const int idx = (i + d_num_ghosts[0]) +
                            (j + d_num_ghosts[1])*d_ghostcell_dims[0] +
                            (k + d_num_ghosts[2])*d_ghostcell_dims[0]*d_ghostcell_dims[1];
                        
                        const int idx_pressure = (i + d_num_subghosts_pressure[0]) +
                            (j + d_num_subghosts_pressure[1])*d_subghostcell_dims_pressure[0] +
                            (k + d_num_subghosts_pressure[2])*d_subghostcell_dims_pressure[0]*d_subghostcell_dims_pressure[1];
                        
                        const int idx_internal_energy = (i + d_num_subghosts_internal_energy[0]) +
                            (j + d_num_subghosts_internal_energy[1])*d_subghostcell_dims_internal_energy[0] +
                            (k + d_num_subghosts_internal_energy[2])*d_subghostcell_dims_internal_energy[0]*
                                d_subghostcell_dims_internal_energy[1];
                        
                        p[idx_pressure] = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getPressure(
                                &rho[idx],
                                &epsilon[idx_internal_energy],
                                thermo_properties_ptr);
                    }
                }
            }
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::computeGlobalCellDataPressureWithInternalEnergy()\n"
            << "Cell data of 'PRESSURE' is not yet registered."
            << std::endl);
    }
}


/*
 * Compute the global cell data of sound speed with pressure in the registered patch.
 */
void
FlowModelSingleSpecies::computeGlobalCellDataSoundSpeedWithPressure()
{
    if (d_num_subghosts_sound_speed > -hier::IntVector::getOne(d_dim))
    {
        // Create the cell data of sound speed.
        d_data_sound_speed.reset(
            new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_sound_speed));
        
        // Get the cell data of the variable density and pressure.
        boost::shared_ptr<pdat::CellData<double> > data_density =
            getGlobalCellDataDensity();
        
        if (!d_data_pressure)
        {
            computeGlobalCellDataPressureWithInternalEnergy();
        }
        
        // Get the pointers to the cell data of sound speed, density and pressure.
        double* c   = d_data_sound_speed->getPointer(0);
        double* rho = data_density->getPointer(0);
        double* p   = d_data_pressure->getPointer(0);
        
        // Get the thermodynamic properties of the species.
        std::vector<const double*> thermo_properties_ptr;
        thermo_properties_ptr.reserve(static_cast<int> (d_thermo_properties.size()));
        for (int ti = 0; ti < static_cast<int> (d_thermo_properties.size()); ti++)
        {
            thermo_properties_ptr.push_back(&d_thermo_properties[ti]);
        }
        
        if (d_dim == tbox::Dimension(1))
        {
            // Compute the sound speed field.
            for (int i = -d_num_subghosts_sound_speed[0];
                 i < d_interior_dims[0] + d_num_subghosts_sound_speed[0];
                 i++)
            {
                // Compute the linear indices.
                const int idx = i + d_num_ghosts[0];
                const int idx_pressure = i + d_num_subghosts_pressure[0];
                const int idx_sound_speed = i + d_num_subghosts_sound_speed[0];
                
                c[idx_sound_speed] = d_equation_of_state_mixing_rules->getEquationOfState()->
                    getSoundSpeed(
                        &rho[idx],
                        &p[idx_pressure],
                        thermo_properties_ptr);
            }
        }
        else if (d_dim == tbox::Dimension(2))
        {
            // Compute the sound speed field.
            for (int j = -d_num_subghosts_sound_speed[1];
                 j < d_interior_dims[1] + d_num_subghosts_sound_speed[1];
                 j++)
            {
                for (int i = -d_num_subghosts_sound_speed[0];
                     i < d_interior_dims[0] + d_num_subghosts_sound_speed[0];
                     i++)
                {
                    // Compute the linear indices.
                    const int idx = (i + d_num_ghosts[0]) +
                        (j + d_num_ghosts[1])*d_ghostcell_dims[0];
                    
                    const int idx_pressure = (i + d_num_subghosts_pressure[0]) +
                        (j + d_num_subghosts_pressure[1])*d_subghostcell_dims_pressure[0];
                    
                    const int idx_sound_speed = (i + d_num_subghosts_sound_speed[0]) +
                        (j + d_num_subghosts_sound_speed[1])*d_subghostcell_dims_sound_speed[0];
                    
                    c[idx_sound_speed] = d_equation_of_state_mixing_rules->getEquationOfState()->
                        getSoundSpeed(
                            &rho[idx],
                            &p[idx_pressure],
                            thermo_properties_ptr);
                }
            }
        }
        else if (d_dim == tbox::Dimension(3))
        {
            // Compute the sound speed field.
            for (int k = -d_num_subghosts_sound_speed[2];
                 k < d_interior_dims[2] + d_num_subghosts_sound_speed[2];
                 k++)
            {
                for (int j = -d_num_subghosts_sound_speed[1];
                     j < d_interior_dims[1] + d_num_subghosts_sound_speed[1];
                     j++)
                {
                    for (int i = -d_num_subghosts_sound_speed[0];
                         i < d_interior_dims[0] + d_num_subghosts_sound_speed[0];
                         i++)
                    {
                        // Compute the linear indices.
                        const int idx = (i + d_num_ghosts[0]) +
                            (j + d_num_ghosts[1])*d_ghostcell_dims[0] +
                            (k + d_num_ghosts[2])*d_ghostcell_dims[0]*d_ghostcell_dims[1];
                        
                        const int idx_pressure = (i + d_num_subghosts_pressure[0]) +
                            (j + d_num_subghosts_pressure[1])*d_subghostcell_dims_pressure[0] +
                            (k + d_num_subghosts_pressure[2])*d_subghostcell_dims_pressure[0]*d_subghostcell_dims_pressure[1];
                        
                        const int idx_sound_speed = (i + d_num_subghosts_sound_speed[0]) +
                            (j + d_num_subghosts_sound_speed[1])*d_subghostcell_dims_sound_speed[0] +
                            (k + d_num_subghosts_sound_speed[2])*d_subghostcell_dims_sound_speed[0]*d_subghostcell_dims_sound_speed[1];
                        
                        c[idx_sound_speed] = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getSoundSpeed(
                                &rho[idx],
                                &p[idx_pressure],
                                thermo_properties_ptr);
                    }
                }
            }
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::computeGlobalCellDataSoundSpeedWithPressure()\n"
            << "Cell data of 'SOUND_SPEED' is not yet registered."
            << std::endl);
    }
}


/*
 * Compute the global cell data of temperature with pressure in the registered patch.
 */
void
FlowModelSingleSpecies::computeGlobalCellDataTemperatureWithPressure()
{
    if (d_num_subghosts_temperature > -hier::IntVector::getOne(d_dim))
    {
        // Create the cell data of temperature.
        d_data_temperature.reset(
            new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_temperature));
        
        // Get the cell data of the variable density and pressure.
        boost::shared_ptr<pdat::CellData<double> > data_density =
            getGlobalCellDataDensity();
        
        if (!d_data_pressure)
        {
            computeGlobalCellDataPressureWithInternalEnergy();
        }
        
        // Get the pointers to the cell data of temperature, density and pressure.
        double* T   = d_data_temperature->getPointer(0);
        double* rho = data_density->getPointer(0);
        double* p   = d_data_pressure->getPointer(0);
        
        // Get the thermodynamic properties of the species.
        std::vector<const double*> thermo_properties_ptr;
        thermo_properties_ptr.reserve(static_cast<int> (d_thermo_properties.size()));
        for (int ti = 0; ti < static_cast<int> (d_thermo_properties.size()); ti++)
        {
            thermo_properties_ptr.push_back(&d_thermo_properties[ti]);
        }
        
        if (d_dim == tbox::Dimension(1))
        {
            // Compute the temperature field.
            for (int i = -d_num_subghosts_temperature[0];
                 i < d_interior_dims[0] + d_num_subghosts_temperature[0];
                 i++)
            {
                // Compute the linear indices.
                const int idx = i + d_num_ghosts[0];
                const int idx_pressure = i + d_num_subghosts_pressure[0];
                const int idx_temperature = i + d_num_subghosts_temperature[0];
                
                T[idx_temperature] = d_equation_of_state_mixing_rules->getEquationOfState()->
                    getTemperature(
                        &rho[idx],
                        &p[idx_pressure],
                        thermo_properties_ptr);
            }
        }
        else if (d_dim == tbox::Dimension(2))
        {
            // Compute the temperature field.
            for (int j = -d_num_subghosts_temperature[1];
                 j < d_interior_dims[1] + d_num_subghosts_temperature[1];
                 j++)
            {
                for (int i = -d_num_subghosts_temperature[0];
                     i < d_interior_dims[0] + d_num_subghosts_temperature[0];
                     i++)
                {
                    // Compute the linear indices.
                    const int idx = (i + d_num_ghosts[0]) +
                        (j + d_num_ghosts[1])*d_ghostcell_dims[0];
                    
                    const int idx_pressure = (i + d_num_subghosts_pressure[0]) +
                        (j + d_num_subghosts_pressure[1])*d_subghostcell_dims_pressure[0];
                    
                    const int idx_temperature = (i + d_num_subghosts_temperature[0]) +
                        (j + d_num_subghosts_temperature[1])*d_subghostcell_dims_temperature[0];
                    
                    T[idx_temperature] = d_equation_of_state_mixing_rules->getEquationOfState()->
                        getTemperature(
                            &rho[idx],
                            &p[idx_pressure],
                            thermo_properties_ptr);
                }
            }
        }
        else if (d_dim == tbox::Dimension(3))
        {
            // Compute the temperature field.
            for (int k = -d_num_subghosts_temperature[2];
                 k < d_interior_dims[2] + d_num_subghosts_temperature[2];
                 k++)
            {
                for (int j = -d_num_subghosts_temperature[1];
                     j < d_interior_dims[1] + d_num_subghosts_temperature[1];
                     j++)
                {
                    for (int i = -d_num_subghosts_temperature[0];
                         i < d_interior_dims[0] + d_num_subghosts_temperature[0];
                         i++)
                    {
                        // Compute the linear indices.
                        const int idx = (i + d_num_ghosts[0]) +
                            (j + d_num_ghosts[1])*d_ghostcell_dims[0] +
                            (k + d_num_ghosts[2])*d_ghostcell_dims[0]*d_ghostcell_dims[1];
                        
                        const int idx_pressure = (i + d_num_subghosts_pressure[0]) +
                            (j + d_num_subghosts_pressure[1])*d_subghostcell_dims_pressure[0] +
                            (k + d_num_subghosts_pressure[2])*d_subghostcell_dims_pressure[0]*d_subghostcell_dims_pressure[1];
                        
                        const int idx_temperature = (i + d_num_subghosts_temperature[0]) +
                            (j + d_num_subghosts_temperature[1])*d_subghostcell_dims_temperature[0] +
                            (k + d_num_subghosts_temperature[2])*d_subghostcell_dims_temperature[0]*d_subghostcell_dims_temperature[1];
                        
                        T[idx_temperature] = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getTemperature(
                                &rho[idx],
                                &p[idx_pressure],
                                thermo_properties_ptr);
                    }
                }
            }
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::computeGlobalCellDataTemperatureWithPressure()\n"
            << "Cell data of 'TEMPERATURE' is not yet registered."
            << std::endl);
    }
}


/*
 * Compute the global cell data of dilatation with velocity in the registered patch.
 */
void
FlowModelSingleSpecies::computeGlobalCellDataDilatationWithVelocity()
{
    if (d_num_subghosts_dilatation > -hier::IntVector::getOne(d_dim))
    {
        // Get the grid spacing.
        const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
            BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                d_patch->getPatchGeometry()));
        
        const double* const dx = patch_geom->getDx();
        
        // Create the cell data of dilatation.
        d_data_dilatation.reset(
            new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_dilatation));
        
        if (!d_data_velocity)
        {
            computeGlobalCellDataVelocity();
        }
        
        // Get the pointer to the cell data of dilatation.
        double* theta = d_data_dilatation->getPointer(0);
        
        if (d_dim == tbox::Dimension(1))
        {
            // Get the pointer to cell data of velocity.
            double* u = d_data_velocity->getPointer(0);
            
            // Compute the dilatation field.
            if (d_num_subghosts_dilatation < d_num_subghosts_velocity)
            {
                for (int i = -d_num_subghosts_dilatation[0];
                     i < d_interior_dims[0] + d_num_subghosts_dilatation[0];
                     i++)
                {
                    // Compute indices of current and neighboring cells.
                    const int idx_x_L = i - 1 + d_num_subghosts_velocity[0];
                    const int idx_x_R = i + 1 + d_num_subghosts_velocity[0];
                    const int idx_dilatation = i + d_num_subghosts_dilatation[0];
                    
                    theta[idx_dilatation] = (u[idx_x_R] - u[idx_x_L])/(2*dx[0]);
                }
            }
            else
            {
                for (int i = -d_num_subghosts_dilatation[0];
                     i < d_interior_dims[0] + d_num_subghosts_dilatation[0];
                     i++)
                {
                    // Compute indices of current and neighboring cells.
                    const int idx = i + d_num_subghosts_velocity[0];
                    const int idx_x_L = i - 1 + d_num_subghosts_velocity[0];
                    const int idx_x_R = i + 1 + d_num_subghosts_velocity[0];
                    const int idx_dilatation = i + d_num_subghosts_dilatation[0];
                    
                    if (i == -d_num_subghosts_velocity[0])
                    {
                        theta[idx_dilatation] = (u[idx_x_R] - u[idx])/(dx[0]);
                    }
                    else if (i == d_interior_dims[0] + d_num_subghosts_velocity[0] - 1)
                    {
                        theta[idx_dilatation] = (u[idx] - u[idx_x_L])/(dx[0]);
                    }
                    else
                    {
                        theta[idx_dilatation] = (u[idx_x_R] - u[idx_x_L])/(2*dx[0]);
                    }
                }
            }
        }
        else if (d_dim == tbox::Dimension(2))
        {
            // Get the pointers to the cell data of velocity.
            double* u = d_data_velocity->getPointer(0);
            double* v = d_data_velocity->getPointer(1);
            
            // Compute the dilatation field.
            if (d_num_subghosts_dilatation < d_num_subghosts_velocity)
            {
                for (int j = -d_num_subghosts_dilatation[1];
                     j < d_interior_dims[1] + d_num_subghosts_dilatation[1];
                     j++)
                {
                    for (int i = -d_num_subghosts_dilatation[0];
                         i < d_interior_dims[0] + d_num_subghosts_dilatation[0];
                         i++)
                    {
                        // Compute indices of current and neighboring cells.
                        const int idx_x_L = (i - 1 + d_num_subghosts_velocity[0]) +
                            (j + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0];
                        
                        const int idx_x_R = (i + 1 + d_num_subghosts_velocity[0]) +
                            (j + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0];
                        
                        const int idx_y_B = (i + d_num_subghosts_velocity[0]) +
                            (j - 1 + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0];
                        
                        const int idx_y_T = (i + d_num_subghosts_velocity[0]) +
                            (j + 1 + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0];
                        
                        const int idx_dilatation = (i + d_num_subghosts_dilatation[0]) +
                            (j + d_num_subghosts_dilatation[1])*d_subghostcell_dims_dilatation[0];
                        
                        double dudx = (u[idx_x_R] - u[idx_x_L])/(2*dx[0]);
                        double dvdy = (v[idx_y_T] - v[idx_y_B])/(2*dx[1]);
                        
                        theta[idx_dilatation] = dudx + dvdy;
                    }
                }
            }
            else
            {
                for (int j = -d_num_subghosts_dilatation[1];
                     j < d_interior_dims[1] + d_num_subghosts_dilatation[1];
                     j++)
                {
                    for (int i = -d_num_subghosts_dilatation[0];
                         i < d_interior_dims[0] + d_num_subghosts_dilatation[0];
                         i++)
                    {
                        // Compute indices of current and neighboring cells.
                        const int idx = (i + d_num_subghosts_velocity[0]) +
                            (j + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0];
                        
                        const int idx_x_L = (i - 1 + d_num_subghosts_velocity[0]) +
                            (j + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0];
                        
                        const int idx_x_R = (i + 1 + d_num_subghosts_velocity[0]) +
                            (j + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0];
                        
                        const int idx_y_B = (i + d_num_subghosts_velocity[0]) +
                            (j - 1 + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0];
                        
                        const int idx_y_T = (i + d_num_subghosts_velocity[0]) +
                            (j + 1 + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0];
                        
                        const int idx_dilatation = (i + d_num_subghosts_dilatation[0]) +
                            (j + d_num_subghosts_dilatation[1])*d_subghostcell_dims_dilatation[0];
                        
                        double dudx, dvdy;
                        
                        if (i == -d_num_subghosts_velocity[0])
                        {
                            dudx = (u[idx_x_R] - u[idx])/(dx[0]);
                        }
                        else if (i == d_interior_dims[0] + d_num_subghosts_velocity[0] - 1)
                        {
                            dudx = (u[idx] - u[idx_x_L])/(dx[0]);
                        }
                        else
                        {
                            dudx = (u[idx_x_R] - u[idx_x_L])/(2*dx[0]);
                        }
                        
                        if (j == -d_num_subghosts_velocity[1])
                        {
                            dvdy = (v[idx_y_T] - v[idx])/(dx[1]);
                        }
                        else if (j == d_interior_dims[1] + d_num_subghosts_velocity[1] - 1)
                        {
                            dvdy = (v[idx] - v[idx_y_B])/(dx[1]);
                        }
                        else
                        {
                            dvdy = (v[idx_y_T] - v[idx_y_B])/(2*dx[1]);
                        }
                        
                        theta[idx_dilatation] = dudx + dvdy;
                    }
                }
            }
        }
        else if (d_dim == tbox::Dimension(3))
        {
            // Get the pointers to the cell data of velocity.
            double* u = d_data_velocity->getPointer(0);
            double* v = d_data_velocity->getPointer(1);
            double* w = d_data_velocity->getPointer(2);
            
            // Compute the dilatation field.
            if (d_num_subghosts_dilatation < d_num_subghosts_velocity)
            {
                for (int k = -d_num_subghosts_dilatation[2];
                     k < d_interior_dims[2] + d_num_subghosts_dilatation[2];
                     k++)
                {
                    for (int j = -d_num_subghosts_dilatation[1];
                         j < d_interior_dims[1] + d_num_subghosts_dilatation[1];
                         j++)
                    {
                        for (int i = -d_num_subghosts_dilatation[0];
                             i < d_interior_dims[0] + d_num_subghosts_dilatation[0];
                             i++)
                        {
                            // Compute indices of current and neighboring cells.
                            const int idx_x_L = (i - 1 + d_num_subghosts_velocity[0]) +
                                (j + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0] +
                                (k + d_num_subghosts_velocity[2])*d_subghostcell_dims_velocity[0]*
                                    d_subghostcell_dims_velocity[1];
                            
                            const int idx_x_R = (i + 1 + d_num_subghosts_velocity[0]) +
                                (j + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0] +
                                (k + d_num_subghosts_velocity[2])*d_subghostcell_dims_velocity[0]*
                                    d_subghostcell_dims_velocity[1];
                            
                            const int idx_y_B = (i + d_num_subghosts_velocity[0]) +
                                (j - 1 + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0] +
                                (k + d_num_subghosts_velocity[2])*d_subghostcell_dims_velocity[0]*
                                    d_subghostcell_dims_velocity[1];
                            
                            const int idx_y_T = (i + d_num_subghosts_velocity[0]) +
                                (j + 1 + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0] +
                                (k + d_num_subghosts_velocity[2])*d_subghostcell_dims_velocity[0]*
                                    d_subghostcell_dims_velocity[1];
                            
                            const int idx_z_B = (i + d_num_subghosts_velocity[0]) +
                                (j + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0] +
                                (k - 1 + d_num_subghosts_velocity[2])*d_subghostcell_dims_velocity[0]*
                                    d_subghostcell_dims_velocity[1];
                            
                            const int idx_z_F = (i + d_num_subghosts_velocity[0]) +
                                (j + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0] +
                                (k + 1 + d_num_subghosts_velocity[2])*d_subghostcell_dims_velocity[0]*
                                    d_subghostcell_dims_velocity[1];
                            
                            const int idx_dilatation = (i + d_num_subghosts_dilatation[0]) +
                                (j + d_num_subghosts_dilatation[1])*d_subghostcell_dims_dilatation[0] +
                                (k + d_num_subghosts_dilatation[2])*d_subghostcell_dims_dilatation[0]*
                                    d_subghostcell_dims_dilatation[1];
                            
                            double dudx = (u[idx_x_R] - u[idx_x_L])/(2*dx[0]);
                            double dvdy = (v[idx_y_T] - v[idx_y_B])/(2*dx[1]);
                            double dwdz = (w[idx_z_F] - w[idx_z_B])/(2*dx[2]);
                            
                            theta[idx_dilatation] = dudx + dvdy + dwdz;
                        }
                    }
                }
            }
            else
            {
                for (int k = -d_num_subghosts_dilatation[2];
                     k < d_interior_dims[2] + d_num_subghosts_dilatation[2];
                     k++)
                {
                    for (int j = -d_num_subghosts_dilatation[1];
                         j < d_interior_dims[1] + d_num_subghosts_dilatation[1];
                         j++)
                    {
                        for (int i = -d_num_subghosts_dilatation[0];
                             i < d_interior_dims[0] + d_num_subghosts_dilatation[0];
                             i++)
                        {
                            // Compute indices of current and neighboring cells.
                            const int idx = (i + d_num_subghosts_velocity[0]) +
                                (j + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0] +
                                (k + d_num_subghosts_velocity[2])*d_subghostcell_dims_velocity[0]*
                                    d_subghostcell_dims_velocity[1];
                            
                            const int idx_x_L = (i - 1 + d_num_subghosts_velocity[0]) +
                                (j + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0] +
                                (k + d_num_subghosts_velocity[2])*d_subghostcell_dims_velocity[0]*
                                    d_subghostcell_dims_velocity[1];
                            
                            const int idx_x_R = (i + 1 + d_num_subghosts_velocity[0]) +
                                (j + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0] +
                                (k + d_num_subghosts_velocity[2])*d_subghostcell_dims_velocity[0]*
                                    d_subghostcell_dims_velocity[1];
                            
                            const int idx_y_B = (i + d_num_subghosts_velocity[0]) +
                                (j - 1 + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0] +
                                (k + d_num_subghosts_velocity[2])*d_subghostcell_dims_velocity[0]*
                                    d_subghostcell_dims_velocity[1];
                            
                            const int idx_y_T = (i + d_num_subghosts_velocity[0]) +
                                (j + 1 + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0] +
                                (k + d_num_subghosts_velocity[2])*d_subghostcell_dims_velocity[0]*
                                    d_subghostcell_dims_velocity[1];
                            
                            const int idx_z_B = (i + d_num_subghosts_velocity[0]) +
                                (j + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0] +
                                (k - 1 + d_num_subghosts_velocity[2])*d_subghostcell_dims_velocity[0]*
                                    d_subghostcell_dims_velocity[1];
                            
                            const int idx_z_F = (i + d_num_subghosts_velocity[0]) +
                                (j + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0] +
                                (k + 1 + d_num_subghosts_velocity[2])*d_subghostcell_dims_velocity[0]*
                                    d_subghostcell_dims_velocity[1];
                            
                            const int idx_dilatation = (i + d_num_subghosts_dilatation[0]) +
                                (j + d_num_subghosts_dilatation[1])*d_subghostcell_dims_dilatation[0] +
                                (k + d_num_subghosts_dilatation[2])*d_subghostcell_dims_dilatation[0]*
                                    d_subghostcell_dims_dilatation[1];
                            
                            double dudx, dvdy, dwdz;
                            
                            if (i == -d_num_subghosts_velocity[0])
                            {
                                dudx = (u[idx_x_R] - u[idx])/(dx[0]);
                            }
                            else if (i == d_interior_dims[0] + d_num_subghosts_velocity[0] - 1)
                            {
                                dudx = (u[idx] - u[idx_x_L])/(dx[0]);
                            }
                            else
                            {
                                dudx = (u[idx_x_R] - u[idx_x_L])/(2*dx[0]);
                            }
                            
                            if (j == -d_num_subghosts_velocity[1])
                            {
                                dvdy = (v[idx_y_T] - v[idx])/(dx[1]);
                            }
                            else if (j == d_interior_dims[1] + d_num_subghosts_velocity[1] - 1)
                            {
                                dvdy = (v[idx] - v[idx_y_B])/(dx[1]);
                            }
                            else
                            {
                                dvdy = (v[idx_y_T] - v[idx_y_B])/(2*dx[1]);
                            }
                            
                            if (k == -d_num_subghosts_velocity[2])
                            {
                                dwdz = (w[idx_z_F] - w[idx])/(dx[2]);
                            }
                            else if (k == d_interior_dims[2] + d_num_subghosts_velocity[2] - 1)
                            {
                                dwdz = (w[idx] - w[idx_z_B])/(dx[2]);
                            }
                            else
                            {
                                dwdz = (w[idx_z_F] - w[idx_z_B])/(2*dx[2]);
                            }
                            
                            theta[idx_dilatation] = dudx + dvdy + dwdz;
                        }
                    }
                }
            }
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::computeGlobalCellDataDilatationWithVelocity()\n"
            << "Cell data of 'DILATATION' is not yet registered."
            << std::endl);
    }
}


/*
 * Compute the global cell data of vorticity with velocity in the registered patch.
 */
void
FlowModelSingleSpecies::computeGlobalCellDataVorticityWithVelocity()
{
    if (d_num_subghosts_vorticity > -hier::IntVector::getOne(d_dim))
    {
        // Get the grid spacing.
        const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
            BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                d_patch->getPatchGeometry()));
        
        const double* const dx = patch_geom->getDx();
        
        if (!d_data_velocity)
        {
            computeGlobalCellDataVelocity();
        }
        
        if (d_dim == tbox::Dimension(1))
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::computeGlobalCellDataVorticityWithVelocity()\n"
                << "Vorticity cannot be found for one-dimensional flow."
                << std::endl);
        }
        else if (d_dim == tbox::Dimension(2))
        {
            d_data_vorticity.reset(
                new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_vorticity));
            
            // Get the pointer to the cell data of vorticity.
            double* omega = d_data_vorticity->getPointer(0);
            
            // Get the pointers to the cell data of velocity.
            double* u = d_data_velocity->getPointer(0);
            double* v = d_data_velocity->getPointer(1);
            
            // Compute the vorticity field.
            if (d_num_subghosts_vorticity < d_num_subghosts_velocity)
            {
                for (int j = -d_num_subghosts_vorticity[1];
                     j < d_interior_dims[1] + d_num_subghosts_vorticity[1];
                     j++)
                {
                    for (int i = -d_num_subghosts_vorticity[0];
                         i < d_interior_dims[0] + d_num_subghosts_vorticity[0];
                         i++)
                    {
                        // Compute indices of current and neighboring cells.
                        const int idx_x_L = (i - 1 + d_num_subghosts_velocity[0]) +
                            (j + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0];
                        
                        const int idx_x_R = (i + 1 + d_num_subghosts_velocity[0]) +
                            (j + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0];
                        
                        const int idx_y_B = (i + d_num_subghosts_velocity[0]) +
                            (j - 1 + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0];
                        
                        const int idx_y_T = (i + d_num_subghosts_velocity[0]) +
                            (j + 1 + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0];
                        
                        const int idx_vorticity = (i + d_num_subghosts_vorticity[0]) +
                            (j + d_num_subghosts_vorticity[1])*d_subghostcell_dims_vorticity[0];
                        
                        double dvdx = (v[idx_x_R] - v[idx_x_L])/(2*dx[0]);
                        double dudy = (u[idx_y_T] - u[idx_y_B])/(2*dx[1]);
                        
                        omega[idx_vorticity] = dvdx - dudy;
                    }
                }
            }
            else
            {
                for (int j = -d_num_subghosts_vorticity[1];
                     j < d_interior_dims[1] + d_num_subghosts_vorticity[1];
                     j++)
                {
                    for (int i = -d_num_subghosts_vorticity[0];
                         i < d_interior_dims[0] + d_num_subghosts_vorticity[0];
                         i++)
                    {
                        // Compute indices of current and neighboring cells.
                        const int idx = (i + d_num_subghosts_velocity[0]) +
                            (j + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0];
                        
                        const int idx_x_L = (i - 1 + d_num_subghosts_velocity[0]) +
                            (j + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0];
                        
                        const int idx_x_R = (i + 1 + d_num_subghosts_velocity[0]) +
                            (j + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0];
                        
                        const int idx_y_B = (i + d_num_subghosts_velocity[0]) +
                            (j - 1 + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0];
                        
                        const int idx_y_T = (i + d_num_subghosts_velocity[0]) +
                            (j + 1 + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0];
                        
                        const int idx_vorticity = (i + d_num_subghosts_vorticity[0]) +
                            (j + d_num_subghosts_vorticity[1])*d_subghostcell_dims_vorticity[0];
                        
                        double dudy, dvdx;
                        
                        if (i == -d_num_subghosts_velocity[0])
                        {
                            dvdx = (v[idx_x_R] - v[idx])/(dx[0]);
                        }
                        else if (i == d_interior_dims[0] + d_num_subghosts_velocity[0] - 1)
                        {
                            dvdx = (v[idx] - v[idx_x_L])/(dx[0]);
                        }
                        else
                        {
                            dvdx = (v[idx_x_R] - v[idx_x_L])/(2*dx[0]);
                        }
                        
                        if (j == -d_num_subghosts_velocity[1])
                        {
                            dudy = (u[idx_y_T] - u[idx])/(dx[1]);
                        }
                        else if (j == d_interior_dims[1] + d_num_subghosts_velocity[1] - 1)
                        {
                            dudy = (u[idx] - u[idx_y_B])/(dx[1]);
                        }
                        else
                        {
                            dudy = (u[idx_y_T] - u[idx_y_B])/(2*dx[1]);
                        }
                        
                        omega[idx_vorticity] = dvdx - dudy;
                    }
                }
            }
        }
        else if (d_dim == tbox::Dimension(3))
        {
            d_data_vorticity.reset(
                new pdat::CellData<double>(d_interior_box, 3, d_num_subghosts_vorticity));
            
            // Get the pointers to the cell data of vorticity.
            double* omega_x = d_data_vorticity->getPointer(0);
            double* omega_y = d_data_vorticity->getPointer(1);
            double* omega_z = d_data_vorticity->getPointer(2);
            
            // Get the pointers to the cell data of velocity.
            double* u = d_data_velocity->getPointer(0);
            double* v = d_data_velocity->getPointer(1);
            double* w = d_data_velocity->getPointer(2);
            
            // Compute the vorticity field.
            if (d_num_subghosts_vorticity < d_num_subghosts_velocity)
            {
                for (int k = -d_num_subghosts_vorticity[2];
                     k < d_interior_dims[2] + d_num_subghosts_vorticity[2];
                     k++)
                {
                    for (int j = -d_num_subghosts_vorticity[1];
                         j < d_interior_dims[1] + d_num_subghosts_vorticity[1];
                         j++)
                    {
                        for (int i = -d_num_subghosts_vorticity[0];
                             i < d_interior_dims[0] + d_num_subghosts_vorticity[0];
                             i++)
                        {
                            // Compute indices of current and neighboring cells.
                            const int idx_x_L = (i - 1 + d_num_subghosts_velocity[0]) +
                                (j + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0] +
                                (k + d_num_subghosts_velocity[2])*d_subghostcell_dims_velocity[0]*
                                    d_subghostcell_dims_velocity[1];
                            
                            const int idx_x_R = (i + 1 + d_num_subghosts_velocity[0]) +
                                (j + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0] +
                                (k + d_num_subghosts_velocity[2])*d_subghostcell_dims_velocity[0]*
                                    d_subghostcell_dims_velocity[1];
                            
                            const int idx_y_B = (i + d_num_subghosts_velocity[0]) +
                                (j - 1 + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0] +
                                (k + d_num_subghosts_velocity[2])*d_subghostcell_dims_velocity[0]*
                                    d_subghostcell_dims_velocity[1];
                            
                            const int idx_y_T = (i + d_num_subghosts_velocity[0]) +
                                (j + 1 + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0] +
                                (k + d_num_subghosts_velocity[2])*d_subghostcell_dims_velocity[0]*
                                    d_subghostcell_dims_velocity[1];
                            
                            const int idx_z_B = (i + d_num_subghosts_velocity[0]) +
                                (j + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0] +
                                (k - 1 + d_num_subghosts_velocity[2])*d_subghostcell_dims_velocity[0]*
                                    d_subghostcell_dims_velocity[1];
                            
                            const int idx_z_F = (i + d_num_subghosts_velocity[0]) +
                                (j + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0] +
                                (k + 1 + d_num_subghosts_velocity[2])*d_subghostcell_dims_velocity[0]*
                                    d_subghostcell_dims_velocity[1];
                            
                            const int idx_vorticity = (i + d_num_subghosts_vorticity[0]) +
                                (j + d_num_subghosts_vorticity[1])*d_subghostcell_dims_vorticity[0] +
                                (k + d_num_subghosts_vorticity[2])*d_subghostcell_dims_vorticity[0]*
                                    d_subghostcell_dims_vorticity[1];
                            
                            double dvdx = (v[idx_x_R] - v[idx_x_L])/(2*dx[0]);
                            double dwdx = (w[idx_x_R] - w[idx_x_L])/(2*dx[0]);
                            double dudy = (u[idx_y_T] - u[idx_y_B])/(2*dx[1]);
                            double dwdy = (w[idx_y_T] - w[idx_y_B])/(2*dx[1]);
                            double dudz = (u[idx_z_F] - u[idx_z_B])/(2*dx[2]);
                            double dvdz = (v[idx_z_F] - v[idx_z_B])/(2*dx[2]);
                            
                            omega_x[idx_vorticity] = dwdy - dvdz;
                            omega_y[idx_vorticity] = dudz - dwdx;
                            omega_z[idx_vorticity] = dvdx - dudy;
                        }
                    }
                }
            }
            else
            {
                for (int k = -d_num_subghosts_vorticity[2]; k < d_interior_dims[2] + d_num_subghosts_vorticity[2]; k++)
                {
                    for (int j = -d_num_subghosts_vorticity[1]; j < d_interior_dims[1] + d_num_subghosts_vorticity[1]; j++)
                    {
                        for (int i = -d_num_subghosts_vorticity[0]; i < d_interior_dims[0] + d_num_subghosts_vorticity[0]; i++)
                        {
                            // Compute indices of current and neighboring cells.
                            const int idx = (i + d_num_subghosts_velocity[0]) +
                                (j + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0] +
                                (k + d_num_subghosts_velocity[2])*d_subghostcell_dims_velocity[0]*
                                    d_subghostcell_dims_velocity[1];
                            
                            const int idx_x_L = (i - 1 + d_num_subghosts_velocity[0]) +
                                (j + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0] +
                                (k + d_num_subghosts_velocity[2])*d_subghostcell_dims_velocity[0]*
                                    d_subghostcell_dims_velocity[1];
                            
                            const int idx_x_R = (i + 1 + d_num_subghosts_velocity[0]) +
                                (j + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0] +
                                (k + d_num_subghosts_velocity[2])*d_subghostcell_dims_velocity[0]*
                                    d_subghostcell_dims_velocity[1];
                            
                            const int idx_y_B = (i + d_num_subghosts_velocity[0]) +
                                (j - 1 + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0] +
                                (k + d_num_subghosts_velocity[2])*d_subghostcell_dims_velocity[0]*
                                    d_subghostcell_dims_velocity[1];
                            
                            const int idx_y_T = (i + d_num_subghosts_velocity[0]) +
                                (j + 1 + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0] +
                                (k + d_num_subghosts_velocity[2])*d_subghostcell_dims_velocity[0]*
                                    d_subghostcell_dims_velocity[1];
                            
                            const int idx_z_B = (i + d_num_subghosts_velocity[0]) +
                                (j + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0] +
                                (k - 1 + d_num_subghosts_velocity[2])*d_subghostcell_dims_velocity[0]*
                                    d_subghostcell_dims_velocity[1];
                            
                            const int idx_z_F = (i + d_num_subghosts_velocity[0]) +
                                (j + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0] +
                                (k + 1 + d_num_subghosts_velocity[2])*d_subghostcell_dims_velocity[0]*
                                    d_subghostcell_dims_velocity[1];
                            
                            const int idx_vorticity = (i + d_num_subghosts_vorticity[0]) +
                                (j + d_num_subghosts_vorticity[1])*d_subghostcell_dims_vorticity[0] +
                                (k + d_num_subghosts_vorticity[2])*d_subghostcell_dims_vorticity[0]*
                                    d_subghostcell_dims_vorticity[1];
                            
                            double dudy, dudz, dvdx, dvdz, dwdx, dwdy;
                            
                            if (i == -d_num_subghosts_velocity[0])
                            {
                                dvdx = (v[idx_x_R] - v[idx])/(dx[0]);
                                dwdx = (w[idx_x_R] - w[idx])/(dx[0]);
                            }
                            else if (i == d_interior_dims[0] + d_num_subghosts_velocity[0] - 1)
                            {
                                dvdx = (v[idx] - v[idx_x_L])/(dx[0]);
                                dwdx = (w[idx] - w[idx_x_L])/(dx[0]);
                            }
                            else
                            {
                                dvdx = (v[idx_x_R] - v[idx_x_L])/(2*dx[0]);
                                dwdx = (w[idx_x_R] - w[idx_x_L])/(2*dx[0]);
                            }
                            
                            if (j == -d_num_subghosts_velocity[1])
                            {
                                dudy = (u[idx_y_T] - u[idx])/(dx[1]);
                                dwdy = (w[idx_y_T] - w[idx])/(dx[1]);
                            }
                            else if (j == d_interior_dims[1] + d_num_subghosts_velocity[1] - 1)
                            {
                                dudy = (u[idx] - u[idx_y_B])/(dx[1]);
                                dwdy = (w[idx] - w[idx_y_B])/(dx[1]);
                            }
                            else
                            {
                                dudy = (u[idx_y_T] - u[idx_y_B])/(2*dx[1]);
                                dwdy = (w[idx_y_T] - w[idx_y_B])/(2*dx[1]);
                            }
                            
                            if (k == -d_num_subghosts_velocity[2])
                            {
                                dudz = (u[idx_z_F] - u[idx])/(dx[2]);
                                dvdz = (v[idx_z_F] - v[idx])/(dx[2]);
                            }
                            else if (k == d_interior_dims[2] + d_num_subghosts_velocity[2] - 1)
                            {
                                dudz = (u[idx] - u[idx_z_B])/(dx[2]);
                                dvdz = (v[idx] - v[idx_z_B])/(dx[2]);
                            }
                            else
                            {
                                dudz = (u[idx_z_F] - u[idx_z_B])/(2*dx[2]);
                                dvdz = (v[idx_z_F] - v[idx_z_B])/(2*dx[2]);
                            }
                            
                            omega_x[idx_vorticity] = dwdy - dvdz;
                            omega_y[idx_vorticity] = dudz - dwdx;
                            omega_z[idx_vorticity] = dvdx - dudy;
                        }
                    }
                }
            }
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::computeGlobalCellDataVorticityWithVelocity()\n"
            << "Cell data of 'VORTICITY' is not yet registered."
            << std::endl);
    }
}


/*
 * Compute the global cell data of enstrophy with vorticity in the registered patch.
 */
void
FlowModelSingleSpecies::computeGlobalCellDataEnstrophyWithVorticity()
{
    if (d_num_subghosts_enstrophy > -hier::IntVector::getOne(d_dim))
    {
        // Create the cell data of enstrophy.
        d_data_enstrophy.reset(
            new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_enstrophy));
        
        // Get the cell data of the vorticity.
        if (!d_data_vorticity) // If the pointer is null.
        {
            computeGlobalCellDataVorticityWithVelocity();
        }
        
        // Get the pointer to the cell data of enstrophy.
        double* Omega = d_data_enstrophy->getPointer(0);
        
        if (d_dim == tbox::Dimension(1))
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Enstrophy cannot be found for one-dimensional flow."
                << std::endl);
        }
        else if (d_dim == tbox::Dimension(2))
        {
            // Get the pointer to the cell data of vorticity.
            double* omega = d_data_vorticity->getPointer(0);
            
            // Compute the enstrophy field.
            for (int j = -d_num_subghosts_enstrophy[1];
                 j < d_interior_dims[1] + d_num_subghosts_enstrophy[1];
                 j++)
            {
                for (int i = -d_num_subghosts_enstrophy[0];
                     i < d_interior_dims[0] + d_num_subghosts_enstrophy[0];
                     i++)
                {
                    // Compute the linear indices.
                    const int idx_vorticity = (i + d_num_subghosts_vorticity[0]) +
                        (j + d_num_subghosts_vorticity[1])*d_subghostcell_dims_vorticity[0];
                    
                    const int idx_enstrophy = (i + d_num_subghosts_enstrophy[0]) +
                        (j + d_num_subghosts_enstrophy[1])*d_subghostcell_dims_enstrophy[0];
                    
                    Omega[idx_enstrophy] = omega[idx_vorticity]*omega[idx_vorticity];
                }
            }
        }
        else if (d_dim == tbox::Dimension(3))
        {
            // Get the pointers to the cell data of vorticity.
            double* omega_x = d_data_vorticity->getPointer(0);
            double* omega_y = d_data_vorticity->getPointer(1);
            double* omega_z = d_data_vorticity->getPointer(2);
            
            // Compute the enstrophy field.
            for (int k = -d_num_subghosts_enstrophy[2];
                 k < d_interior_dims[2] + d_num_subghosts_enstrophy[2];
                 k++)
            {
                for (int j = -d_num_subghosts_enstrophy[1];
                     j < d_interior_dims[1] + d_num_subghosts_enstrophy[1];
                     j++)
                {
                    for (int i = -d_num_subghosts_enstrophy[0];
                         i < d_interior_dims[0] + d_num_subghosts_enstrophy[0];
                         i++)
                    {
                        // Compute the linear indices.
                        const int idx_vorticity = (i + d_num_subghosts_vorticity[0]) +
                                (j + d_num_subghosts_vorticity[1])*d_subghostcell_dims_vorticity[0] +
                                (k + d_num_subghosts_vorticity[2])*d_subghostcell_dims_vorticity[0]*
                                    d_subghostcell_dims_vorticity[1];
                        
                        const int idx_enstrophy = (i + d_num_subghosts_enstrophy[0]) +
                                (j + d_num_subghosts_enstrophy[1])*d_subghostcell_dims_enstrophy[0] +
                                (k + d_num_subghosts_enstrophy[2])*d_subghostcell_dims_enstrophy[0]*
                                    d_subghostcell_dims_enstrophy[1];
                        
                        Omega[idx_enstrophy] = (omega_x[idx_vorticity]*omega_x[idx_vorticity] +
                            omega_y[idx_vorticity]*omega_y[idx_vorticity] +
                            omega_z[idx_vorticity]*omega_z[idx_vorticity]);
                    }
                }
            }
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::computeGlobalCellDataEnstrophyWithVorticity()\n"
            << "Cell data of 'ENSTROPHY' is not yet registered."
            << std::endl);
    }
}


/*
 * Compute the global cell data of convective flux with velocity and pressure in the registered patch.
 */
void
FlowModelSingleSpecies::computeGlobalCellDataConvectiveFluxWithVelocityAndPressure(DIRECTION::TYPE direction)
{
    if (direction == DIRECTION::X_DIRECTION)
    {
        if (d_num_subghosts_convective_flux_x > -hier::IntVector::getOne(d_dim))
        {
            // Create the cell data of convective flux in the x-direction.
            d_data_convective_flux_x.reset(
                new pdat::CellData<double>(d_interior_box, d_num_eqn, d_num_subghosts_convective_flux_x));
            
            // Get the pointers to the components of the convective flux in the x-direction.
            std::vector<double*> F_x;
            F_x.reserve(d_num_eqn);
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_x.push_back(d_data_convective_flux_x->getPointer(ei));
            }
            
            boost::shared_ptr<pdat::CellData<double> > data_momentum =
                getGlobalCellDataMomentum();
            
            boost::shared_ptr<pdat::CellData<double> > data_total_energy =
                getGlobalCellDataTotalEnergy();
            
            if (!d_data_velocity)
            {
                computeGlobalCellDataVelocity();
            }
            
            if (!d_data_pressure)
            {
                computeGlobalCellDataPressureWithInternalEnergy();
            }
            
            // Get the pointers to the cell data of total energy and pressure.
            double* E   = data_total_energy->getPointer(0);
            double* p   = d_data_pressure->getPointer(0);
            
            if (d_dim == tbox::Dimension(1))
            {
                // Get the pointer to the cell data of momentum.
                double* rho_u = data_momentum->getPointer(0);
                
                // Get the pointer to the cell data of velocity.
                double* u = d_data_velocity->getPointer(0);
                
                // Compute the convective flux in the x-direction.
                for (int i = -d_num_subghosts_convective_flux_x[0];
                     i < d_interior_dims[0] + d_num_subghosts_convective_flux_x[0];
                     i++)
                {
                    // Compute the linear indices.
                    const int idx = i + d_num_ghosts[0];
                    const int idx_pressure = i + d_num_subghosts_pressure[0];
                    const int idx_velocity = i + d_num_subghosts_velocity[0];
                    const int idx_convective_flux_x = i + d_num_subghosts_convective_flux_x[0];
                    
                    F_x[0][idx_convective_flux_x] = rho_u[idx];
                    F_x[1][idx_convective_flux_x] = u[idx_velocity]*rho_u[idx] + p[idx_pressure];
                    F_x[2][idx_convective_flux_x] = u[idx_velocity]*(E[idx] + p[idx_pressure]);
                }
            }
            else if (d_dim == tbox::Dimension(2))
            {
                // Get the pointers to the cell data of momentum.
                double* rho_u = data_momentum->getPointer(0);
                double* rho_v = data_momentum->getPointer(1);
                
                // Get the pointer to the cell data of velocity.
                double* u = d_data_velocity->getPointer(0);
                
                // Compute the convective flux in the x-direction.
                for (int j = -d_num_subghosts_convective_flux_x[1];
                     j < d_interior_dims[1] + d_num_subghosts_convective_flux_x[1];
                     j++)
                {
                    for (int i = -d_num_subghosts_convective_flux_x[0];
                         i < d_interior_dims[0] + d_num_subghosts_convective_flux_x[0];
                         i++)
                    {
                        const int idx = (i + d_num_ghosts[0]) +
                        (j + d_num_ghosts[1])*d_ghostcell_dims[0];
                        
                        const int idx_pressure = (i + d_num_subghosts_pressure[0]) +
                            (j + d_num_subghosts_pressure[1])*d_subghostcell_dims_pressure[0];
                        
                        const int idx_velocity = (i + d_num_subghosts_velocity[0]) +
                            (j + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0];
                        
                        const int idx_convective_flux_x = (i + d_num_subghosts_convective_flux_x[0]) +
                            (j + d_num_subghosts_convective_flux_x[1])*d_subghostcell_dims_convective_flux_x[0];
                        
                        F_x[0][idx_convective_flux_x] = rho_u[idx];
                        F_x[1][idx_convective_flux_x] = u[idx_velocity]*rho_u[idx] + p[idx_pressure];
                        F_x[2][idx_convective_flux_x] = u[idx_velocity]*rho_v[idx];
                        F_x[3][idx_convective_flux_x] = u[idx_velocity]*(E[idx] + p[idx_pressure]);
                    }
                }
            }
            else if (d_dim == tbox::Dimension(3))
            {
                // Get the pointers to the cell data of momentum.
                double* rho_u = data_momentum->getPointer(0);
                double* rho_v = data_momentum->getPointer(1);
                double* rho_w = data_momentum->getPointer(2);
                
                // Get the pointer to the cell data of velocity.
                double* u = d_data_velocity->getPointer(0);
                
                // Compute the convective flux in the x-direction.
                for (int k = -d_num_subghosts_convective_flux_x[2];
                     k < d_interior_dims[2] + d_num_subghosts_convective_flux_x[2];
                     k++)
                {
                    for (int j = -d_num_subghosts_convective_flux_x[1];
                         j < d_interior_dims[1] + d_num_subghosts_convective_flux_x[1];
                         j++)
                    {
                        for (int i = -d_num_subghosts_convective_flux_x[0];
                             i < d_interior_dims[0] + d_num_subghosts_convective_flux_x[0];
                             i++)
                        {
                            // Compute the linear indices.
                            const int idx = (i + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*d_ghostcell_dims[0] +
                                (k + d_num_ghosts[2])*d_ghostcell_dims[0]*d_ghostcell_dims[1];
                            
                            const int idx_pressure = (i + d_num_subghosts_pressure[0]) +
                                (j + d_num_subghosts_pressure[1])*d_subghostcell_dims_pressure[0] +
                                (k + d_num_subghosts_pressure[2])*d_subghostcell_dims_pressure[0]*
                                    d_subghostcell_dims_pressure[1];
                            
                            const int idx_velocity = (i + d_num_subghosts_velocity[0]) +
                                (j + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0] +
                                (k + d_num_subghosts_velocity[2])*d_subghostcell_dims_velocity[0]*
                                    d_subghostcell_dims_velocity[1];
                            
                            const int idx_convective_flux_x = (i + d_num_subghosts_convective_flux_x[0]) +
                                (j + d_num_subghosts_convective_flux_x[1])*d_subghostcell_dims_convective_flux_x[0] +
                                (k + d_num_subghosts_convective_flux_x[2])*d_subghostcell_dims_convective_flux_x[0]*
                                    d_subghostcell_dims_convective_flux_x[1];
                            
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
                << ": FlowModelSingleSpecies::computeGlobalCellDataConvectiveFluxWithVelocityAndPressure()\n"
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
            
            // Get the pointers to the components of the convective flux in the y-direction.
            std::vector<double*> F_y;
            F_y.reserve(d_num_eqn);
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_y.push_back(d_data_convective_flux_y->getPointer(ei));
            }
            
            boost::shared_ptr<pdat::CellData<double> > data_momentum =
                getGlobalCellDataMomentum();
            
            boost::shared_ptr<pdat::CellData<double> > data_total_energy =
                getGlobalCellDataTotalEnergy();
            
            if (!d_data_velocity)
            {
                computeGlobalCellDataVelocity();
            }
            
            if (!d_data_pressure)
            {
                computeGlobalCellDataPressureWithInternalEnergy();
            }
            
            // Get the pointers to the cell data of total energy and pressure.
            double* E   = data_total_energy->getPointer(0);
            double* p   = d_data_pressure->getPointer(0);
            
            if (d_dim == tbox::Dimension(1))
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::computeGlobalCellDataConvectiveFluxWithVelocityAndPressure()\n"
                    << "'CONVECTIVE_FLUX_Y' cannot be obtained for problem with dimension less than two."
                    << std::endl);
            }
            else if (d_dim == tbox::Dimension(2))
            {
                // Get the pointers to the cell data of momentum.
                double* rho_u = data_momentum->getPointer(0);
                double* rho_v = data_momentum->getPointer(1);
                
                // Get the pointer to the cell data of velocity.
                double* v = d_data_velocity->getPointer(1);
                
                // Compute the convective flux in the y-direction.
                for (int j = -d_num_subghosts_convective_flux_y[1];
                     j < d_interior_dims[1] + d_num_subghosts_convective_flux_y[1];
                     j++)
                {
                    for (int i = -d_num_subghosts_convective_flux_y[0];
                         i < d_interior_dims[0] + d_num_subghosts_convective_flux_y[0];
                         i++)
                    {
                        // Compute the linear indices.
                        const int idx = (i + d_num_ghosts[0]) +
                        (j + d_num_ghosts[1])*d_ghostcell_dims[0];
                        
                        const int idx_pressure = (i + d_num_subghosts_pressure[0]) +
                            (j + d_num_subghosts_pressure[1])*d_subghostcell_dims_pressure[0];
                        
                        const int idx_velocity = (i + d_num_subghosts_velocity[0]) +
                            (j + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0];
                        
                        const int idx_convective_flux_y = (i + d_num_subghosts_convective_flux_y[0]) +
                            (j + d_num_subghosts_convective_flux_y[1])*d_subghostcell_dims_convective_flux_y[0];
                        
                        F_y[0][idx_convective_flux_y] = rho_v[idx];
                        F_y[1][idx_convective_flux_y] = v[idx_velocity]*rho_u[idx];
                        F_y[2][idx_convective_flux_y] = v[idx_velocity]*rho_v[idx] + p[idx_pressure];
                        F_y[3][idx_convective_flux_y] = v[idx_velocity]*(E[idx] + p[idx_pressure]);
                    }
                }
            }
            else if (d_dim == tbox::Dimension(3))
            {
                // Get the pointers to the cell data of momentum.
                double* rho_u = data_momentum->getPointer(0);
                double* rho_v = data_momentum->getPointer(1);
                double* rho_w = data_momentum->getPointer(2);
                
                // Get the pointer to the cell data of velocity.
                double* v = d_data_velocity->getPointer(1);
                
                // Compute the convective flux in the y-direction.
                for (int k = -d_num_subghosts_convective_flux_y[2];
                     k < d_interior_dims[2] + d_num_subghosts_convective_flux_y[2];
                     k++)
                {
                    for (int j = -d_num_subghosts_convective_flux_y[1];
                         j < d_interior_dims[1] + d_num_subghosts_convective_flux_y[1];
                         j++)
                    {
                        for (int i = -d_num_subghosts_convective_flux_y[0];
                             i < d_interior_dims[0] + d_num_subghosts_convective_flux_y[0];
                             i++)
                        {
                            // Compute the linear indices.
                            const int idx = (i + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*d_ghostcell_dims[0] +
                                (k + d_num_ghosts[2])*d_ghostcell_dims[0]*d_ghostcell_dims[1];
                            
                            const int idx_pressure = (i + d_num_subghosts_pressure[0]) +
                                (j + d_num_subghosts_pressure[1])*d_subghostcell_dims_pressure[0] +
                                (k + d_num_subghosts_pressure[2])*d_subghostcell_dims_pressure[0]*
                                    d_subghostcell_dims_pressure[1];
                            
                            const int idx_velocity = (i + d_num_subghosts_velocity[0]) +
                                (j + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0] +
                                (k + d_num_subghosts_velocity[2])*d_subghostcell_dims_velocity[0]*
                                    d_subghostcell_dims_velocity[1];
                            
                            const int idx_convective_flux_y = (i + d_num_subghosts_convective_flux_y[0]) +
                                (j + d_num_subghosts_convective_flux_y[1])*d_subghostcell_dims_convective_flux_y[0] +
                                (k + d_num_subghosts_convective_flux_y[2])*d_subghostcell_dims_convective_flux_y[0]*
                                    d_subghostcell_dims_convective_flux_y[1];
                            
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
                << ": FlowModelSingleSpecies::computeGlobalCellDataConvectiveFluxWithVelocityAndPressure()\n"
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
            
            // Get the pointers to the components of the convective flux in the z-direction.
            std::vector<double*> F_z;
            F_z.reserve(d_num_eqn);
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_z.push_back(d_data_convective_flux_z->getPointer(ei));
            }
            
            boost::shared_ptr<pdat::CellData<double> > data_momentum =
                getGlobalCellDataMomentum();
            
            boost::shared_ptr<pdat::CellData<double> > data_total_energy =
                getGlobalCellDataTotalEnergy();
            
            if (!d_data_velocity)
            {
                computeGlobalCellDataVelocity();
            }
            
            if (!d_data_pressure)
            {
                computeGlobalCellDataPressureWithInternalEnergy();
            }
            
            // Get the pointers to the cell data of total energy and pressure.
            double* E   = data_total_energy->getPointer(0);
            double* p   = d_data_pressure->getPointer(0);
            
            if (d_dim == tbox::Dimension(1) || d_dim == tbox::Dimension(2))
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::computeGlobalCellDataConvectiveFluxWithVelocityAndPressure()\n"
                    << "'CONVECTIVE_FLUX_Z' cannot be obtained for problem with dimension less than three."
                    << std::endl);
            }
            else if (d_dim == tbox::Dimension(3))
            {
                // Get the pointers to the cell data of momentum.
                double* rho_u = data_momentum->getPointer(0);
                double* rho_v = data_momentum->getPointer(1);
                double* rho_w = data_momentum->getPointer(2);
                
                // Get the pointer to the cell data of velocity.
                double* w = d_data_velocity->getPointer(2);
                
                // Compute the convective flux in the z-direction.
                for (int k = -d_num_subghosts_convective_flux_z[2];
                     k < d_interior_dims[2] + d_num_subghosts_convective_flux_z[2];
                     k++)
                {
                    for (int j = -d_num_subghosts_convective_flux_z[1];
                         j < d_interior_dims[1] + d_num_subghosts_convective_flux_z[1];
                         j++)
                    {
                        for (int i = -d_num_subghosts_convective_flux_z[0];
                             i < d_interior_dims[0] + d_num_subghosts_convective_flux_z[0];
                             i++)
                        {
                            // Compute the linear indices.
                            const int idx = (i + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*d_ghostcell_dims[0] +
                                (k + d_num_ghosts[2])*d_ghostcell_dims[0]*d_ghostcell_dims[1];
                            
                            const int idx_pressure = (i + d_num_subghosts_pressure[0]) +
                                (j + d_num_subghosts_pressure[1])*d_subghostcell_dims_pressure[0] +
                                (k + d_num_subghosts_pressure[2])*d_subghostcell_dims_pressure[0]*
                                    d_subghostcell_dims_pressure[1];
                            
                            const int idx_velocity = (i + d_num_subghosts_velocity[0]) +
                                (j + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0] +
                                (k + d_num_subghosts_velocity[2])*d_subghostcell_dims_velocity[0]*
                                    d_subghostcell_dims_velocity[1];
                            
                            const int idx_convective_flux_z = (i + d_num_subghosts_convective_flux_z[0]) +
                                (j + d_num_subghosts_convective_flux_z[1])*d_subghostcell_dims_convective_flux_z[0] +
                                (k + d_num_subghosts_convective_flux_z[2])*d_subghostcell_dims_convective_flux_z[0]*
                                    d_subghostcell_dims_convective_flux_z[1];
                            
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
                << ": FlowModelSingleSpecies::computeGlobalCellDataConvectiveFluxWithVelocityAndPressure()\n"
                << "Cell data of 'CONVECTIVE_FLUX_Z' is not yet registered."
                << std::endl);
        }
    }
}


/*
 * Compute the global cell data of maximum wave speed with velocity and sound speed in the registered patch.
 */
void
FlowModelSingleSpecies::computeGlobalCellDataMaxWaveSpeedWithVelocityAndSoundSpeed(DIRECTION::TYPE direction)
{
    if (direction == DIRECTION::X_DIRECTION)
    {
        if (d_num_subghosts_max_wave_speed_x > -hier::IntVector::getOne(d_dim))
        {
            // Create the cell data of maximum wave speed in the x-direction.
            d_data_max_wave_speed_x.reset(
                new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_max_wave_speed_x));
            
            if (!d_data_velocity)
            {
                computeGlobalCellDataVelocity();
            }
            
            if (!d_data_sound_speed)
            {
                computeGlobalCellDataSoundSpeedWithPressure();
            }
            
            // Get the pointers to the cell data of maximum wave speed and velocity in x-direction, and sound speed.
            double* lambda_max_x = d_data_max_wave_speed_x->getPointer(0);
            double* u            = d_data_velocity->getPointer(0);
            double* c            = d_data_sound_speed->getPointer(0);
            
            if (d_dim == tbox::Dimension(1))
            {
                // Compute the maximum wave speed in the x-direction.
                for (int i = -d_num_subghosts_max_wave_speed_x[0];
                     i < d_interior_dims[0] + d_num_subghosts_max_wave_speed_x[0];
                     i++)
                {
                    // Compute the linear indices.
                    const int idx_sound_speed = i + d_num_subghosts_sound_speed[0];
                    const int idx_velocity = i + d_num_subghosts_velocity[0];
                    const int idx_max_wave_speed_x = i + d_num_subghosts_max_wave_speed_x[0];
                    
                    lambda_max_x[idx_max_wave_speed_x] = fabs(u[idx_velocity]) + c[idx_sound_speed];
                }
            }
            else if (d_dim == tbox::Dimension(2))
            {
                // Compute the maximum wave speed in the x-direction.
                for (int j = -d_num_subghosts_max_wave_speed_x[1];
                     j < d_interior_dims[1] + d_num_subghosts_max_wave_speed_x[1];
                     j++)
                {
                    for (int i = -d_num_subghosts_max_wave_speed_x[0];
                         i < d_interior_dims[0] + d_num_subghosts_max_wave_speed_x[0];
                         i++)
                    {
                        // Compute the linear indices.
                        const int idx_sound_speed = (i + d_num_subghosts_sound_speed[0]) +
                            (j + d_num_subghosts_sound_speed[1])*d_subghostcell_dims_sound_speed[0];
                        
                        const int idx_velocity = (i + d_num_subghosts_velocity[0]) +
                            (j + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0];
                        
                        const int idx_max_wave_speed_x = (i + d_num_subghosts_max_wave_speed_x[0]) +
                            (j + d_num_subghosts_max_wave_speed_x[1])*d_subghostcell_dims_max_wave_speed_x[0];
                        
                        lambda_max_x[idx_max_wave_speed_x] = fabs(u[idx_velocity]) + c[idx_sound_speed];
                    }
                }
            }
            else if (d_dim == tbox::Dimension(3))
            {
                // Compute the maximum wave speed in the x-direction.
                for (int k = -d_num_subghosts_max_wave_speed_x[2];
                     k < d_interior_dims[2] + d_num_subghosts_max_wave_speed_x[2];
                     k++)
                {
                    for (int j = -d_num_subghosts_max_wave_speed_x[1];
                         j < d_interior_dims[1] + d_num_subghosts_max_wave_speed_x[1];
                         j++)
                    {
                        for (int i = -d_num_subghosts_max_wave_speed_x[0];
                             i < d_interior_dims[0] + d_num_subghosts_max_wave_speed_x[0];
                             i++)
                        {
                            // Compute the linear indices.
                            const int idx_sound_speed = (i + d_num_subghosts_sound_speed[0]) +
                                (j + d_num_subghosts_sound_speed[1])*d_subghostcell_dims_sound_speed[0] +
                                (k + d_num_subghosts_sound_speed[2])*d_subghostcell_dims_sound_speed[0]*
                                    d_subghostcell_dims_sound_speed[1];
                            
                            const int idx_velocity = (i + d_num_subghosts_velocity[0]) +
                                (j + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0] +
                                (k + d_num_subghosts_velocity[2])*d_subghostcell_dims_velocity[0]*
                                    d_subghostcell_dims_velocity[1];
                            
                            const int idx_max_wave_speed_x = (i + d_num_subghosts_max_wave_speed_x[0]) +
                                (j + d_num_subghosts_max_wave_speed_x[1])*d_subghostcell_dims_max_wave_speed_x[0] +
                                (k + d_num_subghosts_max_wave_speed_x[2])*d_subghostcell_dims_max_wave_speed_x[0]*
                                    d_subghostcell_dims_max_wave_speed_x[1];
                            
                            lambda_max_x[idx_max_wave_speed_x] = fabs(u[idx_velocity]) + c[idx_sound_speed];
                        }
                    }
                }
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::computeGlobalCellDataMaxWaveSpeedWithVelocityAndSoundSpeed()\n"
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
            
            if (!d_data_velocity)
            {
                computeGlobalCellDataVelocity();
            }
            
            if (!d_data_sound_speed)
            {
                computeGlobalCellDataSoundSpeedWithPressure();
            }
            
            // Get the pointers to the cell data of maximum wave speed and velocity in y-direction, and sound speed.
            double* lambda_max_y = d_data_max_wave_speed_y->getPointer(0);
            double* v            = d_data_velocity->getPointer(1);
            double* c            = d_data_sound_speed->getPointer(0);
            
            if (d_dim == tbox::Dimension(1))
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::computeGlobalCellDataMaxWaveSpeedWithVelocityAndSoundSpeed()\n"
                    << "'MAX_WAVE_SPEED_Y' cannot be obtained for problem with dimension less than two."
                    << std::endl);
            }
            else if (d_dim == tbox::Dimension(2))
            {
                // Compute the maximum wave speed in the y-direction.
                for (int j = -d_num_subghosts_max_wave_speed_y[1];
                     j < d_interior_dims[1] + d_num_subghosts_max_wave_speed_y[1];
                     j++)
                {
                    for (int i = -d_num_subghosts_max_wave_speed_y[0];
                         i < d_interior_dims[0] + d_num_subghosts_max_wave_speed_y[0];
                         i++)
                    {
                        // Compute the linear indices.
                        const int idx_sound_speed = (i + d_num_subghosts_sound_speed[0]) +
                            (j + d_num_subghosts_sound_speed[1])*d_subghostcell_dims_sound_speed[0];
                        
                        const int idx_velocity = (i + d_num_subghosts_velocity[0]) +
                            (j + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0];
                        
                        const int idx_max_wave_speed_y = (i + d_num_subghosts_max_wave_speed_y[0]) +
                            (j + d_num_subghosts_max_wave_speed_y[1])*d_subghostcell_dims_max_wave_speed_y[0];
                        
                        lambda_max_y[idx_max_wave_speed_y] = fabs(v[idx_velocity]) + c[idx_sound_speed];
                    }
                }
            }
            else if (d_dim == tbox::Dimension(3))
            {
                // Compute the maximum wave speed in the y-direction.
                for (int k = -d_num_subghosts_max_wave_speed_y[2];
                     k < d_interior_dims[2] + d_num_subghosts_max_wave_speed_y[2];
                     k++)
                {
                    for (int j = -d_num_subghosts_max_wave_speed_y[1];
                         j < d_interior_dims[1] + d_num_subghosts_max_wave_speed_y[1];
                         j++)
                    {
                        for (int i = -d_num_subghosts_max_wave_speed_y[0];
                             i < d_interior_dims[0] + d_num_subghosts_max_wave_speed_y[0];
                             i++)
                        {
                            // Compute the linear indices.
                            const int idx_sound_speed = (i + d_num_subghosts_sound_speed[0]) +
                                (j + d_num_subghosts_sound_speed[1])*d_subghostcell_dims_sound_speed[0] +
                                (k + d_num_subghosts_sound_speed[2])*d_subghostcell_dims_sound_speed[0]*
                                    d_subghostcell_dims_sound_speed[1];
                            
                            const int idx_velocity = (i + d_num_subghosts_velocity[0]) +
                                (j + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0] +
                                (k + d_num_subghosts_velocity[2])*d_subghostcell_dims_velocity[0]*
                                    d_subghostcell_dims_velocity[1];
                            
                            const int idx_max_wave_speed_y = (i + d_num_subghosts_max_wave_speed_y[0]) +
                                (j + d_num_subghosts_max_wave_speed_y[1])*d_subghostcell_dims_max_wave_speed_y[0] +
                                (k + d_num_subghosts_max_wave_speed_y[2])*d_subghostcell_dims_max_wave_speed_y[0]*
                                    d_subghostcell_dims_max_wave_speed_y[1];
                            
                            lambda_max_y[idx_max_wave_speed_y] = fabs(v[idx_velocity]) + c[idx_sound_speed];
                        }
                    }
                }
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::computeGlobalCellDataMaxWaveSpeedWithVelocityAndSoundSpeed()\n"
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
            
            if (!d_data_velocity)
            {
                computeGlobalCellDataVelocity();
            }
            
            if (!d_data_sound_speed)
            {
                computeGlobalCellDataSoundSpeedWithPressure();
            }
            
            // Get the pointers to the cell data of maximum wave speed and velocity in z-direction, and sound speed.
            double* lambda_max_z = d_data_max_wave_speed_z->getPointer(0);
            double* w            = d_data_velocity->getPointer(2);
            double* c            = d_data_sound_speed->getPointer(0);
            
            if (d_dim == tbox::Dimension(1) || d_dim == tbox::Dimension(2))
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::computeGlobalCellDataMaxWaveSpeedWithVelocityAndSoundSpeed()\n"
                    << "'MAX_WAVE_SPEED_Z' cannot be obtained for problem with dimension less than three."
                    << std::endl);
            }
            else if (d_dim == tbox::Dimension(3))
            {
                // Compute the maximum wave speed in the z-direction.
                for (int k = -d_num_subghosts_max_wave_speed_z[2];
                     k < d_interior_dims[2] + d_num_subghosts_max_wave_speed_z[2];
                     k++)
                {
                    for (int j = -d_num_subghosts_max_wave_speed_z[1];
                         j < d_interior_dims[1] + d_num_subghosts_max_wave_speed_z[1];
                         j++)
                    {
                        for (int i = -d_num_subghosts_max_wave_speed_z[0];
                             i < d_interior_dims[0] + d_num_subghosts_max_wave_speed_z[0];
                             i++)
                        {
                            // Compute the linear indices.
                            const int idx_sound_speed = (i + d_num_subghosts_sound_speed[0]) +
                                (j + d_num_subghosts_sound_speed[1])*d_subghostcell_dims_sound_speed[0] +
                                (k + d_num_subghosts_sound_speed[2])*d_subghostcell_dims_sound_speed[0]*
                                    d_subghostcell_dims_sound_speed[1];
                            
                            const int idx_velocity = (i + d_num_subghosts_velocity[0]) +
                                (j + d_num_subghosts_velocity[1])*d_subghostcell_dims_velocity[0] +
                                (k + d_num_subghosts_velocity[2])*d_subghostcell_dims_velocity[0]*
                                    d_subghostcell_dims_velocity[1];
                            
                            const int idx_max_wave_speed_z = (i + d_num_subghosts_max_wave_speed_z[0]) +
                                (j + d_num_subghosts_max_wave_speed_z[1])*d_subghostcell_dims_max_wave_speed_z[0] +
                                (k + d_num_subghosts_max_wave_speed_z[2])*d_subghostcell_dims_max_wave_speed_z[0]*
                                    d_subghostcell_dims_max_wave_speed_z[1];
                            
                            lambda_max_z[idx_max_wave_speed_z] = fabs(w[idx_velocity]) + c[idx_sound_speed];
                        }
                    }
                }
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::computeGlobalCellDataMaxWaveSpeedWithVelocityAndSoundSpeed()\n"
                << "Cell data of 'MAX_WAVE_SPEED_Z' is not yet registered."
                << std::endl);
        }
    }
}
