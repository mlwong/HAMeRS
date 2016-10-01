#include "flow/flow_models/FlowModelSingleSpecies.hpp"

FlowModelSingleSpecies::FlowModelSingleSpecies(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
    const hier::IntVector& num_ghosts,
    const int& num_species,
    const boost::shared_ptr<tbox::Database>& flow_model_db):
        FlowModel(
            object_name,
            dim,
            grid_geometry,
            num_ghosts,
            num_species,
            2 + dim.getValue(),
            flow_model_db),
        d_num_subghosts_pressure(-hier::IntVector::getOne(d_dim)),
        d_num_subghosts_velocity(-hier::IntVector::getOne(d_dim)),
        d_num_subghosts_sound_speed(-hier::IntVector::getOne(d_dim)),
        d_num_subghosts_dilatation(-hier::IntVector::getOne(d_dim)),
        d_num_subghosts_vorticity(-hier::IntVector::getOne(d_dim)),
        d_num_subghosts_enstrophy(-hier::IntVector::getOne(d_dim)),
        d_num_subghosts_convective_flux_x(-hier::IntVector::getOne(d_dim)),
        d_num_subghosts_convective_flux_y(-hier::IntVector::getOne(d_dim)),
        d_num_subghosts_convective_flux_z(-hier::IntVector::getOne(d_dim)),
        d_num_subghosts_max_wave_speed_x(-hier::IntVector::getOne(d_dim)),
        d_num_subghosts_max_wave_speed_y(-hier::IntVector::getOne(d_dim)),
        d_num_subghosts_max_wave_speed_z(-hier::IntVector::getOne(d_dim)),
        d_subghost_box_pressure(hier::Box::getEmptyBox(dim)),
        d_subghost_box_velocity(hier::Box::getEmptyBox(dim)),
        d_subghost_box_sound_speed(hier::Box::getEmptyBox(dim)),
        d_subghost_box_dilatation(hier::Box::getEmptyBox(dim)),
        d_subghost_box_vorticity(hier::Box::getEmptyBox(dim)),
        d_subghost_box_enstrophy(hier::Box::getEmptyBox(dim)),
        d_subghost_box_convective_flux_x(hier::Box::getEmptyBox(dim)),
        d_subghost_box_convective_flux_y(hier::Box::getEmptyBox(dim)),
        d_subghost_box_convective_flux_z(hier::Box::getEmptyBox(dim)),
        d_subghost_box_max_wave_speed_x(hier::Box::getEmptyBox(dim)),
        d_subghost_box_max_wave_speed_y(hier::Box::getEmptyBox(dim)),
        d_subghost_box_max_wave_speed_z(hier::Box::getEmptyBox(dim)),
        d_subghostcell_dims_pressure(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_velocity(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_sound_speed(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_dilatation(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_vorticity(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_enstrophy(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_convective_flux_x(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_convective_flux_y(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_convective_flux_z(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_max_wave_speed_x(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_max_wave_speed_y(hier::IntVector::getZero(d_dim)),
        d_subghostcell_dims_max_wave_speed_z(hier::IntVector::getZero(d_dim))
{
    d_eqn_form.reserve(d_num_eqn);
    
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        d_eqn_form.push_back(CONSERVATIVE_EQN);
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
     * Initialize d_equation_of_state_manager and get the equation of state object.
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
    
    d_equation_of_state_manager.reset(new EquationOfStateManager(
        "d_equation_of_state_manager",
        d_dim,
        d_equation_of_state_str));
    
    d_equation_of_state =
        d_equation_of_state_manager->getEquationOfState();
    
    /*
     * Initialize d_equation_of_state_mixing_rules_manager and get the equation of state mixing rules object.
     */
    
    boost::shared_ptr<tbox::Database> species_db;
    
    if (flow_model_db->keyExists("Species"))
    {
        species_db = flow_model_db->getDatabase("Species");
    }
    else if (flow_model_db->keyExists("d_species_db"))
    {
        species_db = flow_model_db->getDatabase("d_species_db");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'Species'/'d_species_db' found in data for flow model"
            << std::endl);
    }
    
    d_equation_of_state_mixing_rules_manager.reset(new EquationOfStateMixingRulesManager(
        "d_equation_of_state_mixing_rules_manager",
        d_dim,
        d_num_species,
        NO_MODEL,
        species_db,
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
     * Initialize the Riemann solvers.
     */
    d_Riemann_solver_HLLC = boost::shared_ptr<RiemannSolverSingleSpeciesHLLC> (
        new RiemannSolverSingleSpeciesHLLC(
            d_object_name,
            d_dim,
            d_num_eqn,
            d_num_species,
            d_equation_of_state,
            d_equation_of_state_mixing_rules));
    
    d_Riemann_solver_HLLC_HLL = boost::shared_ptr<RiemannSolverSingleSpeciesHLLC_HLL> (
        new RiemannSolverSingleSpeciesHLLC_HLL(
            d_object_name,
            d_dim,
            d_num_eqn,
            d_num_species,
            d_equation_of_state,
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
    d_equation_of_state_manager->printClassData(os);
    
    os << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
    d_equation_of_state_mixing_rules_manager->printClassData(os);
}


/*
 * Register the conservative variables.
 */
void
FlowModelSingleSpecies::registerConservativeVariables(RungeKuttaLevelIntegrator* integrator)
{
    integrator->registerVariable(
        d_variable_density,
        d_num_ghosts,
        RungeKuttaLevelIntegrator::TIME_DEP,
        d_grid_geometry,
        "CONSERVATIVE_COARSEN",
        "CONSERVATIVE_LINEAR_REFINE");
    
    integrator->registerVariable(
        d_variable_momentum,
        d_num_ghosts,
        RungeKuttaLevelIntegrator::TIME_DEP,
        d_grid_geometry,
        "CONSERVATIVE_COARSEN",
        "CONSERVATIVE_LINEAR_REFINE");
    
    integrator->registerVariable(
        d_variable_total_energy,
        d_num_ghosts,
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
    
    if (num_subghosts_of_data.find("PRESSURE") != num_subghosts_of_data.end())
    {
        d_num_subghosts_pressure = num_subghosts_of_data.find("PRESSURE")->second;
    }
    
    if (num_subghosts_of_data.find("VELOCITY") != num_subghosts_of_data.end())
    {
        d_num_subghosts_velocity = num_subghosts_of_data.find("VELOCITY")->second;
    }
    
    if (num_subghosts_of_data.find("SOUND_SPEED") != num_subghosts_of_data.end())
    {
        d_num_subghosts_sound_speed = num_subghosts_of_data.find("SOUND_SPEED")->second;
        
        if (d_num_subghosts_pressure > -hier::IntVector::getOne(d_dim))
        {
            if (d_num_subghosts_sound_speed > d_num_subghosts_pressure)
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::registerDerivedCellVariable()\n"
                    << "Number of ghosts of 'SOUND_SPEED' exceeds"
                    << " number of ghosts of 'PRESSURE'."
                    << std::endl);
            }
        }
        else
        {
            d_num_subghosts_pressure = d_num_subghosts_sound_speed;
        }
    }
    
    if (num_subghosts_of_data.find("DILATATION") != num_subghosts_of_data.end())
    {
        d_num_subghosts_dilatation = num_subghosts_of_data.find("DILATATION")->second;
        
        if (d_num_subghosts_velocity > -hier::IntVector::getOne(d_dim))
        {
            if (d_num_subghosts_dilatation > d_num_subghosts_velocity)
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::registerDerivedCellVariable()\n"
                    << "Number of ghosts of 'DILATATION' exceeds"
                    << " number of ghosts of 'VELOCITY'."
                    << std::endl);
            }
        }
        else
        {
            d_num_subghosts_velocity = d_num_subghosts_dilatation;
        }
    }
    
    if (num_subghosts_of_data.find("VORTICITY") != num_subghosts_of_data.end())
    {
        d_num_subghosts_vorticity = num_subghosts_of_data.find("VORTICITY")->second;
        
        if (d_num_subghosts_velocity > -hier::IntVector::getOne(d_dim))
        {
            if (d_num_subghosts_vorticity > d_num_subghosts_velocity)
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::registerDerivedCellVariable()\n"
                    << "Number of ghosts of 'VORTICITY' exceeds"
                    << " number of ghosts of 'VELOCITY'."
                    << std::endl);
            }
        }
        else
        {
            d_num_subghosts_velocity = d_num_subghosts_vorticity;
        }
    }
    
    if (num_subghosts_of_data.find("ENSTROPHY") != num_subghosts_of_data.end())
    {
        d_num_subghosts_enstrophy = num_subghosts_of_data.find("ENSTROPHY")->second;
        
        if (d_num_subghosts_vorticity > -hier::IntVector::getOne(d_dim))
        {
            if (d_num_subghosts_enstrophy > d_num_subghosts_vorticity)
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::registerDerivedCellVariable()\n"
                    << "Number of ghosts of 'ENSTROPHY' exceeds"
                    << " number of ghosts of 'VORTICITY'."
                    << std::endl);
            }
        }
        else
        {
            d_num_subghosts_vorticity = d_num_subghosts_enstrophy;
            
            if (d_num_subghosts_velocity > -hier::IntVector::getOne(d_dim))
            {
                if (d_num_subghosts_enstrophy > d_num_subghosts_velocity)
                {
                    TBOX_ERROR(d_object_name
                        << ": FlowModelSingleSpecies::registerDerivedCellVariable()\n"
                        << "Number of ghosts of 'ENSTROPHY' exceeds"
                        << " number of ghosts of 'VELOCITY'."
                        << std::endl);
                }
            }
            else
            {
                d_num_subghosts_velocity = d_num_subghosts_enstrophy;
            }
        }
    }
    
    if (num_subghosts_of_data.find("CONVECTIVE_FLUX_X") != num_subghosts_of_data.end())
    {
        d_num_subghosts_convective_flux_x = num_subghosts_of_data.find("CONVECTIVE_FLUX_X")->second;
        
        if (d_num_subghosts_pressure > -hier::IntVector::getOne(d_dim))
        {
            if (d_num_subghosts_convective_flux_x > d_num_subghosts_pressure)
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::registerDerivedCellVariable()\n"
                    << "Number of ghosts of 'CONVECTIVE_FLUX_X' exceeds"
                    << " number of ghosts of 'PRESSURE'."
                    << std::endl);
            }
        }
        else
        {
            d_num_subghosts_pressure = d_num_subghosts_convective_flux_x;
        }
        
        if (d_num_subghosts_velocity > -hier::IntVector::getOne(d_dim))
        {
            if (d_num_subghosts_convective_flux_x > d_num_subghosts_velocity)
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::registerDerivedCellVariable()\n"
                    << "Number of ghosts of 'CONVECTIVE_FLUX_X' exceeds"
                    << " number of ghosts of 'VELOCITY'."
                    << std::endl);
            }
        }
        else
        {
            d_num_subghosts_velocity = d_num_subghosts_convective_flux_x;
        }
    }
    
    if (num_subghosts_of_data.find("CONVECTIVE_FLUX_Y") != num_subghosts_of_data.end())
    {
        if (d_dim >= tbox::Dimension(2))
        {
            d_num_subghosts_convective_flux_y = num_subghosts_of_data.find("CONVECTIVE_FLUX_Y")->second;
            
            if (d_num_subghosts_pressure > -hier::IntVector::getOne(d_dim))
            {
                if (d_num_subghosts_convective_flux_y > d_num_subghosts_pressure)
                {
                    TBOX_ERROR(d_object_name
                        << ": FlowModelSingleSpecies::registerDerivedCellVariable()\n"
                        << "Number of ghosts of 'CONVECTIVE_FLUX_Y' exceeds"
                        << " number of ghosts of 'PRESSURE'."
                        << std::endl);
                }
            }
            else
            {
                d_num_subghosts_pressure = d_num_subghosts_convective_flux_y;
            }
            
            if (d_num_subghosts_velocity > -hier::IntVector::getOne(d_dim))
            {
                if (d_num_subghosts_convective_flux_y > d_num_subghosts_velocity)
                {
                    TBOX_ERROR(d_object_name
                        << ": FlowModelSingleSpecies::registerDerivedCellVariable()\n"
                        << "Number of ghosts of 'CONVECTIVE_FLUX_Y' exceeds"
                        << " number of ghosts of 'VELOCITY'."
                        << std::endl);
                }
            }
            else
            {
                d_num_subghosts_velocity = d_num_subghosts_convective_flux_y;
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::registerDerivedCellVariable()\n"
                << "'CONVECTIVE_FLUX_Y' cannot be obtained for problem with dimension less than two."
                << std::endl);
        }
    }
    
    if (num_subghosts_of_data.find("CONVECTIVE_FLUX_Z") != num_subghosts_of_data.end())
    {
        if (d_dim == tbox::Dimension(3))
        {
            d_num_subghosts_convective_flux_z = num_subghosts_of_data.find("CONVECTIVE_FLUX_Z")->second;
            
            if (d_num_subghosts_pressure > -hier::IntVector::getOne(d_dim))
            {
                if (d_num_subghosts_convective_flux_z > d_num_subghosts_pressure)
                {
                    TBOX_ERROR(d_object_name
                        << ": FlowModelSingleSpecies::registerDerivedCellVariable()\n"
                        << "Number of ghosts of 'CONVECTIVE_FLUX_Z' exceeds"
                        << " number of ghosts of 'PRESSURE'."
                        << std::endl);
                }
            }
            else
            {
                d_num_subghosts_pressure = d_num_subghosts_convective_flux_z;
            }
            
            if (d_num_subghosts_velocity > -hier::IntVector::getOne(d_dim))
            {
                if (d_num_subghosts_convective_flux_z > d_num_subghosts_velocity)
                {
                    TBOX_ERROR(d_object_name
                        << ": FlowModelSingleSpecies::registerDerivedCellVariable()\n"
                        << "Number of ghosts of 'CONVECTIVE_FLUX_Z' exceeds"
                        << " number of ghosts of 'VELOCITY'."
                        << std::endl);
                }
            }
            else
            {
                d_num_subghosts_velocity = d_num_subghosts_convective_flux_z;
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::registerDerivedCellVariable()\n"
                << "'CONVECTIVE_FLUX_Z' cannot be obtained for problem with dimension less than three."
                << std::endl);
        }
    }
    
    if (num_subghosts_of_data.find("PRIMITIVE_VARIABLES") != num_subghosts_of_data.end())
    {
        hier::IntVector d_num_subghosts_primitive_variables =
            num_subghosts_of_data.find("PRIMITIVE_VARIABLES")->second;;
        
        if (d_num_subghosts_velocity > -hier::IntVector::getOne(d_dim))
        {
            if (d_num_subghosts_velocity != d_num_subghosts_primitive_variables)
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::registerDerivedCellVariable()\n"
                    << "Number of ghosts of 'VELOCITY' is not equal to"
                    << " number of ghosts of 'PRIMITIVE_VARIABLES'."
                    << std::endl);
            }
        }
        else
        {
            d_num_subghosts_velocity = d_num_subghosts_primitive_variables;
        }
        
        if (d_num_subghosts_pressure > -hier::IntVector::getOne(d_dim))
        {
            if (d_num_subghosts_pressure != d_num_subghosts_primitive_variables)
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::registerDerivedCellVariable()\n"
                    << "Number of ghosts of 'PRESSURE' is not equal to"
                    << " number of ghosts of 'PRIMITIVE_VARIABLES'."
                    << std::endl);
            }
        }
        else
        {
            d_num_subghosts_pressure = d_num_subghosts_primitive_variables;
        }
    }
    
    if (num_subghosts_of_data.find("MAX_WAVE_SPEED_X") != num_subghosts_of_data.end())
    {
        d_num_subghosts_max_wave_speed_x = num_subghosts_of_data.find("MAX_WAVE_SPEED_X")->second;
        
        if (d_num_subghosts_sound_speed > -hier::IntVector::getOne(d_dim))
        {
            if (d_num_subghosts_max_wave_speed_x > d_num_subghosts_sound_speed)
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::registerDerivedCellVariable()\n"
                    << "Number of ghosts of 'MAX_WAVE_SPEED_X' exceeds"
                    << " number of ghosts of 'SOUND_SPEED'."
                    << std::endl);
            }
        }
        else
        {
            d_num_subghosts_sound_speed = d_num_subghosts_max_wave_speed_x;
            
            if (d_num_subghosts_pressure > -hier::IntVector::getOne(d_dim))
            {
                if (d_num_subghosts_max_wave_speed_x > d_num_subghosts_pressure)
                {
                    TBOX_ERROR(d_object_name
                        << ": FlowModelSingleSpecies::registerDerivedCellVariable()\n"
                        << "Number of ghosts of 'MAX_WAVE_SPEED_X' exceeds"
                        << " number of ghosts of 'PRESSURE'."
                        << std::endl);
                }
            }
            else
            {
                d_num_subghosts_pressure = d_num_subghosts_max_wave_speed_x;
            }
        }
        
        if (d_num_subghosts_velocity > -hier::IntVector::getOne(d_dim))
        {
            if (d_num_subghosts_max_wave_speed_x > d_num_subghosts_velocity)
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::registerDerivedCellVariable()\n"
                    << "Number of ghosts of 'MAX_WAVE_SPEED_X' exceeds"
                    << " number of ghosts of 'VELOCITY'."
                    << std::endl);
            }
        }
        else
        {
            d_num_subghosts_velocity = d_num_subghosts_max_wave_speed_x;
        }
    }
    
    if (num_subghosts_of_data.find("MAX_WAVE_SPEED_Y") != num_subghosts_of_data.end())
    {
        d_num_subghosts_max_wave_speed_y = num_subghosts_of_data.find("MAX_WAVE_SPEED_Y")->second;
        
        if (d_num_subghosts_sound_speed > -hier::IntVector::getOne(d_dim))
        {
            if (d_num_subghosts_max_wave_speed_y > d_num_subghosts_sound_speed)
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::registerDerivedCellVariable()\n"
                    << "Number of ghosts of 'MAX_WAVE_SPEED_Y' exceeds"
                    << " number of ghosts of 'SOUND_SPEED'."
                    << std::endl);
            }
        }
        else
        {
            d_num_subghosts_sound_speed = d_num_subghosts_max_wave_speed_y;
            
            if (d_num_subghosts_pressure > -hier::IntVector::getOne(d_dim))
            {
                if (d_num_subghosts_max_wave_speed_y > d_num_subghosts_pressure)
                {
                    TBOX_ERROR(d_object_name
                        << ": FlowModelSingleSpecies::registerDerivedCellVariable()\n"
                        << "Number of ghosts of 'MAX_WAVE_SPEED_Y' exceeds"
                        << " number of ghosts of 'PRESSURE'."
                        << std::endl);
                }
            }
            else
            {
                d_num_subghosts_pressure = d_num_subghosts_max_wave_speed_y;
            }
        }
        
        if (d_num_subghosts_velocity > -hier::IntVector::getOne(d_dim))
        {
            if (d_num_subghosts_max_wave_speed_y > d_num_subghosts_velocity)
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::registerDerivedCellVariable()\n"
                    << "Number of ghosts of 'MAX_WAVE_SPEED_Y' exceeds"
                    << " number of ghosts of 'VELOCITY'."
                    << std::endl);
            }
        }
        else
        {
            d_num_subghosts_velocity = d_num_subghosts_max_wave_speed_y;
        }
    }
    
    if (num_subghosts_of_data.find("MAX_WAVE_SPEED_Z") != num_subghosts_of_data.end())
    {
        d_num_subghosts_max_wave_speed_z = num_subghosts_of_data.find("MAX_WAVE_SPEED_Z")->second;
        
        if (d_num_subghosts_sound_speed > -hier::IntVector::getOne(d_dim))
        {
            if (d_num_subghosts_max_wave_speed_z > d_num_subghosts_sound_speed)
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::registerDerivedCellVariable()\n"
                    << "Number of ghosts of 'MAX_WAVE_SPEED_Z' exceeds"
                    << " number of ghosts of 'SOUND_SPEED'."
                    << std::endl);
            }
        }
        else
        {
            d_num_subghosts_sound_speed = d_num_subghosts_max_wave_speed_z;
            
            if (d_num_subghosts_pressure > -hier::IntVector::getOne(d_dim))
            {
                if (d_num_subghosts_max_wave_speed_z > d_num_subghosts_pressure)
                {
                    TBOX_ERROR(d_object_name
                        << ": FlowModelSingleSpecies::registerDerivedCellVariable()\n"
                        << "Number of ghosts of 'MAX_WAVE_SPEED_Z' exceeds"
                        << " number of ghosts of 'PRESSURE'."
                        << std::endl);
                }
            }
            else
            {
                d_num_subghosts_pressure = d_num_subghosts_max_wave_speed_z;
            }
        }
        
        if (d_num_subghosts_velocity > -hier::IntVector::getOne(d_dim))
        {
            if (d_num_subghosts_max_wave_speed_z > d_num_subghosts_velocity)
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::registerDerivedCellVariable()\n"
                    << "Number of ghosts of 'MAX_WAVE_SPEED_Z' exceeds"
                    << " number of ghosts of 'VELOCITY'."
                    << std::endl);
            }
        }
        else
        {
            d_num_subghosts_velocity = d_num_subghosts_max_wave_speed_z;
        }
    }
}


/*
 * Register the required variables for the computation of projection matrix
 * of conservative variables and its inverse at faces in the registered patch.
 */
void
FlowModelSingleSpecies::registerFaceProjectionMatricesOfConservativeVariables(
    const hier::IntVector& num_subghosts,
    const AVERAGING& averaging)
{
    // Check whether a patch is already registered.
    if (!d_patch)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::registerFaceProjectionMatricesOfConservativeVariables()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    NULL_USE(num_subghosts);
    
    d_proj_mat_conservative_var_averaging = averaging;
    
    d_proj_mat_conservative_var_registered = true;
}


/*
 * Register the required variables for the computation of projection matrix
 * of primitive variables and its inverse at faces in the registered patch.
 */
void
FlowModelSingleSpecies::registerFaceProjectionMatricesOfPrimitiveVariables(
    const hier::IntVector& num_subghosts,
    const AVERAGING& averaging)
{
    // Check whether a patch is already registered.
    if (!d_patch)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::registerFaceProjectionMatricesOfPrimitiveVariables()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    d_proj_mat_primitive_var_averaging = averaging;
    
    if (d_num_subghosts_sound_speed > -hier::IntVector::getOne(d_dim))
    {
        if (num_subghosts > d_num_subghosts_sound_speed)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::registerFaceProjectionMatrices()\n"
                << "Number of ghosts of projection matrices exceeds"
                << " number of ghosts of 'SOUND_SPEED'."
                << std::endl);
        }
    }
    else
    {
        d_num_subghosts_sound_speed = num_subghosts;
        
        if (d_num_subghosts_pressure > -hier::IntVector::getOne(d_dim))
        {
            if (num_subghosts > d_num_subghosts_pressure)
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::registerFaceProjectionMatrices()\n"
                    << "Number of ghosts of projection matrices exceeds"
                    << " number of ghosts of 'PRESSURE'."
                    << std::endl);
            }
        }
        else
        {
            d_num_subghosts_pressure = d_num_subghosts_sound_speed;
        }
    }
    
    d_proj_mat_primitive_var_registered = true;
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
    
    d_num_subghosts_pressure          = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_velocity          = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_sound_speed       = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_dilatation        = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_vorticity         = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_enstrophy         = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_convective_flux_x = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_convective_flux_y = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_convective_flux_z = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_max_wave_speed_x  = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_max_wave_speed_y  = -hier::IntVector::getOne(d_dim);
    d_num_subghosts_max_wave_speed_z  = -hier::IntVector::getOne(d_dim);
    
    d_interior_box                   = hier::Box::getEmptyBox(d_dim);
    d_ghost_box                      = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_pressure          = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_velocity          = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_sound_speed       = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_dilatation        = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_vorticity         = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_enstrophy         = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_convective_flux_x = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_convective_flux_y = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_convective_flux_z = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_max_wave_speed_x  = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_max_wave_speed_y  = hier::Box::getEmptyBox(d_dim);
    d_subghost_box_max_wave_speed_z  = hier::Box::getEmptyBox(d_dim);
    
    d_interior_dims                       = hier::IntVector::getZero(d_dim);
    d_ghostcell_dims                      = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_pressure          = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_velocity          = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_sound_speed       = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_dilatation        = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_vorticity         = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_enstrophy         = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_convective_flux_x = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_convective_flux_y = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_convective_flux_z = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_max_wave_speed_x  = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_max_wave_speed_y  = hier::IntVector::getZero(d_dim);
    d_subghostcell_dims_max_wave_speed_z  = hier::IntVector::getZero(d_dim);
    
    d_data_pressure.reset();
    d_data_velocity.reset();
    d_data_sound_speed.reset();
    d_data_dilatation.reset();
    d_data_vorticity.reset();
    d_data_enstrophy.reset();
    d_data_convective_flux_x.reset();
    d_data_convective_flux_y.reset();
    d_data_convective_flux_z.reset();
    d_data_max_wave_speed_x.reset();
    d_data_max_wave_speed_y.reset();
    d_data_max_wave_speed_z.reset();
    
    d_proj_mat_conservative_var_registered = false;
    d_proj_mat_primitive_var_registered    = false;
    
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
    
    // Compute the pressure cell data.
    if (d_num_subghosts_pressure > -hier::IntVector::getOne(d_dim))
    {
        if (!d_data_pressure)
        {
            computeGlobalCellDataPressure();
        }
    }
    
    // Compute the velocity cell data.
    if (d_num_subghosts_velocity > -hier::IntVector::getOne(d_dim))
    {
        if (!d_data_velocity)
        {
            computeGlobalCellDataVelocity();
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
            computeGlobalCellDataEnstrophyWithVelocityAndVorticity();
        }
    }
    
    // Compute the x-direction convective flux cell data.
    if (d_num_subghosts_convective_flux_x > -hier::IntVector::getOne(d_dim))
    {
        if (!d_data_convective_flux_x)
        {
            computeGlobalCellDataConvectiveFluxWithPressureAndVelocity(X_DIRECTION);
        }
    }
    
    // Compute the y-direction convective flux cell data.
    if (d_num_subghosts_convective_flux_y > -hier::IntVector::getOne(d_dim))
    {
        if (!d_data_convective_flux_y)
        {
            computeGlobalCellDataConvectiveFluxWithPressureAndVelocity(Y_DIRECTION);
        }
    }
    
    // Compute the z-direction convective flux cell data.
    if (d_num_subghosts_convective_flux_z > -hier::IntVector::getOne(d_dim))
    {
        if (!d_data_convective_flux_z)
        {
            computeGlobalCellDataConvectiveFluxWithPressureAndVelocity(Z_DIRECTION);
        }
    }
    
    // Compute the x-direction maximum wave speed cell data.
    if (d_num_subghosts_max_wave_speed_x > -hier::IntVector::getOne(d_dim))
    {
        if (!d_data_max_wave_speed_x)
        {
            computeGlobalCellDataMaxWaveSpeedWithPressureVelocityAndSoundSpeed(X_DIRECTION);
        }
    }
    
    // Compute the y-direction maximum wave speed cell data.
    if (d_num_subghosts_max_wave_speed_y > -hier::IntVector::getOne(d_dim))
    {
        if (!d_data_max_wave_speed_y)
        {
            computeGlobalCellDataMaxWaveSpeedWithPressureVelocityAndSoundSpeed(Y_DIRECTION);
        }
    }
    
    // Compute the z-direction maximum wave speed cell data.
    if (d_num_subghosts_max_wave_speed_z > -hier::IntVector::getOne(d_dim))
    {
        if (!d_data_max_wave_speed_z)
        {
            computeGlobalCellDataMaxWaveSpeedWithPressureVelocityAndSoundSpeed(Z_DIRECTION);
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
        computeGlobalCellDataPressure();
    }
    global_cell_data.push_back(d_data_pressure);
    
    return global_cell_data;
}


/*
 * Compute the local face data of projection matrix of conservative variables in the
 * registered patch.
 */
void
FlowModelSingleSpecies::computeLocalFaceProjectionMatrixOfConservativeVariables(
    boost::multi_array<double, 2>& projection_matrix,
    const hier::Index& cell_index_minus,
    const hier::Index& cell_index_plus,
    const DIRECTION& direction)
{
    NULL_USE(projection_matrix);
    NULL_USE(cell_index_minus);
    NULL_USE(cell_index_plus);
    NULL_USE(direction);
    
    TBOX_ERROR(d_object_name
        << ": FlowModelSingleSpecies::computeLocalFaceProjectionMatrixOfConservativeVariables()\n"
        << "Method computeLocalFaceProjectionMatrixOfConservativeVariables() is not yet implemented."
        << std::endl);
}


/*
 * Compute the local face data of inverse of projection matrix of conservative variables
 * in the registered patch.
 */
void
FlowModelSingleSpecies::computeLocalFaceProjectionMatrixInverseOfConservativeVariables(
    boost::multi_array<double, 2>& projection_matrix_inv,
    const hier::Index& cell_index_minus,
    const hier::Index& cell_index_plus,
    const DIRECTION& direction)
{
    NULL_USE(projection_matrix_inv);
    NULL_USE(cell_index_minus);
    NULL_USE(cell_index_plus);
    NULL_USE(direction);
    
    TBOX_ERROR(d_object_name
        << ": FlowModelSingleSpecies::computeLocalFaceProjectionMatrixInverseOfConservativeVariables()\n"
        << "Method computeLocalFaceProjectionMatrixInverseOfConservativeVariables() is not yet implemented."
        << std::endl);
}


/*
 * Compute the local face data of projection matrix of primitive variables in the
 * registered patch.
 */
void
FlowModelSingleSpecies::computeLocalFaceProjectionMatrixOfPrimitiveVariables(
    boost::multi_array<double, 2>& projection_matrix,
    const hier::Index& cell_index_minus,
    const hier::Index& cell_index_plus,
    const DIRECTION& direction)
{
    if (!d_proj_mat_primitive_var_registered)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::computeLocalFaceProjectionMatrixOfPrimitiveVariables()\n"
            << "Projection matrices is not yet registered."
            << std::endl);
    }
    
    projection_matrix.resize(boost::extents[d_num_eqn][d_num_eqn]);
    
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
    
    // Compute the projection matrix.
    if (d_dim == tbox::Dimension(1))
    {
        // Compute the linear indices.
        const int idx_minus = cell_index_minus[0] + d_num_ghosts[0];
        const int idx_plus = cell_index_plus[0] + d_num_ghosts[0];
        const int idx_sound_speed_minus = cell_index_minus[0] + d_num_subghosts_sound_speed[0];
        const int idx_sound_speed_plus = cell_index_plus[0] + d_num_subghosts_sound_speed[0];
        
        // Compute the average values.
        double rho_average, c_average;
        switch (d_proj_mat_primitive_var_averaging)
        {
            case SIMPLE_AVG:
            {
                rho_average = 0.5*(rho[idx_minus] + rho[idx_plus]);
                c_average = 0.5*(c[idx_sound_speed_minus] + c[idx_sound_speed_plus]);
                
                break;
            }
            case ROE_AVG:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::computeLocalFaceProjectionMatrixOfPrimitiveVariables()\n"
                    << "Roe averaging is not yet implemented."
                    << std::endl);
                
                rho_average = 0.0;
                c_average   = 0.0;
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::computeLocalFaceProjectionMatrixOfPrimitiveVariables()\n"
                    << "Unknown d_proj_mat_primitive_var_averaging given."
                    << std::endl);
                
                rho_average = 0.0;
                c_average   = 0.0;
            }
        }
        
        switch (direction)
        {
            case X_DIRECTION:
            {
                projection_matrix[0][0] = 0.0;
                projection_matrix[0][1] = -0.5*rho_average*c_average;
                projection_matrix[0][2] = 0.5;
                projection_matrix[1][0] = 1.0;
                projection_matrix[1][1] = 0.0;
                projection_matrix[1][2] = -1.0/(c_average*c_average);
                projection_matrix[2][0] = 0.0;
                projection_matrix[2][1] = 0.5*rho_average*c_average;
                projection_matrix[2][2] = 0.5;
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::computeLocalFaceProjectionMatrixOfPrimitiveVariables()\n"
                << "There is only x-direction for one-dimensional problem."
                << std::endl);
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int idx_minus = (cell_index_minus[0] + d_num_ghosts[0]) +
            (cell_index_minus[1] + d_num_ghosts[1])*d_ghostcell_dims[0];
        
        const int idx_plus = (cell_index_plus[0] + d_num_ghosts[0]) +
            (cell_index_plus[1] + d_num_ghosts[1])*d_ghostcell_dims[0];
        
        const int idx_sound_speed_minus = (cell_index_minus[0] + d_num_subghosts_sound_speed[0]) +
            (cell_index_minus[1] + d_num_subghosts_sound_speed[1])*d_subghostcell_dims_sound_speed[0];
        
        const int idx_sound_speed_plus = (cell_index_plus[0] + d_num_subghosts_sound_speed[0]) +
            (cell_index_plus[1] + d_num_subghosts_sound_speed[1])*d_subghostcell_dims_sound_speed[0];
        
        // Compute the average values.
        double rho_average, c_average;
        switch (d_proj_mat_primitive_var_averaging)
        {
            case SIMPLE_AVG:
            {
                rho_average = 0.5*(rho[idx_minus] + rho[idx_plus]);
                c_average = 0.5*(c[idx_sound_speed_minus] + c[idx_sound_speed_plus]);
                
                break;
            }
            case ROE_AVG:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::computeLocalFaceProjectionMatrixOfPrimitiveVariables()\n"
                    << "Roe averaging is not yet implemented."
                    << std::endl);
                
                rho_average = 0.0;
                c_average   = 0.0;
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::computeLocalFaceProjectionMatrixOfPrimitiveVariables()\n"
                    << "Unknown d_proj_mat_primitive_var_averaging given."
                    << std::endl);
                
                rho_average = 0.0;
                c_average   = 0.0;
            }
        }
        
        switch (direction)
        {
            case X_DIRECTION:
            {
                projection_matrix[0][0] = 0.0;
                projection_matrix[0][1] = -0.5*rho_average*c_average;
                projection_matrix[0][2] = 0.0;
                projection_matrix[0][3] = 0.5;
                projection_matrix[1][0] = 1.0;
                projection_matrix[1][1] = 0.0;
                projection_matrix[1][2] = 0.0;
                projection_matrix[1][3] = -1.0/(c_average*c_average);
                projection_matrix[2][0] = 0.0;
                projection_matrix[2][1] = 0.0;
                projection_matrix[2][2] = 1.0;
                projection_matrix[2][3] = 0.0;
                projection_matrix[3][0] = 0.0;
                projection_matrix[3][1] = 0.5*rho_average*c_average;
                projection_matrix[3][2] = 0.0;
                projection_matrix[3][3] = 0.5;
                
                break;
            }
            case Y_DIRECTION:
            {
                projection_matrix[0][0] = 0.0;
                projection_matrix[0][1] = 0.0;
                projection_matrix[0][2] = -0.5*rho_average*c_average;
                projection_matrix[0][3] = 0.5;
                projection_matrix[1][0] = 1.0;
                projection_matrix[1][1] = 0.0;
                projection_matrix[1][2] = 0.0;
                projection_matrix[1][3] = -1.0/(c_average*c_average);
                projection_matrix[2][0] = 0.0;
                projection_matrix[2][1] = 1.0;
                projection_matrix[2][2] = 0.0;
                projection_matrix[2][3] = 0.0;
                projection_matrix[3][0] = 0.0;
                projection_matrix[3][1] = 0.0;
                projection_matrix[3][2] = 0.5*rho_average*c_average;
                projection_matrix[3][3] = 0.5;
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::computeLocalFaceProjectionMatrixOfPrimitiveVariables()\n"
                << "There are only x-direction and y-direction for two-dimensional problem."
                << std::endl);
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int idx_minus = (cell_index_minus[0] + d_num_ghosts[0]) +
            (cell_index_minus[1] + d_num_ghosts[1])*d_ghostcell_dims[0] +
            (cell_index_minus[2] + d_num_ghosts[2])*d_ghostcell_dims[0]*
                d_ghostcell_dims[1];
        
        const int idx_plus = (cell_index_plus[0] + d_num_ghosts[0]) +
            (cell_index_plus[1] + d_num_ghosts[1])*d_ghostcell_dims[0] +
            (cell_index_plus[2] + d_num_ghosts[2])*d_ghostcell_dims[0]*
                d_ghostcell_dims[1];
        
        const int idx_sound_speed_minus = (cell_index_minus[0] + d_num_subghosts_sound_speed[0]) +
            (cell_index_minus[1] + d_num_subghosts_sound_speed[1])*d_subghostcell_dims_sound_speed[0] +
            (cell_index_minus[2] + d_num_subghosts_sound_speed[2])*d_subghostcell_dims_sound_speed[0]*
                d_subghostcell_dims_sound_speed[1];
        
        const int idx_sound_speed_plus = (cell_index_plus[0] + d_num_subghosts_sound_speed[0]) +
            (cell_index_plus[1] + d_num_subghosts_sound_speed[1])*d_subghostcell_dims_sound_speed[0] +
            (cell_index_plus[2] + d_num_subghosts_sound_speed[2])*d_subghostcell_dims_sound_speed[0]*
                d_subghostcell_dims_sound_speed[1];
        
        // Compute the average values.
        double rho_average, c_average;
        switch (d_proj_mat_primitive_var_averaging)
        {
            case SIMPLE_AVG:
            {
                rho_average = 0.5*(rho[idx_minus] + rho[idx_plus]);
                c_average = 0.5*(c[idx_sound_speed_minus] + c[idx_sound_speed_plus]);
                
                break;
            }
            case ROE_AVG:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::computeLocalFaceProjectionMatrixOfPrimitiveVariables()\n"
                    << "Roe averaging is not yet implemented."
                    << std::endl);
                
                rho_average = 0.0;
                c_average   = 0.0;
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::computeLocalFaceProjectionMatrixOfPrimitiveVariables()\n"
                    << "Unknown d_proj_mat_primitive_var_averaging given."
                    << std::endl);
                
                rho_average = 0.0;
                c_average   = 0.0;
            }
        }
        
        switch (direction)
        {
            case X_DIRECTION:
            {
                projection_matrix[0][0] = 0.0;
                projection_matrix[0][1] = -0.5*rho_average*c_average;
                projection_matrix[0][2] = 0.0;
                projection_matrix[0][3] = 0.0;
                projection_matrix[0][4] = 0.5;
                projection_matrix[1][0] = 1.0;
                projection_matrix[1][1] = 0.0;
                projection_matrix[1][2] = 0.0;
                projection_matrix[1][3] = 0.0;
                projection_matrix[1][4] = -1.0/(c_average*c_average);
                projection_matrix[2][0] = 0.0;
                projection_matrix[2][1] = 0.0;
                projection_matrix[2][2] = 1.0;
                projection_matrix[2][3] = 0.0;
                projection_matrix[2][4] = 0.0;
                projection_matrix[3][0] = 0.0;
                projection_matrix[3][1] = 0.0;
                projection_matrix[3][2] = 0.0;
                projection_matrix[3][3] = 1.0;
                projection_matrix[3][4] = 0.0;
                projection_matrix[4][0] = 0.0;
                projection_matrix[4][1] = 0.5*rho_average*c_average;
                projection_matrix[4][2] = 0.0;
                projection_matrix[4][3] = 0.0;
                projection_matrix[4][4] = 0.5;
                
                break;
            }
            case Y_DIRECTION:
            {
                projection_matrix[0][0] = 0.0;
                projection_matrix[0][1] = 0.0;
                projection_matrix[0][2] = -0.5*rho_average*c_average;
                projection_matrix[0][3] = 0.0;
                projection_matrix[0][4] = 0.5;
                projection_matrix[1][0] = 1.0;
                projection_matrix[1][1] = 0.0;
                projection_matrix[1][2] = 0.0;
                projection_matrix[1][3] = 0.0;
                projection_matrix[1][4] = -1.0/(c_average*c_average);
                projection_matrix[2][0] = 0.0;
                projection_matrix[2][1] = 1.0;
                projection_matrix[2][2] = 0.0;
                projection_matrix[2][3] = 0.0;
                projection_matrix[2][4] = 0.0;
                projection_matrix[3][0] = 0.0;
                projection_matrix[3][1] = 0.0;
                projection_matrix[3][2] = 0.0;
                projection_matrix[3][3] = 1.0;
                projection_matrix[3][4] = 0.0;
                projection_matrix[4][0] = 0.0;
                projection_matrix[4][1] = 0.0;
                projection_matrix[4][2] = 0.5*rho_average*c_average;
                projection_matrix[4][3] = 0.0;
                projection_matrix[4][4] = 0.5;
                
                break;
            }
            case Z_DIRECTION:
            {
                projection_matrix[0][0] = 0.0;
                projection_matrix[0][1] = 0.0;
                projection_matrix[0][2] = 0.0;
                projection_matrix[0][3] = -0.5*rho_average*c_average;
                projection_matrix[0][4] = 0.5;
                projection_matrix[1][0] = 1.0;
                projection_matrix[1][1] = 0.0;
                projection_matrix[1][2] = 0.0;
                projection_matrix[1][3] = 0.0;
                projection_matrix[1][4] = -1.0/(c_average*c_average);
                projection_matrix[2][0] = 0.0;
                projection_matrix[2][1] = 1.0;
                projection_matrix[2][2] = 0.0;
                projection_matrix[2][3] = 0.0;
                projection_matrix[2][4] = 0.0;
                projection_matrix[3][0] = 0.0;
                projection_matrix[3][1] = 0.0;
                projection_matrix[3][2] = 1.0;
                projection_matrix[3][3] = 0.0;
                projection_matrix[3][4] = 0.0;
                projection_matrix[4][0] = 0.0;
                projection_matrix[4][1] = 0.0;
                projection_matrix[4][2] = 0.0;
                projection_matrix[4][3] = 0.5*rho_average*c_average;
                projection_matrix[4][4] = 0.5;
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::computeLocalFaceProjectionMatrixOfPrimitiveVariables()\n"
                << "There are only x-direction, y-direction and z-direction for three-dimensional problem."
                << std::endl);
            }
        }
    }
}


/*
 * Compute the local face data of inverse of projection matrix of primitive variables
 * in the registered patch.
 */
void
FlowModelSingleSpecies::computeLocalFaceProjectionMatrixInverseOfPrimitiveVariables(
    boost::multi_array<double, 2>& projection_matrix_inv,
    const hier::Index& cell_index_minus,
    const hier::Index& cell_index_plus,
    const DIRECTION& direction)
{
    if (!d_proj_mat_primitive_var_registered)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::computeLocalFaceProjectionMatrixInverseOfPrimitiveVariables()\n"
            << "Projection matrices is not yet registered."
            << std::endl);
    }
    
    projection_matrix_inv.resize(boost::extents[d_num_eqn][d_num_eqn]);
    
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

    // Compute the inverse of the projection matrix.
    if (d_dim == tbox::Dimension(1))
    {
        // Compute the linear indices.
        const int idx_minus = cell_index_minus[0] + d_num_ghosts[0];
        const int idx_plus = cell_index_plus[0] + d_num_ghosts[0];
        const int idx_sound_speed_minus = cell_index_minus[0] + d_num_subghosts_sound_speed[0];
        const int idx_sound_speed_plus = cell_index_plus[0] + d_num_subghosts_sound_speed[0];
        
        // Compute the average values.
        double rho_average, c_average;
        switch (d_proj_mat_primitive_var_averaging)
        {
            case SIMPLE_AVG:
            {
                rho_average = 0.5*(rho[idx_minus] + rho[idx_plus]);
                c_average = 0.5*(c[idx_sound_speed_minus] + c[idx_sound_speed_plus]);
                
                break;
            }
            case ROE_AVG:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::"
                    << "computeLocalFaceProjectionMatrixInverseOfPrimitiveVariables()\n"
                    << "Roe averaging is not yet implemented."
                    << std::endl);
                
                rho_average = 0.0;
                c_average   = 0.0;
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::"
                    << "computeLocalFaceProjectionMatrixInverseOfPrimitiveVariables()\n"
                    << "Unknown d_proj_mat_primitive_var_averaging given."
                    << std::endl);
                
                rho_average = 0.0;
                c_average   = 0.0;
            }
        }
        
        switch (direction)
        {
            case X_DIRECTION:
            {
                projection_matrix_inv[0][0] = 1.0/(c_average*c_average);
                projection_matrix_inv[0][1] = 1.0;
                projection_matrix_inv[0][2] = 1.0/(c_average*c_average);
                projection_matrix_inv[1][0] = -1.0/(rho_average*c_average);
                projection_matrix_inv[1][1] = 0.0;
                projection_matrix_inv[1][2] = 1.0/(rho_average*c_average);
                projection_matrix_inv[2][0] = 1.0;
                projection_matrix_inv[2][1] = 0.0;
                projection_matrix_inv[2][2] = 1.0;
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::computeLocalFaceProjectionMatrixInverseOfPrimitiveVariables()\n"
                << "There is only x-direction for one-dimensional problem."
                << std::endl);
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int idx_minus = (cell_index_minus[0] + d_num_ghosts[0]) +
            (cell_index_minus[1] + d_num_ghosts[1])*d_ghostcell_dims[0];
        
        const int idx_plus = (cell_index_plus[0] + d_num_ghosts[0]) +
            (cell_index_plus[1] + d_num_ghosts[1])*d_ghostcell_dims[0];
        
        const int idx_sound_speed_minus = (cell_index_minus[0] + d_num_subghosts_sound_speed[0]) +
            (cell_index_minus[1] + d_num_subghosts_sound_speed[1])*d_subghostcell_dims_sound_speed[0];
        
        const int idx_sound_speed_plus = (cell_index_plus[0] + d_num_subghosts_sound_speed[0]) +
            (cell_index_plus[1] + d_num_subghosts_sound_speed[1])*d_subghostcell_dims_sound_speed[0];
        
        // Compute the average values.
        double rho_average, c_average;
        switch (d_proj_mat_primitive_var_averaging)
        {
            case SIMPLE_AVG:
            {
                rho_average = 0.5*(rho[idx_minus] + rho[idx_plus]);
                c_average = 0.5*(c[idx_sound_speed_minus] + c[idx_sound_speed_plus]);
                
                break;
            }
            case ROE_AVG:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::"
                    << "computeLocalFaceProjectionMatrixInverseOfPrimitiveVariables()\n"
                    << "Roe averaging is not yet implemented."
                    << std::endl);
                
                rho_average = 0.0;
                c_average   = 0.0;
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::"
                    << "computeLocalFaceProjectionMatrixInverseOfPrimitiveVariables()\n"
                    << "Unknown d_proj_mat_primitive_var_averaging given."
                    << std::endl);
                
                rho_average = 0.0;
                c_average   = 0.0;
            }
        }
        
        switch (direction)
        {
            case X_DIRECTION:
            {
                projection_matrix_inv[0][0] = 1.0/(c_average*c_average);
                projection_matrix_inv[0][1] = 1.0;
                projection_matrix_inv[0][2] = 0.0;
                projection_matrix_inv[0][3] = 1.0/(c_average*c_average);
                projection_matrix_inv[1][0] = -1.0/(rho_average*c_average);
                projection_matrix_inv[1][1] = 0.0;
                projection_matrix_inv[1][2] = 0.0;
                projection_matrix_inv[1][3] = 1.0/(rho_average*c_average);
                projection_matrix_inv[2][0] = 0.0;
                projection_matrix_inv[2][1] = 0.0;
                projection_matrix_inv[2][2] = 1.0;
                projection_matrix_inv[2][3] = 0.0;
                projection_matrix_inv[3][0] = 1.0;
                projection_matrix_inv[3][1] = 0.0;
                projection_matrix_inv[3][2] = 0.0;
                projection_matrix_inv[3][3] = 1.0;
                
                break;
            }
            case Y_DIRECTION:
            {
                projection_matrix_inv[0][0] = 1.0/(c_average*c_average);
                projection_matrix_inv[0][1] = 1.0;
                projection_matrix_inv[0][2] = 0.0;
                projection_matrix_inv[0][3] = 1.0/(c_average*c_average);
                projection_matrix_inv[1][0] = 0.0;
                projection_matrix_inv[1][1] = 0.0;
                projection_matrix_inv[1][2] = 1.0;
                projection_matrix_inv[1][3] = 0.0;
                projection_matrix_inv[2][0] = -1.0/(rho_average*c_average);
                projection_matrix_inv[2][1] = 0.0;
                projection_matrix_inv[2][2] = 0.0;
                projection_matrix_inv[2][3] = 1.0/(rho_average*c_average);
                projection_matrix_inv[3][0] = 1.0;
                projection_matrix_inv[3][1] = 0.0;
                projection_matrix_inv[3][2] = 0.0;
                projection_matrix_inv[3][3] = 1.0;
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::computeLocalFaceProjectionMatrixInverseOfPrimitiveVariables()\n"
                << "There are only x-direction and y-direction for two-dimensional problem."
                << std::endl);
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int idx_minus = (cell_index_minus[0] + d_num_ghosts[0]) +
            (cell_index_minus[1] + d_num_ghosts[1])*d_ghostcell_dims[0] +
            (cell_index_minus[2] + d_num_ghosts[2])*d_ghostcell_dims[0]*
                d_ghostcell_dims[1];
        
        const int idx_plus = (cell_index_plus[0] + d_num_ghosts[0]) +
            (cell_index_plus[1] + d_num_ghosts[1])*d_ghostcell_dims[0] +
            (cell_index_plus[2] + d_num_ghosts[2])*d_ghostcell_dims[0]*
                d_ghostcell_dims[1];
        
        const int idx_sound_speed_minus = (cell_index_minus[0] + d_num_subghosts_sound_speed[0]) +
            (cell_index_minus[1] + d_num_subghosts_sound_speed[1])*d_subghostcell_dims_sound_speed[0] +
            (cell_index_minus[2] + d_num_subghosts_sound_speed[2])*d_subghostcell_dims_sound_speed[0]*
                d_subghostcell_dims_sound_speed[1];
        
        const int idx_sound_speed_plus = (cell_index_plus[0] + d_num_subghosts_sound_speed[0]) +
            (cell_index_plus[1] + d_num_subghosts_sound_speed[1])*d_subghostcell_dims_sound_speed[0] +
            (cell_index_plus[2] + d_num_subghosts_sound_speed[2])*d_subghostcell_dims_sound_speed[0]*
                d_subghostcell_dims_sound_speed[1];
        
        // Compute the average values.
        double rho_average, c_average;
        switch (d_proj_mat_primitive_var_averaging)
        {
            case SIMPLE_AVG:
            {
                rho_average = 0.5*(rho[idx_minus] + rho[idx_plus]);
                c_average = 0.5*(c[idx_sound_speed_minus] + c[idx_sound_speed_plus]);
                
                break;
            }
            case ROE_AVG:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::"
                    << "computeLocalFaceProjectionMatrixInverseOfPrimitiveVariables()\n"
                    << "Roe averaging is not yet implemented."
                    << std::endl);
                
                rho_average = 0.0;
                c_average   = 0.0;
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::"
                    << "computeLocalFaceProjectionMatrixInverseOfPrimitiveVariables()\n"
                    << "Unknown d_proj_mat_primitive_var_averaging given."
                    << std::endl);
                
                rho_average = 0.0;
                c_average   = 0.0;
            }
        }
        
        switch (direction)
        {
            case X_DIRECTION:
            {
                projection_matrix_inv[0][0] = 1.0/(c_average*c_average);
                projection_matrix_inv[0][1] = 1.0;
                projection_matrix_inv[0][2] = 0.0;
                projection_matrix_inv[0][3] = 0.0;
                projection_matrix_inv[0][4] = 1.0/(c_average*c_average);
                projection_matrix_inv[1][0] = -1.0/(rho_average*c_average);
                projection_matrix_inv[1][1] = 0.0;
                projection_matrix_inv[1][2] = 0.0;
                projection_matrix_inv[1][3] = 0.0;
                projection_matrix_inv[1][4] = 1.0/(rho_average*c_average);
                projection_matrix_inv[2][0] = 0.0;
                projection_matrix_inv[2][1] = 0.0;
                projection_matrix_inv[2][2] = 1.0;
                projection_matrix_inv[2][3] = 0.0;
                projection_matrix_inv[2][4] = 0.0;
                projection_matrix_inv[3][0] = 0.0;
                projection_matrix_inv[3][1] = 0.0;
                projection_matrix_inv[3][2] = 0.0;
                projection_matrix_inv[3][3] = 1.0;
                projection_matrix_inv[3][4] = 0.0;
                projection_matrix_inv[4][0] = 1.0;
                projection_matrix_inv[4][1] = 0.0;
                projection_matrix_inv[4][2] = 0.0;
                projection_matrix_inv[4][3] = 0.0;
                projection_matrix_inv[4][4] = 1.0;
                
                break;
            }
            case Y_DIRECTION:
            {
                projection_matrix_inv[0][0] = 1.0/(c_average*c_average);
                projection_matrix_inv[0][1] = 1.0;
                projection_matrix_inv[0][2] = 0.0;
                projection_matrix_inv[0][3] = 0.0;
                projection_matrix_inv[0][4] = 1.0/(c_average*c_average);
                projection_matrix_inv[1][0] = 0.0;
                projection_matrix_inv[1][1] = 0.0;
                projection_matrix_inv[1][2] = 1.0;
                projection_matrix_inv[1][3] = 0.0;
                projection_matrix_inv[1][4] = 0.0;
                projection_matrix_inv[2][0] = -1.0/(rho_average*c_average);
                projection_matrix_inv[2][1] = 0.0;
                projection_matrix_inv[2][2] = 0.0;
                projection_matrix_inv[2][3] = 0.0;
                projection_matrix_inv[2][4] = 1.0/(rho_average*c_average);
                projection_matrix_inv[3][0] = 0.0;
                projection_matrix_inv[3][1] = 0.0;
                projection_matrix_inv[3][2] = 0.0;
                projection_matrix_inv[3][3] = 1.0;
                projection_matrix_inv[3][4] = 0.0;
                projection_matrix_inv[4][0] = 1.0;
                projection_matrix_inv[4][1] = 0.0;
                projection_matrix_inv[4][2] = 0.0;
                projection_matrix_inv[4][3] = 0.0;
                projection_matrix_inv[4][4] = 1.0;
                
                break;
            }
            case Z_DIRECTION:
            {
                projection_matrix_inv[0][0] = 1.0/(c_average*c_average);
                projection_matrix_inv[0][1] = 1.0;
                projection_matrix_inv[0][2] = 0.0;
                projection_matrix_inv[0][3] = 0.0;
                projection_matrix_inv[0][4] = 1.0/(c_average*c_average);
                projection_matrix_inv[1][0] = 0.0;
                projection_matrix_inv[1][1] = 0.0;
                projection_matrix_inv[1][2] = 1.0;
                projection_matrix_inv[1][3] = 0.0;
                projection_matrix_inv[1][4] = 0.0;
                projection_matrix_inv[2][0] = 0.0;
                projection_matrix_inv[2][1] = 0.0;
                projection_matrix_inv[2][2] = 0.0;
                projection_matrix_inv[2][3] = 1.0;
                projection_matrix_inv[2][4] = 0.0;
                projection_matrix_inv[3][0] = -1.0/(rho_average*c_average);
                projection_matrix_inv[3][1] = 0.0;
                projection_matrix_inv[3][2] = 0.0;
                projection_matrix_inv[3][3] = 0.0;
                projection_matrix_inv[3][4] = 1.0/(rho_average*c_average);
                projection_matrix_inv[4][0] = 1.0;
                projection_matrix_inv[4][1] = 0.0;
                projection_matrix_inv[4][2] = 0.0;
                projection_matrix_inv[4][3] = 0.0;
                projection_matrix_inv[4][4] = 1.0;
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                << ": FlowModelSingleSpecies::computeLocalFaceProjectionMatrixInverseOfPrimitiveVariables()\n"
                << "There are only x-direction, y-direction and z-direction for three-dimensional problem."
                << std::endl);
            }
        }
    }
}


/*
 * Compute the local intercell quantities with conservative variables on each side of the face
 * from Riemann solver at face.
 * fluxes_face: Convective flux at face.
 * velocity_face: Velocity at face.
 * The FlowModelSingleSpecies class modifies nothing for velocity_face.
 */
void
FlowModelSingleSpecies::computeLocalFaceFluxAndVelocityFromRiemannSolverWithConservativeVariables(
    std::vector<boost::reference_wrapper<double> >& flux_face,
    std::vector<boost::reference_wrapper<double> >& velocity_face,
    const std::vector<boost::reference_wrapper<double> >& conservative_variables_minus,
    const std::vector<boost::reference_wrapper<double> >& conservative_variables_plus,
    const DIRECTION& direction,
    const RIEMANN_SOLVER& Riemann_solver)
{
    switch (Riemann_solver)
    {
        case HLLC_RIEMANN_SOLVER:
        {
            d_Riemann_solver_HLLC->computeIntercellFluxFromConservativeVariables(
                flux_face,
                conservative_variables_minus,
                conservative_variables_plus,
                direction);
            
            break;
        }
        case HLLC_HLL_RIEMANN_SOLVER:
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
 * fluxes_face: Convective flux at face.
 * velocity_face: Velocity at face.
 * The FlowModelSingleSpecies class modifies nothing for velocity_face.
 */
void
FlowModelSingleSpecies::computeLocalFaceFluxAndVelocityFromRiemannSolverWithPrimitiveVariables(
    std::vector<boost::reference_wrapper<double> >& flux_face,
    std::vector<boost::reference_wrapper<double> >& velocity_face,
    const std::vector<boost::reference_wrapper<double> >& primitive_variables_minus,
    const std::vector<boost::reference_wrapper<double> >& primitive_variables_plus,
    const DIRECTION& direction,
    const RIEMANN_SOLVER& Riemann_solver)
{
    switch (Riemann_solver)
    {
        case HLLC_RIEMANN_SOLVER:
        {
            d_Riemann_solver_HLLC->computeIntercellFluxFromPrimitiveVariables(
                flux_face,
                primitive_variables_minus,
                primitive_variables_plus,
                direction);
            
            break;
        }
        case HLLC_HLL_RIEMANN_SOLVER:
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
 * Check whether the given conservative variables are within the bounds.
 */
bool
FlowModelSingleSpecies::haveConservativeVariablesBounded(const std::vector<double>& conservative_variables)
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(conservative_variables.size()) == d_num_eqn);
#endif
    
    bool are_bounded = true;
    
    // Check if the density is bounded.
    if (conservative_variables[0] < 0.0)
    {
        are_bounded = false;
    }
    
    // Check if the total energy is bounded.
    if (conservative_variables[d_num_eqn - 1] < 0.0)
    {
        are_bounded = false;
    }
    
    return are_bounded;
}


/*
 * Check whether the given primitive variables are within the bounds.
 */
bool
FlowModelSingleSpecies::havePrimitiveVariablesBounded(const std::vector<double>& primitive_variables)
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(primitive_variables.size()) == d_num_eqn);
#endif
    
    bool are_bounded = true;
    
    // Check if the density is bounded.
    if (primitive_variables[0] < 0.0)
    {
        are_bounded = false;
    }
    
    // Check if the pressure is bounded.
    if (primitive_variables[d_num_eqn - 1] < 0.0)
    {
        are_bounded = false;
    }
    
    return are_bounded;
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
    
    // Compute the pressure.
    const double p = d_equation_of_state->getPressure(
        Q[0],
        m_ptr,
        Q[1 + d_dim.getValue()],
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
    
    // Get the pointers to the velocity components.
    std::vector<const double*> vel_ptr;
    vel_ptr.reserve(d_dim.getValue());
    for (int di = 0; di < d_dim.getValue(); di++)
    {
        vel_ptr.push_back(V[1 + di]);
    }
    
    // Get the thermodynamic properties of the species.
    std::vector<const double*> thermo_properties_ptr;
    thermo_properties_ptr.reserve(static_cast<int> (d_thermo_properties.size()));
    for (int ti = 0; ti < static_cast<int> (d_thermo_properties.size()); ti++)
    {
        thermo_properties_ptr.push_back(&d_thermo_properties[ti]);
    }
    
    // Compute the total energy.
    const double E = d_equation_of_state->getTotalEnergy(
        V[0],
        vel_ptr,
        V[1 + d_dim.getValue()],
        thermo_properties_ptr);
    
    // Convert the primitive variables to conservative variables.
    *Q[0] = *V[0];
    for (int di = 0; di < d_dim.getValue(); di++)
    {
        *Q[1 + di] = (*V[0])*(*V[1 + di]);
    }
    *Q[1 + d_dim.getValue()] = E;
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
    
#ifdef DEBUG_CHECK_ASSERTIONS
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
        
#ifdef DEBUG_CHECK_ASSERTIONS
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
                
                std::vector<const double*> m_ptr;
                m_ptr.push_back(&rho_u[idx_data]);
                
                buffer[idx_region] = d_equation_of_state->getPressure(
                    &rho[idx_data],
                    m_ptr,
                    &E[idx_data],
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
                    
                    std::vector<const double*> m_ptr;
                    m_ptr.push_back(&rho_u[idx_data]);
                    m_ptr.push_back(&rho_v[idx_data]);
                    
                    buffer[idx_region] = d_equation_of_state->getPressure(
                        &rho[idx_data],
                        m_ptr,
                        &E[idx_data],
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
                        
                        std::vector<const double*> m_ptr;
                        m_ptr.push_back(&rho_u[idx_data]);
                        m_ptr.push_back(&rho_v[idx_data]);
                        m_ptr.push_back(&rho_w[idx_data]);
                        
                        buffer[idx_region] = d_equation_of_state->getPressure(
                            &rho[idx_data],
                            m_ptr,
                            &E[idx_data],
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
        
#ifdef DEBUG_CHECK_ASSERTIONS
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
                
                std::vector<const double*> m_ptr;
                m_ptr.push_back(&rho_u[idx_data]);
                
                const double p = d_equation_of_state->getPressure(
                    &rho[idx_data],
                    m_ptr,
                    &E[idx_data],
                    thermo_properties_ptr);
                
                buffer[idx_region] = d_equation_of_state->getSoundSpeed(
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
                    
                    std::vector<const double*> m_ptr;
                    m_ptr.push_back(&rho_u[idx_data]);
                    m_ptr.push_back(&rho_v[idx_data]);
                    
                    const double p = d_equation_of_state->getPressure(
                        &rho[idx_data],
                        m_ptr,
                        &E[idx_data],
                        thermo_properties_ptr);
                    
                    buffer[idx_region] = d_equation_of_state->getSoundSpeed(
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
                        
                        std::vector<const double*> m_ptr;
                        m_ptr.push_back(&rho_u[idx_data]);
                        m_ptr.push_back(&rho_v[idx_data]);
                        m_ptr.push_back(&rho_w[idx_data]);
                        
                        const double p = d_equation_of_state->getPressure(
                            &rho[idx_data],
                            m_ptr,
                            &E[idx_data],
                            thermo_properties_ptr);
                        
                        buffer[idx_region] = d_equation_of_state->getSoundSpeed(
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
        
#ifdef DEBUG_CHECK_ASSERTIONS
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
 * Set the ghost boxes and their dimensions of derived cell variables.
 */
void
FlowModelSingleSpecies::setGhostBoxesAndDimensionsDerivedCellVariables()
{
    if (d_num_subghosts_pressure > -hier::IntVector::getOne(d_dim))
    {
        d_subghost_box_pressure = d_interior_box;
        d_subghost_box_pressure.grow(d_num_subghosts_pressure);
        d_subghostcell_dims_pressure = d_subghost_box_pressure.numberCells();
    }
    
    if (d_num_subghosts_velocity > -hier::IntVector::getOne(d_dim))
    {
        d_subghost_box_velocity = d_interior_box;
        d_subghost_box_velocity.grow(d_num_subghosts_velocity);
        d_subghostcell_dims_velocity = d_subghost_box_velocity.numberCells();
    }
    
    if (d_num_subghosts_sound_speed > -hier::IntVector::getOne(d_dim))
    {
        d_subghost_box_sound_speed = d_interior_box;
        d_subghost_box_sound_speed.grow(d_num_subghosts_sound_speed);
        d_subghostcell_dims_sound_speed = d_subghost_box_sound_speed.numberCells();
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
 * Compute the global cell data of pressure in the registered patch.
 */
void
FlowModelSingleSpecies::computeGlobalCellDataPressure()
{
    if (d_num_subghosts_pressure > -hier::IntVector::getOne(d_dim))
    {
        // Create the cell data of pressure.
        d_data_pressure.reset(
            new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_pressure));
        
        // Get the cell data of the variables density, momentum and total energy.
        boost::shared_ptr<pdat::CellData<double> > data_density =
            getGlobalCellDataDensity();
        
        boost::shared_ptr<pdat::CellData<double> > data_momentum =
            getGlobalCellDataMomentum();
        
        boost::shared_ptr<pdat::CellData<double> > data_total_energy =
            getGlobalCellDataTotalEnergy();
        
        // Get the pointers to the cell data of pressure, density and total energy.
        double* p   = d_data_pressure->getPointer(0);
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
            // Get the pointer to cell data of momentum.
            double* rho_u = data_momentum->getPointer(0);
            
            // Compute the pressure field.
            for (int i = -d_num_subghosts_pressure[0];
                 i < d_interior_dims[0] + d_num_subghosts_pressure[0];
                 i++)
            {
                // Compute the linear indices.
                const int idx = i + d_num_ghosts[0];
                const int idx_pressure = i + d_num_subghosts_pressure[0];
                
                std::vector<const double*> m_ptr;
                m_ptr.reserve(1);
                m_ptr.push_back(&rho_u[idx]);
                
                p[idx_pressure] = d_equation_of_state->
                    getPressure(
                        &rho[idx],
                        m_ptr,
                        &E[idx],
                        thermo_properties_ptr);
            }
        }
        else if (d_dim == tbox::Dimension(2))
        {
            // Get the pointers to the cell data of momentum.
            double* rho_u = data_momentum->getPointer(0);
            double* rho_v = data_momentum->getPointer(1);
            
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
                    
                    const int idx_pressure = (i + d_num_subghosts_pressure[0]) +
                        (j + d_num_subghosts_pressure[1])*d_subghostcell_dims_pressure[0];
                    
                    std::vector<const double*> m_ptr;
                    m_ptr.reserve(2);
                    m_ptr.push_back(&rho_u[idx]);
                    m_ptr.push_back(&rho_v[idx]);
                    
                    p[idx_pressure] = d_equation_of_state->
                        getPressure(
                            &rho[idx],
                            m_ptr,
                            &E[idx],
                            thermo_properties_ptr);
                }
            }
        }
        else if (d_dim == tbox::Dimension(3))
        {
            // Get the pointers to the cell data of momentum.
            double* rho_u = data_momentum->getPointer(0);
            double* rho_v = data_momentum->getPointer(1);
            double* rho_w = data_momentum->getPointer(2);
            
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
                        
                        std::vector<const double*> m_ptr;
                        m_ptr.reserve(3);
                        m_ptr.push_back(&rho_u[idx]);
                        m_ptr.push_back(&rho_v[idx]);
                        m_ptr.push_back(&rho_w[idx]);
                        
                        p[idx_pressure] = d_equation_of_state->
                            getPressure(
                                &rho[idx],
                                m_ptr,
                                &E[idx],
                                thermo_properties_ptr);
                    }
                }
            }
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSingleSpecies::computeGlobalCellDataPressure()\n"
            << "Cell data of 'PRESSURE' is not yet registered."
            << std::endl);
    }
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
            computeGlobalCellDataPressure();
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
                
                c[idx_sound_speed] = d_equation_of_state->
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
                    
                    c[idx_sound_speed] = d_equation_of_state->
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
                        
                        c[idx_sound_speed] = d_equation_of_state->
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
 * Compute the global cell data of enstrophy with velocity and vorticity in the registered patch.
 */
void
FlowModelSingleSpecies::computeGlobalCellDataEnstrophyWithVelocityAndVorticity()
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
                    
                    Omega[idx_enstrophy] = 0.5*omega[idx_vorticity]*omega[idx_vorticity];
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
                        
                        Omega[idx_enstrophy] = 0.5*(omega_x[idx_vorticity]*omega_x[idx_vorticity] +
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
            << ": FlowModelSingleSpecies::computeGlobalCellDataEnstrophyWithVelocityAndVorticity()\n"
            << "Cell data of 'ENSTROPHY' is not yet registered."
            << std::endl);
    }
}


/*
 * Compute the global cell data of convective flux with pressure and velocity in the registered patch.
 */
void
FlowModelSingleSpecies::computeGlobalCellDataConvectiveFluxWithPressureAndVelocity(DIRECTION direction)
{
    if (direction == X_DIRECTION)
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
            
            if (!d_data_pressure)
            {
                computeGlobalCellDataPressure();
            }
            
            if (!d_data_velocity)
            {
                computeGlobalCellDataVelocity();
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
                << ": FlowModelSingleSpecies::computeGlobalCellDataConvectiveFluxWithPressureAndVelocity()\n"
                << "Cell data of 'CONVECTIVE_FLUX_X' is not yet registered."
                << std::endl);
        }
    }
    else if (direction == Y_DIRECTION)
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
            
            if (!d_data_pressure)
            {
                computeGlobalCellDataPressure();
            }
            
            if (!d_data_velocity)
            {
                computeGlobalCellDataVelocity();
            }
            
            // Get the pointers to the cell data of total energy and pressure.
            double* E   = data_total_energy->getPointer(0);
            double* p   = d_data_pressure->getPointer(0);
            
            if (d_dim == tbox::Dimension(1))
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::computeGlobalCellDataConvectiveFluxWithPressureAndVelocity()\n"
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
                << ": FlowModelSingleSpecies::computeGlobalCellDataConvectiveFluxWithPressureAndVelocity()\n"
                << "Cell data of 'CONVECTIVE_FLUX_Y' is not yet registered."
                << std::endl);
        }
    }
    else if (direction == Z_DIRECTION)
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
            
            if (!d_data_pressure)
            {
                computeGlobalCellDataPressure();
            }
            
            if (!d_data_velocity)
            {
                computeGlobalCellDataVelocity();
            }
            
            // Get the pointers to the cell data of total energy and pressure.
            double* E   = data_total_energy->getPointer(0);
            double* p   = d_data_pressure->getPointer(0);
            
            if (d_dim == tbox::Dimension(1) || d_dim == tbox::Dimension(2))
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::computeGlobalCellDataConvectiveFluxWithPressureAndVelocity()\n"
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
                << ": FlowModelSingleSpecies::computeGlobalCellDataConvectiveFluxWithPressureAndVelocity()\n"
                << "Cell data of 'CONVECTIVE_FLUX_Z' is not yet registered."
                << std::endl);
        }
    }
}


/*
 * Compute the global cell data of maximum wave speed with pressure, velocity and sound speed in the registered patch.
 */
void
FlowModelSingleSpecies::computeGlobalCellDataMaxWaveSpeedWithPressureVelocityAndSoundSpeed(DIRECTION direction)
{
    if (direction == X_DIRECTION)
    {
        if (d_num_subghosts_max_wave_speed_x > -hier::IntVector::getOne(d_dim))
        {
            // Create the cell data of maximum wave speed in the x-direction.
            d_data_max_wave_speed_x.reset(
                new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_max_wave_speed_x));
            
            if (!d_data_sound_speed)
            {
                computeGlobalCellDataSoundSpeedWithPressure();
            }
            
            if (!d_data_velocity)
            {
                computeGlobalCellDataVelocity();
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
                << ": FlowModelSingleSpecies::computeGlobalCellDataMaxWaveSpeedWithPressureVelocityAndSoundSpeed()\n"
                << "Cell data of 'MAX_WAVE_SPEED_X' is not yet registered."
                << std::endl);
        }
    }
    else if (direction == Y_DIRECTION)
    {
        if (d_num_subghosts_max_wave_speed_y > -hier::IntVector::getOne(d_dim))
        {
            // Create the cell data of maximum wave speed in the y-direction.
            d_data_max_wave_speed_y.reset(
                new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_max_wave_speed_y));
            
            if (!d_data_sound_speed)
            {
                computeGlobalCellDataSoundSpeedWithPressure();
            }
            
            if (!d_data_velocity)
            {
                computeGlobalCellDataVelocity();
            }
            
            // Get the pointers to the cell data of maximum wave speed and velocity in y-direction, and sound speed.
            double* lambda_max_y = d_data_max_wave_speed_y->getPointer(0);
            double* v            = d_data_velocity->getPointer(1);
            double* c            = d_data_sound_speed->getPointer(0);
            
            if (d_dim == tbox::Dimension(1))
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::computeGlobalCellDataMaxWaveSpeedWithPressureVelocityAndSoundSpeed()\n"
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
                << ": FlowModelSingleSpecies::computeGlobalCellDataMaxWaveSpeedWithPressureVelocityAndSoundSpeed()\n"
                << "Cell data of 'MAX_WAVE_SPEED_Y' is not yet registered."
                << std::endl);
        }
    }
    else if (direction == Z_DIRECTION)
    {
        if (d_num_subghosts_max_wave_speed_z > -hier::IntVector::getOne(d_dim))
        {
            // Create the cell data of maximum wave speed in the z-direction.
            d_data_max_wave_speed_z.reset(
                new pdat::CellData<double>(d_interior_box, 1, d_num_subghosts_max_wave_speed_z));
            
            if (!d_data_sound_speed)
            {
                computeGlobalCellDataSoundSpeedWithPressure();
            }
            
            if (!d_data_velocity)
            {
                computeGlobalCellDataVelocity();
            }
            
            // Get the pointers to the cell data of maximum wave speed and velocity in z-direction, and sound speed.
            double* lambda_max_z = d_data_max_wave_speed_z->getPointer(0);
            double* w            = d_data_velocity->getPointer(2);
            double* c            = d_data_sound_speed->getPointer(0);
            
            if (d_dim == tbox::Dimension(1) || d_dim == tbox::Dimension(2))
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSingleSpecies::computeGlobalCellDataMaxWaveSpeedWithPressureVelocityAndSoundSpeed()\n"
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
                << ": FlowModelSingleSpecies::computeGlobalCellDataMaxWaveSpeedWithPressureVelocityAndSoundSpeed()\n"
                << "Cell data of 'MAX_WAVE_SPEED_Z' is not yet registered."
                << std::endl);
        }
    }
}
