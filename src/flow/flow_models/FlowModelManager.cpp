#include "flow/flow_models/FlowModelManager.hpp"

FlowModelManager::FlowModelManager(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_species,
    const boost::shared_ptr<EquationOfState>& equation_of_state,
    const std::string& flow_model_str):
        d_object_name(object_name),
        d_dim(dim),
        d_grid_geometry(grid_geometry),
        d_num_species(num_species),
        d_equation_of_state(equation_of_state)
{
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(grid_geometry);
    TBOX_ASSERT(equation_of_state);
    
    if (flow_model_str == "SINGLE_SPECIES")
    {
        d_flow_model_label = SINGLE_SPECIES;
        
        d_flow_model.reset(new FlowModelSingleSpecies(
            "d_flow_model",
            d_dim,
            d_grid_geometry,
            hier::IntVector::getZero(d_dim),
            d_num_species,
            d_equation_of_state));
    }
    else if (flow_model_str == "FOUR_EQN_CONSERVATIVE")
    {
        d_flow_model_label = FOUR_EQN_CONSERVATIVE;
        
        d_flow_model.reset(new FlowModelFourEqnConservative(
            "d_flow_model",
            d_dim,
            d_grid_geometry,
            hier::IntVector::getZero(d_dim),
            d_num_species,
            d_equation_of_state));
    }
    else if (flow_model_str == "FIVE_EQN_ALLAIRE")
    {
        d_flow_model_label = FIVE_EQN_ALLAIRE;
        
        d_flow_model.reset(new FlowModelFiveEqnAllaire(
            "d_flow_model",
            d_dim,
            d_grid_geometry,
            hier::IntVector::getZero(d_dim),
            d_num_species,
            d_equation_of_state));
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Unknown flow_model/d_flow_model string = '"
            << flow_model_str
            << "' found in input/restart file."
            << std::endl);        
    }
    
    if (d_num_species > 1 && d_flow_model_label == SINGLE_SPECIES)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Number of species = "
            << d_num_species
            << " shouldn't use single-species model."
            << std::endl); 
    }
}


/*
 * Initialize d_initial_conditions.
 */
void
FlowModelManager::initializeInitialConditions(
    const std::string& project_name,
    boost::shared_ptr<InitialConditions>& initial_conditions)
{
    initial_conditions.reset(new InitialConditions(
        "initial conditions",
        project_name,
        d_dim,
        d_grid_geometry,
        d_flow_model_label,
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
 * Print all characteristics of flow model manager.
 */
void
FlowModelManager::printClassData(std::ostream& os) const
{
    os << "\nPrint FlowModelManager object..."
       << std::endl;
    
    os << std::endl;
    
    os << "FlowModelManager: this = "
       << (FlowModelManager *)this
       << std::endl;
    
    os << "d_object_name = "
       << d_object_name
       << std::endl;
    
    os << "d_flow_model_label = "
       << d_flow_model_label
       << std::endl;
}


/*
 * Initialize boost::shared_ptr of the variables
 * in d_initial_conditions.
 */
void
FlowModelManager::setVariablesForInitialConditions()
{
    std::vector<boost::shared_ptr<pdat::CellVariable<double> > > conservative_variables =
        d_flow_model->getConservativeVariables();
    
    switch (d_flow_model_label)
    {
        case SINGLE_SPECIES:
        {
            d_initial_conditions->setVariablesForSingleSpecies(
                conservative_variables[0],
                conservative_variables[1],
                conservative_variables[2]);
            
            break;
        }
        case FOUR_EQN_CONSERVATIVE:
        {
            d_initial_conditions->setVariablesForFourEqnConservative(
                conservative_variables[0],
                conservative_variables[1],
                conservative_variables[2]);
            
            break;
        }
        case FIVE_EQN_ALLAIRE:
        {
            d_initial_conditions->setVariablesForFiveEqnAllaire(
                conservative_variables[0],
                conservative_variables[1],
                conservative_variables[2],
                conservative_variables[3]);
            
            break;
        }
    }
}
    