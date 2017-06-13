#include "flow/flow_models/FlowModelManager.hpp"

FlowModelManager::FlowModelManager(
    const std::string& object_name,
    const boost::shared_ptr<tbox::Database>& flow_model_db,
    const std::string& flow_model_str):
        d_object_name(object_name),
        d_flow_model_db(flow_model_db)
{
    TBOX_ASSERT(!object_name.empty());
    
    if (flow_model_str == "SINGLE_SPECIES")
    {
        d_flow_model_type = FLOW_MODEL::SINGLE_SPECIES;
    }
    else if (flow_model_str == "FOUR_EQN_CONSERVATIVE")
    {
        d_flow_model_type = FLOW_MODEL::FOUR_EQN_CONSERVATIVE;
    }
    else if (flow_model_str == "FIVE_EQN_ALLAIRE")
    {
        d_flow_model_type = FLOW_MODEL::FIVE_EQN_ALLAIRE;
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
}


/*
 * Create the flow model.
 */
boost::shared_ptr<FlowModel>
FlowModelManager::createFlowModel(
    const tbox::Dimension& dim,
    const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_species)
{
    if (num_species > 1 && d_flow_model_type == FLOW_MODEL::SINGLE_SPECIES)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Number of species = "
            << num_species
            << " shouldn't use single-species model."
            << std::endl); 
    }
    
    boost::shared_ptr<FlowModel> flow_model;
    
    switch (d_flow_model_type)
    {
        case FLOW_MODEL::SINGLE_SPECIES:
        {
            flow_model.reset(new FlowModelSingleSpecies(
                "d_flow_model",
                dim,
                grid_geometry,
                num_species,
                d_flow_model_db));
            
            break;
        }
        case FLOW_MODEL::FOUR_EQN_CONSERVATIVE:
        {
            flow_model.reset(new FlowModelFourEqnConservative(
                "d_flow_model",
                dim,
                grid_geometry,
                num_species,
                d_flow_model_db));
            
            break;
        }
        case FLOW_MODEL::FIVE_EQN_ALLAIRE:
        {
            flow_model.reset(new FlowModelFiveEqnAllaire(
                "d_flow_model",
                dim,
                grid_geometry,
                num_species,
                d_flow_model_db));
            
            break;
        }
        
    }
    
    return flow_model;
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
    
    os << "d_flow_model_type = "
       << d_flow_model_type
       << std::endl;
}
