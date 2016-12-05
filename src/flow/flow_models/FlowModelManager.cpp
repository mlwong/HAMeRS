#include "flow/flow_models/FlowModelManager.hpp"

FlowModelManager::FlowModelManager(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_species,
    const boost::shared_ptr<tbox::Database>& flow_model_db,
    const std::string& flow_model_str):
        d_object_name(object_name),
        d_dim(dim),
        d_grid_geometry(grid_geometry),
        d_num_species(num_species)
{
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(grid_geometry);
    
    if (flow_model_str == "SINGLE_SPECIES")
    {
        d_flow_model_type = FLOW_MODEL::SINGLE_SPECIES;
        
        d_flow_model.reset(new FlowModelSingleSpecies(
            "d_flow_model",
            d_dim,
            d_grid_geometry,
            hier::IntVector::getZero(d_dim),
            d_num_species,
            flow_model_db));
    }
    else if (flow_model_str == "FOUR_EQN_CONSERVATIVE")
    {
        d_flow_model_type = FLOW_MODEL::FOUR_EQN_CONSERVATIVE;
        
        d_flow_model.reset(new FlowModelFourEqnConservative(
            "d_flow_model",
            d_dim,
            d_grid_geometry,
            hier::IntVector::getZero(d_dim),
            d_num_species,
            flow_model_db));
    }
    else if (flow_model_str == "FIVE_EQN_ALLAIRE")
    {
        d_flow_model_type = FLOW_MODEL::FIVE_EQN_ALLAIRE;
        
        d_flow_model.reset(new FlowModelFiveEqnAllaire(
            "d_flow_model",
            d_dim,
            d_grid_geometry,
            hier::IntVector::getZero(d_dim),
            d_num_species,
            flow_model_db));
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
    
    if (d_num_species > 1 && d_flow_model_type == FLOW_MODEL::SINGLE_SPECIES)
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
    
    os << "================================================================================";
    d_flow_model->printClassData(os);
}
