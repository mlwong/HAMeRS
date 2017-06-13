#include "flow/diffusive_flux_reconstructors/DiffusiveFluxReconstructorManager.hpp"

DiffusiveFluxReconstructorManager::DiffusiveFluxReconstructorManager(
    const std::string& object_name,
    const boost::shared_ptr<tbox::Database>& diffusive_flux_reconstructor_db,
    const std::string& diffusive_flux_reconstructor_str):
        d_object_name(object_name),
        d_diffusive_flux_reconstructor_db(diffusive_flux_reconstructor_db)
{
    if (diffusive_flux_reconstructor_str == "SIXTH_ORDER")
    {
        d_diffusive_flux_reconstructor_type = DIFFUSIVE_FLUX_RECONSTRUCTOR::SIXTH_ORDER;
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Unknown diffusive_flux_reconstructor string = '"
            << diffusive_flux_reconstructor_str
            << "' found in input."
            << std::endl);        
    }
}


/*
 * Create the diffusive flux reconstructor.
 */
boost::shared_ptr<DiffusiveFluxReconstructor>
DiffusiveFluxReconstructorManager::createDiffusiveFluxReconstructor(
    const tbox::Dimension& dim,
    const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_species,
    const boost::shared_ptr<FlowModel>& flow_model)
{
    boost::shared_ptr<DiffusiveFluxReconstructor> diff_flux_reconstructor;
    
    switch (d_diffusive_flux_reconstructor_type)
    {
        case (DIFFUSIVE_FLUX_RECONSTRUCTOR::SIXTH_ORDER):
        {
            diff_flux_reconstructor.reset(new DiffusiveFluxReconstructorSixthOrder(
                "d_diffusive_flux_reconstructor",
                dim,
                grid_geometry,
                flow_model->getNumberOfEquations(),
                num_species,
                flow_model,
                d_diffusive_flux_reconstructor_db));
            
            break;
        }
    }
    
    return diff_flux_reconstructor;
}


/*
 * Print all characteristics of diffusive flux reconstructor manager.
 */
void
DiffusiveFluxReconstructorManager::printClassData(std::ostream& os) const
{
    os << "\nPrint DiffusiveFluxReconstructorManager object..."
       << std::endl;
    
    os << std::endl;
    
    os << "DiffusiveFluxReconstructorManager: this = "
       << (DiffusiveFluxReconstructorManager *)this
       << std::endl;
    
    os << "d_object_name = "
       << d_object_name
       << std::endl;
    
    os << "d_diffusive_flux_reconstructor_type = "
       << d_diffusive_flux_reconstructor_type
       << std::endl;
}
