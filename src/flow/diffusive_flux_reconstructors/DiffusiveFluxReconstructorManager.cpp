#include "flow/diffusive_flux_reconstructors/DiffusiveFluxReconstructorManager.hpp"

DiffusiveFluxReconstructorManager::DiffusiveFluxReconstructorManager(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_eqn,
    const boost::shared_ptr<FlowModel>& flow_model,
    const boost::shared_ptr<tbox::Database>& diffusive_flux_reconstructor_db,
    const std::string& diffusive_flux_reconstructor_str):
        d_object_name(object_name)
{
    if (diffusive_flux_reconstructor_str == "SIXTH_ORDER")
    {
        d_diffusive_flux_reconstructor_type = DIFFUSIVE_FLUX_RECONSTRUCTOR::SIXTH_ORDER;
        
        d_diff_flux_reconstructor.reset(new DiffusiveFluxReconstructorSixthOrder(
            "d_diffusive_flux_reconstructor",
            dim,
            grid_geometry,
            flow_model->getNumberOfEquations(),
            flow_model,
            diffusive_flux_reconstructor_db));
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
