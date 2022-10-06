#include "flow/diffusive_flux_reconstructors/DiffusiveFluxReconstructorManager.hpp"

DiffusiveFluxReconstructorManager::DiffusiveFluxReconstructorManager(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_eqn,
    const HAMERS_SHARED_PTR<FlowModel>& flow_model,
    const HAMERS_SHARED_PTR<tbox::Database>& diffusive_flux_reconstructor_db,
    const std::string& diffusive_flux_reconstructor_str):
        d_object_name(object_name)
{
    if (diffusive_flux_reconstructor_str == "MIDPOINT_SECOND_ORDER")
    {
        d_diffusive_flux_reconstructor_type = DIFFUSIVE_FLUX_RECONSTRUCTOR::MIDPOINT_SECOND_ORDER;
        
        d_diff_flux_reconstructor.reset(new DiffusiveFluxReconstructorMidpointSecondOrder(
            "d_diffusive_flux_reconstructor",
            dim,
            grid_geometry,
            flow_model->getNumberOfEquations(),
            flow_model,
            diffusive_flux_reconstructor_db));
    }
    else if (diffusive_flux_reconstructor_str == "MIDPOINT_FOURTH_ORDER")
    {
        d_diffusive_flux_reconstructor_type = DIFFUSIVE_FLUX_RECONSTRUCTOR::MIDPOINT_FOURTH_ORDER;
        
        d_diff_flux_reconstructor.reset(new DiffusiveFluxReconstructorMidpointFourthOrder(
            "d_diffusive_flux_reconstructor",
            dim,
            grid_geometry,
            flow_model->getNumberOfEquations(),
            flow_model,
            diffusive_flux_reconstructor_db));
    }
    else if (diffusive_flux_reconstructor_str == "MIDPOINT_SIXTH_ORDER")
    {
        d_diffusive_flux_reconstructor_type = DIFFUSIVE_FLUX_RECONSTRUCTOR::MIDPOINT_SIXTH_ORDER;
        
        d_diff_flux_reconstructor.reset(new DiffusiveFluxReconstructorMidpointSixthOrder(
            "d_diffusive_flux_reconstructor",
            dim,
            grid_geometry,
            flow_model->getNumberOfEquations(),
            flow_model,
            diffusive_flux_reconstructor_db));
    }
    else if (diffusive_flux_reconstructor_str == "SECOND_ORDER")
    {
        d_diffusive_flux_reconstructor_type = DIFFUSIVE_FLUX_RECONSTRUCTOR::SECOND_ORDER;
        
        d_diff_flux_reconstructor.reset(new DiffusiveFluxReconstructorNodeSecondOrder(
            "d_diffusive_flux_reconstructor",
            dim,
            grid_geometry,
            flow_model->getNumberOfEquations(),
            flow_model,
            diffusive_flux_reconstructor_db));
    }
    else if (diffusive_flux_reconstructor_str == "FOURTH_ORDER")
    {
        d_diffusive_flux_reconstructor_type = DIFFUSIVE_FLUX_RECONSTRUCTOR::FOURTH_ORDER;
        
        d_diff_flux_reconstructor.reset(new DiffusiveFluxReconstructorNodeFourthOrder(
            "d_diffusive_flux_reconstructor",
            dim,
            grid_geometry,
            flow_model->getNumberOfEquations(),
            flow_model,
            diffusive_flux_reconstructor_db));
    }
    else if (diffusive_flux_reconstructor_str == "SIXTH_ORDER")
    {
        d_diffusive_flux_reconstructor_type = DIFFUSIVE_FLUX_RECONSTRUCTOR::SIXTH_ORDER;
        
        d_diff_flux_reconstructor.reset(new DiffusiveFluxReconstructorNodeSixthOrder(
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
