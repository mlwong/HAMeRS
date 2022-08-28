#include "flow/diffusive_flux_reconstructors/midpoint/DiffusiveFluxReconstructorMidpointSixthOrder.hpp"

#include "SAMRAI/geom/CartesianPatchGeometry.h"

#include <map>

DiffusiveFluxReconstructorMidpointSixthOrder::DiffusiveFluxReconstructorMidpointSixthOrder(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_eqn,
    const HAMERS_SHARED_PTR<FlowModel>& flow_model,
    const HAMERS_SHARED_PTR<tbox::Database>& diffusive_flux_reconstructor_db):
        DiffusiveFluxReconstructor(
            object_name,
            dim,
            grid_geometry,
            num_eqn,
            flow_model,
            diffusive_flux_reconstructor_db)
{
    d_num_diff_ghosts = hier::IntVector::getOne(d_dim)*11;
}


/*
 * Print all characteristics of the diffusive flux reconstruction class.
 */
void
DiffusiveFluxReconstructorMidpointSixthOrder::printClassData(
    std::ostream& os) const
{
    os << "\nPrint DiffusiveFluxReconstructorMidpointSixthOrder object..."
       << std::endl;
    
    os << std::endl;
    
    os << "DiffusiveFluxReconstructorMidpointSixthOrder: this = "
       << (DiffusiveFluxReconstructorMidpointSixthOrder *)this
       << std::endl;
    os << "d_object_name = "
       << d_object_name
       << std::endl;
}


/*
 * Put the characteristics of the diffusive flux reconstruction class
 * into the restart database.
 */
void
DiffusiveFluxReconstructorMidpointSixthOrder::putToRestart(
   const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const
{
    restart_db->putString("d_diffusive_flux_reconstructor", "MIDPOINT_SIXTH_ORDER");
}


/*
 * Compute the diffusive flux on a patch.
 */
void
DiffusiveFluxReconstructorMidpointSixthOrder::computeDiffusiveFluxOnPatch(
    hier::Patch& patch,
    const HAMERS_SHARED_PTR<pdat::SideVariable<double> >& variable_diffusive_flux,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double time,
    const double dt,
    const int RK_step_number)
{
    NULL_USE(time);
    NULL_USE(RK_step_number);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(variable_diffusive_flux);
#endif
    
    d_flow_model->setupDiffusiveFluxUtilities();
    
    HAMERS_SHARED_PTR<FlowModelDiffusiveFluxUtilities> diffusive_flux_utilities =
        d_flow_model->getFlowModelDiffusiveFluxUtilities();
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // diffusive ghost cells.
    hier::Box diff_ghost_box = interior_box;
    diff_ghost_box.grow(d_num_diff_ghosts);
    const hier::IntVector diff_ghostcell_dims = diff_ghost_box.numberCells();
    
    // Get the side data of diffusive flux.
    HAMERS_SHARED_PTR<pdat::SideData<double> > diffusive_flux(
        HAMERS_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
            patch.getPatchData(variable_diffusive_flux, data_context)));
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(diffusive_flux);
    TBOX_ASSERT(diffusive_flux->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
#endif
    
    // Initialize the data of diffusive flux to zero.
    diffusive_flux->fillAll(double(0));
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Register the patch and derived cell variables in the flow model and compute the corresponding cell data.
         */
        
        d_flow_model->registerPatchWithDataContext(patch, data_context);
        
        diffusive_flux_utilities->registerDerivedVariablesForDiffusiveFluxes(d_num_diff_ghosts);
        
        diffusive_flux_utilities->allocateMemoryForDerivedCellData();
        
        diffusive_flux_utilities->computeDerivedCellData();
        
        /*
         * Delcare containers for computing diffusivities at midpoints.
         */
        
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > var_cell_data_for_diffusivities;
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > var_side_data_for_diffusivities;
        std::vector<int> var_cell_data_for_diffusivities_component_idx;
        
        diffusive_flux_utilities->getCellDataForInterpolationToSideDataForDiffusiveFluxDiffusivities(
            var_cell_data_for_diffusivities,
            var_cell_data_for_diffusivities_component_idx);
        
        // Interpolate variables from nodes to midpoints for computing diffusivities at midpoints.
        
        // UNFINISHED...
        
        diffusive_flux_utilities->computeSideDataOfDiffusiveFluxDiffusivities(
            var_side_data_for_diffusivities);
        
        /*
         * Delcare containers for computing fluxes in different directions.
         */
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > var_data_x;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > > diffusivities_data_x;
        
        std::vector<std::vector<int> > var_component_idx_x;
        
        std::vector<std::vector<int> > diffusivities_component_idx_x;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > > derivative_x;
        
        std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivative_x_computed;
        
        HAMERS_SHARED_PTR<pdat::SideData<double> > diffusive_flux_midpoint;
        diffusive_flux_midpoint.reset(new pdat::SideData<double>(interior_box, d_num_eqn, d_num_diff_ghosts));
        
        std::vector<double*> F_midpoint_x;
        F_midpoint_x.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            F_midpoint_x.push_back(diffusive_flux_midpoint->getPointer(0, ei));
        }
        
        /*
         * (1) Compute the flux in the x-direction.
         */
        
        // Get the variables for the derivatives in the diffusive flux.
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxVariablesForDerivative(
            var_data_x,
            var_component_idx_x,
            DIRECTION::X_DIRECTION,
            DIRECTION::X_DIRECTION);
        
        // Get the diffusivities in the diffusive flux.
        diffusive_flux_utilities->getSideDataOfDiffusiveFluxDiffusivities(
            diffusivities_data_x,
            diffusivities_component_idx_x,
            DIRECTION::X_DIRECTION,
            DIRECTION::X_DIRECTION);
        
        TBOX_ASSERT(static_cast<int>(var_data_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_x.size()) == d_num_eqn);
        
        /*
         * Compute the derivatives in x-direction for diffusive flux in x-direction.
         */
        
        // computeFirstDerivativesInX(
        //     patch,
        //     derivative_x,
        //     derivative_x_computed,
        //     var_data_x,
        //     var_component_idx_x);
        
        /*
         * Compute diffusive flux in x-direction at nodes.
         */
        
        
        
        /*
         * Unregister the patch and data of all registered derived cell variables in the flow model.
         */
        
        d_flow_model->unregisterPatch();
        
    } // if (d_dim == tbox::Dimension(1))
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Register the patch and derived cell variables in the flow model and compute the corresponding cell data.
         */
        
        d_flow_model->registerPatchWithDataContext(patch, data_context);
        
        diffusive_flux_utilities->registerDerivedVariablesForDiffusiveFluxes(d_num_diff_ghosts);
        
        diffusive_flux_utilities->allocateMemoryForDerivedCellData();
        
        diffusive_flux_utilities->computeDerivedCellData();
        
        /*
         * Delcare containers for computing diffusivities at midpoints.
         */
        
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > var_cell_data_for_diffusivities;
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > var_side_data_for_diffusivities;
        std::vector<int> var_cell_data_for_diffusivities_component_idx;
        
        diffusive_flux_utilities->getCellDataForInterpolationToSideDataForDiffusiveFluxDiffusivities(
            var_cell_data_for_diffusivities,
            var_cell_data_for_diffusivities_component_idx);
        
        // Interpolate variables from nodes to midpoints for computing diffusivities at midpoints.
        
        // UNFINISHED...
        
        diffusive_flux_utilities->computeSideDataOfDiffusiveFluxDiffusivities(
            var_side_data_for_diffusivities);
        
        /*
         * Delcare containers for computing fluxes in different directions.
         */
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > var_data_x;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > var_data_y;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > > diffusivities_data_x;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > > diffusivities_data_y;
        
        std::vector<std::vector<int> > var_component_idx_x;
        std::vector<std::vector<int> > var_component_idx_y;
        
        std::vector<std::vector<int> > diffusivities_component_idx_x;
        std::vector<std::vector<int> > diffusivities_component_idx_y;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > > derivative_x_midpoint;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > > derivative_y_midpoint;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > derivative_x_node;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > derivative_y_node;
        
        std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivative_x_node_computed;
        std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivative_y_node_computed;
        
        HAMERS_SHARED_PTR<pdat::SideData<double> > diffusive_flux_midpoint;
        diffusive_flux_midpoint.reset(new pdat::SideData<double>(interior_box, d_num_eqn, d_num_diff_ghosts));
        
        std::vector<double*> F_midpoint_x;
        std::vector<double*> F_midpoint_y;
        F_midpoint_x.reserve(d_num_eqn);
        F_midpoint_y.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            F_midpoint_x.push_back(diffusive_flux_midpoint->getPointer(0, ei));
            F_midpoint_y.push_back(diffusive_flux_midpoint->getPointer(1, ei));
        }
        
        diffusive_flux_midpoint->fillAll(double(0));
        
        /*
         * (1) Compute the flux in the x-direction.
         */
        
        // Get the variables for the derivatives in the diffusive flux.
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxVariablesForDerivative(
            var_data_x,
            var_component_idx_x,
            DIRECTION::X_DIRECTION,
            DIRECTION::X_DIRECTION);
        
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxVariablesForDerivative(
            var_data_y,
            var_component_idx_y,
            DIRECTION::X_DIRECTION,
            DIRECTION::Y_DIRECTION);
        
        // Get the diffusivities in the diffusive flux.
        diffusive_flux_utilities->getSideDataOfDiffusiveFluxDiffusivities(
            diffusivities_data_x,
            diffusivities_component_idx_x,
            DIRECTION::X_DIRECTION,
            DIRECTION::X_DIRECTION);
        
        diffusive_flux_utilities->getSideDataOfDiffusiveFluxDiffusivities(
            diffusivities_data_y,
            diffusivities_component_idx_y,
            DIRECTION::X_DIRECTION,
            DIRECTION::Y_DIRECTION);
        
        TBOX_ASSERT(static_cast<int>(var_data_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_data_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_y.size()) == d_num_eqn);
        
        /*
         * Compute the derivatives in x-direction for diffusive flux in x-direction.
         */
        
        computeFirstDerivativesInXAtMidpoint(
            patch,
            derivative_x_midpoint,
            var_data_x,
            var_component_idx_x);
        
        /*
         * Compute the derivatives in y-direction for diffusive flux in x-direction.
         */
        
        computeFirstDerivativesInYAtNode(
            patch,
            derivative_y_node,
            derivative_y_node_computed,
            var_data_y,
            var_component_idx_y);
        
        interpolateFromNodeToMidpointX(
            patch,
            derivative_y_midpoint,
            derivative_y_node);
        
        /*
         * Compute diffusive flux in x-direction at midpoints.
         */
        
        /*
         * Reconstruct the flux in x-direction.
         */
        
        
        
        
        var_data_x.clear();
        var_data_y.clear();
        
        diffusivities_data_x.clear();
        diffusivities_data_y.clear();
        
        var_component_idx_x.clear();
        var_component_idx_y.clear();
        
        diffusivities_component_idx_x.clear();
        diffusivities_component_idx_y.clear();
        
        derivative_x_midpoint.clear();
        derivative_y_midpoint.clear();
        derivative_y_node.clear();
        
        /*
         * (2) Compute the flux in the y-direction.
         */
        
        // Get the variables for the derivatives in the diffusive flux.
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxVariablesForDerivative(
            var_data_x,
            var_component_idx_x,
            DIRECTION::Y_DIRECTION,
            DIRECTION::X_DIRECTION);
        
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxVariablesForDerivative(
            var_data_y,
            var_component_idx_y,
            DIRECTION::Y_DIRECTION,
            DIRECTION::Y_DIRECTION);
        
        // Get the diffusivities in the diffusive flux.
        diffusive_flux_utilities->getSideDataOfDiffusiveFluxDiffusivities(
            diffusivities_data_x,
            diffusivities_component_idx_x,
            DIRECTION::Y_DIRECTION,
            DIRECTION::X_DIRECTION);
        
        diffusive_flux_utilities->getSideDataOfDiffusiveFluxDiffusivities(
            diffusivities_data_y,
            diffusivities_component_idx_y,
            DIRECTION::Y_DIRECTION,
            DIRECTION::Y_DIRECTION);
        
        TBOX_ASSERT(static_cast<int>(var_data_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_data_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_y.size()) == d_num_eqn);
        
        /*
         * Compute the derivatives in x-direction for diffusive flux in y-direction.
         */
        
        computeFirstDerivativesInXAtNode(
            patch,
            derivative_x_node,
            derivative_x_node_computed,
            var_data_x,
            var_component_idx_x);
        
        interpolateFromNodeToMidpointY(
            patch,
            derivative_x_midpoint,
            derivative_x_node);
        
        /*
         * Compute the derivatives in y-direction for diffusive flux in y-direction.
         */
        
        computeFirstDerivativesInYAtMidpoint(
            patch,
            derivative_y_midpoint,
            var_data_y,
            var_component_idx_y);
        
        /*
         * Compute diffusive flux in y-direction at midpoints.
         */
        
        
        
        /*
         * Reconstruct the flux in y-direction.
         */
        
        
        
        var_data_x.clear();
        var_data_y.clear();
        
        diffusivities_data_x.clear();
        diffusivities_data_y.clear();
        
        var_component_idx_x.clear();
        var_component_idx_y.clear();
        
        diffusivities_component_idx_x.clear();
        diffusivities_component_idx_y.clear();
        
        derivative_x_midpoint.clear();
        derivative_y_midpoint.clear();
        derivative_x_node.clear();
        
        /*
         * Unregister the patch and data of all registered derived cell variables in the flow model.
         */
        
        d_flow_model->unregisterPatch();
        
    } // if (d_dim == tbox::Dimension(2))
    else if (d_dim == tbox::Dimension(3))
    {
        
    } // if (d_dim == tbox::Dimension(3))
}


/*
 * Compute the derivatives in x-direction at midpoints.
 */
void
DiffusiveFluxReconstructorMidpointSixthOrder::computeFirstDerivativesInXAtMidpoint(
    hier::Patch& patch,
    std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > >& derivative_x_midpoint,
    const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& data_x,
    const std::vector<std::vector<int> >& data_component_idx_x)
{
    
}


/*
 * Compute the derivatives in y-direction at midpoints.
 */
void
DiffusiveFluxReconstructorMidpointSixthOrder::computeFirstDerivativesInYAtMidpoint(
    hier::Patch& patch,
    std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > >& derivative_y_midpoint,
    const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& data_y,
    const std::vector<std::vector<int> >& data_component_idx_y)
{
    
}


/*
 * Compute the derivatives in z-direction at midpoints.
 */
void
DiffusiveFluxReconstructorMidpointSixthOrder::computeFirstDerivativesInZAtMidpoint(
    hier::Patch& patch,
    std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > >& derivative_z_midpoint,
    const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& data_z,
    const std::vector<std::vector<int> >& data_component_idx_z)
{
    
}


/*
 * Compute the derivatives in x-direction at nodes.
 */
void
DiffusiveFluxReconstructorMidpointSixthOrder::computeFirstDerivativesInXAtNode(
    hier::Patch& patch,
    std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& derivative_x_node,
    std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > >& derivative_x_node_computed,
    const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& var_data_x,
    const std::vector<std::vector<int> >& var_component_idx_x)
{
    
}


/*
 * Compute the derivatives in y-direction at nodes.
 */
void
DiffusiveFluxReconstructorMidpointSixthOrder::computeFirstDerivativesInYAtNode(
    hier::Patch& patch,
    std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& derivative_y_node,
    std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > >& derivative_y_node_computed,
    const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& var_data_y,
    const std::vector<std::vector<int> >& var_component_idx_y)
{
    
}


/*
 * Compute the derivatives in z-direction at nodes.
 */
void
DiffusiveFluxReconstructorMidpointSixthOrder::computeFirstDerivativesInZAtNode(
    hier::Patch& patch,
    std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& derivative_z_node,
    std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > >& derivative_z_node_computed,
    const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& var_data_z,
    const std::vector<std::vector<int> >& var_component_idx_z)
{
    
}


/*
 * Interpolate the variables from nodes to midpoints in x-direction.
 */
void
DiffusiveFluxReconstructorMidpointSixthOrder::interpolateFromNodeToMidpointX(
    hier::Patch& patch,
    std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > >& var_midpoint,
    const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& var_node)
{

}


/*
 * Interpolate the variables from nodes to midpoints in y-direction.
 */
void
DiffusiveFluxReconstructorMidpointSixthOrder::interpolateFromNodeToMidpointY(
    hier::Patch& patch,
    std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > >& var_midpoint,
    const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& var_node)
{

}


/*
 * Interpolate the variables from nodes to midpoints in z-direction.
 */
void
DiffusiveFluxReconstructorMidpointSixthOrder::interpolateFromNodeToMidpointZ(
    hier::Patch& patch,
    std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > >& var_midpoint,
    const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& var_node)
{

}
