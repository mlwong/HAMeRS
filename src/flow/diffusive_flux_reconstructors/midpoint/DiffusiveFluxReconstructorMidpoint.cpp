#include "flow/diffusive_flux_reconstructors/midpoint/DiffusiveFluxReconstructorMidpoint.hpp"

#include "SAMRAI/geom/CartesianPatchGeometry.h"

#include <map>

DiffusiveFluxReconstructorMidpoint::DiffusiveFluxReconstructorMidpoint(
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
            diffusive_flux_reconstructor_db),
        d_num_der_midpoint_ghosts(hier::IntVector::getZero(d_dim)),
        d_num_der_node_ghosts(hier::IntVector::getZero(d_dim)),
        d_num_interp_midpoint_ghosts(hier::IntVector::getZero(d_dim)),
        d_num_flux_reconstruct_ghosts(hier::IntVector::getZero(d_dim)),
        d_num_scratch_derivatives_node_used(0),
        d_num_scratch_derivatives_midpoint_x_used(0),
        d_num_scratch_derivatives_midpoint_y_used(0),
        d_num_scratch_derivatives_midpoint_z_used(0),
        d_num_scratch_diffusivities_midpoint_used(0)
{
}


/*
 * Compute the diffusive flux on a patch.
 */
void
DiffusiveFluxReconstructorMidpoint::computeDiffusiveFluxOnPatch(
    hier::Patch& patch,
    const HAMERS_SHARED_PTR<hier::CoarseFineBoundary> coarse_fine_bdry,
    const HAMERS_SHARED_PTR<pdat::SideVariable<Real> >& variable_diffusive_flux,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double time,
    const double dt,
    const int RK_step_number)
{
    NULL_USE(coarse_fine_bdry);
    NULL_USE(time);
    NULL_USE(RK_step_number);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(variable_diffusive_flux);
#endif
    
    const bool allocate_scratch_data_containers = true;
    
    d_flow_model->setupDiffusiveFluxUtilities();
    
    HAMERS_SHARED_PTR<FlowModelDiffusiveFluxUtilities> diffusive_flux_utilities =
        d_flow_model->getFlowModelDiffusiveFluxUtilities();
    
    if (d_dim == tbox::Dimension(1) || d_dim == tbox::Dimension(2))
    {
        if (diffusive_flux_utilities->useSubgridScaleModel())
        {
            TBOX_ERROR(d_object_name
                << ": DiffusiveFluxReconstructorMidpoint::"
                << "computeDiffusiveFluxOnPatch()\n"
                << "Subgrid scale model is not implemented for one-dimensional or two-dimensional problems."
                << std::endl);
        }
    }
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // diffusive ghost cells.
    hier::Box diff_ghost_box = interior_box;
    diff_ghost_box.grow(d_num_diff_ghosts);
    const hier::IntVector diff_ghostcell_dims = diff_ghost_box.numberCells();
    
    // Get the side data of diffusive flux.
    HAMERS_SHARED_PTR<pdat::SideData<Real> > diffusive_flux(
        HAMERS_SHARED_PTR_CAST<pdat::SideData<Real>, hier::PatchData>(
            patch.getPatchData(variable_diffusive_flux, data_context)));
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(diffusive_flux);
    TBOX_ASSERT(diffusive_flux->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
#endif
    
    // Initialize the data of diffusive flux to zero.
    diffusive_flux->fillAll(Real(0));
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the dimension and number of ghost cells.
         */
        
        const int interior_dim_0 = interior_dims[0];
        
        const int num_diff_ghosts_0 = d_num_diff_ghosts[0];
        
        const int num_flux_reconstruct_ghosts_0 = d_num_flux_reconstruct_ghosts[0];
        
        /*
         * Register the patch and derived cell variables in the flow model and compute the corresponding cell data.
         */
        
        d_flow_model->registerPatchWithDataContext(patch, data_context);
        
        diffusive_flux_utilities->registerDerivedVariablesForDiffusiveFluxes(d_num_diff_ghosts, true);
        
        diffusive_flux_utilities->allocateMemoryForDerivedCellData();
        
        diffusive_flux_utilities->computeDerivedCellData();
        
        /*
         * Delcare containers for computing diffusivities at midpoints.
         */
        
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > var_cell_data_for_diffusivities;
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > var_side_data_for_diffusivities;
        std::vector<int> var_cell_data_for_diffusivities_component_idx;
        
        diffusive_flux_utilities->getCellDataForInterpolationToSideDataForDiffusiveFluxDiffusivities(
            var_cell_data_for_diffusivities,
            var_cell_data_for_diffusivities_component_idx);
        
        // Interpolate variables from nodes to midpoints for computing diffusivities at midpoints.
        
        interpolateDiffusivitiesFromNodeToMidpoint(
            var_side_data_for_diffusivities,
            var_cell_data_for_diffusivities,
            var_cell_data_for_diffusivities_component_idx,
            patch,
            allocate_scratch_data_containers);
        
        diffusive_flux_utilities->allocateMemoryForSideDataOfDiffusiveFluxDiffusivities();
        
        diffusive_flux_utilities->computeSideDataOfDiffusiveFluxDiffusivities(
            var_side_data_for_diffusivities);
        
        /*
         * Delcare containers for computing fluxes in different directions.
         */
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > > var_data_x;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > > diffusivities_data_x;
        
        std::vector<std::vector<int> > var_component_idx_x;
        
        std::vector<std::vector<int> > diffusivities_component_idx_x;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > > derivatives_x_midpoint_x;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > > derivatives_x_node;
        
        std::map<Real*, HAMERS_SHARED_PTR<pdat::SideData<Real> > > derivatives_x_midpoint_x_computed;

        std::map<Real*, HAMERS_SHARED_PTR<pdat::CellData<Real> > > derivatives_x_node_computed;
        
        HAMERS_SHARED_PTR<pdat::SideData<Real> > diffusive_flux_midpoint;
        diffusive_flux_midpoint.reset(new pdat::SideData<Real>(interior_box, d_num_eqn, d_num_diff_ghosts));
        
        std::vector<Real*> F_midpoint_x;
        F_midpoint_x.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            F_midpoint_x.push_back(diffusive_flux_midpoint->getPointer(0, ei));
        }
        
        diffusive_flux_midpoint->fillAll(Real(0));
        
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
        
        computeFirstDerivativesInXAtMidpointX(
            derivatives_x_midpoint_x,
            derivatives_x_midpoint_x_computed,
            var_data_x,
            var_component_idx_x,
            patch,
            allocate_scratch_data_containers);
        
        /*
         * Compute diffusive flux in x-direction at midpoints.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ASSERT(static_cast<int>(derivatives_x_midpoint_x[ei].size()) ==
                        static_cast<int>(diffusivities_data_x[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_x[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_x[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_x[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_x[ei][vi];
                
                // Get the pointer to diffusivity.
                Real* mu = diffusivities_data_x[ei][vi]->getPointer(0, mu_idx);
                
                // Get the pointer to derivative.
                Real* dudx = derivatives_x_midpoint_x[ei][vi]->getPointer(0, 0);
                
                /*
                 * Get the sub-ghost cell width of the diffusivity.
                 */
                
                hier::IntVector num_subghosts_diffusivity =
                    diffusivities_data_x[ei][vi]->getGhostCellWidth();
                
                const int num_subghosts_0_diffusivity = num_subghosts_diffusivity[0];
                
                HAMERS_PRAGMA_SIMD
                for (int i = -num_subghosts_0_diffusivity + num_flux_reconstruct_ghosts_0;
                    i < (interior_dim_0 + 1) + num_subghosts_0_diffusivity - num_flux_reconstruct_ghosts_0;
                    i++)
                {
                    // Compute the linear indices.
                    const int idx = i + num_diff_ghosts_0;
                    const int idx_diffusivity = i + num_subghosts_0_diffusivity;
                    
                    F_midpoint_x[ei][idx] += mu[idx_diffusivity]*dudx[idx];
                }
            }
        }
        
        /*
         * Reconstruct the flux in x-direction.
         */
        
        reconstructFluxX(
            diffusive_flux,
            diffusive_flux_midpoint,
            patch,
            dt);
        
        var_data_x.clear();
        diffusivities_data_x.clear();
        var_component_idx_x.clear();
        diffusivities_component_idx_x.clear();
        derivatives_x_midpoint_x.clear();
        
        /*
         * Unregister the patch and data of all registered derived cell variables in the flow model.
         */
        
        d_flow_model->unregisterPatch();
        
    } // if (d_dim == tbox::Dimension(1))
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the dimensions and number of ghost cells.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        const int num_diff_ghosts_0 = d_num_diff_ghosts[0];
        const int num_diff_ghosts_1 = d_num_diff_ghosts[1];
        const int diff_ghostcell_dim_0 = diff_ghostcell_dims[0];
        
        const int num_flux_reconstruct_ghosts_0 = d_num_flux_reconstruct_ghosts[0];
        const int num_flux_reconstruct_ghosts_1 = d_num_flux_reconstruct_ghosts[1];
        
        /*
         * Register the patch and derived cell variables in the flow model and compute the corresponding cell data.
         */
        
        d_flow_model->registerPatchWithDataContext(patch, data_context);
        
        diffusive_flux_utilities->registerDerivedVariablesForDiffusiveFluxes(d_num_diff_ghosts, true);
        
        diffusive_flux_utilities->allocateMemoryForDerivedCellData();
        
        diffusive_flux_utilities->computeDerivedCellData();
        
        /*
         * Delcare containers for computing diffusivities at midpoints.
         */
        
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > var_cell_data_for_diffusivities;
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > var_side_data_for_diffusivities;
        std::vector<int> var_cell_data_for_diffusivities_component_idx;
        
        diffusive_flux_utilities->getCellDataForInterpolationToSideDataForDiffusiveFluxDiffusivities(
            var_cell_data_for_diffusivities,
            var_cell_data_for_diffusivities_component_idx);
        
        // Interpolate variables from nodes to midpoints for computing diffusivities at midpoints.
        
        interpolateDiffusivitiesFromNodeToMidpoint(
            var_side_data_for_diffusivities,
            var_cell_data_for_diffusivities,
            var_cell_data_for_diffusivities_component_idx,
            patch,
            allocate_scratch_data_containers);
        
        diffusive_flux_utilities->allocateMemoryForSideDataOfDiffusiveFluxDiffusivities();
        
        diffusive_flux_utilities->computeSideDataOfDiffusiveFluxDiffusivities(
            var_side_data_for_diffusivities);
        
        /*
         * Delcare containers for computing fluxes in different directions.
         */
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > > var_data_x;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > > var_data_y;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > > diffusivities_data_x;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > > diffusivities_data_y;
        
        std::vector<std::vector<int> > var_component_idx_x;
        std::vector<std::vector<int> > var_component_idx_y;
        
        std::vector<std::vector<int> > diffusivities_component_idx_x;
        std::vector<std::vector<int> > diffusivities_component_idx_y;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > > derivatives_x_midpoint_x;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > > derivatives_y_midpoint_x;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > > derivatives_x_midpoint_y;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > > derivatives_y_midpoint_y;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > > derivatives_x_node;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > > derivatives_y_node;
        
        std::map<Real*, HAMERS_SHARED_PTR<pdat::SideData<Real> > > derivatives_x_midpoint_x_computed;
        std::map<Real*, HAMERS_SHARED_PTR<pdat::SideData<Real> > > derivatives_y_midpoint_x_computed;
        
        std::map<Real*, HAMERS_SHARED_PTR<pdat::SideData<Real> > > derivatives_x_midpoint_y_computed;
        std::map<Real*, HAMERS_SHARED_PTR<pdat::SideData<Real> > > derivatives_y_midpoint_y_computed;
        
        std::map<Real*, HAMERS_SHARED_PTR<pdat::CellData<Real> > > derivatives_x_node_computed;
        std::map<Real*, HAMERS_SHARED_PTR<pdat::CellData<Real> > > derivatives_y_node_computed;
        
        HAMERS_SHARED_PTR<pdat::SideData<Real> > diffusive_flux_midpoint;
        diffusive_flux_midpoint.reset(new pdat::SideData<Real>(interior_box, d_num_eqn, d_num_diff_ghosts));
        
        std::vector<Real*> F_midpoint_x;
        std::vector<Real*> F_midpoint_y;
        F_midpoint_x.reserve(d_num_eqn);
        F_midpoint_y.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            F_midpoint_x.push_back(diffusive_flux_midpoint->getPointer(0, ei));
            F_midpoint_y.push_back(diffusive_flux_midpoint->getPointer(1, ei));
        }
        
        diffusive_flux_midpoint->fillAll(Real(0));
        
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
        
        computeFirstDerivativesInXAtMidpointX(
            derivatives_x_midpoint_x,
            derivatives_x_midpoint_x_computed,
            var_data_x,
            var_component_idx_x,
            patch,
            allocate_scratch_data_containers);
        
        /*
         * Compute the derivatives in y-direction for diffusive flux in x-direction.
         */
        
        computeFirstDerivativesInYAtNode(
            derivatives_y_node,
            derivatives_y_node_computed,
            var_data_y,
            var_component_idx_y,
            patch,
            allocate_scratch_data_containers);
        
        interpolateDerivativesFromNodeToMidpointX(
            derivatives_y_midpoint_x,
            derivatives_y_midpoint_x_computed,
            derivatives_y_node,
            var_data_y,
            var_component_idx_y,
            patch,
            allocate_scratch_data_containers);
        
        /*
         * Compute diffusive flux in x-direction at midpoints.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ASSERT(static_cast<int>(derivatives_x_midpoint_x[ei].size()) ==
                        static_cast<int>(diffusivities_data_x[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_x[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_x[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_x[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_x[ei][vi];
                
                // Get the pointer to diffusivity.
                Real* mu = diffusivities_data_x[ei][vi]->getPointer(0, mu_idx);
                
                // Get the pointer to derivative.
                Real* dudx = derivatives_x_midpoint_x[ei][vi]->getPointer(0, 0);
                
                /*
                 * Get the sub-ghost cell width and ghost box dimensions of the diffusivity.
                 */
                
                hier::IntVector num_subghosts_diffusivity =
                    diffusivities_data_x[ei][vi]->getGhostCellWidth();
                
                hier::IntVector subghostcell_dims_diffusivity =
                    diffusivities_data_x[ei][vi]->getGhostBox().numberCells();
                
                const int num_subghosts_0_diffusivity = num_subghosts_diffusivity[0];
                const int num_subghosts_1_diffusivity = num_subghosts_diffusivity[1];
                const int subghostcell_dim_0_diffusivity = subghostcell_dims_diffusivity[0];
                
                for (int j = 0; j < interior_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = -num_subghosts_0_diffusivity + num_flux_reconstruct_ghosts_0;
                        i < (interior_dim_0 + 1) + num_subghosts_0_diffusivity - num_flux_reconstruct_ghosts_0;
                        i++)
                    {
                        // Compute the linear indices.
                        const int idx = (i + num_diff_ghosts_0) +
                            (j + num_diff_ghosts_1)*(diff_ghostcell_dim_0 + 1);
                        
                        const int idx_diffusivity = (i + num_subghosts_0_diffusivity) +
                            (j + num_subghosts_1_diffusivity)*(subghostcell_dim_0_diffusivity + 1);
                        
                        F_midpoint_x[ei][idx] += mu[idx_diffusivity]*dudx[idx];
                    }
                }
            }
            
            TBOX_ASSERT(static_cast<int>(derivatives_y_midpoint_x[ei].size()) ==
                        static_cast<int>(diffusivities_data_y[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_y[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_y[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_y[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_y[ei][vi];
                
                // Get the pointer to diffusivity.
                Real* mu = diffusivities_data_y[ei][vi]->getPointer(0, mu_idx);
                
                // Get the pointer to derivative.
                Real* dudy = derivatives_y_midpoint_x[ei][vi]->getPointer(0, 0);
                
                /*
                 * Get the sub-ghost cell width and ghost box dimensions of the diffusivity.
                 */
                
                hier::IntVector num_subghosts_diffusivity =
                    diffusivities_data_y[ei][vi]->getGhostCellWidth();
                
                hier::IntVector subghostcell_dims_diffusivity =
                    diffusivities_data_y[ei][vi]->getGhostBox().numberCells();
                
                const int num_subghosts_0_diffusivity = num_subghosts_diffusivity[0];
                const int num_subghosts_1_diffusivity = num_subghosts_diffusivity[1];
                const int subghostcell_dim_0_diffusivity = subghostcell_dims_diffusivity[0];
                
                for (int j = 0; j < interior_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = -num_subghosts_0_diffusivity + num_flux_reconstruct_ghosts_0;
                        i < (interior_dim_0 + 1) + num_subghosts_0_diffusivity - num_flux_reconstruct_ghosts_0;
                        i++)
                    {
                        // Compute the linear indices.
                        const int idx = (i + num_diff_ghosts_0) +
                            (j + num_diff_ghosts_1)*(diff_ghostcell_dim_0 + 1);
                        
                        const int idx_diffusivity = (i + num_subghosts_0_diffusivity) +
                            (j + num_subghosts_1_diffusivity)*(subghostcell_dim_0_diffusivity + 1);
                        
                        F_midpoint_x[ei][idx] += mu[idx_diffusivity]*dudy[idx];
                    }
                }
            }
        }
        
        /*
         * Reconstruct the flux in x-direction.
         */
        
        reconstructFluxX(
            diffusive_flux,
            diffusive_flux_midpoint,
            patch,
            dt);
        
        var_data_x.clear();
        var_data_y.clear();
        
        diffusivities_data_x.clear();
        diffusivities_data_y.clear();
        
        var_component_idx_x.clear();
        var_component_idx_y.clear();
        
        diffusivities_component_idx_x.clear();
        diffusivities_component_idx_y.clear();
        
        derivatives_x_midpoint_x.clear();
        derivatives_y_midpoint_x.clear();
        
        derivatives_y_node.clear();
        
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
            derivatives_x_node,
            derivatives_x_node_computed,
            var_data_x,
            var_component_idx_x,
            patch,
            allocate_scratch_data_containers);
        
        interpolateDerivativesFromNodeToMidpointY(
            derivatives_x_midpoint_y,
            derivatives_x_midpoint_y_computed,
            derivatives_x_node,
            var_data_x,
            var_component_idx_x,
            patch,
            allocate_scratch_data_containers);
        
        /*
         * Compute the derivatives in y-direction for diffusive flux in y-direction.
         */
        
        computeFirstDerivativesInYAtMidpointY(
            derivatives_y_midpoint_y,
            derivatives_y_midpoint_y_computed,
            var_data_y,
            var_component_idx_y,
            patch,
            allocate_scratch_data_containers);
        
        /*
         * Compute diffusive flux in y-direction at midpoints.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ASSERT(static_cast<int>(derivatives_x_midpoint_y[ei].size()) ==
                        static_cast<int>(diffusivities_data_x[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_x[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_x[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_x[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_x[ei][vi];
                
                // Get the pointer to diffusivity.
                Real* mu = diffusivities_data_x[ei][vi]->getPointer(1, mu_idx);
                
                // Get the pointer to derivative.
                Real* dudx = derivatives_x_midpoint_y[ei][vi]->getPointer(1, 0);
                
                /*
                 * Get the sub-ghost cell width and ghost box dimensions of the diffusivity.
                 */
                
                hier::IntVector num_subghosts_diffusivity =
                    diffusivities_data_x[ei][vi]->getGhostCellWidth();
                
                hier::IntVector subghostcell_dims_diffusivity =
                    diffusivities_data_x[ei][vi]->getGhostBox().numberCells();
                
                const int num_subghosts_0_diffusivity = num_subghosts_diffusivity[0];
                const int num_subghosts_1_diffusivity = num_subghosts_diffusivity[1];
                const int subghostcell_dim_0_diffusivity = subghostcell_dims_diffusivity[0];
                
                for (int j = -num_subghosts_1_diffusivity + num_flux_reconstruct_ghosts_1;
                    j < (interior_dim_1 + 1) + num_subghosts_1_diffusivity - num_flux_reconstruct_ghosts_1;
                    j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx = (i + num_diff_ghosts_0) +
                            (j + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                        
                        const int idx_diffusivity = (i + num_subghosts_0_diffusivity) +
                            (j + num_subghosts_1_diffusivity)*subghostcell_dim_0_diffusivity;
                        
                        F_midpoint_y[ei][idx] += mu[idx_diffusivity]*dudx[idx];
                    }
                }
            }
            
            TBOX_ASSERT(static_cast<int>(derivatives_y_midpoint_y[ei].size()) ==
                        static_cast<int>(diffusivities_data_y[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_y[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_y[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_y[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_y[ei][vi];
                
                // Get the pointer to diffusivity.
                Real* mu = diffusivities_data_y[ei][vi]->getPointer(1, mu_idx);
                
                // Get the pointer to derivative.
                Real* dudy = derivatives_y_midpoint_y[ei][vi]->getPointer(1, 0);
                
                /*
                 * Get the sub-ghost cell width and ghost box dimensions of the diffusivity.
                 */
                
                hier::IntVector num_subghosts_diffusivity =
                    diffusivities_data_y[ei][vi]->getGhostCellWidth();
                
                hier::IntVector subghostcell_dims_diffusivity =
                    diffusivities_data_y[ei][vi]->getGhostBox().numberCells();
                
                const int num_subghosts_0_diffusivity = num_subghosts_diffusivity[0];
                const int num_subghosts_1_diffusivity = num_subghosts_diffusivity[1];
                const int subghostcell_dim_0_diffusivity = subghostcell_dims_diffusivity[0];
                
                for (int j = -num_subghosts_1_diffusivity + num_flux_reconstruct_ghosts_1;
                    j < (interior_dim_1 + 1) + num_subghosts_1_diffusivity - num_flux_reconstruct_ghosts_1;
                    j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx = (i + num_diff_ghosts_0) +
                            (j + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                        
                        const int idx_diffusivity = (i + num_subghosts_0_diffusivity) +
                            (j + num_subghosts_1_diffusivity)*subghostcell_dim_0_diffusivity;
                        
                        F_midpoint_y[ei][idx] += mu[idx_diffusivity]*dudy[idx];
                    }
                }
            }
        }
        
        /*
         * Reconstruct the flux in y-direction.
         */
        
        reconstructFluxY(
            diffusive_flux,
            diffusive_flux_midpoint,
            patch,
            dt);
        
        var_data_x.clear();
        var_data_y.clear();
        
        diffusivities_data_x.clear();
        diffusivities_data_y.clear();
        
        var_component_idx_x.clear();
        var_component_idx_y.clear();
        
        diffusivities_component_idx_x.clear();
        diffusivities_component_idx_y.clear();
        
        derivatives_x_midpoint_y.clear();
        derivatives_y_midpoint_y.clear();
        
        derivatives_x_node.clear();
        
        /*
         * Unregister the patch and data of all registered derived cell variables in the flow model.
         */
        
        d_flow_model->unregisterPatch();
        
    } // if (d_dim == tbox::Dimension(2))
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the dimensions and number of ghost cells.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        const int interior_dim_2 = interior_dims[2];
        
        const int num_diff_ghosts_0 = d_num_diff_ghosts[0];
        const int num_diff_ghosts_1 = d_num_diff_ghosts[1];
        const int num_diff_ghosts_2 = d_num_diff_ghosts[2];
        const int diff_ghostcell_dim_0 = diff_ghostcell_dims[0];
        const int diff_ghostcell_dim_1 = diff_ghostcell_dims[1];
        
        const int num_flux_reconstruct_ghosts_0 = d_num_flux_reconstruct_ghosts[0];
        const int num_flux_reconstruct_ghosts_1 = d_num_flux_reconstruct_ghosts[1];
        const int num_flux_reconstruct_ghosts_2 = d_num_flux_reconstruct_ghosts[2];
        
        /*
         * Register the patch and derived cell variables in the flow model and compute the corresponding cell data.
         */
        
        d_flow_model->registerPatchWithDataContext(patch, data_context);
        
        diffusive_flux_utilities->registerDerivedVariablesForDiffusiveFluxes(d_num_diff_ghosts, true);
        
        diffusive_flux_utilities->allocateMemoryForDerivedCellData();
        
        diffusive_flux_utilities->computeDerivedCellData();
        
        /*
         * Delcare containers for computing diffusivities at midpoints.
         */
        
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > var_cell_data_for_diffusivities;
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > var_side_data_for_diffusivities;
        std::vector<int> var_cell_data_for_diffusivities_component_idx;
        
        diffusive_flux_utilities->getCellDataForInterpolationToSideDataForDiffusiveFluxDiffusivities(
            var_cell_data_for_diffusivities,
            var_cell_data_for_diffusivities_component_idx);
        
        // Interpolate variables from nodes to midpoints for computing diffusivities at midpoints.
        
        interpolateDiffusivitiesFromNodeToMidpoint(
            var_side_data_for_diffusivities,
            var_cell_data_for_diffusivities,
            var_cell_data_for_diffusivities_component_idx,
            patch,
            allocate_scratch_data_containers);
        
        /*
         * Delcare containers for computing fluxes in different directions.
         */
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > > var_data_x;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > > var_data_y;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > > var_data_z;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > > diffusivities_data_x;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > > diffusivities_data_y;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > > diffusivities_data_z;
        
        std::vector<std::vector<int> > var_component_idx_x;
        std::vector<std::vector<int> > var_component_idx_y;
        std::vector<std::vector<int> > var_component_idx_z;
        
        std::vector<std::vector<int> > diffusivities_component_idx_x;
        std::vector<std::vector<int> > diffusivities_component_idx_y;
        std::vector<std::vector<int> > diffusivities_component_idx_z;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > > derivatives_x_midpoint_x;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > > derivatives_y_midpoint_x;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > > derivatives_z_midpoint_x;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > > derivatives_x_midpoint_y;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > > derivatives_y_midpoint_y;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > > derivatives_z_midpoint_y;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > > derivatives_x_midpoint_z;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > > derivatives_y_midpoint_z;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > > derivatives_z_midpoint_z;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > > derivatives_x_node;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > > derivatives_y_node;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > > derivatives_z_node;
        
        std::map<Real*, HAMERS_SHARED_PTR<pdat::SideData<Real> > > derivatives_x_midpoint_x_computed;
        std::map<Real*, HAMERS_SHARED_PTR<pdat::SideData<Real> > > derivatives_y_midpoint_x_computed;
        std::map<Real*, HAMERS_SHARED_PTR<pdat::SideData<Real> > > derivatives_z_midpoint_x_computed;
        
        std::map<Real*, HAMERS_SHARED_PTR<pdat::SideData<Real> > > derivatives_x_midpoint_y_computed;
        std::map<Real*, HAMERS_SHARED_PTR<pdat::SideData<Real> > > derivatives_y_midpoint_y_computed;
        std::map<Real*, HAMERS_SHARED_PTR<pdat::SideData<Real> > > derivatives_z_midpoint_y_computed;
        
        std::map<Real*, HAMERS_SHARED_PTR<pdat::SideData<Real> > > derivatives_x_midpoint_z_computed;
        std::map<Real*, HAMERS_SHARED_PTR<pdat::SideData<Real> > > derivatives_y_midpoint_z_computed;
        std::map<Real*, HAMERS_SHARED_PTR<pdat::SideData<Real> > > derivatives_z_midpoint_z_computed;
        
        std::map<Real*, HAMERS_SHARED_PTR<pdat::CellData<Real> > > derivatives_x_node_computed;
        std::map<Real*, HAMERS_SHARED_PTR<pdat::CellData<Real> > > derivatives_y_node_computed;
        std::map<Real*, HAMERS_SHARED_PTR<pdat::CellData<Real> > > derivatives_z_node_computed;
        
        HAMERS_SHARED_PTR<pdat::SideData<Real> > diffusive_flux_midpoint;
        diffusive_flux_midpoint.reset(new pdat::SideData<Real>(interior_box, d_num_eqn, d_num_diff_ghosts));
        
        std::vector<Real*> F_midpoint_x;
        std::vector<Real*> F_midpoint_y;
        std::vector<Real*> F_midpoint_z;
        F_midpoint_x.reserve(d_num_eqn);
        F_midpoint_y.reserve(d_num_eqn);
        F_midpoint_z.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            F_midpoint_x.push_back(diffusive_flux_midpoint->getPointer(0, ei));
            F_midpoint_y.push_back(diffusive_flux_midpoint->getPointer(1, ei));
            F_midpoint_z.push_back(diffusive_flux_midpoint->getPointer(2, ei));
        }
        
        diffusive_flux_midpoint->fillAll(Real(0));
        
        /*
         * Delcare containers for subgrid scale model.
         */
        
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > var_data_sgs_x;
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > var_data_sgs_y;
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > var_data_sgs_z;
        
        std::vector<int> var_component_idx_sgs_x;
        std::vector<int> var_component_idx_sgs_y;
        std::vector<int> var_component_idx_sgs_z;
        
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > derivatives_sgs_x_midpoint_x;
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > derivatives_sgs_y_midpoint_x;
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > derivatives_sgs_z_midpoint_x;
        
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > derivatives_sgs_x_midpoint_y;
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > derivatives_sgs_y_midpoint_y;
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > derivatives_sgs_z_midpoint_y;
        
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > derivatives_sgs_x_midpoint_z;
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > derivatives_sgs_y_midpoint_z;
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > derivatives_sgs_z_midpoint_z;
        
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > derivatives_sgs_x_node;
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > derivatives_sgs_y_node;
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > derivatives_sgs_z_node;
        
        if (diffusive_flux_utilities->useSubgridScaleModel())
        {
            /*
             * (1) Compute the derivatives at midpoints in the x-direction.
             */
            
            diffusive_flux_utilities->getCellDataOfVariablesForSideDerivativeForSubgridScaleViscosity(
                var_data_sgs_x,
                var_component_idx_sgs_x,
                DIRECTION::X_DIRECTION,
                DIRECTION::X_DIRECTION);
            
            diffusive_flux_utilities->getCellDataOfVariablesForSideDerivativeForSubgridScaleViscosity(
                var_data_sgs_y,
                var_component_idx_sgs_y,
                DIRECTION::X_DIRECTION,
                DIRECTION::Y_DIRECTION);
            
            diffusive_flux_utilities->getCellDataOfVariablesForSideDerivativeForSubgridScaleViscosity(
                var_data_sgs_z,
                var_component_idx_sgs_z,
                DIRECTION::X_DIRECTION,
                DIRECTION::Z_DIRECTION);
            
            /*
             * Compute the derivatives in x-direction.
             */
            
            computeFirstDerivativesInXAtMidpointX(
                derivatives_sgs_x_midpoint_x,
                derivatives_x_midpoint_x_computed,
                var_data_sgs_x,
                var_component_idx_sgs_x,
                patch,
                allocate_scratch_data_containers);
            
            /*
             * Compute the derivatives in y-direction.
             */
            
            computeFirstDerivativesInYAtNode(
                derivatives_sgs_y_node,
                derivatives_y_node_computed,
                var_data_sgs_y,
                var_component_idx_sgs_y,
                patch,
                allocate_scratch_data_containers);
            
            interpolateDerivativesFromNodeToMidpointX(
                derivatives_sgs_y_midpoint_x,
                derivatives_y_midpoint_x_computed,
                derivatives_sgs_y_node,
                var_data_sgs_y,
                var_component_idx_sgs_y,
                patch,
                allocate_scratch_data_containers);
            
            /*
             * Compute the derivatives in z-direction.
             */
            
            computeFirstDerivativesInZAtNode(
                derivatives_sgs_z_node,
                derivatives_z_node_computed,
                var_data_sgs_z,
                var_component_idx_sgs_z,
                patch,
                allocate_scratch_data_containers);
            
            interpolateDerivativesFromNodeToMidpointX(
                derivatives_sgs_z_midpoint_x,
                derivatives_z_midpoint_x_computed,
                derivatives_sgs_z_node,
                var_data_sgs_z,
                var_component_idx_sgs_z,
                patch,
                allocate_scratch_data_containers);
            
            var_data_sgs_x.clear();
            var_data_sgs_y.clear();
            var_data_sgs_z.clear();
            
            var_component_idx_sgs_x.clear();
            var_component_idx_sgs_y.clear();
            var_component_idx_sgs_z.clear();
            
            derivatives_sgs_y_node.clear();
            derivatives_sgs_z_node.clear();
            
            /*
             * (2) Compute the derivatives at midpoints in the y-direction.
             */
            
            diffusive_flux_utilities->getCellDataOfVariablesForSideDerivativeForSubgridScaleViscosity(
                var_data_sgs_x,
                var_component_idx_sgs_x,
                DIRECTION::Y_DIRECTION,
                DIRECTION::X_DIRECTION);
            
            diffusive_flux_utilities->getCellDataOfVariablesForSideDerivativeForSubgridScaleViscosity(
                var_data_sgs_y,
                var_component_idx_sgs_y,
                DIRECTION::Y_DIRECTION,
                DIRECTION::Y_DIRECTION);
            
            diffusive_flux_utilities->getCellDataOfVariablesForSideDerivativeForSubgridScaleViscosity(
                var_data_sgs_z,
                var_component_idx_sgs_z,
                DIRECTION::Y_DIRECTION,
                DIRECTION::Z_DIRECTION);
            
            /*
             * Compute the derivatives in x-direction.
             */
            
            computeFirstDerivativesInXAtNode(
                derivatives_sgs_x_node,
                derivatives_x_node_computed,
                var_data_sgs_x,
                var_component_idx_sgs_x,
                patch,
                allocate_scratch_data_containers);
            
            interpolateDerivativesFromNodeToMidpointY(
                derivatives_sgs_x_midpoint_y,
                derivatives_x_midpoint_y_computed,
                derivatives_sgs_x_node,
                var_data_sgs_x,
                var_component_idx_sgs_x,
                patch,
                allocate_scratch_data_containers);
            
            /*
             * Compute the derivatives in y-direction.
             */
            
            computeFirstDerivativesInYAtMidpointY(
                derivatives_sgs_y_midpoint_y,
                derivatives_y_midpoint_y_computed,
                var_data_sgs_y,
                var_component_idx_sgs_y,
                patch,
                allocate_scratch_data_containers);
            
            /*
             * Compute the derivatives in z-direction.
             */
            
            computeFirstDerivativesInZAtNode(
                derivatives_sgs_z_node,
                derivatives_z_node_computed,
                var_data_sgs_z,
                var_component_idx_sgs_z,
                patch,
                allocate_scratch_data_containers);
            
            interpolateDerivativesFromNodeToMidpointY(
                derivatives_sgs_z_midpoint_y,
                derivatives_z_midpoint_y_computed,
                derivatives_sgs_z_node,
                var_data_sgs_z,
                var_component_idx_sgs_z,
                patch,
                allocate_scratch_data_containers);
            
            var_data_sgs_x.clear();
            var_data_sgs_y.clear();
            var_data_sgs_z.clear();
            
            var_component_idx_sgs_x.clear();
            var_component_idx_sgs_y.clear();
            var_component_idx_sgs_z.clear();
            
            derivatives_sgs_x_node.clear();
            derivatives_sgs_z_node.clear();
            
            /*
             * (3) Compute the derivatives at midpoints in the z-direction.
             */
            
            diffusive_flux_utilities->getCellDataOfVariablesForSideDerivativeForSubgridScaleViscosity(
                var_data_sgs_x,
                var_component_idx_sgs_x,
                DIRECTION::Z_DIRECTION,
                DIRECTION::X_DIRECTION);
            
            diffusive_flux_utilities->getCellDataOfVariablesForSideDerivativeForSubgridScaleViscosity(
                var_data_sgs_y,
                var_component_idx_sgs_y,
                DIRECTION::Z_DIRECTION,
                DIRECTION::Y_DIRECTION);
            
            diffusive_flux_utilities->getCellDataOfVariablesForSideDerivativeForSubgridScaleViscosity(
                var_data_sgs_z,
                var_component_idx_sgs_z,
                DIRECTION::Z_DIRECTION,
                DIRECTION::Z_DIRECTION);
            
            /*
             * Compute the derivatives in x-direction.
             */
            
            computeFirstDerivativesInXAtNode(
                derivatives_sgs_x_node,
                derivatives_x_node_computed,
                var_data_sgs_x,
                var_component_idx_sgs_x,
                patch,
                allocate_scratch_data_containers);
            
            interpolateDerivativesFromNodeToMidpointZ(
                derivatives_sgs_x_midpoint_z,
                derivatives_x_midpoint_z_computed,
                derivatives_sgs_x_node,
                var_data_sgs_x,
                var_component_idx_sgs_x,
                patch,
                allocate_scratch_data_containers);
            
            /*
             * Compute the derivatives in y-direction.
             */
            
            computeFirstDerivativesInYAtNode(
                derivatives_sgs_y_node,
                derivatives_y_node_computed,
                var_data_sgs_y,
                var_component_idx_sgs_y,
                patch,
                allocate_scratch_data_containers);
            
            interpolateDerivativesFromNodeToMidpointZ(
                derivatives_sgs_y_midpoint_z,
                derivatives_y_midpoint_z_computed,
                derivatives_sgs_y_node,
                var_data_sgs_y,
                var_component_idx_sgs_y,
                patch,
                allocate_scratch_data_containers);
            
            /*
             * Compute the derivatives in z-direction.
             */
            
            computeFirstDerivativesInZAtMidpointZ(
                derivatives_sgs_z_midpoint_z,
                derivatives_z_midpoint_z_computed,
                var_data_sgs_z,
                var_component_idx_sgs_z,
                patch,
                allocate_scratch_data_containers);
            
            var_data_sgs_x.clear();
            var_data_sgs_y.clear();
            var_data_sgs_z.clear();
            
            var_component_idx_sgs_x.clear();
            var_component_idx_sgs_y.clear();
            var_component_idx_sgs_z.clear();
            
            derivatives_sgs_x_node.clear();
            derivatives_sgs_y_node.clear();
            
            /*
             * (4) Update the diffusivities with subgrid scale model
             */
            
            std::map<DIRECTION::TYPE, std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > > derivatives_midpoint_x;
            std::map<DIRECTION::TYPE, std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > > derivatives_midpoint_y;
            std::map<DIRECTION::TYPE, std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > > derivatives_midpoint_z;
            
            std::pair<DIRECTION::TYPE, std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > >
                derivatives_x_midpoint_x_pair(
                    DIRECTION::X_DIRECTION,
                    derivatives_sgs_x_midpoint_x);
            
            std::pair<DIRECTION::TYPE, std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > >
                derivatives_y_midpoint_x_pair(
                    DIRECTION::Y_DIRECTION,
                    derivatives_sgs_y_midpoint_x);
            
            std::pair<DIRECTION::TYPE, std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > >
                derivatives_z_midpoint_x_pair(
                    DIRECTION::Z_DIRECTION,
                    derivatives_sgs_z_midpoint_x);
            
            std::pair<DIRECTION::TYPE, std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > >
                derivatives_x_midpoint_y_pair(
                    DIRECTION::X_DIRECTION,
                    derivatives_sgs_x_midpoint_y);
            
            std::pair<DIRECTION::TYPE, std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > >
                derivatives_y_midpoint_y_pair(
                    DIRECTION::Y_DIRECTION,
                    derivatives_sgs_y_midpoint_y);
            
            std::pair<DIRECTION::TYPE, std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > >
                derivatives_z_midpoint_y_pair(
                    DIRECTION::Z_DIRECTION,
                    derivatives_sgs_z_midpoint_y);
            
            std::pair<DIRECTION::TYPE, std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > >
                derivatives_x_midpoint_z_pair(
                    DIRECTION::X_DIRECTION,
                    derivatives_sgs_x_midpoint_z);
            
            std::pair<DIRECTION::TYPE, std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > >
                derivatives_y_midpoint_z_pair(
                    DIRECTION::Y_DIRECTION,
                    derivatives_sgs_y_midpoint_z);
            
            std::pair<DIRECTION::TYPE, std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > >
                derivatives_z_midpoint_z_pair(
                    DIRECTION::Z_DIRECTION,
                    derivatives_sgs_z_midpoint_z);
            
            derivatives_midpoint_x.insert(derivatives_x_midpoint_x_pair);
            derivatives_midpoint_x.insert(derivatives_y_midpoint_x_pair);
            derivatives_midpoint_x.insert(derivatives_z_midpoint_x_pair);
            
            derivatives_midpoint_y.insert(derivatives_x_midpoint_y_pair);
            derivatives_midpoint_y.insert(derivatives_y_midpoint_y_pair);
            derivatives_midpoint_y.insert(derivatives_z_midpoint_y_pair);
            
            derivatives_midpoint_z.insert(derivatives_x_midpoint_z_pair);
            derivatives_midpoint_z.insert(derivatives_y_midpoint_z_pair);
            derivatives_midpoint_z.insert(derivatives_z_midpoint_z_pair);
            
            diffusive_flux_utilities->updateSideDataOfDiffusiveFluxDiffusivitiesWithSubgridScaleModel(
                var_side_data_for_diffusivities,
                derivatives_midpoint_x,
                DIRECTION::X_DIRECTION);
            
            diffusive_flux_utilities->updateSideDataOfDiffusiveFluxDiffusivitiesWithSubgridScaleModel(
                var_side_data_for_diffusivities,
                derivatives_midpoint_y,
                DIRECTION::Y_DIRECTION);
            
            diffusive_flux_utilities->updateSideDataOfDiffusiveFluxDiffusivitiesWithSubgridScaleModel(
                var_side_data_for_diffusivities,
                derivatives_midpoint_z,
                DIRECTION::Z_DIRECTION);
        }
        
        
        diffusive_flux_utilities->allocateMemoryForSideDataOfDiffusiveFluxDiffusivities();
        
        diffusive_flux_utilities->computeSideDataOfDiffusiveFluxDiffusivities(
            var_side_data_for_diffusivities);
        
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
        
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxVariablesForDerivative(
            var_data_z,
            var_component_idx_z,
            DIRECTION::X_DIRECTION,
            DIRECTION::Z_DIRECTION);
        
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
        
        diffusive_flux_utilities->getSideDataOfDiffusiveFluxDiffusivities(
            diffusivities_data_z,
            diffusivities_component_idx_z,
            DIRECTION::X_DIRECTION,
            DIRECTION::Z_DIRECTION);
        
        TBOX_ASSERT(static_cast<int>(var_data_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_data_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_data_z.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_z.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_z.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_z.size()) == d_num_eqn);
        
        /*
         * Compute the derivatives in x-direction for diffusive flux in x-direction.
         */
        
        computeFirstDerivativesInXAtMidpointX(
            derivatives_x_midpoint_x,
            derivatives_x_midpoint_x_computed,
            var_data_x,
            var_component_idx_x,
            patch,
            allocate_scratch_data_containers);
        
        /*
         * Compute the derivatives in y-direction for diffusive flux in x-direction.
         */
        
        computeFirstDerivativesInYAtNode(
            derivatives_y_node,
            derivatives_y_node_computed,
            var_data_y,
            var_component_idx_y,
            patch,
            allocate_scratch_data_containers);
        
        interpolateDerivativesFromNodeToMidpointX(
            derivatives_y_midpoint_x,
            derivatives_y_midpoint_x_computed,
            derivatives_y_node,
            var_data_y,
            var_component_idx_y,
            patch,
            allocate_scratch_data_containers);
        
        /*
         * Compute the derivatives in z-direction for diffusive flux in x-direction.
         */
        
        computeFirstDerivativesInZAtNode(
            derivatives_z_node,
            derivatives_z_node_computed,
            var_data_z,
            var_component_idx_z,
            patch,
            allocate_scratch_data_containers);
        
        interpolateDerivativesFromNodeToMidpointX(
            derivatives_z_midpoint_x,
            derivatives_z_midpoint_x_computed,
            derivatives_z_node,
            var_data_z,
            var_component_idx_z,
            patch,
            allocate_scratch_data_containers);
        
        /*
         * Compute diffusive flux in x-direction at midpoints.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ASSERT(static_cast<int>(derivatives_x_midpoint_x[ei].size()) ==
                        static_cast<int>(diffusivities_data_x[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_x[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_x[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_x[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_x[ei][vi];
                
                // Get the pointer to diffusivity.
                Real* mu = diffusivities_data_x[ei][vi]->getPointer(0, mu_idx);
                
                // Get the pointer to derivative.
                Real* dudx = derivatives_x_midpoint_x[ei][vi]->getPointer(0, 0);
                
                /*
                 * Get the sub-ghost cell width and ghost box dimensions of the diffusivity.
                 */
                
                hier::IntVector num_subghosts_diffusivity =
                    diffusivities_data_x[ei][vi]->getGhostCellWidth();
                
                hier::IntVector subghostcell_dims_diffusivity =
                    diffusivities_data_x[ei][vi]->getGhostBox().numberCells();
                
                const int num_subghosts_0_diffusivity = num_subghosts_diffusivity[0];
                const int num_subghosts_1_diffusivity = num_subghosts_diffusivity[1];
                const int num_subghosts_2_diffusivity = num_subghosts_diffusivity[2];
                const int subghostcell_dim_0_diffusivity = subghostcell_dims_diffusivity[0];
                const int subghostcell_dim_1_diffusivity = subghostcell_dims_diffusivity[1];
                
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = -num_subghosts_0_diffusivity + num_flux_reconstruct_ghosts_0;
                            i < (interior_dim_0 + 1) + num_subghosts_0_diffusivity - num_flux_reconstruct_ghosts_0;
                            i++)
                        {
                            // Compute the linear indices.
                            const int idx = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*(diff_ghostcell_dim_0 + 1) +
                                (k + num_diff_ghosts_2)*(diff_ghostcell_dim_0 + 1)*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_diffusivity = (i + num_subghosts_0_diffusivity) +
                                (j + num_subghosts_1_diffusivity)*(subghostcell_dim_0_diffusivity + 1) +
                                (k + num_subghosts_2_diffusivity)*(subghostcell_dim_0_diffusivity + 1)*
                                    subghostcell_dim_1_diffusivity;
                            
                            F_midpoint_x[ei][idx] += mu[idx_diffusivity]*dudx[idx];
                        }
                    }
                }
            }
            
            TBOX_ASSERT(static_cast<int>(derivatives_y_midpoint_x[ei].size()) ==
                        static_cast<int>(diffusivities_data_y[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_y[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_y[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_y[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_y[ei][vi];
                
                // Get the pointer to diffusivity.
                Real* mu = diffusivities_data_y[ei][vi]->getPointer(0, mu_idx);
                
                // Get the pointer to derivative.
                Real* dudy = derivatives_y_midpoint_x[ei][vi]->getPointer(0, 0);
                
                /*
                 * Get the sub-ghost cell width and ghost box dimensions of the diffusivity.
                 */
                
                hier::IntVector num_subghosts_diffusivity =
                    diffusivities_data_y[ei][vi]->getGhostCellWidth();
                
                hier::IntVector subghostcell_dims_diffusivity =
                    diffusivities_data_y[ei][vi]->getGhostBox().numberCells();
                
                const int num_subghosts_0_diffusivity = num_subghosts_diffusivity[0];
                const int num_subghosts_1_diffusivity = num_subghosts_diffusivity[1];
                const int num_subghosts_2_diffusivity = num_subghosts_diffusivity[2];
                const int subghostcell_dim_0_diffusivity = subghostcell_dims_diffusivity[0];
                const int subghostcell_dim_1_diffusivity = subghostcell_dims_diffusivity[1];
                
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = -num_subghosts_0_diffusivity + num_flux_reconstruct_ghosts_0;
                            i < (interior_dim_0 + 1) + num_subghosts_0_diffusivity - num_flux_reconstruct_ghosts_0;
                            i++)
                        {
                            // Compute the linear indices.
                            const int idx = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*(diff_ghostcell_dim_0 + 1) +
                                (k + num_diff_ghosts_2)*(diff_ghostcell_dim_0 + 1)*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_diffusivity = (i + num_subghosts_0_diffusivity) +
                                (j + num_subghosts_1_diffusivity)*(subghostcell_dim_0_diffusivity + 1) +
                                (k + num_subghosts_2_diffusivity)*(subghostcell_dim_0_diffusivity + 1)*
                                    subghostcell_dim_1_diffusivity;
                            
                            F_midpoint_x[ei][idx] += mu[idx_diffusivity]*dudy[idx];
                        }
                    }
                }
            }
            
            TBOX_ASSERT(static_cast<int>(derivatives_z_midpoint_x[ei].size()) ==
                        static_cast<int>(diffusivities_data_z[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_z[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_z[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_z[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_z[ei][vi];
                
                // Get the pointer to diffusivity.
                Real* mu = diffusivities_data_z[ei][vi]->getPointer(0, mu_idx);
                
                // Get the pointer to derivative.
                Real* dudz = derivatives_z_midpoint_x[ei][vi]->getPointer(0, 0);
                
                /*
                 * Get the sub-ghost cell width and ghost box dimensions of the diffusivity.
                 */
                
                hier::IntVector num_subghosts_diffusivity =
                    diffusivities_data_z[ei][vi]->getGhostCellWidth();
                
                hier::IntVector subghostcell_dims_diffusivity =
                    diffusivities_data_z[ei][vi]->getGhostBox().numberCells();
                
                const int num_subghosts_0_diffusivity = num_subghosts_diffusivity[0];
                const int num_subghosts_1_diffusivity = num_subghosts_diffusivity[1];
                const int num_subghosts_2_diffusivity = num_subghosts_diffusivity[2];
                const int subghostcell_dim_0_diffusivity = subghostcell_dims_diffusivity[0];
                const int subghostcell_dim_1_diffusivity = subghostcell_dims_diffusivity[1];
                
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = -num_subghosts_0_diffusivity + num_flux_reconstruct_ghosts_0;
                            i < (interior_dim_0 + 1) + num_subghosts_0_diffusivity - num_flux_reconstruct_ghosts_0;
                            i++)
                        {
                            // Compute the linear indices.
                            const int idx = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*(diff_ghostcell_dim_0 + 1) +
                                (k + num_diff_ghosts_2)*(diff_ghostcell_dim_0 + 1)*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_diffusivity = (i + num_subghosts_0_diffusivity) +
                                (j + num_subghosts_1_diffusivity)*(subghostcell_dim_0_diffusivity + 1) +
                                (k + num_subghosts_2_diffusivity)*(subghostcell_dim_0_diffusivity + 1)*
                                    subghostcell_dim_1_diffusivity;
                            
                            F_midpoint_x[ei][idx] += mu[idx_diffusivity]*dudz[idx];
                        }
                    }
                }
            }
        }
        
        /*
         * Reconstruct the flux in x-direction.
         */
        
        reconstructFluxX(
            diffusive_flux,
            diffusive_flux_midpoint,
            patch,
            dt);
        
        var_data_x.clear();
        var_data_y.clear();
        var_data_z.clear();
        
        diffusivities_data_x.clear();
        diffusivities_data_y.clear();
        diffusivities_data_z.clear();
        
        var_component_idx_x.clear();
        var_component_idx_y.clear();
        var_component_idx_z.clear();
        
        diffusivities_component_idx_x.clear();
        diffusivities_component_idx_y.clear();
        diffusivities_component_idx_z.clear();
        
        derivatives_x_midpoint_x.clear();
        derivatives_y_midpoint_x.clear();
        derivatives_z_midpoint_x.clear();
        
        derivatives_y_node.clear();
        derivatives_z_node.clear();
        
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
        
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxVariablesForDerivative(
            var_data_z,
            var_component_idx_z,
            DIRECTION::Y_DIRECTION,
            DIRECTION::Z_DIRECTION);
        
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
        
        diffusive_flux_utilities->getSideDataOfDiffusiveFluxDiffusivities(
            diffusivities_data_z,
            diffusivities_component_idx_z,
            DIRECTION::Y_DIRECTION,
            DIRECTION::Z_DIRECTION);
        
        TBOX_ASSERT(static_cast<int>(var_data_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_data_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_data_z.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_z.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_z.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_z.size()) == d_num_eqn);
        
        /*
         * Compute the derivatives in x-direction for diffusive flux in y-direction.
         */
        
        computeFirstDerivativesInXAtNode(
            derivatives_x_node,
            derivatives_x_node_computed,
            var_data_x,
            var_component_idx_x,
            patch,
            allocate_scratch_data_containers);
        
        interpolateDerivativesFromNodeToMidpointY(
            derivatives_x_midpoint_y,
            derivatives_x_midpoint_y_computed,
            derivatives_x_node,
            var_data_x,
            var_component_idx_x,
            patch,
            allocate_scratch_data_containers);
        
        /*
         * Compute the derivatives in y-direction for diffusive flux in y-direction.
         */
        
        computeFirstDerivativesInYAtMidpointY(
            derivatives_y_midpoint_y,
            derivatives_y_midpoint_y_computed,
            var_data_y,
            var_component_idx_y,
            patch,
            allocate_scratch_data_containers);
        
        /*
         * Compute the derivatives in z-direction for diffusive flux in y-direction.
         */
        
        computeFirstDerivativesInZAtNode(
            derivatives_z_node,
            derivatives_z_node_computed,
            var_data_z,
            var_component_idx_z,
            patch,
            allocate_scratch_data_containers);
        
        interpolateDerivativesFromNodeToMidpointY(
            derivatives_z_midpoint_y,
            derivatives_z_midpoint_y_computed,
            derivatives_z_node,
            var_data_z,
            var_component_idx_z,
            patch,
            allocate_scratch_data_containers);
        
        /*
         * Compute diffusive flux in y-direction at midpoints.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ASSERT(static_cast<int>(derivatives_x_midpoint_y[ei].size()) ==
                        static_cast<int>(diffusivities_data_x[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_x[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_x[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_x[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_x[ei][vi];
                
                // Get the pointer to diffusivity.
                Real* mu = diffusivities_data_x[ei][vi]->getPointer(1, mu_idx);
                
                // Get the pointer to derivative.
                Real* dudx = derivatives_x_midpoint_y[ei][vi]->getPointer(1, 0);
                
                /*
                 * Get the sub-ghost cell width and ghost box dimensions of the diffusivity.
                 */
                
                hier::IntVector num_subghosts_diffusivity =
                    diffusivities_data_x[ei][vi]->getGhostCellWidth();
                
                hier::IntVector subghostcell_dims_diffusivity =
                    diffusivities_data_x[ei][vi]->getGhostBox().numberCells();
                
                const int num_subghosts_0_diffusivity = num_subghosts_diffusivity[0];
                const int num_subghosts_1_diffusivity = num_subghosts_diffusivity[1];
                const int num_subghosts_2_diffusivity = num_subghosts_diffusivity[2];
                const int subghostcell_dim_0_diffusivity = subghostcell_dims_diffusivity[0];
                const int subghostcell_dim_1_diffusivity = subghostcell_dims_diffusivity[1];
                
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = -num_subghosts_1_diffusivity + num_flux_reconstruct_ghosts_1;
                        j < (interior_dim_1 + 1) + num_subghosts_1_diffusivity - num_flux_reconstruct_ghosts_1;
                        j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    (diff_ghostcell_dim_1 + 1);
                            
                            const int idx_diffusivity = (i + num_subghosts_0_diffusivity) +
                                (j + num_subghosts_1_diffusivity)*subghostcell_dim_0_diffusivity +
                                (k + num_subghosts_2_diffusivity)*subghostcell_dim_0_diffusivity*
                                    (subghostcell_dim_1_diffusivity + 1);
                            
                            F_midpoint_y[ei][idx] += mu[idx_diffusivity]*dudx[idx];
                        }
                    }
                }
            }
            
            TBOX_ASSERT(static_cast<int>(derivatives_y_midpoint_y[ei].size()) ==
                        static_cast<int>(diffusivities_data_y[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_y[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_y[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_y[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_y[ei][vi];
                
                // Get the pointer to diffusivity.
                Real* mu = diffusivities_data_y[ei][vi]->getPointer(1, mu_idx);
                
                // Get the pointer to derivative.
                Real* dudy = derivatives_y_midpoint_y[ei][vi]->getPointer(1, 0);
                
                /*
                 * Get the sub-ghost cell width and ghost box dimensions of the diffusivity.
                 */
                
                hier::IntVector num_subghosts_diffusivity =
                    diffusivities_data_y[ei][vi]->getGhostCellWidth();
                
                hier::IntVector subghostcell_dims_diffusivity =
                    diffusivities_data_y[ei][vi]->getGhostBox().numberCells();
                
                const int num_subghosts_0_diffusivity = num_subghosts_diffusivity[0];
                const int num_subghosts_1_diffusivity = num_subghosts_diffusivity[1];
                const int num_subghosts_2_diffusivity = num_subghosts_diffusivity[2];
                const int subghostcell_dim_0_diffusivity = subghostcell_dims_diffusivity[0];
                const int subghostcell_dim_1_diffusivity = subghostcell_dims_diffusivity[1];
                
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = -num_subghosts_1_diffusivity + num_flux_reconstruct_ghosts_1;
                        j < (interior_dim_1 + 1) + num_subghosts_1_diffusivity - num_flux_reconstruct_ghosts_1;
                        j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    (diff_ghostcell_dim_1 + 1);
                            
                            const int idx_diffusivity = (i + num_subghosts_0_diffusivity) +
                                (j + num_subghosts_1_diffusivity)*subghostcell_dim_0_diffusivity +
                                (k + num_subghosts_2_diffusivity)*subghostcell_dim_0_diffusivity*
                                    (subghostcell_dim_1_diffusivity + 1);
                            
                            F_midpoint_y[ei][idx] += mu[idx_diffusivity]*dudy[idx];
                        }
                    }
                }
            }
            
            TBOX_ASSERT(static_cast<int>(derivatives_z_midpoint_y[ei].size()) ==
                        static_cast<int>(diffusivities_data_z[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_z[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_z[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_z[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_z[ei][vi];
                
                // Get the pointer to diffusivity.
                Real* mu = diffusivities_data_z[ei][vi]->getPointer(1, mu_idx);
                
                // Get the pointer to derivative.
                Real* dudz = derivatives_z_midpoint_y[ei][vi]->getPointer(1, 0);
                
                /*
                 * Get the sub-ghost cell width and ghost box dimensions of the diffusivity.
                 */
                
                hier::IntVector num_subghosts_diffusivity =
                    diffusivities_data_z[ei][vi]->getGhostCellWidth();
                
                hier::IntVector subghostcell_dims_diffusivity =
                    diffusivities_data_z[ei][vi]->getGhostBox().numberCells();
                
                const int num_subghosts_0_diffusivity = num_subghosts_diffusivity[0];
                const int num_subghosts_1_diffusivity = num_subghosts_diffusivity[1];
                const int num_subghosts_2_diffusivity = num_subghosts_diffusivity[2];
                const int subghostcell_dim_0_diffusivity = subghostcell_dims_diffusivity[0];
                const int subghostcell_dim_1_diffusivity = subghostcell_dims_diffusivity[1];
                
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = -num_subghosts_1_diffusivity + num_flux_reconstruct_ghosts_1;
                        j < (interior_dim_1 + 1) + num_subghosts_1_diffusivity - num_flux_reconstruct_ghosts_1;
                        j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    (diff_ghostcell_dim_1 + 1);
                            
                            const int idx_diffusivity = (i + num_subghosts_0_diffusivity) +
                                (j + num_subghosts_1_diffusivity)*subghostcell_dim_0_diffusivity +
                                (k + num_subghosts_2_diffusivity)*subghostcell_dim_0_diffusivity*
                                    (subghostcell_dim_1_diffusivity + 1);
                            
                            F_midpoint_y[ei][idx] += mu[idx_diffusivity]*dudz[idx];
                        }
                    }
                }
            }
        }
        
        /*
         * Reconstruct the flux in y-direction.
         */
        
        reconstructFluxY(
            diffusive_flux,
            diffusive_flux_midpoint,
            patch,
            dt);
        
        var_data_x.clear();
        var_data_y.clear();
        var_data_z.clear();
        
        diffusivities_data_x.clear();
        diffusivities_data_y.clear();
        diffusivities_data_z.clear();
        
        var_component_idx_x.clear();
        var_component_idx_y.clear();
        var_component_idx_z.clear();
        
        diffusivities_component_idx_x.clear();
        diffusivities_component_idx_y.clear();
        diffusivities_component_idx_z.clear();
        
        derivatives_x_midpoint_y.clear();
        derivatives_y_midpoint_y.clear();
        derivatives_z_midpoint_y.clear();
        
        derivatives_x_node.clear();
        derivatives_z_node.clear();
        
        /*
         * (3) Compute the flux in the z-direction.
         */
        
        // Get the variables for the derivatives in the diffusive flux.
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxVariablesForDerivative(
            var_data_x,
            var_component_idx_x,
            DIRECTION::Z_DIRECTION,
            DIRECTION::X_DIRECTION);
        
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxVariablesForDerivative(
            var_data_y,
            var_component_idx_y,
            DIRECTION::Z_DIRECTION,
            DIRECTION::Y_DIRECTION);
        
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxVariablesForDerivative(
            var_data_z,
            var_component_idx_z,
            DIRECTION::Z_DIRECTION,
            DIRECTION::Z_DIRECTION);
        
        // Get the diffusivities in the diffusive flux.
        diffusive_flux_utilities->getSideDataOfDiffusiveFluxDiffusivities(
            diffusivities_data_x,
            diffusivities_component_idx_x,
            DIRECTION::Z_DIRECTION,
            DIRECTION::X_DIRECTION);
        
        diffusive_flux_utilities->getSideDataOfDiffusiveFluxDiffusivities(
            diffusivities_data_y,
            diffusivities_component_idx_y,
            DIRECTION::Z_DIRECTION,
            DIRECTION::Y_DIRECTION);
        
        diffusive_flux_utilities->getSideDataOfDiffusiveFluxDiffusivities(
            diffusivities_data_z,
            diffusivities_component_idx_z,
            DIRECTION::Z_DIRECTION,
            DIRECTION::Z_DIRECTION);
        
        TBOX_ASSERT(static_cast<int>(var_data_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_data_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_data_z.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_z.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_z.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_z.size()) == d_num_eqn);
        
        /*
         * Compute the derivatives in x-direction for diffusive flux in z-direction.
         */
        
        computeFirstDerivativesInXAtNode(
            derivatives_x_node,
            derivatives_x_node_computed,
            var_data_x,
            var_component_idx_x,
            patch,
            allocate_scratch_data_containers);
        
        interpolateDerivativesFromNodeToMidpointZ(
            derivatives_x_midpoint_z,
            derivatives_x_midpoint_z_computed,
            derivatives_x_node,
            var_data_x,
            var_component_idx_x,
            patch,
            allocate_scratch_data_containers);
        
        /*
         * Compute the derivatives in y-direction for diffusive flux in z-direction.
         */
        
        computeFirstDerivativesInYAtNode(
            derivatives_y_node,
            derivatives_y_node_computed,
            var_data_y,
            var_component_idx_y,
            patch,
            allocate_scratch_data_containers);
        
        interpolateDerivativesFromNodeToMidpointZ(
            derivatives_y_midpoint_z,
            derivatives_y_midpoint_z_computed,
            derivatives_y_node,
            var_data_y,
            var_component_idx_y,
            patch,
            allocate_scratch_data_containers);
        
        /*
         * Compute the derivatives in z-direction for diffusive flux in z-direction.
         */
        
        computeFirstDerivativesInZAtMidpointZ(
            derivatives_z_midpoint_z,
            derivatives_z_midpoint_z_computed,
            var_data_z,
            var_component_idx_z,
            patch,
            allocate_scratch_data_containers);
        
        /*
         * Compute diffusive flux in z-direction at midpoints.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ASSERT(static_cast<int>(derivatives_x_midpoint_z[ei].size()) ==
                        static_cast<int>(diffusivities_data_x[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_x[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_x[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_x[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_x[ei][vi];
                
                // Get the pointer to diffusivity.
                Real* mu = diffusivities_data_x[ei][vi]->getPointer(2, mu_idx);
                
                // Get the pointer to derivative.
                Real* dudx = derivatives_x_midpoint_z[ei][vi]->getPointer(2, 0);
                
                /*
                 * Get the sub-ghost cell width and ghost box dimensions of the diffusivity.
                 */
                
                hier::IntVector num_subghosts_diffusivity =
                    diffusivities_data_x[ei][vi]->getGhostCellWidth();
                
                hier::IntVector subghostcell_dims_diffusivity =
                    diffusivities_data_x[ei][vi]->getGhostBox().numberCells();
                
                const int num_subghosts_0_diffusivity = num_subghosts_diffusivity[0];
                const int num_subghosts_1_diffusivity = num_subghosts_diffusivity[1];
                const int num_subghosts_2_diffusivity = num_subghosts_diffusivity[2];
                const int subghostcell_dim_0_diffusivity = subghostcell_dims_diffusivity[0];
                const int subghostcell_dim_1_diffusivity = subghostcell_dims_diffusivity[1];
                
                for (int k = -num_subghosts_2_diffusivity + num_flux_reconstruct_ghosts_2;
                    k < (interior_dim_2 + 1) + num_subghosts_2_diffusivity - num_flux_reconstruct_ghosts_2;
                    k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_diffusivity = (i + num_subghosts_0_diffusivity) +
                                (j + num_subghosts_1_diffusivity)*subghostcell_dim_0_diffusivity +
                                (k + num_subghosts_2_diffusivity)*subghostcell_dim_0_diffusivity*
                                    subghostcell_dim_1_diffusivity;
                            
                            F_midpoint_z[ei][idx] += mu[idx_diffusivity]*dudx[idx];
                        }
                    }
                }
            }
            
            TBOX_ASSERT(static_cast<int>(derivatives_y_midpoint_z[ei].size()) ==
                        static_cast<int>(diffusivities_data_y[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_y[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_y[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_y[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_y[ei][vi];
                
                // Get the pointer to diffusivity.
                Real* mu = diffusivities_data_y[ei][vi]->getPointer(2, mu_idx);
                
                // Get the pointer to derivative.
                Real* dudy = derivatives_y_midpoint_z[ei][vi]->getPointer(2, 0);
                
                /*
                 * Get the sub-ghost cell width and ghost box dimensions of the diffusivity.
                 */
                
                hier::IntVector num_subghosts_diffusivity =
                    diffusivities_data_y[ei][vi]->getGhostCellWidth();
                
                hier::IntVector subghostcell_dims_diffusivity =
                    diffusivities_data_y[ei][vi]->getGhostBox().numberCells();
                
                const int num_subghosts_0_diffusivity = num_subghosts_diffusivity[0];
                const int num_subghosts_1_diffusivity = num_subghosts_diffusivity[1];
                const int num_subghosts_2_diffusivity = num_subghosts_diffusivity[2];
                const int subghostcell_dim_0_diffusivity = subghostcell_dims_diffusivity[0];
                const int subghostcell_dim_1_diffusivity = subghostcell_dims_diffusivity[1];
                
                for (int k = -num_subghosts_2_diffusivity + num_flux_reconstruct_ghosts_2;
                    k < (interior_dim_2 + 1) + num_subghosts_2_diffusivity - num_flux_reconstruct_ghosts_2;
                    k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_diffusivity = (i + num_subghosts_0_diffusivity) +
                                (j + num_subghosts_1_diffusivity)*subghostcell_dim_0_diffusivity +
                                (k + num_subghosts_2_diffusivity)*subghostcell_dim_0_diffusivity*
                                    subghostcell_dim_1_diffusivity;
                            
                            F_midpoint_z[ei][idx] += mu[idx_diffusivity]*dudy[idx];
                        }
                    }
                }
            }
            
            TBOX_ASSERT(static_cast<int>(derivatives_z_midpoint_z[ei].size()) ==
                        static_cast<int>(diffusivities_data_z[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_z[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_z[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_z[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_z[ei][vi];
                
                // Get the pointer to diffusivity.
                Real* mu = diffusivities_data_z[ei][vi]->getPointer(2, mu_idx);
                
                // Get the pointer to derivative.
                Real* dudz = derivatives_z_midpoint_z[ei][vi]->getPointer(2, 0);
                
                /*
                 * Get the sub-ghost cell width and ghost box dimensions of the diffusivity.
                 */
                
                hier::IntVector num_subghosts_diffusivity =
                    diffusivities_data_z[ei][vi]->getGhostCellWidth();
                
                hier::IntVector subghostcell_dims_diffusivity =
                    diffusivities_data_z[ei][vi]->getGhostBox().numberCells();
                
                const int num_subghosts_0_diffusivity = num_subghosts_diffusivity[0];
                const int num_subghosts_1_diffusivity = num_subghosts_diffusivity[1];
                const int num_subghosts_2_diffusivity = num_subghosts_diffusivity[2];
                const int subghostcell_dim_0_diffusivity = subghostcell_dims_diffusivity[0];
                const int subghostcell_dim_1_diffusivity = subghostcell_dims_diffusivity[1];
                
                for (int k = -num_subghosts_2_diffusivity + num_flux_reconstruct_ghosts_2;
                    k < (interior_dim_2 + 1) + num_subghosts_2_diffusivity - num_flux_reconstruct_ghosts_2;
                    k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_diffusivity = (i + num_subghosts_0_diffusivity) +
                                (j + num_subghosts_1_diffusivity)*subghostcell_dim_0_diffusivity +
                                (k + num_subghosts_2_diffusivity)*subghostcell_dim_0_diffusivity*
                                    subghostcell_dim_1_diffusivity;
                            
                            F_midpoint_z[ei][idx] += mu[idx_diffusivity]*dudz[idx];
                        }
                    }
                }
            }
        }
        
        /*
         * Reconstruct the flux in z-direction.
         */
        
        reconstructFluxZ(
            diffusive_flux,
            diffusive_flux_midpoint,
            patch,
            dt);
        
        var_data_x.clear();
        var_data_y.clear();
        var_data_z.clear();
        
        diffusivities_data_x.clear();
        diffusivities_data_y.clear();
        diffusivities_data_z.clear();
        
        var_component_idx_x.clear();
        var_component_idx_y.clear();
        var_component_idx_z.clear();
        
        diffusivities_component_idx_x.clear();
        diffusivities_component_idx_y.clear();
        diffusivities_component_idx_z.clear();
        
        derivatives_x_midpoint_z.clear();
        derivatives_y_midpoint_z.clear();
        derivatives_z_midpoint_z.clear();
        
        derivatives_x_node.clear();
        derivatives_y_node.clear();
        
        /*
         * Unregister the patch and data of all registered derived cell variables in the flow model.
         */
        
        d_flow_model->unregisterPatch();
        
    } // if (d_dim == tbox::Dimension(3))
    
    // Clear the scratch data.
    d_scratch_derivatives_node.clear();
    d_scratch_derivatives_midpoint_x.clear();
    d_scratch_derivatives_midpoint_y.clear();
    d_scratch_derivatives_midpoint_z.clear();
    d_scratch_diffusivities_midpoint.clear();
    
    d_num_scratch_derivatives_node_used = 0;
    d_num_scratch_derivatives_midpoint_x_used = 0;
    d_num_scratch_derivatives_midpoint_y_used = 0;
    d_num_scratch_derivatives_midpoint_z_used = 0;
    d_num_scratch_diffusivities_midpoint_used = 0;
}


/*
 * Compute the derivatives in x-direction at midpoints.
 */
void
DiffusiveFluxReconstructorMidpoint::computeFirstDerivativesInXAtMidpointX(
    std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& derivatives_x_midpoint_x,
    std::map<Real*, HAMERS_SHARED_PTR<pdat::SideData<Real> > >& derivatives_x_midpoint_x_computed,
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& data_x,
    const std::vector<int>& data_component_idx_x,
    const hier::Patch& patch,
    const bool allocate_scratch_derivatives_midpoint)
{
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // diffusive ghost cells.
    hier::Box diff_ghost_box = interior_box;
    diff_ghost_box.grow(d_num_diff_ghosts);
    const hier::IntVector diff_ghostcell_dims = diff_ghost_box.numberCells();
    
    // Get the grid spacing.
    const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
        HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
    const double* const dx = patch_geom->getDx();
    
    const Real dx_0_inv = Real(1)/Real(dx[0]);
    
    hier::IntVector direction = hier::IntVector::getZero(d_dim);
    direction[0] = 1;
    
    derivatives_x_midpoint_x.reserve(static_cast<int>(data_x.size()));
    
    for (int vi = 0; vi < static_cast<int>(data_x.size()); vi++)
    {
        // Get the index of variable for derivative.
        const int u_idx = data_component_idx_x[vi];
        
        if (derivatives_x_midpoint_x_computed.find(data_x[vi]->getPointer(u_idx))
            == derivatives_x_midpoint_x_computed.end())
        {
            // Get the pointer to variable for derivative.
            Real* u = data_x[vi]->getPointer(u_idx);
            
            d_num_scratch_derivatives_midpoint_x_used++;
            // Declare container to store the derivative.
            if (allocate_scratch_derivatives_midpoint)
            {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
                TBOX_ASSERT(static_cast<int>(d_scratch_derivatives_midpoint_x.size()) == d_num_scratch_derivatives_midpoint_x_used - 1);
#endif
                d_scratch_derivatives_midpoint_x.push_back(HAMERS_SHARED_PTR<pdat::SideData<Real> >(
                    new pdat::SideData<Real>(
                        interior_box, 1, d_num_diff_ghosts, direction)));
            }
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
            else
            {
                TBOX_ASSERT(static_cast<int>(d_scratch_derivatives_midpoint_x.size()) >= d_num_scratch_derivatives_midpoint_x_used);
            }
#endif
            
            // Get the pointer to the derivative at midpoint.
            Real* dudx = d_scratch_derivatives_midpoint_x[d_num_scratch_derivatives_midpoint_x_used - 1]->getPointer(0, 0);
            
            /*
             * Get the ghost cell widths and ghost box dimensions of the variables.
             */
            
            const hier::IntVector& num_ghosts_derivative_midpoint = d_num_diff_ghosts;
            
            const hier::IntVector& ghostcell_dims_derivative_midpoint = diff_ghostcell_dims;
            
            const hier::IntVector& num_ghosts_data_node =
                data_x[vi]->getGhostCellWidth();
            
            const hier::IntVector& ghostcell_dims_data_node =
                data_x[vi]->getGhostBox().numberCells();
            
            hier::IntVector domain_lo(d_dim);
            hier::IntVector domain_dims(d_dim);
            
            domain_lo = hier::IntVector::getZero(d_dim);
            domain_dims = interior_dims;
            
            domain_lo[0] += (-d_num_diff_ghosts[0] + d_num_der_midpoint_ghosts[0]);
            domain_dims[0] += 1;
            domain_dims[0] += 2*(d_num_diff_ghosts[0] - d_num_der_midpoint_ghosts[0]);
            
            computeFirstDerivativesInXAtMidpointX(
                dudx,
                u,
                num_ghosts_derivative_midpoint,
                num_ghosts_data_node,
                ghostcell_dims_derivative_midpoint,
                ghostcell_dims_data_node,
                domain_lo,
                domain_dims,
                dx_0_inv);
            
            std::pair<Real*, HAMERS_SHARED_PTR<pdat::SideData<Real> > > derivative_pair(
                u,
                d_scratch_derivatives_midpoint_x[d_num_scratch_derivatives_midpoint_x_used - 1]);
            
            derivatives_x_midpoint_x_computed.insert(derivative_pair);
        }
        
        derivatives_x_midpoint_x.push_back(
            derivatives_x_midpoint_x_computed.find(data_x[vi]->getPointer(u_idx))->
            second);
    }
}


/*
 * Compute the derivatives in x-direction at midpoints.
 */
void
DiffusiveFluxReconstructorMidpoint::computeFirstDerivativesInXAtMidpointX(
    std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > >& derivatives_x_midpoint_x,
    std::map<Real*, HAMERS_SHARED_PTR<pdat::SideData<Real> > >& derivatives_x_midpoint_x_computed,
    const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > >& data_x,
    const std::vector<std::vector<int> >& data_component_idx_x,
    const hier::Patch& patch,
    const bool allocate_scratch_derivatives_midpoint)
{
    const int num_eqn = static_cast<int>(data_x.size());
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(d_num_eqn == num_eqn);
    TBOX_ASSERT(static_cast<int>(data_component_idx_x.size()) == num_eqn);
    
    for (int ei = 0; ei < num_eqn; ei++)
    {
        TBOX_ASSERT(static_cast<int>(data_x[ei].size()) ==
                    static_cast<int>(data_component_idx_x[ei].size()));
    }
#endif
    
    derivatives_x_midpoint_x.resize(num_eqn);
    
    for (int ei = 0; ei < num_eqn; ei++)
    {
        computeFirstDerivativesInXAtMidpointX(
            derivatives_x_midpoint_x[ei],
            derivatives_x_midpoint_x_computed,
            data_x[ei],
            data_component_idx_x[ei],
            patch,
            allocate_scratch_derivatives_midpoint);
    }
}


/*
 * Compute the derivatives in y-direction at midpoints.
 */
void
DiffusiveFluxReconstructorMidpoint::computeFirstDerivativesInYAtMidpointY(
    std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& derivatives_y_midpoint_y,
    std::map<Real*, HAMERS_SHARED_PTR<pdat::SideData<Real> > >& derivatives_y_midpoint_y_computed,
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& data_y,
    const std::vector<int>& data_component_idx_y,
    const hier::Patch& patch,
    const bool allocate_scratch_derivatives_midpoint)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": DiffusiveFluxReconstructorMidpoint::"
            << "computeFirstDerivativesInYAtMidpointY()\n"
            << "There isn't y-direction for one-dimensional problem."
            << std::endl);
    }
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // diffusive ghost cells.
    hier::Box diff_ghost_box = interior_box;
    diff_ghost_box.grow(d_num_diff_ghosts);
    const hier::IntVector diff_ghostcell_dims = diff_ghost_box.numberCells();
    
    // Get the grid spacing.
    const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
        HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
    const double* const dx = patch_geom->getDx();
    
    const Real dx_1_inv = Real(1)/Real(dx[1]);
    
    hier::IntVector direction = hier::IntVector::getZero(d_dim);
    direction[1] = 1;
    
    derivatives_y_midpoint_y.reserve(static_cast<int>(data_y.size()));
    
    for (int vi = 0; vi < static_cast<int>(data_y.size()); vi++)
    {
        // Get the index of variable for derivative.
        const int u_idx = data_component_idx_y[vi];
        
        if (derivatives_y_midpoint_y_computed.find(data_y[vi]->getPointer(u_idx))
            == derivatives_y_midpoint_y_computed.end())
        {
            // Get the pointer to variable for derivative.
            Real* u = data_y[vi]->getPointer(u_idx);
            
            d_num_scratch_derivatives_midpoint_y_used++;
            // Declare container to store the derivative.
            if (allocate_scratch_derivatives_midpoint)
            {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
                TBOX_ASSERT(static_cast<int>(d_scratch_derivatives_midpoint_y.size()) == d_num_scratch_derivatives_midpoint_y_used - 1);
#endif
                d_scratch_derivatives_midpoint_y.push_back(HAMERS_SHARED_PTR<pdat::SideData<Real> >(
                    new pdat::SideData<Real>(
                        interior_box, 1, d_num_diff_ghosts, direction)));
            }
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
            else
            {
                TBOX_ASSERT(static_cast<int>(d_scratch_derivatives_midpoint_y.size()) >= d_num_scratch_derivatives_midpoint_y_used);
            }
#endif
            
            // Get the pointer to the derivative at midpoint.
            Real* dudy = d_scratch_derivatives_midpoint_y[d_num_scratch_derivatives_midpoint_y_used - 1]->getPointer(1, 0);
            
            /*
             * Get the ghost cell widths and ghost box dimensions of the variables.
             */
            
            const hier::IntVector& num_ghosts_derivative_midpoint = d_num_diff_ghosts;
            
            const hier::IntVector& ghostcell_dims_derivative_midpoint = diff_ghostcell_dims;
            
            const hier::IntVector& num_ghosts_data_node =
                data_y[vi]->getGhostCellWidth();
            
            const hier::IntVector& ghostcell_dims_data_node =
                data_y[vi]->getGhostBox().numberCells();
            
            hier::IntVector domain_lo(d_dim);
            hier::IntVector domain_dims(d_dim);
            
            domain_lo = hier::IntVector::getZero(d_dim);
            domain_dims = interior_dims;
            
            domain_lo[1] += (-d_num_diff_ghosts[1] + d_num_der_midpoint_ghosts[1]);
            domain_dims[1] += 1;
            domain_dims[1] += 2*(d_num_diff_ghosts[1] - d_num_der_midpoint_ghosts[1]);
            
            computeFirstDerivativesInYAtMidpointY(
                dudy,
                u,
                num_ghosts_derivative_midpoint,
                num_ghosts_data_node,
                ghostcell_dims_derivative_midpoint,
                ghostcell_dims_data_node,
                domain_lo,
                domain_dims,
                dx_1_inv);
            
            std::pair<Real*, HAMERS_SHARED_PTR<pdat::SideData<Real> > > derivative_pair(
                u,
                d_scratch_derivatives_midpoint_y[d_num_scratch_derivatives_midpoint_y_used - 1]);
            
            derivatives_y_midpoint_y_computed.insert(derivative_pair);
        }
        
        derivatives_y_midpoint_y.push_back(
            derivatives_y_midpoint_y_computed.find(data_y[vi]->getPointer(u_idx))->
            second);
    }
}


/*
 * Compute the derivatives in y-direction at midpoints.
 */
void
DiffusiveFluxReconstructorMidpoint::computeFirstDerivativesInYAtMidpointY(
    std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > >& derivatives_y_midpoint_y,
    std::map<Real*, HAMERS_SHARED_PTR<pdat::SideData<Real> > >& derivatives_y_midpoint_y_computed,
    const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > >& data_y,
    const std::vector<std::vector<int> >& data_component_idx_y,
    const hier::Patch& patch,
    const bool allocate_scratch_derivatives_midpoint)
{
    const int num_eqn = static_cast<int>(data_y.size());
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(d_num_eqn == num_eqn);
    TBOX_ASSERT(static_cast<int>(data_component_idx_y.size()) == num_eqn);
    
    for (int ei = 0; ei < num_eqn; ei++)
    {
        TBOX_ASSERT(static_cast<int>(data_y[ei].size()) ==
                    static_cast<int>(data_component_idx_y[ei].size()));
    }
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": DiffusiveFluxReconstructorMidpoint::"
            << "computeFirstDerivativesInYAtMidpointY()\n"
            << "There isn't y-direction for one-dimensional problem."
            << std::endl);
    }
#endif
    
    derivatives_y_midpoint_y.resize(num_eqn);
    
    for (int ei = 0; ei < num_eqn; ei++)
    {
        computeFirstDerivativesInYAtMidpointY(
            derivatives_y_midpoint_y[ei],
            derivatives_y_midpoint_y_computed,
            data_y[ei],
            data_component_idx_y[ei],
            patch,
            allocate_scratch_derivatives_midpoint);
    }
}


/*
 * Compute the derivatives in z-direction at midpoints.
 */
void
DiffusiveFluxReconstructorMidpoint::computeFirstDerivativesInZAtMidpointZ(
    std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& derivatives_z_midpoint_z,
    std::map<Real*, HAMERS_SHARED_PTR<pdat::SideData<Real> > >& derivatives_z_midpoint_z_computed,
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& data_z,
    const std::vector<int>& data_component_idx_z,
    const hier::Patch& patch,
    const bool allocate_scratch_derivatives_midpoint)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": DiffusiveFluxReconstructorMidpoint::"
            << "computeFirstDerivativesInZAtMidpointZ()\n"
            << "There isn't z-direction for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": DiffusiveFluxReconstructorMidpoint::"
            << "computeFirstDerivativesInZAtMidpointZ()\n"
            << "There isn't z-direction for two-dimensional problem."
            << std::endl);
    }
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // diffusive ghost cells.
    hier::Box diff_ghost_box = interior_box;
    diff_ghost_box.grow(d_num_diff_ghosts);
    const hier::IntVector diff_ghostcell_dims = diff_ghost_box.numberCells();
    
    // Get the grid spacing.
    const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
        HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
    const double* const dx = patch_geom->getDx();
    
    const Real dx_2_inv = Real(1)/Real(dx[2]);
    
    hier::IntVector direction = hier::IntVector::getZero(d_dim);
    direction[2] = 1;
    
    derivatives_z_midpoint_z.reserve(static_cast<int>(data_z.size()));
    
    for (int vi = 0; vi < static_cast<int>(data_z.size()); vi++)
    {
        // Get the index of variable for derivative.
        const int u_idx = data_component_idx_z[vi];
        
        if (derivatives_z_midpoint_z_computed.find(data_z[vi]->getPointer(u_idx))
            == derivatives_z_midpoint_z_computed.end())
        {
            // Get the pointer to variable for derivative.
            Real* u = data_z[vi]->getPointer(u_idx);
            
            d_num_scratch_derivatives_midpoint_z_used++;
            // Declare container to store the derivative.
            if (allocate_scratch_derivatives_midpoint)
            {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
                TBOX_ASSERT(static_cast<int>(d_scratch_derivatives_midpoint_z.size()) == d_num_scratch_derivatives_midpoint_z_used - 1);
#endif
                d_scratch_derivatives_midpoint_z.push_back(HAMERS_SHARED_PTR<pdat::SideData<Real> >(
                    new pdat::SideData<Real>(
                        interior_box, 1, d_num_diff_ghosts, direction)));
            }
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
            else
            {
                TBOX_ASSERT(static_cast<int>(d_scratch_derivatives_midpoint_z.size()) >= d_num_scratch_derivatives_midpoint_z_used);
            }
#endif
            
            // Get the pointer to the derivative at midpoint.
            Real* dudz = d_scratch_derivatives_midpoint_z[d_num_scratch_derivatives_midpoint_z_used - 1]->getPointer(2, 0);
            
            /*
             * Get the ghost cell widths and ghost box dimensions of the variables.
             */
            
            const hier::IntVector& num_ghosts_derivative_midpoint = d_num_diff_ghosts;
            
            const hier::IntVector& ghostcell_dims_derivative_midpoint = diff_ghostcell_dims;
            
            const hier::IntVector& num_ghosts_data_node =
                data_z[vi]->getGhostCellWidth();
            
            const hier::IntVector& ghostcell_dims_data_node =
                data_z[vi]->getGhostBox().numberCells();
            
            hier::IntVector domain_lo(d_dim);
            hier::IntVector domain_dims(d_dim);
            
            domain_lo = hier::IntVector::getZero(d_dim);
            domain_dims = interior_dims;
            
            domain_lo[2] += (-d_num_diff_ghosts[2] + d_num_der_midpoint_ghosts[2]);
            domain_dims[2] += 1;
            domain_dims[2] += 2*(d_num_diff_ghosts[2] - d_num_der_midpoint_ghosts[2]);
            
            computeFirstDerivativesInZAtMidpointZ(
                dudz,
                u,
                num_ghosts_derivative_midpoint,
                num_ghosts_data_node,
                ghostcell_dims_derivative_midpoint,
                ghostcell_dims_data_node,
                domain_lo,
                domain_dims,
                dx_2_inv);
            
            std::pair<Real*, HAMERS_SHARED_PTR<pdat::SideData<Real> > > derivative_pair(
                u,
                d_scratch_derivatives_midpoint_z[d_num_scratch_derivatives_midpoint_z_used - 1]);
            
            derivatives_z_midpoint_z_computed.insert(derivative_pair);
        }
        
        derivatives_z_midpoint_z.push_back(
            derivatives_z_midpoint_z_computed.find(data_z[vi]->getPointer(u_idx))->
            second);
    }
}


/*
 * Compute the derivatives in z-direction at midpoints.
 */
void
DiffusiveFluxReconstructorMidpoint::computeFirstDerivativesInZAtMidpointZ(
    std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > >& derivatives_z_midpoint_z,
    std::map<Real*, HAMERS_SHARED_PTR<pdat::SideData<Real> > >& derivatives_z_midpoint_z_computed,
    const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > >& data_z,
    const std::vector<std::vector<int> >& data_component_idx_z,
    const hier::Patch& patch,
    const bool allocate_scratch_derivatives_midpoint)
{
    const int num_eqn = static_cast<int>(data_z.size());
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(d_num_eqn == num_eqn);
    TBOX_ASSERT(static_cast<int>(data_component_idx_z.size()) == num_eqn);
    
    for (int ei = 0; ei < num_eqn; ei++)
    {
        TBOX_ASSERT(static_cast<int>(data_z[ei].size()) ==
                    static_cast<int>(data_component_idx_z[ei].size()));
    }
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": DiffusiveFluxReconstructorMidpoint::"
            << "computeFirstDerivativesInZAtMidpointZ()\n"
            << "There isn't z-direction for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": DiffusiveFluxReconstructorMidpoint::"
            << "computeFirstDerivativesInZAtMidpointZ()\n"
            << "There isn't z-direction for two-dimensional problem."
            << std::endl);
    }
#endif
    
    derivatives_z_midpoint_z.resize(num_eqn);
    
    for (int ei = 0; ei < num_eqn; ei++)
    {
        computeFirstDerivativesInZAtMidpointZ(
            derivatives_z_midpoint_z[ei],
            derivatives_z_midpoint_z_computed,
            data_z[ei],
            data_component_idx_z[ei],
            patch,
            allocate_scratch_derivatives_midpoint);
    }
}


/*
 * Compute the derivatives in x-direction at nodes.
 */
void
DiffusiveFluxReconstructorMidpoint::computeFirstDerivativesInXAtNode(
    std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& derivatives_x_node,
    std::map<Real*, HAMERS_SHARED_PTR<pdat::CellData<Real> > >& derivatives_x_node_computed,
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& data_x,
    const std::vector<int>& data_component_idx_x,
    const hier::Patch& patch,
    const bool allocate_scratch_derivatives_node)
{
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // diffusive ghost cells.
    hier::Box diff_ghost_box = interior_box;
    diff_ghost_box.grow(d_num_diff_ghosts);
    const hier::IntVector diff_ghostcell_dims = diff_ghost_box.numberCells();
    
    // Get the grid spacing.
    const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
        HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
    const double* const dx = patch_geom->getDx();
    
    const Real dx_0_inv = Real(1)/Real(dx[0]);
    
    derivatives_x_node.reserve(static_cast<int>(data_x.size()));
    
    for (int vi = 0; vi < static_cast<int>(data_x.size()); vi++)
    {
        // Get the index of variable for derivative.
        const int u_idx = data_component_idx_x[vi];
        
        if (derivatives_x_node_computed.find(data_x[vi]->getPointer(u_idx))
            == derivatives_x_node_computed.end())
        {
            // Get the pointer to variable for derivative.
            Real* u = data_x[vi]->getPointer(u_idx);
            
            d_num_scratch_derivatives_node_used++;
            // Declare container to store the derivative.
            if (allocate_scratch_derivatives_node)
            {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
                TBOX_ASSERT(static_cast<int>(d_scratch_derivatives_node.size()) == d_num_scratch_derivatives_node_used - 1);
#endif
                d_scratch_derivatives_node.push_back(HAMERS_SHARED_PTR<pdat::CellData<Real> >(
                    new pdat::CellData<Real>(
                        interior_box, 1, d_num_diff_ghosts)));
            }
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
            else
            {
                TBOX_ASSERT(static_cast<int>(d_scratch_derivatives_node.size()) >= d_num_scratch_derivatives_node_used);
            }
#endif
            
            // Get the pointer to the derivative.
            Real* dudx = d_scratch_derivatives_node[d_num_scratch_derivatives_node_used - 1]->getPointer(0);
            
            /*
             * Get the ghost cell widths and ghost box dimensions of the variables.
             */
            
            const hier::IntVector& num_ghosts_derivative_node = d_num_diff_ghosts;
            
            const hier::IntVector& ghostcell_dims_derivative_node = diff_ghostcell_dims;
            
            const hier::IntVector& num_ghosts_data_node =
                data_x[vi]->getGhostCellWidth();
            
            const hier::IntVector& ghostcell_dims_data_node =
                data_x[vi]->getGhostBox().numberCells();
            
            hier::IntVector domain_lo(d_dim);
            hier::IntVector domain_dims(d_dim);
            
            domain_lo = -d_num_diff_ghosts;
            domain_dims = interior_dims + d_num_diff_ghosts*2;
            
            domain_lo[0] += d_num_der_node_ghosts[0];
            domain_dims[0] -= 2*d_num_der_node_ghosts[0];
            
            computeFirstDerivativesInXAtNode(
                dudx,
                u,
                num_ghosts_derivative_node,
                num_ghosts_data_node,
                ghostcell_dims_derivative_node,
                ghostcell_dims_data_node,
                domain_lo,
                domain_dims,
                dx_0_inv);
            
            std::pair<Real*, HAMERS_SHARED_PTR<pdat::CellData<Real> > > derivative_pair(
                u,
                d_scratch_derivatives_node[d_num_scratch_derivatives_node_used - 1]);
            
            derivatives_x_node_computed.insert(derivative_pair);
        }
        
        derivatives_x_node.push_back(
            derivatives_x_node_computed.find(data_x[vi]->getPointer(u_idx))->
                second);
    }
}


/*
 * Compute the derivatives in x-direction at nodes.
 */
void
DiffusiveFluxReconstructorMidpoint::computeFirstDerivativesInXAtNode(
    std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > >& derivatives_x_node,
    std::map<Real*, HAMERS_SHARED_PTR<pdat::CellData<Real> > >& derivatives_x_node_computed,
    const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > >& data_x,
    const std::vector<std::vector<int> >& data_component_idx_x,
    const hier::Patch& patch,
    const bool allocate_scratch_derivatives_node)
{
    const int num_eqn = static_cast<int>(data_x.size());
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(d_num_eqn == num_eqn);
    TBOX_ASSERT(static_cast<int>(data_component_idx_x.size()) == num_eqn);
    
    for (int ei = 0; ei < num_eqn; ei++)
    {
        TBOX_ASSERT(static_cast<int>(data_x[ei].size()) ==
                    static_cast<int>(data_component_idx_x[ei].size()));
    }
#endif
    
    derivatives_x_node.resize(num_eqn);
    
    for (int ei = 0; ei < num_eqn; ei++)
    {
        computeFirstDerivativesInXAtNode(
            derivatives_x_node[ei],
            derivatives_x_node_computed,
            data_x[ei],
            data_component_idx_x[ei],
            patch,
            allocate_scratch_derivatives_node);
    }
}


/*
 * Compute the derivatives in y-direction at nodes.
 */
void
DiffusiveFluxReconstructorMidpoint::computeFirstDerivativesInYAtNode(
    std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& derivatives_y_node,
    std::map<Real*, HAMERS_SHARED_PTR<pdat::CellData<Real> > >& derivatives_y_node_computed,
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& data_y,
    const std::vector<int>& data_component_idx_y,
    const hier::Patch& patch,
    const bool allocate_scratch_derivatives_node)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": DiffusiveFluxReconstructorMidpoint::"
            << "computeFirstDerivativesInYAtNode()\n"
            << "There isn't y-direction for one-dimensional problem."
            << std::endl);
    }
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // diffusive ghost cells.
    hier::Box diff_ghost_box = interior_box;
    diff_ghost_box.grow(d_num_diff_ghosts);
    const hier::IntVector diff_ghostcell_dims = diff_ghost_box.numberCells();
    
    // Get the grid spacing.
    const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
        HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
    const double* const dx = patch_geom->getDx();
    
    const Real dx_1_inv = Real(1)/Real(dx[1]);
    
    derivatives_y_node.reserve(static_cast<int>(data_y.size()));
    
    for (int vi = 0; vi < static_cast<int>(data_y.size()); vi++)
    {
        // Get the index of variable for derivative.
        const int u_idx = data_component_idx_y[vi];
        
        if (derivatives_y_node_computed.find(data_y[vi]->getPointer(u_idx))
            == derivatives_y_node_computed.end())
        {
            // Get the pointer to variable for derivative.
            Real* u = data_y[vi]->getPointer(u_idx);
            
            d_num_scratch_derivatives_node_used++;
            // Declare container to store the derivative.
            if (allocate_scratch_derivatives_node)
            {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
                TBOX_ASSERT(static_cast<int>(d_scratch_derivatives_node.size()) == d_num_scratch_derivatives_node_used - 1);
#endif
                d_scratch_derivatives_node.push_back(HAMERS_SHARED_PTR<pdat::CellData<Real> >(
                    new pdat::CellData<Real>(
                        interior_box, 1, d_num_diff_ghosts)));
            }
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
            else
            {
                TBOX_ASSERT(static_cast<int>(d_scratch_derivatives_node.size()) >= d_num_scratch_derivatives_node_used);
            }
#endif
            
            // Get the pointer to the derivative.
            Real* dudy = d_scratch_derivatives_node[d_num_scratch_derivatives_node_used - 1]->getPointer(0);
            
            /*
             * Get the ghost cell widths and ghost box dimensions of the variables.
             */
            
            const hier::IntVector& num_ghosts_derivative_node = d_num_diff_ghosts;
            
            const hier::IntVector& ghostcell_dims_derivative_node = diff_ghostcell_dims;
            
            const hier::IntVector& num_ghosts_data_node =
                data_y[vi]->getGhostCellWidth();
            
            const hier::IntVector& ghostcell_dims_data_node =
                data_y[vi]->getGhostBox().numberCells();
            
            hier::IntVector domain_lo(d_dim);
            hier::IntVector domain_dims(d_dim);
            
            domain_lo = -d_num_diff_ghosts;
            domain_dims = interior_dims + d_num_diff_ghosts*2;
            
            domain_lo[1] += d_num_der_node_ghosts[1];
            domain_dims[1] -= 2*d_num_der_node_ghosts[1];
            
            computeFirstDerivativesInYAtNode(
                dudy,
                u,
                num_ghosts_derivative_node,
                num_ghosts_data_node,
                ghostcell_dims_derivative_node,
                ghostcell_dims_data_node,
                domain_lo,
                domain_dims,
                dx_1_inv);
            
            std::pair<Real*, HAMERS_SHARED_PTR<pdat::CellData<Real> > > derivative_pair(
                u,
                d_scratch_derivatives_node[d_num_scratch_derivatives_node_used - 1]);
            
            derivatives_y_node_computed.insert(derivative_pair);
        }
        
        derivatives_y_node.push_back(
            derivatives_y_node_computed.find(data_y[vi]->getPointer(u_idx))->
                second);
    }
}


/*
 * Compute the derivatives in y-direction at nodes.
 */
void
DiffusiveFluxReconstructorMidpoint::computeFirstDerivativesInYAtNode(
    std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > >& derivatives_y_node,
    std::map<Real*, HAMERS_SHARED_PTR<pdat::CellData<Real> > >& derivatives_y_node_computed,
    const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > >& data_y,
    const std::vector<std::vector<int> >& data_component_idx_y,
    const hier::Patch& patch,
    const bool allocate_scratch_derivatives_node)
{
    const int num_eqn = static_cast<int>(data_y.size());
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(d_num_eqn == num_eqn);
    TBOX_ASSERT(static_cast<int>(data_component_idx_y.size()) == num_eqn);
    
    for (int ei = 0; ei < num_eqn; ei++)
    {
        TBOX_ASSERT(static_cast<int>(data_y[ei].size()) ==
                    static_cast<int>(data_component_idx_y[ei].size()));
    }
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": DiffusiveFluxReconstructorMidpoint::"
            << "computeFirstDerivativesInYAtNode()\n"
            << "There isn't y-direction for one-dimensional problem."
            << std::endl);
    }
#endif
    
    derivatives_y_node.resize(num_eqn);
    
    for (int ei = 0; ei < num_eqn; ei++)
    {
        computeFirstDerivativesInYAtNode(
            derivatives_y_node[ei],
            derivatives_y_node_computed,
            data_y[ei],
            data_component_idx_y[ei],
            patch,
            allocate_scratch_derivatives_node);
    }
}


/*
 * Compute the derivatives in z-direction at nodes.
 */
void
DiffusiveFluxReconstructorMidpoint::computeFirstDerivativesInZAtNode(
    std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& derivatives_z_node,
    std::map<Real*, HAMERS_SHARED_PTR<pdat::CellData<Real> > >& derivatives_z_node_computed,
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& data_z,
    const std::vector<int>& data_component_idx_z,
    const hier::Patch& patch,
    const bool allocate_scratch_derivatives_node)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": DiffusiveFluxReconstructorMidpoint::"
            << "computeFirstDerivativesInZAtNode()\n"
            << "There isn't z-direction for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": DiffusiveFluxReconstructorMidpoint::"
            << "computeFirstDerivativesInZAtNode()\n"
            << "There isn't z-direction for two-dimensional problem."
            << std::endl);
    }
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // diffusive ghost cells.
    hier::Box diff_ghost_box = interior_box;
    diff_ghost_box.grow(d_num_diff_ghosts);
    const hier::IntVector diff_ghostcell_dims = diff_ghost_box.numberCells();
    
    // Get the grid spacing.
    const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
        HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
    const double* const dx = patch_geom->getDx();
    
    const Real dx_2_inv = Real(1)/Real(dx[2]);
    
    derivatives_z_node.reserve(static_cast<int>(data_z.size()));
    
    for (int vi = 0; vi < static_cast<int>(data_z.size()); vi++)
    {
        // Get the index of variable for derivative.
        const int u_idx = data_component_idx_z[vi];
        
        if (derivatives_z_node_computed.find(data_z[vi]->getPointer(u_idx))
            == derivatives_z_node_computed.end())
        {
            // Get the pointer to variable for derivative.
            Real* u = data_z[vi]->getPointer(u_idx);
            
            d_num_scratch_derivatives_node_used++;
            // Declare container to store the derivative.
            if (allocate_scratch_derivatives_node)
            {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
                TBOX_ASSERT(static_cast<int>(d_scratch_derivatives_node.size()) == d_num_scratch_derivatives_node_used - 1);
#endif
                d_scratch_derivatives_node.push_back(HAMERS_SHARED_PTR<pdat::CellData<Real> >(
                    new pdat::CellData<Real>(
                        interior_box, 1, d_num_diff_ghosts)));
            }
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
            else
            {
                TBOX_ASSERT(static_cast<int>(d_scratch_derivatives_node.size()) >= d_num_scratch_derivatives_node_used);
            }
#endif
            
            // Get the pointer to the derivative.
            Real* dudz = d_scratch_derivatives_node[d_num_scratch_derivatives_node_used - 1]->getPointer(0);
            
            /*
             * Get the ghost cell widths and ghost box dimensions of the variables.
             */
            
            const hier::IntVector& num_ghosts_derivative_node = d_num_diff_ghosts;
            
            const hier::IntVector& ghostcell_dims_derivative_node = diff_ghostcell_dims;
            
            const hier::IntVector& num_ghosts_data_node =
                data_z[vi]->getGhostCellWidth();
            
            const hier::IntVector& ghostcell_dims_data_node =
                data_z[vi]->getGhostBox().numberCells();
            
            hier::IntVector domain_lo(d_dim);
            hier::IntVector domain_dims(d_dim);
            
            domain_lo = -d_num_diff_ghosts;
            domain_dims = interior_dims + d_num_diff_ghosts*2;
            
            domain_lo[2] += d_num_der_node_ghosts[2];
            domain_dims[2] -= 2*d_num_der_node_ghosts[2];
            
            computeFirstDerivativesInZAtNode(
                dudz,
                u,
                num_ghosts_derivative_node,
                num_ghosts_data_node,
                ghostcell_dims_derivative_node,
                ghostcell_dims_data_node,
                domain_lo,
                domain_dims,
                dx_2_inv);
            
            std::pair<Real*, HAMERS_SHARED_PTR<pdat::CellData<Real> > > derivative_pair(
                u,
                d_scratch_derivatives_node[d_num_scratch_derivatives_node_used - 1]);
            
            derivatives_z_node_computed.insert(derivative_pair);
        }
        
        derivatives_z_node.push_back(
            derivatives_z_node_computed.find(data_z[vi]->getPointer(u_idx))->
                second);
    }
}


/*
 * Compute the derivatives in z-direction at nodes.
 */
void
DiffusiveFluxReconstructorMidpoint::computeFirstDerivativesInZAtNode(
    std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > >& derivatives_z_node,
    std::map<Real*, HAMERS_SHARED_PTR<pdat::CellData<Real> > >& derivatives_z_node_computed,
    const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > >& data_z,
    const std::vector<std::vector<int> >& data_component_idx_z,
    const hier::Patch& patch,
    const bool allocate_scratch_derivatives_node)
{
    const int num_eqn = static_cast<int>(data_z.size());
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(d_num_eqn == num_eqn);
    TBOX_ASSERT(static_cast<int>(data_component_idx_z.size()) == num_eqn);
    
    for (int ei = 0; ei < num_eqn; ei++)
    {
        TBOX_ASSERT(static_cast<int>(data_z[ei].size()) ==
                    static_cast<int>(data_component_idx_z[ei].size()));
    }
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": DiffusiveFluxReconstructorMidpoint::"
            << "computeFirstDerivativesInZAtNode()\n"
            << "There isn't z-direction for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": DiffusiveFluxReconstructorMidpoint::"
            << "computeFirstDerivativesInZAtNode()\n"
            << "There isn't z-direction for two-dimensional problem."
            << std::endl);
    }
#endif
    
    derivatives_z_node.resize(num_eqn);
    
    for (int ei = 0; ei < num_eqn; ei++)
    {
        computeFirstDerivativesInZAtNode(
            derivatives_z_node[ei],
            derivatives_z_node_computed,
            data_z[ei],
            data_component_idx_z[ei],
            patch,
            allocate_scratch_derivatives_node);
    }
}


/*
 * Interpolate the diffusivities from nodes to midpoints.
 */
void
DiffusiveFluxReconstructorMidpoint::interpolateDiffusivitiesFromNodeToMidpoint(
    std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& var_side_data_for_diffusivities,
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& var_cell_data_for_diffusivities,
    const std::vector<int>& var_cell_data_for_diffusivities_component_idx,
    const hier::Patch& patch,
    const bool allocate_scratch_diffusivities_midpoint)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(var_cell_data_for_diffusivities.size()) ==
        static_cast<int>(var_cell_data_for_diffusivities_component_idx.size()));
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // diffusive ghost cells.
    hier::Box diff_ghost_box = interior_box;
    diff_ghost_box.grow(d_num_diff_ghosts);
    const hier::IntVector diff_ghostcell_dims = diff_ghost_box.numberCells();
    
    var_side_data_for_diffusivities.reserve(static_cast<int>(var_cell_data_for_diffusivities.size()));
    
    for (int vi = 0; vi < static_cast<int>(var_cell_data_for_diffusivities.size()); vi++)
    {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
        TBOX_ASSERT(var_cell_data_for_diffusivities[vi]->getGhostCellWidth() >= d_num_diff_ghosts);
#endif
        
        // Get the pointer to variable for interpoation.
        Real* u_node = var_cell_data_for_diffusivities[vi]->getPointer(var_cell_data_for_diffusivities_component_idx[vi]);
        
        d_num_scratch_diffusivities_midpoint_used++;
        // Declare container to store the diffusivities.
        if (allocate_scratch_diffusivities_midpoint)
        {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
            TBOX_ASSERT(static_cast<int>(d_scratch_diffusivities_midpoint.size()) == d_num_scratch_diffusivities_midpoint_used - 1);
#endif
            d_scratch_diffusivities_midpoint.push_back(HAMERS_SHARED_PTR<pdat::SideData<Real> >(
                new pdat::SideData<Real>(
                    interior_box, 1, d_num_diff_ghosts)));
        }
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
        else
        {
            TBOX_ASSERT(static_cast<int>(d_scratch_diffusivities_midpoint.size()) >= d_num_scratch_diffusivities_midpoint_used);
        }
#endif
        
        /*
         * Get the ghost cell widths and ghost box dimensions of the variables.
         */
        
        const hier::IntVector& num_ghosts_data_midpoint = d_num_diff_ghosts;
        
        const hier::IntVector& ghostcell_dims_data_midpoint = diff_ghostcell_dims;
        
        const hier::IntVector& num_ghosts_data_node =
            var_cell_data_for_diffusivities[vi]->getGhostCellWidth();
        
        const hier::IntVector& ghostcell_dims_data_node =
            var_cell_data_for_diffusivities[vi]->getGhostBox().numberCells();
        
        hier::IntVector domain_lo(d_dim);
        hier::IntVector domain_dims(d_dim);
        
        if (d_dim == tbox::Dimension(1))
        {
            // Get the pointer to the diffusivity at midpoint.
            Real* u_midpoint_x =
                d_scratch_diffusivities_midpoint[d_num_scratch_diffusivities_midpoint_used - 1]->getPointer(0, 0);
            
            // Interpolation in x-direction.
            
            domain_lo = hier::IntVector::getZero(d_dim);
            domain_dims = interior_dims;
            
            domain_lo[0] += (-d_num_diff_ghosts[0] + d_num_interp_midpoint_ghosts[0]);
            domain_dims[0] += 1;
            domain_dims[0] += 2*(d_num_diff_ghosts[0] - d_num_interp_midpoint_ghosts[0]);
            
            interpolateDataFromNodeToMidpointX(
                u_midpoint_x,
                u_node,
                num_ghosts_data_midpoint,
                num_ghosts_data_node,
                ghostcell_dims_data_midpoint,
                ghostcell_dims_data_node,
                domain_lo,
                domain_dims);
        }
        else if (d_dim == tbox::Dimension(2))
        {
            // Get the pointers to the diffusivities at midpoints.
            Real* u_midpoint_x =
                d_scratch_diffusivities_midpoint[d_num_scratch_diffusivities_midpoint_used - 1]->getPointer(0, 0);
            Real* u_midpoint_y =
                d_scratch_diffusivities_midpoint[d_num_scratch_diffusivities_midpoint_used - 1]->getPointer(1, 0);
            
            // Interpolation in x-direction.
            
            domain_lo = hier::IntVector::getZero(d_dim);
            domain_dims = interior_dims;
            
            domain_lo[0] += (-d_num_diff_ghosts[0] + d_num_interp_midpoint_ghosts[0]);
            domain_dims[0] += 1;
            domain_dims[0] += 2*(d_num_diff_ghosts[0] - d_num_interp_midpoint_ghosts[0]);
            
            interpolateDataFromNodeToMidpointX(
                u_midpoint_x,
                u_node,
                num_ghosts_data_midpoint,
                num_ghosts_data_node,
                ghostcell_dims_data_midpoint,
                ghostcell_dims_data_node,
                domain_lo,
                domain_dims);
            
            // Interpolation in y-direction.
            
            domain_lo = hier::IntVector::getZero(d_dim);
            domain_dims = interior_dims;
            
            domain_lo[1] += (-d_num_diff_ghosts[1] + d_num_interp_midpoint_ghosts[1]);
            domain_dims[1] += 1;
            domain_dims[1] += 2*(d_num_diff_ghosts[1] - d_num_interp_midpoint_ghosts[1]);
            
            interpolateDataFromNodeToMidpointY(
                u_midpoint_y,
                u_node,
                num_ghosts_data_midpoint,
                num_ghosts_data_node,
                ghostcell_dims_data_midpoint,
                ghostcell_dims_data_node,
                domain_lo,
                domain_dims);
        }
        else if (d_dim == tbox::Dimension(3))
        {
            // Get the pointers to the diffusivities at midpoints.
            Real* u_midpoint_x =
                d_scratch_diffusivities_midpoint[d_num_scratch_diffusivities_midpoint_used - 1]->getPointer(0, 0);
            Real* u_midpoint_y =
                d_scratch_diffusivities_midpoint[d_num_scratch_diffusivities_midpoint_used - 1]->getPointer(1, 0);
            Real* u_midpoint_z =
                d_scratch_diffusivities_midpoint[d_num_scratch_diffusivities_midpoint_used - 1]->getPointer(2, 0);
            
            // Interpolation in x-direction.
            
            domain_lo = hier::IntVector::getZero(d_dim);
            domain_dims = interior_dims;
            
            domain_lo[0] += (-d_num_diff_ghosts[0] + d_num_interp_midpoint_ghosts[0]);
            domain_dims[0] += 1;
            domain_dims[0] += 2*(d_num_diff_ghosts[0] - d_num_interp_midpoint_ghosts[0]);
            
            interpolateDataFromNodeToMidpointX(
                u_midpoint_x,
                u_node,
                num_ghosts_data_midpoint,
                num_ghosts_data_node,
                ghostcell_dims_data_midpoint,
                ghostcell_dims_data_node,
                domain_lo,
                domain_dims);
            
            // Interpolation in y-direction.
            
            domain_lo = hier::IntVector::getZero(d_dim);
            domain_dims = interior_dims;
            
            domain_lo[1] += (-d_num_diff_ghosts[1] + d_num_interp_midpoint_ghosts[1]);
            domain_dims[1] += 1;
            domain_dims[1] += 2*(d_num_diff_ghosts[1] - d_num_interp_midpoint_ghosts[1]);
            
            interpolateDataFromNodeToMidpointY(
                u_midpoint_y,
                u_node,
                num_ghosts_data_midpoint,
                num_ghosts_data_node,
                ghostcell_dims_data_midpoint,
                ghostcell_dims_data_node,
                domain_lo,
                domain_dims);
            
            // Interpolation in z-direction.
            
            domain_lo = hier::IntVector::getZero(d_dim);
            domain_dims = interior_dims;
            
            domain_lo[2] += (-d_num_diff_ghosts[2] + d_num_interp_midpoint_ghosts[2]);
            domain_dims[2] += 1;
            domain_dims[2] += 2*(d_num_diff_ghosts[2] - d_num_interp_midpoint_ghosts[2]);
            
            interpolateDataFromNodeToMidpointZ(
                u_midpoint_z,
                u_node,
                num_ghosts_data_midpoint,
                num_ghosts_data_node,
                ghostcell_dims_data_midpoint,
                ghostcell_dims_data_node,
                domain_lo,
                domain_dims);
        }
        
        var_side_data_for_diffusivities.push_back(d_scratch_diffusivities_midpoint[d_num_scratch_diffusivities_midpoint_used - 1]);
    }
}


/*
 * Interpolate the derivatives from nodes to midpoints in x-direction.
 */
void
DiffusiveFluxReconstructorMidpoint::interpolateDerivativesFromNodeToMidpointX(
    std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& derivatives_midpoint_x,
    std::map<Real*, HAMERS_SHARED_PTR<pdat::SideData<Real> > >& derivatives_midpoint_x_computed,
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& derivative_node,
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& data,
    const std::vector<int>& data_component_idx,
    const hier::Patch& patch,
    const bool allocate_scratch_derivatives_midpoint)
{
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // diffusive ghost cells.
    hier::Box diff_ghost_box = interior_box;
    diff_ghost_box.grow(d_num_diff_ghosts);
    const hier::IntVector diff_ghostcell_dims = diff_ghost_box.numberCells();
    
    // Get the grid spacing.
    const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
        HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
    hier::IntVector direction = hier::IntVector::getZero(d_dim);
    direction[0] = 1;
    
    derivatives_midpoint_x.reserve(static_cast<int>(derivative_node.size()));
    
    for (int vi = 0; vi < static_cast<int>(derivative_node.size()); vi++)
    {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
        TBOX_ASSERT(derivative_node[vi]->getGhostCellWidth() == d_num_diff_ghosts);
#endif
        
        // Get the index of variable for derivative.
        const int u_idx = data_component_idx[vi];
        
        if (derivatives_midpoint_x_computed.find(data[vi]->getPointer(u_idx))
            == derivatives_midpoint_x_computed.end())
        {
            // Get the pointer to variable for derivative.
            Real* u = data[vi]->getPointer(u_idx);
            
            // Get the pointer to variable for interpoation.
            Real* der_node = derivative_node[vi]->getPointer(0);
            
            d_num_scratch_derivatives_midpoint_x_used++;
            // Declare container to store the derivative.
            if (allocate_scratch_derivatives_midpoint)
            {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
                TBOX_ASSERT(static_cast<int>(d_scratch_derivatives_midpoint_x.size()) == d_num_scratch_derivatives_midpoint_x_used - 1);
#endif
                d_scratch_derivatives_midpoint_x.push_back(HAMERS_SHARED_PTR<pdat::SideData<Real> >(
                    new pdat::SideData<Real>(
                        interior_box, 1, d_num_diff_ghosts, direction)));
            }
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
            else
            {
                TBOX_ASSERT(static_cast<int>(d_scratch_derivatives_midpoint_x.size()) >= d_num_scratch_derivatives_midpoint_x_used);
            }
#endif
            
            // Get the pointer to the derivative at midpoint.
            Real* der_midpoint_x = d_scratch_derivatives_midpoint_x[d_num_scratch_derivatives_midpoint_x_used - 1]->getPointer(0, 0);
            
            /*
             * Get the ghost cell widths and ghost box dimensions of the variables.
             */
            
            const hier::IntVector& num_ghosts_derivative_node = d_num_diff_ghosts;
            const hier::IntVector& num_ghosts_derivative_midpoint = d_num_diff_ghosts;
            
            const hier::IntVector& ghostcell_dims_derivative_node = diff_ghostcell_dims;
            const hier::IntVector& ghostcell_dims_derivative_midpoint = diff_ghostcell_dims;
            
            hier::IntVector domain_lo(d_dim);
            hier::IntVector domain_dims(d_dim);
            
            domain_lo = hier::IntVector::getZero(d_dim);
            domain_dims = interior_dims;
            
            domain_lo[0] += (-d_num_diff_ghosts[0] + d_num_interp_midpoint_ghosts[0]);
            domain_dims[0] += 1;
            domain_dims[0] += 2*(d_num_diff_ghosts[0] - d_num_interp_midpoint_ghosts[0]);
            
            interpolateDataFromNodeToMidpointX(
                der_midpoint_x,
                der_node,
                num_ghosts_derivative_midpoint,
                num_ghosts_derivative_node,
                ghostcell_dims_derivative_midpoint,
                ghostcell_dims_derivative_node,
                domain_lo,
                domain_dims);
            
            std::pair<Real*, HAMERS_SHARED_PTR<pdat::SideData<Real> > > derivative_pair(
                u,
                d_scratch_derivatives_midpoint_x[d_num_scratch_derivatives_midpoint_x_used - 1]);
            
            derivatives_midpoint_x_computed.insert(derivative_pair);
        }
        
        derivatives_midpoint_x.push_back(
            derivatives_midpoint_x_computed.find(data[vi]->getPointer(u_idx))->
                second);
    }
}


/*
 * Interpolate the derivatives from nodes to midpoints in x-direction.
 */
void
DiffusiveFluxReconstructorMidpoint::interpolateDerivativesFromNodeToMidpointX(
    std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > >& derivatives_midpoint_x,
    std::map<Real*, HAMERS_SHARED_PTR<pdat::SideData<Real> > >& derivatives_midpoint_x_computed,
    const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > >& derivative_node,
    const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > >& data,
    const std::vector<std::vector<int> >& data_component_idx,
    const hier::Patch& patch,
    const bool allocate_scratch_derivatives_midpoint)
{
    const int num_eqn = static_cast<int>(derivative_node.size());
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(d_num_eqn == num_eqn);
    TBOX_ASSERT(static_cast<int>(data.size()) == num_eqn);
    TBOX_ASSERT(static_cast<int>(data_component_idx.size()) == num_eqn);
#endif
    
    derivatives_midpoint_x.resize(num_eqn);
    
    for (int ei = 0; ei < num_eqn; ei++)
    {
        interpolateDerivativesFromNodeToMidpointX(
            derivatives_midpoint_x[ei],
            derivatives_midpoint_x_computed,
            derivative_node[ei],
            data[ei],
            data_component_idx[ei],
            patch,
            allocate_scratch_derivatives_midpoint);
    }
}


/*
 * Interpolate the derivatives from nodes to midpoints in y-direction.
 */
void
DiffusiveFluxReconstructorMidpoint::interpolateDerivativesFromNodeToMidpointY(
    std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& derivatives_midpoint_y,
    std::map<Real*, HAMERS_SHARED_PTR<pdat::SideData<Real> > >& derivatives_midpoint_y_computed,
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& derivative_node,
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& data,
    const std::vector<int>& data_component_idx,
    const hier::Patch& patch,
    const bool allocate_scratch_derivatives_midpoint)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": DiffusiveFluxReconstructorMidpoint::"
            << "interpolateDerivativesFromNodeToMidpointY()\n"
            << "There isn't y-direction for one-dimensional problem."
            << std::endl);
    }
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // diffusive ghost cells.
    hier::Box diff_ghost_box = interior_box;
    diff_ghost_box.grow(d_num_diff_ghosts);
    const hier::IntVector diff_ghostcell_dims = diff_ghost_box.numberCells();
    
    // Get the grid spacing.
    const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
        HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
    hier::IntVector direction = hier::IntVector::getZero(d_dim);
    direction[1] = 1;
    
    derivatives_midpoint_y.reserve(static_cast<int>(derivative_node.size()));
    
    for (int vi = 0; vi < static_cast<int>(derivative_node.size()); vi++)
    {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
        TBOX_ASSERT(derivative_node[vi]->getGhostCellWidth() == d_num_diff_ghosts);
#endif
        
        // Get the index of variable for derivative.
        const int u_idx = data_component_idx[vi];
        
        if (derivatives_midpoint_y_computed.find(data[vi]->getPointer(u_idx))
            == derivatives_midpoint_y_computed.end())
        {
            // Get the pointer to variable for derivative.
            Real* u = data[vi]->getPointer(u_idx);
            
            // Get the pointer to variable for interpoation.
            Real* der_node = derivative_node[vi]->getPointer(0);
            
            d_num_scratch_derivatives_midpoint_y_used++;
            // Declare container to store the derivative.
            if (allocate_scratch_derivatives_midpoint)
            {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
                TBOX_ASSERT(static_cast<int>(d_scratch_derivatives_midpoint_y.size()) == d_num_scratch_derivatives_midpoint_y_used - 1);
#endif
                d_scratch_derivatives_midpoint_y.push_back(HAMERS_SHARED_PTR<pdat::SideData<Real> >(
                    new pdat::SideData<Real>(
                        interior_box, 1, d_num_diff_ghosts, direction)));
            }
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
            else
            {
                TBOX_ASSERT(static_cast<int>(d_scratch_derivatives_midpoint_y.size()) >= d_num_scratch_derivatives_midpoint_y_used);
            }
#endif
            
            // Get the pointer to the derivative at midpoint.
            Real* der_midpoint_y = d_scratch_derivatives_midpoint_y[d_num_scratch_derivatives_midpoint_y_used - 1]->getPointer(1, 0);
            
            /*
             * Get the ghost cell widths and ghost box dimensions of the variables.
             */
            
            const hier::IntVector& num_ghosts_derivative_node = d_num_diff_ghosts;
            const hier::IntVector& num_ghosts_derivative_midpoint = d_num_diff_ghosts;
            
            const hier::IntVector& ghostcell_dims_derivative_node = diff_ghostcell_dims;
            const hier::IntVector& ghostcell_dims_derivative_midpoint = diff_ghostcell_dims;
            
            hier::IntVector domain_lo(d_dim);
            hier::IntVector domain_dims(d_dim);
            
            domain_lo = hier::IntVector::getZero(d_dim);
            domain_dims = interior_dims;
            
            domain_lo[1] += (-d_num_diff_ghosts[1] + d_num_interp_midpoint_ghosts[1]);
            domain_dims[1] += 1;
            domain_dims[1] += 2*(d_num_diff_ghosts[1] - d_num_interp_midpoint_ghosts[1]);
            
            interpolateDataFromNodeToMidpointY(
                der_midpoint_y,
                der_node,
                num_ghosts_derivative_midpoint,
                num_ghosts_derivative_node,
                ghostcell_dims_derivative_midpoint,
                ghostcell_dims_derivative_node,
                domain_lo,
                domain_dims);
            
            std::pair<Real*, HAMERS_SHARED_PTR<pdat::SideData<Real> > > derivative_pair(
                u,
                d_scratch_derivatives_midpoint_y[d_num_scratch_derivatives_midpoint_y_used - 1]);
            
            derivatives_midpoint_y_computed.insert(derivative_pair);
        }
        
        derivatives_midpoint_y.push_back(
            derivatives_midpoint_y_computed.find(data[vi]->getPointer(u_idx))->
                second);
    }
}


/*
 * Interpolate the derivatives from nodes to midpoints in y-direction.
 */
void
DiffusiveFluxReconstructorMidpoint::interpolateDerivativesFromNodeToMidpointY(
    std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > >& derivatives_midpoint_y,
    std::map<Real*, HAMERS_SHARED_PTR<pdat::SideData<Real> > >& derivatives_midpoint_y_computed,
    const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > >& derivative_node,
    const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > >& data,
    const std::vector<std::vector<int> >& data_component_idx,
    const hier::Patch& patch,
    const bool allocate_scratch_derivatives_midpoint)
{
    const int num_eqn = static_cast<int>(derivative_node.size());
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(d_num_eqn == num_eqn);
    TBOX_ASSERT(static_cast<int>(data.size()) == num_eqn);
    TBOX_ASSERT(static_cast<int>(data_component_idx.size()) == num_eqn);
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": DiffusiveFluxReconstructorMidpoint::"
            << "interpolateDerivativesFromNodeToMidpointY()\n"
            << "There isn't y-direction for one-dimensional problem."
            << std::endl);
    }
#endif
    
    derivatives_midpoint_y.resize(num_eqn);
    
    for (int ei = 0; ei < num_eqn; ei++)
    {
        interpolateDerivativesFromNodeToMidpointY(
            derivatives_midpoint_y[ei],
            derivatives_midpoint_y_computed,
            derivative_node[ei],
            data[ei],
            data_component_idx[ei],
            patch,
            allocate_scratch_derivatives_midpoint);
    }
}


/*
 * Interpolate the derivatives from nodes to midpoints in z-direction.
 */
void
DiffusiveFluxReconstructorMidpoint::interpolateDerivativesFromNodeToMidpointZ(
    std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& derivatives_midpoint_z,
    std::map<Real*, HAMERS_SHARED_PTR<pdat::SideData<Real> > >& derivatives_midpoint_z_computed,
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& derivative_node,
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& data,
    const std::vector<int>& data_component_idx,
    const hier::Patch& patch,
    const bool allocate_scratch_derivatives_midpoint)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": DiffusiveFluxReconstructorMidpoint::"
            << "interpolateDerivativesFromNodeToMidpointZ()\n"
            << "There isn't z-direction for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": DiffusiveFluxReconstructorMidpoint::"
            << "interpolateDerivativesFromNodeToMidpointZ()\n"
            << "There isn't z-direction for two-dimensional problem."
            << std::endl);
    }
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // diffusive ghost cells.
    hier::Box diff_ghost_box = interior_box;
    diff_ghost_box.grow(d_num_diff_ghosts);
    const hier::IntVector diff_ghostcell_dims = diff_ghost_box.numberCells();
    
    // Get the grid spacing.
    const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
        HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
    hier::IntVector direction = hier::IntVector::getZero(d_dim);
    direction[2] = 1;
    
    derivatives_midpoint_z.reserve(static_cast<int>(derivative_node.size()));
    
    for (int vi = 0; vi < static_cast<int>(derivative_node.size()); vi++)
    {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
        TBOX_ASSERT(derivative_node[vi]->getGhostCellWidth() == d_num_diff_ghosts);
#endif
        
        // Get the index of variable for derivative.
        const int u_idx = data_component_idx[vi];
        
        if (derivatives_midpoint_z_computed.find(data[vi]->getPointer(u_idx))
            == derivatives_midpoint_z_computed.end())
        {
            // Get the pointer to variable for derivative.
            Real* u = data[vi]->getPointer(u_idx);
            
            // Get the pointer to variable for interpoation.
            Real* der_node = derivative_node[vi]->getPointer(0);
            
            d_num_scratch_derivatives_midpoint_z_used++;
            // Declare container to store the derivative.
            if (allocate_scratch_derivatives_midpoint)
            {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
                TBOX_ASSERT(static_cast<int>(d_scratch_derivatives_midpoint_z.size()) == d_num_scratch_derivatives_midpoint_z_used - 1);
#endif
                d_scratch_derivatives_midpoint_z.push_back(HAMERS_SHARED_PTR<pdat::SideData<Real> >(
                    new pdat::SideData<Real>(
                        interior_box, 1, d_num_diff_ghosts, direction)));
            }
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
            else
            {
                TBOX_ASSERT(static_cast<int>(d_scratch_derivatives_midpoint_z.size()) >= d_num_scratch_derivatives_midpoint_z_used);
            }
#endif
            
            // Get the pointer to the derivative at midpoint.
            Real* der_midpoint_z = d_scratch_derivatives_midpoint_z[d_num_scratch_derivatives_midpoint_z_used - 1]->getPointer(2, 0);
            
            /*
             * Get the ghost cell widths and ghost box dimensions of the variables.
             */
            
            const hier::IntVector& num_ghosts_derivative_node = d_num_diff_ghosts;
            const hier::IntVector& num_ghosts_derivative_midpoint = d_num_diff_ghosts;
            
            const hier::IntVector& ghostcell_dims_derivative_node = diff_ghostcell_dims;
            const hier::IntVector& ghostcell_dims_derivative_midpoint = diff_ghostcell_dims;
            
            hier::IntVector domain_lo(d_dim);
            hier::IntVector domain_dims(d_dim);
            
            domain_lo = hier::IntVector::getZero(d_dim);
            domain_dims = interior_dims;
            
            domain_lo[2] += (-d_num_diff_ghosts[2] + d_num_interp_midpoint_ghosts[2]);
            domain_dims[2] += 1;
            domain_dims[2] += 2*(d_num_diff_ghosts[2] - d_num_interp_midpoint_ghosts[2]);
            
            interpolateDataFromNodeToMidpointZ(
                der_midpoint_z,
                der_node,
                num_ghosts_derivative_midpoint,
                num_ghosts_derivative_node,
                ghostcell_dims_derivative_midpoint,
                ghostcell_dims_derivative_node,
                domain_lo,
                domain_dims);
            
            std::pair<Real*, HAMERS_SHARED_PTR<pdat::SideData<Real> > > derivative_pair(
                u,
                d_scratch_derivatives_midpoint_z[d_num_scratch_derivatives_midpoint_z_used - 1]);
            
            derivatives_midpoint_z_computed.insert(derivative_pair);
        }
        
        derivatives_midpoint_z.push_back(
            derivatives_midpoint_z_computed.find(data[vi]->getPointer(u_idx))->
                second);
    }
}


/*
 * Interpolate the derivatives from nodes to midpoints in z-direction.
 */
void
DiffusiveFluxReconstructorMidpoint::interpolateDerivativesFromNodeToMidpointZ(
    std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > >& derivatives_midpoint_z,
    std::map<Real*, HAMERS_SHARED_PTR<pdat::SideData<Real> > >& derivatives_midpoint_z_computed,
    const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > >& derivative_node,
    const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > >& data,
    const std::vector<std::vector<int> >& data_component_idx,
    const hier::Patch& patch,
    const bool allocate_scratch_derivatives_midpoint)
{
    const int num_eqn = static_cast<int>(derivative_node.size());
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(d_num_eqn == num_eqn);
    TBOX_ASSERT(static_cast<int>(data.size()) == num_eqn);
    TBOX_ASSERT(static_cast<int>(data_component_idx.size()) == num_eqn);
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": DiffusiveFluxReconstructorMidpoint::"
            << "interpolateDerivativesFromNodeToMidpointZ()\n"
            << "There isn't z-direction for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": DiffusiveFluxReconstructorMidpoint::"
            << "interpolateDerivativesFromNodeToMidpointZ()\n"
            << "There isn't z-direction for two-dimensional problem."
            << std::endl);
    }
#endif
    
    derivatives_midpoint_z.resize(num_eqn);
    
    for (int ei = 0; ei < num_eqn; ei++)
    {
        interpolateDerivativesFromNodeToMidpointZ(
            derivatives_midpoint_z[ei],
            derivatives_midpoint_z_computed,
            derivative_node[ei],
            data[ei],
            data_component_idx[ei],
            patch,
            allocate_scratch_derivatives_midpoint);
    }
}


/*
 * Reconstruct the flux using flux at midpoints in x-direction.
 */
void
DiffusiveFluxReconstructorMidpoint::reconstructFluxX(
    HAMERS_SHARED_PTR<pdat::SideData<Real> >& diffusive_flux,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& diffusive_flux_midpoint,
    const hier::Patch& patch,
    const double dt) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(diffusive_flux->getDepth() == d_num_eqn);
    TBOX_ASSERT(diffusive_flux_midpoint->getDepth() == d_num_eqn);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // diffusive ghost cells.
    hier::Box diff_ghost_box = interior_box;
    diff_ghost_box.grow(d_num_diff_ghosts);
    const hier::IntVector ghostcell_dims_flux_midpoint = diff_ghost_box.numberCells();
    
    const hier::IntVector& num_ghosts_flux_midpoint = d_num_diff_ghosts;
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    domain_lo = hier::IntVector::getZero(d_dim);
    domain_dims = interior_dims;
    domain_dims[0]++;
    
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        Real* F_face_x     = diffusive_flux->getPointer(0, ei);
        Real* F_midpoint_x = diffusive_flux_midpoint->getPointer(0, ei);
        
        reconstructFluxX(
            F_face_x,
            F_midpoint_x,
            num_ghosts_flux_midpoint,
            ghostcell_dims_flux_midpoint,
            domain_lo,
            domain_dims,
            interior_dims,
            dt);
    }
}


/*
 * Reconstruct the flux using flux at midpoints in y-direction.
 */
void
DiffusiveFluxReconstructorMidpoint::reconstructFluxY(
    HAMERS_SHARED_PTR<pdat::SideData<Real> >& diffusive_flux,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& diffusive_flux_midpoint,
    const hier::Patch& patch,
    const double dt) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(diffusive_flux->getDepth() == d_num_eqn);
    TBOX_ASSERT(diffusive_flux_midpoint->getDepth() == d_num_eqn);
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": DiffusiveFluxReconstructorMidpoint::"
            << "reconstructFluxY()\n"
            << "There isn't y-direction for one-dimensional problem."
            << std::endl);
    }
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // diffusive ghost cells.
    hier::Box diff_ghost_box = interior_box;
    diff_ghost_box.grow(d_num_diff_ghosts);
    const hier::IntVector ghostcell_dims_flux_midpoint = diff_ghost_box.numberCells();
    
    const hier::IntVector& num_ghosts_flux_midpoint = d_num_diff_ghosts;
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    domain_lo = hier::IntVector::getZero(d_dim);
    domain_dims = interior_dims;
    domain_dims[1]++;
    
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        Real* F_face_y     = diffusive_flux->getPointer(1, ei);
        Real* F_midpoint_y = diffusive_flux_midpoint->getPointer(1, ei);
        
        reconstructFluxY(
            F_face_y,
            F_midpoint_y,
            num_ghosts_flux_midpoint,
            ghostcell_dims_flux_midpoint,
            domain_lo,
            domain_dims,
            interior_dims,
            dt);
    }
}


/*
 * Reconstruct the flux using flux at midpoints in z-direction.
 */
void
DiffusiveFluxReconstructorMidpoint::reconstructFluxZ(
    HAMERS_SHARED_PTR<pdat::SideData<Real> >& diffusive_flux,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& diffusive_flux_midpoint,
    const hier::Patch& patch,
    const double dt) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(diffusive_flux->getDepth() == d_num_eqn);
    TBOX_ASSERT(diffusive_flux_midpoint->getDepth() == d_num_eqn);
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": DiffusiveFluxReconstructorMidpoint::"
            << "reconstructFluxZ()\n"
            << "There isn't z-direction for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": DiffusiveFluxReconstructorMidpoint::"
            << "reconstructFluxZ()\n"
            << "There isn't z-direction for two-dimensional problem."
            << std::endl);
    }
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // diffusive ghost cells.
    hier::Box diff_ghost_box = interior_box;
    diff_ghost_box.grow(d_num_diff_ghosts);
    const hier::IntVector ghostcell_dims_flux_midpoint = diff_ghost_box.numberCells();
    
    const hier::IntVector& num_ghosts_flux_midpoint = d_num_diff_ghosts;
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    domain_lo = hier::IntVector::getZero(d_dim);
    domain_dims = interior_dims;
    domain_dims[2]++;
    
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        Real* F_face_z     = diffusive_flux->getPointer(2, ei);
        Real* F_midpoint_z = diffusive_flux_midpoint->getPointer(2, ei);
        
        reconstructFluxZ(
            F_face_z,
            F_midpoint_z,
            num_ghosts_flux_midpoint,
            ghostcell_dims_flux_midpoint,
            domain_lo,
            domain_dims,
            interior_dims,
            dt);
    }
}
