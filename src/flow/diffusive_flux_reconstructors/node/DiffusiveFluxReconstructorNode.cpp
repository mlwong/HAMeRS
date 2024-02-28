#include "flow/diffusive_flux_reconstructors/node/DiffusiveFluxReconstructorNode.hpp"

#include "SAMRAI/geom/CartesianPatchGeometry.h"

#include <map>

DiffusiveFluxReconstructorNode::DiffusiveFluxReconstructorNode(
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
        d_num_der_node_ghosts(hier::IntVector::getZero(d_dim)),
        d_num_flux_reconstruct_ghosts(hier::IntVector::getZero(d_dim)),
        d_num_scratch_derivatives_node_used(0)
{
}


/*
 * Compute the diffusive flux on a patch.
 */
void
DiffusiveFluxReconstructorNode::computeDiffusiveFluxOnPatch(
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
    
    if (diffusive_flux_utilities->useSubgridScaleModel())
    {
        TBOX_ERROR(d_object_name
            << ": computeDiffusiveFluxOnPatch::"
            << "computeDiffusiveFluxOnPatch()\n"
            << "Subgrid scale model is not implemented."
            << std::endl);
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
        
        diffusive_flux_utilities->registerDerivedVariablesForDiffusiveFluxes(d_num_diff_ghosts);
        
        diffusive_flux_utilities->allocateMemoryForDerivedCellData();
        
        diffusive_flux_utilities->computeDerivedCellData();
        
        /*
         * Delcare containers for computing fluxes in different directions.
         */
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > > var_data_x;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > > diffusivities_data_x;
        
        std::vector<std::vector<int> > var_component_idx_x;
        
        std::vector<std::vector<int> > diffusivities_component_idx_x;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > > derivatives_x;
        
        std::map<Real*, HAMERS_SHARED_PTR<pdat::CellData<Real> > > derivatives_x_computed;
        
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > diffusive_flux_node(1);
        diffusive_flux_node[0].reset(new pdat::CellData<Real>(interior_box, d_num_eqn, d_num_diff_ghosts));
        
        std::vector<Real*> F_node_x;
        F_node_x.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            F_node_x.push_back(diffusive_flux_node[0]->getPointer(ei));
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
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxDiffusivities(
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
        
        computeFirstDerivativesInX(
            derivatives_x,
            derivatives_x_computed,
            var_data_x,
            var_component_idx_x,
            patch,
            allocate_scratch_data_containers);
        
        /*
         * Compute diffusive flux in x-direction at nodes.
         */
        
        diffusive_flux_node[0]->fillAll(Real(0));
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ASSERT(static_cast<int>(derivatives_x[ei].size()) ==
                        static_cast<int>(diffusivities_data_x[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_x[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_x[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_x[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_x[ei][vi];
                
                // Get the pointer to diffusivity.
                Real* mu = diffusivities_data_x[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivative.
                Real* dudx = derivatives_x[ei][vi]->getPointer(0);
                
                /*
                 * Get the sub-ghost cell width of the diffusivity.
                 */
                
                hier::IntVector num_subghosts_diffusivity =
                    diffusivities_data_x[ei][vi]->getGhostCellWidth();
                
                const int num_subghosts_0_diffusivity = num_subghosts_diffusivity[0];
                
                HAMERS_PRAGMA_SIMD
                for (int i = -num_subghosts_0_diffusivity + num_flux_reconstruct_ghosts_0;
                    i < interior_dim_0 + num_subghosts_0_diffusivity - num_flux_reconstruct_ghosts_0;
                    i++)
                {
                    // Compute the linear indices.
                    const int idx = i + num_diff_ghosts_0;
                    const int idx_diffusivity = i + num_subghosts_0_diffusivity;
                    
                    F_node_x[ei][idx] += mu[idx_diffusivity]*dudx[idx];
                }
            }
        }
        
        /*
         * Reconstruct the flux in x-direction.
         */
        
        reconstructFluxX(
            diffusive_flux,
            diffusive_flux_node[0],
            patch,
            dt);
        
        var_data_x.clear();
        diffusivities_data_x.clear();
        var_component_idx_x.clear();
        diffusivities_component_idx_x.clear();
        derivatives_x.clear();
        
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
        
        diffusive_flux_utilities->registerDerivedVariablesForDiffusiveFluxes(d_num_diff_ghosts);
        
        diffusive_flux_utilities->allocateMemoryForDerivedCellData();
        
        diffusive_flux_utilities->computeDerivedCellData();
        
        /*
         * Delcare containers for computing fluxes in different directions.
         */
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > > var_data_x;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > > var_data_y;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > > diffusivities_data_x;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > > diffusivities_data_y;
        
        std::vector<std::vector<int> > var_component_idx_x;
        std::vector<std::vector<int> > var_component_idx_y;
        
        std::vector<std::vector<int> > diffusivities_component_idx_x;
        std::vector<std::vector<int> > diffusivities_component_idx_y;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > > derivatives_x;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > > derivatives_y;
        
        std::map<Real*, HAMERS_SHARED_PTR<pdat::CellData<Real> > > derivatives_x_computed;
        std::map<Real*, HAMERS_SHARED_PTR<pdat::CellData<Real> > > derivatives_y_computed;
        
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > diffusive_flux_node(2);
        diffusive_flux_node[0].reset(new pdat::CellData<Real>(interior_box, d_num_eqn, d_num_diff_ghosts));
        diffusive_flux_node[1].reset(new pdat::CellData<Real>(interior_box, d_num_eqn, d_num_diff_ghosts));
        
        std::vector<Real*> F_node_x;
        std::vector<Real*> F_node_y;
        F_node_x.reserve(d_num_eqn);
        F_node_y.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            F_node_x.push_back(diffusive_flux_node[0]->getPointer(ei));
            F_node_y.push_back(diffusive_flux_node[1]->getPointer(ei));
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
        
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxVariablesForDerivative(
            var_data_y,
            var_component_idx_y,
            DIRECTION::X_DIRECTION,
            DIRECTION::Y_DIRECTION);
        
        // Get the diffusivities in the diffusive flux.
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxDiffusivities(
            diffusivities_data_x,
            diffusivities_component_idx_x,
            DIRECTION::X_DIRECTION,
            DIRECTION::X_DIRECTION);
        
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxDiffusivities(
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
        
        computeFirstDerivativesInX(
            derivatives_x,
            derivatives_x_computed,
            var_data_x,
            var_component_idx_x,
            patch,
            allocate_scratch_data_containers);
        
        /*
         * Compute the derivatives in y-direction for diffusive flux in x-direction.
         */
        
        computeFirstDerivativesInY(
            derivatives_y,
            derivatives_y_computed,
            var_data_y,
            var_component_idx_y,
            patch,
            allocate_scratch_data_containers);
        
        /*
         * Compute diffusive flux in x-direction at nodes.
         */
        
        diffusive_flux_node[0]->fillAll(Real(0));
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ASSERT(static_cast<int>(derivatives_x[ei].size()) ==
                        static_cast<int>(diffusivities_data_x[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_x[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_x[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_x[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_x[ei][vi];
                
                // Get the pointer to diffusivity.
                Real* mu = diffusivities_data_x[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivative.
                Real* dudx = derivatives_x[ei][vi]->getPointer(0);
                
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
                    for (int i = -3; i < interior_dim_0 + 3; i++)
                    {
                        // Compute the linear indices.
                        const int idx = (i + num_diff_ghosts_0) +
                            (j + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                        
                        const int idx_diffusivity = (i + num_subghosts_0_diffusivity) +
                            (j + num_subghosts_1_diffusivity)*subghostcell_dim_0_diffusivity;
                        
                        F_node_x[ei][idx] += mu[idx_diffusivity]*dudx[idx];
                    }
                }
            }
            
            TBOX_ASSERT(static_cast<int>(derivatives_y[ei].size()) ==
                        static_cast<int>(diffusivities_data_y[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_y[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_y[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_y[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_y[ei][vi];
                
                // Get the pointer to diffusivity.
                Real* mu = diffusivities_data_y[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivative.
                Real* dudy = derivatives_y[ei][vi]->getPointer(0);
                
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
                        i < interior_dim_0 + num_subghosts_0_diffusivity - num_flux_reconstruct_ghosts_0;
                        i++)
                    {
                        // Compute the linear indices.
                        const int idx = (i + num_diff_ghosts_0) +
                            (j + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                        
                        const int idx_diffusivity = (i + num_subghosts_0_diffusivity) +
                            (j + num_subghosts_1_diffusivity)*subghostcell_dim_0_diffusivity;
                        
                        F_node_x[ei][idx] += mu[idx_diffusivity]*dudy[idx];
                    }
                }
            }
        }
        
        /*
         * Reconstruct the flux in x-direction.
         */
        
        reconstructFluxX(
            diffusive_flux,
            diffusive_flux_node[0],
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
        
        derivatives_x.clear();
        derivatives_y.clear();
        
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
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxDiffusivities(
            diffusivities_data_x,
            diffusivities_component_idx_x,
            DIRECTION::Y_DIRECTION,
            DIRECTION::X_DIRECTION);
        
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxDiffusivities(
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
        
        computeFirstDerivativesInX(
            derivatives_x,
            derivatives_x_computed,
            var_data_x,
            var_component_idx_x,
            patch,
            allocate_scratch_data_containers);
        
        /*
         * Compute the derivatives in y-direction for diffusive flux in y-direction.
         */
        
        computeFirstDerivativesInY(
            derivatives_y,
            derivatives_y_computed,
            var_data_y,
            var_component_idx_y,
            patch,
            allocate_scratch_data_containers);
        
        /*
         * Compute diffusive flux in y-direction at nodes.
         */
        
        diffusive_flux_node[1]->fillAll(Real(0));
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ASSERT(static_cast<int>(derivatives_x[ei].size()) ==
                        static_cast<int>(diffusivities_data_x[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_x[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_x[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_x[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_x[ei][vi];
                
                // Get the pointer to diffusivity.
                Real* mu = diffusivities_data_x[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivative.
                Real* dudx = derivatives_x[ei][vi]->getPointer(0);
                
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
                
                for (int j = -3; j < interior_dim_1 + 3; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx = (i + num_diff_ghosts_0) +
                            (j + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                        
                        const int idx_diffusivity = (i + num_subghosts_0_diffusivity) +
                            (j + num_subghosts_1_diffusivity)*subghostcell_dim_0_diffusivity;
                        
                        F_node_y[ei][idx] += mu[idx_diffusivity]*dudx[idx];
                    }
                }
            }
            
            TBOX_ASSERT(static_cast<int>(derivatives_y[ei].size()) ==
                        static_cast<int>(diffusivities_data_y[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_y[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_y[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_y[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_y[ei][vi];
                
                // Get the pointer to diffusivity.
                Real* mu = diffusivities_data_y[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivative.
                Real* dudy = derivatives_y[ei][vi]->getPointer(0);
                
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
                    j < interior_dim_1 + num_subghosts_1_diffusivity - num_flux_reconstruct_ghosts_1;
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
                        
                        F_node_y[ei][idx] += mu[idx_diffusivity]*dudy[idx];
                    }
                }
            }
        }
        
        /*
         * Reconstruct the flux in y-direction.
         */
        
        reconstructFluxY(
            diffusive_flux,
            diffusive_flux_node[1],
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
        
        derivatives_x.clear();
        derivatives_y.clear();
        
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
        
        diffusive_flux_utilities->registerDerivedVariablesForDiffusiveFluxes(d_num_diff_ghosts);
        
        diffusive_flux_utilities->allocateMemoryForDerivedCellData();
        
        diffusive_flux_utilities->computeDerivedCellData();
        
        /*
         * Delcare containers for computing fluxes in different directions.
         */
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > > var_data_x;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > > var_data_y;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > > var_data_z;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > > diffusivities_data_x;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > > diffusivities_data_y;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > > diffusivities_data_z;
        
        std::vector<std::vector<int> > var_component_idx_x;
        std::vector<std::vector<int> > var_component_idx_y;
        std::vector<std::vector<int> > var_component_idx_z;
        
        std::vector<std::vector<int> > diffusivities_component_idx_x;
        std::vector<std::vector<int> > diffusivities_component_idx_y;
        std::vector<std::vector<int> > diffusivities_component_idx_z;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > > derivatives_x;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > > derivatives_y;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > > derivatives_z;
        
        std::map<Real*, HAMERS_SHARED_PTR<pdat::CellData<Real> > > derivatives_x_computed;
        std::map<Real*, HAMERS_SHARED_PTR<pdat::CellData<Real> > > derivatives_y_computed;
        std::map<Real*, HAMERS_SHARED_PTR<pdat::CellData<Real> > > derivatives_z_computed;
        
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > diffusive_flux_node(3);
        diffusive_flux_node[0].reset(new pdat::CellData<Real>(interior_box, d_num_eqn, d_num_diff_ghosts));
        diffusive_flux_node[1].reset(new pdat::CellData<Real>(interior_box, d_num_eqn, d_num_diff_ghosts));
        diffusive_flux_node[2].reset(new pdat::CellData<Real>(interior_box, d_num_eqn, d_num_diff_ghosts));
        
        std::vector<Real*> F_node_x;
        std::vector<Real*> F_node_y;
        std::vector<Real*> F_node_z;
        F_node_x.reserve(d_num_eqn);
        F_node_y.reserve(d_num_eqn);
        F_node_z.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            F_node_x.push_back(diffusive_flux_node[0]->getPointer(ei));
            F_node_y.push_back(diffusive_flux_node[1]->getPointer(ei));
            F_node_z.push_back(diffusive_flux_node[2]->getPointer(ei));
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
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxDiffusivities(
            diffusivities_data_x,
            diffusivities_component_idx_x,
            DIRECTION::X_DIRECTION,
            DIRECTION::X_DIRECTION);
        
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxDiffusivities(
            diffusivities_data_y,
            diffusivities_component_idx_y,
            DIRECTION::X_DIRECTION,
            DIRECTION::Y_DIRECTION);
        
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxDiffusivities(
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
        
        computeFirstDerivativesInX(
            derivatives_x,
            derivatives_x_computed,
            var_data_x,
            var_component_idx_x,
            patch,
            allocate_scratch_data_containers);
        
        /*
         * Compute the derivatives in y-direction for diffusive flux in x-direction.
         */
        
        computeFirstDerivativesInY(
            derivatives_y,
            derivatives_y_computed,
            var_data_y,
            var_component_idx_y,
            patch,
            allocate_scratch_data_containers);
        
        /*
         * Compute the derivatives in z-direction for diffusive flux in x-direction.
         */
        
        computeFirstDerivativesInZ(
            derivatives_z,
            derivatives_z_computed,
            var_data_z,
            var_component_idx_z,
            patch,
            allocate_scratch_data_containers);
        
        /*
         * Compute diffusive flux in x-direction at nodes.
         */
        
        diffusive_flux_node[0]->fillAll(Real(0));
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ASSERT(static_cast<int>(derivatives_x[ei].size()) ==
                        static_cast<int>(diffusivities_data_x[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_x[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_x[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_x[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_x[ei][vi];
                
                // Get the pointer to diffusivity.
                Real* mu = diffusivities_data_x[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivative.
                Real* dudx = derivatives_x[ei][vi]->getPointer(0);
                
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
                        for (int i = -3; i < interior_dim_0 + 3; i++)
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
                            
                            F_node_x[ei][idx] += mu[idx_diffusivity]*dudx[idx];
                        }
                    }
                }
            }
            
            TBOX_ASSERT(static_cast<int>(derivatives_y[ei].size()) ==
                        static_cast<int>(diffusivities_data_y[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_y[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_y[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_y[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_y[ei][vi];
                
                // Get the pointer to diffusivity.
                Real* mu = diffusivities_data_y[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivative.
                Real* dudy = derivatives_y[ei][vi]->getPointer(0);
                
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
                        for (int i = -3; i < interior_dim_0 + 3; i++)
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
                            
                            F_node_x[ei][idx] += mu[idx_diffusivity]*dudy[idx];
                        }
                    }
                }
            }
            
            TBOX_ASSERT(static_cast<int>(derivatives_z[ei].size()) ==
                        static_cast<int>(diffusivities_data_z[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_z[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_z[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_z[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_z[ei][vi];
                
                // Get the pointer to diffusivity.
                Real* mu = diffusivities_data_z[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivative.
                Real* dudz = derivatives_z[ei][vi]->getPointer(0);
                
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
                            i < interior_dim_0 + num_subghosts_0_diffusivity - num_flux_reconstruct_ghosts_0;
                            i++)
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
                            
                            F_node_x[ei][idx] += mu[idx_diffusivity]*dudz[idx];
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
            diffusive_flux_node[0],
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
        
        derivatives_x.clear();
        derivatives_y.clear();
        derivatives_z.clear();
        
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
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxDiffusivities(
            diffusivities_data_x,
            diffusivities_component_idx_x,
            DIRECTION::Y_DIRECTION,
            DIRECTION::X_DIRECTION);
        
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxDiffusivities(
            diffusivities_data_y,
            diffusivities_component_idx_y,
            DIRECTION::Y_DIRECTION,
            DIRECTION::Y_DIRECTION);
        
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxDiffusivities(
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
        
        computeFirstDerivativesInX(
            derivatives_x,
            derivatives_x_computed,
            var_data_x,
            var_component_idx_x,
            patch,
            allocate_scratch_data_containers);
        
        /*
         * Compute the derivatives in y-direction for diffusive flux in y-direction.
         */
        
        computeFirstDerivativesInY(
            derivatives_y,
            derivatives_y_computed,
            var_data_y,
            var_component_idx_y,
            patch,
            allocate_scratch_data_containers);
        
        /*
         * Compute the derivatives in z-direction for diffusive flux in y-direction.
         */
        
        computeFirstDerivativesInZ(
            derivatives_z,
            derivatives_z_computed,
            var_data_z,
            var_component_idx_z,
            patch,
            allocate_scratch_data_containers);
        
        /*
         * Compute diffusive flux in y-direction at nodes.
         */
        
        diffusive_flux_node[1]->fillAll(Real(0));
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ASSERT(static_cast<int>(derivatives_x[ei].size()) ==
                        static_cast<int>(diffusivities_data_x[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_x[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_x[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_x[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_x[ei][vi];
                
                // Get the pointer to diffusivity.
                Real* mu = diffusivities_data_x[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivative.
                Real* dudx = derivatives_x[ei][vi]->getPointer(0);
                
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
                    for (int j = -3; j < interior_dim_1 + 3; j++)
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
                            
                            F_node_y[ei][idx] += mu[idx_diffusivity]*dudx[idx];
                        }
                    }
                }
            }
            
            TBOX_ASSERT(static_cast<int>(derivatives_y[ei].size()) ==
                        static_cast<int>(diffusivities_data_y[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_y[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_y[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_y[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_y[ei][vi];
                
                // Get the pointer to diffusivity.
                Real* mu = diffusivities_data_y[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivative.
                Real* dudy = derivatives_y[ei][vi]->getPointer(0);
                
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
                    for (int j = -3; j < interior_dim_1 + 3; j++)
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
                            
                            F_node_y[ei][idx] += mu[idx_diffusivity]*dudy[idx];
                        }
                    }
                }
            }
            
            TBOX_ASSERT(static_cast<int>(derivatives_z[ei].size()) ==
                        static_cast<int>(diffusivities_data_z[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_z[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_z[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_z[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_z[ei][vi];
                
                // Get the pointer to diffusivity.
                Real* mu = diffusivities_data_z[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivative.
                Real* dudz = derivatives_z[ei][vi]->getPointer(0);
                
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
                        j < interior_dim_1 + num_subghosts_1_diffusivity - num_flux_reconstruct_ghosts_1;
                        j++)
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
                            
                            F_node_y[ei][idx] += mu[idx_diffusivity]*dudz[idx];
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
            diffusive_flux_node[1],
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
        
        derivatives_x.clear();
        derivatives_y.clear();
        derivatives_z.clear();
        
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
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxDiffusivities(
            diffusivities_data_x,
            diffusivities_component_idx_x,
            DIRECTION::Z_DIRECTION,
            DIRECTION::X_DIRECTION);
        
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxDiffusivities(
            diffusivities_data_y,
            diffusivities_component_idx_y,
            DIRECTION::Z_DIRECTION,
            DIRECTION::Y_DIRECTION);
        
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxDiffusivities(
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
        
        computeFirstDerivativesInX(
            derivatives_x,
            derivatives_x_computed,
            var_data_x,
            var_component_idx_x,
            patch,
            allocate_scratch_data_containers);
        
        /*
         * Compute the derivatives in y-direction for diffusive flux in z-direction.
         */
        
        computeFirstDerivativesInY(
            derivatives_y,
            derivatives_y_computed,
            var_data_y,
            var_component_idx_y,
            patch,
            allocate_scratch_data_containers);
        
        /*
         * Compute the derivatives in z-direction for diffusive flux in z-direction.
         */
        
        computeFirstDerivativesInZ(
            derivatives_z,
            derivatives_z_computed,
            var_data_z,
            var_component_idx_z,
            patch,
            allocate_scratch_data_containers);
        
        /*
         * Compute diffusive flux in z-direction at nodes.
         */
        
        diffusive_flux_node[2]->fillAll(Real(0));
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ASSERT(static_cast<int>(derivatives_x[ei].size()) ==
                        static_cast<int>(diffusivities_data_x[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_x[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_x[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_x[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_x[ei][vi];
                
                // Get the pointer to diffusivity.
                Real* mu = diffusivities_data_x[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivative.
                Real* dudx = derivatives_x[ei][vi]->getPointer(0);
                
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
                
                for (int k = -3; k < interior_dim_2 + 3; k++)
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
                            
                            F_node_z[ei][idx] += mu[idx_diffusivity]*dudx[idx];
                        }
                    }
                }
            }
            
            TBOX_ASSERT(static_cast<int>(derivatives_y[ei].size()) ==
                        static_cast<int>(diffusivities_data_y[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_y[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_y[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_y[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_y[ei][vi];
                
                // Get the pointer to diffusivity.
                Real* mu = diffusivities_data_y[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivative.
                Real* dudy = derivatives_y[ei][vi]->getPointer(0);
                
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
                
                for (int k = -3; k < interior_dim_2 + 3; k++)
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
                            
                            F_node_z[ei][idx] += mu[idx_diffusivity]*dudy[idx];
                        }
                    }
                }
            }
            
            TBOX_ASSERT(static_cast<int>(derivatives_z[ei].size()) ==
                        static_cast<int>(diffusivities_data_z[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_z[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_z[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_z[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_z[ei][vi];
                
                // Get the pointer to diffusivity.
                Real* mu = diffusivities_data_z[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivative.
                Real* dudz = derivatives_z[ei][vi]->getPointer(0);
                
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
                    k < interior_dim_2 + num_subghosts_2_diffusivity - num_flux_reconstruct_ghosts_2;
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
                            
                            F_node_z[ei][idx] += mu[idx_diffusivity]*dudz[idx];
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
            diffusive_flux_node[2],
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
        
        derivatives_x.clear();
        derivatives_y.clear();
        derivatives_z.clear();
        
        /*
         * Unregister the patch and data of all registered derived cell variables in the flow model.
         */
        
        d_flow_model->unregisterPatch();
        
    } // if (d_dim == tbox::Dimension(3))
    
    // Clear the scratch data.
    d_scratch_derivatives_node.clear();
    
    d_num_scratch_derivatives_node_used = 0;
}


/*
 * Compute the first derivatives in the x-direction.
 */
void
DiffusiveFluxReconstructorNode::computeFirstDerivativesInX(
    std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& derivatives_x,
    std::map<Real*, HAMERS_SHARED_PTR<pdat::CellData<Real> > >& derivatives_x_computed,
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
    
    derivatives_x.reserve(static_cast<int>(data_x.size()));
    
    for (int vi = 0; vi < static_cast<int>(data_x.size()); vi++)
    {
        // Get the index of variable for derivative.
        const int u_idx = data_component_idx_x[vi];
        
        if (derivatives_x_computed.find(data_x[vi]->getPointer(u_idx))
            == derivatives_x_computed.end())
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
            
            computeFirstDerivativesInX(
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
            
            derivatives_x_computed.insert(derivative_pair);
        }
        
        derivatives_x.push_back(
            derivatives_x_computed.find(data_x[vi]->getPointer(u_idx))->
                second);
    }
}


/*
 * Compute the first derivatives in the x-direction.
 */
void
DiffusiveFluxReconstructorNode::computeFirstDerivativesInX(
    std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > >& derivatives_x,
    std::map<Real*, HAMERS_SHARED_PTR<pdat::CellData<Real> > >& derivatives_x_computed,
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
    
    derivatives_x.resize(num_eqn);
    
    for (int ei = 0; ei < num_eqn; ei++)
    {
        computeFirstDerivativesInX(
            derivatives_x[ei],
            derivatives_x_computed,
            data_x[ei],
            data_component_idx_x[ei],
            patch,
            allocate_scratch_derivatives_node);
    }
}


/*
 * Compute the first derivatives in the y-direction.
 */
void
DiffusiveFluxReconstructorNode::computeFirstDerivativesInY(
    std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& derivatives_y,
    std::map<Real*, HAMERS_SHARED_PTR<pdat::CellData<Real> > >& derivatives_y_computed,
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& data_y,
    const std::vector<int>& data_component_idx_y,
    const hier::Patch& patch,
    const bool allocate_scratch_derivatives_node)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": DiffusiveFluxReconstructorNode::"
            << "computeFirstDerivativesInY()\n"
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
    
    derivatives_y.reserve(static_cast<int>(data_y.size()));
    
    for (int vi = 0; vi < static_cast<int>(data_y.size()); vi++)
    {
        // Get the index of variable for derivative.
        const int u_idx = data_component_idx_y[vi];
        
        if (derivatives_y_computed.find(data_y[vi]->getPointer(u_idx))
            == derivatives_y_computed.end())
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
            
            computeFirstDerivativesInY(
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
            
            derivatives_y_computed.insert(derivative_pair);
        }
        
        derivatives_y.push_back(
            derivatives_y_computed.find(data_y[vi]->getPointer(u_idx))->
                second);
    }
}


/*
 * Compute the first derivatives in the y-direction.
 */
void
DiffusiveFluxReconstructorNode::computeFirstDerivativesInY(
    std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > >& derivatives_y,
    std::map<Real*, HAMERS_SHARED_PTR<pdat::CellData<Real> > >& derivatives_y_computed,
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
            << ": DiffusiveFluxReconstructorNode::"
            << "computeFirstDerivativesInY()\n"
            << "There isn't y-direction for one-dimensional problem."
            << std::endl);
    }
#endif
    
    derivatives_y.resize(num_eqn);
    
    for (int ei = 0; ei < num_eqn; ei++)
    {
        computeFirstDerivativesInY(
            derivatives_y[ei],
            derivatives_y_computed,
            data_y[ei],
            data_component_idx_y[ei],
            patch,
            allocate_scratch_derivatives_node);
    }
}


/*
 * Compute the first derivatives in the z-direction.
 */
void
DiffusiveFluxReconstructorNode::computeFirstDerivativesInZ(
    std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& derivatives_z,
    std::map<Real*, HAMERS_SHARED_PTR<pdat::CellData<Real> > >& derivatives_z_computed,
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& data_z,
    const std::vector<int>& data_component_idx_z,
    const hier::Patch& patch,
    const bool allocate_scratch_derivatives_node)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": DiffusiveFluxReconstructorNode::"
            << "computeFirstDerivativesInZ()\n"
            << "There isn't z-direction for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": DiffusiveFluxReconstructorNode::"
            << "computeFirstDerivativesInZ()\n"
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
    
    derivatives_z.reserve(static_cast<int>(data_z.size()));
    
    for (int vi = 0; vi < static_cast<int>(data_z.size()); vi++)
    {
        // Get the index of variable for derivative.
        const int u_idx = data_component_idx_z[vi];
        
        if (derivatives_z_computed.find(data_z[vi]->getPointer(u_idx))
            == derivatives_z_computed.end())
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
            
            computeFirstDerivativesInZ(
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
            
            derivatives_z_computed.insert(derivative_pair);
        }
        
        derivatives_z.push_back(
            derivatives_z_computed.find(data_z[vi]->getPointer(u_idx))->
                second);
    }
}


/*
 * Compute the first derivatives in the z-direction.
 */
void
DiffusiveFluxReconstructorNode::computeFirstDerivativesInZ(
    std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > >& derivatives_z,
    std::map<Real*, HAMERS_SHARED_PTR<pdat::CellData<Real> > >& derivatives_z_computed,
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
            << ": DiffusiveFluxReconstructorNode::"
            << "computeFirstDerivativesInZ()\n"
            << "There isn't z-direction for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": DiffusiveFluxReconstructorNode::"
            << "computeFirstDerivativesInZ()\n"
            << "There isn't z-direction for two-dimensional problem."
            << std::endl);
    }
#endif
    
    derivatives_z.resize(num_eqn);
    
    for (int ei = 0; ei < num_eqn; ei++)
    {
        computeFirstDerivativesInZ(
            derivatives_z[ei],
            derivatives_z_computed,
            data_z[ei],
            data_component_idx_z[ei],
            patch,
            allocate_scratch_derivatives_node);
    }
}


/*
 * Reconstruct the flux using flux at nodes in x-direction.
 */
void
DiffusiveFluxReconstructorNode::reconstructFluxX(
    HAMERS_SHARED_PTR<pdat::SideData<Real> >& diffusive_flux,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& diffusive_flux_node,
    const hier::Patch& patch,
    const double dt) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(diffusive_flux->getDepth() == d_num_eqn);
    TBOX_ASSERT(diffusive_flux_node->getDepth() == d_num_eqn);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // diffusive ghost cells.
    hier::Box diff_ghost_box = interior_box;
    diff_ghost_box.grow(d_num_diff_ghosts);
    const hier::IntVector ghostcell_dims_flux_node = diff_ghost_box.numberCells();
    
    const hier::IntVector& num_ghosts_flux_node = d_num_diff_ghosts;
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    domain_lo = hier::IntVector::getZero(d_dim);
    domain_dims = interior_dims;
    domain_dims[0]++;
    
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        Real* F_face_x = diffusive_flux->getPointer(0, ei);
        Real* F_node_x = diffusive_flux_node->getPointer(ei);
        
        reconstructFluxX(
            F_face_x,
            F_node_x,
            num_ghosts_flux_node,
            ghostcell_dims_flux_node,
            domain_lo,
            domain_dims,
            interior_dims,
            dt);
    }
}


/*
 * Reconstruct the flux using flux at nodes in y-direction.
 */
void
DiffusiveFluxReconstructorNode::reconstructFluxY(
    HAMERS_SHARED_PTR<pdat::SideData<Real> >& diffusive_flux,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& diffusive_flux_node,
    const hier::Patch& patch,
    const double dt) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(diffusive_flux->getDepth() == d_num_eqn);
    TBOX_ASSERT(diffusive_flux_node->getDepth() == d_num_eqn);
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": DiffusiveFluxReconstructorNode::"
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
    const hier::IntVector ghostcell_dims_flux_node = diff_ghost_box.numberCells();
    
    const hier::IntVector& num_ghosts_flux_node = d_num_diff_ghosts;
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    domain_lo = hier::IntVector::getZero(d_dim);
    domain_dims = interior_dims;
    domain_dims[1]++;
    
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        Real* F_face_y = diffusive_flux->getPointer(1, ei);
        Real* F_node_y = diffusive_flux_node->getPointer(ei);
        
        reconstructFluxY(
            F_face_y,
            F_node_y,
            num_ghosts_flux_node,
            ghostcell_dims_flux_node,
            domain_lo,
            domain_dims,
            interior_dims,
            dt);
    }
}


/*
 * Reconstruct the flux using flux at nodes in z-direction.
 */
void
DiffusiveFluxReconstructorNode::reconstructFluxZ(
    HAMERS_SHARED_PTR<pdat::SideData<Real> >& diffusive_flux,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& diffusive_flux_node,
    const hier::Patch& patch,
    const double dt) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(diffusive_flux->getDepth() == d_num_eqn);
    TBOX_ASSERT(diffusive_flux_node->getDepth() == d_num_eqn);
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": DiffusiveFluxReconstructorNode::"
            << "reconstructFluxZ()\n"
            << "There isn't z-direction for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": DiffusiveFluxReconstructorNode::"
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
    const hier::IntVector ghostcell_dims_flux_node = diff_ghost_box.numberCells();
    
    const hier::IntVector& num_ghosts_flux_node = d_num_diff_ghosts;
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    domain_lo = hier::IntVector::getZero(d_dim);
    domain_dims = interior_dims;
    domain_dims[2]++;
    
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        Real* F_face_z = diffusive_flux->getPointer(2, ei);
        Real* F_node_z = diffusive_flux_node->getPointer(ei);
        
        reconstructFluxZ(
            F_face_z,
            F_node_z,
            num_ghosts_flux_node,
            ghostcell_dims_flux_node,
            domain_lo,
            domain_dims,
            interior_dims,
            dt);
    }
}
