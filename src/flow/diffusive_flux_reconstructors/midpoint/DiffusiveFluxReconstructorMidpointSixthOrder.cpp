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
    d_num_diff_ghosts = hier::IntVector::getOne(d_dim)*5;
    
    // The number of ghost cells used for finite difference and interpolation schemes at midpoints are the same.
    d_num_der_midpoint_ghosts     = hier::IntVector::getOne(d_dim)*3;
    d_num_der_node_ghosts         = hier::IntVector::getOne(d_dim)*3;
    d_num_interp_midpoint_ghosts  = hier::IntVector::getOne(d_dim)*3;
    d_num_flux_reconstruct_ghosts = hier::IntVector::getOne(d_dim)*3;
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
        
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > var_cell_data_for_diffusivities;
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > var_side_data_for_diffusivities;
        std::vector<int> var_cell_data_for_diffusivities_component_idx;
        
        diffusive_flux_utilities->getCellDataForInterpolationToSideDataForDiffusiveFluxDiffusivities(
            var_cell_data_for_diffusivities,
            var_cell_data_for_diffusivities_component_idx);
        
        // Interpolate variables from nodes to midpoints for computing diffusivities at midpoints.
        
        interpolateDiffusivitiesFromNodeToMidpoint(
            patch,
            var_side_data_for_diffusivities,
            var_cell_data_for_diffusivities,
            var_cell_data_for_diffusivities_component_idx);
        
        diffusive_flux_utilities->allocateMemoryForSideDataOfDiffusiveFluxDiffusivities();
        
        diffusive_flux_utilities->computeSideDataOfDiffusiveFluxDiffusivities(
            var_side_data_for_diffusivities);
        
        /*
         * Delcare containers for computing fluxes in different directions.
         */
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > var_data_x;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > > diffusivities_data_x;
        
        std::vector<std::vector<int> > var_component_idx_x;
        
        std::vector<std::vector<int> > diffusivities_component_idx_x;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > > derivative_x_midpoint;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > derivative_x_node;
        
        std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivative_x_node_computed;
        
        HAMERS_SHARED_PTR<pdat::SideData<double> > diffusive_flux_midpoint;
        diffusive_flux_midpoint.reset(new pdat::SideData<double>(interior_box, d_num_eqn, d_num_diff_ghosts));
        
        std::vector<double*> F_midpoint_x;
        F_midpoint_x.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            F_midpoint_x.push_back(diffusive_flux_midpoint->getPointer(0, ei));
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
        
        computeFirstDerivativesInXAtMidpoint(
            patch,
            derivative_x_midpoint,
            var_data_x,
            var_component_idx_x);
        
        /*
         * Compute diffusive flux in x-direction at midpoints.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ASSERT(static_cast<int>(derivative_x_midpoint[ei].size()) ==
                        static_cast<int>(diffusivities_data_x[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_x[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_x[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_x[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_x[ei][vi];
                
                // Get the pointer to diffusivity.
                double* mu = diffusivities_data_x[ei][vi]->getPointer(0, mu_idx);
                
                // Get the pointer to derivative.
                double* dudx = derivative_x_midpoint[ei][vi]->getPointer(0, 0);
                
                /*
                 * Get the sub-ghost cell width of the diffusivity.
                 */
                
                hier::IntVector num_subghosts_diffusivity =
                    diffusivities_data_x[ei][vi]->getGhostCellWidth();
                
                const int num_subghosts_0_diffusivity = num_subghosts_diffusivity[0];
                
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
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
            patch,
            diffusive_flux,
            diffusive_flux_midpoint,
            dt);
        
        var_data_x.clear();
        diffusivities_data_x.clear();
        var_component_idx_x.clear();
        diffusivities_component_idx_x.clear();
        derivative_x_midpoint.clear();
        
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
        
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > var_cell_data_for_diffusivities;
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > var_side_data_for_diffusivities;
        std::vector<int> var_cell_data_for_diffusivities_component_idx;
        
        diffusive_flux_utilities->getCellDataForInterpolationToSideDataForDiffusiveFluxDiffusivities(
            var_cell_data_for_diffusivities,
            var_cell_data_for_diffusivities_component_idx);
        
        // Interpolate variables from nodes to midpoints for computing diffusivities at midpoints.
        
        interpolateDiffusivitiesFromNodeToMidpoint(
            patch,
            var_side_data_for_diffusivities,
            var_cell_data_for_diffusivities,
            var_cell_data_for_diffusivities_component_idx);
        
        diffusive_flux_utilities->allocateMemoryForSideDataOfDiffusiveFluxDiffusivities();
        
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
        
        interpolateDerivativesFromNodeToMidpointX(
            patch,
            derivative_y_midpoint,
            derivative_y_node);
        
        /*
         * Compute diffusive flux in x-direction at midpoints.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ASSERT(static_cast<int>(derivative_x_midpoint[ei].size()) ==
                        static_cast<int>(diffusivities_data_x[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_x[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_x[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_x[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_x[ei][vi];
                
                // Get the pointer to diffusivity.
                double* mu = diffusivities_data_x[ei][vi]->getPointer(0, mu_idx);
                
                // Get the pointer to derivative.
                double* dudx = derivative_x_midpoint[ei][vi]->getPointer(0, 0);
                
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
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
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
            
            TBOX_ASSERT(static_cast<int>(derivative_y_midpoint[ei].size()) ==
                        static_cast<int>(diffusivities_data_y[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_y[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_y[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_y[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_y[ei][vi];
                
                // Get the pointer to diffusivity.
                double* mu = diffusivities_data_y[ei][vi]->getPointer(0, mu_idx);
                
                // Get the pointer to derivative.
                double* dudy = derivative_y_midpoint[ei][vi]->getPointer(0, 0);
                
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
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
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
            patch,
            diffusive_flux,
            diffusive_flux_midpoint,
            dt);
        
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
        
        interpolateDerivativesFromNodeToMidpointY(
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
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ASSERT(static_cast<int>(derivative_x_midpoint[ei].size()) ==
                        static_cast<int>(diffusivities_data_x[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_x[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_x[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_x[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_x[ei][vi];
                
                // Get the pointer to diffusivity.
                double* mu = diffusivities_data_x[ei][vi]->getPointer(1, mu_idx);
                
                // Get the pointer to derivative.
                double* dudx = derivative_x_midpoint[ei][vi]->getPointer(1, 0);
                
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
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
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
            
            TBOX_ASSERT(static_cast<int>(derivative_y_midpoint[ei].size()) ==
                        static_cast<int>(diffusivities_data_y[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_y[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_y[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_y[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_y[ei][vi];
                
                // Get the pointer to diffusivity.
                double* mu = diffusivities_data_y[ei][vi]->getPointer(1, mu_idx);
                
                // Get the pointer to derivative.
                double* dudy = derivative_y_midpoint[ei][vi]->getPointer(1, 0);
                
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
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
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
            patch,
            diffusive_flux,
            diffusive_flux_midpoint,
            dt);
        
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
        
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > var_cell_data_for_diffusivities;
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > var_side_data_for_diffusivities;
        std::vector<int> var_cell_data_for_diffusivities_component_idx;
        
        diffusive_flux_utilities->getCellDataForInterpolationToSideDataForDiffusiveFluxDiffusivities(
            var_cell_data_for_diffusivities,
            var_cell_data_for_diffusivities_component_idx);
        
        // Interpolate variables from nodes to midpoints for computing diffusivities at midpoints.
        
        interpolateDiffusivitiesFromNodeToMidpoint(
            patch,
            var_side_data_for_diffusivities,
            var_cell_data_for_diffusivities,
            var_cell_data_for_diffusivities_component_idx);
        
        diffusive_flux_utilities->allocateMemoryForSideDataOfDiffusiveFluxDiffusivities();
        
        diffusive_flux_utilities->computeSideDataOfDiffusiveFluxDiffusivities(
            var_side_data_for_diffusivities);
        
        /*
         * Delcare containers for computing fluxes in different directions.
         */
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > var_data_x;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > var_data_y;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > var_data_z;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > > diffusivities_data_x;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > > diffusivities_data_y;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > > diffusivities_data_z;
        
        std::vector<std::vector<int> > var_component_idx_x;
        std::vector<std::vector<int> > var_component_idx_y;
        std::vector<std::vector<int> > var_component_idx_z;
        
        std::vector<std::vector<int> > diffusivities_component_idx_x;
        std::vector<std::vector<int> > diffusivities_component_idx_y;
        std::vector<std::vector<int> > diffusivities_component_idx_z;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > > derivative_x_midpoint;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > > derivative_y_midpoint;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > > derivative_z_midpoint;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > derivative_x_node;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > derivative_y_node;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > derivative_z_node;
        
        std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivative_x_node_computed;
        std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivative_y_node_computed;
        std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivative_z_node_computed;
        
        HAMERS_SHARED_PTR<pdat::SideData<double> > diffusive_flux_midpoint;
        diffusive_flux_midpoint.reset(new pdat::SideData<double>(interior_box, d_num_eqn, d_num_diff_ghosts));
        
        std::vector<double*> F_midpoint_x;
        std::vector<double*> F_midpoint_y;
        std::vector<double*> F_midpoint_z;
        F_midpoint_x.reserve(d_num_eqn);
        F_midpoint_y.reserve(d_num_eqn);
        F_midpoint_z.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            F_midpoint_x.push_back(diffusive_flux_midpoint->getPointer(0, ei));
            F_midpoint_y.push_back(diffusive_flux_midpoint->getPointer(1, ei));
            F_midpoint_z.push_back(diffusive_flux_midpoint->getPointer(2, ei));
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
        
        interpolateDerivativesFromNodeToMidpointX(
            patch,
            derivative_y_midpoint,
            derivative_y_node);
        
        /*
         * Compute the derivatives in z-direction for diffusive flux in x-direction.
         */
        
        computeFirstDerivativesInZAtNode(
            patch,
            derivative_z_node,
            derivative_z_node_computed,
            var_data_z,
            var_component_idx_z);
        
        interpolateDerivativesFromNodeToMidpointX(
            patch,
            derivative_z_midpoint,
            derivative_z_node);
        
        /*
         * Compute diffusive flux in x-direction at midpoints.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ASSERT(static_cast<int>(derivative_x_midpoint[ei].size()) ==
                        static_cast<int>(diffusivities_data_x[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_x[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_x[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_x[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_x[ei][vi];
                
                // Get the pointer to diffusivity.
                double* mu = diffusivities_data_x[ei][vi]->getPointer(0, mu_idx);
                
                // Get the pointer to derivative.
                double* dudx = derivative_x_midpoint[ei][vi]->getPointer(0, 0);
                
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
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
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
            
            TBOX_ASSERT(static_cast<int>(derivative_y_midpoint[ei].size()) ==
                        static_cast<int>(diffusivities_data_y[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_y[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_y[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_y[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_y[ei][vi];
                
                // Get the pointer to diffusivity.
                double* mu = diffusivities_data_y[ei][vi]->getPointer(0, mu_idx);
                
                // Get the pointer to derivative.
                double* dudy = derivative_y_midpoint[ei][vi]->getPointer(0, 0);
                
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
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
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
            
            TBOX_ASSERT(static_cast<int>(derivative_z_midpoint[ei].size()) ==
                        static_cast<int>(diffusivities_data_z[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_z[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_z[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_z[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_z[ei][vi];
                
                // Get the pointer to diffusivity.
                double* mu = diffusivities_data_z[ei][vi]->getPointer(0, mu_idx);
                
                // Get the pointer to derivative.
                double* dudz = derivative_z_midpoint[ei][vi]->getPointer(0, 0);
                
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
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
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
            patch,
            diffusive_flux,
            diffusive_flux_midpoint,
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
        
        derivative_x_midpoint.clear();
        derivative_y_midpoint.clear();
        derivative_z_midpoint.clear();
        
        derivative_y_node.clear();
        derivative_z_node.clear();
        
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
            patch,
            derivative_x_node,
            derivative_x_node_computed,
            var_data_x,
            var_component_idx_x);
        
        interpolateDerivativesFromNodeToMidpointY(
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
         * Compute the derivatives in z-direction for diffusive flux in y-direction.
         */
        
        computeFirstDerivativesInZAtNode(
            patch,
            derivative_z_node,
            derivative_z_node_computed,
            var_data_z,
            var_component_idx_z);
        
        interpolateDerivativesFromNodeToMidpointY(
            patch,
            derivative_z_midpoint,
            derivative_z_node);
        
        /*
         * Compute diffusive flux in y-direction at midpoints.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ASSERT(static_cast<int>(derivative_x_midpoint[ei].size()) ==
                        static_cast<int>(diffusivities_data_x[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_x[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_x[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_x[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_x[ei][vi];
                
                // Get the pointer to diffusivity.
                double* mu = diffusivities_data_x[ei][vi]->getPointer(1, mu_idx);
                
                // Get the pointer to derivative.
                double* dudx = derivative_x_midpoint[ei][vi]->getPointer(1, 0);
                
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
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
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
            
            TBOX_ASSERT(static_cast<int>(derivative_y_midpoint[ei].size()) ==
                        static_cast<int>(diffusivities_data_y[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_y[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_y[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_y[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_y[ei][vi];
                
                // Get the pointer to diffusivity.
                double* mu = diffusivities_data_y[ei][vi]->getPointer(1, mu_idx);
                
                // Get the pointer to derivative.
                double* dudy = derivative_y_midpoint[ei][vi]->getPointer(1, 0);
                
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
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
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
            
            TBOX_ASSERT(static_cast<int>(derivative_z_midpoint[ei].size()) ==
                        static_cast<int>(diffusivities_data_z[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_z[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_z[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_z[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_z[ei][vi];
                
                // Get the pointer to diffusivity.
                double* mu = diffusivities_data_z[ei][vi]->getPointer(1, mu_idx);
                
                // Get the pointer to derivative.
                double* dudz = derivative_z_midpoint[ei][vi]->getPointer(1, 0);
                
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
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
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
            patch,
            diffusive_flux,
            diffusive_flux_midpoint,
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
        
        derivative_x_midpoint.clear();
        derivative_y_midpoint.clear();
        derivative_z_midpoint.clear();
        
        derivative_x_node.clear();
        derivative_z_node.clear();
        
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
            patch,
            derivative_x_node,
            derivative_x_node_computed,
            var_data_x,
            var_component_idx_x);
        
        interpolateDerivativesFromNodeToMidpointZ(
            patch,
            derivative_x_midpoint,
            derivative_x_node);
        
        /*
         * Compute the derivatives in y-direction for diffusive flux in z-direction.
         */
        
        computeFirstDerivativesInYAtNode(
            patch,
            derivative_y_node,
            derivative_y_node_computed,
            var_data_y,
            var_component_idx_y);
        
        interpolateDerivativesFromNodeToMidpointZ(
            patch,
            derivative_y_midpoint,
            derivative_y_node);
        
        /*
         * Compute the derivatives in z-direction for diffusive flux in z-direction.
         */
        
        computeFirstDerivativesInZAtMidpoint(
            patch,
            derivative_z_midpoint,
            var_data_z,
            var_component_idx_z);
        
        
        /*
         * Compute diffusive flux in z-direction at midpoints.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ASSERT(static_cast<int>(derivative_x_midpoint[ei].size()) ==
                        static_cast<int>(diffusivities_data_x[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_x[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_x[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_x[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_x[ei][vi];
                
                // Get the pointer to diffusivity.
                double* mu = diffusivities_data_x[ei][vi]->getPointer(2, mu_idx);
                
                // Get the pointer to derivative.
                double* dudx = derivative_x_midpoint[ei][vi]->getPointer(2, 0);
                
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
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
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
            
            TBOX_ASSERT(static_cast<int>(derivative_y_midpoint[ei].size()) ==
                        static_cast<int>(diffusivities_data_y[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_y[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_y[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_y[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_y[ei][vi];
                
                // Get the pointer to diffusivity.
                double* mu = diffusivities_data_y[ei][vi]->getPointer(2, mu_idx);
                
                // Get the pointer to derivative.
                double* dudy = derivative_y_midpoint[ei][vi]->getPointer(2, 0);
                
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
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
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
            
            TBOX_ASSERT(static_cast<int>(derivative_z_midpoint[ei].size()) ==
                        static_cast<int>(diffusivities_data_z[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_z[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_z[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_z[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_z[ei][vi];
                
                // Get the pointer to diffusivity.
                double* mu = diffusivities_data_z[ei][vi]->getPointer(2, mu_idx);
                
                // Get the pointer to derivative.
                double* dudz = derivative_z_midpoint[ei][vi]->getPointer(2, 0);
                
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
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
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
            patch,
            diffusive_flux,
            diffusive_flux_midpoint,
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
        
        derivative_x_midpoint.clear();
        derivative_y_midpoint.clear();
        derivative_z_midpoint.clear();
        
        derivative_x_node.clear();
        derivative_y_node.clear();
        
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
DiffusiveFluxReconstructorMidpointSixthOrder::computeFirstDerivativesInXAtMidpoint(
    hier::Patch& patch,
    std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > >& derivative_x_midpoint,
    const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& data_x,
    const std::vector<std::vector<int> >& data_component_idx_x,
    const bool allocate_scratch_derivatives_midpoint)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        TBOX_ASSERT(static_cast<int>(data_x[ei].size()) ==
                    static_cast<int>(data_component_idx_x[ei].size()));
    }
#endif
    
    const double a_n =  double(75)/double(64);
    const double b_n = -double(25)/double(384);
    const double c_n =  double(3)/double(640);
    
    derivative_x_midpoint.resize(d_num_eqn);
    
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
    
    const double inv_dx_0 = double(1)/dx[0];
    
    hier::IntVector direction = hier::IntVector::getZero(d_dim);
    direction[0] = 1;
    
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        derivative_x_midpoint[ei].reserve(static_cast<int>(data_x[ei].size()));
        
        for (int vi = 0; vi < static_cast<int>(data_x[ei].size()); vi++)
        {
            // Get the index of variable for derivative.
            const int u_idx = data_component_idx_x[ei][vi];
            
            // Get the pointer to variable for derivative.
            double* u = data_x[ei][vi]->getPointer(u_idx);
            
            d_num_scratch_derivatives_midpoint_x_used++;
            // Declare container to store the derivative.
            if (allocate_scratch_derivatives_midpoint)
            {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
                TBOX_ASSERT(static_cast<int>(d_scratch_derivatives_midpoint_x.size()) == d_num_scratch_derivatives_midpoint_x_used - 1);
#endif
                d_scratch_derivatives_midpoint_x.push_back(HAMERS_SHARED_PTR<pdat::SideData<double> >(
                    new pdat::SideData<double>(
                        interior_box, 1, d_num_diff_ghosts, direction)));
            }
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
            else
            {
                TBOX_ASSERT(static_cast<int>(d_scratch_derivatives_midpoint_x.size()) >= d_num_scratch_derivatives_midpoint_x_used);
            }
#endif
            
            // Get the pointer to the derivative at midpoint.
            double* der_midpoint_x = d_scratch_derivatives_midpoint_x[d_num_scratch_derivatives_midpoint_x_used - 1]->getPointer(0, 0);
            
            /*
             * Get the sub-ghost cell width of the variable.
             */
            
            hier::IntVector num_subghosts_data =
                data_x[ei][vi]->getGhostCellWidth();
            
            hier::IntVector subghostcell_dims_data =
                data_x[ei][vi]->getGhostBox().numberCells();
            
            if (d_dim == tbox::Dimension(1))
            {
                /*
                 * Get the dimensions and number of ghost cells.
                 */
                
                const int interior_dim_0 = interior_dims[0];
                
                const int num_diff_ghosts_0 = d_num_diff_ghosts[0];
                
                const int num_subghosts_0_data = num_subghosts_data[0];
                
                const int num_der_midpoint_ghosts_0 = d_num_der_midpoint_ghosts[0];
                
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = -num_diff_ghosts_0 + num_der_midpoint_ghosts_0;
                    i < (interior_dim_0 + 1) + num_diff_ghosts_0 - num_der_midpoint_ghosts_0;
                    i++)
                {
                    // Compute the linear indices.
                    const int idx = i + num_diff_ghosts_0;
                    
                    const int idx_data_LLL = i - 3 + num_subghosts_0_data;
                    const int idx_data_LL  = i - 2 + num_subghosts_0_data;
                    const int idx_data_L   = i - 1 + num_subghosts_0_data;
                    const int idx_data_R   = i + 0 + num_subghosts_0_data;
                    const int idx_data_RR  = i + 1 + num_subghosts_0_data;
                    const int idx_data_RRR = i + 2 + num_subghosts_0_data;
                    
                    der_midpoint_x[idx] = (a_n*(u[idx_data_R]   - u[idx_data_L]) +
                                           b_n*(u[idx_data_RR]  - u[idx_data_LL]) +
                                           c_n*(u[idx_data_RRR] - u[idx_data_LLL])
                                          )*inv_dx_0;
                }
            }
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
                
                const int num_subghosts_0_data = num_subghosts_data[0];
                const int num_subghosts_1_data = num_subghosts_data[1];
                const int subghostcell_dim_0_data = subghostcell_dims_data[0];
                
                const int num_der_midpoint_ghosts_0 = d_num_der_midpoint_ghosts[0];
                
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = -num_diff_ghosts_0 + num_der_midpoint_ghosts_0;
                        i < (interior_dim_0 + 1) + num_diff_ghosts_0 - num_der_midpoint_ghosts_0;
                        i++)
                    {
                        // Compute the linear indices.
                        const int idx = (i + num_diff_ghosts_0) +
                            (j + num_diff_ghosts_1)*(diff_ghostcell_dim_0 + 1);
                        
                        const int idx_data_LLL = (i - 3 + num_subghosts_0_data) +
                            (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                        
                        const int idx_data_LL = (i - 2 + num_subghosts_0_data) +
                            (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                        
                        const int idx_data_L = (i - 1 + num_subghosts_0_data) +
                            (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                        
                        const int idx_data_R = (i + num_subghosts_0_data) +
                            (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                        
                        const int idx_data_RR = (i + 1 + num_subghosts_0_data) +
                            (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                        
                        const int idx_data_RRR = (i + 2 + num_subghosts_0_data) +
                            (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                        
                        der_midpoint_x[idx] = (a_n*(u[idx_data_R]   - u[idx_data_L]) +
                                               b_n*(u[idx_data_RR]  - u[idx_data_LL]) +
                                               c_n*(u[idx_data_RRR] - u[idx_data_LLL])
                                              )*inv_dx_0;
                    }
                }
            }
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
                
                const int num_subghosts_0_data = num_subghosts_data[0];
                const int num_subghosts_1_data = num_subghosts_data[1];
                const int num_subghosts_2_data = num_subghosts_data[2];
                const int subghostcell_dim_0_data = subghostcell_dims_data[0];
                const int subghostcell_dim_1_data = subghostcell_dims_data[1];
                
                const int num_der_midpoint_ghosts_0 = d_num_der_midpoint_ghosts[0];
                
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = -num_diff_ghosts_0 + num_der_midpoint_ghosts_0;
                            i < (interior_dim_0 + 1) + num_diff_ghosts_0 - num_der_midpoint_ghosts_0;
                            i++)
                        {
                            // Compute the linear indices.
                            const int idx = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*(diff_ghostcell_dim_0 + 1) +
                                (k + num_diff_ghosts_2)*(diff_ghostcell_dim_0 + 1)*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_data_LLL = (i - 3 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                    subghostcell_dim_1_data;
                            
                            const int idx_data_LL = (i - 2 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                    subghostcell_dim_1_data;
                            
                            const int idx_data_L = (i - 1 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                    subghostcell_dim_1_data;
                            
                            const int idx_data_R = (i + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                    subghostcell_dim_1_data;
                            
                            const int idx_data_RR = (i + 1 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                    subghostcell_dim_1_data;
                            
                            const int idx_data_RRR = (i + 2 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                    subghostcell_dim_1_data;
                            
                            der_midpoint_x[idx] = (a_n*(u[idx_data_R]   - u[idx_data_L]) +
                                                   b_n*(u[idx_data_RR]  - u[idx_data_LL]) +
                                                   c_n*(u[idx_data_RRR] - u[idx_data_LLL])
                                                  )*inv_dx_0;
                        }
                    }
                }
            }
            
            derivative_x_midpoint[ei].push_back(d_scratch_derivatives_midpoint_x[d_num_scratch_derivatives_midpoint_x_used - 1]);
        }
    }
}


/*
 * Compute the derivatives in y-direction at midpoints.
 */
void
DiffusiveFluxReconstructorMidpointSixthOrder::computeFirstDerivativesInYAtMidpoint(
    hier::Patch& patch,
    std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > >& derivative_y_midpoint,
    const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& data_y,
    const std::vector<std::vector<int> >& data_component_idx_y,
    const bool allocate_scratch_derivatives_midpoint)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        TBOX_ASSERT(static_cast<int>(data_y[ei].size()) ==
                    static_cast<int>(data_component_idx_y[ei].size()));
    }
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": DiffusiveFluxReconstructorMidpointSixthOrder::"
            << "computeFirstDerivativesInYAtMidpoint()\n"
            << "There isn't y-direction for 1D problem."
            << std::endl);
    }
#endif
    
    const double a_n =  double(75)/double(64);
    const double b_n = -double(25)/double(384);
    const double c_n =  double(3)/double(640);
    
    derivative_y_midpoint.resize(d_num_eqn);
    
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
    
    const double inv_dx_1 = double(1)/dx[1];
    
    hier::IntVector direction = hier::IntVector::getZero(d_dim);
    direction[1] = 1;
    
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        derivative_y_midpoint[ei].reserve(static_cast<int>(data_y[ei].size()));
        
        for (int vi = 0; vi < static_cast<int>(data_y[ei].size()); vi++)
        {
            // Get the index of variable for derivative.
            const int u_idx = data_component_idx_y[ei][vi];
            
            // Get the pointer to variable for derivative.
            double* u = data_y[ei][vi]->getPointer(u_idx);
            
            d_num_scratch_derivatives_midpoint_y_used++;
            // Declare container to store the derivative.
            if (allocate_scratch_derivatives_midpoint)
            {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
                TBOX_ASSERT(static_cast<int>(d_scratch_derivatives_midpoint_y.size()) == d_num_scratch_derivatives_midpoint_y_used - 1);
#endif
                d_scratch_derivatives_midpoint_y.push_back(HAMERS_SHARED_PTR<pdat::SideData<double> >(
                    new pdat::SideData<double>(
                        interior_box, 1, d_num_diff_ghosts, direction)));
            }
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
            else
            {
                TBOX_ASSERT(static_cast<int>(d_scratch_derivatives_midpoint_y.size()) >= d_num_scratch_derivatives_midpoint_y_used);
            }
#endif
            
            // Get the pointer to the derivative at midpoint.
            double* der_midpoint_y = d_scratch_derivatives_midpoint_y[d_num_scratch_derivatives_midpoint_y_used - 1]->getPointer(1, 0);
            
            /*
             * Get the sub-ghost cell width and ghost box dimensions of the variable.
             */
            
            hier::IntVector num_subghosts_data =
                data_y[ei][vi]->getGhostCellWidth();
            
            hier::IntVector subghostcell_dims_data =
                data_y[ei][vi]->getGhostBox().numberCells();
            
            if (d_dim == tbox::Dimension(2))
            {
                /*
                 * Get the dimensions and number of ghost cells.
                 */
                
                const int interior_dim_0 = interior_dims[0];
                const int interior_dim_1 = interior_dims[1];
                
                const int num_diff_ghosts_0 = d_num_diff_ghosts[0];
                const int num_diff_ghosts_1 = d_num_diff_ghosts[1];
                const int diff_ghostcell_dim_0 = diff_ghostcell_dims[0];
                
                const int num_subghosts_0_data = num_subghosts_data[0];
                const int num_subghosts_1_data = num_subghosts_data[1];
                const int subghostcell_dim_0_data = subghostcell_dims_data[0];
                
                const int num_der_midpoint_ghosts_1 = d_num_der_midpoint_ghosts[1];
                
                for (int j = -num_diff_ghosts_1 + num_der_midpoint_ghosts_1;
                    j < (interior_dim_1 + 1) + num_diff_ghosts_1 - num_der_midpoint_ghosts_1;
                    j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx = (i + num_diff_ghosts_0) +
                            (j + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                        
                        const int idx_data_BBB = (i + num_subghosts_0_data) +
                            (j - 3 + num_subghosts_1_data)*subghostcell_dim_0_data;
                        
                        const int idx_data_BB = (i + num_subghosts_0_data) +
                            (j - 2 + num_subghosts_1_data)*subghostcell_dim_0_data;
                        
                        const int idx_data_B = (i + num_subghosts_0_data) +
                            (j - 1 + num_subghosts_1_data)*subghostcell_dim_0_data;
                        
                        const int idx_data_T = (i + num_subghosts_0_data) +
                            (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                        
                        const int idx_data_TT = (i + num_subghosts_0_data) +
                            (j + 1 + num_subghosts_1_data)*subghostcell_dim_0_data;
                        
                        const int idx_data_TTT = (i + num_subghosts_0_data) +
                            (j + 2 + num_subghosts_1_data)*subghostcell_dim_0_data;
                        
                        der_midpoint_y[idx] = (a_n*(u[idx_data_T]   - u[idx_data_B]) +
                                               b_n*(u[idx_data_TT]  - u[idx_data_BB]) +
                                               c_n*(u[idx_data_TTT] - u[idx_data_BBB])
                                              )*inv_dx_1;
                    }
                }
            }
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
                
                const int num_subghosts_0_data = num_subghosts_data[0];
                const int num_subghosts_1_data = num_subghosts_data[1];
                const int num_subghosts_2_data = num_subghosts_data[2];
                const int subghostcell_dim_0_data = subghostcell_dims_data[0];
                const int subghostcell_dim_1_data = subghostcell_dims_data[1];
                
                const int num_der_midpoint_ghosts_1 = d_num_der_midpoint_ghosts[1];
                
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = -num_diff_ghosts_1 + num_der_midpoint_ghosts_1;
                        j < (interior_dim_1 + 1) + num_diff_ghosts_1 - num_der_midpoint_ghosts_1;
                        j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    (diff_ghostcell_dim_1 + 1);
                            
                            const int idx_data_BBB = (i + num_subghosts_0_data) +
                                (j - 3 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                    subghostcell_dim_1_data;
                            
                            const int idx_data_BB = (i + num_subghosts_0_data) +
                                (j - 2 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                    subghostcell_dim_1_data;
                            
                            const int idx_data_B = (i + num_subghosts_0_data) +
                                (j - 1 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                    subghostcell_dim_1_data;
                            
                            const int idx_data_T = (i + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                    subghostcell_dim_1_data;
                            
                            const int idx_data_TT = (i + num_subghosts_0_data) +
                                (j + 1 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                    subghostcell_dim_1_data;
                            
                            const int idx_data_TTT = (i + num_subghosts_0_data) +
                                (j + 2 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                    subghostcell_dim_1_data;
                            
                            der_midpoint_y[idx] = (a_n*(u[idx_data_T]   - u[idx_data_B]) +
                                                   b_n*(u[idx_data_TT]  - u[idx_data_BB]) +
                                                   c_n*(u[idx_data_TTT] - u[idx_data_BBB])
                                                  )*inv_dx_1;
                        }
                    }
                }
            }
            
            derivative_y_midpoint[ei].push_back(d_scratch_derivatives_midpoint_y[d_num_scratch_derivatives_midpoint_y_used - 1]);
        }
    }
}


/*
 * Compute the derivatives in z-direction at midpoints.
 */
void
DiffusiveFluxReconstructorMidpointSixthOrder::computeFirstDerivativesInZAtMidpoint(
    hier::Patch& patch,
    std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > >& derivative_z_midpoint,
    const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& data_z,
    const std::vector<std::vector<int> >& data_component_idx_z,
    const bool allocate_scratch_derivatives_midpoint)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        TBOX_ASSERT(static_cast<int>(data_y[ei].size()) ==
                    static_cast<int>(data_component_idx_y[ei].size()));
    }
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": DiffusiveFluxReconstructorMidpointSixthOrder::"
            << "computeFirstDerivativesInZAtMidpoint()\n"
            << "There isn't z-direction for 1D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": DiffusiveFluxReconstructorMidpointSixthOrder::"
            << "computeFirstDerivativesInZAtMidpoint()\n"
            << "There isn't z-direction for 2D problem."
            << std::endl);
    }
#endif
    
    const double a_n =  double(75)/double(64);
    const double b_n = -double(25)/double(384);
    const double c_n =  double(3)/double(640);
    
    derivative_z_midpoint.resize(d_num_eqn);
    
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
    
    const double inv_dx_2 = double(1)/dx[2];
    
    hier::IntVector direction = hier::IntVector::getZero(d_dim);
    direction[2] = 1;
    
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        derivative_z_midpoint[ei].reserve(static_cast<int>(data_z[ei].size()));
        
        for (int vi = 0; vi < static_cast<int>(data_z[ei].size()); vi++)
        {
            // Get the index of variable for derivative.
            const int u_idx = data_component_idx_z[ei][vi];
            
            // Get the pointer to variable for derivative.
            double* u = data_z[ei][vi]->getPointer(u_idx);
            
            d_num_scratch_derivatives_midpoint_z_used++;
            // Declare container to store the derivative.
            if (allocate_scratch_derivatives_midpoint)
            {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
                TBOX_ASSERT(static_cast<int>(d_scratch_derivatives_midpoint_z.size()) == d_num_scratch_derivatives_midpoint_z_used - 1);
#endif
                d_scratch_derivatives_midpoint_z.push_back(HAMERS_SHARED_PTR<pdat::SideData<double> >(
                    new pdat::SideData<double>(
                        interior_box, 1, d_num_diff_ghosts, direction)));
            }
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
            else
            {
                TBOX_ASSERT(static_cast<int>(d_scratch_derivatives_midpoint_z.size()) >= d_num_scratch_derivatives_midpoint_z_used);
            }
#endif
            
            // Get the pointer to the derivative at midpoint.
            double* der_midpoint_z = d_scratch_derivatives_midpoint_z[d_num_scratch_derivatives_midpoint_z_used - 1]->getPointer(2, 0);
            
            /*
             * Get the sub-ghost cell width and ghost box dimensions of the variable.
             */
            
            hier::IntVector num_subghosts_data =
                data_z[ei][vi]->getGhostCellWidth();
            
            hier::IntVector subghostcell_dims_data =
                data_z[ei][vi]->getGhostBox().numberCells();
            
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
            
            const int num_subghosts_0_data = num_subghosts_data[0];
            const int num_subghosts_1_data = num_subghosts_data[1];
            const int num_subghosts_2_data = num_subghosts_data[2];
            const int subghostcell_dim_0_data = subghostcell_dims_data[0];
            const int subghostcell_dim_1_data = subghostcell_dims_data[1];
            
            const int num_der_midpoint_ghosts_2 = d_num_der_midpoint_ghosts[2];
            
            for (int k = -num_diff_ghosts_2 + num_der_midpoint_ghosts_2;
                k < (interior_dim_2 + 1) + num_diff_ghosts_2 - num_der_midpoint_ghosts_2;
                k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx = (i + num_diff_ghosts_0) +
                            (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                            (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                diff_ghostcell_dim_1;
                        
                        const int idx_data_BBB = (i + num_subghosts_0_data) +
                            (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                            (k - 3 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                subghostcell_dim_1_data;
                        
                        const int idx_data_BB = (i + num_subghosts_0_data) +
                            (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                            (k - 2 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                subghostcell_dim_1_data;
                        
                        const int idx_data_B = (i + num_subghosts_0_data) +
                            (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                            (k - 1 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                subghostcell_dim_1_data;
                        
                        const int idx_data_F = (i + num_subghosts_0_data) +
                            (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                            (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                subghostcell_dim_1_data;
                        
                        const int idx_data_FF = (i + num_subghosts_0_data) +
                            (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                            (k + 1 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                subghostcell_dim_1_data;
                        
                        const int idx_data_FFF = (i + num_subghosts_0_data) +
                            (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                            (k + 2 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                subghostcell_dim_1_data;
                        
                        der_midpoint_z[idx] = (a_n*(u[idx_data_F]   - u[idx_data_B]) +
                                               b_n*(u[idx_data_FF]  - u[idx_data_BB]) +
                                               c_n*(u[idx_data_FFF] - u[idx_data_BBB])
                                              )*inv_dx_2;
                    }
                }
            }
            
            derivative_z_midpoint[ei].push_back(d_scratch_derivatives_midpoint_z[d_num_scratch_derivatives_midpoint_z_used - 1]);
        }
    }
}


/*
 * Compute the derivatives in x-direction at nodes.
 */
void
DiffusiveFluxReconstructorMidpointSixthOrder::computeFirstDerivativesInXAtNode(
    hier::Patch& patch,
    std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& derivative_x_node,
    std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > >& derivative_x_node_computed,
    const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& data_x,
    const std::vector<std::vector<int> >& data_component_idx_x,
    const bool allocate_scratch_derivatives_node)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        TBOX_ASSERT(static_cast<int>(data_x[ei].size()) ==
                    static_cast<int>(data_component_idx_x[ei].size()));
    }
#endif
    
    const double a_n =  double(3)/double(4);
    const double b_n = -double(3)/double(20);
    const double c_n =  double(1)/double(60);
    
    derivative_x_node.resize(d_num_eqn);
    
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
    
    const double inv_dx_0 = double(1)/dx[0];
    
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        derivative_x_node[ei].reserve(static_cast<int>(data_x[ei].size()));
        
        for (int vi = 0; vi < static_cast<int>(data_x[ei].size()); vi++)
        {
            // Get the index of variable for derivative.
            const int u_idx = data_component_idx_x[ei][vi];
            
            if (derivative_x_node_computed.find(data_x[ei][vi]->getPointer(u_idx))
                == derivative_x_node_computed.end())
            {
                // Get the pointer to variable for derivative.
                double* u = data_x[ei][vi]->getPointer(u_idx);
                
                d_num_scratch_derivatives_node_used++;
                // Declare container to store the derivative.
                if (allocate_scratch_derivatives_node)
                {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
                    TBOX_ASSERT(static_cast<int>(d_scratch_derivatives_node.size()) == d_num_scratch_derivatives_node_used - 1);
#endif
                    d_scratch_derivatives_node.push_back(HAMERS_SHARED_PTR<pdat::CellData<double> >(
                        new pdat::CellData<double>(
                            interior_box, 1, d_num_diff_ghosts)));
                }
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
                else
                {
                    TBOX_ASSERT(static_cast<int>(d_scratch_derivatives_node.size()) >= d_num_scratch_derivatives_node_used);
                }
#endif
                
                // Get the pointer to the derivative.
                double* dudx = d_scratch_derivatives_node[d_num_scratch_derivatives_node_used - 1]->getPointer(0);
                
                /*
                 * Get the sub-ghost cell width of the variable.
                 */
                
                hier::IntVector num_subghosts_data =
                    data_x[ei][vi]->getGhostCellWidth();
                
                hier::IntVector subghostcell_dims_data =
                    data_x[ei][vi]->getGhostBox().numberCells();
                
                if (d_dim == tbox::Dimension(1))
                {
                    /*
                     * Get the dimensions and number of ghost cells.
                     */
                    
                    const int interior_dim_0 = interior_dims[0];
                    
                    const int num_diff_ghosts_0 = d_num_diff_ghosts[0];
                    
                    const int num_subghosts_0_data = num_subghosts_data[0];
                    
                    const int num_der_node_ghosts_0 = d_num_der_node_ghosts[0];
                    
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = -num_diff_ghosts_0 + num_der_node_ghosts_0;
                        i < interior_dim_0 + num_diff_ghosts_0 - num_der_node_ghosts_0;
                        i++)
                    {
                        // Compute the linear indices.
                        const int idx = i + num_diff_ghosts_0;
                        
                        const int idx_data_LLL = i - 3 + num_subghosts_0_data;
                        const int idx_data_LL  = i - 2 + num_subghosts_0_data;
                        const int idx_data_L   = i - 1 + num_subghosts_0_data;
                        const int idx_data_R   = i + 1 + num_subghosts_0_data;
                        const int idx_data_RR  = i + 2 + num_subghosts_0_data;
                        const int idx_data_RRR = i + 3 + num_subghosts_0_data;
                        
                        dudx[idx] = (a_n*(u[idx_data_R]   - u[idx_data_L]) +
                                     b_n*(u[idx_data_RR]  - u[idx_data_LL]) +
                                     c_n*(u[idx_data_RRR] - u[idx_data_LLL])
                                    )*inv_dx_0;
                    }
                }
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
                    
                    const int num_subghosts_0_data = num_subghosts_data[0];
                    const int num_subghosts_1_data = num_subghosts_data[1];
                    const int subghostcell_dim_0_data = subghostcell_dims_data[0];
                    
                    const int num_der_node_ghosts_0 = d_num_der_node_ghosts[0];
                    
                    for (int j = -num_diff_ghosts_1; j < interior_dim_1 + num_diff_ghosts_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = -num_diff_ghosts_0 + num_der_node_ghosts_0;
                            i < interior_dim_0 + num_diff_ghosts_0 - num_der_node_ghosts_0;
                            i++)
                        {
                            // Compute the linear indices.
                            const int idx = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                            
                            const int idx_data_LLL = (i - 3 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_LL = (i - 2 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_L = (i - 1 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_R = (i + 1 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_RR = (i + 2 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_RRR = (i + 3 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            dudx[idx] = (a_n*(u[idx_data_R]   - u[idx_data_L]) +
                                         b_n*(u[idx_data_RR]  - u[idx_data_LL]) +
                                         c_n*(u[idx_data_RRR] - u[idx_data_LLL])
                                        )*inv_dx_0;
                        }
                    }
                }
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
                    
                    const int num_subghosts_0_data = num_subghosts_data[0];
                    const int num_subghosts_1_data = num_subghosts_data[1];
                    const int num_subghosts_2_data = num_subghosts_data[2];
                    const int subghostcell_dim_0_data = subghostcell_dims_data[0];
                    const int subghostcell_dim_1_data = subghostcell_dims_data[1];
                    
                    const int num_der_node_ghosts_0 = d_num_der_node_ghosts[0];
                    
                    for (int k = -num_diff_ghosts_2; k < interior_dim_2 + num_diff_ghosts_2; k++)
                    {
                        for (int j = -num_diff_ghosts_1; j < interior_dim_1 + num_diff_ghosts_1; j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = -num_diff_ghosts_0 + num_der_node_ghosts_0;
                                i < interior_dim_0 + num_diff_ghosts_0 - num_der_node_ghosts_0;
                                i++)
                            {
                                // Compute the linear indices.
                                const int idx = (i + num_diff_ghosts_0) +
                                    (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                    (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                        diff_ghostcell_dim_1;
                                
                                const int idx_data_LLL = (i - 3 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_LL = (i - 2 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_L = (i - 1 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_R = (i + 1 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_RR = (i + 2 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_RRR = (i + 3 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                dudx[idx] = (a_n*(u[idx_data_R]   - u[idx_data_L]) +
                                             b_n*(u[idx_data_RR]  - u[idx_data_LL]) +
                                             c_n*(u[idx_data_RRR] - u[idx_data_LLL])
                                            )*inv_dx_0;
                            }
                        }
                    }
                }
                
                std::pair<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivative_pair(
                    u,
                    d_scratch_derivatives_node[d_num_scratch_derivatives_node_used - 1]);
                
                derivative_x_node_computed.insert(derivative_pair);
            }
            
            derivative_x_node[ei].push_back(
                derivative_x_node_computed.find(data_x[ei][vi]->getPointer(u_idx))->
                    second);
        }
    }
}


/*
 * Compute the derivatives in y-direction at nodes.
 */
void
DiffusiveFluxReconstructorMidpointSixthOrder::computeFirstDerivativesInYAtNode(
    hier::Patch& patch,
    std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& derivative_y_node,
    std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > >& derivative_y_node_computed,
    const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& data_y,
    const std::vector<std::vector<int> >& data_component_idx_y,
    const bool allocate_scratch_derivatives_node)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        TBOX_ASSERT(static_cast<int>(data_y[ei].size()) ==
                    static_cast<int>(data_component_idx_y[ei].size()));
    }
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": DiffusiveFluxReconstructorMidpointSixthOrder::"
            << "computeFirstDerivativesInYAtNode()\n"
            << "There isn't y-direction for 1D problem."
            << std::endl);
    }
#endif
    
    const double a_n =  double(3)/double(4);
    const double b_n = -double(3)/double(20);
    const double c_n =  double(1)/double(60);
    
    derivative_y_node.resize(d_num_eqn);
    
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
    
    const double inv_dx_1 = double(1)/dx[1];
    
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        derivative_y_node[ei].reserve(static_cast<int>(data_y[ei].size()));
        
        for (int vi = 0; vi < static_cast<int>(data_y[ei].size()); vi++)
        {
            // Get the index of variable for derivative.
            const int u_idx = data_component_idx_y[ei][vi];
            
            if (derivative_y_node_computed.find(data_y[ei][vi]->getPointer(u_idx))
                == derivative_y_node_computed.end())
            {
                // Get the pointer to variable for derivative.
                double* u = data_y[ei][vi]->getPointer(u_idx);
                
                d_num_scratch_derivatives_node_used++;
                // Declare container to store the derivative.
                if (allocate_scratch_derivatives_node)
                {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
                    TBOX_ASSERT(static_cast<int>(d_scratch_derivatives_node.size()) == d_num_scratch_derivatives_node_used - 1);
#endif
                    d_scratch_derivatives_node.push_back(HAMERS_SHARED_PTR<pdat::CellData<double> >(
                        new pdat::CellData<double>(
                            interior_box, 1, d_num_diff_ghosts)));
                }
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
                else
                {
                    TBOX_ASSERT(static_cast<int>(d_scratch_derivatives_node.size()) >= d_num_scratch_derivatives_node_used);
                }
#endif
                
                // Get the pointer to the derivative.
                double* dudy = d_scratch_derivatives_node[d_num_scratch_derivatives_node_used - 1]->getPointer(0);
                
                /*
                 * Get the sub-ghost cell width and ghost box dimensions of the variable.
                 */
                
                hier::IntVector num_subghosts_data =
                    data_y[ei][vi]->getGhostCellWidth();
                
                hier::IntVector subghostcell_dims_data =
                    data_y[ei][vi]->getGhostBox().numberCells();
                
                if (d_dim == tbox::Dimension(2))
                {
                    /*
                     * Get the dimensions and number of ghost cells.
                     */
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    
                    const int num_diff_ghosts_0 = d_num_diff_ghosts[0];
                    const int num_diff_ghosts_1 = d_num_diff_ghosts[1];
                    const int diff_ghostcell_dim_0 = diff_ghostcell_dims[0];
                    
                    const int num_subghosts_0_data = num_subghosts_data[0];
                    const int num_subghosts_1_data = num_subghosts_data[1];
                    const int subghostcell_dim_0_data = subghostcell_dims_data[0];
                    
                    const int num_der_node_ghosts_1 = d_num_der_node_ghosts[1];
                    
                    for (int j = -num_diff_ghosts_1 + num_der_node_ghosts_1;
                        j < interior_dim_1 + num_diff_ghosts_1 - num_der_node_ghosts_1;
                        j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = -num_diff_ghosts_0; i < interior_dim_0 + num_diff_ghosts_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                            
                            const int idx_data_BBB = (i + num_subghosts_0_data) +
                                (j - 3 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_BB = (i + num_subghosts_0_data) +
                                (j - 2 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_B = (i + num_subghosts_0_data) +
                                (j - 1 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_T = (i + num_subghosts_0_data) +
                                (j + 1 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_TT = (i + num_subghosts_0_data) +
                                (j + 2 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_TTT = (i + num_subghosts_0_data) +
                                (j + 3 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            dudy[idx] = (a_n*(u[idx_data_T]   - u[idx_data_B]) +
                                         b_n*(u[idx_data_TT]  - u[idx_data_BB]) +
                                         c_n*(u[idx_data_TTT] - u[idx_data_BBB])
                                        )*inv_dx_1;
                        }
                    }
                }
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
                    
                    const int num_subghosts_0_data = num_subghosts_data[0];
                    const int num_subghosts_1_data = num_subghosts_data[1];
                    const int num_subghosts_2_data = num_subghosts_data[2];
                    const int subghostcell_dim_0_data = subghostcell_dims_data[0];
                    const int subghostcell_dim_1_data = subghostcell_dims_data[1];
                    
                    const int num_der_node_ghosts_1 = d_num_der_node_ghosts[1];
                    
                    for (int k = -num_diff_ghosts_2; k < interior_dim_2 + num_diff_ghosts_2; k++)
                    {
                        for (int j = -num_diff_ghosts_1 + num_der_node_ghosts_1;
                            j < interior_dim_1 + num_diff_ghosts_1 - num_der_node_ghosts_1;
                            j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = -num_diff_ghosts_0; i < interior_dim_0 + num_diff_ghosts_0; i++)
                            {
                                // Compute the linear indices.
                                const int idx = (i + num_diff_ghosts_0) +
                                    (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                    (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                        diff_ghostcell_dim_1;
                                
                                const int idx_data_BBB = (i + num_subghosts_0_data) +
                                    (j - 3 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_BB = (i + num_subghosts_0_data) +
                                    (j - 2 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_B = (i + num_subghosts_0_data) +
                                    (j - 1 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_T = (i + num_subghosts_0_data) +
                                    (j + 1 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_TT = (i + num_subghosts_0_data) +
                                    (j + 2 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_TTT = (i + num_subghosts_0_data) +
                                    (j + 3 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                dudy[idx] = (a_n*(u[idx_data_T]   - u[idx_data_B]) +
                                             b_n*(u[idx_data_TT]  - u[idx_data_BB]) +
                                             c_n*(u[idx_data_TTT] - u[idx_data_BBB])
                                            )*inv_dx_1;
                            }
                        }
                    }
                }
                
                std::pair<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivative_pair(
                    u,
                    d_scratch_derivatives_node[d_num_scratch_derivatives_node_used - 1]);
                
                derivative_y_node_computed.insert(derivative_pair);
            }
            
            derivative_y_node[ei].push_back(
                derivative_y_node_computed.find(data_y[ei][vi]->getPointer(u_idx))->
                    second);
        }
    }
}


/*
 * Compute the derivatives in z-direction at nodes.
 */
void
DiffusiveFluxReconstructorMidpointSixthOrder::computeFirstDerivativesInZAtNode(
    hier::Patch& patch,
    std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& derivative_z_node,
    std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > >& derivative_z_node_computed,
    const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& data_z,
    const std::vector<std::vector<int> >& data_component_idx_z,
    const bool allocate_scratch_derivatives_node)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        TBOX_ASSERT(static_cast<int>(data_z[ei].size()) ==
                    static_cast<int>(data_component_idx_z[ei].size()));
    }
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": DiffusiveFluxReconstructorMidpointSixthOrder::"
            << "computeFirstDerivativesInZAtNode()\n"
            << "There isn't z-direction for 1D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": DiffusiveFluxReconstructorMidpointSixthOrder::"
            << "computeFirstDerivativesInZAtNode()\n"
            << "There isn't z-direction for 2D problem."
            << std::endl);
    }
#endif
    
    const double a_n =  double(3)/double(4);
    const double b_n = -double(3)/double(20);
    const double c_n =  double(1)/double(60);
    
    derivative_z_node.resize(d_num_eqn);
    
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
    
    const double inv_dx_2 = double(1)/dx[2];
    
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        derivative_z_node[ei].reserve(static_cast<int>(data_z[ei].size()));
        
        for (int vi = 0; vi < static_cast<int>(data_z[ei].size()); vi++)
        {
            // Get the index of variable for derivative.
            const int u_idx = data_component_idx_z[ei][vi];
            
            if (derivative_z_node_computed.find(data_z[ei][vi]->getPointer(u_idx))
                == derivative_z_node_computed.end())
            {
                // Get the pointer to variable for derivative.
                double* u = data_z[ei][vi]->getPointer(u_idx);
                
                d_num_scratch_derivatives_node_used++;
                // Declare container to store the derivative.
                if (allocate_scratch_derivatives_node)
                {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
                    TBOX_ASSERT(static_cast<int>(d_scratch_derivatives_node.size()) == d_num_scratch_derivatives_node_used - 1);
#endif
                    d_scratch_derivatives_node.push_back(HAMERS_SHARED_PTR<pdat::CellData<double> >(
                        new pdat::CellData<double>(
                            interior_box, 1, d_num_diff_ghosts)));
                }
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
                else
                {
                    TBOX_ASSERT(static_cast<int>(d_scratch_derivatives_node.size()) >= d_num_scratch_derivatives_node_used);
                }
#endif
                
                // Get the pointer to the derivative.
                double* dudz = d_scratch_derivatives_node[d_num_scratch_derivatives_node_used - 1]->getPointer(0);
                
                /*
                 * Get the sub-ghost cell width and ghost box dimensions of the variable.
                 */
                
                hier::IntVector num_subghosts_data =
                    data_z[ei][vi]->getGhostCellWidth();
                
                hier::IntVector subghostcell_dims_data =
                    data_z[ei][vi]->getGhostBox().numberCells();
                
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
                
                const int num_subghosts_0_data = num_subghosts_data[0];
                const int num_subghosts_1_data = num_subghosts_data[1];
                const int num_subghosts_2_data = num_subghosts_data[2];
                const int subghostcell_dim_0_data = subghostcell_dims_data[0];
                const int subghostcell_dim_1_data = subghostcell_dims_data[1];
                
                const int num_der_node_ghosts_2 = d_num_der_node_ghosts[2];
                
                for (int k = -num_diff_ghosts_2 + num_der_node_ghosts_2;
                    k < interior_dim_2 + num_diff_ghosts_2 - num_der_node_ghosts_2;
                    k++)
                {
                    for (int j = -num_diff_ghosts_1; j < interior_dim_1 + num_diff_ghosts_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = -num_diff_ghosts_0; i < interior_dim_0 + num_diff_ghosts_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_data_BBB = (i + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                (k - 3 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                    subghostcell_dim_1_data;
                            
                            const int idx_data_BB = (i + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                (k - 2 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                    subghostcell_dim_1_data;
                            
                            const int idx_data_B = (i + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                (k - 1 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                    subghostcell_dim_1_data;
                            
                            const int idx_data_F = (i + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                (k + 1 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                    subghostcell_dim_1_data;
                            
                            const int idx_data_FF = (i + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                (k + 2 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                    subghostcell_dim_1_data;
                            
                            const int idx_data_FFF = (i + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                (k + 3 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                    subghostcell_dim_1_data;
                            
                            dudz[idx] = (a_n*(u[idx_data_F]   - u[idx_data_B]) +
                                         b_n*(u[idx_data_FF]  - u[idx_data_BB]) +
                                         c_n*(u[idx_data_FFF] - u[idx_data_BBB])
                                        )*inv_dx_2;
                        }
                    }
                }
                
                std::pair<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivative_pair(
                    u,
                    d_scratch_derivatives_node[d_num_scratch_derivatives_node_used - 1]);
                
                derivative_z_node_computed.insert(derivative_pair);
            }
            
            derivative_z_node[ei].push_back(
                derivative_z_node_computed.find(data_z[ei][vi]->getPointer(u_idx))->
                    second);
        }
    }
}


/*
 * Interpolate the diffusivities from nodes to midpoints.
 */
void
DiffusiveFluxReconstructorMidpointSixthOrder::interpolateDiffusivitiesFromNodeToMidpoint(
    hier::Patch& patch,
    std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& var_side_data_for_diffusivities,
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& var_cell_data_for_diffusivities,
    const std::vector<int>& var_cell_data_for_diffusivities_component_idx,
    const bool allocate_scratch_diffusivities_midpoint)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(var_cell_data_for_diffusivities.size()) ==
        static_cast<int>(var_cell_data_for_diffusivities_component_idx.size()));
#endif
    
    const double a_n =  double(75)/double(128);
    const double b_n = -double(25)/double(256);
    const double c_n =  double(3)/double(256);
    
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
        TBOX_ASSERT(var_cell_data_for_diffusivities[vi]->->getGhostCellWidth() >= d_num_diff_ghosts);
#endif
        
        // Get the pointer to variable for interpoation.
        double* u = var_cell_data_for_diffusivities[vi]->getPointer(var_cell_data_for_diffusivities_component_idx[vi]);
        
        d_num_scratch_diffusivities_midpoint_used++;
        // Declare container to store the diffusivities.
        if (allocate_scratch_diffusivities_midpoint)
        {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
            TBOX_ASSERT(static_cast<int>(d_scratch_diffusivities_midpoint.size()) == d_num_scratch_diffusivities_midpoint_used - 1);
#endif
            d_scratch_diffusivities_midpoint.push_back(HAMERS_SHARED_PTR<pdat::SideData<double> >(
                new pdat::SideData<double>(
                    interior_box, 1, d_num_diff_ghosts)));
        }
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
        else
        {
            TBOX_ASSERT(static_cast<int>(d_scratch_diffusivities_midpoint.size()) >= d_num_scratch_diffusivities_midpoint_used);
        }
#endif
        
        /*
         * Get the sub-ghost cell width of the variable.
         */
        
        hier::IntVector num_subghosts_data =
            var_cell_data_for_diffusivities[vi]->getGhostCellWidth();
        
        hier::IntVector subghostcell_dims_data =
            var_cell_data_for_diffusivities[vi]->getGhostBox().numberCells();
            
        if (d_dim == tbox::Dimension(1))
        {
            // Get the pointer to the diffusivity at midpoint.
            double* diffus_midpoint_x =
                d_scratch_diffusivities_midpoint[d_num_scratch_diffusivities_midpoint_used - 1]->getPointer(0, 0);
            
            /*
             * Get the dimensions and number of ghost cells.
             */
            
            const int interior_dim_0 = interior_dims[0];
            
            const int num_diff_ghosts_0 = d_num_diff_ghosts[0];
            
            const int num_subghosts_0_data = num_subghosts_data[0];
            
            const int num_interp_midpoint_ghosts_0 = d_num_interp_midpoint_ghosts[0];
            
            // Interpolation in x-direction.
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = -num_diff_ghosts_0 + num_interp_midpoint_ghosts_0;
                i < (interior_dim_0 + 1) + num_diff_ghosts_0 - num_interp_midpoint_ghosts_0;
                i++)
            {
                // Compute the linear indices.
                const int idx = i + num_diff_ghosts_0;
                
                const int idx_data_LLL = i - 3 + num_subghosts_0_data;
                const int idx_data_LL  = i - 2 + num_subghosts_0_data;
                const int idx_data_L   = i - 1 + num_subghosts_0_data;
                const int idx_data_R   = i + 0 + num_subghosts_0_data;
                const int idx_data_RR  = i + 1 + num_subghosts_0_data;
                const int idx_data_RRR = i + 2 + num_subghosts_0_data;
                
                diffus_midpoint_x[idx] = (a_n*(u[idx_data_R]   + u[idx_data_L]) +
                                          b_n*(u[idx_data_RR]  + u[idx_data_LL]) +
                                          c_n*(u[idx_data_RRR] + u[idx_data_LLL])
                                         );
            }
            
        }
        else if (d_dim == tbox::Dimension(2))
        {
            // Get the pointers to the diffusivities at midpoints.
            double* diffus_midpoint_x =
                d_scratch_diffusivities_midpoint[d_num_scratch_diffusivities_midpoint_used - 1]->getPointer(0, 0);
            double* diffus_midpoint_y =
                d_scratch_diffusivities_midpoint[d_num_scratch_diffusivities_midpoint_used - 1]->getPointer(1, 0);
            
            /*
             * Get the dimensions and number of ghost cells.
             */
            
            const int interior_dim_0 = interior_dims[0];
            const int interior_dim_1 = interior_dims[1];
            
            const int num_diff_ghosts_0 = d_num_diff_ghosts[0];
            const int num_diff_ghosts_1 = d_num_diff_ghosts[1];
            const int diff_ghostcell_dim_0 = diff_ghostcell_dims[0];
            
            const int num_subghosts_0_data = num_subghosts_data[0];
            const int num_subghosts_1_data = num_subghosts_data[1];
            const int subghostcell_dim_0_data = subghostcell_dims_data[0];
            
            const int num_interp_midpoint_ghosts_0 = d_num_interp_midpoint_ghosts[0];
            const int num_interp_midpoint_ghosts_1 = d_num_interp_midpoint_ghosts[1];
            
            // Interpolation in x-direction.
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = -num_diff_ghosts_0 + num_interp_midpoint_ghosts_0;
                    i < (interior_dim_0 + 1) + num_diff_ghosts_0 - num_interp_midpoint_ghosts_0;
                    i++)
                {
                    // Compute the linear indices.
                    const int idx = (i + num_diff_ghosts_0) +
                        (j + num_diff_ghosts_1)*(diff_ghostcell_dim_0 + 1);
                    
                    const int idx_data_LLL = (i - 3 + num_subghosts_0_data) +
                        (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                    
                    const int idx_data_LL = (i - 2 + num_subghosts_0_data) +
                        (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                    
                    const int idx_data_L = (i - 1 + num_subghosts_0_data) +
                        (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                    
                    const int idx_data_R = (i + num_subghosts_0_data) +
                        (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                    
                    const int idx_data_RR = (i + 1 + num_subghosts_0_data) +
                        (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                    
                    const int idx_data_RRR = (i + 2 + num_subghosts_0_data) +
                        (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                    
                    diffus_midpoint_x[idx] = (a_n*(u[idx_data_R]   + u[idx_data_L]) +
                                              b_n*(u[idx_data_RR]  + u[idx_data_LL]) +
                                              c_n*(u[idx_data_RRR] + u[idx_data_LLL])
                                             );
                }
            }
            
            // Interpolation in y-direction.
            for (int j = -num_diff_ghosts_1 + num_interp_midpoint_ghosts_1;
                j < (interior_dim_1 + 1) + num_diff_ghosts_1 - num_interp_midpoint_ghosts_1;
                j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx = (i + num_diff_ghosts_0) +
                        (j + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                    
                    const int idx_data_BBB = (i + num_subghosts_0_data) +
                        (j - 3 + num_subghosts_1_data)*subghostcell_dim_0_data;
                    
                    const int idx_data_BB = (i + num_subghosts_0_data) +
                        (j - 2 + num_subghosts_1_data)*subghostcell_dim_0_data;
                    
                    const int idx_data_B = (i + num_subghosts_0_data) +
                        (j - 1 + num_subghosts_1_data)*subghostcell_dim_0_data;
                    
                    const int idx_data_T = (i + num_subghosts_0_data) +
                        (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                    
                    const int idx_data_TT = (i + num_subghosts_0_data) +
                        (j + 1 + num_subghosts_1_data)*subghostcell_dim_0_data;
                    
                    const int idx_data_TTT = (i + num_subghosts_0_data) +
                        (j + 2 + num_subghosts_1_data)*subghostcell_dim_0_data;
                    
                    diffus_midpoint_y[idx] = (a_n*(u[idx_data_T]   + u[idx_data_B]) +
                                              b_n*(u[idx_data_TT]  + u[idx_data_BB]) +
                                              c_n*(u[idx_data_TTT] + u[idx_data_BBB])
                                             );
                }
            }
        }
        else if (d_dim == tbox::Dimension(3))
        {
            // Get the pointers to the diffusivities at midpoints.
            double* diffus_midpoint_x =
                d_scratch_diffusivities_midpoint[d_num_scratch_diffusivities_midpoint_used - 1]->getPointer(0, 0);
            double* diffus_midpoint_y =
                d_scratch_diffusivities_midpoint[d_num_scratch_diffusivities_midpoint_used - 1]->getPointer(1, 0);
            double* diffus_midpoint_z =
                d_scratch_diffusivities_midpoint[d_num_scratch_diffusivities_midpoint_used - 1]->getPointer(2, 0);
            
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
            
            const int num_subghosts_0_data = num_subghosts_data[0];
            const int num_subghosts_1_data = num_subghosts_data[1];
            const int num_subghosts_2_data = num_subghosts_data[2];
            const int subghostcell_dim_0_data = subghostcell_dims_data[0];
            const int subghostcell_dim_1_data = subghostcell_dims_data[1];
            
            const int num_interp_midpoint_ghosts_0 = d_num_interp_midpoint_ghosts[0];
            const int num_interp_midpoint_ghosts_1 = d_num_interp_midpoint_ghosts[1];
            const int num_interp_midpoint_ghosts_2 = d_num_interp_midpoint_ghosts[2];
            
            // Interpolation in x-direction.
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = -num_diff_ghosts_0 + num_interp_midpoint_ghosts_0;
                        i < (interior_dim_0 + 1) + num_diff_ghosts_0 - num_interp_midpoint_ghosts_0;
                        i++)
                    {
                        // Compute the linear indices.
                        const int idx = (i + num_diff_ghosts_0) +
                            (j + num_diff_ghosts_1)*(diff_ghostcell_dim_0 + 1) +
                            (k + num_diff_ghosts_2)*(diff_ghostcell_dim_0 + 1)*
                                diff_ghostcell_dim_1;
                        
                        const int idx_data_LLL = (i - 3 + num_subghosts_0_data) +
                            (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                            (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                subghostcell_dim_1_data;
                        
                        const int idx_data_LL = (i - 2 + num_subghosts_0_data) +
                            (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                            (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                subghostcell_dim_1_data;
                        
                        const int idx_data_L = (i - 1 + num_subghosts_0_data) +
                            (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                            (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                subghostcell_dim_1_data;
                        
                        const int idx_data_R = (i + num_subghosts_0_data) +
                            (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                            (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                subghostcell_dim_1_data;
                        
                        const int idx_data_RR = (i + 1 + num_subghosts_0_data) +
                            (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                            (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                subghostcell_dim_1_data;
                        
                        const int idx_data_RRR = (i + 2 + num_subghosts_0_data) +
                            (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                            (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                subghostcell_dim_1_data;
                        
                        diffus_midpoint_x[idx] = (a_n*(u[idx_data_R]   + u[idx_data_L]) +
                                                  b_n*(u[idx_data_RR]  + u[idx_data_LL]) +
                                                  c_n*(u[idx_data_RRR] + u[idx_data_LLL])
                                                 );
                    }
                }
            }
            
            // Interpolation in y-direction.
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = -num_diff_ghosts_1 + num_interp_midpoint_ghosts_1;
                    j < (interior_dim_1 + 1) + num_diff_ghosts_1 - num_interp_midpoint_ghosts_1;
                    j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx = (i + num_diff_ghosts_0) +
                            (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                            (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                (diff_ghostcell_dim_1 + 1);
                        
                        const int idx_data_BBB = (i + num_subghosts_0_data) +
                            (j - 3 + num_subghosts_1_data)*subghostcell_dim_0_data +
                            (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                subghostcell_dim_1_data;
                        
                        const int idx_data_BB = (i + num_subghosts_0_data) +
                            (j - 2 + num_subghosts_1_data)*subghostcell_dim_0_data +
                            (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                subghostcell_dim_1_data;
                        
                        const int idx_data_B = (i + num_subghosts_0_data) +
                            (j - 1 + num_subghosts_1_data)*subghostcell_dim_0_data +
                            (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                subghostcell_dim_1_data;
                        
                        const int idx_data_T = (i + num_subghosts_0_data) +
                            (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                            (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                subghostcell_dim_1_data;
                        
                        const int idx_data_TT = (i + num_subghosts_0_data) +
                            (j + 1 + num_subghosts_1_data)*subghostcell_dim_0_data +
                            (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                subghostcell_dim_1_data;
                        
                        const int idx_data_TTT = (i + num_subghosts_0_data) +
                            (j + 2 + num_subghosts_1_data)*subghostcell_dim_0_data +
                            (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                subghostcell_dim_1_data;
                        
                        diffus_midpoint_y[idx] = (a_n*(u[idx_data_T]   + u[idx_data_B]) +
                                                  b_n*(u[idx_data_TT]  + u[idx_data_BB]) +
                                                  c_n*(u[idx_data_TTT] + u[idx_data_BBB])
                                                 );
                    }
                }
            }
            
            // Interpolation in z-direction.
            for (int k = -num_diff_ghosts_2 + num_interp_midpoint_ghosts_2;
                k < (interior_dim_2 + 1) + num_diff_ghosts_2 - num_interp_midpoint_ghosts_2;
                k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx = (i + num_diff_ghosts_0) +
                            (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                            (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                diff_ghostcell_dim_1;
                        
                        const int idx_data_BBB = (i + num_subghosts_0_data) +
                            (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                            (k - 3 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                subghostcell_dim_1_data;
                        
                        const int idx_data_BB = (i + num_subghosts_0_data) +
                            (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                            (k - 2 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                subghostcell_dim_1_data;
                        
                        const int idx_data_B = (i + num_subghosts_0_data) +
                            (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                            (k - 1 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                subghostcell_dim_1_data;
                        
                        const int idx_data_F = (i + num_subghosts_0_data) +
                            (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                            (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                subghostcell_dim_1_data;
                        
                        const int idx_data_FF = (i + num_subghosts_0_data) +
                            (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                            (k + 1 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                subghostcell_dim_1_data;
                        
                        const int idx_data_FFF = (i + num_subghosts_0_data) +
                            (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                            (k + 2 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                subghostcell_dim_1_data;
                        
                        diffus_midpoint_z[idx] = (a_n*(u[idx_data_F]   + u[idx_data_B]) +
                                                  b_n*(u[idx_data_FF]  + u[idx_data_BB]) +
                                                  c_n*(u[idx_data_FFF] + u[idx_data_BBB])
                                                 );
                    }
                }
            }
        }
        
        var_side_data_for_diffusivities.push_back(d_scratch_diffusivities_midpoint[d_num_scratch_diffusivities_midpoint_used - 1]);
    }
}


/*
 * Interpolate the derivatives from nodes to midpoints in x-direction.
 */
void
DiffusiveFluxReconstructorMidpointSixthOrder::interpolateDerivativesFromNodeToMidpointX(
    hier::Patch& patch,
    std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > >& derivative_midpoint_x,
    const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& derivative_node,
    const bool allocate_scratch_derivatives_midpoint)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(derivative_node.size()) == d_num_eqn);
#endif
    
    const double a_n =  double(75)/double(128);
    const double b_n = -double(25)/double(256);
    const double c_n =  double(3)/double(256);
    
    derivative_midpoint_x.resize(d_num_eqn);
    
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
    
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        derivative_midpoint_x[ei].reserve(static_cast<int>(derivative_node[ei].size()));
        
        for (int vi = 0; vi < static_cast<int>(derivative_node[ei].size()); vi++)
        {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
            TBOX_ASSERT(derivative_node[ei][vi]->getGhostCellWidth() == d_num_diff_ghosts);
#endif
            
            // Get the pointer to variable for interpoation.
            double* der_node = derivative_node[ei][vi]->getPointer(0);
            
            d_num_scratch_derivatives_midpoint_x_used++;
            // Declare container to store the derivative.
            if (allocate_scratch_derivatives_midpoint)
            {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
                TBOX_ASSERT(static_cast<int>(d_scratch_derivatives_midpoint_x.size()) == d_num_scratch_derivatives_midpoint_x_used - 1);
#endif
                d_scratch_derivatives_midpoint_x.push_back(HAMERS_SHARED_PTR<pdat::SideData<double> >(
                    new pdat::SideData<double>(
                        interior_box, 1, d_num_diff_ghosts, direction)));
            }
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
            else
            {
                TBOX_ASSERT(static_cast<int>(d_scratch_derivatives_midpoint_x.size()) >= d_num_scratch_derivatives_midpoint_x_used);
            }
#endif
            
            // Get the pointer to the derivative at midpoint.
            double* der_midpoint_x = d_scratch_derivatives_midpoint_x[d_num_scratch_derivatives_midpoint_x_used - 1]->getPointer(0, 0);
            
            if (d_dim == tbox::Dimension(1))
            {
                /*
                 * Get the dimensions and number of ghost cells.
                 */
                
                const int interior_dim_0 = interior_dims[0];
                
                const int num_diff_ghosts_0 = d_num_diff_ghosts[0];
                
                const int num_interp_midpoint_ghosts_0 = d_num_interp_midpoint_ghosts[0];
                
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = -num_diff_ghosts_0 + num_interp_midpoint_ghosts_0;
                    i < (interior_dim_0 + 1) + num_diff_ghosts_0 - num_interp_midpoint_ghosts_0;
                    i++)
                {
                    // Compute the linear indices.
                    const int idx = i + num_diff_ghosts_0;
                    
                    const int idx_node_LLL = i - 3 + num_diff_ghosts_0;
                    const int idx_node_LL  = i - 2 + num_diff_ghosts_0;
                    const int idx_node_L   = i - 1 + num_diff_ghosts_0;
                    const int idx_node_R   = i + 0 + num_diff_ghosts_0;
                    const int idx_node_RR  = i + 1 + num_diff_ghosts_0;
                    const int idx_node_RRR = i + 2 + num_diff_ghosts_0;
                    
                    der_midpoint_x[idx] = (a_n*(der_node[idx_node_R]   + der_node[idx_node_L]) +
                                           b_n*(der_node[idx_node_RR]  + der_node[idx_node_LL]) +
                                           c_n*(der_node[idx_node_RRR] + der_node[idx_node_LLL])
                                          );
                }
            }
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
                
                const int num_interp_midpoint_ghosts_0 = d_num_interp_midpoint_ghosts[0];
                
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = -num_diff_ghosts_0 + num_interp_midpoint_ghosts_0;
                        i < (interior_dim_0 + 1) + num_diff_ghosts_0 - num_interp_midpoint_ghosts_0;
                        i++)
                    {
                        // Compute the linear indices.
                        const int idx = (i + num_diff_ghosts_0) +
                            (j + num_diff_ghosts_1)*(diff_ghostcell_dim_0 + 1);
                        
                        const int idx_node_LLL = (i - 3 + num_diff_ghosts_0) +
                            (j + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                        
                        const int idx_node_LL = (i - 2 + num_diff_ghosts_0) +
                            (j + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                        
                        const int idx_node_L = (i - 1 + num_diff_ghosts_0) +
                            (j + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                        
                        const int idx_node_R = (i + num_diff_ghosts_0) +
                            (j + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                        
                        const int idx_node_RR = (i + 1 + num_diff_ghosts_0) +
                            (j + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                        
                        const int idx_node_RRR = (i + 2 + num_diff_ghosts_0) +
                            (j + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                        
                        der_midpoint_x[idx] = (a_n*(der_node[idx_node_R]   + der_node[idx_node_L]) +
                                               b_n*(der_node[idx_node_RR]  + der_node[idx_node_LL]) +
                                               c_n*(der_node[idx_node_RRR] + der_node[idx_node_LLL])
                                              );
                    }
                }
            }
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
                
                const int num_interp_midpoint_ghosts_0 = d_num_interp_midpoint_ghosts[0];
                
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = -num_diff_ghosts_0 + num_interp_midpoint_ghosts_0;
                            i < (interior_dim_0 + 1) + num_diff_ghosts_0 - num_interp_midpoint_ghosts_0;
                            i++)
                        {
                            // Compute the linear indices.
                            const int idx = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*(diff_ghostcell_dim_0 + 1) +
                                (k + num_diff_ghosts_2)*(diff_ghostcell_dim_0 + 1)*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_LLL = (i - 3 + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_LL = (i - 2 + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_L = (i - 1 + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_R = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_RR = (i + 1 + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_RRR = (i + 2 + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            der_midpoint_x[idx] = (a_n*(der_node[idx_node_R]   + der_node[idx_node_L]) +
                                                   b_n*(der_node[idx_node_RR]  + der_node[idx_node_LL]) +
                                                   c_n*(der_node[idx_node_RRR] + der_node[idx_node_LLL])
                                                  );
                        }
                    }
                }
            }
            
            derivative_midpoint_x[ei].push_back(d_scratch_derivatives_midpoint_x[d_num_scratch_derivatives_midpoint_x_used - 1]);
        }
    }
}


/*
 * Interpolate the derivatives from nodes to midpoints in y-direction.
 */
void
DiffusiveFluxReconstructorMidpointSixthOrder::interpolateDerivativesFromNodeToMidpointY(
    hier::Patch& patch,
    std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > >& derivative_midpoint_y,
    const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& derivative_node,
    const bool allocate_scratch_derivatives_midpoint)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(derivative_node.size()) == d_num_eqn);
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": DiffusiveFluxReconstructorMidpointSixthOrder::"
            << "interpolateDerivativesFromNodeToMidpointY()\n"
            << "There isn't y-direction for 1D problem."
            << std::endl);
    }
#endif
    
    const double a_n =  double(75)/double(128);
    const double b_n = -double(25)/double(256);
    const double c_n =  double(3)/double(256);
    
    derivative_midpoint_y.resize(d_num_eqn);
    
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
    
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        derivative_midpoint_y[ei].reserve(static_cast<int>(derivative_node[ei].size()));
        
        for (int vi = 0; vi < static_cast<int>(derivative_node[ei].size()); vi++)
        {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
            TBOX_ASSERT(derivative_node[ei][vi]->getGhostCellWidth() == d_num_diff_ghosts);
#endif
            
            // Get the pointer to variable for interpoation.
            double* der_node = derivative_node[ei][vi]->getPointer(0);
            
            d_num_scratch_derivatives_midpoint_y_used++;
            // Declare container to store the derivative.
            if (allocate_scratch_derivatives_midpoint)
            {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
                TBOX_ASSERT(static_cast<int>(d_scratch_derivatives_midpoint_y.size()) == d_num_scratch_derivatives_midpoint_y_used - 1);
#endif
                d_scratch_derivatives_midpoint_y.push_back(HAMERS_SHARED_PTR<pdat::SideData<double> >(
                    new pdat::SideData<double>(
                        interior_box, 1, d_num_diff_ghosts, direction)));
            }
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
            else
            {
                TBOX_ASSERT(static_cast<int>(d_scratch_derivatives_midpoint_y.size()) >= d_num_scratch_derivatives_midpoint_y_used);
            }
#endif
            
            // Get the pointer to the derivative at midpoint.
            double* der_midpoint_y = d_scratch_derivatives_midpoint_y[d_num_scratch_derivatives_midpoint_y_used - 1]->getPointer(1, 0);
            
            if (d_dim == tbox::Dimension(2))
            {
                /*
                 * Get the dimensions and number of ghost cells.
                 */
                
                const int interior_dim_0 = interior_dims[0];
                const int interior_dim_1 = interior_dims[1];
                
                const int num_diff_ghosts_0 = d_num_diff_ghosts[0];
                const int num_diff_ghosts_1 = d_num_diff_ghosts[1];
                const int diff_ghostcell_dim_0 = diff_ghostcell_dims[0];
                
                const int num_interp_midpoint_ghosts_1 = d_num_interp_midpoint_ghosts[1];
                
                for (int j = -num_diff_ghosts_1 + num_interp_midpoint_ghosts_1;
                    j < (interior_dim_1 + 1) + num_diff_ghosts_1 - num_interp_midpoint_ghosts_1;
                    j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx = (i + num_diff_ghosts_0) +
                            (j + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                        
                        const int idx_node_BBB = (i + num_diff_ghosts_0) +
                            (j - 3 + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                        
                        const int idx_node_BB = (i + num_diff_ghosts_0) +
                            (j - 2 + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                        
                        const int idx_node_B = (i + num_diff_ghosts_0) +
                            (j - 1 + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                        
                        const int idx_node_T = (i + num_diff_ghosts_0) +
                            (j + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                        
                        const int idx_node_TT = (i + num_diff_ghosts_0) +
                            (j + 1 + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                        
                        const int idx_node_TTT = (i + num_diff_ghosts_0) +
                            (j + 2 + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                        
                        der_midpoint_y[idx] = (a_n*(der_node[idx_node_T]   + der_node[idx_node_B]) +
                                               b_n*(der_node[idx_node_TT]  + der_node[idx_node_BB]) +
                                               c_n*(der_node[idx_node_TTT] + der_node[idx_node_BBB])
                                              );
                    }
                }
            }
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
                
                const int num_interp_midpoint_ghosts_1 = d_num_interp_midpoint_ghosts[1];
                
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = -num_diff_ghosts_1 + num_interp_midpoint_ghosts_1;
                        j < (interior_dim_1 + 1) + num_diff_ghosts_1 - num_interp_midpoint_ghosts_1;
                        j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    (diff_ghostcell_dim_1 + 1);
                            
                            const int idx_node_BBB = (i + num_diff_ghosts_0) +
                                (j - 3 + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_BB = (i + num_diff_ghosts_0) +
                                (j - 2 + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_B = (i + num_diff_ghosts_0) +
                                (j - 1 + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_T = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_TT = (i + num_diff_ghosts_0) +
                                (j + 1 + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_TTT = (i + num_diff_ghosts_0) +
                                (j + 2 + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            der_midpoint_y[idx] = (a_n*(der_node[idx_node_T]   + der_node[idx_node_B]) +
                                                   b_n*(der_node[idx_node_TT]  + der_node[idx_node_BB]) +
                                                   c_n*(der_node[idx_node_TTT] + der_node[idx_node_BBB])
                                                  );
                        }
                    }
                }
                
            }
            
            derivative_midpoint_y[ei].push_back(d_scratch_derivatives_midpoint_y[d_num_scratch_derivatives_midpoint_y_used - 1]);
        }
    }
}


/*
 * Interpolate the derivatives from nodes to midpoints in z-direction.
 */
void
DiffusiveFluxReconstructorMidpointSixthOrder::interpolateDerivativesFromNodeToMidpointZ(
    hier::Patch& patch,
    std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > >& derivative_midpoint_z,
    const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& derivative_node,
    const bool allocate_scratch_derivatives_midpoint)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(derivative_node.size()) == d_num_eqn);
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": DiffusiveFluxReconstructorMidpointSixthOrder::"
            << "interpolateDerivativesFromNodeToMidpointZ()\n"
            << "There isn't z-direction for 1D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": DiffusiveFluxReconstructorMidpointSixthOrder::"
            << "interpolateDerivativesFromNodeToMidpointZ()\n"
            << "There isn't z-direction for 2D problem."
            << std::endl);
    }
#endif
    
    const double a_n =  double(75)/double(128);
    const double b_n = -double(25)/double(256);
    const double c_n =  double(3)/double(256);
    
    derivative_midpoint_z.resize(d_num_eqn);
    
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
    
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        derivative_midpoint_z[ei].reserve(static_cast<int>(derivative_node[ei].size()));
        
        for (int vi = 0; vi < static_cast<int>(derivative_node[ei].size()); vi++)
        {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
            TBOX_ASSERT(derivative_node[ei][vi]->getGhostCellWidth() == d_num_diff_ghosts);
#endif
            
            // Get the pointer to variable for interpoation.
            double* der_node = derivative_node[ei][vi]->getPointer(0);
            
            d_num_scratch_derivatives_midpoint_z_used++;
            // Declare container to store the derivative.
            if (allocate_scratch_derivatives_midpoint)
            {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
                TBOX_ASSERT(static_cast<int>(d_scratch_derivatives_midpoint_z.size()) == d_num_scratch_derivatives_midpoint_z_used - 1);
#endif
                d_scratch_derivatives_midpoint_z.push_back(HAMERS_SHARED_PTR<pdat::SideData<double> >(
                    new pdat::SideData<double>(
                        interior_box, 1, d_num_diff_ghosts, direction)));
            }
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
            else
            {
                TBOX_ASSERT(static_cast<int>(d_scratch_derivatives_midpoint_z.size()) >= d_num_scratch_derivatives_midpoint_z_used);
            }
#endif
            
            // Get the pointer to the derivative at midpoint.
            double* der_midpoint_z = d_scratch_derivatives_midpoint_z[d_num_scratch_derivatives_midpoint_z_used - 1]->getPointer(2, 0);
            
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
            
            const int num_interp_midpoint_ghosts_2 = d_num_interp_midpoint_ghosts[2];
            
            for (int k = -num_diff_ghosts_2 + num_interp_midpoint_ghosts_2;
                k < (interior_dim_2 + 1) + num_diff_ghosts_2 - num_interp_midpoint_ghosts_2;
                k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx = (i + num_diff_ghosts_0) +
                            (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                            (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                diff_ghostcell_dim_1;
                        
                        const int idx_node_BBB = (i + num_diff_ghosts_0) +
                            (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                            (k - 3 + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                diff_ghostcell_dim_1;
                        
                        const int idx_node_BB = (i + num_diff_ghosts_0) +
                            (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                            (k - 2 + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                diff_ghostcell_dim_1;
                        
                        const int idx_node_B = (i + num_diff_ghosts_0) +
                            (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                            (k - 1 + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                diff_ghostcell_dim_1;
                        
                        const int idx_node_F = (i + num_diff_ghosts_0) +
                            (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                            (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                diff_ghostcell_dim_1;
                        
                        const int idx_node_FF = (i + num_diff_ghosts_0) +
                            (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                            (k + 1 + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                diff_ghostcell_dim_1;
                        
                        const int idx_node_FFF = (i + num_diff_ghosts_0) +
                            (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                            (k + 2 + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                diff_ghostcell_dim_1;
                        
                        der_midpoint_z[idx] = (a_n*(der_node[idx_node_F]   + der_node[idx_node_B]) +
                                               b_n*(der_node[idx_node_FF]  + der_node[idx_node_BB]) +
                                               c_n*(der_node[idx_node_FFF] + der_node[idx_node_BBB])
                                              );
                    }
                }
            }
            
            derivative_midpoint_z[ei].push_back(d_scratch_derivatives_midpoint_z[d_num_scratch_derivatives_midpoint_z_used - 1]);
        }
    }
}


/*
 * Reconstruct the flux using flux at midpoints in x-direction.
 */
void
DiffusiveFluxReconstructorMidpointSixthOrder::reconstructFluxX(
    hier::Patch& patch,
    HAMERS_SHARED_PTR<pdat::SideData<double> >& diffusive_flux,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& diffusive_flux_midpoint,
    const double dt)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(diffusive_flux.getDepth() == d_num_eqn);
    TBOX_ASSERT(diffusive_flux_midpoint.getDepth() == d_num_eqn);
#endif
    
    const double a_m =  double(75)/double(64);
    const double b_m = -double(25)/double(384);
    const double c_m =  double(3)/double(640);
    
    const double a_r = a_m + b_m + c_m;
    const double b_r = b_m + c_m;
    const double c_r = c_m;
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // diffusive ghost cells.
    hier::Box diff_ghost_box = interior_box;
    diff_ghost_box.grow(d_num_diff_ghosts);
    const hier::IntVector diff_ghostcell_dims = diff_ghost_box.numberCells();
    
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        double* F_face_x     = diffusive_flux->getPointer(0, ei);
        double* F_midpoint_x = diffusive_flux_midpoint->getPointer(0, ei);
        
        if (d_dim == tbox::Dimension(1))
        {
            /*
             * Get the dimensions and number of ghost cells.
             */
            
            const int interior_dim_0 = interior_dims[0];
            
            const int num_diff_ghosts_0 = d_num_diff_ghosts[0];
            
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0 + 1; i++)
            {
                // Compute the linear indices.
                const int idx_face_x = i;
                
                const int idx_midpoint_LL  = i - 2 + num_diff_ghosts_0;
                const int idx_midpoint_L   = i - 1 + num_diff_ghosts_0;
                const int idx_midpoint     = i + 0 + num_diff_ghosts_0;
                const int idx_midpoint_R   = i + 1 + num_diff_ghosts_0;
                const int idx_midpoint_RR  = i + 2 + num_diff_ghosts_0;
                
                F_face_x[idx_face_x] += dt*(
                    a_r*(F_midpoint_x[idx_midpoint]) +
                    b_r*(F_midpoint_x[idx_midpoint_L]  + F_midpoint_x[idx_midpoint_R]) +
                    c_r*(F_midpoint_x[idx_midpoint_LL] + F_midpoint_x[idx_midpoint_RR])
                    );
            }
        }
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
            
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0 + 1; i++)
                {
                    // Compute the linear indices.
                    const int idx_face_x = i +
                        j*(interior_dim_0 + 1);
                    
                    const int idx_midpoint_LL = (i - 2 + num_diff_ghosts_0) +
                        (j + num_diff_ghosts_1)*(diff_ghostcell_dim_0 + 1);
                    
                    const int idx_midpoint_L = (i - 1 + num_diff_ghosts_0) +
                        (j + num_diff_ghosts_1)*(diff_ghostcell_dim_0 + 1);
                    
                    const int idx_midpoint = (i + num_diff_ghosts_0) +
                        (j + num_diff_ghosts_1)*(diff_ghostcell_dim_0 + 1);
                    
                    const int idx_midpoint_R = (i + 1 + num_diff_ghosts_0) +
                        (j + num_diff_ghosts_1)*(diff_ghostcell_dim_0 + 1);
                    
                    const int idx_midpoint_RR = (i + 2 + num_diff_ghosts_0) +
                        (j + num_diff_ghosts_1)*(diff_ghostcell_dim_0 + 1);
                    
                    F_face_x[idx_face_x] += dt*(
                        a_r*(F_midpoint_x[idx_midpoint]) +
                        b_r*(F_midpoint_x[idx_midpoint_L]  + F_midpoint_x[idx_midpoint_R]) +
                        c_r*(F_midpoint_x[idx_midpoint_LL] + F_midpoint_x[idx_midpoint_RR])
                        );
                }
            }
        }
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
            
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0 + 1; i++)
                    {
                        // Compute the linear indices.
                        const int idx_face_x = i +
                            j*(interior_dim_0 + 1) +
                            k*(interior_dim_0 + 1)*interior_dim_1;
                        
                        const int idx_midpoint_LL = (i - 2 + num_diff_ghosts_0) +
                            (j + num_diff_ghosts_1)*(diff_ghostcell_dim_0 + 1) +
                            (k + num_diff_ghosts_2)*(diff_ghostcell_dim_0 + 1)*
                                diff_ghostcell_dim_1;
                        
                        const int idx_midpoint_L = (i - 1 + num_diff_ghosts_0) +
                            (j + num_diff_ghosts_1)*(diff_ghostcell_dim_0 + 1) +
                            (k + num_diff_ghosts_2)*(diff_ghostcell_dim_0 + 1)*
                                diff_ghostcell_dim_1;
                        
                        const int idx_midpoint = (i + num_diff_ghosts_0) +
                            (j + num_diff_ghosts_1)*(diff_ghostcell_dim_0 + 1) +
                            (k + num_diff_ghosts_2)*(diff_ghostcell_dim_0 + 1)*
                                diff_ghostcell_dim_1;
                        
                        const int idx_midpoint_R = (i + 1 + num_diff_ghosts_0) +
                            (j + num_diff_ghosts_1)*(diff_ghostcell_dim_0 + 1) +
                            (k + num_diff_ghosts_2)*(diff_ghostcell_dim_0 + 1)*
                                diff_ghostcell_dim_1;
                        
                        const int idx_midpoint_RR = (i + 2 + num_diff_ghosts_0) +
                            (j + num_diff_ghosts_1)*(diff_ghostcell_dim_0 + 1) +
                            (k + num_diff_ghosts_2)*(diff_ghostcell_dim_0 + 1)*
                                diff_ghostcell_dim_1;
                        
                        F_face_x[idx_face_x] += dt*(
                            a_r*(F_midpoint_x[idx_midpoint]) +
                            b_r*(F_midpoint_x[idx_midpoint_L]  + F_midpoint_x[idx_midpoint_R]) +
                            c_r*(F_midpoint_x[idx_midpoint_LL] + F_midpoint_x[idx_midpoint_RR])
                            );
                    }
                }
            }
        }
    }
}


/*
 * Reconstruct the flux using flux at midpoints in y-direction.
 */
void
DiffusiveFluxReconstructorMidpointSixthOrder::reconstructFluxY(
    hier::Patch& patch,
    HAMERS_SHARED_PTR<pdat::SideData<double> >& diffusive_flux,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& diffusive_flux_midpoint,
    const double dt)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(diffusive_flux.getDepth() == d_num_eqn);
    TBOX_ASSERT(diffusive_flux_midpoint.getDepth() == d_num_eqn);
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": DiffusiveFluxReconstructorMidpointSixthOrder::"
            << "reconstructFluxY()\n"
            << "There isn't y-direction for 1D problem."
            << std::endl);
    }
#endif
    
    const double a_m =  double(75)/double(64);
    const double b_m = -double(25)/double(384);
    const double c_m =  double(3)/double(640);
    
    const double a_r = a_m + b_m + c_m;
    const double b_r = b_m + c_m;
    const double c_r = c_m;
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // diffusive ghost cells.
    hier::Box diff_ghost_box = interior_box;
    diff_ghost_box.grow(d_num_diff_ghosts);
    const hier::IntVector diff_ghostcell_dims = diff_ghost_box.numberCells();
    
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        double* F_face_y     = diffusive_flux->getPointer(1, ei);
        double* F_midpoint_y = diffusive_flux_midpoint->getPointer(1, ei);
        
        if (d_dim == tbox::Dimension(2))
        {
            /*
             * Get the dimensions and number of ghost cells.
             */
            
            const int interior_dim_0 = interior_dims[0];
            const int interior_dim_1 = interior_dims[1];
            
            const int num_diff_ghosts_0 = d_num_diff_ghosts[0];
            const int num_diff_ghosts_1 = d_num_diff_ghosts[1];
            const int diff_ghostcell_dim_0 = diff_ghostcell_dims[0];
            
            for (int j = 0; j < interior_dim_1 + 1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_face_y = i +
                        j*interior_dim_0;
                    
                    const int idx_midpoint_BB = (i + num_diff_ghosts_0) +
                        (j - 2 + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                    
                    const int idx_midpoint_B = (i + num_diff_ghosts_0) +
                        (j - 1 + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                    
                    const int idx_midpoint = (i + num_diff_ghosts_0) +
                        (j + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                    
                    const int idx_midpoint_T = (i + num_diff_ghosts_0) +
                        (j + 1 + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                    
                    const int idx_midpoint_TT = (i + num_diff_ghosts_0) +
                        (j + 2 + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                    
                    F_face_y[idx_face_y] += dt*(
                        a_r*(F_midpoint_y[idx_midpoint]) +
                        b_r*(F_midpoint_y[idx_midpoint_B]  + F_midpoint_y[idx_midpoint_T]) +
                        c_r*(F_midpoint_y[idx_midpoint_BB] + F_midpoint_y[idx_midpoint_TT])
                        );
                }
            }
        }
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
            
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1 + 1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_face_y = i +
                            j*interior_dim_0 +
                            k*interior_dim_0*(interior_dim_1 + 1);
                        
                        const int idx_midpoint_BB = (i + num_diff_ghosts_0) +
                            (j - 2 + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                            (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                (diff_ghostcell_dim_1 + 1);
                        
                        const int idx_midpoint_B = (i + num_diff_ghosts_0) +
                            (j - 1 + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                            (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                (diff_ghostcell_dim_1 + 1);
                        
                        const int idx_midpoint = (i + num_diff_ghosts_0) +
                            (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                            (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                (diff_ghostcell_dim_1 + 1);
                        
                        const int idx_midpoint_T = (i + num_diff_ghosts_0) +
                            (j + 1 + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                            (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                (diff_ghostcell_dim_1 + 1);
                        
                        const int idx_midpoint_TT = (i + num_diff_ghosts_0) +
                            (j + 2 + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                            (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                (diff_ghostcell_dim_1 + 1);
                        
                        F_face_y[idx_face_y] += dt*(
                            a_r*(F_midpoint_y[idx_midpoint]) +
                            b_r*(F_midpoint_y[idx_midpoint_B]  + F_midpoint_y[idx_midpoint_T]) +
                            c_r*(F_midpoint_y[idx_midpoint_BB] + F_midpoint_y[idx_midpoint_TT])
                            );
                    }
                }
            }
        }
    }
}


/*
 * Reconstruct the flux using flux at midpoints in z-direction.
 */
void
DiffusiveFluxReconstructorMidpointSixthOrder::reconstructFluxZ(
    hier::Patch& patch,
    HAMERS_SHARED_PTR<pdat::SideData<double> >& diffusive_flux,
    const HAMERS_SHARED_PTR<pdat::SideData<double> >& diffusive_flux_midpoint,
    const double dt)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(diffusive_flux.getDepth() == d_num_eqn);
    TBOX_ASSERT(diffusive_flux_midpoint.getDepth() == d_num_eqn);
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": DiffusiveFluxReconstructorMidpointSixthOrder::"
            << "reconstructFluxZ()\n"
            << "There isn't z-direction for 1D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": DiffusiveFluxReconstructorMidpointSixthOrder::"
            << "reconstructFluxZ()\n"
            << "There isn't z-direction for 2D problem."
            << std::endl);
    }
#endif
    
    const double a_m =  double(75)/double(64);
    const double b_m = -double(25)/double(384);
    const double c_m =  double(3)/double(640);
    
    const double a_r = a_m + b_m + c_m;
    const double b_r = b_m + c_m;
    const double c_r = c_m;
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // diffusive ghost cells.
    hier::Box diff_ghost_box = interior_box;
    diff_ghost_box.grow(d_num_diff_ghosts);
    const hier::IntVector diff_ghostcell_dims = diff_ghost_box.numberCells();
    
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        double* F_face_z     = diffusive_flux->getPointer(2, ei);
        double* F_midpoint_z = diffusive_flux_midpoint->getPointer(2, ei);
        
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
        
        for (int k = 0; k < interior_dim_2 + 1; k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_face_z = i +
                        j*interior_dim_0 +
                        k*interior_dim_0*interior_dim_1;
                    
                    const int idx_midpoint_BB = (i + num_diff_ghosts_0) +
                        (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                        (k - 2 + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                            diff_ghostcell_dim_1;
                    
                    const int idx_midpoint_B = (i + num_diff_ghosts_0) +
                        (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                        (k - 1 + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                            diff_ghostcell_dim_1;
                    
                    const int idx_midpoint = (i + num_diff_ghosts_0) +
                        (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                        (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                            diff_ghostcell_dim_1;
                    
                    const int idx_midpoint_F = (i + num_diff_ghosts_0) +
                        (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                        (k + 1 + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                            diff_ghostcell_dim_1;
                    
                    const int idx_midpoint_FF = (i + num_diff_ghosts_0) +
                        (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                        (k + 2 + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                            diff_ghostcell_dim_1;
                    
                    F_face_z[idx_face_z] += dt*(
                        a_r*(F_midpoint_z[idx_midpoint]) +
                        b_r*(F_midpoint_z[idx_midpoint_B]  + F_midpoint_z[idx_midpoint_F]) +
                        c_r*(F_midpoint_z[idx_midpoint_BB] + F_midpoint_z[idx_midpoint_FF])
                        );
                }
            }
        }
    }
}
