#include "flow/diffusive_flux_reconstructors/sixth_order/DiffusiveFluxReconstructorSixthOrder.hpp"

#include "SAMRAI/geom/CartesianPatchGeometry.h"

#include <map>

DiffusiveFluxReconstructorSixthOrder::DiffusiveFluxReconstructorSixthOrder(
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
    d_num_diff_ghosts = hier::IntVector::getOne(d_dim)*6;
    
    d_num_der_node_ghosts         = hier::IntVector::getOne(d_dim)*3;
    d_num_flux_reconstruct_ghosts = hier::IntVector::getOne(d_dim)*3;
}


/*
 * Print all characteristics of the diffusive flux reconstruction class.
 */
void
DiffusiveFluxReconstructorSixthOrder::printClassData(
    std::ostream& os) const
{
    os << "\nPrint DiffusiveFluxReconstructorSixthOrder object..."
       << std::endl;
    
    os << std::endl;
    
    os << "DiffusiveFluxReconstructorSixthOrder: this = "
       << (DiffusiveFluxReconstructorSixthOrder *)this
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
DiffusiveFluxReconstructorSixthOrder::putToRestart(
   const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const
{
    restart_db->putString("d_diffusive_flux_reconstructor", "SIXTH_ORDER");
}


/*
 * Compute the diffusive flux on a patch.
 */
void
DiffusiveFluxReconstructorSixthOrder::computeDiffusiveFluxOnPatch(
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
    
    const bool allocate_scratch_data_containers = true;
    
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
        
        diffusive_flux_utilities->registerDerivedVariablesForDiffusiveFluxes(d_num_diff_ghosts);
        
        diffusive_flux_utilities->allocateMemoryForDerivedCellData();
        
        diffusive_flux_utilities->computeDerivedCellData();
        
        /*
         * Delcare containers for computing fluxes in different directions.
         */
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > var_data_x;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > diffusivities_data_x;
        
        std::vector<std::vector<int> > var_component_idx_x;
        
        std::vector<std::vector<int> > diffusivities_component_idx_x;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > derivatives_x;
        
        std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivatives_x_computed;
        
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > diffusive_flux_node(1);
        diffusive_flux_node[0].reset(new pdat::CellData<double>(interior_box, d_num_eqn, d_num_diff_ghosts));
        
        std::vector<double*> F_node_x;
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
            patch,
            derivatives_x,
            derivatives_x_computed,
            var_data_x,
            var_component_idx_x,
            allocate_scratch_data_containers);
        
        /*
         * Compute diffusive flux in x-direction at nodes.
         */
        
        diffusive_flux_node[0]->fillAll(double(0));
        
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
                double* mu = diffusivities_data_x[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivative.
                double* dudx = derivatives_x[ei][vi]->getPointer(0);
                
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
            patch,
            diffusive_flux,
            diffusive_flux_node[0],
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
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > var_data_x;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > var_data_y;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > diffusivities_data_x;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > diffusivities_data_y;
        
        std::vector<std::vector<int> > var_component_idx_x;
        std::vector<std::vector<int> > var_component_idx_y;
        
        std::vector<std::vector<int> > diffusivities_component_idx_x;
        std::vector<std::vector<int> > diffusivities_component_idx_y;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > derivatives_x;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > derivatives_y;
        
        std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivatives_x_computed;
        std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivatives_y_computed;
        
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > diffusive_flux_node(2);
        diffusive_flux_node[0].reset(new pdat::CellData<double>(interior_box, d_num_eqn, d_num_diff_ghosts));
        diffusive_flux_node[1].reset(new pdat::CellData<double>(interior_box, d_num_eqn, d_num_diff_ghosts));
        
        std::vector<double*> F_node_x;
        std::vector<double*> F_node_y;
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
            patch,
            derivatives_x,
            derivatives_x_computed,
            var_data_x,
            var_component_idx_x,
            allocate_scratch_data_containers);
        
        /*
         * Compute the derivatives in y-direction for diffusive flux in x-direction.
         */
        
        computeFirstDerivativesInY(
            patch,
            derivatives_y,
            derivatives_y_computed,
            var_data_y,
            var_component_idx_y,
            allocate_scratch_data_containers);
        
        /*
         * Compute diffusive flux in x-direction at nodes.
         */
        
        diffusive_flux_node[0]->fillAll(double(0));
        
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
                double* mu = diffusivities_data_x[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivative.
                double* dudx = derivatives_x[ei][vi]->getPointer(0);
                
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
                double* mu = diffusivities_data_y[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivative.
                double* dudy = derivatives_y[ei][vi]->getPointer(0);
                
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
            patch,
            diffusive_flux,
            diffusive_flux_node[0],
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
            patch,
            derivatives_x,
            derivatives_x_computed,
            var_data_x,
            var_component_idx_x,
            allocate_scratch_data_containers);
        
        /*
         * Compute the derivatives in y-direction for diffusive flux in y-direction.
         */
        
        computeFirstDerivativesInY(
            patch,
            derivatives_y,
            derivatives_y_computed,
            var_data_y,
            var_component_idx_y,
            allocate_scratch_data_containers);
        
        /*
         * Compute diffusive flux in y-direction at nodes.
         */
        
        diffusive_flux_node[1]->fillAll(double(0));
        
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
                double* mu = diffusivities_data_x[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivative.
                double* dudx = derivatives_x[ei][vi]->getPointer(0);
                
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
                double* mu = diffusivities_data_y[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivative.
                double* dudy = derivatives_y[ei][vi]->getPointer(0);
                
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
                        
                        F_node_y[ei][idx] += mu[idx_diffusivity]*dudy[idx];
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
            diffusive_flux_node[1],
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
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > var_data_x;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > var_data_y;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > var_data_z;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > diffusivities_data_x;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > diffusivities_data_y;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > diffusivities_data_z;
        
        std::vector<std::vector<int> > var_component_idx_x;
        std::vector<std::vector<int> > var_component_idx_y;
        std::vector<std::vector<int> > var_component_idx_z;
        
        std::vector<std::vector<int> > diffusivities_component_idx_x;
        std::vector<std::vector<int> > diffusivities_component_idx_y;
        std::vector<std::vector<int> > diffusivities_component_idx_z;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > derivatives_x;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > derivatives_y;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > derivatives_z;
        
        std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivatives_x_computed;
        std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivatives_y_computed;
        std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivatives_z_computed;
        
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > diffusive_flux_node(3);
        diffusive_flux_node[0].reset(new pdat::CellData<double>(interior_box, d_num_eqn, d_num_diff_ghosts));
        diffusive_flux_node[1].reset(new pdat::CellData<double>(interior_box, d_num_eqn, d_num_diff_ghosts));
        diffusive_flux_node[2].reset(new pdat::CellData<double>(interior_box, d_num_eqn, d_num_diff_ghosts));
        
        std::vector<double*> F_node_x;
        std::vector<double*> F_node_y;
        std::vector<double*> F_node_z;
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
            patch,
            derivatives_x,
            derivatives_x_computed,
            var_data_x,
            var_component_idx_x,
            allocate_scratch_data_containers);
        
        /*
         * Compute the derivatives in y-direction for diffusive flux in x-direction.
         */
        
        computeFirstDerivativesInY(
            patch,
            derivatives_y,
            derivatives_y_computed,
            var_data_y,
            var_component_idx_y,
            allocate_scratch_data_containers);
        
        /*
         * Compute the derivatives in z-direction for diffusive flux in x-direction.
         */
        
        computeFirstDerivativesInZ(
            patch,
            derivatives_z,
            derivatives_z_computed,
            var_data_z,
            var_component_idx_z,
            allocate_scratch_data_containers);
        
        /*
         * Compute diffusive flux in x-direction at nodes.
         */
        
        diffusive_flux_node[0]->fillAll(double(0));
        
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
                double* mu = diffusivities_data_x[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivative.
                double* dudx = derivatives_x[ei][vi]->getPointer(0);
                
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
                double* mu = diffusivities_data_y[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivative.
                double* dudy = derivatives_y[ei][vi]->getPointer(0);
                
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
                double* mu = diffusivities_data_z[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivative.
                double* dudz = derivatives_z[ei][vi]->getPointer(0);
                
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
            patch,
            diffusive_flux,
            diffusive_flux_node[0],
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
            patch,
            derivatives_x,
            derivatives_x_computed,
            var_data_x,
            var_component_idx_x,
            allocate_scratch_data_containers);
        
        /*
         * Compute the derivatives in y-direction for diffusive flux in y-direction.
         */
        
        computeFirstDerivativesInY(
            patch,
            derivatives_y,
            derivatives_y_computed,
            var_data_y,
            var_component_idx_y,
            allocate_scratch_data_containers);
        
        /*
         * Compute the derivatives in z-direction for diffusive flux in y-direction.
         */
        
        computeFirstDerivativesInZ(
            patch,
            derivatives_z,
            derivatives_z_computed,
            var_data_z,
            var_component_idx_z,
            allocate_scratch_data_containers);
        
        /*
         * Compute diffusive flux in y-direction at nodes.
         */
        
        diffusive_flux_node[1]->fillAll(double(0));
        
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
                double* mu = diffusivities_data_x[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivative.
                double* dudx = derivatives_x[ei][vi]->getPointer(0);
                
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
                double* mu = diffusivities_data_y[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivative.
                double* dudy = derivatives_y[ei][vi]->getPointer(0);
                
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
                double* mu = diffusivities_data_z[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivative.
                double* dudz = derivatives_z[ei][vi]->getPointer(0);
                
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
            patch,
            diffusive_flux,
            diffusive_flux_node[1],
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
            patch,
            derivatives_x,
            derivatives_x_computed,
            var_data_x,
            var_component_idx_x,
            allocate_scratch_data_containers);
        
        /*
         * Compute the derivatives in y-direction for diffusive flux in z-direction.
         */
        
        computeFirstDerivativesInY(
            patch,
            derivatives_y,
            derivatives_y_computed,
            var_data_y,
            var_component_idx_y,
            allocate_scratch_data_containers);
        
        /*
         * Compute the derivatives in z-direction for diffusive flux in z-direction.
         */
        
        computeFirstDerivativesInZ(
            patch,
            derivatives_z,
            derivatives_z_computed,
            var_data_z,
            var_component_idx_z,
            allocate_scratch_data_containers);
        
        /*
         * Compute diffusive flux in z-direction at nodes.
         */
        
        diffusive_flux_node[2]->fillAll(double(0));
        
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
                double* mu = diffusivities_data_x[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivative.
                double* dudx = derivatives_x[ei][vi]->getPointer(0);
                
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
                double* mu = diffusivities_data_y[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivative.
                double* dudy = derivatives_y[ei][vi]->getPointer(0);
                
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
                double* mu = diffusivities_data_z[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivative.
                double* dudz = derivatives_z[ei][vi]->getPointer(0);
                
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
            patch,
            diffusive_flux,
            diffusive_flux_node[2],
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
DiffusiveFluxReconstructorSixthOrder::computeFirstDerivativesInX(
    hier::Patch& patch,
    std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& derivatives_x,
    std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > >& derivatives_x_computed,
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
    
    derivatives_x.resize(d_num_eqn);
    
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
    
    const double dx_0_inv = double(1)/dx[0];
    
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        derivatives_x[ei].reserve(static_cast<int>(data_x[ei].size()));
        
        for (int vi = 0; vi < static_cast<int>(data_x[ei].size()); vi++)
        {
            // Get the index of variable for derivative.
            const int u_idx = data_component_idx_x[ei][vi];
            
            if (derivatives_x_computed.find(data_x[ei][vi]->getPointer(u_idx))
                == derivatives_x_computed.end())
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
                 * Get the ghost cell widths and ghost box dimensions of the variables.
                 */
                
                const hier::IntVector& num_ghosts_derivative_node = d_num_diff_ghosts;
                
                const hier::IntVector& ghostcell_dims_derivative_node = diff_ghostcell_dims;
                
                const hier::IntVector& num_ghosts_data_node =
                    data_x[ei][vi]->getGhostCellWidth();
                
                const hier::IntVector& ghostcell_dims_data_node =
                    data_x[ei][vi]->getGhostBox().numberCells();
                
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
                
                std::pair<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivative_pair(
                    u,
                    d_scratch_derivatives_node[d_num_scratch_derivatives_node_used - 1]);
                
                derivatives_x_computed.insert(derivative_pair);
            }
            
            derivatives_x[ei].push_back(
                derivatives_x_computed.find(data_x[ei][vi]->getPointer(u_idx))->
                    second);
        }
    }
}


/*
 * Compute the first derivatives in the y-direction.
 */
void
DiffusiveFluxReconstructorSixthOrder::computeFirstDerivativesInY(
    hier::Patch& patch,
    std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& derivatives_y,
    std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > >& derivatives_y_computed,
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
            << ": DiffusiveFluxReconstructorSixthOrder::"
            << "computeFirstDerivativesInY()\n"
            << "There isn't y-direction for one-dimensional problem."
            << std::endl);
    }
#endif
    
    derivatives_y.resize(d_num_eqn);
    
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
    
    const double dx_1_inv = double(1)/dx[1];
    
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        derivatives_y[ei].reserve(static_cast<int>(data_y[ei].size()));
        
        for (int vi = 0; vi < static_cast<int>(data_y[ei].size()); vi++)
        {
            // Get the index of variable for derivative.
            const int u_idx = data_component_idx_y[ei][vi];
            
            if (derivatives_y_computed.find(data_y[ei][vi]->getPointer(u_idx))
                == derivatives_y_computed.end())
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
                 * Get the ghost cell widths and ghost box dimensions of the variables.
                 */
                
                const hier::IntVector& num_ghosts_derivative_node = d_num_diff_ghosts;
                
                const hier::IntVector& ghostcell_dims_derivative_node = diff_ghostcell_dims;
                
                const hier::IntVector& num_ghosts_data_node =
                    data_y[ei][vi]->getGhostCellWidth();
                
                const hier::IntVector& ghostcell_dims_data_node =
                    data_y[ei][vi]->getGhostBox().numberCells();
                
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
                
                std::pair<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivative_pair(
                    u,
                    d_scratch_derivatives_node[d_num_scratch_derivatives_node_used - 1]);
                
                derivatives_y_computed.insert(derivative_pair);
            }
            
            derivatives_y[ei].push_back(
                derivatives_y_computed.find(data_y[ei][vi]->getPointer(u_idx))->
                    second);
        }
    }
}


/*
 * Compute the first derivatives in the z-direction.
 */
void
DiffusiveFluxReconstructorSixthOrder::computeFirstDerivativesInZ(
    hier::Patch& patch,
    std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& derivatives_z,
    std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > >& derivatives_z_computed,
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
            << ": DiffusiveFluxReconstructorSixthOrder::"
            << "computeFirstDerivativesInZ()\n"
            << "There isn't z-direction for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": DiffusiveFluxReconstructorSixthOrder::"
            << "computeFirstDerivativesInZ()\n"
            << "There isn't z-direction for two-dimensional problem."
            << std::endl);
    }
#endif
    
    derivatives_z.resize(d_num_eqn);
    
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
    
    const double dx_2_inv = double(1)/dx[2];
    
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        derivatives_z[ei].reserve(static_cast<int>(data_z[ei].size()));
        
        for (int vi = 0; vi < static_cast<int>(data_z[ei].size()); vi++)
        {
            // Get the index of variable for derivative.
            const int u_idx = data_component_idx_z[ei][vi];
            
            if (derivatives_z_computed.find(data_z[ei][vi]->getPointer(u_idx))
                == derivatives_z_computed.end())
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
                 * Get the ghost cell widths and ghost box dimensions of the variables.
                 */
                
                const hier::IntVector& num_ghosts_derivative_node = d_num_diff_ghosts;
                
                const hier::IntVector& ghostcell_dims_derivative_node = diff_ghostcell_dims;
                
                const hier::IntVector& num_ghosts_data_node =
                    data_z[ei][vi]->getGhostCellWidth();
                
                const hier::IntVector& ghostcell_dims_data_node =
                    data_z[ei][vi]->getGhostBox().numberCells();
                
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
                
                std::pair<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivative_pair(
                    u,
                    d_scratch_derivatives_node[d_num_scratch_derivatives_node_used - 1]);
                
                derivatives_z_computed.insert(derivative_pair);
            }
            
            derivatives_z[ei].push_back(
                derivatives_z_computed.find(data_z[ei][vi]->getPointer(u_idx))->
                    second);
        }
    }
}


/*
 * Reconstruct the flux using flux at nodes in x-direction.
 */
void
DiffusiveFluxReconstructorSixthOrder::reconstructFluxX(
    hier::Patch& patch,
    HAMERS_SHARED_PTR<pdat::SideData<double> >& diffusive_flux,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& diffusive_flux_node,
    const double dt)
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
        double* F_face_x = diffusive_flux->getPointer(0, ei);
        double* F_node_x = diffusive_flux_node->getPointer(ei);
        
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
DiffusiveFluxReconstructorSixthOrder::reconstructFluxY(
    hier::Patch& patch,
    HAMERS_SHARED_PTR<pdat::SideData<double> >& diffusive_flux,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& diffusive_flux_node,
    const double dt)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(diffusive_flux->getDepth() == d_num_eqn);
    TBOX_ASSERT(diffusive_flux_node->getDepth() == d_num_eqn);
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": DiffusiveFluxReconstructorSixthOrder::"
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
        double* F_face_y = diffusive_flux->getPointer(1, ei);
        double* F_node_y = diffusive_flux_node->getPointer(ei);
        
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
DiffusiveFluxReconstructorSixthOrder::reconstructFluxZ(
    hier::Patch& patch,
    HAMERS_SHARED_PTR<pdat::SideData<double> >& diffusive_flux,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& diffusive_flux_node,
    const double dt)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(diffusive_flux->getDepth() == d_num_eqn);
    TBOX_ASSERT(diffusive_flux_node->getDepth() == d_num_eqn);
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": DiffusiveFluxReconstructorSixthOrder::"
            << "reconstructFluxZ()\n"
            << "There isn't z-direction for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": DiffusiveFluxReconstructorSixthOrder::"
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
        double* F_face_z = diffusive_flux->getPointer(2, ei);
        double* F_node_z = diffusive_flux_node->getPointer(ei);
        
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


///


/*
 * Kernel to compute the first derivatives in the x-direction.
 */
void
DiffusiveFluxReconstructorSixthOrder::computeFirstDerivativesInX(
    double* dudx,
    const double* const u,
    const hier::IntVector& num_ghosts_derivative_node,
    const hier::IntVector& num_ghosts_data_node,
    const hier::IntVector& ghostcell_dims_derivative_node,
    const hier::IntVector& ghostcell_dims_data_node,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims,
    const double& dx_0_inv)
{
    const double a_n =  double(3)/double(4);
    const double b_n = -double(3)/double(20);
    const double c_n =  double(1)/double(60);
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and number of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int num_ghosts_0_derivative_node = num_ghosts_derivative_node[0];
        
        const int num_ghosts_0_data = num_ghosts_data_node[0];
        
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx = i + num_ghosts_0_derivative_node;
            
            const int idx_data_LLL = i - 3 + num_ghosts_0_data;
            const int idx_data_LL  = i - 2 + num_ghosts_0_data;
            const int idx_data_L   = i - 1 + num_ghosts_0_data;
            const int idx_data_R   = i + 1 + num_ghosts_0_data;
            const int idx_data_RR  = i + 2 + num_ghosts_0_data;
            const int idx_data_RRR = i + 3 + num_ghosts_0_data;
            
            dudx[idx] = (a_n*(u[idx_data_R]   - u[idx_data_L]) +
                         b_n*(u[idx_data_RR]  - u[idx_data_LL]) +
                         c_n*(u[idx_data_RRR] - u[idx_data_LLL])
                        )*dx_0_inv;
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and number of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        
        const int num_ghosts_0_derivative_node = num_ghosts_derivative_node[0];
        const int num_ghosts_1_derivative_node = num_ghosts_derivative_node[1];
        const int ghostcell_dim_0_derivative_node = ghostcell_dims_derivative_node[0];
        
        const int num_ghosts_0_data = num_ghosts_data_node[0];
        const int num_ghosts_1_data = num_ghosts_data_node[1];
        const int ghostcell_dim_0_data = ghostcell_dims_data_node[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx = (i + num_ghosts_0_derivative_node) +
                    (j + num_ghosts_1_derivative_node)*ghostcell_dim_0_derivative_node;
                
                const int idx_data_LLL = (i - 3 + num_ghosts_0_data) +
                    (j + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                const int idx_data_LL = (i - 2 + num_ghosts_0_data) +
                    (j + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                const int idx_data_L = (i - 1 + num_ghosts_0_data) +
                    (j + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                const int idx_data_R = (i + 1 + num_ghosts_0_data) +
                    (j + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                const int idx_data_RR = (i + 2 + num_ghosts_0_data) +
                    (j + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                const int idx_data_RRR = (i + 3 + num_ghosts_0_data) +
                    (j + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                dudx[idx] = (a_n*(u[idx_data_R]   - u[idx_data_L]) +
                             b_n*(u[idx_data_RR]  - u[idx_data_LL]) +
                             c_n*(u[idx_data_RRR] - u[idx_data_LLL])
                            )*dx_0_inv;
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and number of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_lo_2 = domain_lo[2];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        const int domain_dim_2 = domain_dims[2];
        
        const int num_ghosts_0_derivative_node = num_ghosts_derivative_node[0];
        const int num_ghosts_1_derivative_node = num_ghosts_derivative_node[1];
        const int num_ghosts_2_derivative_node = num_ghosts_derivative_node[2];
        const int ghostcell_dim_0_derivative_node = ghostcell_dims_derivative_node[0];
        const int ghostcell_dim_1_derivative_node = ghostcell_dims_derivative_node[1];
        
        const int num_ghosts_0_data = num_ghosts_data_node[0];
        const int num_ghosts_1_data = num_ghosts_data_node[1];
        const int num_ghosts_2_data = num_ghosts_data_node[2];
        const int ghostcell_dim_0_data = ghostcell_dims_data_node[0];
        const int ghostcell_dim_1_data = ghostcell_dims_data_node[1];
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx = (i + num_ghosts_0_derivative_node) +
                        (j + num_ghosts_1_derivative_node)*ghostcell_dim_0_derivative_node +
                        (k + num_ghosts_2_derivative_node)*ghostcell_dim_0_derivative_node*
                            ghostcell_dim_1_derivative_node;
                    
                    const int idx_data_LLL = (i - 3 + num_ghosts_0_data) +
                        (j + num_ghosts_1_data)*ghostcell_dim_0_data +
                        (k + num_ghosts_2_data)*ghostcell_dim_0_data*
                            ghostcell_dim_1_data;
                    
                    const int idx_data_LL = (i - 2 + num_ghosts_0_data) +
                        (j + num_ghosts_1_data)*ghostcell_dim_0_data +
                        (k + num_ghosts_2_data)*ghostcell_dim_0_data*
                            ghostcell_dim_1_data;
                    
                    const int idx_data_L = (i - 1 + num_ghosts_0_data) +
                        (j + num_ghosts_1_data)*ghostcell_dim_0_data +
                        (k + num_ghosts_2_data)*ghostcell_dim_0_data*
                            ghostcell_dim_1_data;
                    
                    const int idx_data_R = (i + 1 + num_ghosts_0_data) +
                        (j + num_ghosts_1_data)*ghostcell_dim_0_data +
                        (k + num_ghosts_2_data)*ghostcell_dim_0_data*
                            ghostcell_dim_1_data;
                    
                    const int idx_data_RR = (i + 2 + num_ghosts_0_data) +
                        (j + num_ghosts_1_data)*ghostcell_dim_0_data +
                        (k + num_ghosts_2_data)*ghostcell_dim_0_data*
                            ghostcell_dim_1_data;
                    
                    const int idx_data_RRR = (i + 3 + num_ghosts_0_data) +
                        (j + num_ghosts_1_data)*ghostcell_dim_0_data +
                        (k + num_ghosts_2_data)*ghostcell_dim_0_data*
                            ghostcell_dim_1_data;
                    
                    dudx[idx] = (a_n*(u[idx_data_R]   - u[idx_data_L]) +
                                 b_n*(u[idx_data_RR]  - u[idx_data_LL]) +
                                 c_n*(u[idx_data_RRR] - u[idx_data_LLL])
                                )*dx_0_inv;
                }
            }
        }
    }
}


/*
 * Kernel to compute the first derivatives in the y-direction.
 */
void
DiffusiveFluxReconstructorSixthOrder::computeFirstDerivativesInY(
    double* dudy,
    const double* const u,
    const hier::IntVector& num_ghosts_derivative_node,
    const hier::IntVector& num_ghosts_data_node,
    const hier::IntVector& ghostcell_dims_derivative_node,
    const hier::IntVector& ghostcell_dims_data_node,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims,
    const double& dx_1_inv)
{
    const double a_n =  double(3)/double(4);
    const double b_n = -double(3)/double(20);
    const double c_n =  double(1)/double(60);
    
    if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and number of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        
        const int num_ghosts_0_derivative_node = num_ghosts_derivative_node[0];
        const int num_ghosts_1_derivative_node = num_ghosts_derivative_node[1];
        const int ghostcell_dim_0_derivative_node = ghostcell_dims_derivative_node[0];
        
        const int num_ghosts_0_data = num_ghosts_data_node[0];
        const int num_ghosts_1_data = num_ghosts_data_node[1];
        const int ghostcell_dim_0_data = ghostcell_dims_data_node[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx = (i + num_ghosts_0_derivative_node) +
                    (j + num_ghosts_1_derivative_node)*ghostcell_dim_0_derivative_node;
                
                const int idx_data_BBB = (i + num_ghosts_0_data) +
                    (j - 3 + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                const int idx_data_BB = (i + num_ghosts_0_data) +
                    (j - 2 + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                const int idx_data_B = (i + num_ghosts_0_data) +
                    (j - 1 + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                const int idx_data_T = (i + num_ghosts_0_data) +
                    (j + 1 + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                const int idx_data_TT = (i + num_ghosts_0_data) +
                    (j + 2 + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                const int idx_data_TTT = (i + num_ghosts_0_data) +
                    (j + 3 + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                dudy[idx] = (a_n*(u[idx_data_T]   - u[idx_data_B]) +
                             b_n*(u[idx_data_TT]  - u[idx_data_BB]) +
                             c_n*(u[idx_data_TTT] - u[idx_data_BBB])
                            )*dx_1_inv;
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and number of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_lo_2 = domain_lo[2];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        const int domain_dim_2 = domain_dims[2];
        
        const int num_ghosts_0_derivative_node = num_ghosts_derivative_node[0];
        const int num_ghosts_1_derivative_node = num_ghosts_derivative_node[1];
        const int num_ghosts_2_derivative_node = num_ghosts_derivative_node[2];
        const int ghostcell_dim_0_derivative_node = ghostcell_dims_derivative_node[0];
        const int ghostcell_dim_1_derivative_node = ghostcell_dims_derivative_node[1];
        
        const int num_ghosts_0_data = num_ghosts_data_node[0];
        const int num_ghosts_1_data = num_ghosts_data_node[1];
        const int num_ghosts_2_data = num_ghosts_data_node[2];
        const int ghostcell_dim_0_data = ghostcell_dims_data_node[0];
        const int ghostcell_dim_1_data = ghostcell_dims_data_node[1];
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx = (i + num_ghosts_0_derivative_node) +
                        (j + num_ghosts_1_derivative_node)*ghostcell_dim_0_derivative_node +
                        (k + num_ghosts_2_derivative_node)*ghostcell_dim_0_derivative_node*
                            ghostcell_dim_1_derivative_node;
                    
                    const int idx_data_BBB = (i + num_ghosts_0_data) +
                        (j - 3 + num_ghosts_1_data)*ghostcell_dim_0_data +
                        (k + num_ghosts_2_data)*ghostcell_dim_0_data*
                            ghostcell_dim_1_data;
                    
                    const int idx_data_BB = (i + num_ghosts_0_data) +
                        (j - 2 + num_ghosts_1_data)*ghostcell_dim_0_data +
                        (k + num_ghosts_2_data)*ghostcell_dim_0_data*
                            ghostcell_dim_1_data;
                    
                    const int idx_data_B = (i + num_ghosts_0_data) +
                        (j - 1 + num_ghosts_1_data)*ghostcell_dim_0_data +
                        (k + num_ghosts_2_data)*ghostcell_dim_0_data*
                            ghostcell_dim_1_data;
                    
                    const int idx_data_T = (i + num_ghosts_0_data) +
                        (j + 1 + num_ghosts_1_data)*ghostcell_dim_0_data +
                        (k + num_ghosts_2_data)*ghostcell_dim_0_data*
                            ghostcell_dim_1_data;
                    
                    const int idx_data_TT = (i + num_ghosts_0_data) +
                        (j + 2 + num_ghosts_1_data)*ghostcell_dim_0_data +
                        (k + num_ghosts_2_data)*ghostcell_dim_0_data*
                            ghostcell_dim_1_data;
                    
                    const int idx_data_TTT = (i + num_ghosts_0_data) +
                        (j + 3 + num_ghosts_1_data)*ghostcell_dim_0_data +
                        (k + num_ghosts_2_data)*ghostcell_dim_0_data*
                            ghostcell_dim_1_data;
                    
                    dudy[idx] = (a_n*(u[idx_data_T]   - u[idx_data_B]) +
                                 b_n*(u[idx_data_TT]  - u[idx_data_BB]) +
                                 c_n*(u[idx_data_TTT] - u[idx_data_BBB])
                                )*dx_1_inv;
                }
            }
        }
    }
}


/*
 * Kernel to compute the first derivatives in the z-direction.
 */
void
DiffusiveFluxReconstructorSixthOrder::computeFirstDerivativesInZ(
    double* dudz,
    const double* const u,
    const hier::IntVector& num_ghosts_derivative_node,
    const hier::IntVector& num_ghosts_data_node,
    const hier::IntVector& ghostcell_dims_derivative_node,
    const hier::IntVector& ghostcell_dims_data_node,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims,
    const double& dx_2_inv)
{
    const double a_n =  double(3)/double(4);
    const double b_n = -double(3)/double(20);
    const double c_n =  double(1)/double(60);
    
    /*
     * Get the local lower index, numbers of cells in each dimension and number of ghost cells.
     */
    
    const int domain_lo_0 = domain_lo[0];
    const int domain_lo_1 = domain_lo[1];
    const int domain_lo_2 = domain_lo[2];
    const int domain_dim_0 = domain_dims[0];
    const int domain_dim_1 = domain_dims[1];
    const int domain_dim_2 = domain_dims[2];
    
    const int num_ghosts_0_derivative_node = num_ghosts_derivative_node[0];
    const int num_ghosts_1_derivative_node = num_ghosts_derivative_node[1];
    const int num_ghosts_2_derivative_node = num_ghosts_derivative_node[2];
    const int ghostcell_dim_0_derivative_node = ghostcell_dims_derivative_node[0];
    const int ghostcell_dim_1_derivative_node = ghostcell_dims_derivative_node[1];
    
    const int num_ghosts_0_data = num_ghosts_data_node[0];
    const int num_ghosts_1_data = num_ghosts_data_node[1];
    const int num_ghosts_2_data = num_ghosts_data_node[2];
    const int ghostcell_dim_0_data = ghostcell_dims_data_node[0];
    const int ghostcell_dim_1_data = ghostcell_dims_data_node[1];
    
    for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
    {
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx = (i + num_ghosts_0_derivative_node) +
                    (j + num_ghosts_1_derivative_node)*ghostcell_dim_0_derivative_node +
                    (k + num_ghosts_2_derivative_node)*ghostcell_dim_0_derivative_node*
                        ghostcell_dim_1_derivative_node;
                
                const int idx_data_BBB = (i + num_ghosts_0_data) +
                    (j + num_ghosts_1_data)*ghostcell_dim_0_data +
                    (k - 3 + num_ghosts_2_data)*ghostcell_dim_0_data*
                        ghostcell_dim_1_data;
                
                const int idx_data_BB = (i + num_ghosts_0_data) +
                    (j + num_ghosts_1_data)*ghostcell_dim_0_data +
                    (k - 2 + num_ghosts_2_data)*ghostcell_dim_0_data*
                        ghostcell_dim_1_data;
                
                const int idx_data_B = (i + num_ghosts_0_data) +
                    (j + num_ghosts_1_data)*ghostcell_dim_0_data +
                    (k - 1 + num_ghosts_2_data)*ghostcell_dim_0_data*
                        ghostcell_dim_1_data;
                
                const int idx_data_F = (i + num_ghosts_0_data) +
                    (j + num_ghosts_1_data)*ghostcell_dim_0_data +
                    (k + 1 + num_ghosts_2_data)*ghostcell_dim_0_data*
                        ghostcell_dim_1_data;
                
                const int idx_data_FF = (i + num_ghosts_0_data) +
                    (j + num_ghosts_1_data)*ghostcell_dim_0_data +
                    (k + 2 + num_ghosts_2_data)*ghostcell_dim_0_data*
                        ghostcell_dim_1_data;
                
                const int idx_data_FFF = (i + num_ghosts_0_data) +
                    (j + num_ghosts_1_data)*ghostcell_dim_0_data +
                    (k + 3 + num_ghosts_2_data)*ghostcell_dim_0_data*
                        ghostcell_dim_1_data;
                
                dudz[idx] = (a_n*(u[idx_data_F]   - u[idx_data_B]) +
                             b_n*(u[idx_data_FF]  - u[idx_data_BB]) +
                             c_n*(u[idx_data_FFF] - u[idx_data_BBB])
                            )*dx_2_inv;
            }
        }
    }
}


/*
 * Kernel to reconstruct the flux using flux at nodes in x-direction.
 */
void
DiffusiveFluxReconstructorSixthOrder::reconstructFluxX(
    double* F_face_x,
    const double* const F_node_x,
    const hier::IntVector& num_ghosts_flux_node,
    const hier::IntVector& ghostcell_dims_flux_node,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims,
    const hier::IntVector& interior_dims,
    const double& dt)
{
    const double a_n =  double(3)/double(4);
    const double b_n = -double(3)/double(20);
    const double c_n =  double(1)/double(60);
    
    const double a_r = a_n + b_n + c_n;
    const double b_r = b_n + c_n;
    const double c_r = c_n;
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and number of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int num_ghosts_0_flux_node = num_ghosts_flux_node[0];
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_face_x = i;
            
            const int idx_node_LLL = i - 3 + num_ghosts_0_flux_node;
            const int idx_node_LL  = i - 2 + num_ghosts_0_flux_node;
            const int idx_node_L   = i - 1 + num_ghosts_0_flux_node;
            const int idx_node_R   = i + 0 + num_ghosts_0_flux_node;
            const int idx_node_RR  = i + 1 + num_ghosts_0_flux_node;
            const int idx_node_RRR = i + 2 + num_ghosts_0_flux_node;
            
            F_face_x[idx_face_x] += dt*(
                a_r*(F_node_x[idx_node_L]   + F_node_x[idx_node_R]) +
                b_r*(F_node_x[idx_node_LL]  + F_node_x[idx_node_RR]) +
                c_r*(F_node_x[idx_node_LLL] + F_node_x[idx_node_RRR])
                );
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and number of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        
        const int interior_dim_0 = interior_dims[0];
        
        const int num_ghosts_0_flux_node = num_ghosts_flux_node[0];
        const int num_ghosts_1_flux_node = num_ghosts_flux_node[1];
        const int ghostcell_dim_0_flux_node = ghostcell_dims_flux_node[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_face_x = i +
                    j*(interior_dim_0 + 1);
                
                const int idx_node_LLL = (i - 3 + num_ghosts_0_flux_node) +
                    (j + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node;
                
                const int idx_node_LL = (i - 2 + num_ghosts_0_flux_node) +
                    (j + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node;
                
                const int idx_node_L = (i - 1 + num_ghosts_0_flux_node) +
                    (j + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node;
                
                const int idx_node_R = (i + num_ghosts_0_flux_node) +
                    (j + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node;
                
                const int idx_node_RR = (i + 1 + num_ghosts_0_flux_node) +
                    (j + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node;
                
                const int idx_node_RRR = (i + 2 + num_ghosts_0_flux_node) +
                    (j + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node;
                
                F_face_x[idx_face_x] += dt*(
                    a_r*(F_node_x[idx_node_L]   + F_node_x[idx_node_R]) +
                    b_r*(F_node_x[idx_node_LL]  + F_node_x[idx_node_RR]) +
                    c_r*(F_node_x[idx_node_LLL] + F_node_x[idx_node_RRR])
                    );
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and number of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_lo_2 = domain_lo[2];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        const int domain_dim_2 = domain_dims[2];
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        const int num_ghosts_0_flux_node = num_ghosts_flux_node[0];
        const int num_ghosts_1_flux_node = num_ghosts_flux_node[1];
        const int num_ghosts_2_flux_node = num_ghosts_flux_node[2];
        const int ghostcell_dim_0_flux_node = ghostcell_dims_flux_node[0];
        const int ghostcell_dim_1_flux_node = ghostcell_dims_flux_node[1];
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_face_x = i +
                        j*(interior_dim_0 + 1) +
                        k*(interior_dim_0 + 1)*interior_dim_1;
                    
                    const int idx_node_LLL = (i - 3 + num_ghosts_0_flux_node) +
                        (j + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node +
                        (k + num_ghosts_2_flux_node)*ghostcell_dim_0_flux_node*
                            ghostcell_dim_1_flux_node;
                    
                    const int idx_node_LL = (i - 2 + num_ghosts_0_flux_node) +
                        (j + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node +
                        (k + num_ghosts_2_flux_node)*ghostcell_dim_0_flux_node*
                            ghostcell_dim_1_flux_node;
                    
                    const int idx_node_L = (i - 1 + num_ghosts_0_flux_node) +
                        (j + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node +
                        (k + num_ghosts_2_flux_node)*ghostcell_dim_0_flux_node*
                            ghostcell_dim_1_flux_node;
                    
                    const int idx_node_R = (i + num_ghosts_0_flux_node) +
                        (j + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node +
                        (k + num_ghosts_2_flux_node)*ghostcell_dim_0_flux_node*
                            ghostcell_dim_1_flux_node;
                    
                    const int idx_node_RR = (i + 1 + num_ghosts_0_flux_node) +
                        (j + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node +
                        (k + num_ghosts_2_flux_node)*ghostcell_dim_0_flux_node*
                            ghostcell_dim_1_flux_node;
                    
                    const int idx_node_RRR = (i + 2 + num_ghosts_0_flux_node) +
                        (j + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node +
                        (k + num_ghosts_2_flux_node)*ghostcell_dim_0_flux_node*
                            ghostcell_dim_1_flux_node;
                    
                    F_face_x[idx_face_x] += dt*(
                        a_r*(F_node_x[idx_node_L]   + F_node_x[idx_node_R]) +
                        b_r*(F_node_x[idx_node_LL]  + F_node_x[idx_node_RR]) +
                        c_r*(F_node_x[idx_node_LLL] + F_node_x[idx_node_RRR])
                        );
                }
            }
        }
    }
}


/*
 * Kernel to reconstruct the flux using flux at nodes in y-direction.
 */
void
DiffusiveFluxReconstructorSixthOrder::reconstructFluxY(
    double* F_face_y,
    const double* const F_node_y,
    const hier::IntVector& num_ghosts_flux_node,
    const hier::IntVector& ghostcell_dims_flux_node,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims,
    const hier::IntVector& interior_dims,
    const double& dt)
{
    const double a_n =  double(3)/double(4);
    const double b_n = -double(3)/double(20);
    const double c_n =  double(1)/double(60);
    
    const double a_r = a_n + b_n + c_n;
    const double b_r = b_n + c_n;
    const double c_r = c_n;
    
    if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and number of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        
        const int interior_dim_0 = interior_dims[0];
        
        const int num_ghosts_0_flux_node = num_ghosts_flux_node[0];
        const int num_ghosts_1_flux_node = num_ghosts_flux_node[1];
        const int ghostcell_dim_0_flux_node = ghostcell_dims_flux_node[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_face_y = i +
                    j*interior_dim_0;
                
                const int idx_node_BBB = (i + num_ghosts_0_flux_node) +
                    (j - 3 + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node;
                
                const int idx_node_BB = (i + num_ghosts_0_flux_node) +
                    (j - 2 + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node;
                
                const int idx_node_B = (i + num_ghosts_0_flux_node) +
                    (j - 1 + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node;
                
                const int idx_node_T = (i + num_ghosts_0_flux_node) +
                    (j + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node;
                
                const int idx_node_TT = (i + num_ghosts_0_flux_node) +
                    (j + 1 + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node;
                
                const int idx_node_TTT = (i + num_ghosts_0_flux_node) +
                    (j + 2 + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node;
                
                F_face_y[idx_face_y] += dt*(
                    a_r*(F_node_y[idx_node_B]   + F_node_y[idx_node_T]) +
                    b_r*(F_node_y[idx_node_BB]  + F_node_y[idx_node_TT]) +
                    c_r*(F_node_y[idx_node_BBB] + F_node_y[idx_node_TTT])
                    );
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and number of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_lo_2 = domain_lo[2];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        const int domain_dim_2 = domain_dims[2];
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        const int num_ghosts_0_flux_node = num_ghosts_flux_node[0];
        const int num_ghosts_1_flux_node = num_ghosts_flux_node[1];
        const int num_ghosts_2_flux_node = num_ghosts_flux_node[2];
        const int ghostcell_dim_0_flux_node = ghostcell_dims_flux_node[0];
        const int ghostcell_dim_1_flux_node = ghostcell_dims_flux_node[1];
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_face_y = i +
                        j*interior_dim_0 +
                        k*interior_dim_0*(interior_dim_1 + 1);
                    
                    const int idx_node_BBB = (i + num_ghosts_0_flux_node) +
                        (j - 3 + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node +
                        (k + num_ghosts_2_flux_node)*ghostcell_dim_0_flux_node*
                            ghostcell_dim_1_flux_node;
                    
                    const int idx_node_BB = (i + num_ghosts_0_flux_node) +
                        (j - 2 + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node +
                        (k + num_ghosts_2_flux_node)*ghostcell_dim_0_flux_node*
                            ghostcell_dim_1_flux_node;
                    
                    const int idx_node_B = (i + num_ghosts_0_flux_node) +
                        (j - 1 + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node +
                        (k + num_ghosts_2_flux_node)*ghostcell_dim_0_flux_node*
                            ghostcell_dim_1_flux_node;
                    
                    const int idx_node_T = (i + num_ghosts_0_flux_node) +
                        (j + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node +
                        (k + num_ghosts_2_flux_node)*ghostcell_dim_0_flux_node*
                            ghostcell_dim_1_flux_node;
                    
                    const int idx_node_TT = (i + num_ghosts_0_flux_node) +
                        (j + 1 + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node +
                        (k + num_ghosts_2_flux_node)*ghostcell_dim_0_flux_node*
                            ghostcell_dim_1_flux_node;
                    
                    const int idx_node_TTT = (i + num_ghosts_0_flux_node) +
                        (j + 2 + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node +
                        (k + num_ghosts_2_flux_node)*ghostcell_dim_0_flux_node*
                            ghostcell_dim_1_flux_node;
                    
                    F_face_y[idx_face_y] += dt*(
                        a_r*(F_node_y[idx_node_B]   + F_node_y[idx_node_T]) +
                        b_r*(F_node_y[idx_node_BB]  + F_node_y[idx_node_TT]) +
                        c_r*(F_node_y[idx_node_BBB] + F_node_y[idx_node_TTT])
                        );
                }
            }
        }
    }
}


/*
 * Kernel to reconstruct the flux using flux at nodes in z-direction.
 */
void
DiffusiveFluxReconstructorSixthOrder::reconstructFluxZ(
    double* F_face_z,
    const double* const F_node_z,
    const hier::IntVector& num_ghosts_flux_node,
    const hier::IntVector& ghostcell_dims_flux_node,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims,
    const hier::IntVector& interior_dims,
    const double& dt)
{
    const double a_n =  double(3)/double(4);
    const double b_n = -double(3)/double(20);
    const double c_n =  double(1)/double(60);
    
    const double a_r = a_n + b_n + c_n;
    const double b_r = b_n + c_n;
    const double c_r = c_n;
    
    /*
     * Get the local lower index, numbers of cells in each dimension and number of ghost cells.
     */
    
    const int domain_lo_0 = domain_lo[0];
    const int domain_lo_1 = domain_lo[1];
    const int domain_lo_2 = domain_lo[2];
    const int domain_dim_0 = domain_dims[0];
    const int domain_dim_1 = domain_dims[1];
    const int domain_dim_2 = domain_dims[2];
    
    const int interior_dim_0 = interior_dims[0];
    const int interior_dim_1 = interior_dims[1];
    
    const int num_ghosts_0_flux_node = num_ghosts_flux_node[0];
    const int num_ghosts_1_flux_node = num_ghosts_flux_node[1];
    const int num_ghosts_2_flux_node = num_ghosts_flux_node[2];
    const int ghostcell_dim_0_flux_node = ghostcell_dims_flux_node[0];
    const int ghostcell_dim_1_flux_node = ghostcell_dims_flux_node[1];
    
    for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
    {
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_face_z = i +
                    j*interior_dim_0 +
                    k*interior_dim_0*interior_dim_1;
                
                const int idx_node_BBB = (i + num_ghosts_0_flux_node) +
                    (j + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node +
                    (k - 3 + num_ghosts_2_flux_node)*ghostcell_dim_0_flux_node*
                        ghostcell_dim_1_flux_node;
                
                const int idx_node_BB = (i + num_ghosts_0_flux_node) +
                    (j + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node +
                    (k - 2 + num_ghosts_2_flux_node)*ghostcell_dim_0_flux_node*
                        ghostcell_dim_1_flux_node;
                
                const int idx_node_B = (i + num_ghosts_0_flux_node) +
                    (j + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node +
                    (k - 1 + num_ghosts_2_flux_node)*ghostcell_dim_0_flux_node*
                        ghostcell_dim_1_flux_node;
                
                const int idx_node_F = (i + num_ghosts_0_flux_node) +
                    (j + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node +
                    (k + num_ghosts_2_flux_node)*ghostcell_dim_0_flux_node*
                        ghostcell_dim_1_flux_node;
                
                const int idx_node_FF = (i + num_ghosts_0_flux_node) +
                    (j + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node +
                    (k + 1 + num_ghosts_2_flux_node)*ghostcell_dim_0_flux_node*
                        ghostcell_dim_1_flux_node;
                
                const int idx_node_FFF = (i + num_ghosts_0_flux_node) +
                    (j + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node +
                    (k + 2 + num_ghosts_2_flux_node)*ghostcell_dim_0_flux_node*
                        ghostcell_dim_1_flux_node;
                
                F_face_z[idx_face_z] += dt*(
                    a_r*(F_node_z[idx_node_B]   + F_node_z[idx_node_F]) +
                    b_r*(F_node_z[idx_node_BB]  + F_node_z[idx_node_FF]) +
                    c_r*(F_node_z[idx_node_BBB] + F_node_z[idx_node_FFF])
                    );
            }
        }
    }
}
