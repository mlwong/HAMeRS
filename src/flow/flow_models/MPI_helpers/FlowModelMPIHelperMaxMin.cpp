#include "flow/flow_models/MPI_helpers/FlowModelMPIHelperMaxMin.hpp"

#include "extn/patch_hierarchies/ExtendedFlattenedHierarchy.hpp"
#include "util/derivatives/DerivativeFirstOrder.hpp"

#include <limits>

/*
 * Compute maximum value with only x-direction as inhomogeneous direction.
 */
std::vector<double>
FlowModelMPIHelperMaxMin::getMaxQuantityWithInhomogeneousXDirection(
    const std::string quantity_name,
    const int component_idx,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
    std::vector<double> max_quantity;
    
    const int num_levels = d_patch_hierarchy->getNumberOfLevels();
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    HAMERS_SHARED_PTR<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *d_patch_hierarchy,
            0,
            num_levels - 1));
    
    if (d_dim == tbox::Dimension(1))
    {
        const int finest_level_dim_0 = d_finest_level_dims[0];
        
        double* u_max_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        max_quantity.resize(finest_level_dim_0);
        double* u_max_global = max_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_max_local[i]  = std::numeric_limits<double>::min();
            u_max_global[i] = std::numeric_limits<double>::min();
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                d_patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratio_to_coarest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarest_level/ratio_to_coarest_level;
            
            const int ratio_to_finest_level_0 = ratio_to_finest_level[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower indices.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                d_flow_model->registerDerivedVariables(num_subghosts_of_data);
                
                d_flow_model->allocateMemoryForDerivedCellData();
                
                d_flow_model->computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_quantity =
                    d_flow_model->getCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        /*
                         * Compute the linear index and update the max.
                         */
                        
                        const int idx = relative_idx_lo_0 + i + num_ghosts_0_quantity;
                        
                        for (int ii = 0; ii < ratio_to_finest_level_0; ii++)
                        {
                            const int idx_fine = (idx_lo_0 + i)*ratio_to_finest_level_0 + ii;
                            
                            u_max_local[idx_fine] = fmax(u_max_local[idx_fine], u[idx]);
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global max.
         */
        
        d_mpi.Allreduce(
            u_max_local,
            u_max_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_MAX);
        
        std::free(u_max_local);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_0 = d_finest_level_dims[0];
        
        double* u_max_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        max_quantity.resize(finest_level_dim_0);
        double* u_max_global = max_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_max_local[i]  = std::numeric_limits<double>::min();
            u_max_global[i] = std::numeric_limits<double>::min();
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                d_patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratio_to_coarest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarest_level/ratio_to_coarest_level;
            
            const int ratio_to_finest_level_0 = ratio_to_finest_level[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower indices.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                d_flow_model->registerDerivedVariables(num_subghosts_of_data);
                
                d_flow_model->allocateMemoryForDerivedCellData();
                
                d_flow_model->computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_quantity =
                    d_flow_model->getCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            /*
                             * Compute the linear index and update the max.
                             */
                            
                            const int idx = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                            
                            for (int ii = 0; ii < ratio_to_finest_level_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratio_to_finest_level_0 + ii;
                                
                                u_max_local[idx_fine] = fmax(u_max_local[idx_fine], u[idx]);
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global max.
         */
        
        d_mpi.Allreduce(
            u_max_local,
            u_max_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_MAX);
        
        std::free(u_max_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_0 = d_finest_level_dims[0];
        
        double* u_max_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        max_quantity.resize(finest_level_dim_0);
        double* u_max_global = max_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_max_local[i]  = std::numeric_limits<double>::min();
            u_max_global[i] = std::numeric_limits<double>::min();
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                d_patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratio_to_coarest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarest_level/ratio_to_coarest_level;
            
            const int ratio_to_finest_level_0 = ratio_to_finest_level[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower indices.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                d_flow_model->registerDerivedVariables(num_subghosts_of_data);
                
                d_flow_model->allocateMemoryForDerivedCellData();
                
                d_flow_model->computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_quantity =
                    d_flow_model->getCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int num_ghosts_2_quantity = num_ghosts_quantity[2];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                const int ghostcell_dim_1_quantity = ghostcell_dims_quantity[1];
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    const int interior_dim_2 = interior_dims[2];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int idx_lo_2 = index_lo[2];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    const int relative_idx_lo_2 = relative_index_lo[2];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                /*
                                 * Compute the linear index and update the max.
                                 */
                                
                                const int idx = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                    (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                        ghostcell_dim_1_quantity;
                                
                                for (int ii = 0; ii < ratio_to_finest_level_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratio_to_finest_level_0 + ii;
                                    
                                    u_max_local[idx_fine] = fmax(u_max_local[idx_fine], u[idx]);
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global max.
         */
        
        d_mpi.Allreduce(
            u_max_local,
            u_max_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_MAX);
        
        std::free(u_max_local);
    }
    
    return max_quantity;
}


/*
 * Compute minimum value with only x-direction as inhomogeneous direction.
 */
std::vector<double>
FlowModelMPIHelperMaxMin::getMinQuantityWithInhomogeneousXDirection(
    const std::string quantity_name,
    const int component_idx,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
    std::vector<double> min_quantity;
    
    const int num_levels = d_patch_hierarchy->getNumberOfLevels();
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    HAMERS_SHARED_PTR<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *d_patch_hierarchy,
            0,
            num_levels - 1));
    
    if (d_dim == tbox::Dimension(1))
    {
        const int finest_level_dim_0 = d_finest_level_dims[0];
        
        double* u_min_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        min_quantity.resize(finest_level_dim_0);
        double* u_min_global = min_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_min_local[i]  = std::numeric_limits<double>::max();
            u_min_global[i] = std::numeric_limits<double>::max();
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                d_patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratio_to_coarest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarest_level/ratio_to_coarest_level;
            
            const int ratio_to_finest_level_0 = ratio_to_finest_level[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower indices.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                d_flow_model->registerDerivedVariables(num_subghosts_of_data);
                
                d_flow_model->allocateMemoryForDerivedCellData();
                
                d_flow_model->computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_quantity =
                    d_flow_model->getCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        /*
                         * Compute the linear index and update the min.
                         */
                        
                        const int idx = relative_idx_lo_0 + i + num_ghosts_0_quantity;
                        
                        for (int ii = 0; ii < ratio_to_finest_level_0; ii++)
                        {
                            const int idx_fine = (idx_lo_0 + i)*ratio_to_finest_level_0 + ii;
                            
                            u_min_local[idx_fine] = fmin(u_min_local[idx_fine], u[idx]);
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global min.
         */
        
        d_mpi.Allreduce(
            u_min_local,
            u_min_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_MIN);
        
        std::free(u_min_local);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_0 = d_finest_level_dims[0];
        
        double* u_min_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        min_quantity.resize(finest_level_dim_0);
        double* u_min_global = min_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_min_local[i]  = std::numeric_limits<double>::max();
            u_min_global[i] = std::numeric_limits<double>::max();
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                d_patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratio_to_coarest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarest_level/ratio_to_coarest_level;
            
            const int ratio_to_finest_level_0 = ratio_to_finest_level[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower indices.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                d_flow_model->registerDerivedVariables(num_subghosts_of_data);
                
                d_flow_model->allocateMemoryForDerivedCellData();
                
                d_flow_model->computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_quantity =
                    d_flow_model->getCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            /*
                             * Compute the linear index and update the min.
                             */
                            
                            const int idx = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                            
                            for (int ii = 0; ii < ratio_to_finest_level_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratio_to_finest_level_0 + ii;
                                
                                u_min_local[idx_fine] = fmin(u_min_local[idx_fine], u[idx]);
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global min.
         */
        
        d_mpi.Allreduce(
            u_min_local,
            u_min_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_MIN);
        
        std::free(u_min_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_0 = d_finest_level_dims[0];
        
        double* u_min_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        min_quantity.resize(finest_level_dim_0);
        double* u_min_global = min_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_min_local[i]  = std::numeric_limits<double>::max();
            u_min_global[i] = std::numeric_limits<double>::max();
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                d_patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratio_to_coarest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarest_level/ratio_to_coarest_level;
            
            const int ratio_to_finest_level_0 = ratio_to_finest_level[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower indices.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                d_flow_model->registerDerivedVariables(num_subghosts_of_data);
                
                d_flow_model->allocateMemoryForDerivedCellData();
                
                d_flow_model->computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_quantity =
                    d_flow_model->getCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int num_ghosts_2_quantity = num_ghosts_quantity[2];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                const int ghostcell_dim_1_quantity = ghostcell_dims_quantity[1];
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    const int interior_dim_2 = interior_dims[2];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int idx_lo_2 = index_lo[2];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    const int relative_idx_lo_2 = relative_index_lo[2];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                /*
                                 * Compute the linear index and update the min.
                                 */
                                
                                const int idx = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                    (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                        ghostcell_dim_1_quantity;
                                
                                for (int ii = 0; ii < ratio_to_finest_level_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratio_to_finest_level_0 + ii;
                                    
                                    u_min_local[idx_fine] = fmin(u_min_local[idx_fine], u[idx]);
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global min.
         */
        
        d_mpi.Allreduce(
            u_min_local,
            u_min_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_MIN);
        
        std::free(u_min_local);
    }
    
    return min_quantity;
}


/*
 * Compute maximum value with only y-direction as inhomogeneous direction.
 */
std::vector<double>
FlowModelMPIHelperMaxMin::getMaxQuantityWithInhomogeneousYDirection(
    const std::string quantity_name,
    const int component_idx,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
    std::vector<double> max_quantity;
    
    const int num_levels = d_patch_hierarchy->getNumberOfLevels();
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    HAMERS_SHARED_PTR<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *d_patch_hierarchy,
            0,
            num_levels - 1));
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperMaxMin::getMaxQuantityWithInhomogeneousYDirection():\n"
            << "Cannot compute maximum value with only y-direction as inhomogeneous direction for"
            << " one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_1 = d_finest_level_dims[1];
        
        double* u_max_local = (double*)std::malloc(finest_level_dim_1*sizeof(double));
        
        max_quantity.resize(finest_level_dim_1);
        double* u_max_global = max_quantity.data();
        
        for (int j = 0; j < finest_level_dim_1; j++)
        {
            u_max_local[j]  = std::numeric_limits<double>::min();
            u_max_global[j] = std::numeric_limits<double>::min();
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                d_patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratio_to_coarest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarest_level/ratio_to_coarest_level;
            
            const int ratio_to_finest_level_1 = ratio_to_finest_level[1];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower indices.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                d_flow_model->registerDerivedVariables(num_subghosts_of_data);
                
                d_flow_model->allocateMemoryForDerivedCellData();
                
                d_flow_model->computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_quantity =
                    d_flow_model->getCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            /*
                             * Compute the linear index and update the max.
                             */
                            
                            const int idx = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                            
                            for (int jj = 0; jj < ratio_to_finest_level_1; jj++)
                            {
                                const int idx_fine = (idx_lo_1 + j)*ratio_to_finest_level_1 + jj;
                                
                                u_max_local[idx_fine] = fmax(u_max_local[idx_fine], u[idx]);
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global max.
         */
        
        d_mpi.Allreduce(
            u_max_local,
            u_max_global,
            finest_level_dim_1,
            MPI_DOUBLE,
            MPI_MAX);
        
        std::free(u_max_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_1 = d_finest_level_dims[1];
        
        double* u_max_local = (double*)std::malloc(finest_level_dim_1*sizeof(double));
        
        max_quantity.resize(finest_level_dim_1);
        double* u_max_global = max_quantity.data();
        
        for (int j = 0; j < finest_level_dim_1; j++)
        {
            u_max_local[j]  = std::numeric_limits<double>::min();
            u_max_global[j] = std::numeric_limits<double>::min();
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                d_patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratio_to_coarest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarest_level/ratio_to_coarest_level;
            
            const int ratio_to_finest_level_1 = ratio_to_finest_level[1];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower indices.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                d_flow_model->registerDerivedVariables(num_subghosts_of_data);
                
                d_flow_model->allocateMemoryForDerivedCellData();
                
                d_flow_model->computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_quantity =
                    d_flow_model->getCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int num_ghosts_2_quantity = num_ghosts_quantity[2];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                const int ghostcell_dim_1_quantity = ghostcell_dims_quantity[1];
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    const int interior_dim_2 = interior_dims[2];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int idx_lo_2 = index_lo[2];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    const int relative_idx_lo_2 = relative_index_lo[2];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                /*
                                 * Compute the linear index and update the max.
                                 */
                                
                                const int idx = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                    (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                        ghostcell_dim_1_quantity;
                                
                                for (int jj = 0; jj < ratio_to_finest_level_1; jj++)
                                {
                                    const int idx_fine = (idx_lo_1 + j)*ratio_to_finest_level_1 + jj;
                                    
                                    u_max_local[idx_fine] = fmax(u_max_local[idx_fine], u[idx]);
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global max.
         */
        
        d_mpi.Allreduce(
            u_max_local,
            u_max_global,
            finest_level_dim_1,
            MPI_DOUBLE,
            MPI_MAX);
        
        std::free(u_max_local);
    }
    
    return max_quantity;
}


/*
 * Compute minimum value with only y-direction as inhomogeneous direction.
 */
std::vector<double>
FlowModelMPIHelperMaxMin::getMinQuantityWithInhomogeneousYDirection(
    const std::string quantity_name,
    const int component_idx,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
    std::vector<double> min_quantity;
    
    const int num_levels = d_patch_hierarchy->getNumberOfLevels();
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    HAMERS_SHARED_PTR<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *d_patch_hierarchy,
            0,
            num_levels - 1));
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperMaxMin::getMinQuantityWithInhomogeneousYDirection():\n"
            << "Cannot compute minimum value with only y-direction as inhomogeneous direction for"
            << " one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_1 = d_finest_level_dims[1];
        
        double* u_min_local = (double*)std::malloc(finest_level_dim_1*sizeof(double));
        
        min_quantity.resize(finest_level_dim_1);
        double* u_min_global = min_quantity.data();
        
        for (int j = 0; j < finest_level_dim_1; j++)
        {
            u_min_local[j]  = std::numeric_limits<double>::max();
            u_min_global[j] = std::numeric_limits<double>::max();
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                d_patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratio_to_coarest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarest_level/ratio_to_coarest_level;
            
            const int ratio_to_finest_level_1 = ratio_to_finest_level[1];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower indices.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                d_flow_model->registerDerivedVariables(num_subghosts_of_data);
                
                d_flow_model->allocateMemoryForDerivedCellData();
                
                d_flow_model->computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_quantity =
                    d_flow_model->getCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            /*
                             * Compute the linear index and update the min.
                             */
                            
                            const int idx = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                            
                            for (int jj = 0; jj < ratio_to_finest_level_1; jj++)
                            {
                                const int idx_fine = (idx_lo_1 + j)*ratio_to_finest_level_1 + jj;
                                
                                u_min_local[idx_fine] = fmin(u_min_local[idx_fine], u[idx]);
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global min.
         */
        
        d_mpi.Allreduce(
            u_min_local,
            u_min_global,
            finest_level_dim_1,
            MPI_DOUBLE,
            MPI_MIN);
        
        std::free(u_min_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_1 = d_finest_level_dims[1];
        
        double* u_min_local = (double*)std::malloc(finest_level_dim_1*sizeof(double));
        
        min_quantity.resize(finest_level_dim_1);
        double* u_min_global = min_quantity.data();
        
        for (int j = 0; j < finest_level_dim_1; j++)
        {
            u_min_local[j]  = std::numeric_limits<double>::max();
            u_min_global[j] = std::numeric_limits<double>::max();
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                d_patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratio_to_coarest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarest_level/ratio_to_coarest_level;
            
            const int ratio_to_finest_level_1 = ratio_to_finest_level[1];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower indices.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                d_flow_model->registerDerivedVariables(num_subghosts_of_data);
                
                d_flow_model->allocateMemoryForDerivedCellData();
                
                d_flow_model->computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_quantity =
                    d_flow_model->getCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int num_ghosts_2_quantity = num_ghosts_quantity[2];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                const int ghostcell_dim_1_quantity = ghostcell_dims_quantity[1];
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    const int interior_dim_2 = interior_dims[2];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int idx_lo_2 = index_lo[2];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    const int relative_idx_lo_2 = relative_index_lo[2];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                /*
                                 * Compute the linear index and update the min.
                                 */
                                
                                const int idx = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                    (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                        ghostcell_dim_1_quantity;
                                
                                for (int jj = 0; jj < ratio_to_finest_level_1; jj++)
                                {
                                    const int idx_fine = (idx_lo_1 + j)*ratio_to_finest_level_1 + jj;
                                    
                                    u_min_local[idx_fine] = fmin(u_min_local[idx_fine], u[idx]);
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global min.
         */
        
        d_mpi.Allreduce(
            u_min_local,
            u_min_global,
            finest_level_dim_1,
            MPI_DOUBLE,
            MPI_MIN);
        
        std::free(u_min_local);
    }
    
    return min_quantity;
}


/*
 * Compute maximum value with only z-direction as inhomogeneous direction.
 */
std::vector<double>
FlowModelMPIHelperMaxMin::getMaxQuantityWithInhomogeneousZDirection(
    const std::string quantity_name,
    const int component_idx,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
    std::vector<double> max_quantity;
    
    const int num_levels = d_patch_hierarchy->getNumberOfLevels();
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    HAMERS_SHARED_PTR<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *d_patch_hierarchy,
            0,
            num_levels - 1));
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperMaxMin::getMaxQuantityWithInhomogeneousZDirection():\n"
            << "Cannot compute maximum value with only z-direction as inhomogeneous direction for"
            << " one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperMaxMin::getMaxQuantityWithInhomogeneousZDirection():\n"
            << "Cannot compute maximum value with only z-direction as inhomogeneous direction for"
            << " two-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_2 = d_finest_level_dims[2];
        
        double* u_max_local = (double*)std::malloc(finest_level_dim_2*sizeof(double));
        
        max_quantity.resize(finest_level_dim_2);
        double* u_max_global = max_quantity.data();
        
        for (int k = 0; k < finest_level_dim_2; k++)
        {
            u_max_local[k]  = std::numeric_limits<double>::min();
            u_max_global[k] = std::numeric_limits<double>::min();
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                d_patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratio_to_coarest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarest_level/ratio_to_coarest_level;
            
            const int ratio_to_finest_level_2 = ratio_to_finest_level[2];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower indices.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                d_flow_model->registerDerivedVariables(num_subghosts_of_data);
                
                d_flow_model->allocateMemoryForDerivedCellData();
                
                d_flow_model->computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_quantity =
                    d_flow_model->getCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int num_ghosts_2_quantity = num_ghosts_quantity[2];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                const int ghostcell_dim_1_quantity = ghostcell_dims_quantity[1];
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    const int interior_dim_2 = interior_dims[2];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int idx_lo_2 = index_lo[2];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    const int relative_idx_lo_2 = relative_index_lo[2];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                /*
                                 * Compute the linear index and update the max.
                                 */
                                
                                const int idx = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                    (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                        ghostcell_dim_1_quantity;
                                
                                for (int kk = 0; kk < ratio_to_finest_level_2; kk++)
                                {
                                    const int idx_fine = (idx_lo_2 + k)*ratio_to_finest_level_2 + kk;
                                    
                                    u_max_local[idx_fine] = fmax(u_max_local[idx_fine], u[idx]);
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global max.
         */
        
        d_mpi.Allreduce(
            u_max_local,
            u_max_global,
            finest_level_dim_2,
            MPI_DOUBLE,
            MPI_MAX);
        
        std::free(u_max_local);
    }
    
    return max_quantity;
}


/*
 * Compute minimum value with only z-direction as inhomogeneous direction.
 */
std::vector<double>
FlowModelMPIHelperMaxMin::getMinQuantityWithInhomogeneousZDirection(
    const std::string quantity_name,
    const int component_idx,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
    std::vector<double> min_quantity;
    
    const int num_levels = d_patch_hierarchy->getNumberOfLevels();
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    HAMERS_SHARED_PTR<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *d_patch_hierarchy,
            0,
            num_levels - 1));
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperMaxMin::getMinQuantityWithInhomogeneousZDirection():\n"
            << "Cannot compute minimum value with only z-direction as inhomogeneous direction for"
            << " one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperMaxMin::getMinQuantityWithInhomogeneousZDirection():\n"
            << "Cannot compute minimum value with only z-direction as inhomogeneous direction for"
            << " two-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_2 = d_finest_level_dims[2];
        
        double* u_min_local = (double*)std::malloc(finest_level_dim_2*sizeof(double));
        
        min_quantity.resize(finest_level_dim_2);
        double* u_min_global = min_quantity.data();
        
        for (int k = 0; k < finest_level_dim_2; k++)
        {
            u_min_local[k]  = std::numeric_limits<double>::max();
            u_min_global[k] = std::numeric_limits<double>::max();
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                d_patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratio_to_coarest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarest_level/ratio_to_coarest_level;
            
            const int ratio_to_finest_level_2 = ratio_to_finest_level[2];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower indices.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                d_flow_model->registerDerivedVariables(num_subghosts_of_data);
                
                d_flow_model->allocateMemoryForDerivedCellData();
                
                d_flow_model->computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_quantity =
                    d_flow_model->getCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int num_ghosts_2_quantity = num_ghosts_quantity[2];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                const int ghostcell_dim_1_quantity = ghostcell_dims_quantity[1];
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    const int interior_dim_2 = interior_dims[2];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int idx_lo_2 = index_lo[2];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    const int relative_idx_lo_2 = relative_index_lo[2];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                /*
                                 * Compute the linear index and update the min.
                                 */
                                
                                const int idx = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                    (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                        ghostcell_dim_1_quantity;
                                
                                for (int kk = 0; kk < ratio_to_finest_level_2; kk++)
                                {
                                    const int idx_fine = (idx_lo_2 + k)*ratio_to_finest_level_2 + kk;
                                    
                                    u_min_local[idx_fine] = fmin(u_min_local[idx_fine], u[idx]);
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global min.
         */
        
        d_mpi.Allreduce(
            u_min_local,
            u_min_global,
            finest_level_dim_2,
            MPI_DOUBLE,
            MPI_MIN);
        
        std::free(u_min_local);
    }
    
    return min_quantity;
}


/*
 * Compute maximum location within quantity bounds in x-direction.
 */
double
FlowModelMPIHelperMaxMin::getMaxLocationWithinQuantityBoundsInXDirection(
    const std::string quantity_name,
    const int component_idx,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double bound_lo,
    const double bound_hi) const
{
    double location_x_max_global;
    
    const int num_levels = d_patch_hierarchy->getNumberOfLevels();
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    HAMERS_SHARED_PTR<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *d_patch_hierarchy,
            0,
            num_levels - 1));
    
    /*
     * Get the lower indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    
    if (d_dim == tbox::Dimension(1))
    {
        double location_x_max_local = x_lo[0];
        location_x_max_global       = x_lo[0];
        
        const double half = double(1)/double(2);
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                d_patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower indices, grid spacings and the lower spatial coordinates.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                const double* const x_lo_patch = patch_geom->getXLower();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                d_flow_model->registerDerivedVariables(num_subghosts_of_data);
                
                d_flow_model->allocateMemoryForDerivedCellData();
                
                d_flow_model->computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_quantity =
                    d_flow_model->getCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        /*
                         * Compute the linear index and update the max.
                         */
                        
                        const int idx = relative_idx_lo_0 + i + num_ghosts_0_quantity;
                        
                        const double x = (relative_idx_lo_0 + i + half)*dx[0] + x_lo_patch[0];
                        
                        if (u[idx] > bound_lo && u[idx] < bound_hi)
                        {
                            if (x > location_x_max_local)
                            {
                                location_x_max_local = x;
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global max.
         */
        
        d_mpi.Allreduce(
            &location_x_max_local,
            &location_x_max_global,
            1,
            MPI_DOUBLE,
            MPI_MAX);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        double location_x_max_local = x_lo[0];
        location_x_max_global       = x_lo[0];
        
        const double half = double(1)/double(2);
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                d_patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower indices, grid spacings and the lower spatial coordinates.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                const double* const x_lo_patch = patch_geom->getXLower();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                d_flow_model->registerDerivedVariables(num_subghosts_of_data);
                
                d_flow_model->allocateMemoryForDerivedCellData();
                
                d_flow_model->computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_quantity =
                    d_flow_model->getCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            /*
                             * Compute the linear index and update the max.
                             */
                            
                            const int idx = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                            
                            const double x = (relative_idx_lo_0 + i + half)*dx[0] + x_lo_patch[0];
                            
                            if (u[idx] > bound_lo && u[idx] < bound_hi)
                            {
                                if (x > location_x_max_local)
                                {
                                    location_x_max_local = x;
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global max.
         */
        
        d_mpi.Allreduce(
            &location_x_max_local,
            &location_x_max_global,
            1,
            MPI_DOUBLE,
            MPI_MAX);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        double location_x_max_local = x_lo[0];
        location_x_max_global       = x_lo[0];
        
        const double half = double(1)/double(2);
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                d_patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower indices, grid spacings and the lower spatial coordinates.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                const double* const x_lo_patch = patch_geom->getXLower();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                d_flow_model->registerDerivedVariables(num_subghosts_of_data);
                
                d_flow_model->allocateMemoryForDerivedCellData();
                
                d_flow_model->computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_quantity =
                    d_flow_model->getCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int num_ghosts_2_quantity = num_ghosts_quantity[2];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                const int ghostcell_dim_1_quantity = ghostcell_dims_quantity[1];
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    const int interior_dim_2 = interior_dims[2];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int idx_lo_2 = index_lo[2];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    const int relative_idx_lo_2 = relative_index_lo[2];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                /*
                                 * Compute the linear index and update the max.
                                 */
                                
                                const int idx = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                    (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                        ghostcell_dim_1_quantity;
                                
                                const double x = (relative_idx_lo_0 + i + half)*dx[0] + x_lo_patch[0];
                                
                                if (u[idx] > bound_lo && u[idx] < bound_hi)
                                {
                                    if (x > location_x_max_local)
                                    {
                                        location_x_max_local = x;
                                    }
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global max.
         */
        
        d_mpi.Allreduce(
            &location_x_max_local,
            &location_x_max_global,
            1,
            MPI_DOUBLE,
            MPI_MAX);
    }
    
    return location_x_max_global;
}


/*
 * Compute minimum location within quantity bounds in x-direction.
 */
double
FlowModelMPIHelperMaxMin::getMinLocationWithinQuantityBoundsInXDirection(
    const std::string quantity_name,
    const int component_idx,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double bound_lo,
    const double bound_hi) const
{
    double location_x_min_global;
    
    const int num_levels = d_patch_hierarchy->getNumberOfLevels();
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    HAMERS_SHARED_PTR<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *d_patch_hierarchy,
            0,
            num_levels - 1));
    
    /*
     * Get the upper indices of the physical domain.
     */
    
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        double location_x_min_local = x_hi[0];
        location_x_min_global       = x_hi[0];
        
        const double half = double(1)/double(2);
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                d_patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower indices, grid spacings and the lower spatial coordinates.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                const double* const x_lo_patch = patch_geom->getXLower();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                d_flow_model->registerDerivedVariables(num_subghosts_of_data);
                
                d_flow_model->allocateMemoryForDerivedCellData();
                
                d_flow_model->computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_quantity =
                    d_flow_model->getCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        /*
                         * Compute the linear index and update the min.
                         */
                        
                        const int idx = relative_idx_lo_0 + i + num_ghosts_0_quantity;
                        
                        const double x = (relative_idx_lo_0 + i + half)*dx[0] + x_lo_patch[0];
                        
                        if (u[idx] > bound_lo && u[idx] < bound_hi)
                        {
                            if (x < location_x_min_local)
                            {
                                location_x_min_local = x;
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global min.
         */
        
        d_mpi.Allreduce(
            &location_x_min_local,
            &location_x_min_global,
            1,
            MPI_DOUBLE,
            MPI_MIN);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        double location_x_min_local = x_hi[0];
        location_x_min_global       = x_hi[0];
        
        const double half = double(1)/double(2);
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                d_patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower indices, grid spacings and the lower spatial coordinates.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                const double* const x_lo_patch = patch_geom->getXLower();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                d_flow_model->registerDerivedVariables(num_subghosts_of_data);
                
                d_flow_model->allocateMemoryForDerivedCellData();
                
                d_flow_model->computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_quantity =
                    d_flow_model->getCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            /*
                             * Compute the linear index and update the min.
                             */
                            
                            const int idx = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                            
                            const double x = (relative_idx_lo_0 + i + half)*dx[0] + x_lo_patch[0];
                            
                            if (u[idx] > bound_lo && u[idx] < bound_hi)
                            {
                                if (x < location_x_min_local)
                                {
                                    location_x_min_local = x;
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global min.
         */
        
        d_mpi.Allreduce(
            &location_x_min_local,
            &location_x_min_global,
            1,
            MPI_DOUBLE,
            MPI_MIN);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        double location_x_min_local = x_hi[0];
        location_x_min_global       = x_hi[0];
        
        const double half = double(1)/double(2);
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                d_patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower indices, grid spacings and the lower spatial coordinates.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                const double* const x_lo_patch = patch_geom->getXLower();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                d_flow_model->registerDerivedVariables(num_subghosts_of_data);
                
                d_flow_model->allocateMemoryForDerivedCellData();
                
                d_flow_model->computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_quantity =
                    d_flow_model->getCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int num_ghosts_2_quantity = num_ghosts_quantity[2];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                const int ghostcell_dim_1_quantity = ghostcell_dims_quantity[1];
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    const int interior_dim_2 = interior_dims[2];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int idx_lo_2 = index_lo[2];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    const int relative_idx_lo_2 = relative_index_lo[2];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                /*
                                 * Compute the linear index and update the min.
                                 */
                                
                                const int idx = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                    (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                        ghostcell_dim_1_quantity;
                                
                                const double x = (relative_idx_lo_0 + i + half)*dx[0] + x_lo_patch[0];
                                
                                if (u[idx] > bound_lo && u[idx] < bound_hi)
                                {
                                    if (x < location_x_min_local)
                                    {
                                        location_x_min_local = x;
                                    }
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global min.
         */
        
        d_mpi.Allreduce(
            &location_x_min_local,
            &location_x_min_global,
            1,
            MPI_DOUBLE,
            MPI_MIN);
    }
    
    return location_x_min_global;
}


/*
 * Compute maximum location within quantity bounds in y-direction.
 */
double
FlowModelMPIHelperMaxMin::getMaxLocationWithinQuantityBoundsInYDirection(
    const std::string quantity_name,
    const int component_idx,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double bound_lo,
    const double bound_hi) const
{
    double location_y_max_global;
    
    const int num_levels = d_patch_hierarchy->getNumberOfLevels();
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    HAMERS_SHARED_PTR<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *d_patch_hierarchy,
            0,
            num_levels - 1));
    
    /*
     * Get the lower indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperMaxMin::getMaxLocationWithinQuantityBoundsInYDirection():\n"
            << "Cannot compute maximum location within quantity bounds in y-direction for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        double location_y_max_local = x_lo[1];
        location_y_max_global       = x_lo[1];
        
        const double half = double(1)/double(2);
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                d_patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower indices, grid spacings and the lower spatial coordinates.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                const double* const x_lo_patch = patch_geom->getXLower();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                d_flow_model->registerDerivedVariables(num_subghosts_of_data);
                
                d_flow_model->allocateMemoryForDerivedCellData();
                
                d_flow_model->computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_quantity =
                    d_flow_model->getCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            /*
                             * Compute the linear index and update the max.
                             */
                            
                            const int idx = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                            
                            const double y = (relative_idx_lo_1 + j + half)*dx[1] + x_lo_patch[1];
                            
                            if (u[idx] > bound_lo && u[idx] < bound_hi)
                            {
                                if (y > location_y_max_local)
                                {
                                    location_y_max_local = y;
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global max.
         */
        
        d_mpi.Allreduce(
            &location_y_max_local,
            &location_y_max_global,
            1,
            MPI_DOUBLE,
            MPI_MAX);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        double location_y_max_local = x_lo[1];
        location_y_max_global       = x_lo[1];
        
        const double half = double(1)/double(2);
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                d_patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower indices, grid spacings and the lower spatial coordinates.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                const double* const x_lo_patch = patch_geom->getXLower();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                d_flow_model->registerDerivedVariables(num_subghosts_of_data);
                
                d_flow_model->allocateMemoryForDerivedCellData();
                
                d_flow_model->computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_quantity =
                    d_flow_model->getCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int num_ghosts_2_quantity = num_ghosts_quantity[2];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                const int ghostcell_dim_1_quantity = ghostcell_dims_quantity[1];
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    const int interior_dim_2 = interior_dims[2];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int idx_lo_2 = index_lo[2];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    const int relative_idx_lo_2 = relative_index_lo[2];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                /*
                                 * Compute the linear index and update the max.
                                 */
                                
                                const int idx = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                    (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                        ghostcell_dim_1_quantity;
                                
                                const double y = (relative_idx_lo_1 + j + half)*dx[1] + x_lo_patch[1];
                                
                                if (u[idx] > bound_lo && u[idx] < bound_hi)
                                {
                                    if (y > location_y_max_local)
                                    {
                                        location_y_max_local = y;
                                    }
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global max.
         */
        
        d_mpi.Allreduce(
            &location_y_max_local,
            &location_y_max_global,
            1,
            MPI_DOUBLE,
            MPI_MAX);
    }
    
    return location_y_max_global;
}


/*
 * Compute minimum location within quantity bounds in y-direction.
 */
double
FlowModelMPIHelperMaxMin::getMinLocationWithinQuantityBoundsInYDirection(
    const std::string quantity_name,
    const int component_idx,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double bound_lo,
    const double bound_hi) const
{
    double location_y_min_global;
    
    const int num_levels = d_patch_hierarchy->getNumberOfLevels();
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    HAMERS_SHARED_PTR<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *d_patch_hierarchy,
            0,
            num_levels - 1));
    
    /*
     * Get the upper indices of the physical domain.
     */
    
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperMaxMin::getMinLocationWithinQuantityBoundsInYDirection():\n"
            << "Cannot compute maximum location within quantity bounds in y-direction for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        double location_y_min_local = x_hi[1];
        location_y_min_global       = x_hi[1];
        
        const double half = double(1)/double(2);
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                d_patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower indices, grid spacings and the lower spatial coordinates.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                const double* const x_lo_patch = patch_geom->getXLower();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                d_flow_model->registerDerivedVariables(num_subghosts_of_data);
                
                d_flow_model->allocateMemoryForDerivedCellData();
                
                d_flow_model->computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_quantity =
                    d_flow_model->getCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            /*
                             * Compute the linear index and update the max.
                             */
                            
                            const int idx = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                            
                            const double y = (relative_idx_lo_1 + j + half)*dx[1] + x_lo_patch[1];
                            
                            if (u[idx] > bound_lo && u[idx] < bound_hi)
                            {
                                if (y < location_y_min_local)
                                {
                                    location_y_min_local = y;
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global min.
         */
        
        d_mpi.Allreduce(
            &location_y_min_local,
            &location_y_min_global,
            1,
            MPI_DOUBLE,
            MPI_MIN);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        double location_y_min_local = x_hi[1];
        location_y_min_global       = x_hi[1];
        
        const double half = double(1)/double(2);
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                d_patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower indices, grid spacings and the lower spatial coordinates.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                const double* const x_lo_patch = patch_geom->getXLower();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                d_flow_model->registerDerivedVariables(num_subghosts_of_data);
                
                d_flow_model->allocateMemoryForDerivedCellData();
                
                d_flow_model->computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_quantity =
                    d_flow_model->getCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int num_ghosts_2_quantity = num_ghosts_quantity[2];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                const int ghostcell_dim_1_quantity = ghostcell_dims_quantity[1];
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    const int interior_dim_2 = interior_dims[2];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int idx_lo_2 = index_lo[2];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    const int relative_idx_lo_2 = relative_index_lo[2];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                /*
                                 * Compute the linear index and update the max.
                                 */
                                
                                const int idx = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                    (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                        ghostcell_dim_1_quantity;
                                
                                const double y = (relative_idx_lo_1 + j + half)*dx[1] + x_lo_patch[1];
                                
                                if (u[idx] > bound_lo && u[idx] < bound_hi)
                                {
                                    if (y < location_y_min_local)
                                    {
                                        location_y_min_local = y;
                                    }
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global min.
         */
        
        d_mpi.Allreduce(
            &location_y_min_local,
            &location_y_min_global,
            1,
            MPI_DOUBLE,
            MPI_MIN);
    }
    
    return location_y_min_global;
}


/*
 * Compute maximum location within quantity bounds in z-direction.
 */
double
FlowModelMPIHelperMaxMin::getMaxLocationWithinQuantityBoundsInZDirection(
    const std::string quantity_name,
    const int component_idx,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double bound_lo,
    const double bound_hi) const
{
    double location_z_max_global;
    
    const int num_levels = d_patch_hierarchy->getNumberOfLevels();
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    HAMERS_SHARED_PTR<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *d_patch_hierarchy,
            0,
            num_levels - 1));
    
    /*
     * Get the lower indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperMaxMin::getMaxLocationWithinQuantityBoundsInZDirection():\n"
            << "Cannot compute maximum location within quantity bounds in z-direction for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperMaxMin::getMaxLocationWithinQuantityBoundsInZDirection():\n"
            << "Cannot compute maximum location within quantity bounds in z-direction for two-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        double location_z_max_local = x_lo[2];
        location_z_max_global       = x_lo[2];
        
        const double half = double(1)/double(2);
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                d_patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower indices, grid spacings and the lower spatial coordinates.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                const double* const x_lo_patch = patch_geom->getXLower();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                d_flow_model->registerDerivedVariables(num_subghosts_of_data);
                
                d_flow_model->allocateMemoryForDerivedCellData();
                
                d_flow_model->computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_quantity =
                    d_flow_model->getCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int num_ghosts_2_quantity = num_ghosts_quantity[2];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                const int ghostcell_dim_1_quantity = ghostcell_dims_quantity[1];
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    const int interior_dim_2 = interior_dims[2];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int idx_lo_2 = index_lo[2];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    const int relative_idx_lo_2 = relative_index_lo[2];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                /*
                                 * Compute the linear index and update the max.
                                 */
                                
                                const int idx = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                    (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                        ghostcell_dim_1_quantity;
                                
                                const double z = (relative_idx_lo_2 + k + half)*dx[2] + x_lo_patch[2];
                                
                                if (u[idx] > bound_lo && u[idx] < bound_hi)
                                {
                                    if (z > location_z_max_local)
                                    {
                                        location_z_max_local = z;
                                    }
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global max.
         */
        
        d_mpi.Allreduce(
            &location_z_max_local,
            &location_z_max_global,
            1,
            MPI_DOUBLE,
            MPI_MAX);
    }
    
    return location_z_max_global;
}


/*
 * Compute minimum location within quantity bounds in z-direction.
 */
double
FlowModelMPIHelperMaxMin::getMinLocationWithinQuantityBoundsInZDirection(
    const std::string quantity_name,
    const int component_idx,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double bound_lo,
    const double bound_hi) const
{
    double location_z_min_global;
    
    const int num_levels = d_patch_hierarchy->getNumberOfLevels();
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    HAMERS_SHARED_PTR<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *d_patch_hierarchy,
            0,
            num_levels - 1));
    
    /*
     * Get the upper indices of the physical domain.
     */
    
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperMaxMin::getMinLocationWithinQuantityBoundsInZDirection():\n"
            << "Cannot compute maximum location within quantity bounds in z-direction for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperMaxMin::getMinLocationWithinQuantityBoundsInZDirection():\n"
            << "Cannot compute maximum location within quantity bounds in z-direction for two-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        double location_z_min_local = x_hi[2];
        location_z_min_global       = x_hi[2];
        
        const double half = double(1)/double(2);
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                d_patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower indices, grid spacings and the lower spatial coordinates.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                const double* const x_lo_patch = patch_geom->getXLower();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model->registerPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                d_flow_model->registerDerivedVariables(num_subghosts_of_data);
                
                d_flow_model->allocateMemoryForDerivedCellData();
                
                d_flow_model->computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_quantity =
                    d_flow_model->getCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int num_ghosts_2_quantity = num_ghosts_quantity[2];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                const int ghostcell_dim_1_quantity = ghostcell_dims_quantity[1];
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    const int interior_dim_2 = interior_dims[2];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int idx_lo_2 = index_lo[2];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    const int relative_idx_lo_2 = relative_index_lo[2];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                /*
                                 * Compute the linear index and update the max.
                                 */
                                
                                const int idx = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                    (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                        ghostcell_dim_1_quantity;
                                
                                const double z = (relative_idx_lo_2 + k + half)*dx[2] + x_lo_patch[2];
                                
                                if (u[idx] > bound_lo && u[idx] < bound_hi)
                                {
                                    if (z < location_z_min_local)
                                    {
                                        location_z_min_local = z;
                                    }
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global min.
         */
        
        d_mpi.Allreduce(
            &location_z_min_local,
            &location_z_min_global,
            1,
            MPI_DOUBLE,
            MPI_MIN);
    }
    
    return location_z_min_global;
}


/*
 * Compute maximum value of absolute value of gradient with only x-direction as inhomogeneous direction.
 */
std::vector<double>
FlowModelMPIHelperMaxMin::getMaxAbsoluteGradientWithInhomogeneousXDirection(
    const std::string quantity_name,
    const int component_idx,
    const int derivative_direction,
    const int num_ghosts_derivative,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
    if (d_dim == tbox::Dimension(1) && derivative_direction > 0)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperMaxMin::getMaxAbsoluteGradientWithInhomogeneousXDirection():\n"
            << "Cannot take derivative for one-dimensional problem!\n"
            << "derivative_direction = " << derivative_direction << " given!\n"
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2) && derivative_direction > 1)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperMaxMin::getMaxAbsoluteGradientWithInhomogeneousXDirection():\n"
            << "Cannot take derivative for two-dimensional problem!\n"
            << "derivative_direction = " << derivative_direction << " given!\n"
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(3) && derivative_direction > 2)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperMaxMin::getMaxAbsoluteGradientWithInhomogeneousXDirection():\n"
            << "Cannot take derivative for three-dimensional problem!\n"
            << "derivative_direction = " << derivative_direction << " given!\n"
            << std::endl);
    }
    
    hier::IntVector num_ghosts_der = hier::IntVector::getZero(d_dim);
    num_ghosts_der[derivative_direction] = num_ghosts_derivative;
    
    std::vector<double> max_abs_derivative;
    
    const int num_levels = d_patch_hierarchy->getNumberOfLevels();
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    HAMERS_SHARED_PTR<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *d_patch_hierarchy,
            0,
            num_levels - 1));
    
    if (d_dim == tbox::Dimension(1))
    {
        HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_x(
            new DerivativeFirstOrder(
                "first order derivative in x-direction",
                d_dim,
                DIRECTION::X_DIRECTION,
                num_ghosts_derivative));
        
        const int finest_level_dim_0 = d_finest_level_dims[0];
        
        double* abs_der_max_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        max_abs_derivative.resize(finest_level_dim_0);
        double* abs_der_max_global = max_abs_derivative.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            abs_der_max_local[i]  = double(0);
            abs_der_max_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                d_patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratio_to_coarest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarest_level/ratio_to_coarest_level;
            
            const int ratio_to_finest_level_0 = ratio_to_finest_level[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower indices.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model->getNumberOfGhostCells();
                TBOX_ASSERT(num_ghosts >= num_ghosts_der);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, num_ghosts_der));
                
                d_flow_model->registerDerivedVariables(num_subghosts_of_data);
                
                d_flow_model->allocateMemoryForDerivedCellData();
                
                d_flow_model->computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_quantity =
                    d_flow_model->getCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                
                /*
                 * Initialize cell data for the derivative and get pointer to the cell data.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_derivative(
                    new pdat::CellData<double>(patch_box, 1, hier::IntVector::getZero(d_dim)));
                
                double* der = data_derivative->getPointer(0);
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    if (derivative_direction == 0)
                    {
                        derivative_first_order_x->computeDerivative(
                            data_derivative,
                            data_quantity,
                            dx[0],
                            patch_visible_box,
                            0,
                            0);
                    }
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        /*
                         * Compute the linear index and update the max.
                         */
                        
                        const int idx = relative_idx_lo_0 + i;
                        
                        for (int ii = 0; ii < ratio_to_finest_level_0; ii++)
                        {
                            const int idx_fine = (idx_lo_0 + i)*ratio_to_finest_level_0 + ii;
                            
                            abs_der_max_local[idx_fine] = fmax(abs_der_max_local[idx_fine], fabs(der[idx]));
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global max.
         */
        
        d_mpi.Allreduce(
            abs_der_max_local,
            abs_der_max_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_MAX);
        
        std::free(abs_der_max_local);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_x(
            new DerivativeFirstOrder(
                "first order derivative in x-direction",
                d_dim,
                DIRECTION::X_DIRECTION,
                num_ghosts_derivative));
        
        HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_y(
            new DerivativeFirstOrder(
                "first order derivative in y-direction",
                d_dim,
                DIRECTION::Y_DIRECTION,
                num_ghosts_derivative));
        
        const int finest_level_dim_0 = d_finest_level_dims[0];
        
        double* abs_der_max_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        max_abs_derivative.resize(finest_level_dim_0);
        double* abs_der_max_global = max_abs_derivative.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            abs_der_max_local[i]  = double(0);
            abs_der_max_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                d_patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratio_to_coarest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarest_level/ratio_to_coarest_level;
            
            const int ratio_to_finest_level_0 = ratio_to_finest_level[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower indices and grid spacings.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model->getNumberOfGhostCells();
                TBOX_ASSERT(num_ghosts >= num_ghosts_der);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, num_ghosts_der));
                
                d_flow_model->registerDerivedVariables(num_subghosts_of_data);
                
                d_flow_model->allocateMemoryForDerivedCellData();
                
                d_flow_model->computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_quantity =
                    d_flow_model->getCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                
                /*
                 * Initialize cell data for the derivative and get pointer to the cell data.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_derivative(
                    new pdat::CellData<double>(patch_box, 1, hier::IntVector::getZero(d_dim)));
                
                double* der = data_derivative->getPointer(0);
                
                const hier::IntVector patch_interior_dims = patch_box.numberCells();
                const int patch_interior_dim_0 = patch_interior_dims[0];
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    if (derivative_direction == 0)
                    {
                        derivative_first_order_x->computeDerivative(
                            data_derivative,
                            data_quantity,
                            dx[0],
                            patch_visible_box,
                            0,
                            0);
                    }
                    else if (derivative_direction == 1)
                    {
                        derivative_first_order_y->computeDerivative(
                            data_derivative,
                            data_quantity,
                            dx[1],
                            patch_visible_box,
                            0,
                            0);
                    }
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            /*
                             * Compute the linear index and update the max.
                             */
                            
                            const int idx = (relative_idx_lo_0 + i) +
                                (relative_idx_lo_1 + j)*patch_interior_dim_0;
                            
                            for (int ii = 0; ii < ratio_to_finest_level_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratio_to_finest_level_0 + ii;
                                
                                abs_der_max_local[idx_fine] = fmax(abs_der_max_local[idx_fine], fabs(der[idx]));
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global max.
         */
        
        d_mpi.Allreduce(
            abs_der_max_local,
            abs_der_max_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_MAX);
        
        std::free(abs_der_max_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_x(
            new DerivativeFirstOrder(
                "first order derivative in x-direction",
                d_dim,
                DIRECTION::X_DIRECTION,
                num_ghosts_derivative));
        
        HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_y(
            new DerivativeFirstOrder(
                "first order derivative in y-direction",
                d_dim,
                DIRECTION::Y_DIRECTION,
                num_ghosts_derivative));
        
        HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_z(
            new DerivativeFirstOrder(
                "first order derivative in z-direction",
                d_dim,
                DIRECTION::Z_DIRECTION,
                num_ghosts_derivative));
        
        const int finest_level_dim_0 = d_finest_level_dims[0];
        
        double* abs_der_max_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        max_abs_derivative.resize(finest_level_dim_0);
        double* abs_der_max_global = max_abs_derivative.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            abs_der_max_local[i]  = double(0);
            abs_der_max_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                d_patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratio_to_coarest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarest_level/ratio_to_coarest_level;
            
            const int ratio_to_finest_level_0 = ratio_to_finest_level[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower indices and grid spacings.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model->getNumberOfGhostCells();
                TBOX_ASSERT(num_ghosts >= num_ghosts_der);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, num_ghosts));
                
                d_flow_model->registerDerivedVariables(num_subghosts_of_data);
                
                d_flow_model->allocateMemoryForDerivedCellData();
                
                d_flow_model->computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_quantity =
                    d_flow_model->getCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int num_ghosts_2_quantity = num_ghosts_quantity[2];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                const int ghostcell_dim_1_quantity = ghostcell_dims_quantity[1];
                
                /*
                 * Initialize cell data for the derivative and get pointer to the cell data.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_derivative(
                    new pdat::CellData<double>(patch_box, 1, hier::IntVector::getZero(d_dim)));
                
                double* der = data_derivative->getPointer(0);
                
                const hier::IntVector patch_interior_dims = patch_box.numberCells();
                const int patch_interior_dim_0 = patch_interior_dims[0];
                const int patch_interior_dim_1 = patch_interior_dims[1];
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    if (derivative_direction == 0)
                    {
                        derivative_first_order_x->computeDerivative(
                            data_derivative,
                            data_quantity,
                            dx[0],
                            patch_visible_box,
                            0,
                            0);
                    }
                    else if (derivative_direction == 1)
                    {
                        derivative_first_order_y->computeDerivative(
                            data_derivative,
                            data_quantity,
                            dx[1],
                            patch_visible_box,
                            0,
                            0);
                    }
                    else if (derivative_direction == 2)
                    {
                         derivative_first_order_z->computeDerivative(
                            data_derivative,
                            data_quantity,
                            dx[2],
                            patch_visible_box,
                            0,
                            0);
                    }
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    const int interior_dim_2 = interior_dims[2];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int idx_lo_2 = index_lo[2];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    const int relative_idx_lo_2 = relative_index_lo[2];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                /*
                                 * Compute the linear index and update the max.
                                 */
                                
                                const int idx = (relative_idx_lo_0 + i) +
                                    (relative_idx_lo_1 + j)*patch_interior_dim_0 +
                                    (relative_idx_lo_2 + k)*patch_interior_dim_0*
                                        patch_interior_dim_1;
                                
                                for (int ii = 0; ii < ratio_to_finest_level_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratio_to_finest_level_0 + ii;
                                    
                                    abs_der_max_local[idx_fine] = fmax(abs_der_max_local[idx_fine], fabs(der[idx]));
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global max.
         */
        
        d_mpi.Allreduce(
            abs_der_max_local,
            abs_der_max_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_MAX);
        
        std::free(abs_der_max_local);
    }
    
    return max_abs_derivative;
}



/*
 * Compute maximum value of magnitude of gradient with only x-direction as inhomogeneous direction.
 */
std::vector<double>
FlowModelMPIHelperMaxMin::getMaxMagnitudeGradientWithInhomogeneousXDirection(
    const std::string quantity_name,
    const int component_idx,
    const int num_ghosts_derivative,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
    hier::IntVector num_ghosts_der = hier::IntVector::getOne(d_dim)*num_ghosts_derivative;
    
    std::vector<double> max_mag_gradient;
    
    const int num_levels = d_patch_hierarchy->getNumberOfLevels();
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    HAMERS_SHARED_PTR<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *d_patch_hierarchy,
            0,
            num_levels - 1));
    
    if (d_dim == tbox::Dimension(1))
    {
        HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_x(
            new DerivativeFirstOrder(
                "first order derivative in x-direction",
                d_dim,
                DIRECTION::X_DIRECTION,
                num_ghosts_derivative));
        
        const int finest_level_dim_0 = d_finest_level_dims[0];
        
        double* mag_grad_max_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        max_mag_gradient.resize(finest_level_dim_0);
        double* mag_grad_max_global = max_mag_gradient.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            mag_grad_max_local[i]  = double(0);
            mag_grad_max_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                d_patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratio_to_coarest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarest_level/ratio_to_coarest_level;
            
            const int ratio_to_finest_level_0 = ratio_to_finest_level[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower indices.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model->getNumberOfGhostCells();
                TBOX_ASSERT(num_ghosts >= num_ghosts_der);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, num_ghosts_der));
                
                d_flow_model->registerDerivedVariables(num_subghosts_of_data);
                
                d_flow_model->allocateMemoryForDerivedCellData();
                
                d_flow_model->computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_quantity =
                    d_flow_model->getCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                
                /*
                 * Initialize cell data for the derivative and get pointer to the cell data.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_derivative(
                    new pdat::CellData<double>(patch_box, d_dim.getValue(), hier::IntVector::getZero(d_dim)));
                
                double* der_x = data_derivative->getPointer(0);
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    derivative_first_order_x->computeDerivative(
                        data_derivative,
                        data_quantity,
                        dx[0],
                        patch_visible_box,
                        0,
                        0);
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        /*
                         * Compute the linear index and update the max.
                         */
                        
                        const int idx = relative_idx_lo_0 + i;
                        
                        const double mag_grad = sqrt(der_x[idx]*der_x[idx]);
                        
                        for (int ii = 0; ii < ratio_to_finest_level_0; ii++)
                        {
                            const int idx_fine = (idx_lo_0 + i)*ratio_to_finest_level_0 + ii;
                            
                            mag_grad_max_local[idx_fine] = fmax(mag_grad_max_local[idx_fine], mag_grad);
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global max.
         */
        
        d_mpi.Allreduce(
            mag_grad_max_local,
            mag_grad_max_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_MAX);
        
        std::free(mag_grad_max_local);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_x(
            new DerivativeFirstOrder(
                "first order derivative in x-direction",
                d_dim,
                DIRECTION::X_DIRECTION,
                num_ghosts_derivative));
        
        HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_y(
            new DerivativeFirstOrder(
                "first order derivative in y-direction",
                d_dim,
                DIRECTION::Y_DIRECTION,
                num_ghosts_derivative));
        
        const int finest_level_dim_0 = d_finest_level_dims[0];
        
        double* mag_grad_max_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        max_mag_gradient.resize(finest_level_dim_0);
        double* mag_grad_max_global = max_mag_gradient.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            mag_grad_max_local[i]  = double(0);
            mag_grad_max_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                d_patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratio_to_coarest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarest_level/ratio_to_coarest_level;
            
            const int ratio_to_finest_level_0 = ratio_to_finest_level[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower indices and grid spacings.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model->getNumberOfGhostCells();
                TBOX_ASSERT(num_ghosts >= num_ghosts_der);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, num_ghosts_der));
                
                d_flow_model->registerDerivedVariables(num_subghosts_of_data);
                
                d_flow_model->allocateMemoryForDerivedCellData();
                
                d_flow_model->computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_quantity =
                    d_flow_model->getCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                
                /*
                 * Initialize cell data for the derivatives and get pointers to the cell data.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_derivatives(
                    new pdat::CellData<double>(patch_box, d_dim.getValue(), hier::IntVector::getZero(d_dim)));
                
                double* der_x = data_derivatives->getPointer(0);
                double* der_y = data_derivatives->getPointer(1);
                
                const hier::IntVector patch_interior_dims = patch_box.numberCells();
                const int patch_interior_dim_0 = patch_interior_dims[0];
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    derivative_first_order_x->computeDerivative(
                        data_derivatives,
                        data_quantity,
                        dx[0],
                        patch_visible_box,
                        0,
                        0);
                    
                    derivative_first_order_y->computeDerivative(
                        data_derivatives,
                        data_quantity,
                        dx[1],
                        patch_visible_box,
                        1,
                        0);
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            /*
                             * Compute the linear index and update the max.
                             */
                            
                            const int idx = (relative_idx_lo_0 + i) +
                                (relative_idx_lo_1 + j)*patch_interior_dim_0;
                            
                            const double mag_grad = sqrt(der_x[idx]*der_x[idx] + der_y[idx]*der_y[idx]);
                            
                            for (int ii = 0; ii < ratio_to_finest_level_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratio_to_finest_level_0 + ii;
                                
                                mag_grad_max_local[idx_fine] = fmax(mag_grad_max_local[idx_fine], mag_grad);
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global max.
         */
        
        d_mpi.Allreduce(
            mag_grad_max_local,
            mag_grad_max_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_MAX);
        
        std::free(mag_grad_max_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_x(
            new DerivativeFirstOrder(
                "first order derivative in x-direction",
                d_dim,
                DIRECTION::X_DIRECTION,
                num_ghosts_derivative));
        
        HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_y(
            new DerivativeFirstOrder(
                "first order derivative in y-direction",
                d_dim,
                DIRECTION::Y_DIRECTION,
                num_ghosts_derivative));
        
        HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_z(
            new DerivativeFirstOrder(
                "first order derivative in z-direction",
                d_dim,
                DIRECTION::Z_DIRECTION,
                num_ghosts_derivative));
        
        const int finest_level_dim_0 = d_finest_level_dims[0];
        
        double* mag_grad_max_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        max_mag_gradient.resize(finest_level_dim_0);
        double* mag_grad_max_global = max_mag_gradient.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            mag_grad_max_local[i]  = double(0);
            mag_grad_max_global[i] = double(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                d_patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratio_to_coarest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarest_level/ratio_to_coarest_level;
            
            const int ratio_to_finest_level_0 = ratio_to_finest_level[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower indices and grid spacings.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model->getNumberOfGhostCells();
                TBOX_ASSERT(num_ghosts >= num_ghosts_der);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, num_ghosts));
                
                d_flow_model->registerDerivedVariables(num_subghosts_of_data);
                
                d_flow_model->allocateMemoryForDerivedCellData();
                
                d_flow_model->computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_quantity =
                    d_flow_model->getCellData(quantity_name);
                
                double* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int num_ghosts_2_quantity = num_ghosts_quantity[2];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                const int ghostcell_dim_1_quantity = ghostcell_dims_quantity[1];
                
                /*
                 * Initialize cell data for the derivatives and get pointers to the cell data.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_derivatives(
                    new pdat::CellData<double>(patch_box, d_dim.getValue(), hier::IntVector::getZero(d_dim)));
                
                double* der_x = data_derivatives->getPointer(0);
                double* der_y = data_derivatives->getPointer(1);
                double* der_z = data_derivatives->getPointer(2);
                
                const hier::IntVector patch_interior_dims = patch_box.numberCells();
                const int patch_interior_dim_0 = patch_interior_dims[0];
                const int patch_interior_dim_1 = patch_interior_dims[1];
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    derivative_first_order_x->computeDerivative(
                        data_derivatives,
                        data_quantity,
                        dx[0],
                        patch_visible_box,
                        0,
                        0);
                    
                    derivative_first_order_y->computeDerivative(
                        data_derivatives,
                        data_quantity,
                        dx[1],
                        patch_visible_box,
                        1,
                        0);
                    
                     derivative_first_order_z->computeDerivative(
                        data_derivatives,
                        data_quantity,
                        dx[2],
                        patch_visible_box,
                        2,
                        0);
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    const int interior_dim_2 = interior_dims[2];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int idx_lo_2 = index_lo[2];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    const int relative_idx_lo_2 = relative_index_lo[2];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                /*
                                 * Compute the linear index and update the max.
                                 */
                                
                                const int idx = (relative_idx_lo_0 + i) +
                                    (relative_idx_lo_1 + j)*patch_interior_dim_0 +
                                    (relative_idx_lo_2 + k)*patch_interior_dim_0*
                                        patch_interior_dim_1;
                                
                                const double mag_grad = sqrt(der_x[idx]*der_x[idx] + der_y[idx]*der_y[idx] + der_z[idx]*der_z[idx]);
                                
                                for (int ii = 0; ii < ratio_to_finest_level_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratio_to_finest_level_0 + ii;
                                    
                                    mag_grad_max_local[idx_fine] = fmax(mag_grad_max_local[idx_fine], mag_grad);
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global max.
         */
        
        d_mpi.Allreduce(
            mag_grad_max_local,
            mag_grad_max_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_MAX);
        
        std::free(mag_grad_max_local);
    }
    
    return max_mag_gradient;
}
