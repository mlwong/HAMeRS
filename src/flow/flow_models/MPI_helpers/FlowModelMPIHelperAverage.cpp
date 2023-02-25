#include "flow/flow_models/MPI_helpers/FlowModelMPIHelperAverage.hpp"

#include "extn/patch_hierarchies/ExtendedFlattenedHierarchy.hpp"
#include "util/derivatives/DerivativeFirstOrder.hpp"

/*
 * Compute averaged value over the entire domain.
 */
Real
FlowModelMPIHelperAverage::getAveragedQuantity(
    const std::string quantity_name,
    const int component_idx,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    Real u_avg_local  = Real(0);
    Real u_avg_global = Real(0);
    
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
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_x = Real(x_hi[0] - x_lo[0]);
        
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
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity =
                    getCellData(quantity_name);
                
                Real* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                
                const Real weight = Real(dx[0])/L_x;
                
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
                         * Compute the index of the data point and count how many times the data is repeated.
                         */
                        
                        const hier::Index idx_pt(tbox::Dimension(1), idx_lo_0 + i);
                        
                        int n_overlapped = 1;
                        
                        for (hier::BoxContainer::BoxContainerConstIterator iob(
                                patch_overlapped_visible_boxes.begin());
                             iob != patch_overlapped_visible_boxes.end();
                             iob++)
                        {
                            const hier::Box& patch_overlapped_visible_box = *iob;
                            
                            if (patch_overlapped_visible_box.contains(idx_pt))
                            {
                                n_overlapped++;
                            }
                        }
                        
                        /*
                         * Compute the linear index and the data to add.
                         */
                        
                        const int idx = relative_idx_lo_0 + i + num_ghosts_0_quantity;
                        
                        u_avg_local += (u[idx]*weight/((Real) n_overlapped));
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        d_mpi.Allreduce(
            &u_avg_local,
            &u_avg_global,
            1,
            MPI_DOUBLE,
            MPI_SUM);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_x = Real(x_hi[0] - x_lo[0]);
        const Real L_y = Real(x_hi[1] - x_lo[1]);
        
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
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity =
                    getCellData(quantity_name);
                
                Real* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                
                const Real weight = Real(dx[0]*dx[1])/(L_x*L_y);
                
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
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j);
                            
                            int n_overlapped = 1;
                            
                            for (hier::BoxContainer::BoxContainerConstIterator iob(
                                    patch_overlapped_visible_boxes.begin());
                                 iob != patch_overlapped_visible_boxes.end();
                                 iob++)
                            {
                                const hier::Box& patch_overlapped_visible_box = *iob;
                                
                                if (patch_overlapped_visible_box.contains(idx_pt))
                                {
                                    n_overlapped++;
                                }
                            }
                            
                            /*
                             * Compute the linear index and the data to add.
                             */
                            
                            const int idx = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                            
                            u_avg_local += u[idx]*weight/((Real) n_overlapped);
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        d_mpi.Allreduce(
            &u_avg_local,
            &u_avg_global,
            1,
            MPI_DOUBLE,
            MPI_SUM);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_x = Real(x_hi[0] - x_lo[0]);
        const Real L_y = Real(x_hi[1] - x_lo[1]);
        const Real L_z = Real(x_hi[2] - x_lo[2]);
        
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
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity =
                    getCellData(quantity_name);
                
                Real* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int num_ghosts_2_quantity = num_ghosts_quantity[2];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                const int ghostcell_dim_1_quantity = ghostcell_dims_quantity[1];
                
                const Real weight = Real(dx[0]*dx[1]*dx[2])/(L_x*L_y*L_z);
                
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
                                 * Compute the index of the data point and count how many times the data is repeated.
                                 */
                                
                                const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                                
                                int n_overlapped = 1;
                                
                                for (hier::BoxContainer::BoxContainerConstIterator iob(
                                        patch_overlapped_visible_boxes.begin());
                                     iob != patch_overlapped_visible_boxes.end();
                                     iob++)
                                {
                                    const hier::Box& patch_overlapped_visible_box = *iob;
                                    
                                    if (patch_overlapped_visible_box.contains(idx_pt))
                                    {
                                        n_overlapped++;
                                    }
                                }
                                
                                /*
                                 * Compute the linear index and the data to add.
                                 */
                                
                                const int idx = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                    (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                        ghostcell_dim_1_quantity;
                                
                                u_avg_local += u[idx]*weight/((Real) n_overlapped);
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        d_mpi.Allreduce(
            &u_avg_local,
            &u_avg_global,
            1,
            MPI_DOUBLE,
            MPI_SUM);
    }
    
    return u_avg_global;
}


/*
 * Compute averaged reciprocal of value over the entire domain.
 */
Real
FlowModelMPIHelperAverage::getAveragedReciprocalOfQuantity(
    const std::string quantity_name,
    const int component_idx,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    Real u_inv_avg_local  = Real(0);
    Real u_inv_avg_global = Real(0);
    
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
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_x = Real(x_hi[0] - x_lo[0]);
        
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
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity =
                    getCellData(quantity_name);
                
                Real* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                
                const Real weight = Real(dx[0])/L_x;
                
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
                         * Compute the index of the data point and count how many times the data is repeated.
                         */
                        
                        const hier::Index idx_pt(tbox::Dimension(1), idx_lo_0 + i);
                        
                        int n_overlapped = 1;
                        
                        for (hier::BoxContainer::BoxContainerConstIterator iob(
                                patch_overlapped_visible_boxes.begin());
                             iob != patch_overlapped_visible_boxes.end();
                             iob++)
                        {
                            const hier::Box& patch_overlapped_visible_box = *iob;
                            
                            if (patch_overlapped_visible_box.contains(idx_pt))
                            {
                                n_overlapped++;
                            }
                        }
                        
                        /*
                         * Compute the linear index and the data to add.
                         */
                        
                        const int idx = relative_idx_lo_0 + i + num_ghosts_0_quantity;
                        
                        u_inv_avg_local += (weight/u[idx])/((Real) n_overlapped);
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average of reciprocal.
         */
        
        d_mpi.Allreduce(
            &u_inv_avg_local,
            &u_inv_avg_global,
            1,
            MPI_DOUBLE,
            MPI_SUM);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_x = Real(x_hi[0] - x_lo[0]);
        const Real L_y = Real(x_hi[1] - x_lo[1]);
        
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
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity =
                    getCellData(quantity_name);
                
                Real* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                
                const Real weight = Real(dx[0]*dx[1])/(L_x*L_y);
                
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
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j);
                            
                            int n_overlapped = 1;
                            
                            for (hier::BoxContainer::BoxContainerConstIterator iob(
                                    patch_overlapped_visible_boxes.begin());
                                 iob != patch_overlapped_visible_boxes.end();
                                 iob++)
                            {
                                const hier::Box& patch_overlapped_visible_box = *iob;
                                
                                if (patch_overlapped_visible_box.contains(idx_pt))
                                {
                                    n_overlapped++;
                                }
                            }
                            
                            /*
                             * Compute the linear index and the data to add.
                             */
                            
                            const int idx = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                            
                            u_inv_avg_local += (weight/u[idx])/((Real) n_overlapped);
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average of reciprocal.
         */
        
        d_mpi.Allreduce(
            &u_inv_avg_local,
            &u_inv_avg_global,
            1,
            MPI_DOUBLE,
            MPI_SUM);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_x = Real(x_hi[0] - x_lo[0]);
        const Real L_y = Real(x_hi[1] - x_lo[1]);
        const Real L_z = Real(x_hi[2] - x_lo[2]);
        
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
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity =
                    getCellData(quantity_name);
                
                Real* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int num_ghosts_2_quantity = num_ghosts_quantity[2];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                const int ghostcell_dim_1_quantity = ghostcell_dims_quantity[1];
                
                const Real weight = Real(dx[0]*dx[1]*dx[2])/(L_x*L_y*L_z);
                
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
                                 * Compute the index of the data point and count how many times the data is repeated.
                                 */
                                
                                const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                                
                                int n_overlapped = 1;
                                
                                for (hier::BoxContainer::BoxContainerConstIterator iob(
                                        patch_overlapped_visible_boxes.begin());
                                     iob != patch_overlapped_visible_boxes.end();
                                     iob++)
                                {
                                    const hier::Box& patch_overlapped_visible_box = *iob;
                                    
                                    if (patch_overlapped_visible_box.contains(idx_pt))
                                    {
                                        n_overlapped++;
                                    }
                                }
                                
                                /*
                                 * Compute the linear index and the data to add.
                                 */
                                
                                const int idx = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                    (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                        ghostcell_dim_1_quantity;
                                
                                u_inv_avg_local += (weight/u[idx])/((Real) n_overlapped);
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average of reciprocal.
         */
        
        d_mpi.Allreduce(
            &u_inv_avg_local,
            &u_inv_avg_global,
            1,
            MPI_DOUBLE,
            MPI_SUM);
    }
    
    return u_inv_avg_global;
}


/*
 * Compute averaged value (on product of variables) over the entire domain.
 */
Real
FlowModelMPIHelperAverage::getAveragedQuantity(
    const std::vector<std::string>& quantity_names,
    const std::vector<int>& component_indices,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    int num_quantities = static_cast<int>(quantity_names.size());
    
    TBOX_ASSERT(static_cast<int>(component_indices.size()) == num_quantities);
    
    Real avg_local  = Real(0);
    Real avg_global = Real(0);
    
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
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_x = Real(x_hi[0] - x_lo[0]);
        
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
                 * Get the patch lower indices.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and data in the flow model and compute the corresponding
                 * average.
                 */
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                }
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = getCellData(quantity_names[qi]);
                }
                
                std::vector<Real*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                }
                
                const Real weight = Real(dx[0])/L_x;
                
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
                         * Compute the index of the data point and count how many times the data is repeated.
                         */
                        
                        const hier::Index idx_pt(tbox::Dimension(1), idx_lo_0 + i);
                        
                        int n_overlapped = 1;
                        
                        for (hier::BoxContainer::BoxContainerConstIterator iob(
                                patch_overlapped_visible_boxes.begin());
                             iob != patch_overlapped_visible_boxes.end();
                             iob++)
                        {
                            const hier::Box& patch_overlapped_visible_box = *iob;
                            
                            if (patch_overlapped_visible_box.contains(idx_pt))
                            {
                                n_overlapped++;
                            }
                        }
                        
                        /*
                         * Compute the linear indices and the data to add.
                         */
                        
                        const int idx_q0 = relative_idx_lo_0 + i + num_ghosts_0_u_qi[0];
                        
                        Real avg = u_qi[0][idx_q0];
                        
                        for (int qi = 1; qi < num_quantities; qi++)
                        {
                            const int idx_qi = relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi];
                            
                            avg *= u_qi[qi][idx_qi];
                        }
                        
                        avg_local += (avg*weight/((Real) n_overlapped));
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        d_mpi.Allreduce(
            &avg_local,
            &avg_global,
            1,
            MPI_DOUBLE,
            MPI_SUM);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_x = Real(x_hi[0] - x_lo[0]);
        const Real L_y = Real(x_hi[1] - x_lo[1]);
        
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
                 * Get the patch lower indices and grid spacings.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and data in the flow model and compute the corresponding
                 * average.
                 */
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                }
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = getCellData(quantity_names[qi]);
                }
                
                std::vector<Real*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                }
                
                const Real weight = Real(dx[0]*dx[1])/(L_x*L_y);
                
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
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j);
                            
                            int n_overlapped = 1;
                            
                            for (hier::BoxContainer::BoxContainerConstIterator iob(
                                    patch_overlapped_visible_boxes.begin());
                                 iob != patch_overlapped_visible_boxes.end();
                                 iob++)
                            {
                                const hier::Box& patch_overlapped_visible_box = *iob;
                                
                                if (patch_overlapped_visible_box.contains(idx_pt))
                                {
                                    n_overlapped++;
                                }
                            }
                            
                            /*
                             * Compute the linear indices and the data to add.
                             */
                            
                            const int idx_q0 = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0];
                            
                            Real avg = u_qi[0][idx_q0];
                            
                            for (int qi = 1; qi < num_quantities; qi++)
                            {
                                const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                
                                avg *= u_qi[qi][idx_qi];
                            }
                            
                            avg_local += (avg*weight/((Real) n_overlapped));
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        d_mpi.Allreduce(
            &avg_local,
            &avg_global,
            1,
            MPI_DOUBLE,
            MPI_SUM);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_x = Real(x_hi[0] - x_lo[0]);
        const Real L_y = Real(x_hi[1] - x_lo[1]);
        const Real L_z = Real(x_hi[2] - x_lo[2]);
        
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
                 * Get the patch lower indices and grid spacings.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and data in the flow model and compute the corresponding
                 * average.
                 */
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                }
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = getCellData(quantity_names[qi]);
                }
                
                std::vector<Real*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> num_ghosts_2_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                std::vector<int> ghostcell_dim_1_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                num_ghosts_2_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                ghostcell_dim_1_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    num_ghosts_2_u_qi.push_back(num_ghosts_u_qi[2]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                    ghostcell_dim_1_u_qi.push_back(ghostcell_dims_u_qi[1]);
                }
                
                const Real weight = Real(dx[0]*dx[1]*dx[2])/(L_x*L_y*L_z);
                
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
                                 * Compute the index of the data point and count how many times the data is repeated.
                                 */
                                
                                const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                                
                                int n_overlapped = 1;
                                
                                for (hier::BoxContainer::BoxContainerConstIterator iob(
                                        patch_overlapped_visible_boxes.begin());
                                     iob != patch_overlapped_visible_boxes.end();
                                     iob++)
                                {
                                    const hier::Box& patch_overlapped_visible_box = *iob;
                                    
                                    if (patch_overlapped_visible_box.contains(idx_pt))
                                    {
                                        n_overlapped++;
                                    }
                                }
                                
                                /*
                                 * Compute the linear index and the data to add.
                                 */
                                
                                const int idx_q0 = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0] +
                                    (relative_idx_lo_2 + k + num_ghosts_2_u_qi[0])*ghostcell_dim_0_u_qi[0]*
                                        ghostcell_dim_1_u_qi[0];
                                
                                Real avg = u_qi[0][idx_q0];
                                
                                for (int qi = 1; qi < num_quantities; qi++)
                                {
                                    const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                        (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                            ghostcell_dim_1_u_qi[qi];
                                    
                                    avg *= u_qi[qi][idx_qi];
                                }
                                
                                avg_local += (avg*weight/((Real) n_overlapped));
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        d_mpi.Allreduce(
            &avg_local,
            &avg_global,
            1,
            MPI_DOUBLE,
            MPI_SUM);
    }
    
    return avg_global;
}


/*
 * Compute averaged value (on product of variables) over the entire domain.
 */
Real
FlowModelMPIHelperAverage::getAveragedQuantity(
    const std::vector<std::string>& quantity_names,
    const std::vector<int>& component_indices,
    const std::vector<bool>& use_reciprocal,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    int num_quantities = static_cast<int>(quantity_names.size());
    
    TBOX_ASSERT(static_cast<int>(component_indices.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(use_reciprocal.size()) == num_quantities);
    
    Real avg_local  = Real(0);
    Real avg_global = Real(0);
    
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
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_x = Real(x_hi[0] - x_lo[0]);
        
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
                 * Get the patch lower indices.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and data in the flow model and compute the corresponding
                 * average.
                 */
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                }
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = getCellData(quantity_names[qi]);
                }
                
                std::vector<Real*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                }
                
                const Real weight = Real(dx[0])/L_x;
                
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
                         * Compute the index of the data point and count how many times the data is repeated.
                         */
                        
                        const hier::Index idx_pt(tbox::Dimension(1), idx_lo_0 + i);
                        
                        int n_overlapped = 1;
                        
                        for (hier::BoxContainer::BoxContainerConstIterator iob(
                                patch_overlapped_visible_boxes.begin());
                             iob != patch_overlapped_visible_boxes.end();
                             iob++)
                        {
                            const hier::Box& patch_overlapped_visible_box = *iob;
                            
                            if (patch_overlapped_visible_box.contains(idx_pt))
                            {
                                n_overlapped++;
                            }
                        }
                        
                        /*
                         * Compute the linear indices and the data to add.
                         */
                        
                        Real avg = Real(1);
                        
                        for (int qi = 0; qi < num_quantities; qi++)
                        {
                            const int idx_qi = relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi];
                            
                            if (use_reciprocal[qi])
                            {
                                avg /= u_qi[qi][idx_qi];
                            }
                            else
                            {
                                avg *= u_qi[qi][idx_qi];
                            }
                        }
                        
                        avg_local += (avg*weight/((Real) n_overlapped));
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        d_mpi.Allreduce(
            &avg_local,
            &avg_global,
            1,
            MPI_DOUBLE,
            MPI_SUM);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_x = Real(x_hi[0] - x_lo[0]);
        const Real L_y = Real(x_hi[1] - x_lo[1]);
        
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
                 * Get the patch lower indices and grid spacings.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and data in the flow model and compute the corresponding
                 * average.
                 */
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                }
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = getCellData(quantity_names[qi]);
                }
                
                std::vector<Real*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                }
                
                const Real weight = Real(dx[0]*dx[1])/(L_x*L_y);
                
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
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j);
                            
                            int n_overlapped = 1;
                            
                            for (hier::BoxContainer::BoxContainerConstIterator iob(
                                    patch_overlapped_visible_boxes.begin());
                                 iob != patch_overlapped_visible_boxes.end();
                                 iob++)
                            {
                                const hier::Box& patch_overlapped_visible_box = *iob;
                                
                                if (patch_overlapped_visible_box.contains(idx_pt))
                                {
                                    n_overlapped++;
                                }
                            }
                            
                            /*
                             * Compute the linear indices and the data to add.
                             */
                            
                            Real avg = Real(1);
                            
                            for (int qi = 0; qi < num_quantities; qi++)
                            {
                                const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                
                                if (use_reciprocal[qi])
                                {
                                    avg /= u_qi[qi][idx_qi];
                                }
                                else
                                {
                                    avg *= u_qi[qi][idx_qi];
                                }
                            }
                            
                            avg_local += (avg*weight/((Real) n_overlapped));
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        d_mpi.Allreduce(
            &avg_local,
            &avg_global,
            1,
            MPI_DOUBLE,
            MPI_SUM);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_x = Real(x_hi[0] - x_lo[0]);
        const Real L_y = Real(x_hi[1] - x_lo[1]);
        const Real L_z = Real(x_hi[2] - x_lo[2]);
        
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
                 * Get the patch lower indices and grid spacings.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and data in the flow model and compute the corresponding
                 * average.
                 */
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                }
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = getCellData(quantity_names[qi]);
                }
                
                std::vector<Real*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> num_ghosts_2_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                std::vector<int> ghostcell_dim_1_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                num_ghosts_2_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                ghostcell_dim_1_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    num_ghosts_2_u_qi.push_back(num_ghosts_u_qi[2]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                    ghostcell_dim_1_u_qi.push_back(ghostcell_dims_u_qi[1]);
                }
                
                const Real weight = Real(dx[0]*dx[1]*dx[2])/(L_x*L_y*L_z);
                
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
                                 * Compute the index of the data point and count how many times the data is repeated.
                                 */
                                
                                const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                                
                                int n_overlapped = 1;
                                
                                for (hier::BoxContainer::BoxContainerConstIterator iob(
                                        patch_overlapped_visible_boxes.begin());
                                     iob != patch_overlapped_visible_boxes.end();
                                     iob++)
                                {
                                    const hier::Box& patch_overlapped_visible_box = *iob;
                                    
                                    if (patch_overlapped_visible_box.contains(idx_pt))
                                    {
                                        n_overlapped++;
                                    }
                                }
                                
                                /*
                                 * Compute the linear index and the data to add.
                                 */
                                
                                Real avg = Real(1);
                                
                                for (int qi = 0; qi < num_quantities; qi++)
                                {
                                    const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                        (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                            ghostcell_dim_1_u_qi[qi];
                                    
                                    if (use_reciprocal[qi])
                                    {
                                        avg /= u_qi[qi][idx_qi];
                                    }
                                    else
                                    {
                                        avg *= u_qi[qi][idx_qi];
                                    }
                                }
                                
                                avg_local += (avg*weight/((Real) n_overlapped));
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        d_mpi.Allreduce(
            &avg_local,
            &avg_global,
            1,
            MPI_DOUBLE,
            MPI_SUM);
    }
    
    return avg_global;
}


/*
 * Compute averaged value (on product of variable derivatives) over the entire domain.
 */
Real
FlowModelMPIHelperAverage::getAveragedQuantity(
    const std::vector<std::string>& quantity_names,
    const std::vector<int>& component_indices,
    const std::vector<bool>& use_derivative,
    const std::vector<int>& derivative_directions,
    const int num_ghosts_derivative,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    int num_quantities = static_cast<int>(quantity_names.size());
    
    TBOX_ASSERT(static_cast<int>(component_indices.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(use_derivative.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(derivative_directions.size()) == num_quantities);
    
    std::vector<bool> use_reciprocal(num_quantities, false);
    
    return getAveragedQuantity(
        quantity_names,
        component_indices,
        use_derivative,
        derivative_directions,
        use_reciprocal,
        num_ghosts_derivative,
        data_context);
}


/*
 * Compute averaged value (on product of variable derivatives) over the entire domain.
 */
Real
FlowModelMPIHelperAverage::getAveragedQuantity(
    const std::vector<std::string>& quantity_names,
    const std::vector<int>& component_indices,
    const std::vector<bool>& use_derivative,
    const std::vector<int>& derivative_directions,
    const std::vector<bool>& use_reciprocal,
    const int num_ghosts_derivative,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    int num_quantities = static_cast<int>(quantity_names.size());
    
    TBOX_ASSERT(static_cast<int>(component_indices.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(use_derivative.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(derivative_directions.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(use_reciprocal.size()) == num_quantities);
    
    Real avg_local  = Real(0);
    Real avg_global = Real(0);
    
    int num_use_derivative = 0;
    
    if (d_dim == tbox::Dimension(1))
    {
        for (int qi = 0; qi < num_quantities; qi++)
        {
            if (use_derivative[qi])
            {
                num_use_derivative++;
                if (derivative_directions[qi] < 0 || derivative_directions[qi] > 0)
                {
                    TBOX_ERROR(d_object_name
                        << ": FlowModelMPIHelperAverage::getAveragedQuantity():\n"
                        << "Cannot take derivative for one-dimensional problem!\n"
                        << "derivative_directions[" << qi << "] = " << derivative_directions[qi] << " given!\n"
                        << std::endl);
                }
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        for (int qi = 0; qi < num_quantities; qi++)
        {
            if (use_derivative[qi])
            {
                num_use_derivative++;
                if (derivative_directions[qi] < 0 || derivative_directions[qi] > 1)
                {
                    TBOX_ERROR(d_object_name
                        << ": FlowModelMPIHelperAverage::getAveragedQuantity():\n"
                        << "Cannot take derivative for two-dimensional problem!\n"
                        << "derivative_directions[" << qi << "] = " << derivative_directions[qi] << " given!\n"
                        << std::endl);
                }
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        for (int qi = 0; qi < num_quantities; qi++)
        {
            if (use_derivative[qi])
            {
                num_use_derivative++;
                if (derivative_directions[qi] < 0 || derivative_directions[qi] > 2)
                {
                    TBOX_ERROR(d_object_name
                        << ": FlowModelMPIHelperAverage::getAveragedQuantity():\n"
                        << "Cannot take derivative for three-dimensional problem!\n"
                        << "derivative_directions[" << qi << "] = " << derivative_directions[qi] << " given!\n"
                        << std::endl);
                }
            }
        }
    }
    
    std::vector<Real> averaged_quantity;
    
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
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_x(
            new DerivativeFirstOrder(
                "first order derivative in x-direction",
                d_dim,
                DIRECTION::X_DIRECTION,
                num_ghosts_derivative));
        
        hier::IntVector num_ghosts_der = hier::IntVector::getOne(d_dim)*num_ghosts_derivative;
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_x = Real(x_hi[0] - x_lo[0]);
        
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
                 * Get the patch lower indices and grid spacings.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and data in the flow model and compute the corresponding
                 * average.
                 */
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    if (num_use_derivative > 0)
                    {
                        num_subghosts_of_data.insert(
                            std::pair<std::string, hier::IntVector>(quantity_names[qi], num_ghosts_der));
                    }
                    else
                    {
                        num_subghosts_of_data.insert(
                            std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                    }
                }
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = getCellData(quantity_names[qi]);
                }
                
                std::vector<Real*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                }
                
                /*
                 * Initialize cell data for the derivatives and get pointers to the cell data.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_derivative;
                std::vector<Real*> der_qi;
                
                if (num_use_derivative > 0)
                {
                    data_derivative = HAMERS_MAKE_SHARED<pdat::CellData<Real> >(patch_box, num_use_derivative, hier::IntVector::getZero(d_dim));
                    
                    der_qi.resize(num_use_derivative);
                    for (int qi = 0; qi < num_use_derivative; qi++)
                    {
                        der_qi[qi] = data_derivative->getPointer(qi);
                    }
                }
                
                const Real weight = Real(dx[0])/L_x;
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    int count_derivative = 0;
                    for (int qi = 0; qi < num_quantities; qi++)
                    {
                        if (use_derivative[qi] && derivative_directions[qi] == 0)
                        {
                            derivative_first_order_x->computeDerivative(
                                data_derivative,
                                data_quantities[qi],
                                Real(dx[0]),
                                patch_visible_box,
                                count_derivative,
                                component_indices[qi]);
                            
                            count_derivative++;
                        }
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
                         * Compute the index of the data point and count how many times the data is repeated.
                         */
                        
                        const hier::Index idx_pt(tbox::Dimension(1), idx_lo_0 + i);
                        
                        int n_overlapped = 1;
                        
                        for (hier::BoxContainer::BoxContainerConstIterator iob(
                                patch_overlapped_visible_boxes.begin());
                             iob != patch_overlapped_visible_boxes.end();
                             iob++)
                        {
                            const hier::Box& patch_overlapped_visible_box = *iob;
                            
                            if (patch_overlapped_visible_box.contains(idx_pt))
                            {
                                n_overlapped++;
                            }
                        }
                        
                        /*
                         * Compute the linear indices and the data to add.
                         */
                        
                        Real avg = Real(1);
                        
                        count_derivative = 0;
                        for (int qi = 0; qi < num_quantities; qi++)
                        {
                            if (use_reciprocal[qi])
                            {
                                if (use_derivative[qi])
                                {
                                    const int idx_der = relative_idx_lo_0 + i;
                                    
                                    avg /= der_qi[count_derivative][idx_der];
                                    count_derivative++;
                                }
                                else
                                {
                                    const int idx_qi = relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi];
                                    
                                    avg /= u_qi[qi][idx_qi];
                                }
                            }
                            else
                            {
                                if (use_derivative[qi])
                                {
                                    const int idx_der = relative_idx_lo_0 + i;
                                    
                                    avg *= der_qi[count_derivative][idx_der];
                                    count_derivative++;
                                }
                                else
                                {
                                    const int idx_qi = relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi];
                                    
                                    avg *= u_qi[qi][idx_qi];
                                }
                            }
                        }
                        
                        /*
                         * Add the data.
                         */
                        
                        avg_local += (avg*weight/((Real) n_overlapped));
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        d_mpi.Allreduce(
            &avg_local,
            &avg_global,
            1,
            MPI_DOUBLE,
            MPI_SUM);
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
        
        hier::IntVector num_ghosts_der = hier::IntVector::getOne(d_dim)*num_ghosts_derivative;
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_x = Real(x_hi[0] - x_lo[0]);
        const Real L_y = Real(x_hi[1] - x_lo[1]);
        
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
                 * Get the patch lower indices and grid spacings.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and data in the flow model and compute the corresponding
                 * average.
                 */
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    if (num_use_derivative > 0)
                    {
                        num_subghosts_of_data.insert(
                            std::pair<std::string, hier::IntVector>(quantity_names[qi], num_ghosts_der));
                    }
                    else
                    {
                        num_subghosts_of_data.insert(
                            std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                    }
                }
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = getCellData(quantity_names[qi]);
                }
                
                std::vector<Real*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                }
                
                /*
                 * Initialize cell data for the derivatives and get pointers to the cell data.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_derivative;
                std::vector<Real*> der_qi;
                
                if (num_use_derivative > 0)
                {
                    data_derivative = HAMERS_MAKE_SHARED<pdat::CellData<Real> >(patch_box, num_use_derivative, hier::IntVector::getZero(d_dim));
                    
                    der_qi.resize(num_use_derivative);
                    for (int qi = 0; qi < num_use_derivative; qi++)
                    {
                        der_qi[qi] = data_derivative->getPointer(qi);
                    }
                }
                
                const hier::IntVector patch_interior_dims = patch_box.numberCells();
                const int patch_interior_dim_0 = patch_interior_dims[0];
                
                const Real weight = Real(dx[0]*dx[1])/(L_x*L_y);
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    int count_derivative = 0;
                    for (int qi = 0; qi < num_quantities; qi++)
                    {
                        if (use_derivative[qi] && derivative_directions[qi] == 0)
                        {
                            derivative_first_order_x->computeDerivative(
                                data_derivative,
                                data_quantities[qi],
                                Real(dx[0]),
                                patch_visible_box,
                                count_derivative,
                                component_indices[qi]);
                            
                            count_derivative++;
                        }
                        else if (use_derivative[qi] && derivative_directions[qi] == 1)
                        {
                            derivative_first_order_y->computeDerivative(
                                data_derivative,
                                data_quantities[qi],
                                Real(dx[1]),
                                patch_visible_box,
                                count_derivative,
                                component_indices[qi]);
                            
                            count_derivative++;
                        }
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
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j);
                            
                            int n_overlapped = 1;
                            
                            for (hier::BoxContainer::BoxContainerConstIterator iob(
                                    patch_overlapped_visible_boxes.begin());
                                 iob != patch_overlapped_visible_boxes.end();
                                 iob++)
                            {
                                const hier::Box& patch_overlapped_visible_box = *iob;
                                
                                if (patch_overlapped_visible_box.contains(idx_pt))
                                {
                                    n_overlapped++;
                                }
                            }
                            
                            /*
                             * Compute the linear indices and the data to add.
                             */
                            
                            Real avg = Real(1);
                            
                            count_derivative = 0;
                            for (int qi = 0; qi < num_quantities; qi++)
                            {
                                if (use_reciprocal[qi])
                                {
                                    if (use_derivative[qi])
                                    {
                                        const int idx_der = (relative_idx_lo_0 + i) +
                                            (relative_idx_lo_1 + j)*patch_interior_dim_0;
                                        
                                        avg /= der_qi[count_derivative][idx_der];
                                        
                                        count_derivative++;
                                    }
                                    else
                                    {
                                        const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                        
                                        avg /= u_qi[qi][idx_qi];
                                    }
                                }
                                else
                                {
                                    if (use_derivative[qi])
                                    {
                                        const int idx_der = (relative_idx_lo_0 + i) +
                                            (relative_idx_lo_1 + j)*patch_interior_dim_0;
                                        
                                        avg *= der_qi[count_derivative][idx_der];
                                        
                                        count_derivative++;
                                    }
                                    else
                                    {
                                        const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                        
                                        avg *= u_qi[qi][idx_qi];
                                    }
                                }
                            }
                            
                            /*
                             * Add the data.
                             */
                            
                            avg_local += (avg*weight/((Real) n_overlapped));
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        d_mpi.Allreduce(
            &avg_local,
            &avg_global,
            1,
            MPI_DOUBLE,
            MPI_SUM);
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
        
        hier::IntVector num_ghosts_der = hier::IntVector::getOne(d_dim)*num_ghosts_derivative;
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_x = Real(x_hi[0] - x_lo[0]);
        const Real L_y = Real(x_hi[1] - x_lo[1]);
        const Real L_z = Real(x_hi[2] - x_lo[2]);
        
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
                 * Get the patch lower indices and grid spacings.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and data in the flow model and compute the corresponding
                 * average.
                 */
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    if (num_use_derivative > 0)
                    {
                        num_subghosts_of_data.insert(
                            std::pair<std::string, hier::IntVector>(quantity_names[qi], num_ghosts_der));
                    }
                    else
                    {
                        num_subghosts_of_data.insert(
                            std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                    }
                }
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = getCellData(quantity_names[qi]);
                }
                
                std::vector<Real*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> num_ghosts_2_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                std::vector<int> ghostcell_dim_1_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                num_ghosts_2_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                ghostcell_dim_1_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    num_ghosts_2_u_qi.push_back(num_ghosts_u_qi[2]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                    ghostcell_dim_1_u_qi.push_back(ghostcell_dims_u_qi[1]);
                }
                
                /*
                 * Initialize cell data for the derivatives and get pointers to the cell data.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_derivative;
                std::vector<Real*> der_qi;
                
                if (num_use_derivative > 0)
                {
                    data_derivative = HAMERS_MAKE_SHARED<pdat::CellData<Real> >(patch_box, num_use_derivative, hier::IntVector::getZero(d_dim));
                    
                    der_qi.resize(num_use_derivative);
                    for (int qi = 0; qi < num_use_derivative; qi++)
                    {
                        der_qi[qi] = data_derivative->getPointer(qi);
                    }
                }
                
                const hier::IntVector patch_interior_dims = patch_box.numberCells();
                const int patch_interior_dim_0 = patch_interior_dims[0];
                const int patch_interior_dim_1 = patch_interior_dims[1];
                
                const Real weight = Real(dx[0]*dx[1]*dx[2])/(L_x*L_y*L_z);
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    int count_derivative = 0;
                    for (int qi = 0; qi < num_quantities; qi++)
                    {
                        if (use_derivative[qi] && derivative_directions[qi] == 0)
                        {
                            derivative_first_order_x->computeDerivative(
                                data_derivative,
                                data_quantities[qi],
                                Real(dx[0]),
                                patch_visible_box,
                                count_derivative,
                                component_indices[qi]);
                            
                            count_derivative++;
                        }
                        else if (use_derivative[qi] && derivative_directions[qi] == 1)
                        {
                            derivative_first_order_y->computeDerivative(
                                data_derivative,
                                data_quantities[qi],
                                Real(dx[1]),
                                patch_visible_box,
                                count_derivative,
                                component_indices[qi]);
                            
                            count_derivative++;
                        }
                        else if (use_derivative[qi] && derivative_directions[qi] == 2)
                        {
                            derivative_first_order_z->computeDerivative(
                                data_derivative,
                                data_quantities[qi],
                                Real(dx[2]),
                                patch_visible_box,
                                count_derivative,
                                component_indices[qi]);
                            
                            count_derivative++;
                        }
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
                                 * Compute the index of the data point and count how many times the data is repeated.
                                 */
                                
                                const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                                
                                int n_overlapped = 1;
                                
                                for (hier::BoxContainer::BoxContainerConstIterator iob(
                                        patch_overlapped_visible_boxes.begin());
                                     iob != patch_overlapped_visible_boxes.end();
                                     iob++)
                                {
                                    const hier::Box& patch_overlapped_visible_box = *iob;
                                    
                                    if (patch_overlapped_visible_box.contains(idx_pt))
                                    {
                                        n_overlapped++;
                                    }
                                }
                                
                                /*
                                 * Compute the linear index and the data to add.
                                 */
                                
                                Real avg = Real(1);
                                
                                count_derivative = 0;
                                for (int qi = 0; qi < num_quantities; qi++)
                                {
                                    if (use_reciprocal[qi])
                                    {
                                        if (use_derivative[qi])
                                        {
                                            const int idx_der = (relative_idx_lo_0 + i) +
                                                (relative_idx_lo_1 + j)*patch_interior_dim_0 +
                                                (relative_idx_lo_2 + k)*patch_interior_dim_0*
                                                    patch_interior_dim_1;
                                            
                                            avg /= der_qi[count_derivative][idx_der];
                                            
                                            count_derivative++;
                                        }
                                        else
                                        {
                                            const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                                (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                                (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                    ghostcell_dim_1_u_qi[qi];
                                            
                                            avg /= u_qi[qi][idx_qi];
                                        }
                                    }
                                    else
                                    {
                                        if (use_derivative[qi])
                                        {
                                            const int idx_der = (relative_idx_lo_0 + i) +
                                                (relative_idx_lo_1 + j)*patch_interior_dim_0 +
                                                (relative_idx_lo_2 + k)*patch_interior_dim_0*
                                                    patch_interior_dim_1;
                                            
                                            avg *= der_qi[count_derivative][idx_der];
                                            
                                            count_derivative++;
                                        }
                                        else
                                        {
                                            const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                                (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                                (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                    ghostcell_dim_1_u_qi[qi];
                                            
                                            avg *= u_qi[qi][idx_qi];
                                        }
                                    }
                                }
                                
                                /*
                                 * Add the data.
                                 */
                                
                                avg_local += (avg*weight/((Real) n_overlapped));
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        d_mpi.Allreduce(
            &avg_local,
            &avg_global,
            1,
            MPI_DOUBLE,
            MPI_SUM);
    }
    
    return avg_global;
}


/*
 * Compute averaged derivative of value (on product of variables) over the entire domain.
 */
Real
FlowModelMPIHelperAverage::getAveragedDerivativeOfQuantity(
    const std::vector<std::string>& quantity_names,
    const std::vector<int>& component_indices,
    const std::vector<bool>& use_reciprocal,
    const int derivative_direction,
    const int num_ghosts_derivative,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    int num_quantities = static_cast<int>(quantity_names.size());
    
    TBOX_ASSERT(static_cast<int>(component_indices.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(use_reciprocal.size()) == num_quantities);
    
    Real der_avg_local  = Real(0);
    Real der_avg_global = Real(0);
    
    if (d_dim == tbox::Dimension(1) && (derivative_direction < 0 || derivative_direction > 0))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperAverage::getAveragedDerivativeOfQuantity():\n"
            << "Cannot take derivative for one-dimensional problem!\n"
            << "derivative_direction = " << derivative_direction << " given!\n"
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2) && (derivative_direction < 0 || derivative_direction > 1))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperAverage::getAveragedDerivativeOfQuantity():\n"
            << "Cannot take derivative for two-dimensional problem!\n"
            << "derivative_direction = " << derivative_direction << " given!\n"
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(3) && (derivative_direction < 0 || derivative_direction > 2))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperAverage::getAveragedDerivativeOfQuantity():\n"
            << "Cannot take derivative for three-dimensional problem!\n"
            << "derivative_direction = " << derivative_direction << " given!\n"
            << std::endl);
    }
    
    hier::IntVector num_ghosts_der = hier::IntVector::getZero(d_dim);
    num_ghosts_der[derivative_direction] = num_ghosts_derivative;
    
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
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_x(
            new DerivativeFirstOrder(
                "first order derivative in x-direction",
                d_dim,
                DIRECTION::X_DIRECTION,
                num_ghosts_derivative));
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_x = Real(x_hi[0] - x_lo[0]);
        
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
                 * Get the patch lower indices and grid spacings.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantities in the flow model and compute the
                 * corresponding product.
                 */
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model->getNumberOfGhostCells();
                TBOX_ASSERT(num_ghosts >= num_ghosts_der);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], num_ghosts_der));
                }
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = getCellData(quantity_names[qi]);
                }
                
                std::vector<Real*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                }
                
                /*
                 * Initialize cell data for the derivative and the product of variables.
                 * Also, get pointers to the cell data containers.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_derivative(
                    new pdat::CellData<Real>(patch_box, 1, hier::IntVector::getZero(d_dim)));
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_product(
                    new pdat::CellData<Real>(patch_box, 1, num_ghosts_der));
                
                Real* der     = data_derivative->getPointer(0);
                Real* product = data_product->getPointer(0);
                
                const hier::IntVector num_ghosts_product = data_product->getGhostCellWidth();
                
                const int num_ghosts_0_product = num_ghosts_product[0];
                
                const hier::IntVector patch_interior_dims = patch_box.numberCells();
                const int patch_interior_dim_0 = patch_interior_dims[0];
                
                data_product->fillAll(Real(1));
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    if (use_reciprocal[qi])
                    {
                        for (int i = -num_ghosts_0_product; i < patch_interior_dim_0 + num_ghosts_0_product; i++)
                        {
                            const int idx_product = i + num_ghosts_0_product;
                            const int idx_u_qi    = i + num_ghosts_0_u_qi[qi];
                            
                            product[idx_product] /= u_qi[qi][idx_u_qi];
                        }
                    }
                    else
                    {
                        for (int i = -num_ghosts_0_product; i < patch_interior_dim_0 + num_ghosts_0_product; i++)
                        {
                            const int idx_product = i + num_ghosts_0_product;
                            const int idx_u_qi    = i + num_ghosts_0_u_qi[qi];
                            
                            product[idx_product] *= u_qi[qi][idx_u_qi];
                        }
                    }
                }
                
                const Real weight = Real(dx[0])/L_x;
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    if (derivative_direction == 0)
                    {
                        derivative_first_order_x->computeDerivative(
                            data_derivative,
                            data_product,
                            Real(dx[0]),
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
                         * Compute the index of the data point and count how many times the data is repeated.
                         */
                        
                        const hier::Index idx_pt(tbox::Dimension(1), idx_lo_0 + i);
                        
                        int n_overlapped = 1;
                        
                        for (hier::BoxContainer::BoxContainerConstIterator iob(
                                patch_overlapped_visible_boxes.begin());
                             iob != patch_overlapped_visible_boxes.end();
                             iob++)
                        {
                            const hier::Box& patch_overlapped_visible_box = *iob;
                            
                            if (patch_overlapped_visible_box.contains(idx_pt))
                            {
                                n_overlapped++;
                            }
                        }
                        
                        /*
                         * Compute the linear index and add the data.
                         */
                        
                        const int idx = relative_idx_lo_0 + i;
                        
                        der_avg_local += (der[idx]*weight)/((Real) n_overlapped);
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average of derivative.
         */
        
        d_mpi.Allreduce(
            &der_avg_local,
            &der_avg_global,
            1,
            MPI_DOUBLE,
            MPI_SUM);
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
        
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_x = Real(x_hi[0] - x_lo[0]);
        const Real L_y = Real(x_hi[1] - x_lo[1]);
        
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
                 * Get the patch lower indices and grid spacings.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantities in the flow model and compute the
                 * corresponding product.
                 */
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model->getNumberOfGhostCells();
                TBOX_ASSERT(num_ghosts >= num_ghosts_der);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], num_ghosts_der));
                }
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = getCellData(quantity_names[qi]);
                }
                
                std::vector<Real*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                }
                
                /*
                 * Initialize cell data for the derivative and the product of variables.
                 * Also, get pointers to the cell data containers.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_derivative(
                    new pdat::CellData<Real>(patch_box, 1, hier::IntVector::getZero(d_dim)));
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_product(
                    new pdat::CellData<Real>(patch_box, 1, num_ghosts_der));
                
                Real* der     = data_derivative->getPointer(0);
                Real* product = data_product->getPointer(0);
                
                const hier::IntVector num_ghosts_product = data_product->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_product = data_product->getGhostBox().numberCells();
                
                const int num_ghosts_0_product = num_ghosts_product[0];
                const int num_ghosts_1_product = num_ghosts_product[1];
                const int ghostcell_dim_0_product = ghostcell_dims_product[0];
                
                const hier::IntVector patch_interior_dims = patch_box.numberCells();
                const int patch_interior_dim_0 = patch_interior_dims[0];
                const int patch_interior_dim_1 = patch_interior_dims[1];
                
                data_product->fillAll(Real(1));
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    if (use_reciprocal[qi])
                    {
                        for (int j = -num_ghosts_1_product; j < patch_interior_dim_1 + num_ghosts_1_product; j++)
                        {
                            for (int i = -num_ghosts_0_product; i < patch_interior_dim_0 + num_ghosts_0_product; i++)
                            {
                                const int idx_product = (i + num_ghosts_0_product) +
                                    (j + num_ghosts_1_product)*ghostcell_dim_0_product;
                                
                                const int idx_u_qi = (i + num_ghosts_0_u_qi[qi]) +
                                    (j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                product[idx_product] /= u_qi[qi][idx_u_qi];
                            }
                        }
                    }
                    else
                    {
                        for (int j = -num_ghosts_1_product; j < patch_interior_dim_1 + num_ghosts_1_product; j++)
                        {
                            for (int i = -num_ghosts_0_product; i < patch_interior_dim_0 + num_ghosts_0_product; i++)
                            {
                                const int idx_product = (i + num_ghosts_0_product) +
                                    (j + num_ghosts_1_product)*ghostcell_dim_0_product;
                                
                                const int idx_u_qi = (i + num_ghosts_0_u_qi[qi]) +
                                    (j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                product[idx_product] *= u_qi[qi][idx_u_qi];
                            }
                        }
                    }
                }
                
                const Real weight = Real(dx[0]*dx[1])/(L_x*L_y);
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    if (derivative_direction == 0)
                    {
                        derivative_first_order_x->computeDerivative(
                            data_derivative,
                            data_product,
                            Real(dx[0]),
                            patch_visible_box,
                            0,
                            0);
                    }
                    else if (derivative_direction == 1)
                    {
                        derivative_first_order_y->computeDerivative(
                            data_derivative,
                            data_product,
                            Real(dx[1]),
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
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j);
                            
                            int n_overlapped = 1;
                            
                            for (hier::BoxContainer::BoxContainerConstIterator iob(
                                    patch_overlapped_visible_boxes.begin());
                                 iob != patch_overlapped_visible_boxes.end();
                                 iob++)
                            {
                                const hier::Box& patch_overlapped_visible_box = *iob;
                                
                                if (patch_overlapped_visible_box.contains(idx_pt))
                                {
                                    n_overlapped++;
                                }
                            }
                            
                            /*
                             * Compute the linear index and add the data.
                             */
                            
                            const int idx = (relative_idx_lo_0 + i) +
                                (relative_idx_lo_1 + j)*patch_interior_dim_0;
                            
                            der_avg_local += (der[idx]*weight)/((Real) n_overlapped);
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average of derivative.
         */
        
        d_mpi.Allreduce(
            &der_avg_local,
            &der_avg_global,
            1,
            MPI_DOUBLE,
            MPI_SUM);
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
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_x = Real(x_hi[0] - x_lo[0]);
        const Real L_y = Real(x_hi[1] - x_lo[1]);
        const Real L_z = Real(x_hi[2] - x_lo[2]);
        
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
                 * Get the patch lower indices and grid spacings.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and the quantities in the flow model and compute the
                 * corresponding product.
                 */
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model->getNumberOfGhostCells();
                TBOX_ASSERT(num_ghosts >= num_ghosts_der);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], num_ghosts));
                }
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = getCellData(quantity_names[qi]);
                }
                
                std::vector<Real*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> num_ghosts_2_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                std::vector<int> ghostcell_dim_1_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                num_ghosts_2_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                ghostcell_dim_1_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    num_ghosts_2_u_qi.push_back(num_ghosts_u_qi[2]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                    ghostcell_dim_1_u_qi.push_back(ghostcell_dims_u_qi[1]);
                }
                
                /*
                 * Initialize cell data for the derivative and the product of variables.
                 * Also, get pointers to the cell data containers.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_derivative(
                    new pdat::CellData<Real>(patch_box, 1, hier::IntVector::getZero(d_dim)));
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_product(
                    new pdat::CellData<Real>(patch_box, 1, num_ghosts_der));
                
                Real* der     = data_derivative->getPointer(0);
                Real* product = data_product->getPointer(0);
                
                const hier::IntVector num_ghosts_product = data_product->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_product = data_product->getGhostBox().numberCells();
                
                const int num_ghosts_0_product = num_ghosts_product[0];
                const int num_ghosts_1_product = num_ghosts_product[1];
                const int num_ghosts_2_product = num_ghosts_product[2];
                const int ghostcell_dim_0_product = ghostcell_dims_product[0];
                const int ghostcell_dim_1_product = ghostcell_dims_product[1];
                
                const hier::IntVector patch_interior_dims = patch_box.numberCells();
                const int patch_interior_dim_0 = patch_interior_dims[0];
                const int patch_interior_dim_1 = patch_interior_dims[1];
                const int patch_interior_dim_2 = patch_interior_dims[2];
                
                data_product->fillAll(Real(1));
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    if (use_reciprocal[qi])
                    {
                        for (int k = -num_ghosts_2_product; k < patch_interior_dim_2 + num_ghosts_2_product; k++)
                        {
                            for (int j = -num_ghosts_1_product; j < patch_interior_dim_1 + num_ghosts_1_product; j++)
                            {
                                for (int i = -num_ghosts_0_product; i < patch_interior_dim_0 + num_ghosts_0_product; i++)
                                {
                                    const int idx_product = (i + num_ghosts_0_product) +
                                        (j + num_ghosts_1_product)*ghostcell_dim_0_product +
                                        (k + num_ghosts_2_product)*ghostcell_dim_0_product*
                                            ghostcell_dim_1_product;
                                    
                                    const int idx_u_qi = (i + num_ghosts_0_u_qi[qi]) +
                                        (j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                        (k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                            ghostcell_dim_1_u_qi[qi];
                                        
                                    product[idx_product] /= u_qi[qi][idx_u_qi];
                                }
                            }
                        }
                    }
                    else
                    {
                        for (int k = -num_ghosts_2_product; k < patch_interior_dim_2 + num_ghosts_2_product; k++)
                        {
                            for (int j = -num_ghosts_1_product; j < patch_interior_dim_1 + num_ghosts_1_product; j++)
                            {
                                for (int i = -num_ghosts_0_product; i < patch_interior_dim_0 + num_ghosts_0_product; i++)
                                {
                                    const int idx_product = (i + num_ghosts_0_product) +
                                        (j + num_ghosts_1_product)*ghostcell_dim_0_product +
                                        (k + num_ghosts_2_product)*ghostcell_dim_0_product*
                                            ghostcell_dim_1_product;
                                    
                                    const int idx_u_qi = (i + num_ghosts_0_u_qi[qi]) +
                                        (j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                        (k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                            ghostcell_dim_1_u_qi[qi];
                                        
                                    product[idx_product] *= u_qi[qi][idx_u_qi];
                                }
                            }
                        }
                    }
                }
                
                const Real weight = Real(dx[0]*dx[1]*dx[2])/(L_x*L_y*L_z);
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    if (derivative_direction == 0)
                    {
                        derivative_first_order_x->computeDerivative(
                            data_derivative,
                            data_product,
                            Real(dx[0]),
                            patch_visible_box,
                            0,
                            0);
                    }
                    else if (derivative_direction == 1)
                    {
                        derivative_first_order_y->computeDerivative(
                            data_derivative,
                            data_product,
                            Real(dx[1]),
                            patch_visible_box,
                            0,
                            0);
                    }
                    else if (derivative_direction == 2)
                    {
                         derivative_first_order_z->computeDerivative(
                            data_derivative,
                            data_product,
                            Real(dx[2]),
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
                                 * Compute the index of the data point and count how many times the data is repeated.
                                 */
                                
                                const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                                
                                int n_overlapped = 1;
                                
                                for (hier::BoxContainer::BoxContainerConstIterator iob(
                                        patch_overlapped_visible_boxes.begin());
                                     iob != patch_overlapped_visible_boxes.end();
                                     iob++)
                                {
                                    const hier::Box& patch_overlapped_visible_box = *iob;
                                    
                                    if (patch_overlapped_visible_box.contains(idx_pt))
                                    {
                                        n_overlapped++;
                                    }
                                }
                                
                                /*
                                 * Compute the linear index and the data to add.
                                 */
                                
                                const int idx = (relative_idx_lo_0 + i) +
                                    (relative_idx_lo_1 + j)*patch_interior_dim_0 +
                                    (relative_idx_lo_2 + k)*patch_interior_dim_0*
                                        patch_interior_dim_1;
                                
                                der_avg_local += (der[idx]*weight)/((Real) n_overlapped);
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average of derivative.
         */
        
        d_mpi.Allreduce(
            &der_avg_local,
            &der_avg_global,
            1,
            MPI_DOUBLE,
            MPI_SUM);
    }
    
    return der_avg_global;
}

/*
 * Compute averaged value with only x-direction as inhomogeneous direction on the coarsest level.
 */
std::vector<Real>
FlowModelMPIHelperAverage::getAveragedQuantityWithInhomogeneousXDirectionOnCoarsestLevel(
    const std::string quantity_name,
    const int component_idx,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    std::vector<Real> averaged_quantity;
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space. Coarsest and finest levels are both set to zero.
     */
    
    HAMERS_SHARED_PTR<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *d_patch_hierarchy,
            0,
            0));
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        const int coarsest_level_dim_0 = d_coarsest_level_dims[0];
        
        Real* u_avg_local = (Real*)std::malloc(coarsest_level_dim_0*sizeof(Real));
        
        averaged_quantity.resize(coarsest_level_dim_0);
        Real* u_avg_global = averaged_quantity.data();
        
        for (int i = 0; i < coarsest_level_dim_0; i++)
        {
            u_avg_local[i]  = Real(0);
            u_avg_global[i] = Real(0);
        }
        
        /*
         * Get the coarsest patch level.
         */
        
        HAMERS_SHARED_PTR<hier::PatchLevel> patch_level_coarsest(
            d_patch_hierarchy->getPatchLevel(0));
        
        for (hier::PatchLevel::iterator ip(patch_level_coarsest->begin());
             ip != patch_level_coarsest->end();
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
            
            setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
            
            std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
            
            num_subghosts_of_data.insert(
                std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
            
            registerDerivedVariables(num_subghosts_of_data);
            
            allocateMemoryForDerivedCellData();
            
            computeDerivedCellData();
            
            /*
             * Get the pointer to data inside the flow model.
             */
            
            HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity =
                getCellData(quantity_name);
            
            Real* u = data_quantity->getPointer(component_idx);
            
            const hier::BoxContainer& patch_visible_boxes =
                flattened_hierarchy->getVisibleBoxes(
                    patch_box,
                    0);
            
            const hier::BoxContainer& patch_overlapped_visible_boxes =
                flattened_hierarchy->getOverlappedVisibleBoxes(
                    patch_box,
                    0);
            
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
                     * Compute the index of the data point and count how many times the data is repeated.
                     */
                    
                    const hier::Index idx_pt(tbox::Dimension(1), idx_lo_0 + i);
                    
                    int n_overlapped = 1;
                    
                    for (hier::BoxContainer::BoxContainerConstIterator iob(
                            patch_overlapped_visible_boxes.begin());
                         iob != patch_overlapped_visible_boxes.end();
                         iob++)
                    {
                        const hier::Box& patch_overlapped_visible_box = *iob;
                        
                        if (patch_overlapped_visible_box.contains(idx_pt))
                        {
                            n_overlapped++;
                        }
                    }
                    
                    /*
                     * Compute the linear indices and the data to add.
                     */
                    
                    const int idx     = relative_idx_lo_0 + i + num_ghosts_0_quantity;
                    const int idx_avg = idx_lo_0 + i;
                    
                    const Real value_to_add = u[idx]/((Real) n_overlapped);
                    
                    u_avg_local[idx_avg] += value_to_add;
                }
            }
            
            /*
             * Unregister the patch and data of all registered derived cell variables in the flow model.
             */
            
            unregisterPatch();
        }
        
        /*
         * Reduction to get the global average.
         */
        
        d_mpi.Allreduce(
            u_avg_local,
            u_avg_global,
            coarsest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_avg_local);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int coarsest_level_dim_0 = d_coarsest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_y = Real(x_hi[1] - x_lo[1]);
        
        Real* u_avg_local = (Real*)std::malloc(coarsest_level_dim_0*sizeof(Real));
        
        averaged_quantity.resize(coarsest_level_dim_0);
        Real* u_avg_global = averaged_quantity.data();
        
        for (int i = 0; i < coarsest_level_dim_0; i++)
        {
            u_avg_local[i]  = Real(0);
            u_avg_global[i] = Real(0);
        }
        
        /*
         * Get the coarsest patch level.
         */
        
        HAMERS_SHARED_PTR<hier::PatchLevel> patch_level_coarsest(
            d_patch_hierarchy->getPatchLevel(0));
        
        for (hier::PatchLevel::iterator ip(patch_level_coarsest->begin());
             ip != patch_level_coarsest->end();
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
            
            setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
            
            std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
            
            num_subghosts_of_data.insert(
                std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
            
            registerDerivedVariables(num_subghosts_of_data);
            
            allocateMemoryForDerivedCellData();
            
            computeDerivedCellData();
            
            /*
             * Get the pointer to data inside the flow model.
             */
            
            HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity =
                getCellData(quantity_name);
            
            Real* u = data_quantity->getPointer(component_idx);
            
            const hier::BoxContainer& patch_visible_boxes =
                flattened_hierarchy->getVisibleBoxes(
                    patch_box,
                    0);
            
            const hier::BoxContainer& patch_overlapped_visible_boxes =
                flattened_hierarchy->getOverlappedVisibleBoxes(
                    patch_box,
                    0);
            
            const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
            const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
            
            const int num_ghosts_0_quantity = num_ghosts_quantity[0];
            const int num_ghosts_1_quantity = num_ghosts_quantity[1];
            const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
            
            const Real weight = Real(dx[1])/L_y;
            
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
                         * Compute the index of the data point and count how many times the data is repeated.
                         */
                        
                        const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j);
                        
                        int n_overlapped = 1;
                        
                        for (hier::BoxContainer::BoxContainerConstIterator iob(
                                patch_overlapped_visible_boxes.begin());
                             iob != patch_overlapped_visible_boxes.end();
                             iob++)
                        {
                            const hier::Box& patch_overlapped_visible_box = *iob;
                            
                            if (patch_overlapped_visible_box.contains(idx_pt))
                            {
                                n_overlapped++;
                            }
                        }
                        
                        /*
                         * Compute the linear indices and the data to add.
                         */
                        
                        const int idx = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                            (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                        
                        const int idx_avg = idx_lo_0 + i;
                        
                        const Real value_to_add = u[idx]*weight/((Real) n_overlapped);
                        
                        u_avg_local[idx_avg] += value_to_add;
                    }
                }
            }
            
            /*
             * Unregister the patch and data of all registered derived cell variables in the flow model.
             */
            
            unregisterPatch();
        }
        
        /*
         * Reduction to get the global average.
         */
        
        d_mpi.Allreduce(
            u_avg_local,
            u_avg_global,
            coarsest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_avg_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int coarsest_level_dim_0 = d_coarsest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_y = Real(x_hi[1] - x_lo[1]);
        const Real L_z = Real(x_hi[2] - x_lo[2]);
        
        Real* u_avg_local = (Real*)std::malloc(coarsest_level_dim_0*sizeof(Real));
        
        averaged_quantity.resize(coarsest_level_dim_0);
        Real* u_avg_global = averaged_quantity.data();
        
        for (int i = 0; i < coarsest_level_dim_0; i++)
        {
            u_avg_local[i]  = Real(0);
            u_avg_global[i] = Real(0);
        }
        
        /*
         * Get the coarsest patch level.
         */
        
        HAMERS_SHARED_PTR<hier::PatchLevel> patch_level_coarsest(
            d_patch_hierarchy->getPatchLevel(0));
        
        for (hier::PatchLevel::iterator ip(patch_level_coarsest->begin());
             ip != patch_level_coarsest->end();
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
            
            setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
            
            std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
            
            num_subghosts_of_data.insert(
                std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
            
            registerDerivedVariables(num_subghosts_of_data);
            
            allocateMemoryForDerivedCellData();
            
            computeDerivedCellData();
            
            /*
             * Get the pointer to data inside the flow model.
             */
            
            HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity =
                getCellData(quantity_name);
            
            Real* u = data_quantity->getPointer(component_idx);
            
            const hier::BoxContainer& patch_visible_boxes =
                flattened_hierarchy->getVisibleBoxes(
                    patch_box,
                    0);
            
            const hier::BoxContainer& patch_overlapped_visible_boxes =
                flattened_hierarchy->getOverlappedVisibleBoxes(
                    patch_box,
                    0);
            
            const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
            const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
            
            const int num_ghosts_0_quantity = num_ghosts_quantity[0];
            const int num_ghosts_1_quantity = num_ghosts_quantity[1];
            const int num_ghosts_2_quantity = num_ghosts_quantity[2];
            const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
            const int ghostcell_dim_1_quantity = ghostcell_dims_quantity[1];
            
            const Real weight = Real(dx[1]*dx[2])/(L_y*L_z);
            
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
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                            
                            int n_overlapped = 1;
                            
                            for (hier::BoxContainer::BoxContainerConstIterator iob(
                                    patch_overlapped_visible_boxes.begin());
                                 iob != patch_overlapped_visible_boxes.end();
                                 iob++)
                            {
                                const hier::Box& patch_overlapped_visible_box = *iob;
                                
                                if (patch_overlapped_visible_box.contains(idx_pt))
                                {
                                    n_overlapped++;
                                }
                            }
                            
                            /*
                             * Compute the linear indices and the data to add.
                             */
                            
                            const int idx = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                    ghostcell_dim_1_quantity;
                            
                            const int idx_avg = idx_lo_0 + i;
                            
                            const Real value_to_add = u[idx]*weight/((Real) n_overlapped);
                            
                            u_avg_local[idx_avg] += value_to_add;
                        }
                    }
                }
            }
            
            /*
             * Unregister the patch and data of all registered derived cell variables in the flow model.
             */
            
            unregisterPatch();
        }
        
        /*
         * Reduction to get the global average.
         */
        
        d_mpi.Allreduce(
            u_avg_local,
            u_avg_global,
            coarsest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_avg_local);
    }
    
    return averaged_quantity;
}


/*
 * Compute averaged reciprocal of value with only x-direction as inhomogeneous direction on the coarsest level.
 */
std::vector<Real>
FlowModelMPIHelperAverage::getAveragedReciprocalOfQuantityWithInhomogeneousXDirectionOnCoarsestLevel(
    const std::string quantity_name,
    const int component_idx,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    std::vector<Real> averaged_reciprocal_quantity;
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space. Coarsest and finest levels are both set to zero.
     */
    
    HAMERS_SHARED_PTR<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *d_patch_hierarchy,
            0,
            0));
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        const int coarsest_level_dim_0 = d_finest_level_dims[0];
        
        Real* u_inv_avg_local = (Real*)std::malloc(coarsest_level_dim_0*sizeof(Real));
        
        averaged_reciprocal_quantity.resize(coarsest_level_dim_0);
        Real* u_inv_avg_global = averaged_reciprocal_quantity.data();
        
        for (int i = 0; i < coarsest_level_dim_0; i++)
        {
            u_inv_avg_local[i]  = Real(0);
            u_inv_avg_global[i] = Real(0);
        }
        
        /*
         * Get the coarsest patch level.
         */
        
        HAMERS_SHARED_PTR<hier::PatchLevel> patch_level_coarsest(
            d_patch_hierarchy->getPatchLevel(0));
        
        for (hier::PatchLevel::iterator ip(patch_level_coarsest->begin());
             ip != patch_level_coarsest->end();
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
            
            setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
            
            std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
            
            num_subghosts_of_data.insert(
                std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
            
            registerDerivedVariables(num_subghosts_of_data);
            
            allocateMemoryForDerivedCellData();
            
            computeDerivedCellData();
            
            /*
             * Get the pointer to data inside the flow model.
             */
            
            HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity =
                getCellData(quantity_name);
            
            Real* u = data_quantity->getPointer(component_idx);
            
            const hier::BoxContainer& patch_visible_boxes =
                flattened_hierarchy->getVisibleBoxes(
                    patch_box,
                    0);
            
            const hier::BoxContainer& patch_overlapped_visible_boxes =
                flattened_hierarchy->getOverlappedVisibleBoxes(
                    patch_box,
                    0);
            
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
                     * Compute the index of the data point and count how many times the data is repeated.
                     */
                    
                    const hier::Index idx_pt(tbox::Dimension(1), idx_lo_0 + i);
                    
                    int n_overlapped = 1;
                    
                    for (hier::BoxContainer::BoxContainerConstIterator iob(
                            patch_overlapped_visible_boxes.begin());
                         iob != patch_overlapped_visible_boxes.end();
                         iob++)
                    {
                        const hier::Box& patch_overlapped_visible_box = *iob;
                        
                        if (patch_overlapped_visible_box.contains(idx_pt))
                        {
                            n_overlapped++;
                        }
                    }
                    
                    /*
                     * Compute the linear indices and the data to add.
                     */
                    
                    const int idx     = relative_idx_lo_0 + i + num_ghosts_0_quantity;
                    const int idx_avg = idx_lo_0 + i;
                    
                    const Real value_to_add = (Real(1)/u[idx])/((Real) n_overlapped);
                    
                    u_inv_avg_local[idx_avg] += value_to_add;
                }
            }
            
            /*
             * Unregister the patch and data of all registered derived cell variables in the flow model.
             */
            
            unregisterPatch();
        }
        
        /*
         * Reduction to get the global average of reciprocal.
         */
        
        d_mpi.Allreduce(
            u_inv_avg_local,
            u_inv_avg_global,
            coarsest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_inv_avg_local);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int coarsest_level_dim_0 = d_finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_y = Real(x_hi[1] - x_lo[1]);
        
        Real* u_inv_avg_local = (Real*)std::malloc(coarsest_level_dim_0*sizeof(Real));
        
        averaged_reciprocal_quantity.resize(coarsest_level_dim_0);
        Real* u_inv_avg_global = averaged_reciprocal_quantity.data();
        
        for (int i = 0; i < coarsest_level_dim_0; i++)
        {
            u_inv_avg_local[i]  = Real(0);
            u_inv_avg_global[i] = Real(0);
        }
        
        /*
         * Get the coarsest patch level.
         */
        
        HAMERS_SHARED_PTR<hier::PatchLevel> patch_level_coarsest(
            d_patch_hierarchy->getPatchLevel(0));
        
        for (hier::PatchLevel::iterator ip(patch_level_coarsest->begin());
             ip != patch_level_coarsest->end();
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
            
            setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
            
            std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
            
            num_subghosts_of_data.insert(
                std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
            
            registerDerivedVariables(num_subghosts_of_data);
            
            allocateMemoryForDerivedCellData();
            
            computeDerivedCellData();
            
            /*
             * Get the pointer to data inside the flow model.
             */
            
            HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity =
                getCellData(quantity_name);
            
            Real* u = data_quantity->getPointer(component_idx);
            
            const hier::BoxContainer& patch_visible_boxes =
                flattened_hierarchy->getVisibleBoxes(
                    patch_box,
                    0);
            
            const hier::BoxContainer& patch_overlapped_visible_boxes =
                flattened_hierarchy->getOverlappedVisibleBoxes(
                    patch_box,
                    0);
            
            const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
            const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
            
            const int num_ghosts_0_quantity = num_ghosts_quantity[0];
            const int num_ghosts_1_quantity = num_ghosts_quantity[1];
            const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
            
            const Real weight = Real(dx[1])/L_y;
            
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
                         * Compute the index of the data point and count how many times the data is repeated.
                         */
                        
                        const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j);
                        
                        int n_overlapped = 1;
                        
                        for (hier::BoxContainer::BoxContainerConstIterator iob(
                                patch_overlapped_visible_boxes.begin());
                             iob != patch_overlapped_visible_boxes.end();
                             iob++)
                        {
                            const hier::Box& patch_overlapped_visible_box = *iob;
                            
                            if (patch_overlapped_visible_box.contains(idx_pt))
                            {
                                n_overlapped++;
                            }
                        }
                        
                        /*
                         * Compute the linear indices and the data to add.
                         */
                        
                        const int idx = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                            (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                        
                        const int idx_avg = idx_lo_0 + i;
                        
                        const Real value_to_add = (Real(1)/u[idx])*weight/((Real) n_overlapped);
                        
                        u_inv_avg_local[idx_avg] += value_to_add;
                    }
                }
            }
            
            /*
             * Unregister the patch and data of all registered derived cell variables in the flow model.
             */
            
            unregisterPatch();
        }
        
        /*
         * Reduction to get the global average of reciprocal.
         */
        
        d_mpi.Allreduce(
            u_inv_avg_local,
            u_inv_avg_global,
            coarsest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_inv_avg_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int coarsest_level_dim_0 = d_finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_y = Real(x_hi[1] - x_lo[1]);
        const Real L_z = Real(x_hi[2] - x_lo[2]);
        
        Real* u_inv_avg_local = (Real*)std::malloc(coarsest_level_dim_0*sizeof(Real));
        
        averaged_reciprocal_quantity.resize(coarsest_level_dim_0);
        Real* u_inv_avg_global = averaged_reciprocal_quantity.data();
        
        for (int i = 0; i < coarsest_level_dim_0; i++)
        {
            u_inv_avg_local[i]  = Real(0);
            u_inv_avg_global[i] = Real(0);
        }
        
        /*
         * Get the coarsest patch level.
         */
        
        HAMERS_SHARED_PTR<hier::PatchLevel> patch_level_coarsest(
            d_patch_hierarchy->getPatchLevel(0));
        
        for (hier::PatchLevel::iterator ip(patch_level_coarsest->begin());
             ip != patch_level_coarsest->end();
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
            
            setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
            
            std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
            
            num_subghosts_of_data.insert(
                std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
            
            registerDerivedVariables(num_subghosts_of_data);
            
            allocateMemoryForDerivedCellData();
            
            computeDerivedCellData();
            
            /*
             * Get the pointer to data inside the flow model.
             */
            
            HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity =
                getCellData(quantity_name);
            
            Real* u = data_quantity->getPointer(component_idx);
            
            const hier::BoxContainer& patch_visible_boxes =
                flattened_hierarchy->getVisibleBoxes(
                    patch_box,
                    0);
            
            const hier::BoxContainer& patch_overlapped_visible_boxes =
                flattened_hierarchy->getOverlappedVisibleBoxes(
                    patch_box,
                    0);
            
            const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
            const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
            
            const int num_ghosts_0_quantity = num_ghosts_quantity[0];
            const int num_ghosts_1_quantity = num_ghosts_quantity[1];
            const int num_ghosts_2_quantity = num_ghosts_quantity[2];
            const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
            const int ghostcell_dim_1_quantity = ghostcell_dims_quantity[1];
            
            const Real weight = Real(dx[1]*dx[2])/(L_y*L_z);
            
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
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                            
                            int n_overlapped = 1;
                            
                            for (hier::BoxContainer::BoxContainerConstIterator iob(
                                    patch_overlapped_visible_boxes.begin());
                                 iob != patch_overlapped_visible_boxes.end();
                                 iob++)
                            {
                                const hier::Box& patch_overlapped_visible_box = *iob;
                                
                                if (patch_overlapped_visible_box.contains(idx_pt))
                                {
                                    n_overlapped++;
                                }
                            }
                            
                            /*
                             * Compute the linear indices and the data to add.
                             */
                            
                            const int idx = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                    ghostcell_dim_1_quantity;
                            
                            const int idx_avg = idx_lo_0 + i;
                            
                            const Real value_to_add = (Real(1)/u[idx])*weight/((Real) n_overlapped);
                            
                            u_inv_avg_local[idx_avg] += value_to_add;
                        }
                    }
                }
            }
            
            /*
             * Unregister the patch and data of all registered derived cell variables in the flow model.
             */
            
            unregisterPatch();
        }
        
        /*
         * Reduction to get the global average of reciprocal.
         */
        
        d_mpi.Allreduce(
            u_inv_avg_local,
            u_inv_avg_global,
            coarsest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_inv_avg_local);
    }
    
    return averaged_reciprocal_quantity;
}


/*
 * Compute averaged value (on product of variables) with only x-direction as inhomogeneous direction on
 * the coarsest level.
 */
std::vector<Real>
FlowModelMPIHelperAverage::getAveragedQuantityWithInhomogeneousXDirectionOnCoarsestLevel(
    const std::vector<std::string>& quantity_names,
    const std::vector<int>& component_indices,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    int num_quantities = static_cast<int>(quantity_names.size());
    
    TBOX_ASSERT(static_cast<int>(component_indices.size()) == num_quantities);
    
    std::vector<bool> use_reciprocal(num_quantities, false);
    
    return getAveragedQuantityWithInhomogeneousXDirectionOnCoarsestLevel(
        quantity_names,
        component_indices,
        use_reciprocal,
        data_context);
}


/*
 * Compute averaged value (on product of variables) with only x-direction as inhomogeneous direction on
 * the coarsest level.
 */
std::vector<Real>
FlowModelMPIHelperAverage::getAveragedQuantityWithInhomogeneousXDirectionOnCoarsestLevel(
    const std::vector<std::string>& quantity_names,
    const std::vector<int>& component_indices,
    const std::vector<bool>& use_reciprocal,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    int num_quantities = static_cast<int>(quantity_names.size());
    
    TBOX_ASSERT(static_cast<int>(component_indices.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(use_reciprocal.size()) == num_quantities);
    
    std::vector<Real> averaged_quantity;
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space. Coarsest and finest levels are both set to zero.
     */
    
    HAMERS_SHARED_PTR<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *d_patch_hierarchy,
            0,
            0));
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        const int coarsest_level_dim_0 = d_coarsest_level_dims[0];
        
        Real* avg_local = (Real*)std::malloc(coarsest_level_dim_0*sizeof(Real));
        
        averaged_quantity.resize(coarsest_level_dim_0);
        Real* avg_global = averaged_quantity.data();
        
        for (int i = 0; i < coarsest_level_dim_0; i++)
        {
            avg_local[i]  = Real(0);
            avg_global[i] = Real(0);
        }
        
        /*
         * Get the coarsest patch level.
         */
        
        HAMERS_SHARED_PTR<hier::PatchLevel> patch_level_coarsest(
            d_patch_hierarchy->getPatchLevel(0));
        
        for (hier::PatchLevel::iterator ip(patch_level_coarsest->begin());
             ip != patch_level_coarsest->end();
             ip++)
        {
            const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
            
            /*
             * Get the patch lower indices.
             */
            
            const hier::Box& patch_box = patch->getBox();
            
            const hier::Index& patch_index_lo = patch_box.lower();
            
            /*
             * Register the patch and data in the flow model and compute the corresponding
             * average.
             */
            
            setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
            
            std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
            
            for (int qi = 0; qi < num_quantities; qi++)
            {
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
            }
            
            registerDerivedVariables(num_subghosts_of_data);
            
            allocateMemoryForDerivedCellData();
            
            computeDerivedCellData();
            
            /*
             * Get the pointers to data inside the flow model.
             */
            
            std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > data_quantities;
            data_quantities.resize(num_quantities);
            for (int qi = 0; qi < num_quantities; qi++)
            {
                data_quantities[qi] = getCellData(quantity_names[qi]);
            }
            
            std::vector<Real*> u_qi;
            u_qi.resize(num_quantities);
            for (int qi = 0; qi < num_quantities; qi++)
            {
                u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
            }
            
            const hier::BoxContainer& patch_visible_boxes =
                flattened_hierarchy->getVisibleBoxes(
                    patch_box,
                    0);
            
            const hier::BoxContainer& patch_overlapped_visible_boxes =
                flattened_hierarchy->getOverlappedVisibleBoxes(
                    patch_box,
                    0);
            
            std::vector<int> num_ghosts_0_u_qi;
            num_ghosts_0_u_qi.reserve(num_quantities);
            
            for (int qi = 0; qi < num_quantities; qi++)
            {
                const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
            }
            
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
                     * Compute the index of the data point and count how many times the data is repeated.
                     */
                    
                    const hier::Index idx_pt(tbox::Dimension(1), idx_lo_0 + i);
                    
                    int n_overlapped = 1;
                    
                    for (hier::BoxContainer::BoxContainerConstIterator iob(
                            patch_overlapped_visible_boxes.begin());
                         iob != patch_overlapped_visible_boxes.end();
                         iob++)
                    {
                        const hier::Box& patch_overlapped_visible_box = *iob;
                        
                        if (patch_overlapped_visible_box.contains(idx_pt))
                        {
                            n_overlapped++;
                        }
                    }
                    
                    /*
                     * Compute the linear indices and the data to add.
                     */
                    
                    const int idx_avg = idx_lo_0 + i;
                    
                    Real avg = Real(1);
                    
                    for (int qi = 0; qi < num_quantities; qi++)
                    {
                        const int idx_qi = relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi];
                        
                        if (use_reciprocal[qi])
                        {
                            avg /= u_qi[qi][idx_qi];
                        }
                        else
                        {
                            avg *= u_qi[qi][idx_qi];
                        }
                    }
                    
                    avg_local[idx_avg] += (avg/((Real) n_overlapped));
                }
            }
            
            /*
             * Unregister the patch and data of all registered derived cell variables in the flow model.
             */
            
            unregisterPatch();
        }
        
        /*
         * Reduction to get the global average.
         */
        
        d_mpi.Allreduce(
            avg_local,
            avg_global,
            coarsest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(avg_local);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int coarsest_level_dim_0 = d_coarsest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_y = Real(x_hi[1] - x_lo[1]);
        
        Real* avg_local = (Real*)std::malloc(coarsest_level_dim_0*sizeof(Real));
        
        averaged_quantity.resize(coarsest_level_dim_0);
        Real* avg_global = averaged_quantity.data();
        
        for (int i = 0; i < coarsest_level_dim_0; i++)
        {
            avg_local[i]  = Real(0);
            avg_global[i] = Real(0);
        }
        
        /*
         * Get the coarsest patch level.
         */
        
        HAMERS_SHARED_PTR<hier::PatchLevel> patch_level_coarsest(
            d_patch_hierarchy->getPatchLevel(0));
        
        for (hier::PatchLevel::iterator ip(patch_level_coarsest->begin());
             ip != patch_level_coarsest->end();
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
             * Register the patch and data in the flow model and compute the corresponding
             * average.
             */
            
            setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
            
            std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
            
            for (int qi = 0; qi < num_quantities; qi++)
            {
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
            }
            
            registerDerivedVariables(num_subghosts_of_data);
            
            allocateMemoryForDerivedCellData();
            
            computeDerivedCellData();
            
            /*
             * Get the pointers to data inside the flow model.
             */
            
            std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > data_quantities;
            data_quantities.resize(num_quantities);
            for (int qi = 0; qi < num_quantities; qi++)
            {
                data_quantities[qi] = getCellData(quantity_names[qi]);
            }
            
            std::vector<Real*> u_qi;
            u_qi.resize(num_quantities);
            for (int qi = 0; qi < num_quantities; qi++)
            {
                u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
            }
            
            const hier::BoxContainer& patch_visible_boxes =
                flattened_hierarchy->getVisibleBoxes(
                    patch_box,
                    0);
            
            const hier::BoxContainer& patch_overlapped_visible_boxes =
                flattened_hierarchy->getOverlappedVisibleBoxes(
                    patch_box,
                    0);
            
            std::vector<int> num_ghosts_0_u_qi;
            std::vector<int> num_ghosts_1_u_qi;
            std::vector<int> ghostcell_dim_0_u_qi;
            num_ghosts_0_u_qi.reserve(num_quantities);
            num_ghosts_1_u_qi.reserve(num_quantities);
            ghostcell_dim_0_u_qi.reserve(num_quantities);
            
            for (int qi = 0; qi < num_quantities; qi++)
            {
                const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                
                num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
            }
            
            const Real weight = Real(dx[1])/L_y;
            
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
                         * Compute the index of the data point and count how many times the data is repeated.
                         */
                        
                        const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j);
                        
                        int n_overlapped = 1;
                        
                        for (hier::BoxContainer::BoxContainerConstIterator iob(
                                patch_overlapped_visible_boxes.begin());
                             iob != patch_overlapped_visible_boxes.end();
                             iob++)
                        {
                            const hier::Box& patch_overlapped_visible_box = *iob;
                            
                            if (patch_overlapped_visible_box.contains(idx_pt))
                            {
                                n_overlapped++;
                            }
                        }
                        
                        /*
                         * Compute the linear indices and the data to add.
                         */
                        
                        const int idx_avg = idx_lo_0 + i;
                        
                        Real avg = Real(1);
                        
                        for (int qi = 0; qi < num_quantities; qi++)
                        {
                            const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                            
                            if (use_reciprocal[qi])
                            {
                                avg /= u_qi[qi][idx_qi];
                            }
                            else
                            {
                                avg *= u_qi[qi][idx_qi];
                            }
                        }
                        
                        avg_local[idx_avg] += (avg*weight/((Real) n_overlapped));
                    }
                }
            }
            
            /*
             * Unregister the patch and data of all registered derived cell variables in the flow model.
             */
            
            unregisterPatch();
        }
        
        /*
         * Reduction to get the global average.
         */
        
        d_mpi.Allreduce(
            avg_local,
            avg_global,
            coarsest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(avg_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int coarsest_level_dim_0 = d_coarsest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_y = Real(x_hi[1] - x_lo[1]);
        const Real L_z = Real(x_hi[2] - x_lo[2]);
        
        Real* avg_local = (Real*)std::malloc(coarsest_level_dim_0*sizeof(Real));
        
        averaged_quantity.resize(coarsest_level_dim_0);
        Real* avg_global = averaged_quantity.data();
        
        for (int i = 0; i < coarsest_level_dim_0; i++)
        {
            avg_local[i]  = Real(0);
            avg_global[i] = Real(0);
        }
        
        /*
         * Get the coarsest patch level.
         */
        
        HAMERS_SHARED_PTR<hier::PatchLevel> patch_level_coarsest(
            d_patch_hierarchy->getPatchLevel(0));
        
        for (hier::PatchLevel::iterator ip(patch_level_coarsest->begin());
             ip != patch_level_coarsest->end();
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
             * Register the patch and data in the flow model and compute the corresponding
             * average.
             */
            
            setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
            
            std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
            
            for (int qi = 0; qi < num_quantities; qi++)
            {
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
            }
            
            registerDerivedVariables(num_subghosts_of_data);
            
            allocateMemoryForDerivedCellData();
            
            computeDerivedCellData();
            
            /*
             * Get the pointers to data inside the flow model.
             */
            
            std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > data_quantities;
            data_quantities.resize(num_quantities);
            for (int qi = 0; qi < num_quantities; qi++)
            {
                data_quantities[qi] = getCellData(quantity_names[qi]);
            }
            
            std::vector<Real*> u_qi;
            u_qi.resize(num_quantities);
            for (int qi = 0; qi < num_quantities; qi++)
            {
                u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
            }
            
            const hier::BoxContainer& patch_visible_boxes =
                flattened_hierarchy->getVisibleBoxes(
                    patch_box,
                    0);
            
            const hier::BoxContainer& patch_overlapped_visible_boxes =
                flattened_hierarchy->getOverlappedVisibleBoxes(
                    patch_box,
                    0);
            
            std::vector<int> num_ghosts_0_u_qi;
            std::vector<int> num_ghosts_1_u_qi;
            std::vector<int> num_ghosts_2_u_qi;
            std::vector<int> ghostcell_dim_0_u_qi;
            std::vector<int> ghostcell_dim_1_u_qi;
            num_ghosts_0_u_qi.reserve(num_quantities);
            num_ghosts_1_u_qi.reserve(num_quantities);
            num_ghosts_2_u_qi.reserve(num_quantities);
            ghostcell_dim_0_u_qi.reserve(num_quantities);
            ghostcell_dim_1_u_qi.reserve(num_quantities);
            
            for (int qi = 0; qi < num_quantities; qi++)
            {
                const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                
                num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                num_ghosts_2_u_qi.push_back(num_ghosts_u_qi[2]);
                ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                ghostcell_dim_1_u_qi.push_back(ghostcell_dims_u_qi[1]);
            }
            
            const Real weight = Real(dx[1]*dx[2])/(L_y*L_z);
            
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
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                            
                            int n_overlapped = 1;
                            
                            for (hier::BoxContainer::BoxContainerConstIterator iob(
                                    patch_overlapped_visible_boxes.begin());
                                 iob != patch_overlapped_visible_boxes.end();
                                 iob++)
                            {
                                const hier::Box& patch_overlapped_visible_box = *iob;
                                
                                if (patch_overlapped_visible_box.contains(idx_pt))
                                {
                                    n_overlapped++;
                                }
                            }
                            
                            /*
                             * Compute the linear index and the data to add.
                             */
                            
                            const int idx_avg = idx_lo_0 + i;
                            
                            Real avg = Real(1);
                            
                            for (int qi = 0; qi < num_quantities; qi++)
                            {
                                const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                    (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                        ghostcell_dim_1_u_qi[qi];
                                
                                if (use_reciprocal[qi])
                                {
                                    avg /= u_qi[qi][idx_qi];
                                }
                                else
                                {
                                    avg *= u_qi[qi][idx_qi];
                                }
                            }
                            
                            avg_local[idx_avg] += (avg*weight/((Real) n_overlapped));
                        }
                    }
                }
            }
            
            /*
             * Unregister the patch and data of all registered derived cell variables in the flow model.
             */
            
            unregisterPatch();
        }
        
        /*
         * Reduction to get the global average.
         */
        
        d_mpi.Allreduce(
            avg_local,
            avg_global,
            coarsest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(avg_local);
    }
    
    return averaged_quantity;
}


/*
 * Compute averaged value (on product of variable derivatives) with only x direction as inhomogeneous direction
 * on the coarsest level.
 */
std::vector<Real>
FlowModelMPIHelperAverage::getAveragedQuantityWithInhomogeneousXDirectionOnCoarsestLevel(
    const std::vector<std::string>& quantity_names,
    const std::vector<int>& component_indices,
    const std::vector<bool>& use_derivative,
    const std::vector<int>& derivative_directions,
    const int num_ghosts_derivative,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    int num_quantities = static_cast<int>(quantity_names.size());
    
    TBOX_ASSERT(static_cast<int>(component_indices.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(use_derivative.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(derivative_directions.size()) == num_quantities);
    
    std::vector<bool> use_reciprocal(num_quantities, false);
    
    return getAveragedQuantityWithInhomogeneousXDirectionOnCoarsestLevel(
        quantity_names,
        component_indices,
        use_derivative,
        derivative_directions,
        use_reciprocal,
        num_ghosts_derivative,
        data_context);
}


/*
 * Compute averaged value (on product of variable derivatives) with only x direction as inhomogeneous direction
 * on the coarsest level.
 */
std::vector<Real>
FlowModelMPIHelperAverage::getAveragedQuantityWithInhomogeneousXDirectionOnCoarsestLevel(
    const std::vector<std::string>& quantity_names,
    const std::vector<int>& component_indices,
    const std::vector<bool>& use_derivative,
    const std::vector<int>& derivative_directions,
    const std::vector<bool>& use_reciprocal,
    const int num_ghosts_derivative,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    int num_quantities = static_cast<int>(quantity_names.size());
    
    TBOX_ASSERT(static_cast<int>(component_indices.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(use_derivative.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(derivative_directions.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(use_reciprocal.size()) == num_quantities);
    
    int num_use_derivative = 0;
    
    if (d_dim == tbox::Dimension(1))
    {
        for (int qi = 0; qi < num_quantities; qi++)
        {
            if (use_derivative[qi])
            {
                num_use_derivative++;
                if (derivative_directions[qi] < 0 || derivative_directions[qi] > 0)
                {
                    TBOX_ERROR(d_object_name
                        << ": FlowModelMPIHelperAverage::getAveragedQuantityWithInhomogeneousXDirectionOnCoarsestLevel():\n"
                        << "Cannot take derivative for one-dimensional problem!\n"
                        << "derivative_directions[" << qi << "] = " << derivative_directions[qi] << " given!\n"
                        << std::endl);
                }
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        for (int qi = 0; qi < num_quantities; qi++)
        {
            if (use_derivative[qi])
            {
                num_use_derivative++;
                if (derivative_directions[qi] < 0 || derivative_directions[qi] > 1)
                {
                    TBOX_ERROR(d_object_name
                        << ": FlowModelMPIHelperAverage::getAveragedQuantityWithInhomogeneousXDirectionOnCoarsestLevel():\n"
                        << "Cannot take derivative for two-dimensional problem!\n"
                        << "derivative_directions[" << qi << "] = " << derivative_directions[qi] << " given!\n"
                        << std::endl);
                }
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        for (int qi = 0; qi < num_quantities; qi++)
        {
            if (use_derivative[qi])
            {
                num_use_derivative++;
                if (derivative_directions[qi] < 0 || derivative_directions[qi] > 2)
                {
                    TBOX_ERROR(d_object_name
                        << ": FlowModelMPIHelperAverage::getAveragedQuantityWithInhomogeneousXDirectionOnCoarsestLevel():\n"
                        << "Cannot take derivative for three-dimensional problem!\n"
                        << "derivative_directions[" << qi << "] = " << derivative_directions[qi] << " given!\n"
                        << std::endl);
                }
            }
        }
    }
    
    std::vector<Real> averaged_quantity;
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space. Coarsest and finest levels are both set to zero.
     */
    
    HAMERS_SHARED_PTR<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *d_patch_hierarchy,
            0,
            0));
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_x(
            new DerivativeFirstOrder(
                "first order derivative in x-direction",
                d_dim,
                DIRECTION::X_DIRECTION,
                num_ghosts_derivative));
        
        hier::IntVector num_ghosts_der = hier::IntVector::getOne(d_dim)*num_ghosts_derivative;
        
        const int coarsest_level_dim_0 = d_coarsest_level_dims[0];
        
        Real* avg_local = (Real*)std::malloc(coarsest_level_dim_0*sizeof(Real));
        
        averaged_quantity.resize(coarsest_level_dim_0);
        Real* avg_global = averaged_quantity.data();
        
        for (int i = 0; i < coarsest_level_dim_0; i++)
        {
            avg_local[i]  = Real(0);
            avg_global[i] = Real(0);
        }
        
        /*
         * Get the coarsest patch level.
         */
        
        HAMERS_SHARED_PTR<hier::PatchLevel> patch_level_coarsest(
            d_patch_hierarchy->getPatchLevel(0));
        
        for (hier::PatchLevel::iterator ip(patch_level_coarsest->begin());
             ip != patch_level_coarsest->end();
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
             * Register the patch and data in the flow model and compute the corresponding
             * average.
             */
            
            setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
            
            std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
            
            for (int qi = 0; qi < num_quantities; qi++)
            {
                if (num_use_derivative > 0)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], num_ghosts_der));
                }
                else
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                }
            }
            
            registerDerivedVariables(num_subghosts_of_data);
            
            allocateMemoryForDerivedCellData();
            
            computeDerivedCellData();
            
            /*
             * Get the pointers to data inside the flow model.
             */
            
            std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > data_quantities;
            data_quantities.resize(num_quantities);
            for (int qi = 0; qi < num_quantities; qi++)
            {
                data_quantities[qi] = getCellData(quantity_names[qi]);
            }
            
            std::vector<Real*> u_qi;
            u_qi.resize(num_quantities);
            for (int qi = 0; qi < num_quantities; qi++)
            {
                u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
            }
            
            const hier::BoxContainer& patch_visible_boxes =
                flattened_hierarchy->getVisibleBoxes(
                    patch_box,
                    0);
            
            const hier::BoxContainer& patch_overlapped_visible_boxes =
                flattened_hierarchy->getOverlappedVisibleBoxes(
                    patch_box,
                    0);
            
            std::vector<int> num_ghosts_0_u_qi;
            num_ghosts_0_u_qi.reserve(num_quantities);
            
            for (int qi = 0; qi < num_quantities; qi++)
            {
                const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
            }
            
            /*
             * Initialize cell data for the derivatives and get pointers to the cell data.
             */
            
            HAMERS_SHARED_PTR<pdat::CellData<Real> > data_derivative;
            std::vector<Real*> der_qi;
            
            if (num_use_derivative > 0)
            {
                data_derivative = HAMERS_MAKE_SHARED<pdat::CellData<Real> >(patch_box, num_use_derivative, hier::IntVector::getZero(d_dim));
                
                der_qi.resize(num_use_derivative);
                for (int qi = 0; qi < num_use_derivative; qi++)
                {
                    der_qi[qi] = data_derivative->getPointer(qi);
                }
            }
            
            for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                 ib != patch_visible_boxes.end();
                 ib++)
            {
                const hier::Box& patch_visible_box = *ib;
                
                int count_derivative = 0;
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    if (use_derivative[qi] && derivative_directions[qi] == 0)
                    {
                        derivative_first_order_x->computeDerivative(
                            data_derivative,
                            data_quantities[qi],
                            Real(dx[0]),
                            patch_visible_box,
                            count_derivative,
                            component_indices[qi]);
                        
                        count_derivative++;
                    }
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
                     * Compute the index of the data point and count how many times the data is repeated.
                     */
                    
                    const hier::Index idx_pt(tbox::Dimension(1), idx_lo_0 + i);
                    
                    int n_overlapped = 1;
                    
                    for (hier::BoxContainer::BoxContainerConstIterator iob(
                            patch_overlapped_visible_boxes.begin());
                         iob != patch_overlapped_visible_boxes.end();
                         iob++)
                    {
                        const hier::Box& patch_overlapped_visible_box = *iob;
                        
                        if (patch_overlapped_visible_box.contains(idx_pt))
                        {
                            n_overlapped++;
                        }
                    }
                    
                    /*
                     * Compute the linear indices and the data to add.
                     */
                    
                    const int idx_avg = idx_lo_0 + i;
                    
                    Real avg = Real(1);
                    
                    count_derivative = 0;
                    for (int qi = 0; qi < num_quantities; qi++)
                    {
                        if (use_reciprocal[qi])
                        {
                            if (use_derivative[qi])
                            {
                                const int idx_der = relative_idx_lo_0 + i;
                                
                                avg /= der_qi[count_derivative][idx_der];
                                count_derivative++;
                            }
                            else
                            {
                                const int idx_qi = relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi];
                                
                                avg /= u_qi[qi][idx_qi];
                            }
                        }
                        else
                        {
                            if (use_derivative[qi])
                            {
                                const int idx_der = relative_idx_lo_0 + i;
                                
                                avg *= der_qi[count_derivative][idx_der];
                                count_derivative++;
                            }
                            else
                            {
                                const int idx_qi = relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi];
                                
                                avg *= u_qi[qi][idx_qi];
                            }
                        }
                    }
                    
                    avg_local[idx_avg] += (avg/((Real) n_overlapped));
                }
            }
            
            /*
             * Unregister the patch and data of all registered derived cell variables in the flow model.
             */
            
            unregisterPatch();
        }
        
        /*
         * Reduction to get the global average.
         */
        
        d_mpi.Allreduce(
            avg_local,
            avg_global,
            coarsest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(avg_local);
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
        
        hier::IntVector num_ghosts_der = hier::IntVector::getOne(d_dim)*num_ghosts_derivative;
        
        const int coarsest_level_dim_0 = d_coarsest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_y = Real(x_hi[1] - x_lo[1]);
        
        Real* avg_local = (Real*)std::malloc(coarsest_level_dim_0*sizeof(Real));
        
        averaged_quantity.resize(coarsest_level_dim_0);
        Real* avg_global = averaged_quantity.data();
        
        for (int i = 0; i < coarsest_level_dim_0; i++)
        {
            avg_local[i]  = Real(0);
            avg_global[i] = Real(0);
        }
        
        /*
         * Get the coarsest patch level.
         */
        
        HAMERS_SHARED_PTR<hier::PatchLevel> patch_level_coarsest(
            d_patch_hierarchy->getPatchLevel(0));
        
        for (hier::PatchLevel::iterator ip(patch_level_coarsest->begin());
             ip != patch_level_coarsest->end();
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
             * Register the patch and data in the flow model and compute the corresponding
             * average.
             */
            
            setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
            
            std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
            
            for (int qi = 0; qi < num_quantities; qi++)
            {
                if (num_use_derivative > 0)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], num_ghosts_der));
                }
                else
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                }
            }
            
            registerDerivedVariables(num_subghosts_of_data);
            
            allocateMemoryForDerivedCellData();
            
            computeDerivedCellData();
            
            /*
             * Get the pointers to data inside the flow model.
             */
            
            std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > data_quantities;
            data_quantities.resize(num_quantities);
            for (int qi = 0; qi < num_quantities; qi++)
            {
                data_quantities[qi] = getCellData(quantity_names[qi]);
            }
            
            std::vector<Real*> u_qi;
            u_qi.resize(num_quantities);
            for (int qi = 0; qi < num_quantities; qi++)
            {
                u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
            }
            
            const hier::BoxContainer& patch_visible_boxes =
                flattened_hierarchy->getVisibleBoxes(
                    patch_box,
                    0);
            
            const hier::BoxContainer& patch_overlapped_visible_boxes =
                flattened_hierarchy->getOverlappedVisibleBoxes(
                    patch_box,
                    0);
            
            std::vector<int> num_ghosts_0_u_qi;
            std::vector<int> num_ghosts_1_u_qi;
            std::vector<int> ghostcell_dim_0_u_qi;
            num_ghosts_0_u_qi.reserve(num_quantities);
            num_ghosts_1_u_qi.reserve(num_quantities);
            ghostcell_dim_0_u_qi.reserve(num_quantities);
            
            for (int qi = 0; qi < num_quantities; qi++)
            {
                const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                
                num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
            }
            
            /*
             * Initialize cell data for the derivatives and get pointers to the cell data.
             */
            
            HAMERS_SHARED_PTR<pdat::CellData<Real> > data_derivative;
            std::vector<Real*> der_qi;
            
            if (num_use_derivative > 0)
            {
                data_derivative = HAMERS_MAKE_SHARED<pdat::CellData<Real> >(patch_box, num_use_derivative, hier::IntVector::getZero(d_dim));
                
                der_qi.resize(num_use_derivative);
                for (int qi = 0; qi < num_use_derivative; qi++)
                {
                    der_qi[qi] = data_derivative->getPointer(qi);
                }
            }
            
            const hier::IntVector patch_interior_dims = patch_box.numberCells();
            const int patch_interior_dim_0 = patch_interior_dims[0];
            
            const Real weight = Real(dx[1])/L_y;
            
            for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                 ib != patch_visible_boxes.end();
                 ib++)
            {
                const hier::Box& patch_visible_box = *ib;
                
                int count_derivative = 0;
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    if (use_derivative[qi] && derivative_directions[qi] == 0)
                    {
                        derivative_first_order_x->computeDerivative(
                            data_derivative,
                            data_quantities[qi],
                            Real(dx[0]),
                            patch_visible_box,
                            count_derivative,
                            component_indices[qi]);
                        
                        count_derivative++;
                    }
                    else if (use_derivative[qi] && derivative_directions[qi] == 1)
                    {
                        derivative_first_order_y->computeDerivative(
                            data_derivative,
                            data_quantities[qi],
                            Real(dx[1]),
                            patch_visible_box,
                            count_derivative,
                            component_indices[qi]);
                        
                        count_derivative++;
                    }
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
                         * Compute the index of the data point and count how many times the data is repeated.
                         */
                        
                        const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j);
                        
                        int n_overlapped = 1;
                        
                        for (hier::BoxContainer::BoxContainerConstIterator iob(
                                patch_overlapped_visible_boxes.begin());
                             iob != patch_overlapped_visible_boxes.end();
                             iob++)
                        {
                            const hier::Box& patch_overlapped_visible_box = *iob;
                            
                            if (patch_overlapped_visible_box.contains(idx_pt))
                            {
                                n_overlapped++;
                            }
                        }
                        
                        /*
                         * Compute the linear indices and the data to add.
                         */
                        
                        const int idx_avg = idx_lo_0 + i;
                        
                        Real avg = Real(1);
                        
                        count_derivative = 0;
                        for (int qi = 0; qi < num_quantities; qi++)
                        {
                            if (use_reciprocal[qi])
                            {
                                if (use_derivative[qi])
                                {
                                    const int idx_der = (relative_idx_lo_0 + i) +
                                        (relative_idx_lo_1 + j)*patch_interior_dim_0;
                                    
                                    avg /= der_qi[count_derivative][idx_der];
                                    
                                    count_derivative++;
                                }
                                else
                                {
                                    const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    avg /= u_qi[qi][idx_qi];
                                }
                            }
                            else
                            {
                                if (use_derivative[qi])
                                {
                                    const int idx_der = (relative_idx_lo_0 + i) +
                                        (relative_idx_lo_1 + j)*patch_interior_dim_0;
                                    
                                    avg *= der_qi[count_derivative][idx_der];
                                    
                                    count_derivative++;
                                }
                                else
                                {
                                    const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    avg *= u_qi[qi][idx_qi];
                                }
                            }
                        }
                        
                        avg_local[idx_avg] += (avg*weight/((Real) n_overlapped));
                    }
                }
            }
            
            /*
             * Unregister the patch and data of all registered derived cell variables in the flow model.
             */
            
            unregisterPatch();
        }
        
        /*
         * Reduction to get the global average.
         */
        
        d_mpi.Allreduce(
            avg_local,
            avg_global,
            coarsest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(avg_local);
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
        
        hier::IntVector num_ghosts_der = hier::IntVector::getOne(d_dim)*num_ghosts_derivative;
        
        const int coarsest_level_dim_0 = d_coarsest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_y = Real(x_hi[1] - x_lo[1]);
        const Real L_z = Real(x_hi[2] - x_lo[2]);
        
        Real* avg_local = (Real*)std::malloc(coarsest_level_dim_0*sizeof(Real));
        
        averaged_quantity.resize(coarsest_level_dim_0);
        Real* avg_global = averaged_quantity.data();
        
        for (int i = 0; i < coarsest_level_dim_0; i++)
        {
            avg_local[i]  = Real(0);
            avg_global[i] = Real(0);
        }
        
        /*
         * Get the coarsest patch level.
         */
        
        HAMERS_SHARED_PTR<hier::PatchLevel> patch_level_coarsest(
            d_patch_hierarchy->getPatchLevel(0));
        
        for (hier::PatchLevel::iterator ip(patch_level_coarsest->begin());
             ip != patch_level_coarsest->end();
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
             * Register the patch and data in the flow model and compute the corresponding
             * average.
             */
            
            setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
            
            std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
            
            for (int qi = 0; qi < num_quantities; qi++)
            {
                if (num_use_derivative > 0)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], num_ghosts_der));
                }
                else
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                }
            }
            
            registerDerivedVariables(num_subghosts_of_data);
            
            allocateMemoryForDerivedCellData();
            
            computeDerivedCellData();
            
            /*
             * Get the pointers to data inside the flow model.
             */
            
            std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > data_quantities;
            data_quantities.resize(num_quantities);
            for (int qi = 0; qi < num_quantities; qi++)
            {
                data_quantities[qi] = getCellData(quantity_names[qi]);
            }
            
            std::vector<Real*> u_qi;
            u_qi.resize(num_quantities);
            for (int qi = 0; qi < num_quantities; qi++)
            {
                u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
            }
            
            const hier::BoxContainer& patch_visible_boxes =
                flattened_hierarchy->getVisibleBoxes(
                    patch_box,
                    0);
            
            const hier::BoxContainer& patch_overlapped_visible_boxes =
                flattened_hierarchy->getOverlappedVisibleBoxes(
                    patch_box,
                    0);
            
            std::vector<int> num_ghosts_0_u_qi;
            std::vector<int> num_ghosts_1_u_qi;
            std::vector<int> num_ghosts_2_u_qi;
            std::vector<int> ghostcell_dim_0_u_qi;
            std::vector<int> ghostcell_dim_1_u_qi;
            num_ghosts_0_u_qi.reserve(num_quantities);
            num_ghosts_1_u_qi.reserve(num_quantities);
            num_ghosts_2_u_qi.reserve(num_quantities);
            ghostcell_dim_0_u_qi.reserve(num_quantities);
            ghostcell_dim_1_u_qi.reserve(num_quantities);
            
            for (int qi = 0; qi < num_quantities; qi++)
            {
                const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                
                num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                num_ghosts_2_u_qi.push_back(num_ghosts_u_qi[2]);
                ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                ghostcell_dim_1_u_qi.push_back(ghostcell_dims_u_qi[1]);
            }
            
            /*
             * Initialize cell data for the derivatives and get pointers to the cell data.
             */
            
            HAMERS_SHARED_PTR<pdat::CellData<Real> > data_derivative;
            std::vector<Real*> der_qi;
            
            if (num_use_derivative > 0)
            {
                data_derivative = HAMERS_MAKE_SHARED<pdat::CellData<Real> >(patch_box, num_use_derivative, hier::IntVector::getZero(d_dim));
                
                der_qi.resize(num_use_derivative);
                for (int qi = 0; qi < num_use_derivative; qi++)
                {
                    der_qi[qi] = data_derivative->getPointer(qi);
                }
            }
            
            const hier::IntVector patch_interior_dims = patch_box.numberCells();
            const int patch_interior_dim_0 = patch_interior_dims[0];
            const int patch_interior_dim_1 = patch_interior_dims[1];
            
            const Real weight = Real(dx[1]*dx[2])/(L_y*L_z);
            
            for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                 ib != patch_visible_boxes.end();
                 ib++)
            {
                const hier::Box& patch_visible_box = *ib;
                
                int count_derivative = 0;
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    if (use_derivative[qi] && derivative_directions[qi] == 0)
                    {
                        derivative_first_order_x->computeDerivative(
                            data_derivative,
                            data_quantities[qi],
                            Real(dx[0]),
                            patch_visible_box,
                            count_derivative,
                            component_indices[qi]);
                        
                        count_derivative++;
                    }
                    else if (use_derivative[qi] && derivative_directions[qi] == 1)
                    {
                        derivative_first_order_y->computeDerivative(
                            data_derivative,
                            data_quantities[qi],
                            Real(dx[1]),
                            patch_visible_box,
                            count_derivative,
                            component_indices[qi]);
                        
                        count_derivative++;
                    }
                    else if (use_derivative[qi] && derivative_directions[qi] == 2)
                    {
                        derivative_first_order_z->computeDerivative(
                            data_derivative,
                            data_quantities[qi],
                            Real(dx[2]),
                            patch_visible_box,
                            count_derivative,
                            component_indices[qi]);
                        
                        count_derivative++;
                    }
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
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                            
                            int n_overlapped = 1;
                            
                            for (hier::BoxContainer::BoxContainerConstIterator iob(
                                    patch_overlapped_visible_boxes.begin());
                                 iob != patch_overlapped_visible_boxes.end();
                                 iob++)
                            {
                                const hier::Box& patch_overlapped_visible_box = *iob;
                                
                                if (patch_overlapped_visible_box.contains(idx_pt))
                                {
                                    n_overlapped++;
                                }
                            }
                            
                            /*
                             * Compute the linear index and the data to add.
                             */
                            
                            const int idx_avg = idx_lo_0 + i;
                            
                            Real avg = Real(1);
                            
                            count_derivative = 0;
                            for (int qi = 0; qi < num_quantities; qi++)
                            {
                                if (use_reciprocal[qi])
                                {
                                    if (use_derivative[qi])
                                    {
                                        const int idx_der = (relative_idx_lo_0 + i) +
                                            (relative_idx_lo_1 + j)*patch_interior_dim_0 +
                                            (relative_idx_lo_2 + k)*patch_interior_dim_0*
                                                patch_interior_dim_1;
                                        
                                        avg /= der_qi[count_derivative][idx_der];
                                        
                                        count_derivative++;
                                    }
                                    else
                                    {
                                        const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        avg /= u_qi[qi][idx_qi];
                                    }
                                }
                                else
                                {
                                    if (use_derivative[qi])
                                    {
                                        const int idx_der = (relative_idx_lo_0 + i) +
                                            (relative_idx_lo_1 + j)*patch_interior_dim_0 +
                                            (relative_idx_lo_2 + k)*patch_interior_dim_0*
                                                patch_interior_dim_1;
                                        
                                        avg *= der_qi[count_derivative][idx_der];
                                        
                                        count_derivative++;
                                    }
                                    else
                                    {
                                        const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        avg *= u_qi[qi][idx_qi];
                                    }
                                }
                            }
                            
                            avg_local[idx_avg] += (avg*weight/((Real) n_overlapped));
                        }
                    }
                }
            }
            
            /*
             * Unregister the patch and data of all registered derived cell variables in the flow model.
             */
            
            unregisterPatch();
        }
        
        /*
         * Reduction to get the global average.
         */
        
        d_mpi.Allreduce(
            avg_local,
            avg_global,
            coarsest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(avg_local);
    }
    
    return averaged_quantity;
}


/*
 * Compute averaged derivative of value (on product of variables) with only x direction as inhomogeneous direction
 * on the coarsest level.
 */
std::vector<Real>
FlowModelMPIHelperAverage::getAveragedDerivativeOfQuantityWithInhomogeneousXDirectionOnCoarsestLevel(
    const std::vector<std::string>& quantity_names,
    const std::vector<int>& component_indices,
    const std::vector<bool>& use_reciprocal,
    const int derivative_direction,
    const int num_ghosts_derivative,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    int num_quantities = static_cast<int>(quantity_names.size());
    
    TBOX_ASSERT(static_cast<int>(component_indices.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(use_reciprocal.size()) == num_quantities);
    
    if (d_dim == tbox::Dimension(1) && (derivative_direction < 0 || derivative_direction > 0))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperAverage::getAveragedDerivativeOfQuantityWithInhomogeneousXDirectionOnCoarsestLevel():\n"
            << "Cannot take derivative for one-dimensional problem!\n"
            << "derivative_direction = " << derivative_direction << " given!\n"
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2) && (derivative_direction < 0 || derivative_direction > 1))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperAverage::getAveragedDerivativeOfQuantityWithInhomogeneousXDirectionOnCoarsestLevel():\n"
            << "Cannot take derivative for two-dimensional problem!\n"
            << "derivative_direction = " << derivative_direction << " given!\n"
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(3) && (derivative_direction < 0 || derivative_direction > 2))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperAverage::getAveragedDerivativeOfQuantityWithInhomogeneousXDirectionOnCoarsestLevel():\n"
            << "Cannot take derivative for three-dimensional problem!\n"
            << "derivative_direction = " << derivative_direction << " given!\n"
            << std::endl);
    }
    
    hier::IntVector num_ghosts_der = hier::IntVector::getZero(d_dim);
    num_ghosts_der[derivative_direction] = num_ghosts_derivative;
    
    std::vector<Real> averaged_derivative;
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space. Coarsest and finest levels are both set to zero.
     */
    
    HAMERS_SHARED_PTR<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *d_patch_hierarchy,
            0,
            0));
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_x(
            new DerivativeFirstOrder(
                "first order derivative in x-direction",
                d_dim,
                DIRECTION::X_DIRECTION,
                num_ghosts_derivative));
        
        const int coarsest_level_dim_0 = d_coarsest_level_dims[0];
        
        Real* der_avg_local = (Real*)std::malloc(coarsest_level_dim_0*sizeof(Real));
        
        averaged_derivative.resize(coarsest_level_dim_0);
        Real* der_avg_global = averaged_derivative.data();
        
        for (int i = 0; i < coarsest_level_dim_0; i++)
        {
            der_avg_local[i]  = Real(0);
            der_avg_global[i] = Real(0);
        }
        
        /*
         * Get the coarsest patch level.
         */
        
        HAMERS_SHARED_PTR<hier::PatchLevel> patch_level_coarsest(
            d_patch_hierarchy->getPatchLevel(0));
        
        for (hier::PatchLevel::iterator ip(patch_level_coarsest->begin());
             ip != patch_level_coarsest->end();
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
             * Register the patch and the quantities in the flow model and compute the
             * corresponding product.
             */
            
            setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
            
            hier::IntVector num_ghosts = d_flow_model->getNumberOfGhostCells();
            TBOX_ASSERT(num_ghosts >= num_ghosts_der);
            
            std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
            
            for (int qi = 0; qi < num_quantities; qi++)
            {
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_names[qi], num_ghosts_der));
            }
            
            registerDerivedVariables(num_subghosts_of_data);
            
            allocateMemoryForDerivedCellData();
            
            computeDerivedCellData();
            
            /*
             * Get the pointers to data inside the flow model.
             */
            
            std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > data_quantities;
            data_quantities.resize(num_quantities);
            for (int qi = 0; qi < num_quantities; qi++)
            {
                data_quantities[qi] = getCellData(quantity_names[qi]);
            }
            
            std::vector<Real*> u_qi;
            u_qi.resize(num_quantities);
            for (int qi = 0; qi < num_quantities; qi++)
            {
                u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
            }
            
            const hier::BoxContainer& patch_visible_boxes =
                flattened_hierarchy->getVisibleBoxes(
                    patch_box,
                    0);
            
            const hier::BoxContainer& patch_overlapped_visible_boxes =
                flattened_hierarchy->getOverlappedVisibleBoxes(
                    patch_box,
                    0);
            
            std::vector<int> num_ghosts_0_u_qi;
            num_ghosts_0_u_qi.reserve(num_quantities);
            
            for (int qi = 0; qi < num_quantities; qi++)
            {
                const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
            }
            
            /*
             * Initialize cell data for the derivative and the product of variables.
             * Also, get pointers to the cell data containers.
             */
            
            HAMERS_SHARED_PTR<pdat::CellData<Real> > data_derivative(
                new pdat::CellData<Real>(patch_box, 1, hier::IntVector::getZero(d_dim)));
            
            HAMERS_SHARED_PTR<pdat::CellData<Real> > data_product(
                new pdat::CellData<Real>(patch_box, 1, num_ghosts_der));
            
            Real* der     = data_derivative->getPointer(0);
            Real* product = data_product->getPointer(0);
            
            const hier::IntVector num_ghosts_product = data_product->getGhostCellWidth();
            
            const int num_ghosts_0_product = num_ghosts_product[0];
            
            const hier::IntVector patch_interior_dims = patch_box.numberCells();
            const int patch_interior_dim_0 = patch_interior_dims[0];
            
            data_product->fillAll(Real(1));
            
            for (int qi = 0; qi < num_quantities; qi++)
            {
                if (use_reciprocal[qi])
                {
                    for (int i = -num_ghosts_0_product; i < patch_interior_dim_0 + num_ghosts_0_product; i++)
                    {
                        const int idx_product = i + num_ghosts_0_product;
                        const int idx_u_qi    = i + num_ghosts_0_u_qi[qi];
                        
                        product[idx_product] /= u_qi[qi][idx_u_qi];
                    }
                }
                else
                {
                    for (int i = -num_ghosts_0_product; i < patch_interior_dim_0 + num_ghosts_0_product; i++)
                    {
                        const int idx_product = i + num_ghosts_0_product;
                        const int idx_u_qi    = i + num_ghosts_0_u_qi[qi];
                        
                        product[idx_product] *= u_qi[qi][idx_u_qi];
                    }
                }
            }
            
            for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                 ib != patch_visible_boxes.end();
                 ib++)
            {
                const hier::Box& patch_visible_box = *ib;
                
                if (derivative_direction == 0)
                {
                    derivative_first_order_x->computeDerivative(
                        data_derivative,
                        data_product,
                        Real(dx[0]),
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
                     * Compute the index of the data point and count how many times the data is repeated.
                     */
                    
                    const hier::Index idx_pt(tbox::Dimension(1), idx_lo_0 + i);
                    
                    int n_overlapped = 1;
                    
                    for (hier::BoxContainer::BoxContainerConstIterator iob(
                            patch_overlapped_visible_boxes.begin());
                         iob != patch_overlapped_visible_boxes.end();
                         iob++)
                    {
                        const hier::Box& patch_overlapped_visible_box = *iob;
                        
                        if (patch_overlapped_visible_box.contains(idx_pt))
                        {
                            n_overlapped++;
                        }
                    }
                    
                    /*
                     * Compute the linear indices and the data to add.
                     */
                    
                    const int idx_avg = idx_lo_0 + i;
                    
                    const int idx = relative_idx_lo_0 + i;
                    
                    const Real value_to_add = der[idx]/((Real) n_overlapped);
                    
                    der_avg_local[idx_avg] += value_to_add;
                }
            }
            
            /*
             * Unregister the patch and data of all registered derived cell variables in the flow model.
             */
            
            unregisterPatch();
        }
        
        /*
         * Reduction to get the global average of derivative.
         */
        
        d_mpi.Allreduce(
            der_avg_local,
            der_avg_global,
            coarsest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(der_avg_local);
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
        
        const int coarsest_level_dim_0 = d_coarsest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_y = Real(x_hi[1] - x_lo[1]);
        
        Real* der_avg_local = (Real*)std::malloc(coarsest_level_dim_0*sizeof(Real));
        
        averaged_derivative.resize(coarsest_level_dim_0);
        Real* der_avg_global = averaged_derivative.data();
        
        for (int i = 0; i < coarsest_level_dim_0; i++)
        {
            der_avg_local[i]  = Real(0);
            der_avg_global[i] = Real(0);
        }
        
        /*
         * Get the coarsest patch level.
         */
        
        HAMERS_SHARED_PTR<hier::PatchLevel> patch_level_coarsest(
            d_patch_hierarchy->getPatchLevel(0));
        
        for (hier::PatchLevel::iterator ip(patch_level_coarsest->begin());
             ip != patch_level_coarsest->end();
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
             * Register the patch and the quantities in the flow model and compute the
             * corresponding product.
             */
            
            setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
            
            hier::IntVector num_ghosts = d_flow_model->getNumberOfGhostCells();
            TBOX_ASSERT(num_ghosts >= num_ghosts_der);
            
            std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
            
            for (int qi = 0; qi < num_quantities; qi++)
            {
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_names[qi], num_ghosts_der));
            }
            
            registerDerivedVariables(num_subghosts_of_data);
            
            allocateMemoryForDerivedCellData();
            
            computeDerivedCellData();
            
            /*
             * Get the pointers to data inside the flow model.
             */
            
            std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > data_quantities;
            data_quantities.resize(num_quantities);
            for (int qi = 0; qi < num_quantities; qi++)
            {
                data_quantities[qi] = getCellData(quantity_names[qi]);
            }
            
            std::vector<Real*> u_qi;
            u_qi.resize(num_quantities);
            for (int qi = 0; qi < num_quantities; qi++)
            {
                u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
            }
            
            const hier::BoxContainer& patch_visible_boxes =
                flattened_hierarchy->getVisibleBoxes(
                    patch_box,
                    0);
            
            const hier::BoxContainer& patch_overlapped_visible_boxes =
                flattened_hierarchy->getOverlappedVisibleBoxes(
                    patch_box,
                    0);
            
            std::vector<int> num_ghosts_0_u_qi;
            std::vector<int> num_ghosts_1_u_qi;
            std::vector<int> ghostcell_dim_0_u_qi;
            num_ghosts_0_u_qi.reserve(num_quantities);
            num_ghosts_1_u_qi.reserve(num_quantities);
            ghostcell_dim_0_u_qi.reserve(num_quantities);
            
            for (int qi = 0; qi < num_quantities; qi++)
            {
                const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                
                num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
            }
            
            /*
             * Initialize cell data for the derivative and the product of variables.
             * Also, get pointers to the cell data containers.
             */
            
            HAMERS_SHARED_PTR<pdat::CellData<Real> > data_derivative(
                new pdat::CellData<Real>(patch_box, 1, hier::IntVector::getZero(d_dim)));
            
            HAMERS_SHARED_PTR<pdat::CellData<Real> > data_product(
                new pdat::CellData<Real>(patch_box, 1, num_ghosts_der));
            
            Real* der     = data_derivative->getPointer(0);
            Real* product = data_product->getPointer(0);
            
            const hier::IntVector num_ghosts_product = data_product->getGhostCellWidth();
            const hier::IntVector ghostcell_dims_product = data_product->getGhostBox().numberCells();
            
            const int num_ghosts_0_product = num_ghosts_product[0];
            const int num_ghosts_1_product = num_ghosts_product[1];
            const int ghostcell_dim_0_product = ghostcell_dims_product[0];
            
            const hier::IntVector patch_interior_dims = patch_box.numberCells();
            const int patch_interior_dim_0 = patch_interior_dims[0];
            const int patch_interior_dim_1 = patch_interior_dims[1];
            
            data_product->fillAll(Real(1));
            
            for (int qi = 0; qi < num_quantities; qi++)
            {
                if (use_reciprocal[qi])
                {
                    for (int j = -num_ghosts_1_product; j < patch_interior_dim_1 + num_ghosts_1_product; j++)
                    {
                        for (int i = -num_ghosts_0_product; i < patch_interior_dim_0 + num_ghosts_0_product; i++)
                        {
                            const int idx_product = (i + num_ghosts_0_product) +
                                (j + num_ghosts_1_product)*ghostcell_dim_0_product;
                            
                            const int idx_u_qi = (i + num_ghosts_0_u_qi[qi]) +
                                (j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                
                            product[idx_product] /= u_qi[qi][idx_u_qi];
                        }
                    }
                }
                else
                {
                    for (int j = -num_ghosts_1_product; j < patch_interior_dim_1 + num_ghosts_1_product; j++)
                    {
                        for (int i = -num_ghosts_0_product; i < patch_interior_dim_0 + num_ghosts_0_product; i++)
                        {
                            const int idx_product = (i + num_ghosts_0_product) +
                                (j + num_ghosts_1_product)*ghostcell_dim_0_product;
                            
                            const int idx_u_qi = (i + num_ghosts_0_u_qi[qi]) +
                                (j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                
                            product[idx_product] *= u_qi[qi][idx_u_qi];
                        }
                    }
                }
            }
            
            const Real weight = Real(dx[1])/L_y;
            
            for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                 ib != patch_visible_boxes.end();
                 ib++)
            {
                const hier::Box& patch_visible_box = *ib;
                
                if (derivative_direction == 0)
                {
                    derivative_first_order_x->computeDerivative(
                        data_derivative,
                        data_product,
                        Real(dx[0]),
                        patch_visible_box,
                        0,
                        0);
                }
                else if (derivative_direction == 1)
                {
                    derivative_first_order_y->computeDerivative(
                        data_derivative,
                        data_product,
                        Real(dx[1]),
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
                         * Compute the index of the data point and count how many times the data is repeated.
                         */
                        
                        const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j);
                        
                        int n_overlapped = 1;
                        
                        for (hier::BoxContainer::BoxContainerConstIterator iob(
                                patch_overlapped_visible_boxes.begin());
                             iob != patch_overlapped_visible_boxes.end();
                             iob++)
                        {
                            const hier::Box& patch_overlapped_visible_box = *iob;
                            
                            if (patch_overlapped_visible_box.contains(idx_pt))
                            {
                                n_overlapped++;
                            }
                        }
                        
                        /*
                         * Compute the linear indices and the data to add.
                         */
                        
                        const int idx_avg = idx_lo_0 + i;
                        
                        const int idx = (relative_idx_lo_0 + i) +
                            (relative_idx_lo_1 + j)*patch_interior_dim_0;
                        
                        const Real value_to_add = der[idx]*weight/((Real) n_overlapped);
                        
                        der_avg_local[idx_avg] += value_to_add;
                    }
                }
            }
            
            /*
             * Unregister the patch and data of all registered derived cell variables in the flow model.
             */
            
            unregisterPatch();
        }
        
        /*
         * Reduction to get the global average of derivative.
         */
        
        d_mpi.Allreduce(
            der_avg_local,
            der_avg_global,
            coarsest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(der_avg_local);
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
        
        const int coarsest_level_dim_0 = d_coarsest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_y = Real(x_hi[1] - x_lo[1]);
        const Real L_z = Real(x_hi[2] - x_lo[2]);
        
        Real* der_avg_local = (Real*)std::malloc(coarsest_level_dim_0*sizeof(Real));
        
        averaged_derivative.resize(coarsest_level_dim_0);
        Real* der_avg_global = averaged_derivative.data();
        
        for (int i = 0; i < coarsest_level_dim_0; i++)
        {
            der_avg_local[i]  = Real(0);
            der_avg_global[i] = Real(0);
        }
        
        /*
         * Get the coarsest patch level.
         */
        
        HAMERS_SHARED_PTR<hier::PatchLevel> patch_level_coarsest(
            d_patch_hierarchy->getPatchLevel(0));
        
        for (hier::PatchLevel::iterator ip(patch_level_coarsest->begin());
             ip != patch_level_coarsest->end();
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
             * Register the patch and the quantities in the flow model and compute the
             * corresponding product.
             */
            
            setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
            
            hier::IntVector num_ghosts = d_flow_model->getNumberOfGhostCells();
            TBOX_ASSERT(num_ghosts >= num_ghosts_der);
            
            std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
            
            for (int qi = 0; qi < num_quantities; qi++)
            {
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_names[qi], num_ghosts));
            }
            
            registerDerivedVariables(num_subghosts_of_data);
            
            allocateMemoryForDerivedCellData();
            
            computeDerivedCellData();
            
            /*
             * Get the pointers to data inside the flow model.
             */
            
            std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > data_quantities;
            data_quantities.resize(num_quantities);
            for (int qi = 0; qi < num_quantities; qi++)
            {
                data_quantities[qi] = getCellData(quantity_names[qi]);
            }
            
            std::vector<Real*> u_qi;
            u_qi.resize(num_quantities);
            for (int qi = 0; qi < num_quantities; qi++)
            {
                u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
            }
            
            const hier::BoxContainer& patch_visible_boxes =
                flattened_hierarchy->getVisibleBoxes(
                    patch_box,
                    0);
            
            const hier::BoxContainer& patch_overlapped_visible_boxes =
                flattened_hierarchy->getOverlappedVisibleBoxes(
                    patch_box,
                    0);
            
            std::vector<int> num_ghosts_0_u_qi;
            std::vector<int> num_ghosts_1_u_qi;
            std::vector<int> num_ghosts_2_u_qi;
            std::vector<int> ghostcell_dim_0_u_qi;
            std::vector<int> ghostcell_dim_1_u_qi;
            num_ghosts_0_u_qi.reserve(num_quantities);
            num_ghosts_1_u_qi.reserve(num_quantities);
            num_ghosts_2_u_qi.reserve(num_quantities);
            ghostcell_dim_0_u_qi.reserve(num_quantities);
            ghostcell_dim_1_u_qi.reserve(num_quantities);
            
            for (int qi = 0; qi < num_quantities; qi++)
            {
                const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                
                num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                num_ghosts_2_u_qi.push_back(num_ghosts_u_qi[2]);
                ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                ghostcell_dim_1_u_qi.push_back(ghostcell_dims_u_qi[1]);
            }
            
            /*
             * Initialize cell data for the derivative and the product of variables.
             * Also, get pointers to the cell data containers.
             */
            
            HAMERS_SHARED_PTR<pdat::CellData<Real> > data_derivative(
                new pdat::CellData<Real>(patch_box, 1, hier::IntVector::getZero(d_dim)));
            
            HAMERS_SHARED_PTR<pdat::CellData<Real> > data_product(
                new pdat::CellData<Real>(patch_box, 1, num_ghosts_der));
            
            Real* der     = data_derivative->getPointer(0);
            Real* product = data_product->getPointer(0);
            
            const hier::IntVector num_ghosts_product = data_product->getGhostCellWidth();
            const hier::IntVector ghostcell_dims_product = data_product->getGhostBox().numberCells();
            
            const int num_ghosts_0_product = num_ghosts_product[0];
            const int num_ghosts_1_product = num_ghosts_product[1];
            const int num_ghosts_2_product = num_ghosts_product[2];
            const int ghostcell_dim_0_product = ghostcell_dims_product[0];
            const int ghostcell_dim_1_product = ghostcell_dims_product[1];
            
            const hier::IntVector patch_interior_dims = patch_box.numberCells();
            const int patch_interior_dim_0 = patch_interior_dims[0];
            const int patch_interior_dim_1 = patch_interior_dims[1];
            const int patch_interior_dim_2 = patch_interior_dims[2];
            
            data_product->fillAll(Real(1));
            
            for (int qi = 0; qi < num_quantities; qi++)
            {
                if (use_reciprocal[qi])
                {
                    for (int k = -num_ghosts_2_product; k < patch_interior_dim_2 + num_ghosts_2_product; k++)
                    {
                        for (int j = -num_ghosts_1_product; j < patch_interior_dim_1 + num_ghosts_1_product; j++)
                        {
                            for (int i = -num_ghosts_0_product; i < patch_interior_dim_0 + num_ghosts_0_product; i++)
                            {
                                const int idx_product = (i + num_ghosts_0_product) +
                                    (j + num_ghosts_1_product)*ghostcell_dim_0_product +
                                    (k + num_ghosts_2_product)*ghostcell_dim_0_product*
                                        ghostcell_dim_1_product;
                                
                                const int idx_u_qi = (i + num_ghosts_0_u_qi[qi]) +
                                    (j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                    (k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                        ghostcell_dim_1_u_qi[qi];
                                    
                                product[idx_product] /= u_qi[qi][idx_u_qi];
                            }
                        }
                    }
                }
                else
                {
                    for (int k = -num_ghosts_2_product; k < patch_interior_dim_2 + num_ghosts_2_product; k++)
                    {
                        for (int j = -num_ghosts_1_product; j < patch_interior_dim_1 + num_ghosts_1_product; j++)
                        {
                            for (int i = -num_ghosts_0_product; i < patch_interior_dim_0 + num_ghosts_0_product; i++)
                            {
                                const int idx_product = (i + num_ghosts_0_product) +
                                    (j + num_ghosts_1_product)*ghostcell_dim_0_product +
                                    (k + num_ghosts_2_product)*ghostcell_dim_0_product*
                                        ghostcell_dim_1_product;
                                
                                const int idx_u_qi = (i + num_ghosts_0_u_qi[qi]) +
                                    (j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                    (k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                        ghostcell_dim_1_u_qi[qi];
                                    
                                product[idx_product] *= u_qi[qi][idx_u_qi];
                            }
                        }
                    }
                }
            }
            
            const Real weight = Real(dx[1]*dx[2])/(L_y*L_z);
            
            for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                 ib != patch_visible_boxes.end();
                 ib++)
            {
                const hier::Box& patch_visible_box = *ib;
                
                if (derivative_direction == 0)
                {
                    derivative_first_order_x->computeDerivative(
                        data_derivative,
                        data_product,
                        Real(dx[0]),
                        patch_visible_box,
                        0,
                        0);
                }
                else if (derivative_direction == 1)
                {
                    derivative_first_order_y->computeDerivative(
                        data_derivative,
                        data_product,
                        Real(dx[1]),
                        patch_visible_box,
                        0,
                        0);
                }
                else if (derivative_direction == 2)
                {
                     derivative_first_order_z->computeDerivative(
                        data_derivative,
                        data_product,
                        Real(dx[2]),
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
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                            
                            int n_overlapped = 1;
                            
                            for (hier::BoxContainer::BoxContainerConstIterator iob(
                                    patch_overlapped_visible_boxes.begin());
                                 iob != patch_overlapped_visible_boxes.end();
                                 iob++)
                            {
                                const hier::Box& patch_overlapped_visible_box = *iob;
                                
                                if (patch_overlapped_visible_box.contains(idx_pt))
                                {
                                    n_overlapped++;
                                }
                            }
                            
                            /*
                             * Compute the linear indices and the data to add.
                             */
                            
                            const int idx_avg = idx_lo_0 + i;
                        
                            const int idx = (relative_idx_lo_0 + i) +
                                (relative_idx_lo_1 + j)*patch_interior_dim_0 +
                                (relative_idx_lo_2 + k)*patch_interior_dim_0*
                                    patch_interior_dim_1;
                            
                            const Real value_to_add = der[idx]*weight/((Real) n_overlapped);
                            
                            der_avg_local[idx_avg] += value_to_add;
                        }
                    }
                }
            }
            
            /*
             * Unregister the patch and data of all registered derived cell variables in the flow model.
             */
            
            unregisterPatch();
        }
        
        /*
         * Reduction to get the global average of derivative.
         */
        
        d_mpi.Allreduce(
            der_avg_local,
            der_avg_global,
            coarsest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(der_avg_local);
    }
    
    return averaged_derivative;
}


/*
 * Compute averaged value with only x-direction as inhomogeneous direction.
 */
std::vector<Real>
FlowModelMPIHelperAverage::getAveragedQuantityWithInhomogeneousXDirection(
    const std::string quantity_name,
    const int component_idx,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    std::vector<Real> averaged_quantity;
    
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
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        const int finest_level_dim_0 = d_finest_level_dims[0];
        
        Real* u_avg_local = (Real*)std::malloc(finest_level_dim_0*sizeof(Real));
        
        averaged_quantity.resize(finest_level_dim_0);
        Real* u_avg_global = averaged_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_avg_local[i]  = Real(0);
            u_avg_global[i] = Real(0);
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
            
            hier::IntVector ratio_to_coarsest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarsest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarsest_level/ratio_to_coarsest_level;
            
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
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity =
                    getCellData(quantity_name);
                
                Real* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
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
                         * Compute the index of the data point and count how many times the data is repeated.
                         */
                        
                        const hier::Index idx_pt(tbox::Dimension(1), idx_lo_0 + i);
                        
                        int n_overlapped = 1;
                        
                        for (hier::BoxContainer::BoxContainerConstIterator iob(
                                patch_overlapped_visible_boxes.begin());
                             iob != patch_overlapped_visible_boxes.end();
                             iob++)
                        {
                            const hier::Box& patch_overlapped_visible_box = *iob;
                            
                            if (patch_overlapped_visible_box.contains(idx_pt))
                            {
                                n_overlapped++;
                            }
                        }
                        
                        /*
                         * Compute the linear index and the data to add.
                         */
                        
                        const int idx = relative_idx_lo_0 + i + num_ghosts_0_quantity;
                        
                        const Real value_to_add = u[idx]/((Real) n_overlapped);
                        
                        for (int ii = 0; ii < ratio_to_finest_level_0; ii++)
                        {
                            const int idx_fine = (idx_lo_0 + i)*ratio_to_finest_level_0 + ii;
                            
                            u_avg_local[idx_fine] += value_to_add;
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        d_mpi.Allreduce(
            u_avg_local,
            u_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_avg_local);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_0 = d_finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_y = Real(x_hi[1] - x_lo[1]);
        
        Real* u_avg_local = (Real*)std::malloc(finest_level_dim_0*sizeof(Real));
        
        averaged_quantity.resize(finest_level_dim_0);
        Real* u_avg_global = averaged_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_avg_local[i]  = Real(0);
            u_avg_global[i] = Real(0);
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
            
            hier::IntVector ratio_to_coarsest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarsest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarsest_level/ratio_to_coarsest_level;
            
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
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity =
                    getCellData(quantity_name);
                
                Real* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                
                const Real weight = Real(dx[1])/L_y;
                
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
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j);
                            
                            int n_overlapped = 1;
                            
                            for (hier::BoxContainer::BoxContainerConstIterator iob(
                                    patch_overlapped_visible_boxes.begin());
                                 iob != patch_overlapped_visible_boxes.end();
                                 iob++)
                            {
                                const hier::Box& patch_overlapped_visible_box = *iob;
                                
                                if (patch_overlapped_visible_box.contains(idx_pt))
                                {
                                    n_overlapped++;
                                }
                            }
                            
                            /*
                             * Compute the linear index and the data to add.
                             */
                            
                            const int idx = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                            
                            const Real value_to_add = u[idx]*weight/((Real) n_overlapped);
                            
                            for (int ii = 0; ii < ratio_to_finest_level_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratio_to_finest_level_0 + ii;
                                
                                u_avg_local[idx_fine] += value_to_add;
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        d_mpi.Allreduce(
            u_avg_local,
            u_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_avg_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_0 = d_finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_y = Real(x_hi[1] - x_lo[1]);
        const Real L_z = Real(x_hi[2] - x_lo[2]);
        
        Real* u_avg_local = (Real*)std::malloc(finest_level_dim_0*sizeof(Real));
        
        averaged_quantity.resize(finest_level_dim_0);
        Real* u_avg_global = averaged_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_avg_local[i]  = Real(0);
            u_avg_global[i] = Real(0);
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
            
            hier::IntVector ratio_to_coarsest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarsest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarsest_level/ratio_to_coarsest_level;
            
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
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity =
                    getCellData(quantity_name);
                
                Real* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int num_ghosts_2_quantity = num_ghosts_quantity[2];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                const int ghostcell_dim_1_quantity = ghostcell_dims_quantity[1];
                
                const Real weight = Real(dx[1]*dx[2])/(L_y*L_z);
                
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
                                 * Compute the index of the data point and count how many times the data is repeated.
                                 */
                                
                                const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                                
                                int n_overlapped = 1;
                                
                                for (hier::BoxContainer::BoxContainerConstIterator iob(
                                        patch_overlapped_visible_boxes.begin());
                                     iob != patch_overlapped_visible_boxes.end();
                                     iob++)
                                {
                                    const hier::Box& patch_overlapped_visible_box = *iob;
                                    
                                    if (patch_overlapped_visible_box.contains(idx_pt))
                                    {
                                        n_overlapped++;
                                    }
                                }
                                
                                /*
                                 * Compute the linear index and the data to add.
                                 */
                                
                                const int idx = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                    (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                        ghostcell_dim_1_quantity;
                                
                                const Real value_to_add = u[idx]*weight/((Real) n_overlapped);
                                
                                for (int ii = 0; ii < ratio_to_finest_level_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratio_to_finest_level_0 + ii;
                                    
                                    u_avg_local[idx_fine] += value_to_add;
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        d_mpi.Allreduce(
            u_avg_local,
            u_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_avg_local);
    }
    
    return averaged_quantity;
}


/*
 * Compute averaged reciprocal of value with only x-direction as inhomogeneous direction.
 */
std::vector<Real>
FlowModelMPIHelperAverage::getAveragedReciprocalOfQuantityWithInhomogeneousXDirection(
    const std::string quantity_name,
    const int component_idx,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    std::vector<Real> averaged_reciprocal_quantity;
    
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
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        const int finest_level_dim_0 = d_finest_level_dims[0];
        
        Real* u_inv_avg_local = (Real*)std::malloc(finest_level_dim_0*sizeof(Real));
        
        averaged_reciprocal_quantity.resize(finest_level_dim_0);
        Real* u_inv_avg_global = averaged_reciprocal_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_inv_avg_local[i]  = Real(0);
            u_inv_avg_global[i] = Real(0);
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
            
            hier::IntVector ratio_to_coarsest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarsest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarsest_level/ratio_to_coarsest_level;
            
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
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity =
                    getCellData(quantity_name);
                
                Real* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
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
                         * Compute the index of the data point and count how many times the data is repeated.
                         */
                        
                        const hier::Index idx_pt(tbox::Dimension(1), idx_lo_0 + i);
                        
                        int n_overlapped = 1;
                        
                        for (hier::BoxContainer::BoxContainerConstIterator iob(
                                patch_overlapped_visible_boxes.begin());
                             iob != patch_overlapped_visible_boxes.end();
                             iob++)
                        {
                            const hier::Box& patch_overlapped_visible_box = *iob;
                            
                            if (patch_overlapped_visible_box.contains(idx_pt))
                            {
                                n_overlapped++;
                            }
                        }
                        
                        /*
                         * Compute the linear index and the data to add.
                         */
                        
                        const int idx = relative_idx_lo_0 + i + num_ghosts_0_quantity;
                        
                        const Real value_to_add = (Real(1)/u[idx])/((Real) n_overlapped);
                        
                        for (int ii = 0; ii < ratio_to_finest_level_0; ii++)
                        {
                            const int idx_fine = (idx_lo_0 + i)*ratio_to_finest_level_0 + ii;
                            
                            u_inv_avg_local[idx_fine] += value_to_add;
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average of reciprocal.
         */
        
        d_mpi.Allreduce(
            u_inv_avg_local,
            u_inv_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_inv_avg_local);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_0 = d_finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_y = Real(x_hi[1] - x_lo[1]);
        
        Real* u_inv_avg_local = (Real*)std::malloc(finest_level_dim_0*sizeof(Real));
        
        averaged_reciprocal_quantity.resize(finest_level_dim_0);
        Real* u_inv_avg_global = averaged_reciprocal_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_inv_avg_local[i]  = Real(0);
            u_inv_avg_global[i] = Real(0);
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
            
            hier::IntVector ratio_to_coarsest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarsest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarsest_level/ratio_to_coarsest_level;
            
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
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity =
                    getCellData(quantity_name);
                
                Real* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                
                const Real weight = Real(dx[1])/L_y;
                
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
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j);
                            
                            int n_overlapped = 1;
                            
                            for (hier::BoxContainer::BoxContainerConstIterator iob(
                                    patch_overlapped_visible_boxes.begin());
                                 iob != patch_overlapped_visible_boxes.end();
                                 iob++)
                            {
                                const hier::Box& patch_overlapped_visible_box = *iob;
                                
                                if (patch_overlapped_visible_box.contains(idx_pt))
                                {
                                    n_overlapped++;
                                }
                            }
                            
                            /*
                             * Compute the linear index and the data to add.
                             */
                            
                            const int idx = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                            
                            const Real value_to_add = (Real(1)/u[idx])*weight/((Real) n_overlapped);
                            
                            for (int ii = 0; ii < ratio_to_finest_level_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratio_to_finest_level_0 + ii;
                                
                                u_inv_avg_local[idx_fine] += value_to_add;
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average of reciprocal.
         */
        
        d_mpi.Allreduce(
            u_inv_avg_local,
            u_inv_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_inv_avg_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_0 = d_finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_y = Real(x_hi[1] - x_lo[1]);
        const Real L_z = Real(x_hi[2] - x_lo[2]);
        
        Real* u_inv_avg_local = (Real*)std::malloc(finest_level_dim_0*sizeof(Real));
        
        averaged_reciprocal_quantity.resize(finest_level_dim_0);
        Real* u_inv_avg_global = averaged_reciprocal_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_inv_avg_local[i]  = Real(0);
            u_inv_avg_global[i] = Real(0);
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
            
            hier::IntVector ratio_to_coarsest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarsest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarsest_level/ratio_to_coarsest_level;
            
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
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity =
                    getCellData(quantity_name);
                
                Real* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int num_ghosts_2_quantity = num_ghosts_quantity[2];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                const int ghostcell_dim_1_quantity = ghostcell_dims_quantity[1];
                
                const Real weight = Real(dx[1]*dx[2])/(L_y*L_z);
                
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
                                 * Compute the index of the data point and count how many times the data is repeated.
                                 */
                                
                                const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                                
                                int n_overlapped = 1;
                                
                                for (hier::BoxContainer::BoxContainerConstIterator iob(
                                        patch_overlapped_visible_boxes.begin());
                                     iob != patch_overlapped_visible_boxes.end();
                                     iob++)
                                {
                                    const hier::Box& patch_overlapped_visible_box = *iob;
                                    
                                    if (patch_overlapped_visible_box.contains(idx_pt))
                                    {
                                        n_overlapped++;
                                    }
                                }
                                
                                /*
                                 * Compute the linear index and the data to add.
                                 */
                                
                                const int idx = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                    (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                        ghostcell_dim_1_quantity;
                                
                                const Real value_to_add = (Real(1)/u[idx])*weight/((Real) n_overlapped);
                                
                                for (int ii = 0; ii < ratio_to_finest_level_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratio_to_finest_level_0 + ii;
                                    
                                    u_inv_avg_local[idx_fine] += value_to_add;
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average of reciprocal.
         */
        
        d_mpi.Allreduce(
            u_inv_avg_local,
            u_inv_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_inv_avg_local);
    }
    
    return averaged_reciprocal_quantity;
}


/*
 * Compute averaged value (on product of variables) with only x-direction as inhomogeneous direction.
 */
std::vector<Real>
FlowModelMPIHelperAverage::getAveragedQuantityWithInhomogeneousXDirection(
    const std::vector<std::string>& quantity_names,
    const std::vector<int>& component_indices,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    int num_quantities = static_cast<int>(quantity_names.size());
    
    TBOX_ASSERT(static_cast<int>(component_indices.size()) == num_quantities);
    
    std::vector<bool> use_reciprocal(num_quantities, false);
    
    return getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        use_reciprocal,
        data_context);
}


/*
 * Compute averaged value (on product of variables) with only x-direction as inhomogeneous direction.
 */
std::vector<Real>
FlowModelMPIHelperAverage::getAveragedQuantityWithInhomogeneousXDirection(
    const std::vector<std::string>& quantity_names,
    const std::vector<int>& component_indices,
    const std::vector<bool>& use_reciprocal,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    int num_quantities = static_cast<int>(quantity_names.size());
    
    TBOX_ASSERT(static_cast<int>(component_indices.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(use_reciprocal.size()) == num_quantities);
    
    std::vector<Real> averaged_quantity;
    
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
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        const int finest_level_dim_0 = d_finest_level_dims[0];
        
        Real* avg_local = (Real*)std::malloc(finest_level_dim_0*sizeof(Real));
        
        averaged_quantity.resize(finest_level_dim_0);
        Real* avg_global = averaged_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            avg_local[i]  = Real(0);
            avg_global[i] = Real(0);
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
            
            hier::IntVector ratio_to_coarsest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarsest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarsest_level/ratio_to_coarsest_level;
            
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
                 * Register the patch and data in the flow model and compute the corresponding
                 * average.
                 */
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                }
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = getCellData(quantity_names[qi]);
                }
                
                std::vector<Real*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                }
                
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
                         * Compute the index of the data point and count how many times the data is repeated.
                         */
                        
                        const hier::Index idx_pt(tbox::Dimension(1), idx_lo_0 + i);
                        
                        int n_overlapped = 1;
                        
                        for (hier::BoxContainer::BoxContainerConstIterator iob(
                                patch_overlapped_visible_boxes.begin());
                             iob != patch_overlapped_visible_boxes.end();
                             iob++)
                        {
                            const hier::Box& patch_overlapped_visible_box = *iob;
                            
                            if (patch_overlapped_visible_box.contains(idx_pt))
                            {
                                n_overlapped++;
                            }
                        }
                        
                        /*
                         * Compute the linear indices and the data to add.
                         */
                        
                        for (int ii = 0; ii < ratio_to_finest_level_0; ii++)
                        {
                            const int idx_fine = (idx_lo_0 + i)*ratio_to_finest_level_0 + ii;
                            
                            Real avg = Real(1);
                            
                            for (int qi = 0; qi < num_quantities; qi++)
                            {
                                const int idx_qi = relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi];
                                
                                if (use_reciprocal[qi])
                                {
                                    avg /= u_qi[qi][idx_qi];
                                }
                                else
                                {
                                    avg *= u_qi[qi][idx_qi];
                                }
                            }
                            
                            avg_local[idx_fine] += (avg/((Real) n_overlapped));
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        d_mpi.Allreduce(
            avg_local,
            avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(avg_local);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_0 = d_finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_y = Real(x_hi[1] - x_lo[1]);
        
        Real* avg_local = (Real*)std::malloc(finest_level_dim_0*sizeof(Real));
        
        averaged_quantity.resize(finest_level_dim_0);
        Real* avg_global = averaged_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            avg_local[i]  = Real(0);
            avg_global[i] = Real(0);
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
            
            hier::IntVector ratio_to_coarsest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarsest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarsest_level/ratio_to_coarsest_level;
            
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
                 * Register the patch and data in the flow model and compute the corresponding
                 * average.
                 */
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                }
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = getCellData(quantity_names[qi]);
                }
                
                std::vector<Real*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                }
                
                const Real weight = Real(dx[1])/L_y;
                
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
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j);
                            
                            int n_overlapped = 1;
                            
                            for (hier::BoxContainer::BoxContainerConstIterator iob(
                                    patch_overlapped_visible_boxes.begin());
                                 iob != patch_overlapped_visible_boxes.end();
                                 iob++)
                            {
                                const hier::Box& patch_overlapped_visible_box = *iob;
                                
                                if (patch_overlapped_visible_box.contains(idx_pt))
                                {
                                    n_overlapped++;
                                }
                            }
                            
                            /*
                             * Compute the linear indices and the data to add.
                             */
                            
                            for (int ii = 0; ii < ratio_to_finest_level_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratio_to_finest_level_0 + ii;
                                
                                Real avg = Real(1);
                                
                                for (int qi = 0; qi < num_quantities; qi++)
                                {
                                    const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    if (use_reciprocal[qi])
                                    {
                                        avg /= u_qi[qi][idx_qi];
                                    }
                                    else
                                    {
                                        avg *= u_qi[qi][idx_qi];
                                    }
                                }
                                
                                avg_local[idx_fine] += (avg*weight/((Real) n_overlapped));
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        d_mpi.Allreduce(
            avg_local,
            avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(avg_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_0 = d_finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_y = Real(x_hi[1] - x_lo[1]);
        const Real L_z = Real(x_hi[2] - x_lo[2]);
        
        Real* avg_local = (Real*)std::malloc(finest_level_dim_0*sizeof(Real));
        
        averaged_quantity.resize(finest_level_dim_0);
        Real* avg_global = averaged_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            avg_local[i]  = Real(0);
            avg_global[i] = Real(0);
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
            
            hier::IntVector ratio_to_coarsest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarsest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarsest_level/ratio_to_coarsest_level;
            
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
                 * Register the patch and data in the flow model and compute the corresponding
                 * average.
                 */
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                }
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = getCellData(quantity_names[qi]);
                }
                
                std::vector<Real*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> num_ghosts_2_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                std::vector<int> ghostcell_dim_1_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                num_ghosts_2_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                ghostcell_dim_1_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    num_ghosts_2_u_qi.push_back(num_ghosts_u_qi[2]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                    ghostcell_dim_1_u_qi.push_back(ghostcell_dims_u_qi[1]);
                }
                
                const Real weight = Real(dx[1]*dx[2])/(L_y*L_z);
                
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
                                 * Compute the index of the data point and count how many times the data is repeated.
                                 */
                                
                                const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                                
                                int n_overlapped = 1;
                                
                                for (hier::BoxContainer::BoxContainerConstIterator iob(
                                        patch_overlapped_visible_boxes.begin());
                                     iob != patch_overlapped_visible_boxes.end();
                                     iob++)
                                {
                                    const hier::Box& patch_overlapped_visible_box = *iob;
                                    
                                    if (patch_overlapped_visible_box.contains(idx_pt))
                                    {
                                        n_overlapped++;
                                    }
                                }
                                
                                /*
                                 * Compute the linear index and the data to add.
                                 */
                                
                                for (int ii = 0; ii < ratio_to_finest_level_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratio_to_finest_level_0 + ii;
                                    
                                    Real avg = Real(1);
                                    
                                    for (int qi = 0; qi < num_quantities; qi++)
                                    {
                                        const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        if (use_reciprocal[qi])
                                        {
                                            avg /= u_qi[qi][idx_qi];
                                        }
                                        else
                                        {
                                            avg *= u_qi[qi][idx_qi];
                                        }
                                    }
                                    
                                    avg_local[idx_fine] += (avg*weight/((Real) n_overlapped));
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        d_mpi.Allreduce(
            avg_local,
            avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(avg_local);
    }
    
    return averaged_quantity;
}


/*
 * Compute averaged value (on product of variable derivatives) with only x direction as inhomogeneous direction.
 */
std::vector<Real>
FlowModelMPIHelperAverage::getAveragedQuantityWithInhomogeneousXDirection(
    const std::vector<std::string>& quantity_names,
    const std::vector<int>& component_indices,
    const std::vector<bool>& use_derivative,
    const std::vector<int>& derivative_directions,
    const int num_ghosts_derivative,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    int num_quantities = static_cast<int>(quantity_names.size());
    
    TBOX_ASSERT(static_cast<int>(component_indices.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(use_derivative.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(derivative_directions.size()) == num_quantities);
    
    std::vector<bool> use_reciprocal(num_quantities, false);
    
    return getAveragedQuantityWithInhomogeneousXDirection(
        quantity_names,
        component_indices,
        use_derivative,
        derivative_directions,
        use_reciprocal,
        num_ghosts_derivative,
        data_context);
}


/*
 * Compute averaged value (on product of variable derivatives) with only x direction as inhomogeneous direction.
 */
std::vector<Real>
FlowModelMPIHelperAverage::getAveragedQuantityWithInhomogeneousXDirection(
    const std::vector<std::string>& quantity_names,
    const std::vector<int>& component_indices,
    const std::vector<bool>& use_derivative,
    const std::vector<int>& derivative_directions,
    const std::vector<bool>& use_reciprocal,
    const int num_ghosts_derivative,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    int num_quantities = static_cast<int>(quantity_names.size());
    
    TBOX_ASSERT(static_cast<int>(component_indices.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(use_derivative.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(derivative_directions.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(use_reciprocal.size()) == num_quantities);
    
    int num_use_derivative = 0;
    
    if (d_dim == tbox::Dimension(1))
    {
        for (int qi = 0; qi < num_quantities; qi++)
        {
            if (use_derivative[qi])
            {
                num_use_derivative++;
                if (derivative_directions[qi] < 0 || derivative_directions[qi] > 0)
                {
                    TBOX_ERROR(d_object_name
                        << ": FlowModelMPIHelperAverage::getAveragedQuantityWithInhomogeneousXDirection():\n"
                        << "Cannot take derivative for one-dimensional problem!\n"
                        << "derivative_directions[" << qi << "] = " << derivative_directions[qi] << " given!\n"
                        << std::endl);
                }
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        for (int qi = 0; qi < num_quantities; qi++)
        {
            if (use_derivative[qi])
            {
                num_use_derivative++;
                if (derivative_directions[qi] < 0 || derivative_directions[qi] > 1)
                {
                    TBOX_ERROR(d_object_name
                        << ": FlowModelMPIHelperAverage::getAveragedQuantityWithInhomogeneousXDirection():\n"
                        << "Cannot take derivative for two-dimensional problem!\n"
                        << "derivative_directions[" << qi << "] = " << derivative_directions[qi] << " given!\n"
                        << std::endl);
                }
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        for (int qi = 0; qi < num_quantities; qi++)
        {
            if (use_derivative[qi])
            {
                num_use_derivative++;
                if (derivative_directions[qi] < 0 || derivative_directions[qi] > 2)
                {
                    TBOX_ERROR(d_object_name
                        << ": FlowModelMPIHelperAverage::getAveragedQuantityWithInhomogeneousXDirection():\n"
                        << "Cannot take derivative for three-dimensional problem!\n"
                        << "derivative_directions[" << qi << "] = " << derivative_directions[qi] << " given!\n"
                        << std::endl);
                }
            }
        }
    }
    
    std::vector<Real> averaged_quantity;
    
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
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_x(
            new DerivativeFirstOrder(
                "first order derivative in x-direction",
                d_dim,
                DIRECTION::X_DIRECTION,
                num_ghosts_derivative));
        
        hier::IntVector num_ghosts_der = hier::IntVector::getOne(d_dim)*num_ghosts_derivative;
        
        const int finest_level_dim_0 = d_finest_level_dims[0];
        
        Real* avg_local = (Real*)std::malloc(finest_level_dim_0*sizeof(Real));
        
        averaged_quantity.resize(finest_level_dim_0);
        Real* avg_global = averaged_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            avg_local[i]  = Real(0);
            avg_global[i] = Real(0);
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
            
            hier::IntVector ratio_to_coarsest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarsest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarsest_level/ratio_to_coarsest_level;
            
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
                 * Register the patch and data in the flow model and compute the corresponding
                 * average.
                 */
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    if (num_use_derivative > 0)
                    {
                        num_subghosts_of_data.insert(
                            std::pair<std::string, hier::IntVector>(quantity_names[qi], num_ghosts_der));
                    }
                    else
                    {
                        num_subghosts_of_data.insert(
                            std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                    }
                }
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = getCellData(quantity_names[qi]);
                }
                
                std::vector<Real*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                }
                
                /*
                 * Initialize cell data for the derivatives and get pointers to the cell data.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_derivative;
                std::vector<Real*> der_qi;
                
                if (num_use_derivative > 0)
                {
                    data_derivative = HAMERS_MAKE_SHARED<pdat::CellData<Real> >(patch_box, num_use_derivative, hier::IntVector::getZero(d_dim));
                    
                    der_qi.resize(num_use_derivative);
                    for (int qi = 0; qi < num_use_derivative; qi++)
                    {
                        der_qi[qi] = data_derivative->getPointer(qi);
                    }
                }
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    int count_derivative = 0;
                    for (int qi = 0; qi < num_quantities; qi++)
                    {
                        if (use_derivative[qi] && derivative_directions[qi] == 0)
                        {
                            derivative_first_order_x->computeDerivative(
                                data_derivative,
                                data_quantities[qi],
                                Real(dx[0]),
                                patch_visible_box,
                                count_derivative,
                                component_indices[qi]);
                            
                            count_derivative++;
                        }
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
                         * Compute the index of the data point and count how many times the data is repeated.
                         */
                        
                        const hier::Index idx_pt(tbox::Dimension(1), idx_lo_0 + i);
                        
                        int n_overlapped = 1;
                        
                        for (hier::BoxContainer::BoxContainerConstIterator iob(
                                patch_overlapped_visible_boxes.begin());
                             iob != patch_overlapped_visible_boxes.end();
                             iob++)
                        {
                            const hier::Box& patch_overlapped_visible_box = *iob;
                            
                            if (patch_overlapped_visible_box.contains(idx_pt))
                            {
                                n_overlapped++;
                            }
                        }
                        
                        /*
                         * Compute the linear indices and the data to add.
                         */
                        
                        Real avg = Real(1);
                        
                        count_derivative = 0;
                        for (int qi = 0; qi < num_quantities; qi++)
                        {
                            if (use_reciprocal[qi])
                            {
                                if (use_derivative[qi])
                                {
                                    const int idx_der = relative_idx_lo_0 + i;
                                    
                                    avg /= der_qi[count_derivative][idx_der];
                                    count_derivative++;
                                }
                                else
                                {
                                    const int idx_qi = relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi];
                                    
                                    avg /= u_qi[qi][idx_qi];
                                }
                            }
                            else
                            {
                                if (use_derivative[qi])
                                {
                                    const int idx_der = relative_idx_lo_0 + i;
                                    
                                    avg *= der_qi[count_derivative][idx_der];
                                    count_derivative++;
                                }
                                else
                                {
                                    const int idx_qi = relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi];
                                    
                                    avg *= u_qi[qi][idx_qi];
                                }
                            }
                        }
                        
                        /*
                         * Add the data.
                         */
                        
                        for (int ii = 0; ii < ratio_to_finest_level_0; ii++)
                        {
                            const int idx_fine = (idx_lo_0 + i)*ratio_to_finest_level_0 + ii;
                            
                            avg_local[idx_fine] += (avg/((Real) n_overlapped));
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        d_mpi.Allreduce(
            avg_local,
            avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(avg_local);
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
        
        hier::IntVector num_ghosts_der = hier::IntVector::getOne(d_dim)*num_ghosts_derivative;
        
        const int finest_level_dim_0 = d_finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_y = Real(x_hi[1] - x_lo[1]);
        
        Real* avg_local = (Real*)std::malloc(finest_level_dim_0*sizeof(Real));
        
        averaged_quantity.resize(finest_level_dim_0);
        Real* avg_global = averaged_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            avg_local[i]  = Real(0);
            avg_global[i] = Real(0);
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
            
            hier::IntVector ratio_to_coarsest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarsest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarsest_level/ratio_to_coarsest_level;
            
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
                 * Register the patch and data in the flow model and compute the corresponding
                 * average.
                 */
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    if (num_use_derivative > 0)
                    {
                        num_subghosts_of_data.insert(
                            std::pair<std::string, hier::IntVector>(quantity_names[qi], num_ghosts_der));
                    }
                    else
                    {
                        num_subghosts_of_data.insert(
                            std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                    }
                }
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = getCellData(quantity_names[qi]);
                }
                
                std::vector<Real*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                }
                
                /*
                 * Initialize cell data for the derivatives and get pointers to the cell data.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_derivative;
                std::vector<Real*> der_qi;
                
                if (num_use_derivative > 0)
                {
                    data_derivative = HAMERS_MAKE_SHARED<pdat::CellData<Real> >(patch_box, num_use_derivative, hier::IntVector::getZero(d_dim));
                    
                    der_qi.resize(num_use_derivative);
                    for (int qi = 0; qi < num_use_derivative; qi++)
                    {
                        der_qi[qi] = data_derivative->getPointer(qi);
                    }
                }
                
                const hier::IntVector patch_interior_dims = patch_box.numberCells();
                const int patch_interior_dim_0 = patch_interior_dims[0];
                
                const Real weight = Real(dx[1])/L_y;
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    int count_derivative = 0;
                    for (int qi = 0; qi < num_quantities; qi++)
                    {
                        if (use_derivative[qi] && derivative_directions[qi] == 0)
                        {
                            derivative_first_order_x->computeDerivative(
                                data_derivative,
                                data_quantities[qi],
                                Real(dx[0]),
                                patch_visible_box,
                                count_derivative,
                                component_indices[qi]);
                            
                            count_derivative++;
                        }
                        else if (use_derivative[qi] && derivative_directions[qi] == 1)
                        {
                            derivative_first_order_y->computeDerivative(
                                data_derivative,
                                data_quantities[qi],
                                Real(dx[1]),
                                patch_visible_box,
                                count_derivative,
                                component_indices[qi]);
                            
                            count_derivative++;
                        }
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
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j);
                            
                            int n_overlapped = 1;
                            
                            for (hier::BoxContainer::BoxContainerConstIterator iob(
                                    patch_overlapped_visible_boxes.begin());
                                 iob != patch_overlapped_visible_boxes.end();
                                 iob++)
                            {
                                const hier::Box& patch_overlapped_visible_box = *iob;
                                
                                if (patch_overlapped_visible_box.contains(idx_pt))
                                {
                                    n_overlapped++;
                                }
                            }
                            
                            /*
                             * Compute the linear indices and the data to add.
                             */
                            
                            Real avg = Real(1);
                            
                            count_derivative = 0;
                            for (int qi = 0; qi < num_quantities; qi++)
                            {
                                if (use_reciprocal[qi])
                                {
                                    if (use_derivative[qi])
                                    {
                                        const int idx_der = (relative_idx_lo_0 + i) +
                                            (relative_idx_lo_1 + j)*patch_interior_dim_0;
                                        
                                        avg /= der_qi[count_derivative][idx_der];
                                        
                                        count_derivative++;
                                    }
                                    else
                                    {
                                        const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                        
                                        avg /= u_qi[qi][idx_qi];
                                    }
                                }
                                else
                                {
                                    if (use_derivative[qi])
                                    {
                                        const int idx_der = (relative_idx_lo_0 + i) +
                                            (relative_idx_lo_1 + j)*patch_interior_dim_0;
                                        
                                        avg *= der_qi[count_derivative][idx_der];
                                        
                                        count_derivative++;
                                    }
                                    else
                                    {
                                        const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                        
                                        avg *= u_qi[qi][idx_qi];
                                    }
                                }
                            }
                            
                            /*
                             * Add the data.
                             */
                            
                            for (int ii = 0; ii < ratio_to_finest_level_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratio_to_finest_level_0 + ii;
                                
                                avg_local[idx_fine] += (avg*weight/((Real) n_overlapped));
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        d_mpi.Allreduce(
            avg_local,
            avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(avg_local);
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
        
        hier::IntVector num_ghosts_der = hier::IntVector::getOne(d_dim)*num_ghosts_derivative;
        
        const int finest_level_dim_0 = d_finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_y = Real(x_hi[1] - x_lo[1]);
        const Real L_z = Real(x_hi[2] - x_lo[2]);
        
        Real* avg_local = (Real*)std::malloc(finest_level_dim_0*sizeof(Real));
        
        averaged_quantity.resize(finest_level_dim_0);
        Real* avg_global = averaged_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            avg_local[i]  = Real(0);
            avg_global[i] = Real(0);
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
            
            hier::IntVector ratio_to_coarsest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarsest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarsest_level/ratio_to_coarsest_level;
            
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
                 * Register the patch and data in the flow model and compute the corresponding
                 * average.
                 */
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    if (num_use_derivative > 0)
                    {
                        num_subghosts_of_data.insert(
                            std::pair<std::string, hier::IntVector>(quantity_names[qi], num_ghosts_der));
                    }
                    else
                    {
                        num_subghosts_of_data.insert(
                            std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                    }
                }
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = getCellData(quantity_names[qi]);
                }
                
                std::vector<Real*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> num_ghosts_2_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                std::vector<int> ghostcell_dim_1_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                num_ghosts_2_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                ghostcell_dim_1_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    num_ghosts_2_u_qi.push_back(num_ghosts_u_qi[2]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                    ghostcell_dim_1_u_qi.push_back(ghostcell_dims_u_qi[1]);
                }
                
                /*
                 * Initialize cell data for the derivatives and get pointers to the cell data.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_derivative;
                std::vector<Real*> der_qi;
                
                if (num_use_derivative > 0)
                {
                    data_derivative = HAMERS_MAKE_SHARED<pdat::CellData<Real> >(patch_box, num_use_derivative, hier::IntVector::getZero(d_dim));
                    
                    der_qi.resize(num_use_derivative);
                    for (int qi = 0; qi < num_use_derivative; qi++)
                    {
                        der_qi[qi] = data_derivative->getPointer(qi);
                    }
                }
                
                const hier::IntVector patch_interior_dims = patch_box.numberCells();
                const int patch_interior_dim_0 = patch_interior_dims[0];
                const int patch_interior_dim_1 = patch_interior_dims[1];
                
                const Real weight = Real(dx[1]*dx[2])/(L_y*L_z);
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    int count_derivative = 0;
                    for (int qi = 0; qi < num_quantities; qi++)
                    {
                        if (use_derivative[qi] && derivative_directions[qi] == 0)
                        {
                            derivative_first_order_x->computeDerivative(
                                data_derivative,
                                data_quantities[qi],
                                Real(dx[0]),
                                patch_visible_box,
                                count_derivative,
                                component_indices[qi]);
                            
                            count_derivative++;
                        }
                        else if (use_derivative[qi] && derivative_directions[qi] == 1)
                        {
                            derivative_first_order_y->computeDerivative(
                                data_derivative,
                                data_quantities[qi],
                                Real(dx[1]),
                                patch_visible_box,
                                count_derivative,
                                component_indices[qi]);
                            
                            count_derivative++;
                        }
                        else if (use_derivative[qi] && derivative_directions[qi] == 2)
                        {
                            derivative_first_order_z->computeDerivative(
                                data_derivative,
                                data_quantities[qi],
                                Real(dx[2]),
                                patch_visible_box,
                                count_derivative,
                                component_indices[qi]);
                            
                            count_derivative++;
                        }
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
                                 * Compute the index of the data point and count how many times the data is repeated.
                                 */
                                
                                const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                                
                                int n_overlapped = 1;
                                
                                for (hier::BoxContainer::BoxContainerConstIterator iob(
                                        patch_overlapped_visible_boxes.begin());
                                     iob != patch_overlapped_visible_boxes.end();
                                     iob++)
                                {
                                    const hier::Box& patch_overlapped_visible_box = *iob;
                                    
                                    if (patch_overlapped_visible_box.contains(idx_pt))
                                    {
                                        n_overlapped++;
                                    }
                                }
                                
                                /*
                                 * Compute the linear index and the data to add.
                                 */
                                
                                Real avg = Real(1);
                                
                                count_derivative = 0;
                                for (int qi = 0; qi < num_quantities; qi++)
                                {
                                    if (use_reciprocal[qi])
                                    {
                                        if (use_derivative[qi])
                                        {
                                            const int idx_der = (relative_idx_lo_0 + i) +
                                                (relative_idx_lo_1 + j)*patch_interior_dim_0 +
                                                (relative_idx_lo_2 + k)*patch_interior_dim_0*
                                                    patch_interior_dim_1;
                                            
                                            avg /= der_qi[count_derivative][idx_der];
                                            
                                            count_derivative++;
                                        }
                                        else
                                        {
                                            const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                                (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                                (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                    ghostcell_dim_1_u_qi[qi];
                                            
                                            avg /= u_qi[qi][idx_qi];
                                        }
                                    }
                                    else
                                    {
                                        if (use_derivative[qi])
                                        {
                                            const int idx_der = (relative_idx_lo_0 + i) +
                                                (relative_idx_lo_1 + j)*patch_interior_dim_0 +
                                                (relative_idx_lo_2 + k)*patch_interior_dim_0*
                                                    patch_interior_dim_1;
                                            
                                            avg *= der_qi[count_derivative][idx_der];
                                            
                                            count_derivative++;
                                        }
                                        else
                                        {
                                            const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                                (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                                (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                    ghostcell_dim_1_u_qi[qi];
                                            
                                            avg *= u_qi[qi][idx_qi];
                                        }
                                    }
                                }
                                
                                /*
                                 * Add the data.
                                 */
                                
                                for (int ii = 0; ii < ratio_to_finest_level_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratio_to_finest_level_0 + ii;
                                    
                                    avg_local[idx_fine] += (avg*weight/((Real) n_overlapped));
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        d_mpi.Allreduce(
            avg_local,
            avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(avg_local);
    }
    
    return averaged_quantity;
}


/*
 * Compute averaged derivative of value (on product of variables) with only x direction as inhomogeneous direction.
 */
std::vector<Real>
FlowModelMPIHelperAverage::getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
    const std::vector<std::string>& quantity_names,
    const std::vector<int>& component_indices,
    const std::vector<bool>& use_reciprocal,
    const int derivative_direction,
    const int num_ghosts_derivative,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    int num_quantities = static_cast<int>(quantity_names.size());
    
    TBOX_ASSERT(static_cast<int>(component_indices.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(use_reciprocal.size()) == num_quantities);
    
    if (d_dim == tbox::Dimension(1) && (derivative_direction < 0 || derivative_direction > 0))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperAverage::getAveragedDerivativeOfQuantityWithInhomogeneousXDirection():\n"
            << "Cannot take derivative for one-dimensional problem!\n"
            << "derivative_direction = " << derivative_direction << " given!\n"
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2) && (derivative_direction < 0 || derivative_direction > 1))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperAverage::getAveragedDerivativeOfQuantityWithInhomogeneousXDirection():\n"
            << "Cannot take derivative for two-dimensional problem!\n"
            << "derivative_direction = " << derivative_direction << " given!\n"
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(3) && (derivative_direction < 0 || derivative_direction > 2))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperAverage::getAveragedDerivativeOfQuantityWithInhomogeneousXDirection():\n"
            << "Cannot take derivative for three-dimensional problem!\n"
            << "derivative_direction = " << derivative_direction << " given!\n"
            << std::endl);
    }
    
    hier::IntVector num_ghosts_der = hier::IntVector::getZero(d_dim);
    num_ghosts_der[derivative_direction] = num_ghosts_derivative;
    
    std::vector<Real> averaged_derivative;
    
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
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_x(
            new DerivativeFirstOrder(
                "first order derivative in x-direction",
                d_dim,
                DIRECTION::X_DIRECTION,
                num_ghosts_derivative));
        
        const int finest_level_dim_0 = d_finest_level_dims[0];
        
        Real* der_avg_local = (Real*)std::malloc(finest_level_dim_0*sizeof(Real));
        
        averaged_derivative.resize(finest_level_dim_0);
        Real* der_avg_global = averaged_derivative.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            der_avg_local[i]  = Real(0);
            der_avg_global[i] = Real(0);
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
            
            hier::IntVector ratio_to_coarsest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarsest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarsest_level/ratio_to_coarsest_level;
            
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
                 * Register the patch and the quantities in the flow model and compute the
                 * corresponding product.
                 */
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model->getNumberOfGhostCells();
                TBOX_ASSERT(num_ghosts >= num_ghosts_der);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], num_ghosts_der));
                }
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = getCellData(quantity_names[qi]);
                }
                
                std::vector<Real*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                }
                
                /*
                 * Initialize cell data for the derivative and the product of variables.
                 * Also, get pointers to the cell data containers.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_derivative(
                    new pdat::CellData<Real>(patch_box, 1, hier::IntVector::getZero(d_dim)));
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_product(
                    new pdat::CellData<Real>(patch_box, 1, num_ghosts_der));
                
                Real* der     = data_derivative->getPointer(0);
                Real* product = data_product->getPointer(0);
                
                const hier::IntVector num_ghosts_product = data_product->getGhostCellWidth();
                
                const int num_ghosts_0_product = num_ghosts_product[0];
                
                const hier::IntVector patch_interior_dims = patch_box.numberCells();
                const int patch_interior_dim_0 = patch_interior_dims[0];
                
                data_product->fillAll(Real(1));
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    if (use_reciprocal[qi])
                    {
                        for (int i = -num_ghosts_0_product; i < patch_interior_dim_0 + num_ghosts_0_product; i++)
                        {
                            const int idx_product = i + num_ghosts_0_product;
                            const int idx_u_qi    = i + num_ghosts_0_u_qi[qi];
                            
                            product[idx_product] /= u_qi[qi][idx_u_qi];
                        }
                    }
                    else
                    {
                        for (int i = -num_ghosts_0_product; i < patch_interior_dim_0 + num_ghosts_0_product; i++)
                        {
                            const int idx_product = i + num_ghosts_0_product;
                            const int idx_u_qi    = i + num_ghosts_0_u_qi[qi];
                            
                            product[idx_product] *= u_qi[qi][idx_u_qi];
                        }
                    }
                }
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    if (derivative_direction == 0)
                    {
                        derivative_first_order_x->computeDerivative(
                            data_derivative,
                            data_product,
                            Real(dx[0]),
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
                         * Compute the index of the data point and count how many times the data is repeated.
                         */
                        
                        const hier::Index idx_pt(tbox::Dimension(1), idx_lo_0 + i);
                        
                        int n_overlapped = 1;
                        
                        for (hier::BoxContainer::BoxContainerConstIterator iob(
                                patch_overlapped_visible_boxes.begin());
                             iob != patch_overlapped_visible_boxes.end();
                             iob++)
                        {
                            const hier::Box& patch_overlapped_visible_box = *iob;
                            
                            if (patch_overlapped_visible_box.contains(idx_pt))
                            {
                                n_overlapped++;
                            }
                        }
                        
                        /*
                         * Compute the linear index and the data to add.
                         */
                        
                        const int idx = relative_idx_lo_0 + i;
                        
                        const Real value_to_add = der[idx]/((Real) n_overlapped);
                        
                        /*
                         * Add the data.
                         */
                        
                        for (int ii = 0; ii < ratio_to_finest_level_0; ii++)
                        {
                            const int idx_fine = (idx_lo_0 + i)*ratio_to_finest_level_0 + ii;
                            
                            der_avg_local[idx_fine] += value_to_add;
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average of derivative.
         */
        
        d_mpi.Allreduce(
            der_avg_local,
            der_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(der_avg_local);
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
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_y = Real(x_hi[1] - x_lo[1]);
        
        Real* der_avg_local = (Real*)std::malloc(finest_level_dim_0*sizeof(Real));
        
        averaged_derivative.resize(finest_level_dim_0);
        Real* der_avg_global = averaged_derivative.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            der_avg_local[i]  = Real(0);
            der_avg_global[i] = Real(0);
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
            
            hier::IntVector ratio_to_coarsest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarsest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarsest_level/ratio_to_coarsest_level;
            
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
                 * Register the patch and the quantities in the flow model and compute the
                 * corresponding product.
                 */
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model->getNumberOfGhostCells();
                TBOX_ASSERT(num_ghosts >= num_ghosts_der);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], num_ghosts_der));
                }
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = getCellData(quantity_names[qi]);
                }
                
                std::vector<Real*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                }
                
                /*
                 * Initialize cell data for the derivative and the product of variables.
                 * Also, get pointers to the cell data containers.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_derivative(
                    new pdat::CellData<Real>(patch_box, 1, hier::IntVector::getZero(d_dim)));
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_product(
                    new pdat::CellData<Real>(patch_box, 1, num_ghosts_der));
                
                Real* der     = data_derivative->getPointer(0);
                Real* product = data_product->getPointer(0);
                
                const hier::IntVector num_ghosts_product = data_product->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_product = data_product->getGhostBox().numberCells();
                
                const int num_ghosts_0_product = num_ghosts_product[0];
                const int num_ghosts_1_product = num_ghosts_product[1];
                const int ghostcell_dim_0_product = ghostcell_dims_product[0];
                
                const hier::IntVector patch_interior_dims = patch_box.numberCells();
                const int patch_interior_dim_0 = patch_interior_dims[0];
                const int patch_interior_dim_1 = patch_interior_dims[1];
                
                data_product->fillAll(Real(1));
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    if (use_reciprocal[qi])
                    {
                        for (int j = -num_ghosts_1_product; j < patch_interior_dim_1 + num_ghosts_1_product; j++)
                        {
                            for (int i = -num_ghosts_0_product; i < patch_interior_dim_0 + num_ghosts_0_product; i++)
                            {
                                const int idx_product = (i + num_ghosts_0_product) +
                                    (j + num_ghosts_1_product)*ghostcell_dim_0_product;
                                
                                const int idx_u_qi = (i + num_ghosts_0_u_qi[qi]) +
                                    (j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                product[idx_product] /= u_qi[qi][idx_u_qi];
                            }
                        }
                    }
                    else
                    {
                        for (int j = -num_ghosts_1_product; j < patch_interior_dim_1 + num_ghosts_1_product; j++)
                        {
                            for (int i = -num_ghosts_0_product; i < patch_interior_dim_0 + num_ghosts_0_product; i++)
                            {
                                const int idx_product = (i + num_ghosts_0_product) +
                                    (j + num_ghosts_1_product)*ghostcell_dim_0_product;
                                
                                const int idx_u_qi = (i + num_ghosts_0_u_qi[qi]) +
                                    (j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                product[idx_product] *= u_qi[qi][idx_u_qi];
                            }
                        }
                    }
                }
                
                const Real weight = Real(dx[1])/L_y;
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    if (derivative_direction == 0)
                    {
                        derivative_first_order_x->computeDerivative(
                            data_derivative,
                            data_product,
                            Real(dx[0]),
                            patch_visible_box,
                            0,
                            0);
                    }
                    else if (derivative_direction == 1)
                    {
                        derivative_first_order_y->computeDerivative(
                            data_derivative,
                            data_product,
                            Real(dx[1]),
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
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j);
                            
                            int n_overlapped = 1;
                            
                            for (hier::BoxContainer::BoxContainerConstIterator iob(
                                    patch_overlapped_visible_boxes.begin());
                                 iob != patch_overlapped_visible_boxes.end();
                                 iob++)
                            {
                                const hier::Box& patch_overlapped_visible_box = *iob;
                                
                                if (patch_overlapped_visible_box.contains(idx_pt))
                                {
                                    n_overlapped++;
                                }
                            }
                            
                            /*
                             * Compute the linear index and the data to add.
                             */
                            
                            const int idx = (relative_idx_lo_0 + i) +
                                (relative_idx_lo_1 + j)*patch_interior_dim_0;
                            
                            const Real value_to_add = der[idx]*weight/((Real) n_overlapped);
                            
                            /*
                             * Add the data.
                             */
                            
                            for (int ii = 0; ii < ratio_to_finest_level_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratio_to_finest_level_0 + ii;
                                
                                der_avg_local[idx_fine] += value_to_add;
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average of derivative.
         */
        
        d_mpi.Allreduce(
            der_avg_local,
            der_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(der_avg_local);
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
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_y = Real(x_hi[1] - x_lo[1]);
        const Real L_z = Real(x_hi[2] - x_lo[2]);
        
        Real* der_avg_local = (Real*)std::malloc(finest_level_dim_0*sizeof(Real));
        
        averaged_derivative.resize(finest_level_dim_0);
        Real* der_avg_global = averaged_derivative.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            der_avg_local[i]  = Real(0);
            der_avg_global[i] = Real(0);
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
            
            hier::IntVector ratio_to_coarsest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarsest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarsest_level/ratio_to_coarsest_level;
            
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
                 * Register the patch and the quantities in the flow model and compute the
                 * corresponding product.
                 */
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model->getNumberOfGhostCells();
                TBOX_ASSERT(num_ghosts >= num_ghosts_der);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], num_ghosts));
                }
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = getCellData(quantity_names[qi]);
                }
                
                std::vector<Real*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> num_ghosts_2_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                std::vector<int> ghostcell_dim_1_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                num_ghosts_2_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                ghostcell_dim_1_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    num_ghosts_2_u_qi.push_back(num_ghosts_u_qi[2]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                    ghostcell_dim_1_u_qi.push_back(ghostcell_dims_u_qi[1]);
                }
                
                /*
                 * Initialize cell data for the derivative and the product of variables.
                 * Also, get pointers to the cell data containers.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_derivative(
                    new pdat::CellData<Real>(patch_box, 1, hier::IntVector::getZero(d_dim)));
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_product(
                    new pdat::CellData<Real>(patch_box, 1, num_ghosts_der));
                
                Real* der     = data_derivative->getPointer(0);
                Real* product = data_product->getPointer(0);
                
                const hier::IntVector num_ghosts_product = data_product->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_product = data_product->getGhostBox().numberCells();
                
                const int num_ghosts_0_product = num_ghosts_product[0];
                const int num_ghosts_1_product = num_ghosts_product[1];
                const int num_ghosts_2_product = num_ghosts_product[2];
                const int ghostcell_dim_0_product = ghostcell_dims_product[0];
                const int ghostcell_dim_1_product = ghostcell_dims_product[1];
                
                const hier::IntVector patch_interior_dims = patch_box.numberCells();
                const int patch_interior_dim_0 = patch_interior_dims[0];
                const int patch_interior_dim_1 = patch_interior_dims[1];
                const int patch_interior_dim_2 = patch_interior_dims[2];
                
                data_product->fillAll(Real(1));
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    if (use_reciprocal[qi])
                    {
                        for (int k = -num_ghosts_2_product; k < patch_interior_dim_2 + num_ghosts_2_product; k++)
                        {
                            for (int j = -num_ghosts_1_product; j < patch_interior_dim_1 + num_ghosts_1_product; j++)
                            {
                                for (int i = -num_ghosts_0_product; i < patch_interior_dim_0 + num_ghosts_0_product; i++)
                                {
                                    const int idx_product = (i + num_ghosts_0_product) +
                                        (j + num_ghosts_1_product)*ghostcell_dim_0_product +
                                        (k + num_ghosts_2_product)*ghostcell_dim_0_product*
                                            ghostcell_dim_1_product;
                                    
                                    const int idx_u_qi = (i + num_ghosts_0_u_qi[qi]) +
                                        (j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                        (k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                            ghostcell_dim_1_u_qi[qi];
                                        
                                    product[idx_product] /= u_qi[qi][idx_u_qi];
                                }
                            }
                        }
                    }
                    else
                    {
                        for (int k = -num_ghosts_2_product; k < patch_interior_dim_2 + num_ghosts_2_product; k++)
                        {
                            for (int j = -num_ghosts_1_product; j < patch_interior_dim_1 + num_ghosts_1_product; j++)
                            {
                                for (int i = -num_ghosts_0_product; i < patch_interior_dim_0 + num_ghosts_0_product; i++)
                                {
                                    const int idx_product = (i + num_ghosts_0_product) +
                                        (j + num_ghosts_1_product)*ghostcell_dim_0_product +
                                        (k + num_ghosts_2_product)*ghostcell_dim_0_product*
                                            ghostcell_dim_1_product;
                                    
                                    const int idx_u_qi = (i + num_ghosts_0_u_qi[qi]) +
                                        (j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                        (k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                            ghostcell_dim_1_u_qi[qi];
                                        
                                    product[idx_product] *= u_qi[qi][idx_u_qi];
                                }
                            }
                        }
                    }
                }
                
                const Real weight = Real(dx[1]*dx[2])/(L_y*L_z);
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    if (derivative_direction == 0)
                    {
                        derivative_first_order_x->computeDerivative(
                            data_derivative,
                            data_product,
                            Real(dx[0]),
                            patch_visible_box,
                            0,
                            0);
                    }
                    else if (derivative_direction == 1)
                    {
                        derivative_first_order_y->computeDerivative(
                            data_derivative,
                            data_product,
                            Real(dx[1]),
                            patch_visible_box,
                            0,
                            0);
                    }
                    else if (derivative_direction == 2)
                    {
                         derivative_first_order_z->computeDerivative(
                            data_derivative,
                            data_product,
                            Real(dx[2]),
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
                                 * Compute the index of the data point and count how many times the data is repeated.
                                 */
                                
                                const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                                
                                int n_overlapped = 1;
                                
                                for (hier::BoxContainer::BoxContainerConstIterator iob(
                                        patch_overlapped_visible_boxes.begin());
                                     iob != patch_overlapped_visible_boxes.end();
                                     iob++)
                                {
                                    const hier::Box& patch_overlapped_visible_box = *iob;
                                    
                                    if (patch_overlapped_visible_box.contains(idx_pt))
                                    {
                                        n_overlapped++;
                                    }
                                }
                                
                                /*
                                 * Compute the linear index and the data to add.
                                 */
                                
                                const int idx = (relative_idx_lo_0 + i) +
                                    (relative_idx_lo_1 + j)*patch_interior_dim_0 +
                                    (relative_idx_lo_2 + k)*patch_interior_dim_0*
                                        patch_interior_dim_1;
                                
                                const Real value_to_add = der[idx]*weight/((Real) n_overlapped);
                                
                                /*
                                 * Add the data.
                                 */
                                
                                for (int ii = 0; ii < ratio_to_finest_level_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratio_to_finest_level_0 + ii;
                                    
                                    der_avg_local[idx_fine] += value_to_add;
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average of derivative.
         */
        
        d_mpi.Allreduce(
            der_avg_local,
            der_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(der_avg_local);
    }
    
    return averaged_derivative;
}


/*
 * Compute averaged value with only y-direction as inhomogeneous direction.
 */
std::vector<Real> FlowModelMPIHelperAverage::getAveragedQuantityWithInhomogeneousYDirection(
    const std::string quantity_name,
    const int component_idx,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    std::vector<Real> averaged_quantity;
    
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
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperAverage::getAveragedQuantityWithInhomogeneousYDirection():\n"
            << "Cannot compute averaged value with only y-direction as inhomogeneous direction for"
            << " one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_1 = d_finest_level_dims[1];
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_x = Real(x_hi[0] - x_lo[0]);
        
        Real* u_avg_local = (Real*)std::malloc(finest_level_dim_1*sizeof(Real));
        
        averaged_quantity.resize(finest_level_dim_1);
        Real* u_avg_global = averaged_quantity.data();
        
        for (int j = 0; j < finest_level_dim_1; j++)
        {
            u_avg_local[j]  = Real(0);
            u_avg_global[j] = Real(0);
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
            
            hier::IntVector ratio_to_coarsest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarsest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarsest_level/ratio_to_coarsest_level;
            
            const int ratio_to_finest_level_1 = ratio_to_finest_level[1];
            
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
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity =
                    getCellData(quantity_name);
                
                Real* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                
                const Real weight = Real(dx[0])/L_x;
                
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
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j);
                            
                            int n_overlapped = 1;
                            
                            for (hier::BoxContainer::BoxContainerConstIterator iob(
                                    patch_overlapped_visible_boxes.begin());
                                 iob != patch_overlapped_visible_boxes.end();
                                 iob++)
                            {
                                const hier::Box& patch_overlapped_visible_box = *iob;
                                
                                if (patch_overlapped_visible_box.contains(idx_pt))
                                {
                                    n_overlapped++;
                                }
                            }
                            
                            /*
                             * Compute the linear index and the data to add.
                             */
                            
                            const int idx = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                            
                            const Real value_to_add = u[idx]*weight/((Real) n_overlapped);
                            
                            for (int jj = 0; jj < ratio_to_finest_level_1; jj++)
                            {
                                const int idx_fine = (idx_lo_1 + j)*ratio_to_finest_level_1 + jj;
                                
                                u_avg_local[idx_fine] += value_to_add;
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        d_mpi.Allreduce(
            u_avg_local,
            u_avg_global,
            finest_level_dim_1,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_avg_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_1 = d_finest_level_dims[1];
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_x = Real(x_hi[0] - x_lo[0]);
        const Real L_z = Real(x_hi[2] - x_lo[2]);
        
        Real* u_avg_local = (Real*)std::malloc(finest_level_dim_1*sizeof(Real));
        
        averaged_quantity.resize(finest_level_dim_1);
        Real* u_avg_global = averaged_quantity.data();
        
        for (int j = 0; j < finest_level_dim_1; j++)
        {
            u_avg_local[j]  = Real(0);
            u_avg_global[j] = Real(0);
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
            
            hier::IntVector ratio_to_coarsest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarsest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarsest_level/ratio_to_coarsest_level;
            
            const int ratio_to_finest_level_1 = ratio_to_finest_level[1];
            
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
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity =
                    getCellData(quantity_name);
                
                Real* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int num_ghosts_2_quantity = num_ghosts_quantity[2];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                const int ghostcell_dim_1_quantity = ghostcell_dims_quantity[1];
                
                const Real weight = Real(dx[0]*dx[2])/(L_x*L_z);
                
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
                                 * Compute the index of the data point and count how many times the data is repeated.
                                 */
                                
                                const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                                
                                int n_overlapped = 1;
                                
                                for (hier::BoxContainer::BoxContainerConstIterator iob(
                                        patch_overlapped_visible_boxes.begin());
                                     iob != patch_overlapped_visible_boxes.end();
                                     iob++)
                                {
                                    const hier::Box& patch_overlapped_visible_box = *iob;
                                    
                                    if (patch_overlapped_visible_box.contains(idx_pt))
                                    {
                                        n_overlapped++;
                                    }
                                }
                                
                                /*
                                 * Compute the linear index and the data to add.
                                 */
                                
                                const int idx = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                    (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                        ghostcell_dim_1_quantity;
                                
                                const Real value_to_add = u[idx]*weight/((Real) n_overlapped);
                                
                                for (int jj = 0; jj < ratio_to_finest_level_1; jj++)
                                {
                                    const int idx_fine = (idx_lo_1 + j)*ratio_to_finest_level_1 + jj;
                                    
                                    u_avg_local[idx_fine] += value_to_add;
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        d_mpi.Allreduce(
            u_avg_local,
            u_avg_global,
            finest_level_dim_1,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_avg_local);
    }
    
    return averaged_quantity;
}


/*
 * Compute averaged reciprocal of value with only y-direction as inhomogeneous direction.
 */
std::vector<Real>
FlowModelMPIHelperAverage::getAveragedReciprocalOfQuantityWithInhomogeneousYDirection(
    const std::string quantity_name,
    const int component_idx,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    std::vector<Real> averaged_reciprocal_quantity;
    
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
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperAverage::getAveragedReciprocalOfQuantityWithInhomogeneousYDirection():\n"
            << "Cannot compute averaged value with only y-direction as inhomogeneous direction for"
            << " one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_1 = d_finest_level_dims[1];
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_x = Real(x_hi[0] - x_lo[0]);
        
        Real* u_inv_avg_local = (Real*)std::malloc(finest_level_dim_1*sizeof(Real));
        
        averaged_reciprocal_quantity.resize(finest_level_dim_1);
        Real* u_inv_avg_global = averaged_reciprocal_quantity.data();
        
        for (int j = 0; j < finest_level_dim_1; j++)
        {
            u_inv_avg_local[j]  = Real(0);
            u_inv_avg_global[j] = Real(0);
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
            
            hier::IntVector ratio_to_coarsest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarsest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarsest_level/ratio_to_coarsest_level;
            
            const int ratio_to_finest_level_1 = ratio_to_finest_level[1];
            
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
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity =
                    getCellData(quantity_name);
                
                Real* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                
                const Real weight = Real(dx[0])/L_x;
                
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
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j);
                            
                            int n_overlapped = 1;
                            
                            for (hier::BoxContainer::BoxContainerConstIterator iob(
                                    patch_overlapped_visible_boxes.begin());
                                 iob != patch_overlapped_visible_boxes.end();
                                 iob++)
                            {
                                const hier::Box& patch_overlapped_visible_box = *iob;
                                
                                if (patch_overlapped_visible_box.contains(idx_pt))
                                {
                                    n_overlapped++;
                                }
                            }
                            
                            /*
                             * Compute the linear index and the data to add.
                             */
                            
                            const int idx = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity;
                            
                            const Real value_to_add = (Real(1)/u[idx])*weight/((Real) n_overlapped);
                            
                            for (int jj = 0; jj < ratio_to_finest_level_1; jj++)
                            {
                                const int idx_fine = (idx_lo_1 + j)*ratio_to_finest_level_1 + jj;
                                
                                u_inv_avg_local[idx_fine] += value_to_add;
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average of reciprocal.
         */
        
        d_mpi.Allreduce(
            u_inv_avg_local,
            u_inv_avg_global,
            finest_level_dim_1,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_inv_avg_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_1 = d_finest_level_dims[1];
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_x = Real(x_hi[0] - x_lo[0]);
        const Real L_z = Real(x_hi[2] - x_lo[2]);
        
        Real* u_inv_avg_local = (Real*)std::malloc(finest_level_dim_1*sizeof(Real));
        
        averaged_reciprocal_quantity.resize(finest_level_dim_1);
        Real* u_inv_avg_global = averaged_reciprocal_quantity.data();
        
        for (int j = 0; j < finest_level_dim_1; j++)
        {
            u_inv_avg_local[j]  = Real(0);
            u_inv_avg_global[j] = Real(0);
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
            
            hier::IntVector ratio_to_coarsest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarsest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarsest_level/ratio_to_coarsest_level;
            
            const int ratio_to_finest_level_1 = ratio_to_finest_level[1];
            
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
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity =
                    getCellData(quantity_name);
                
                Real* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int num_ghosts_2_quantity = num_ghosts_quantity[2];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                const int ghostcell_dim_1_quantity = ghostcell_dims_quantity[1];
                
                const Real weight = Real(dx[0]*dx[2])/(L_x*L_z);
                
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
                                 * Compute the index of the data point and count how many times the data is repeated.
                                 */
                                
                                const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                                
                                int n_overlapped = 1;
                                
                                for (hier::BoxContainer::BoxContainerConstIterator iob(
                                        patch_overlapped_visible_boxes.begin());
                                     iob != patch_overlapped_visible_boxes.end();
                                     iob++)
                                {
                                    const hier::Box& patch_overlapped_visible_box = *iob;
                                    
                                    if (patch_overlapped_visible_box.contains(idx_pt))
                                    {
                                        n_overlapped++;
                                    }
                                }
                                
                                /*
                                 * Compute the linear index and the data to add.
                                 */
                                
                                const int idx = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                    (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                        ghostcell_dim_1_quantity;
                                
                                const Real value_to_add = (Real(1)/u[idx])*weight/((Real) n_overlapped);
                                
                                for (int jj = 0; jj < ratio_to_finest_level_1; jj++)
                                {
                                    const int idx_fine = (idx_lo_1 + j)*ratio_to_finest_level_1 + jj;
                                    
                                    u_inv_avg_local[idx_fine] += value_to_add;
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average of reciprocal.
         */
        
        d_mpi.Allreduce(
            u_inv_avg_local,
            u_inv_avg_global,
            finest_level_dim_1,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_inv_avg_local);
    }
    
    return averaged_reciprocal_quantity;
}


/*
 * Compute averaged value (on product of variables) with only y-direction as inhomogeneous direction.
 */
std::vector<Real>
FlowModelMPIHelperAverage::getAveragedQuantityWithInhomogeneousYDirection(
    const std::vector<std::string>& quantity_names,
    const std::vector<int>& component_indices,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    int num_quantities = static_cast<int>(quantity_names.size());
    
    TBOX_ASSERT(static_cast<int>(component_indices.size()) == num_quantities);
    
    std::vector<bool> use_reciprocal(num_quantities, false);
    
    return getAveragedQuantityWithInhomogeneousYDirection(
        quantity_names,
        component_indices,
        use_reciprocal,
        data_context);
}


/*
 * Compute averaged value (on product of variables) with only y-direction as inhomogeneous direction.
 */
std::vector<Real>
FlowModelMPIHelperAverage::getAveragedQuantityWithInhomogeneousYDirection(
    const std::vector<std::string>& quantity_names,
    const std::vector<int>& component_indices,
    const std::vector<bool>& use_reciprocal,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    int num_quantities = static_cast<int>(quantity_names.size());
    
    TBOX_ASSERT(static_cast<int>(component_indices.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(use_reciprocal.size()) == num_quantities);
    
    std::vector<Real> averaged_quantity;
    
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
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperAverage::getAveragedQuantityWithInhomogeneousYDirection():\n"
            << "Cannot compute averaged value with only y-direction as inhomogeneous direction for"
            << " one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_1 = d_finest_level_dims[1];
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_x = Real(x_hi[0] - x_lo[0]);
        
        Real* avg_local = (Real*)std::malloc(finest_level_dim_1*sizeof(Real));
        
        averaged_quantity.resize(finest_level_dim_1);
        Real* avg_global = averaged_quantity.data();
        
        for (int j = 0; j < finest_level_dim_1; j++)
        {
            avg_local[j]  = Real(0);
            avg_global[j] = Real(0);
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
            
            hier::IntVector ratio_to_coarsest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarsest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarsest_level/ratio_to_coarsest_level;
            
            const int ratio_to_finest_level_1 = ratio_to_finest_level[1];
            
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
                 * Register the patch and data in the flow model and compute the corresponding
                 * average.
                 */
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                }
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = getCellData(quantity_names[qi]);
                }
                
                std::vector<Real*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                }
                
                const Real weight = Real(dx[0])/L_x;
                
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
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j);
                            
                            int n_overlapped = 1;
                            
                            for (hier::BoxContainer::BoxContainerConstIterator iob(
                                    patch_overlapped_visible_boxes.begin());
                                 iob != patch_overlapped_visible_boxes.end();
                                 iob++)
                            {
                                const hier::Box& patch_overlapped_visible_box = *iob;
                                
                                if (patch_overlapped_visible_box.contains(idx_pt))
                                {
                                    n_overlapped++;
                                }
                            }
                            
                            /*
                             * Compute the linear indices and the data to add.
                             */
                            
                            for (int jj = 0; jj < ratio_to_finest_level_1; jj++)
                            {
                                const int idx_fine = (idx_lo_1 + j)*ratio_to_finest_level_1 + jj;
                                
                                Real avg = Real(1);
                                
                                for (int qi = 0; qi < num_quantities; qi++)
                                {
                                    const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    if (use_reciprocal[qi])
                                    {
                                        avg /= u_qi[qi][idx_qi];
                                    }
                                    else
                                    {
                                        avg *= u_qi[qi][idx_qi];
                                    }
                                }
                                
                                avg_local[idx_fine] += (avg*weight/((Real) n_overlapped));
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        d_mpi.Allreduce(
            avg_local,
            avg_global,
            finest_level_dim_1,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(avg_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_1 = d_finest_level_dims[1];
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_x = Real(x_hi[0] - x_lo[0]);
        const Real L_z = Real(x_hi[2] - x_lo[2]);
        
        Real* avg_local = (Real*)std::malloc(finest_level_dim_1*sizeof(Real));
        
        averaged_quantity.resize(finest_level_dim_1);
        Real* avg_global = averaged_quantity.data();
        
        for (int j = 0; j < finest_level_dim_1; j++)
        {
            avg_local[j]  = Real(0);
            avg_global[j] = Real(0);
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
            
            hier::IntVector ratio_to_coarsest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarsest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarsest_level/ratio_to_coarsest_level;
            
            const int ratio_to_finest_level_1 = ratio_to_finest_level[1];
            
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
                 * Register the patch and data in the flow model and compute the corresponding
                 * average.
                 */
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                }
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = getCellData(quantity_names[qi]);
                }
                
                std::vector<Real*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> num_ghosts_2_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                std::vector<int> ghostcell_dim_1_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                num_ghosts_2_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                ghostcell_dim_1_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    num_ghosts_2_u_qi.push_back(num_ghosts_u_qi[2]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                    ghostcell_dim_1_u_qi.push_back(ghostcell_dims_u_qi[1]);
                }
                
                const Real weight = Real(dx[0]*dx[2])/(L_x*L_z);
                
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
                                 * Compute the index of the data point and count how many times the data is repeated.
                                 */
                                
                                const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                                
                                int n_overlapped = 1;
                                
                                for (hier::BoxContainer::BoxContainerConstIterator iob(
                                        patch_overlapped_visible_boxes.begin());
                                     iob != patch_overlapped_visible_boxes.end();
                                     iob++)
                                {
                                    const hier::Box& patch_overlapped_visible_box = *iob;
                                    
                                    if (patch_overlapped_visible_box.contains(idx_pt))
                                    {
                                        n_overlapped++;
                                    }
                                }
                                
                                /*
                                 * Compute the linear index and the data to add.
                                 */
                                
                                for (int jj = 0; jj < ratio_to_finest_level_1; jj++)
                                {
                                    const int idx_fine = (idx_lo_1 + j)*ratio_to_finest_level_1 + jj;
                                    
                                    Real avg = Real(1);
                                    
                                    for (int qi = 0; qi < num_quantities; qi++)
                                    {
                                        const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        if (use_reciprocal[qi])
                                        {
                                            avg /= u_qi[qi][idx_qi];
                                        }
                                        else
                                        {
                                            avg *= u_qi[qi][idx_qi];
                                        }
                                    }
                                    
                                    avg_local[idx_fine] += (avg*weight/((Real) n_overlapped));
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        d_mpi.Allreduce(
            avg_local,
            avg_global,
            finest_level_dim_1,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(avg_local);
    }
    
    return averaged_quantity;
}


/*
 * Compute averaged value (on product of variable derivatives) with only y direction as inhomogeneous direction.
 */
std::vector<Real>
FlowModelMPIHelperAverage::getAveragedQuantityWithInhomogeneousYDirection(
    const std::vector<std::string>& quantity_names,
    const std::vector<int>& component_indices,
    const std::vector<bool>& use_derivative,
    const std::vector<int>& derivative_directions,
    const int num_ghosts_derivative,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    int num_quantities = static_cast<int>(quantity_names.size());
    
    TBOX_ASSERT(static_cast<int>(component_indices.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(use_derivative.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(derivative_directions.size()) == num_quantities);
    
    std::vector<bool> use_reciprocal(num_quantities, false);
    
    return getAveragedQuantityWithInhomogeneousYDirection(
        quantity_names,
        component_indices,
        use_derivative,
        derivative_directions,
        use_reciprocal,
        num_ghosts_derivative,
        data_context);
}


/*
 * Compute averaged value (on product of variable derivatives) with only y direction as inhomogeneous direction.
 */
std::vector<Real>
FlowModelMPIHelperAverage::getAveragedQuantityWithInhomogeneousYDirection(
    const std::vector<std::string>& quantity_names,
    const std::vector<int>& component_indices,
    const std::vector<bool>& use_derivative,
    const std::vector<int>& derivative_directions,
    const std::vector<bool>& use_reciprocal,
    const int num_ghosts_derivative,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    int num_quantities = static_cast<int>(quantity_names.size());
    
    TBOX_ASSERT(static_cast<int>(component_indices.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(use_derivative.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(derivative_directions.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(use_reciprocal.size()) == num_quantities);
    
    int num_use_derivative = 0;
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperAverage::getAveragedQuantityWithInhomogeneousYDirection():\n"
            << "Cannot compute averaged value with only y-direction as inhomogeneous direction for"
            << " one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        for (int qi = 0; qi < num_quantities; qi++)
        {
            if (use_derivative[qi])
            {
                num_use_derivative++;
                if (derivative_directions[qi] < 0 || derivative_directions[qi] > 1)
                {
                    TBOX_ERROR(d_object_name
                        << ": FlowModelMPIHelperAverage::getAveragedQuantityWithInhomogeneousYDirection():\n"
                        << "Cannot take derivative for two-dimensional problem!\n"
                        << "derivative_directions[" << qi << "] = " << derivative_directions[qi] << " given!\n"
                        << std::endl);
                }
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        for (int qi = 0; qi < num_quantities; qi++)
        {
            if (use_derivative[qi])
            {
                num_use_derivative++;
                if (derivative_directions[qi] < 0 || derivative_directions[qi] > 2)
                {
                    TBOX_ERROR(d_object_name
                        << ": FlowModelMPIHelperAverage::getAveragedQuantityWithInhomogeneousYDirection():\n"
                        << "Cannot take derivative for three-dimensional problem!\n"
                        << "derivative_directions[" << qi << "] = " << derivative_directions[qi] << " given!\n"
                        << std::endl);
                }
            }
        }
    }
    
    std::vector<Real> averaged_quantity;
    
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
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(2))
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
        
        hier::IntVector num_ghosts_der = hier::IntVector::getOne(d_dim)*num_ghosts_derivative;
        
        const int finest_level_dim_1 = d_finest_level_dims[1];
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_x = Real(x_hi[0] - x_lo[0]);
        
        Real* avg_local = (Real*)std::malloc(finest_level_dim_1*sizeof(Real));
        
        averaged_quantity.resize(finest_level_dim_1);
        Real* avg_global = averaged_quantity.data();
        
        for (int j = 0; j < finest_level_dim_1; j++)
        {
            avg_local[j]  = Real(0);
            avg_global[j] = Real(0);
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
            
            hier::IntVector ratio_to_coarsest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarsest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarsest_level/ratio_to_coarsest_level;
            
            const int ratio_to_finest_level_1 = ratio_to_finest_level[1];
            
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
                 * Register the patch and data in the flow model and compute the corresponding
                 * average.
                 */
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    if (num_use_derivative > 0)
                    {
                        num_subghosts_of_data.insert(
                            std::pair<std::string, hier::IntVector>(quantity_names[qi], num_ghosts_der));
                    }
                    else
                    {
                        num_subghosts_of_data.insert(
                            std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                    }
                }
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = getCellData(quantity_names[qi]);
                }
                
                std::vector<Real*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                }
                
                /*
                 * Initialize cell data for the derivatives and get pointers to the cell data.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_derivative;
                std::vector<Real*> der_qi;
                
                if (num_use_derivative > 0)
                {
                    data_derivative = HAMERS_MAKE_SHARED<pdat::CellData<Real> >(patch_box, num_use_derivative, hier::IntVector::getZero(d_dim));
                    
                    der_qi.resize(num_use_derivative);
                    for (int qi = 0; qi < num_use_derivative; qi++)
                    {
                        der_qi[qi] = data_derivative->getPointer(qi);
                    }
                }
                
                const hier::IntVector patch_interior_dims = patch_box.numberCells();
                const int patch_interior_dim_0 = patch_interior_dims[0];
                
                const Real weight = Real(dx[0])/L_x;
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    int count_derivative = 0;
                    for (int qi = 0; qi < num_quantities; qi++)
                    {
                        if (use_derivative[qi] && derivative_directions[qi] == 0)
                        {
                            derivative_first_order_x->computeDerivative(
                                data_derivative,
                                data_quantities[qi],
                                Real(dx[0]),
                                patch_visible_box,
                                count_derivative,
                                component_indices[qi]);
                            
                            count_derivative++;
                        }
                        else if (use_derivative[qi] && derivative_directions[qi] == 1)
                        {
                            derivative_first_order_y->computeDerivative(
                                data_derivative,
                                data_quantities[qi],
                                Real(dx[1]),
                                patch_visible_box,
                                count_derivative,
                                component_indices[qi]);
                            
                            count_derivative++;
                        }
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
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j);
                            
                            int n_overlapped = 1;
                            
                            for (hier::BoxContainer::BoxContainerConstIterator iob(
                                    patch_overlapped_visible_boxes.begin());
                                 iob != patch_overlapped_visible_boxes.end();
                                 iob++)
                            {
                                const hier::Box& patch_overlapped_visible_box = *iob;
                                
                                if (patch_overlapped_visible_box.contains(idx_pt))
                                {
                                    n_overlapped++;
                                }
                            }
                            
                            /*
                             * Compute the linear indices and the data to add.
                             */
                            
                            Real avg = Real(1);
                            
                            count_derivative = 0;
                            for (int qi = 0; qi < num_quantities; qi++)
                            {
                                if (use_reciprocal[qi])
                                {
                                    if (use_derivative[qi])
                                    {
                                        const int idx_der = (relative_idx_lo_0 + i) +
                                            (relative_idx_lo_1 + j)*patch_interior_dim_0;
                                        
                                        avg /= der_qi[count_derivative][idx_der];
                                        
                                        count_derivative++;
                                    }
                                    else
                                    {
                                        const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                        
                                        avg /= u_qi[qi][idx_qi];
                                    }
                                }
                                else
                                {
                                    if (use_derivative[qi])
                                    {
                                        const int idx_der = (relative_idx_lo_0 + i) +
                                            (relative_idx_lo_1 + j)*patch_interior_dim_0;
                                        
                                        avg *= der_qi[count_derivative][idx_der];
                                        
                                        count_derivative++;
                                    }
                                    else
                                    {
                                        const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                        
                                        avg *= u_qi[qi][idx_qi];
                                    }
                                }
                            }
                            
                            /*
                             * Add the data.
                             */
                            
                            for (int jj = 0; jj < ratio_to_finest_level_1; jj++)
                            {
                                const int idx_fine = (idx_lo_1 + j)*ratio_to_finest_level_1 + jj;
                                
                                avg_local[idx_fine] += (avg*weight/((Real) n_overlapped));
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        d_mpi.Allreduce(
            avg_local,
            avg_global,
            finest_level_dim_1,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(avg_local);
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
        
        hier::IntVector num_ghosts_der = hier::IntVector::getOne(d_dim)*num_ghosts_derivative;
        
        const int finest_level_dim_1 = d_finest_level_dims[1];
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_x = Real(x_hi[0] - x_lo[0]);
        const Real L_z = Real(x_hi[2] - x_lo[2]);
        
        Real* avg_local = (Real*)std::malloc(finest_level_dim_1*sizeof(Real));
        
        averaged_quantity.resize(finest_level_dim_1);
        Real* avg_global = averaged_quantity.data();
        
        for (int j = 0; j < finest_level_dim_1; j++)
        {
            avg_local[j]  = Real(0);
            avg_global[j] = Real(0);
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
            
            hier::IntVector ratio_to_coarsest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarsest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarsest_level/ratio_to_coarsest_level;
            
            const int ratio_to_finest_level_1 = ratio_to_finest_level[1];
            
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
                 * Register the patch and data in the flow model and compute the corresponding
                 * average.
                 */
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    if (num_use_derivative > 0)
                    {
                        num_subghosts_of_data.insert(
                            std::pair<std::string, hier::IntVector>(quantity_names[qi], num_ghosts_der));
                    }
                    else
                    {
                        num_subghosts_of_data.insert(
                            std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                    }
                }
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = getCellData(quantity_names[qi]);
                }
                
                std::vector<Real*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> num_ghosts_2_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                std::vector<int> ghostcell_dim_1_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                num_ghosts_2_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                ghostcell_dim_1_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    num_ghosts_2_u_qi.push_back(num_ghosts_u_qi[2]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                    ghostcell_dim_1_u_qi.push_back(ghostcell_dims_u_qi[1]);
                }
                
                /*
                 * Initialize cell data for the derivatives and get pointers to the cell data.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_derivative;
                std::vector<Real*> der_qi;
                
                if (num_use_derivative > 0)
                {
                    data_derivative = HAMERS_MAKE_SHARED<pdat::CellData<Real> >(patch_box, num_use_derivative, hier::IntVector::getZero(d_dim));
                    
                    der_qi.resize(num_use_derivative);
                    for (int qi = 0; qi < num_use_derivative; qi++)
                    {
                        der_qi[qi] = data_derivative->getPointer(qi);
                    }
                }
                
                const hier::IntVector patch_interior_dims = patch_box.numberCells();
                const int patch_interior_dim_0 = patch_interior_dims[0];
                const int patch_interior_dim_1 = patch_interior_dims[1];
                
                const Real weight = Real(dx[0]*dx[2])/(L_x*L_z);
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    int count_derivative = 0;
                    for (int qi = 0; qi < num_quantities; qi++)
                    {
                        if (use_derivative[qi] && derivative_directions[qi] == 0)
                        {
                            derivative_first_order_x->computeDerivative(
                                data_derivative,
                                data_quantities[qi],
                                Real(dx[0]),
                                patch_visible_box,
                                count_derivative,
                                component_indices[qi]);
                            
                            count_derivative++;
                        }
                        else if (use_derivative[qi] && derivative_directions[qi] == 1)
                        {
                            derivative_first_order_y->computeDerivative(
                                data_derivative,
                                data_quantities[qi],
                                Real(dx[1]),
                                patch_visible_box,
                                count_derivative,
                                component_indices[qi]);
                            
                            count_derivative++;
                        }
                        else if (use_derivative[qi] && derivative_directions[qi] == 2)
                        {
                            derivative_first_order_z->computeDerivative(
                                data_derivative,
                                data_quantities[qi],
                                Real(dx[2]),
                                patch_visible_box,
                                count_derivative,
                                component_indices[qi]);
                            
                            count_derivative++;
                        }
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
                                 * Compute the index of the data point and count how many times the data is repeated.
                                 */
                                
                                const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                                
                                int n_overlapped = 1;
                                
                                for (hier::BoxContainer::BoxContainerConstIterator iob(
                                        patch_overlapped_visible_boxes.begin());
                                     iob != patch_overlapped_visible_boxes.end();
                                     iob++)
                                {
                                    const hier::Box& patch_overlapped_visible_box = *iob;
                                    
                                    if (patch_overlapped_visible_box.contains(idx_pt))
                                    {
                                        n_overlapped++;
                                    }
                                }
                                
                                /*
                                 * Compute the linear index and the data to add.
                                 */
                                
                                Real avg = Real(1);
                                
                                count_derivative = 0;
                                for (int qi = 0; qi < num_quantities; qi++)
                                {
                                    if (use_reciprocal[qi])
                                    {
                                        if (use_derivative[qi])
                                        {
                                            const int idx_der = (relative_idx_lo_0 + i) +
                                                (relative_idx_lo_1 + j)*patch_interior_dim_0 +
                                                (relative_idx_lo_2 + k)*patch_interior_dim_0*
                                                    patch_interior_dim_1;
                                            
                                            avg /= der_qi[count_derivative][idx_der];
                                            
                                            count_derivative++;
                                        }
                                        else
                                        {
                                            const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                                (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                                (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                    ghostcell_dim_1_u_qi[qi];
                                            
                                            avg /= u_qi[qi][idx_qi];
                                        }
                                    }
                                    else
                                    {
                                        if (use_derivative[qi])
                                        {
                                            const int idx_der = (relative_idx_lo_0 + i) +
                                                (relative_idx_lo_1 + j)*patch_interior_dim_0 +
                                                (relative_idx_lo_2 + k)*patch_interior_dim_0*
                                                    patch_interior_dim_1;
                                            
                                            avg *= der_qi[count_derivative][idx_der];
                                            
                                            count_derivative++;
                                        }
                                        else
                                        {
                                            const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                                (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                                (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                    ghostcell_dim_1_u_qi[qi];
                                            
                                            avg *= u_qi[qi][idx_qi];
                                        }
                                    }
                                }
                                
                                /*
                                 * Add the data.
                                 */
                                
                                for (int jj = 0; jj < ratio_to_finest_level_1; jj++)
                                {
                                    const int idx_fine = (idx_lo_1 + j)*ratio_to_finest_level_1 + jj;
                                    
                                    avg_local[idx_fine] += (avg*weight/((Real) n_overlapped));
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        d_mpi.Allreduce(
            avg_local,
            avg_global,
            finest_level_dim_1,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(avg_local);
    }
    
    return averaged_quantity;
}


/*
 * Compute averaged derivative of value (on product of variables) with only y direction as inhomogeneous direction.
 */
std::vector<Real>
FlowModelMPIHelperAverage::getAveragedDerivativeOfQuantityWithInhomogeneousYDirection(
    const std::vector<std::string>& quantity_names,
    const std::vector<int>& component_indices,
    const std::vector<bool>& use_reciprocal,
    const int derivative_direction,
    const int num_ghosts_derivative,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    int num_quantities = static_cast<int>(quantity_names.size());
    
    TBOX_ASSERT(static_cast<int>(component_indices.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(use_reciprocal.size()) == num_quantities);
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperAverage::getAveragedDerivativeOfQuantityWithInhomogeneousYDirection():\n"
            << "Cannot compute averaged value with only y-direction as inhomogeneous direction for"
            << " one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2) && (derivative_direction < 0 || derivative_direction > 1))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperAverage::getAveragedDerivativeOfQuantityWithInhomogeneousYDirection():\n"
            << "Cannot take derivative for two-dimensional problem!\n"
            << "derivative_direction = " << derivative_direction << " given!\n"
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(3) && (derivative_direction < 0 || derivative_direction > 2))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperAverage::getAveragedDerivativeOfQuantityWithInhomogeneousYDirection():\n"
            << "Cannot take derivative for three-dimensional problem!\n"
            << "derivative_direction = " << derivative_direction << " given!\n"
            << std::endl);
    }
    
    hier::IntVector num_ghosts_der = hier::IntVector::getZero(d_dim);
    num_ghosts_der[derivative_direction] = num_ghosts_derivative;
    
    std::vector<Real> averaged_derivative;
    
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
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(2))
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
        
        const int finest_level_dim_1 = d_finest_level_dims[1];
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_x = Real(x_hi[0] - x_lo[0]);
        
        Real* der_avg_local = (Real*)std::malloc(finest_level_dim_1*sizeof(Real));
        
        averaged_derivative.resize(finest_level_dim_1);
        Real* der_avg_global = averaged_derivative.data();
        
        for (int j = 0; j < finest_level_dim_1; j++)
        {
            der_avg_local[j]  = Real(0);
            der_avg_global[j] = Real(0);
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
            
            hier::IntVector ratio_to_coarsest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarsest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarsest_level/ratio_to_coarsest_level;
            
            const int ratio_to_finest_level_1 = ratio_to_finest_level[1];
            
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
                 * Register the patch and the quantities in the flow model and compute the
                 * corresponding product.
                 */
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model->getNumberOfGhostCells();
                TBOX_ASSERT(num_ghosts >= num_ghosts_der);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], num_ghosts_der));
                }
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = getCellData(quantity_names[qi]);
                }
                
                std::vector<Real*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                }
                
                /*
                 * Initialize cell data for the derivative and the product of variables.
                 * Also, get pointers to the cell data containers.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_derivative(
                    new pdat::CellData<Real>(patch_box, 1, hier::IntVector::getZero(d_dim)));
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_product(
                    new pdat::CellData<Real>(patch_box, 1, num_ghosts_der));
                
                Real* der     = data_derivative->getPointer(0);
                Real* product = data_product->getPointer(0);
                
                const hier::IntVector num_ghosts_product = data_product->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_product = data_product->getGhostBox().numberCells();
                
                const int num_ghosts_0_product = num_ghosts_product[0];
                const int num_ghosts_1_product = num_ghosts_product[1];
                const int ghostcell_dim_0_product = ghostcell_dims_product[0];
                
                const hier::IntVector patch_interior_dims = patch_box.numberCells();
                const int patch_interior_dim_0 = patch_interior_dims[0];
                const int patch_interior_dim_1 = patch_interior_dims[1];
                
                data_product->fillAll(Real(1));
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    if (use_reciprocal[qi])
                    {
                        for (int j = -num_ghosts_1_product; j < patch_interior_dim_1 + num_ghosts_1_product; j++)
                        {
                            for (int i = -num_ghosts_0_product; i < patch_interior_dim_0 + num_ghosts_0_product; i++)
                            {
                                const int idx_product = (i + num_ghosts_0_product) +
                                    (j + num_ghosts_1_product)*ghostcell_dim_0_product;
                                
                                const int idx_u_qi = (i + num_ghosts_0_u_qi[qi]) +
                                    (j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                product[idx_product] /= u_qi[qi][idx_u_qi];
                            }
                        }
                    }
                    else
                    {
                        for (int j = -num_ghosts_1_product; j < patch_interior_dim_1 + num_ghosts_1_product; j++)
                        {
                            for (int i = -num_ghosts_0_product; i < patch_interior_dim_0 + num_ghosts_0_product; i++)
                            {
                                const int idx_product = (i + num_ghosts_0_product) +
                                    (j + num_ghosts_1_product)*ghostcell_dim_0_product;
                                
                                const int idx_u_qi = (i + num_ghosts_0_u_qi[qi]) +
                                    (j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                product[idx_product] *= u_qi[qi][idx_u_qi];
                            }
                        }
                    }
                }
                
                const Real weight = Real(dx[0])/L_x;
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    if (derivative_direction == 0)
                    {
                        derivative_first_order_x->computeDerivative(
                            data_derivative,
                            data_product,
                            Real(dx[0]),
                            patch_visible_box,
                            0,
                            0);
                    }
                    else if (derivative_direction == 1)
                    {
                        derivative_first_order_y->computeDerivative(
                            data_derivative,
                            data_product,
                            Real(dx[1]),
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
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j);
                            
                            int n_overlapped = 1;
                            
                            for (hier::BoxContainer::BoxContainerConstIterator iob(
                                    patch_overlapped_visible_boxes.begin());
                                 iob != patch_overlapped_visible_boxes.end();
                                 iob++)
                            {
                                const hier::Box& patch_overlapped_visible_box = *iob;
                                
                                if (patch_overlapped_visible_box.contains(idx_pt))
                                {
                                    n_overlapped++;
                                }
                            }
                            
                            /*
                             * Compute the linear index and the data to add.
                             */
                            
                            const int idx = (relative_idx_lo_0 + i) +
                                (relative_idx_lo_1 + j)*patch_interior_dim_0;
                            
                            const Real value_to_add = der[idx]*weight/((Real) n_overlapped);
                            
                            /*
                             * Add the data.
                             */
                            
                            for (int jj = 0; jj < ratio_to_finest_level_1; jj++)
                            {
                                const int idx_fine = (idx_lo_1 + j)*ratio_to_finest_level_1 + jj;
                                
                                der_avg_local[idx_fine] += value_to_add;
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average of derivative.
         */
        
        d_mpi.Allreduce(
            der_avg_local,
            der_avg_global,
            finest_level_dim_1,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(der_avg_local);
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
        
        const int finest_level_dim_1 = d_finest_level_dims[1];
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_x = Real(x_hi[0] - x_lo[0]);
        const Real L_z = Real(x_hi[2] - x_lo[2]);
        
        Real* der_avg_local = (Real*)std::malloc(finest_level_dim_1*sizeof(Real));
        
        averaged_derivative.resize(finest_level_dim_1);
        Real* der_avg_global = averaged_derivative.data();
        
        for (int j = 0; j < finest_level_dim_1; j++)
        {
            der_avg_local[j]  = Real(0);
            der_avg_global[j] = Real(0);
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
            
            hier::IntVector ratio_to_coarsest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarsest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarsest_level/ratio_to_coarsest_level;
            
            const int ratio_to_finest_level_1 = ratio_to_finest_level[1];
            
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
                 * Register the patch and the quantities in the flow model and compute the
                 * corresponding product.
                 */
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model->getNumberOfGhostCells();
                TBOX_ASSERT(num_ghosts >= num_ghosts_der);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], num_ghosts));
                }
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = getCellData(quantity_names[qi]);
                }
                
                std::vector<Real*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> num_ghosts_2_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                std::vector<int> ghostcell_dim_1_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                num_ghosts_2_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                ghostcell_dim_1_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    num_ghosts_2_u_qi.push_back(num_ghosts_u_qi[2]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                    ghostcell_dim_1_u_qi.push_back(ghostcell_dims_u_qi[1]);
                }
                
                /*
                 * Initialize cell data for the derivative and the product of variables.
                 * Also, get pointers to the cell data containers.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_derivative(
                    new pdat::CellData<Real>(patch_box, 1, hier::IntVector::getZero(d_dim)));
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_product(
                    new pdat::CellData<Real>(patch_box, 1, num_ghosts_der));
                
                Real* der     = data_derivative->getPointer(0);
                Real* product = data_product->getPointer(0);
                
                const hier::IntVector num_ghosts_product = data_product->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_product = data_product->getGhostBox().numberCells();
                
                const int num_ghosts_0_product = num_ghosts_product[0];
                const int num_ghosts_1_product = num_ghosts_product[1];
                const int num_ghosts_2_product = num_ghosts_product[2];
                const int ghostcell_dim_0_product = ghostcell_dims_product[0];
                const int ghostcell_dim_1_product = ghostcell_dims_product[1];
                
                const hier::IntVector patch_interior_dims = patch_box.numberCells();
                const int patch_interior_dim_0 = patch_interior_dims[0];
                const int patch_interior_dim_1 = patch_interior_dims[1];
                const int patch_interior_dim_2 = patch_interior_dims[2];
                
                data_product->fillAll(Real(1));
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    if (use_reciprocal[qi])
                    {
                        for (int k = -num_ghosts_2_product; k < patch_interior_dim_2 + num_ghosts_2_product; k++)
                        {
                            for (int j = -num_ghosts_1_product; j < patch_interior_dim_1 + num_ghosts_1_product; j++)
                            {
                                for (int i = -num_ghosts_0_product; i < patch_interior_dim_0 + num_ghosts_0_product; i++)
                                {
                                    const int idx_product = (i + num_ghosts_0_product) +
                                        (j + num_ghosts_1_product)*ghostcell_dim_0_product +
                                        (k + num_ghosts_2_product)*ghostcell_dim_0_product*
                                            ghostcell_dim_1_product;
                                    
                                    const int idx_u_qi = (i + num_ghosts_0_u_qi[qi]) +
                                        (j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                        (k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                            ghostcell_dim_1_u_qi[qi];
                                        
                                    product[idx_product] /= u_qi[qi][idx_u_qi];
                                }
                            }
                        }
                    }
                    else
                    {
                        for (int k = -num_ghosts_2_product; k < patch_interior_dim_2 + num_ghosts_2_product; k++)
                        {
                            for (int j = -num_ghosts_1_product; j < patch_interior_dim_1 + num_ghosts_1_product; j++)
                            {
                                for (int i = -num_ghosts_0_product; i < patch_interior_dim_0 + num_ghosts_0_product; i++)
                                {
                                    const int idx_product = (i + num_ghosts_0_product) +
                                        (j + num_ghosts_1_product)*ghostcell_dim_0_product +
                                        (k + num_ghosts_2_product)*ghostcell_dim_0_product*
                                            ghostcell_dim_1_product;
                                    
                                    const int idx_u_qi = (i + num_ghosts_0_u_qi[qi]) +
                                        (j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                        (k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                            ghostcell_dim_1_u_qi[qi];
                                        
                                    product[idx_product] *= u_qi[qi][idx_u_qi];
                                }
                            }
                        }
                    }
                }
                
                const Real weight = Real(dx[0]*dx[2])/(L_x*L_z);
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    if (derivative_direction == 0)
                    {
                        derivative_first_order_x->computeDerivative(
                            data_derivative,
                            data_product,
                            Real(dx[0]),
                            patch_visible_box,
                            0,
                            0);
                    }
                    else if (derivative_direction == 1)
                    {
                        derivative_first_order_y->computeDerivative(
                            data_derivative,
                            data_product,
                            Real(dx[1]),
                            patch_visible_box,
                            0,
                            0);
                    }
                    else if (derivative_direction == 2)
                    {
                         derivative_first_order_z->computeDerivative(
                            data_derivative,
                            data_product,
                            Real(dx[2]),
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
                                 * Compute the index of the data point and count how many times the data is repeated.
                                 */
                                
                                const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                                
                                int n_overlapped = 1;
                                
                                for (hier::BoxContainer::BoxContainerConstIterator iob(
                                        patch_overlapped_visible_boxes.begin());
                                     iob != patch_overlapped_visible_boxes.end();
                                     iob++)
                                {
                                    const hier::Box& patch_overlapped_visible_box = *iob;
                                    
                                    if (patch_overlapped_visible_box.contains(idx_pt))
                                    {
                                        n_overlapped++;
                                    }
                                }
                                
                                /*
                                 * Compute the linear index and the data to add.
                                 */
                                
                                const int idx = (relative_idx_lo_0 + i) +
                                    (relative_idx_lo_1 + j)*patch_interior_dim_0 +
                                    (relative_idx_lo_2 + k)*patch_interior_dim_0*
                                        patch_interior_dim_1;
                                
                                const Real value_to_add = der[idx]*weight/((Real) n_overlapped);
                                
                                /*
                                 * Add the data.
                                 */
                                
                                for (int jj = 0; jj < ratio_to_finest_level_1; jj++)
                                {
                                    const int idx_fine = (idx_lo_1 + j)*ratio_to_finest_level_1 + jj;
                                    
                                    der_avg_local[idx_fine] += value_to_add;
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average of derivative.
         */
        
        d_mpi.Allreduce(
            der_avg_local,
            der_avg_global,
            finest_level_dim_1,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(der_avg_local);
    }
    
    return averaged_derivative;
}


/*
 * Compute averaged value with only z-direction as inhomogeneous direction.
 */
std::vector<Real> FlowModelMPIHelperAverage::getAveragedQuantityWithInhomogeneousZDirection(
    const std::string quantity_name,
    const int component_idx,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    std::vector<Real> averaged_quantity;
    
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
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperAverage::getAveragedQuantityWithInhomogeneousZDirection():\n"
            << "Cannot compute averaged value with only z-direction as inhomogeneous direction for"
            << " one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperAverage::getAveragedQuantityWithInhomogeneousZDirection():\n"
            << "Cannot compute averaged value with only z-direction as inhomogeneous direction for"
            << " two-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_2 = d_finest_level_dims[2];
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_x = Real(x_hi[0] - x_lo[0]);
        const Real L_y = Real(x_hi[1] - x_lo[1]);
        
        Real* u_avg_local = (Real*)std::malloc(finest_level_dim_2*sizeof(Real));
        
        averaged_quantity.resize(finest_level_dim_2);
        Real* u_avg_global = averaged_quantity.data();
        
        for (int k = 0; k < finest_level_dim_2; k++)
        {
            u_avg_local[k]  = Real(0);
            u_avg_global[k] = Real(0);
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
            
            hier::IntVector ratio_to_coarsest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarsest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarsest_level/ratio_to_coarsest_level;
            
            const int ratio_to_finest_level_2 = ratio_to_finest_level[2];
            
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
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity =
                    getCellData(quantity_name);
                
                Real* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int num_ghosts_2_quantity = num_ghosts_quantity[2];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                const int ghostcell_dim_1_quantity = ghostcell_dims_quantity[1];
                
                const Real weight = Real(dx[0]*dx[1])/(L_x*L_y);
                
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
                                 * Compute the index of the data point and count how many times the data is repeated.
                                 */
                                
                                const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                                
                                int n_overlapped = 1;
                                
                                for (hier::BoxContainer::BoxContainerConstIterator iob(
                                        patch_overlapped_visible_boxes.begin());
                                     iob != patch_overlapped_visible_boxes.end();
                                     iob++)
                                {
                                    const hier::Box& patch_overlapped_visible_box = *iob;
                                    
                                    if (patch_overlapped_visible_box.contains(idx_pt))
                                    {
                                        n_overlapped++;
                                    }
                                }
                                
                                /*
                                 * Compute the linear index and the data to add.
                                 */
                                
                                const int idx = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                    (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                        ghostcell_dim_1_quantity;
                                
                                const Real value_to_add = u[idx]*weight/((Real) n_overlapped);
                                
                                for (int kk = 0; kk < ratio_to_finest_level_2; kk++)
                                {
                                    const int idx_fine = (idx_lo_2 + k)*ratio_to_finest_level_2 + kk;
                                    
                                    u_avg_local[idx_fine] += value_to_add;
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        d_mpi.Allreduce(
            u_avg_local,
            u_avg_global,
            finest_level_dim_2,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_avg_local);
    }
    
    return averaged_quantity;
}


/*
 * Compute averaged reciprocal of value with only z-direction as inhomogeneous direction.
 */
std::vector<Real>
FlowModelMPIHelperAverage::getAveragedReciprocalOfQuantityWithInhomogeneousZDirection(
    const std::string quantity_name,
    const int component_idx,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    std::vector<Real> averaged_reciprocal_quantity;
    
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
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperAverage::getAveragedReciprocalOfQuantityWithInhomogeneousZDirection():\n"
            << "Cannot compute averaged value with only z-direction as inhomogeneous direction for"
            << " one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperAverage::getAveragedReciprocalOfQuantityWithInhomogeneousZDirection():\n"
            << "Cannot compute averaged value with only z-direction as inhomogeneous direction for"
            << " two-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_2 = d_finest_level_dims[2];
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_x = Real(x_hi[0] - x_lo[0]);
        const Real L_y = Real(x_hi[1] - x_lo[1]);
        
        Real* u_inv_avg_local = (Real*)std::malloc(finest_level_dim_2*sizeof(Real));
        
        averaged_reciprocal_quantity.resize(finest_level_dim_2);
        Real* u_inv_avg_global = averaged_reciprocal_quantity.data();
        
        for (int k = 0; k < finest_level_dim_2; k++)
        {
            u_inv_avg_local[k]  = Real(0);
            u_inv_avg_global[k] = Real(0);
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
            
            hier::IntVector ratio_to_coarsest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarsest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarsest_level/ratio_to_coarsest_level;
            
            const int ratio_to_finest_level_2 = ratio_to_finest_level[2];
            
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
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>(quantity_name, hier::IntVector::getZero(d_dim)));
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointer to data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity =
                    getCellData(quantity_name);
                
                Real* u = data_quantity->getPointer(component_idx);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_quantity = data_quantity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_quantity = data_quantity->getGhostBox().numberCells();
                
                const int num_ghosts_0_quantity = num_ghosts_quantity[0];
                const int num_ghosts_1_quantity = num_ghosts_quantity[1];
                const int num_ghosts_2_quantity = num_ghosts_quantity[2];
                const int ghostcell_dim_0_quantity = ghostcell_dims_quantity[0];
                const int ghostcell_dim_1_quantity = ghostcell_dims_quantity[1];
                
                const Real weight = Real(dx[0]*dx[1])/(L_x*L_y);
                
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
                                 * Compute the index of the data point and count how many times the data is repeated.
                                 */
                                
                                const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                                
                                int n_overlapped = 1;
                                
                                for (hier::BoxContainer::BoxContainerConstIterator iob(
                                        patch_overlapped_visible_boxes.begin());
                                     iob != patch_overlapped_visible_boxes.end();
                                     iob++)
                                {
                                    const hier::Box& patch_overlapped_visible_box = *iob;
                                    
                                    if (patch_overlapped_visible_box.contains(idx_pt))
                                    {
                                        n_overlapped++;
                                    }
                                }
                                
                                /*
                                 * Compute the linear index and the data to add.
                                 */
                                
                                const int idx = (relative_idx_lo_0 + i + num_ghosts_0_quantity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_quantity)*ghostcell_dim_0_quantity +
                                    (relative_idx_lo_2 + k + num_ghosts_2_quantity)*ghostcell_dim_0_quantity*
                                        ghostcell_dim_1_quantity;
                                
                                const Real value_to_add = (Real(1)/u[idx])*weight/((Real) n_overlapped);
                                
                                for (int kk = 0; kk < ratio_to_finest_level_2; kk++)
                                {
                                    const int idx_fine = (idx_lo_2 + k)*ratio_to_finest_level_2 + kk;
                                    
                                    u_inv_avg_local[idx_fine] += value_to_add;
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average of reciprocal.
         */
        
        d_mpi.Allreduce(
            u_inv_avg_local,
            u_inv_avg_global,
            finest_level_dim_2,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(u_inv_avg_local);
    }
    
    return averaged_reciprocal_quantity;
}


/*
 * Compute averaged value (on product of variables) with only z-direction as inhomogeneous direction.
 */
std::vector<Real>
FlowModelMPIHelperAverage::getAveragedQuantityWithInhomogeneousZDirection(
    const std::vector<std::string>& quantity_names,
    const std::vector<int>& component_indices,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    int num_quantities = static_cast<int>(quantity_names.size());
    
    TBOX_ASSERT(static_cast<int>(component_indices.size()) == num_quantities);
    
    std::vector<bool> use_reciprocal(num_quantities, false);
    
    return getAveragedQuantityWithInhomogeneousZDirection(
        quantity_names,
        component_indices,
        use_reciprocal,
        data_context);
}


/*
 * Compute averaged value (on product of variables) with only z-direction as inhomogeneous direction.
 */
std::vector<Real>
FlowModelMPIHelperAverage::getAveragedQuantityWithInhomogeneousZDirection(
    const std::vector<std::string>& quantity_names,
    const std::vector<int>& component_indices,
    const std::vector<bool>& use_reciprocal,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    int num_quantities = static_cast<int>(quantity_names.size());
    
    TBOX_ASSERT(static_cast<int>(component_indices.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(use_reciprocal.size()) == num_quantities);
    
    std::vector<Real> averaged_quantity;
    
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
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperAverage::getAveragedQuantityWithInhomogeneousZDirection():\n"
            << "Cannot compute averaged value with only z-direction as inhomogeneous direction for"
            << " one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperAverage::getAveragedQuantityWithInhomogeneousZDirection():\n"
            << "Cannot compute averaged value with only z-direction as inhomogeneous direction for"
            << " two-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_2 = d_finest_level_dims[2];
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_x = Real(x_hi[0] - x_lo[0]);
        const Real L_y = Real(x_hi[1] - x_lo[1]);
        
        Real* avg_local = (Real*)std::malloc(finest_level_dim_2*sizeof(Real));
        
        averaged_quantity.resize(finest_level_dim_2);
        Real* avg_global = averaged_quantity.data();
        
        for (int k = 0; k < finest_level_dim_2; k++)
        {
            avg_local[k]  = Real(0);
            avg_global[k] = Real(0);
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
            
            hier::IntVector ratio_to_coarsest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarsest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarsest_level/ratio_to_coarsest_level;
            
            const int ratio_to_finest_level_2 = ratio_to_finest_level[2];
            
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
                 * Register the patch and data in the flow model and compute the corresponding
                 * average.
                 */
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                }
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = getCellData(quantity_names[qi]);
                }
                
                std::vector<Real*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> num_ghosts_2_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                std::vector<int> ghostcell_dim_1_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                num_ghosts_2_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                ghostcell_dim_1_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    num_ghosts_2_u_qi.push_back(num_ghosts_u_qi[2]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                    ghostcell_dim_1_u_qi.push_back(ghostcell_dims_u_qi[1]);
                }
                
                const Real weight = Real(dx[0]*dx[1])/(L_x*L_y);
                
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
                                 * Compute the index of the data point and count how many times the data is repeated.
                                 */
                                
                                const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                                
                                int n_overlapped = 1;
                                
                                for (hier::BoxContainer::BoxContainerConstIterator iob(
                                        patch_overlapped_visible_boxes.begin());
                                     iob != patch_overlapped_visible_boxes.end();
                                     iob++)
                                {
                                    const hier::Box& patch_overlapped_visible_box = *iob;
                                    
                                    if (patch_overlapped_visible_box.contains(idx_pt))
                                    {
                                        n_overlapped++;
                                    }
                                }
                                
                                /*
                                 * Compute the linear index and the data to add.
                                 */
                                
                                for (int kk = 0; kk < ratio_to_finest_level_2; kk++)
                                {
                                    const int idx_fine = (idx_lo_2 + k)*ratio_to_finest_level_2 + kk;
                                    
                                    Real avg = Real(1);
                                    
                                    for (int qi = 0; qi < num_quantities; qi++)
                                    {
                                        const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        if (use_reciprocal[qi])
                                        {
                                            avg /= u_qi[qi][idx_qi];
                                        }
                                        else
                                        {
                                            avg *= u_qi[qi][idx_qi];
                                        }
                                    }
                                    
                                    avg_local[idx_fine] += (avg*weight/((Real) n_overlapped));
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        d_mpi.Allreduce(
            avg_local,
            avg_global,
            finest_level_dim_2,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(avg_local);
    }
    
    return averaged_quantity;
}


/*
 * Compute averaged value (on product of variable derivatives) with only z direction as inhomogeneous direction.
 */
std::vector<Real>
FlowModelMPIHelperAverage::getAveragedQuantityWithInhomogeneousZDirection(
    const std::vector<std::string>& quantity_names,
    const std::vector<int>& component_indices,
    const std::vector<bool>& use_derivative,
    const std::vector<int>& derivative_directions,
    const int num_ghosts_derivative,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    int num_quantities = static_cast<int>(quantity_names.size());
    
    TBOX_ASSERT(static_cast<int>(component_indices.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(use_derivative.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(derivative_directions.size()) == num_quantities);
    
    std::vector<bool> use_reciprocal(num_quantities, false);
    
    return getAveragedQuantityWithInhomogeneousZDirection(
        quantity_names,
        component_indices,
        use_derivative,
        derivative_directions,
        use_reciprocal,
        num_ghosts_derivative,
        data_context);
}


/*
 * Compute averaged value (on product of variable derivatives) with only z direction as inhomogeneous direction.
 */
std::vector<Real>
FlowModelMPIHelperAverage::getAveragedQuantityWithInhomogeneousZDirection(
    const std::vector<std::string>& quantity_names,
    const std::vector<int>& component_indices,
    const std::vector<bool>& use_derivative,
    const std::vector<int>& derivative_directions,
    const std::vector<bool>& use_reciprocal,
    const int num_ghosts_derivative,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    int num_quantities = static_cast<int>(quantity_names.size());
    
    TBOX_ASSERT(static_cast<int>(component_indices.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(use_derivative.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(derivative_directions.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(use_reciprocal.size()) == num_quantities);
    
    int num_use_derivative = 0;
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperAverage::getAveragedQuantityWithInhomogeneousZDirection():\n"
            << "Cannot compute averaged value with only z-direction as inhomogeneous direction for"
            << " one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperAverage::getAveragedQuantityWithInhomogeneousZDirection():\n"
            << "Cannot compute averaged value with only z-direction as inhomogeneous direction for"
            << " two-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        for (int qi = 0; qi < num_quantities; qi++)
        {
            if (use_derivative[qi])
            {
                num_use_derivative++;
                if (derivative_directions[qi] < 0 || derivative_directions[qi] > 2)
                {
                    TBOX_ERROR(d_object_name
                        << ": FlowModelMPIHelperAverage::getAveragedQuantityWithInhomogeneousZDirection():\n"
                        << "Cannot take derivative for three-dimensional problem!\n"
                        << "derivative_directions[" << qi << "] = " << derivative_directions[qi] << " given!\n"
                        << std::endl);
                }
            }
        }
    }
    
    std::vector<Real> averaged_quantity;
    
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
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(3))
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
        
        hier::IntVector num_ghosts_der = hier::IntVector::getOne(d_dim)*num_ghosts_derivative;
        
        const int finest_level_dim_2 = d_finest_level_dims[2];
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_x = Real(x_hi[0] - x_lo[0]);
        const Real L_y = Real(x_hi[1] - x_lo[1]);
        
        Real* avg_local = (Real*)std::malloc(finest_level_dim_2*sizeof(Real));
        
        averaged_quantity.resize(finest_level_dim_2);
        Real* avg_global = averaged_quantity.data();
        
        for (int k = 0; k < finest_level_dim_2; k++)
        {
            avg_local[k]  = Real(0);
            avg_global[k] = Real(0);
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
            
            hier::IntVector ratio_to_coarsest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarsest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarsest_level/ratio_to_coarsest_level;
            
            const int ratio_to_finest_level_2 = ratio_to_finest_level[2];
            
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
                 * Register the patch and data in the flow model and compute the corresponding
                 * average.
                 */
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    if (num_use_derivative > 0)
                    {
                        num_subghosts_of_data.insert(
                            std::pair<std::string, hier::IntVector>(quantity_names[qi], num_ghosts_der));
                    }
                    else
                    {
                        num_subghosts_of_data.insert(
                            std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                    }
                }
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = getCellData(quantity_names[qi]);
                }
                
                std::vector<Real*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> num_ghosts_2_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                std::vector<int> ghostcell_dim_1_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                num_ghosts_2_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                ghostcell_dim_1_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    num_ghosts_2_u_qi.push_back(num_ghosts_u_qi[2]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                    ghostcell_dim_1_u_qi.push_back(ghostcell_dims_u_qi[1]);
                }
                
                /*
                 * Initialize cell data for the derivatives and get pointers to the cell data.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_derivative;
                std::vector<Real*> der_qi;
                
                if (num_use_derivative > 0)
                {
                    data_derivative = HAMERS_MAKE_SHARED<pdat::CellData<Real> >(patch_box, num_use_derivative, hier::IntVector::getZero(d_dim));
                    
                    der_qi.resize(num_use_derivative);
                    for (int qi = 0; qi < num_use_derivative; qi++)
                    {
                        der_qi[qi] = data_derivative->getPointer(qi);
                    }
                }
                
                const hier::IntVector patch_interior_dims = patch_box.numberCells();
                const int patch_interior_dim_0 = patch_interior_dims[0];
                const int patch_interior_dim_1 = patch_interior_dims[1];
                
                const Real weight = Real(dx[0]*dx[1])/(L_x*L_y);
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    int count_derivative = 0;
                    for (int qi = 0; qi < num_quantities; qi++)
                    {
                        if (use_derivative[qi] && derivative_directions[qi] == 0)
                        {
                            derivative_first_order_x->computeDerivative(
                                data_derivative,
                                data_quantities[qi],
                                Real(dx[0]),
                                patch_visible_box,
                                count_derivative,
                                component_indices[qi]);
                            
                            count_derivative++;
                        }
                        else if (use_derivative[qi] && derivative_directions[qi] == 1)
                        {
                            derivative_first_order_y->computeDerivative(
                                data_derivative,
                                data_quantities[qi],
                                Real(dx[1]),
                                patch_visible_box,
                                count_derivative,
                                component_indices[qi]);
                            
                            count_derivative++;
                        }
                        else if (use_derivative[qi] && derivative_directions[qi] == 2)
                        {
                            derivative_first_order_z->computeDerivative(
                                data_derivative,
                                data_quantities[qi],
                                Real(dx[2]),
                                patch_visible_box,
                                count_derivative,
                                component_indices[qi]);
                            
                            count_derivative++;
                        }
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
                                 * Compute the index of the data point and count how many times the data is repeated.
                                 */
                                
                                const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                                
                                int n_overlapped = 1;
                                
                                for (hier::BoxContainer::BoxContainerConstIterator iob(
                                        patch_overlapped_visible_boxes.begin());
                                     iob != patch_overlapped_visible_boxes.end();
                                     iob++)
                                {
                                    const hier::Box& patch_overlapped_visible_box = *iob;
                                    
                                    if (patch_overlapped_visible_box.contains(idx_pt))
                                    {
                                        n_overlapped++;
                                    }
                                }
                                
                                /*
                                 * Compute the linear index and the data to add.
                                 */
                                
                                Real avg = Real(1);
                                
                                count_derivative = 0;
                                for (int qi = 0; qi < num_quantities; qi++)
                                {
                                    if (use_reciprocal[qi])
                                    {
                                        if (use_derivative[qi])
                                        {
                                            const int idx_der = (relative_idx_lo_0 + i) +
                                                (relative_idx_lo_1 + j)*patch_interior_dim_0 +
                                                (relative_idx_lo_2 + k)*patch_interior_dim_0*
                                                    patch_interior_dim_1;
                                            
                                            avg /= der_qi[count_derivative][idx_der];
                                            
                                            count_derivative++;
                                        }
                                        else
                                        {
                                            const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                                (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                                (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                    ghostcell_dim_1_u_qi[qi];
                                            
                                            avg /= u_qi[qi][idx_qi];
                                        }
                                    }
                                    else
                                    {
                                        if (use_derivative[qi])
                                        {
                                            const int idx_der = (relative_idx_lo_0 + i) +
                                                (relative_idx_lo_1 + j)*patch_interior_dim_0 +
                                                (relative_idx_lo_2 + k)*patch_interior_dim_0*
                                                    patch_interior_dim_1;
                                            
                                            avg *= der_qi[count_derivative][idx_der];
                                            
                                            count_derivative++;
                                        }
                                        else
                                        {
                                            const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                                (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                                (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                    ghostcell_dim_1_u_qi[qi];
                                            
                                            avg *= u_qi[qi][idx_qi];
                                        }
                                    }
                                }
                                
                                /*
                                 * Add the data.
                                 */
                                
                                for (int kk = 0; kk < ratio_to_finest_level_2; kk++)
                                {
                                    const int idx_fine = (idx_lo_2 + k)*ratio_to_finest_level_2 + kk;
                                    
                                    avg_local[idx_fine] += (avg*weight/((Real) n_overlapped));
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        d_mpi.Allreduce(
            avg_local,
            avg_global,
            finest_level_dim_2,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(avg_local);
    }
    
    return averaged_quantity;
}


/*
 * Compute averaged derivative of value (on product of variables) with only z direction as inhomogeneous direction.
 */
std::vector<Real>
FlowModelMPIHelperAverage::getAveragedDerivativeOfQuantityWithInhomogeneousZDirection(
    const std::vector<std::string>& quantity_names,
    const std::vector<int>& component_indices,
    const std::vector<bool>& use_reciprocal,
    const int derivative_direction,
    const int num_ghosts_derivative,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    int num_quantities = static_cast<int>(quantity_names.size());
    
    TBOX_ASSERT(static_cast<int>(component_indices.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(use_reciprocal.size()) == num_quantities);
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperAverage::getAveragedDerivativeOfQuantityWithInhomogeneousZDirection():\n"
            << "Cannot compute averaged value with only z-direction as inhomogeneous direction for"
            << " one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperAverage::getAveragedDerivativeOfQuantityWithInhomogeneousZDirection():\n"
            << "Cannot compute averaged value with only z-direction as inhomogeneous direction for"
            << " two-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(3) && (derivative_direction < 0 || derivative_direction > 2))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperAverage::getAveragedDerivativeOfQuantityWithInhomogeneousZDirection():\n"
            << "Cannot take derivative for three-dimensional problem!\n"
            << "derivative_direction = " << derivative_direction << " given!\n"
            << std::endl);
    }
    
    hier::IntVector num_ghosts_der = hier::IntVector::getZero(d_dim);
    num_ghosts_der[derivative_direction] = num_ghosts_derivative;
    
    std::vector<Real> averaged_derivative;
    
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
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(3))
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
        
        const int finest_level_dim_2 = d_finest_level_dims[2];
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_x = Real(x_hi[0] - x_lo[0]);
        const Real L_y = Real(x_hi[1] - x_lo[1]);
        
        Real* der_avg_local = (Real*)std::malloc(finest_level_dim_2*sizeof(Real));
        
        averaged_derivative.resize(finest_level_dim_2);
        Real* der_avg_global = averaged_derivative.data();
        
        for (int k = 0; k < finest_level_dim_2; k++)
        {
            der_avg_local[k]  = Real(0);
            der_avg_global[k] = Real(0);
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
            
            hier::IntVector ratio_to_coarsest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarsest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarsest_level/ratio_to_coarsest_level;
            
            const int ratio_to_finest_level_2 = ratio_to_finest_level[2];
            
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
                 * Register the patch and the quantities in the flow model and compute the
                 * corresponding product.
                 */
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model->getNumberOfGhostCells();
                TBOX_ASSERT(num_ghosts >= num_ghosts_der);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], num_ghosts));
                }
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = getCellData(quantity_names[qi]);
                }
                
                std::vector<Real*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> num_ghosts_2_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                std::vector<int> ghostcell_dim_1_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                num_ghosts_2_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                ghostcell_dim_1_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    num_ghosts_2_u_qi.push_back(num_ghosts_u_qi[2]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                    ghostcell_dim_1_u_qi.push_back(ghostcell_dims_u_qi[1]);
                }
                
                /*
                 * Initialize cell data for the derivative and the product of variables.
                 * Also, get pointers to the cell data containers.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_derivative(
                    new pdat::CellData<Real>(patch_box, 1, hier::IntVector::getZero(d_dim)));
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_product(
                    new pdat::CellData<Real>(patch_box, 1, num_ghosts_der));
                
                Real* der     = data_derivative->getPointer(0);
                Real* product = data_product->getPointer(0);
                
                const hier::IntVector num_ghosts_product = data_product->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_product = data_product->getGhostBox().numberCells();
                
                const int num_ghosts_0_product = num_ghosts_product[0];
                const int num_ghosts_1_product = num_ghosts_product[1];
                const int num_ghosts_2_product = num_ghosts_product[2];
                const int ghostcell_dim_0_product = ghostcell_dims_product[0];
                const int ghostcell_dim_1_product = ghostcell_dims_product[1];
                
                const hier::IntVector patch_interior_dims = patch_box.numberCells();
                const int patch_interior_dim_0 = patch_interior_dims[0];
                const int patch_interior_dim_1 = patch_interior_dims[1];
                const int patch_interior_dim_2 = patch_interior_dims[2];
                
                data_product->fillAll(Real(1));
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    if (use_reciprocal[qi])
                    {
                        for (int k = -num_ghosts_2_product; k < patch_interior_dim_2 + num_ghosts_2_product; k++)
                        {
                            for (int j = -num_ghosts_1_product; j < patch_interior_dim_1 + num_ghosts_1_product; j++)
                            {
                                for (int i = -num_ghosts_0_product; i < patch_interior_dim_0 + num_ghosts_0_product; i++)
                                {
                                    const int idx_product = (i + num_ghosts_0_product) +
                                        (j + num_ghosts_1_product)*ghostcell_dim_0_product +
                                        (k + num_ghosts_2_product)*ghostcell_dim_0_product*
                                            ghostcell_dim_1_product;
                                    
                                    const int idx_u_qi = (i + num_ghosts_0_u_qi[qi]) +
                                        (j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                        (k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                            ghostcell_dim_1_u_qi[qi];
                                        
                                    product[idx_product] /= u_qi[qi][idx_u_qi];
                                }
                            }
                        }
                    }
                    else
                    {
                        for (int k = -num_ghosts_2_product; k < patch_interior_dim_2 + num_ghosts_2_product; k++)
                        {
                            for (int j = -num_ghosts_1_product; j < patch_interior_dim_1 + num_ghosts_1_product; j++)
                            {
                                for (int i = -num_ghosts_0_product; i < patch_interior_dim_0 + num_ghosts_0_product; i++)
                                {
                                    const int idx_product = (i + num_ghosts_0_product) +
                                        (j + num_ghosts_1_product)*ghostcell_dim_0_product +
                                        (k + num_ghosts_2_product)*ghostcell_dim_0_product*
                                            ghostcell_dim_1_product;
                                    
                                    const int idx_u_qi = (i + num_ghosts_0_u_qi[qi]) +
                                        (j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                        (k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                            ghostcell_dim_1_u_qi[qi];
                                        
                                    product[idx_product] *= u_qi[qi][idx_u_qi];
                                }
                            }
                        }
                    }
                }
                
                const Real weight = Real(dx[0]*dx[1])/(L_x*L_y);
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    if (derivative_direction == 0)
                    {
                        derivative_first_order_x->computeDerivative(
                            data_derivative,
                            data_product,
                            Real(dx[0]),
                            patch_visible_box,
                            0,
                            0);
                    }
                    else if (derivative_direction == 1)
                    {
                        derivative_first_order_y->computeDerivative(
                            data_derivative,
                            data_product,
                            Real(dx[1]),
                            patch_visible_box,
                            0,
                            0);
                    }
                    else if (derivative_direction == 2)
                    {
                         derivative_first_order_z->computeDerivative(
                            data_derivative,
                            data_product,
                            Real(dx[2]),
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
                                 * Compute the index of the data point and count how many times the data is repeated.
                                 */
                                
                                const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                                
                                int n_overlapped = 1;
                                
                                for (hier::BoxContainer::BoxContainerConstIterator iob(
                                        patch_overlapped_visible_boxes.begin());
                                     iob != patch_overlapped_visible_boxes.end();
                                     iob++)
                                {
                                    const hier::Box& patch_overlapped_visible_box = *iob;
                                    
                                    if (patch_overlapped_visible_box.contains(idx_pt))
                                    {
                                        n_overlapped++;
                                    }
                                }
                                
                                /*
                                 * Compute the linear index and the data to add.
                                 */
                                
                                const int idx = (relative_idx_lo_0 + i) +
                                    (relative_idx_lo_1 + j)*patch_interior_dim_0 +
                                    (relative_idx_lo_2 + k)*patch_interior_dim_0*
                                        patch_interior_dim_1;
                                
                                const Real value_to_add = der[idx]*weight/((Real) n_overlapped);
                                
                                /*
                                 * Add the data.
                                 */
                                
                                for (int kk = 0; kk < ratio_to_finest_level_2; kk++)
                                {
                                    const int idx_fine = (idx_lo_2 + k)*ratio_to_finest_level_2 + kk;
                                    
                                    der_avg_local[idx_fine] += value_to_add;
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average of derivative.
         */
        
        d_mpi.Allreduce(
            der_avg_local,
            der_avg_global,
            finest_level_dim_2,
            MPI_DOUBLE,
            MPI_SUM);
        
        std::free(der_avg_local);
    }
    
    return averaged_derivative;
}
