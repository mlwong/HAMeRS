#include "util/MPI_helpers/MPIHelperMaxMin.hpp"

#include "extn/patch_hierarchies/ExtendedFlattenedHierarchy.hpp"

#include "SAMRAI/geom/CartesianPatchGeometry.h"

#include <limits>

/*
 * Compute maximum value over the entire domain.
 */
Real
MPIHelperMaxMin::getMaxQuantity(
    HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_quantity,
    const int component_idx,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
    Real u_max_local  = std::numeric_limits<Real>::min();
    Real u_max_global = std::numeric_limits<Real>::min();
    
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
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity(
                    HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
                        patch->getPatchData(variable_quantity, data_context)));
                
                Real* u = data_quantity->getPointer(component_idx);
                
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
                    
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        /*
                         * Compute the linear index and update the max.
                         */
                        
                        const int idx = relative_idx_lo_0 + i + num_ghosts_0_quantity;
                        
                        u_max_local = std::max(u_max_local, u[idx]);
                    }
                }
            }
        }
        
        /*
         * Reduction to get the global max.
         */
        
        d_mpi.Allreduce(
            &u_max_local,
            &u_max_global,
            1,
            HAMERS_MPI_REAL,
            MPI_MAX);
    }
    else if (d_dim == tbox::Dimension(2))
    {
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
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity(
                    HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
                        patch->getPatchData(variable_quantity, data_context)));
                
                Real* u = data_quantity->getPointer(component_idx);
                
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
                            
                            u_max_local = std::max(u_max_local, u[idx]);
                        }
                    }
                }
            }
        }
        
        /*
         * Reduction to get the global max.
         */
        
        d_mpi.Allreduce(
            &u_max_local,
            &u_max_global,
            1,
            HAMERS_MPI_REAL,
            MPI_MAX);
    }
    else if (d_dim == tbox::Dimension(3))
    {
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
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity(
                    HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
                        patch->getPatchData(variable_quantity, data_context)));
                
                Real* u = data_quantity->getPointer(component_idx);
                
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
                                
                                u_max_local = std::max(u_max_local, u[idx]);
                            }
                        }
                    }
                }
            }
        }
        
        /*
         * Reduction to get the global max.
         */
        
        d_mpi.Allreduce(
            &u_max_local,
            &u_max_global,
            1,
            HAMERS_MPI_REAL,
            MPI_MAX);
    }
    
    return u_max_global;
}


/*
 * Compute minimum value over the entire domain.
 */
Real
MPIHelperMaxMin::getMinQuantity(
    HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_quantity,
    const int component_idx,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
    Real u_min_local  = std::numeric_limits<Real>::max();
    Real u_min_global = std::numeric_limits<Real>::max();
    
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
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity(
                    HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
                        patch->getPatchData(variable_quantity, data_context)));
                
                Real* u = data_quantity->getPointer(component_idx);
                
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
                    
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        /*
                         * Compute the linear index and update the min.
                         */
                        
                        const int idx = relative_idx_lo_0 + i + num_ghosts_0_quantity;
                        
                        u_min_local = std::min(u_min_local, u[idx]);
                    }
                }
            }
        }
        
        /*
         * Reduction to get the global min.
         */
        
        d_mpi.Allreduce(
            &u_min_local,
            &u_min_global,
            1,
            HAMERS_MPI_REAL,
            MPI_MIN);
    }
    else if (d_dim == tbox::Dimension(2))
    {
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
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity(
                    HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
                        patch->getPatchData(variable_quantity, data_context)));
                
                Real* u = data_quantity->getPointer(component_idx);
                
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
                            
                            u_min_local = std::min(u_min_local, u[idx]);
                        }
                    }
                }
            }
        }
        
        /*
         * Reduction to get the global min.
         */
        
        d_mpi.Allreduce(
            &u_min_local,
            &u_min_global,
            1,
            HAMERS_MPI_REAL,
            MPI_MIN);
    }
    else if (d_dim == tbox::Dimension(3))
    {
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
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity(
                    HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
                        patch->getPatchData(variable_quantity, data_context)));
                
                Real* u = data_quantity->getPointer(component_idx);
                
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
                                
                                u_min_local = std::min(u_min_local, u[idx]);
                            }
                        }
                    }
                }
            }
        }
        
        /*
         * Reduction to get the global min.
         */
        
        d_mpi.Allreduce(
            &u_min_local,
            &u_min_global,
            1,
            HAMERS_MPI_REAL,
            MPI_MIN);
    }
    
    return u_min_global;
}


/*
 * Compute maximum value with only x-direction as inhomogeneous direction.
 */
std::vector<Real>
MPIHelperMaxMin::getMaxQuantityWithInhomogeneousXDirection(
    HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_quantity,
    const int component_idx,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
    std::vector<Real> max_quantity;
    
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
        
        Real* u_max_local = (Real*)std::malloc(finest_level_dim_0*sizeof(Real));
        
        max_quantity.resize(finest_level_dim_0);
        Real* u_max_global = max_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_max_local[i]  = std::numeric_limits<Real>::min();
            u_max_global[i] = std::numeric_limits<Real>::min();
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
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity(
                    HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
                        patch->getPatchData(variable_quantity, data_context)));
                
                Real* u = data_quantity->getPointer(component_idx);
                
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
                            
                            u_max_local[idx_fine] = std::max(u_max_local[idx_fine], u[idx]);
                        }
                    }
                }
            }
        }
        
        /*
         * Reduction to get the global max.
         */
        
        d_mpi.Allreduce(
            u_max_local,
            u_max_global,
            finest_level_dim_0,
            HAMERS_MPI_REAL,
            MPI_MAX);
        
        std::free(u_max_local);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_0 = d_finest_level_dims[0];
        
        Real* u_max_local = (Real*)std::malloc(finest_level_dim_0*sizeof(Real));
        
        max_quantity.resize(finest_level_dim_0);
        Real* u_max_global = max_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_max_local[i]  = std::numeric_limits<Real>::min();
            u_max_global[i] = std::numeric_limits<Real>::min();
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
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity(
                    HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
                        patch->getPatchData(variable_quantity, data_context)));
                
                Real* u = data_quantity->getPointer(component_idx);
                
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
                    // const int idx_lo_1 = index_lo[1];
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
                                
                                u_max_local[idx_fine] = std::max(u_max_local[idx_fine], u[idx]);
                            }
                        }
                    }
                }
            }
        }
        
        /*
         * Reduction to get the global max.
         */
        
        d_mpi.Allreduce(
            u_max_local,
            u_max_global,
            finest_level_dim_0,
            HAMERS_MPI_REAL,
            MPI_MAX);
        
        std::free(u_max_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_0 = d_finest_level_dims[0];
        
        Real* u_max_local = (Real*)std::malloc(finest_level_dim_0*sizeof(Real));
        
        max_quantity.resize(finest_level_dim_0);
        Real* u_max_global = max_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_max_local[i]  = std::numeric_limits<Real>::min();
            u_max_global[i] = std::numeric_limits<Real>::min();
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
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity(
                    HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
                        patch->getPatchData(variable_quantity, data_context)));
                
                Real* u = data_quantity->getPointer(component_idx);
                
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
                    // const int idx_lo_1 = index_lo[1];
                    // const int idx_lo_2 = index_lo[2];
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
                                    
                                    u_max_local[idx_fine] = std::max(u_max_local[idx_fine], u[idx]);
                                }
                            }
                        }
                    }
                }
            }
        }
        
        /*
         * Reduction to get the global max.
         */
        
        d_mpi.Allreduce(
            u_max_local,
            u_max_global,
            finest_level_dim_0,
            HAMERS_MPI_REAL,
            MPI_MAX);
        
        std::free(u_max_local);
    }
    
    return max_quantity;
}


/*
 * Compute minimum value with only x-direction as inhomogeneous direction.
 */
std::vector<Real>
MPIHelperMaxMin::getMinQuantityWithInhomogeneousXDirection(
    HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_quantity,
    const int component_idx,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
    std::vector<Real> min_quantity;
    
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
        
        Real* u_min_local = (Real*)std::malloc(finest_level_dim_0*sizeof(Real));
        
        min_quantity.resize(finest_level_dim_0);
        Real* u_min_global = min_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_min_local[i]  = std::numeric_limits<Real>::max();
            u_min_global[i] = std::numeric_limits<Real>::max();
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
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity(
                    HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
                        patch->getPatchData(variable_quantity, data_context)));
                
                Real* u = data_quantity->getPointer(component_idx);
                
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
                            
                            u_min_local[idx_fine] = std::min(u_min_local[idx_fine], u[idx]);
                        }
                    }
                }
            }
        }
        
        /*
         * Reduction to get the global min.
         */
        
        d_mpi.Allreduce(
            u_min_local,
            u_min_global,
            finest_level_dim_0,
            HAMERS_MPI_REAL,
            MPI_MIN);
        
        std::free(u_min_local);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_0 = d_finest_level_dims[0];
        
        Real* u_min_local = (Real*)std::malloc(finest_level_dim_0*sizeof(Real));
        
        min_quantity.resize(finest_level_dim_0);
        Real* u_min_global = min_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_min_local[i]  = std::numeric_limits<Real>::max();
            u_min_global[i] = std::numeric_limits<Real>::max();
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
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity(
                    HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
                        patch->getPatchData(variable_quantity, data_context)));
                
                Real* u = data_quantity->getPointer(component_idx);
                
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
                    // const int idx_lo_1 = index_lo[1];
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
                                
                                u_min_local[idx_fine] = std::min(u_min_local[idx_fine], u[idx]);
                            }
                        }
                    }
                }
            }
        }
        
        /*
         * Reduction to get the global min.
         */
        
        d_mpi.Allreduce(
            u_min_local,
            u_min_global,
            finest_level_dim_0,
            HAMERS_MPI_REAL,
            MPI_MIN);
        
        std::free(u_min_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_0 = d_finest_level_dims[0];
        
        Real* u_min_local = (Real*)std::malloc(finest_level_dim_0*sizeof(Real));
        
        min_quantity.resize(finest_level_dim_0);
        Real* u_min_global = min_quantity.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            u_min_local[i]  = std::numeric_limits<Real>::max();
            u_min_global[i] = std::numeric_limits<Real>::max();
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
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity(
                    HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
                        patch->getPatchData(variable_quantity, data_context)));
                
                Real* u = data_quantity->getPointer(component_idx);
                
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
                    // const int idx_lo_1 = index_lo[1];
                    // const int idx_lo_2 = index_lo[2];
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
                                    
                                    u_min_local[idx_fine] = std::min(u_min_local[idx_fine], u[idx]);
                                }
                            }
                        }
                    }
                }
            }
        }
        
        /*
         * Reduction to get the global min.
         */
        
        d_mpi.Allreduce(
            u_min_local,
            u_min_global,
            finest_level_dim_0,
            HAMERS_MPI_REAL,
            MPI_MIN);
        
        std::free(u_min_local);
    }
    
    return min_quantity;
}


/*
 * Compute maximum value with only y-direction as inhomogeneous direction.
 */
std::vector<Real>
MPIHelperMaxMin::getMaxQuantityWithInhomogeneousYDirection(
    HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_quantity,
    const int component_idx,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
    std::vector<Real> max_quantity;
    
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
            << ": MPIHelperMaxMin::getMaxQuantityWithInhomogeneousYDirection():\n"
            << "Cannot compute maximum value with only y-direction as inhomogeneous direction for"
            << " one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_1 = d_finest_level_dims[1];
        
        Real* u_max_local = (Real*)std::malloc(finest_level_dim_1*sizeof(Real));
        
        max_quantity.resize(finest_level_dim_1);
        Real* u_max_global = max_quantity.data();
        
        for (int j = 0; j < finest_level_dim_1; j++)
        {
            u_max_local[j]  = std::numeric_limits<Real>::min();
            u_max_global[j] = std::numeric_limits<Real>::min();
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
                 * Get the patch lower indices.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity(
                    HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
                        patch->getPatchData(variable_quantity, data_context)));
                
                Real* u = data_quantity->getPointer(component_idx);
                
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
                    
                    // const int idx_lo_0 = index_lo[0];
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
                                
                                u_max_local[idx_fine] = std::max(u_max_local[idx_fine], u[idx]);
                            }
                        }
                    }
                }
            }
        }
        
        /*
         * Reduction to get the global max.
         */
        
        d_mpi.Allreduce(
            u_max_local,
            u_max_global,
            finest_level_dim_1,
            HAMERS_MPI_REAL,
            MPI_MAX);
        
        std::free(u_max_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_1 = d_finest_level_dims[1];
        
        Real* u_max_local = (Real*)std::malloc(finest_level_dim_1*sizeof(Real));
        
        max_quantity.resize(finest_level_dim_1);
        Real* u_max_global = max_quantity.data();
        
        for (int j = 0; j < finest_level_dim_1; j++)
        {
            u_max_local[j]  = std::numeric_limits<Real>::min();
            u_max_global[j] = std::numeric_limits<Real>::min();
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
                 * Get the patch lower indices.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity(
                    HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
                        patch->getPatchData(variable_quantity, data_context)));
                
                Real* u = data_quantity->getPointer(component_idx);
                
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
                    
                    // const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    // const int idx_lo_2 = index_lo[2];
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
                                    
                                    u_max_local[idx_fine] = std::max(u_max_local[idx_fine], u[idx]);
                                }
                            }
                        }
                    }
                }
            }
        }
        
        /*
         * Reduction to get the global max.
         */
        
        d_mpi.Allreduce(
            u_max_local,
            u_max_global,
            finest_level_dim_1,
            HAMERS_MPI_REAL,
            MPI_MAX);
        
        std::free(u_max_local);
    }
    
    return max_quantity;
}


/*
 * Compute minimum value with only y-direction as inhomogeneous direction.
 */
std::vector<Real>
MPIHelperMaxMin::getMinQuantityWithInhomogeneousYDirection(
    HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_quantity,
    const int component_idx,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
    std::vector<Real> min_quantity;
    
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
            << ": MPIHelperMaxMin::getMinQuantityWithInhomogeneousYDirection():\n"
            << "Cannot compute minimum value with only y-direction as inhomogeneous direction for"
            << " one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_1 = d_finest_level_dims[1];
        
        Real* u_min_local = (Real*)std::malloc(finest_level_dim_1*sizeof(Real));
        
        min_quantity.resize(finest_level_dim_1);
        Real* u_min_global = min_quantity.data();
        
        for (int j = 0; j < finest_level_dim_1; j++)
        {
            u_min_local[j]  = std::numeric_limits<Real>::max();
            u_min_global[j] = std::numeric_limits<Real>::max();
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
                 * Get the patch lower indices.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity(
                    HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
                        patch->getPatchData(variable_quantity, data_context)));
                
                Real* u = data_quantity->getPointer(component_idx);
                
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
                    
                    // const int idx_lo_0 = index_lo[0];
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
                                
                                u_min_local[idx_fine] = std::min(u_min_local[idx_fine], u[idx]);
                            }
                        }
                    }
                }
            }
        }
        
        /*
         * Reduction to get the global min.
         */
        
        d_mpi.Allreduce(
            u_min_local,
            u_min_global,
            finest_level_dim_1,
            HAMERS_MPI_REAL,
            MPI_MIN);
        
        std::free(u_min_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_1 = d_finest_level_dims[1];
        
        Real* u_min_local = (Real*)std::malloc(finest_level_dim_1*sizeof(Real));
        
        min_quantity.resize(finest_level_dim_1);
        Real* u_min_global = min_quantity.data();
        
        for (int j = 0; j < finest_level_dim_1; j++)
        {
            u_min_local[j]  = std::numeric_limits<Real>::max();
            u_min_global[j] = std::numeric_limits<Real>::max();
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
                 * Get the patch lower indices.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity(
                    HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
                        patch->getPatchData(variable_quantity, data_context)));
                
                Real* u = data_quantity->getPointer(component_idx);
                
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
                    
                    // const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    // const int idx_lo_2 = index_lo[2];
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
                                    
                                    u_min_local[idx_fine] = std::min(u_min_local[idx_fine], u[idx]);
                                }
                            }
                        }
                    }
                }
            }
        }
        
        /*
         * Reduction to get the global min.
         */
        
        d_mpi.Allreduce(
            u_min_local,
            u_min_global,
            finest_level_dim_1,
            HAMERS_MPI_REAL,
            MPI_MIN);
        
        std::free(u_min_local);
    }
    
    return min_quantity;
}


/*
 * Compute maximum value with only z-direction as inhomogeneous direction.
 */
std::vector<Real>
MPIHelperMaxMin::getMaxQuantityWithInhomogeneousZDirection(
    HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_quantity,
    const int component_idx,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
    std::vector<Real> max_quantity;
    
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
            << ": MPIHelperMaxMin::getMaxQuantityWithInhomogeneousZDirection():\n"
            << "Cannot compute maximum value with only z-direction as inhomogeneous direction for"
            << " one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": MPIHelperMaxMin::getMaxQuantityWithInhomogeneousZDirection():\n"
            << "Cannot compute maximum value with only z-direction as inhomogeneous direction for"
            << " two-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_2 = d_finest_level_dims[2];
        
        Real* u_max_local = (Real*)std::malloc(finest_level_dim_2*sizeof(Real));
        
        max_quantity.resize(finest_level_dim_2);
        Real* u_max_global = max_quantity.data();
        
        for (int k = 0; k < finest_level_dim_2; k++)
        {
            u_max_local[k]  = std::numeric_limits<Real>::min();
            u_max_global[k] = std::numeric_limits<Real>::min();
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
                 * Get the patch lower indices.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity(
                    HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
                        patch->getPatchData(variable_quantity, data_context)));
                
                Real* u = data_quantity->getPointer(component_idx);
                
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
                    
                    // const int idx_lo_0 = index_lo[0];
                    // const int idx_lo_1 = index_lo[1];
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
                                    
                                    u_max_local[idx_fine] = std::max(u_max_local[idx_fine], u[idx]);
                                }
                            }
                        }
                    }
                }
            }
        }
        
        /*
         * Reduction to get the global max.
         */
        
        d_mpi.Allreduce(
            u_max_local,
            u_max_global,
            finest_level_dim_2,
            HAMERS_MPI_REAL,
            MPI_MAX);
        
        std::free(u_max_local);
    }
    
    return max_quantity;
}


/*
 * Compute minimum value with only z-direction as inhomogeneous direction.
 */
std::vector<Real>
MPIHelperMaxMin::getMinQuantityWithInhomogeneousZDirection(
    HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_quantity,
    const int component_idx,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
    std::vector<Real> min_quantity;
    
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
            << ": MPIHelperMaxMin::getMinQuantityWithInhomogeneousZDirection():\n"
            << "Cannot compute minimum value with only z-direction as inhomogeneous direction for"
            << " one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": MPIHelperMaxMin::getMinQuantityWithInhomogeneousZDirection():\n"
            << "Cannot compute minimum value with only z-direction as inhomogeneous direction for"
            << " two-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_2 = d_finest_level_dims[2];
        
        Real* u_min_local = (Real*)std::malloc(finest_level_dim_2*sizeof(Real));
        
        min_quantity.resize(finest_level_dim_2);
        Real* u_min_global = min_quantity.data();
        
        for (int k = 0; k < finest_level_dim_2; k++)
        {
            u_min_local[k]  = std::numeric_limits<Real>::max();
            u_min_global[k] = std::numeric_limits<Real>::max();
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
                 * Get the patch lower indices.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity(
                    HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
                        patch->getPatchData(variable_quantity, data_context)));
                
                Real* u = data_quantity->getPointer(component_idx);
                
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
                    
                    // const int idx_lo_0 = index_lo[0];
                    // const int idx_lo_1 = index_lo[1];
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
                                    
                                    u_min_local[idx_fine] = std::min(u_min_local[idx_fine], u[idx]);
                                }
                            }
                        }
                    }
                }
            }
        }
        
        /*
         * Reduction to get the global min.
         */
        
        d_mpi.Allreduce(
            u_min_local,
            u_min_global,
            finest_level_dim_2,
            HAMERS_MPI_REAL,
            MPI_MIN);
        
        std::free(u_min_local);
    }
    
    return min_quantity;
}


/*
 * Compute maximum location within quantity bounds in x-direction.
 */
Real
MPIHelperMaxMin::getMaxLocationWithinQuantityBoundsInXDirection(
    HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_quantity,
    const int component_idx,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const Real bound_lo,
    const Real bound_hi) const
{
    Real location_x_max_global;
    
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
        Real location_x_max_local = Real(x_lo[0]);
        location_x_max_global     = Real(x_lo[0]);
        
        const Real half = Real(1)/Real(2);
        
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
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity(
                    HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
                        patch->getPatchData(variable_quantity, data_context)));
                
                Real* u = data_quantity->getPointer(component_idx);
                
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
                    
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        /*
                         * Compute the linear index and update the max.
                         */
                        
                        const int idx = relative_idx_lo_0 + i + num_ghosts_0_quantity;
                        
                        const Real x = (Real(relative_idx_lo_0 + i) + half)*Real(dx[0]) + Real(x_lo_patch[0]);
                        
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
         * Reduction to get the global max.
         */
        
        d_mpi.Allreduce(
            &location_x_max_local,
            &location_x_max_global,
            1,
            HAMERS_MPI_REAL,
            MPI_MAX);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        Real location_x_max_local = Real(x_lo[0]);
        location_x_max_global     = Real(x_lo[0]);
        
        const Real half = Real(1)/Real(2);
        
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
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity(
                    HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
                        patch->getPatchData(variable_quantity, data_context)));
                
                Real* u = data_quantity->getPointer(component_idx);
                
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
                            
                            const Real x = (Real(relative_idx_lo_0 + i) + half)*Real(dx[0]) + Real(x_lo_patch[0]);
                            
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
        }
        
        /*
         * Reduction to get the global max.
         */
        
        d_mpi.Allreduce(
            &location_x_max_local,
            &location_x_max_global,
            1,
            HAMERS_MPI_REAL,
            MPI_MAX);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        Real location_x_max_local = Real(x_lo[0]);
        location_x_max_global     = Real(x_lo[0]);
        
        const Real half = Real(1)/Real(2);
        
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
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity(
                    HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
                        patch->getPatchData(variable_quantity, data_context)));
                
                Real* u = data_quantity->getPointer(component_idx);
                
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
                                
                                const Real x = (Real(relative_idx_lo_0 + i) + half)*Real(dx[0]) + Real(x_lo_patch[0]);
                                
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
            }
        }
        
        /*
         * Reduction to get the global max.
         */
        
        d_mpi.Allreduce(
            &location_x_max_local,
            &location_x_max_global,
            1,
            HAMERS_MPI_REAL,
            MPI_MAX);
    }
    
    return location_x_max_global;
}


/*
 * Compute minimum location within quantity bounds in x-direction.
 */
Real
MPIHelperMaxMin::getMinLocationWithinQuantityBoundsInXDirection(
    HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_quantity,
    const int component_idx,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const Real bound_lo,
    const Real bound_hi) const
{
    Real location_x_min_global;
    
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
        Real location_x_min_local = Real(x_hi[0]);
        location_x_min_global     = Real(x_hi[0]);
        
        const Real half = Real(1)/Real(2);
        
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
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity(
                    HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
                        patch->getPatchData(variable_quantity, data_context)));
                
                Real* u = data_quantity->getPointer(component_idx);
                
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
                    
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        /*
                         * Compute the linear index and update the max.
                         */
                        
                        const int idx = relative_idx_lo_0 + i + num_ghosts_0_quantity;
                        
                        const Real x = (Real(relative_idx_lo_0 + i) + half)*Real(dx[0]) + Real(x_lo_patch[0]);
                        
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
         * Reduction to get the global min.
         */
        
        d_mpi.Allreduce(
            &location_x_min_local,
            &location_x_min_global,
            1,
            HAMERS_MPI_REAL,
            MPI_MIN);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        Real location_x_min_local = Real(x_hi[0]);
        location_x_min_global     = Real(x_hi[0]);
        
        const Real half = Real(1)/Real(2);
        
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
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity(
                    HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
                        patch->getPatchData(variable_quantity, data_context)));
                
                Real* u = data_quantity->getPointer(component_idx);
                
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
                            
                            const Real x = (Real(relative_idx_lo_0 + i) + half)*Real(dx[0]) + Real(x_lo_patch[0]);
                            
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
        }
        
        /*
         * Reduction to get the global min.
         */
        
        d_mpi.Allreduce(
            &location_x_min_local,
            &location_x_min_global,
            1,
            HAMERS_MPI_REAL,
            MPI_MIN);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        Real location_x_min_local = Real(x_hi[0]);
        location_x_min_global     = Real(x_hi[0]);
        
        const Real half = Real(1)/Real(2);
        
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
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity(
                    HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
                        patch->getPatchData(variable_quantity, data_context)));
                
                Real* u = data_quantity->getPointer(component_idx);
                
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
                                
                                const Real x = (Real(relative_idx_lo_0 + i) + half)*Real(dx[0]) + Real(x_lo_patch[0]);
                                
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
            }
        }
        
        /*
         * Reduction to get the global min.
         */
        
        d_mpi.Allreduce(
            &location_x_min_local,
            &location_x_min_global,
            1,
            HAMERS_MPI_REAL,
            MPI_MIN);
    }
    
    return location_x_min_global;
}


/*
 * Compute maximum location within quantity bounds in y-direction.
 */
Real
MPIHelperMaxMin::getMaxLocationWithinQuantityBoundsInYDirection(
    HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_quantity,
    const int component_idx,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const Real bound_lo,
    const Real bound_hi) const
{
    Real location_y_max_global;
    
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
            << ": MPIHelperMaxMin::getMaxLocationWithinQuantityBoundsInYDirection():\n"
            << "Cannot compute maximum location within quantity bounds in y-direction for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        Real location_y_max_local = Real(x_lo[1]);
        location_y_max_global     = Real(x_lo[1]);
        
        const Real half = Real(1)/Real(2);
        
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
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity(
                    HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
                        patch->getPatchData(variable_quantity, data_context)));
                
                Real* u = data_quantity->getPointer(component_idx);
                
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
                            
                            const Real y = (Real(relative_idx_lo_1 + j) + half)*Real(dx[1]) + Real(x_lo_patch[1]);
                            
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
        }
        
        /*
         * Reduction to get the global max.
         */
        
        d_mpi.Allreduce(
            &location_y_max_local,
            &location_y_max_global,
            1,
            HAMERS_MPI_REAL,
            MPI_MAX);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        Real location_y_max_local = Real(x_lo[1]);
        location_y_max_global     = Real(x_lo[1]);
        
        const Real half = Real(1)/Real(2);
        
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
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity(
                    HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
                        patch->getPatchData(variable_quantity, data_context)));
                
                Real* u = data_quantity->getPointer(component_idx);
                
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
                                
                                const Real y = (Real(relative_idx_lo_1 + j) + half)*Real(dx[1]) + Real(x_lo_patch[1]);
                                
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
            }
        }
        
        /*
         * Reduction to get the global max.
         */
        
        d_mpi.Allreduce(
            &location_y_max_local,
            &location_y_max_global,
            1,
            HAMERS_MPI_REAL,
            MPI_MAX);
    }
    
    return location_y_max_global;
}


/*
 * Compute minimum location within quantity bounds in y-direction.
 */
Real
MPIHelperMaxMin::getMinLocationWithinQuantityBoundsInYDirection(
    HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_quantity,
    const int component_idx,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const Real bound_lo,
    const Real bound_hi) const
{
    Real location_y_min_global;
    
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
            << ": MPIHelperMaxMin::getMinLocationWithinQuantityBoundsInYDirection():\n"
            << "Cannot compute maximum location within quantity bounds in y-direction for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        Real location_y_min_local = Real(x_hi[1]);
        location_y_min_global     = Real(x_hi[1]);
        
        const Real half = Real(1)/Real(2);
        
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
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity(
                    HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
                        patch->getPatchData(variable_quantity, data_context)));
                
                Real* u = data_quantity->getPointer(component_idx);
                
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
                            
                            const Real y = (Real(relative_idx_lo_1 + j) + half)*Real(dx[1]) + Real(x_lo_patch[1]);
                            
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
        }
        
        /*
         * Reduction to get the global min.
         */
        
        d_mpi.Allreduce(
            &location_y_min_local,
            &location_y_min_global,
            1,
            HAMERS_MPI_REAL,
            MPI_MIN);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        Real location_y_min_local = Real(x_hi[1]);
        location_y_min_global     = Real(x_hi[1]);
        
        const Real half = Real(1)/Real(2);
        
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
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity(
                    HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
                        patch->getPatchData(variable_quantity, data_context)));
                
                Real* u = data_quantity->getPointer(component_idx);
                
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
                                
                                const Real y = (Real(relative_idx_lo_1 + j) + half)*Real(dx[1]) + Real(x_lo_patch[1]);
                                
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
            }
        }
        
        /*
         * Reduction to get the global min.
         */
        
        d_mpi.Allreduce(
            &location_y_min_local,
            &location_y_min_global,
            1,
            HAMERS_MPI_REAL,
            MPI_MIN);
    }
    
    return location_y_min_global;
}


/*
 * Compute maximum location within quantity bounds in z-direction.
 */
Real
MPIHelperMaxMin::getMaxLocationWithinQuantityBoundsInZDirection(
    HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_quantity,
    const int component_idx,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const Real bound_lo,
    const Real bound_hi) const
{
    Real location_z_max_global;
    
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
            << ": MPIHelperMaxMin::getMaxLocationWithinQuantityBoundsInZDirection():\n"
            << "Cannot compute maximum location within quantity bounds in z-direction for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": MPIHelperMaxMin::getMaxLocationWithinQuantityBoundsInZDirection():\n"
            << "Cannot compute maximum location within quantity bounds in z-direction for two-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        Real location_z_max_local = Real(x_lo[2]);
        location_z_max_global     = Real(x_lo[2]);
        
        const Real half = Real(1)/Real(2);
        
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
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity(
                    HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
                        patch->getPatchData(variable_quantity, data_context)));
                
                Real* u = data_quantity->getPointer(component_idx);
                
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
                                
                                const Real z = (Real(relative_idx_lo_2 + k) + half)*Real(dx[2]) + Real(x_lo_patch[2]);
                                
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
            }
        }
        
        /*
         * Reduction to get the global max.
         */
        
        d_mpi.Allreduce(
            &location_z_max_local,
            &location_z_max_global,
            1,
            HAMERS_MPI_REAL,
            MPI_MAX);
    }
    
    return location_z_max_global;
}


/*
 * Compute minimum location within quantity bounds in z-direction.
 */
Real
MPIHelperMaxMin::getMinLocationWithinQuantityBoundsInZDirection(
    HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_quantity,
    const int component_idx,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const Real bound_lo,
    const Real bound_hi) const
{
    Real location_z_min_global;
    
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
            << ": MPIHelperMaxMin::getMinLocationWithinQuantityBoundsInZDirection():\n"
            << "Cannot compute maximum location within quantity bounds in z-direction for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": MPIHelperMaxMin::getMinLocationWithinQuantityBoundsInZDirection():\n"
            << "Cannot compute maximum location within quantity bounds in z-direction for two-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        Real location_z_min_local = Real(x_hi[2]);
        location_z_min_global     = Real(x_hi[2]);
        
        const Real half = Real(1)/Real(2);
        
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
                
                HAMERS_SHARED_PTR<pdat::CellData<Real> > data_quantity(
                    HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
                        patch->getPatchData(variable_quantity, data_context)));
                
                Real* u = data_quantity->getPointer(component_idx);
                
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
                                
                                const Real z = (Real(relative_idx_lo_2 + k) + half)*Real(dx[2]) + Real(x_lo_patch[2]);
                                
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
            }
        }
        
        /*
         * Reduction to get the global min.
         */
        
        d_mpi.Allreduce(
            &location_z_min_local,
            &location_z_min_global,
            1,
            HAMERS_MPI_REAL,
            MPI_MIN);
    }
    
    return location_z_min_global;
}
