#include "util/MPI_helpers/MPIHelperGrid.hpp"

#include "extn/patch_hierarchies/ExtendedFlattenedHierarchy.hpp"

#include "SAMRAI/geom/CartesianPatchGeometry.h"

/*
 * Compute number of cells.
 */
Real
MPIHelperGrid::getNumberOfCells() const
{
    Real num_cells_global;
    
    const int num_levels = d_patch_hierarchy->getNumberOfLevels();
    
    if (d_dim == tbox::Dimension(1))
    {
        Real num_cells_local = Real(0);
        num_cells_global       = Real(0);
        
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
                 * Add the number of cells in current patch.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::IntVector interior_dims = patch_box.numberCells();
                
                num_cells_local += Real(interior_dims[0]);
            }
        }
        
        /*
         * Reduction to get the global number of cells.
         */
        
        d_mpi.Allreduce(
            &num_cells_local,
            &num_cells_global,
            1,
            HAMERS_MPI_REAL,
            MPI_SUM);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        Real num_cells_local = Real(0);
        num_cells_global       = Real(0);
        
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
                 * Add the number of cells in current patch.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::IntVector interior_dims = patch_box.numberCells();
                
                num_cells_local += Real(interior_dims[0])*Real(interior_dims[1]);
            }
        }
        
        /*
         * Reduction to get the global integrals.
         */
        
        d_mpi.Allreduce(
            &num_cells_local,
            &num_cells_global,
            1,
            HAMERS_MPI_REAL,
            MPI_SUM);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        Real num_cells_local = Real(0);
        num_cells_global       = Real(0);
        
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
                 * Add the number of cells in current patch.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::IntVector interior_dims = patch_box.numberCells();
                
                num_cells_local += Real(interior_dims[0])*Real(interior_dims[1])*Real(interior_dims[2]);
            }
        }
        
        /*
         * Reduction to get the global number of cells.
         */
        
        d_mpi.Allreduce(
            &num_cells_local,
            &num_cells_global,
            1,
            HAMERS_MPI_REAL,
            MPI_SUM);
    }
    
    return num_cells_global;
}


/*
 * Compute weighted number of cells.
 */
Real
MPIHelperGrid::getWeightedNumberOfCells() const
{
    Real weighted_num_cells_global;
    
    const int num_levels = d_patch_hierarchy->getNumberOfLevels();
    
    if (d_dim == tbox::Dimension(1))
    {
        Real weighted_num_cells_local = Real(0);
        weighted_num_cells_global       = Real(0);
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the refinement ratio from the current level to the coarsest level.
             */
            
            hier::IntVector ratioCurrentLevelToCoarsestLevel =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioCurrentLevelToCoarsestLevel *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            if (li == 0)
            {
                ratioCurrentLevelToCoarsestLevel = hier::IntVector::getOne(d_dim);
            }
            
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
                 * Add the number of cells in current patch.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::IntVector interior_dims = patch_box.numberCells();
                
                weighted_num_cells_local += Real(interior_dims[0])*Real(ratioCurrentLevelToCoarsestLevel[0]);
            }
        }
        
        /*
         * Reduction to get the global number of cells.
         */
        
        d_mpi.Allreduce(
            &weighted_num_cells_local,
            &weighted_num_cells_global,
            1,
            HAMERS_MPI_REAL,
            MPI_SUM);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        Real weighted_num_cells_local = Real(0);
        weighted_num_cells_global       = Real(0);
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the refinement ratio from the current level to the coarsest level.
             */
            
            hier::IntVector ratioCurrentLevelToCoarsestLevel =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioCurrentLevelToCoarsestLevel *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            if (li == 0)
            {
                ratioCurrentLevelToCoarsestLevel = hier::IntVector::getOne(d_dim);
            }
            
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
                 * Add the number of cells in current patch.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::IntVector interior_dims = patch_box.numberCells();
                
                weighted_num_cells_local += Real(interior_dims[0])*Real(interior_dims[1])*
                    Real(ratioCurrentLevelToCoarsestLevel[0]);
            }
        }
        
        /*
         * Reduction to get the global integrals.
         */
        
        d_mpi.Allreduce(
            &weighted_num_cells_local,
            &weighted_num_cells_global,
            1,
            HAMERS_MPI_REAL,
            MPI_SUM);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        Real weighted_num_cells_local = Real(0);
        weighted_num_cells_global       = Real(0);
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the refinement ratio from the current level to the coarsest level.
             */
            
            hier::IntVector ratioCurrentLevelToCoarsestLevel =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioCurrentLevelToCoarsestLevel *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            if (li == 0)
            {
                ratioCurrentLevelToCoarsestLevel = hier::IntVector::getOne(d_dim);
            }
            
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
                 * Add the number of cells in current patch.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::IntVector interior_dims = patch_box.numberCells();
                
                weighted_num_cells_local += Real(interior_dims[0])*Real(interior_dims[1])*
                    Real(ratioCurrentLevelToCoarsestLevel[0]);
            }
        }
        
        /*
         * Reduction to get the global number of cells.
         */
        
        d_mpi.Allreduce(
            &weighted_num_cells_local,
            &weighted_num_cells_global,
            1,
            HAMERS_MPI_REAL,
            MPI_SUM);
    }
    
    return weighted_num_cells_global;
}


/*
 * Compute averaged grid level number with only x-direction as inhomogeneous direction.
 */
std::vector<Real>
MPIHelperGrid::getAveragedGridLevelNumberWithInhomogeneousXDirection() const
{
    std::vector<Real> averaged_grid_level_num;
    
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
        
        Real* ln_avg_local = (Real*)std::malloc(finest_level_dim_0*sizeof(Real));
        
        averaged_grid_level_num.resize(finest_level_dim_0);
        Real* ln_avg_global = averaged_grid_level_num.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            ln_avg_local[i]  = Real(0);
            ln_avg_global[i] = Real(0);
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
                 * Get the patch box.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    
                    const int idx_lo_0 = index_lo[0];
                    
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
                         * Compute the linear index and the grid level to add.
                         */
                        
                        const Real value_to_add = Real(li)/((Real) n_overlapped);
                        
                        for (int ii = 0; ii < ratio_to_finest_level_0; ii++)
                        {
                            const int idx_fine = (idx_lo_0 + i)*ratio_to_finest_level_0 + ii;
                            
                            ln_avg_local[idx_fine] += value_to_add;
                        }
                    }
                }
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        d_mpi.Allreduce(
            ln_avg_local,
            ln_avg_global,
            finest_level_dim_0,
            HAMERS_MPI_REAL,
            MPI_SUM);
        
        std::free(ln_avg_local);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_0 = d_finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_y = Real(x_hi[1] - x_lo[1]);
        
        Real* ln_avg_local = (Real*)std::malloc(finest_level_dim_0*sizeof(Real));
        
        averaged_grid_level_num.resize(finest_level_dim_0);
        Real* ln_avg_global = averaged_grid_level_num.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            ln_avg_local[i]  = Real(0);
            ln_avg_global[i] = Real(0);
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
                 * Get the patch box and grid spacings.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
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
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    
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
                             * Compute the linear index and the grid level to add.
                             */
                            
                            const Real value_to_add = Real(li)*weight/((Real) n_overlapped);
                            
                            for (int ii = 0; ii < ratio_to_finest_level_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratio_to_finest_level_0 + ii;
                                
                                ln_avg_local[idx_fine] += value_to_add;
                            }
                        }
                    }
                }
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        d_mpi.Allreduce(
            ln_avg_local,
            ln_avg_global,
            finest_level_dim_0,
            HAMERS_MPI_REAL,
            MPI_SUM);
        
        std::free(ln_avg_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_0 = d_finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_y = Real(x_hi[1] - x_lo[1]);
        const Real L_z = Real(x_hi[2] - x_lo[2]);
        
        Real* ln_avg_local = (Real*)std::malloc(finest_level_dim_0*sizeof(Real));
        
        averaged_grid_level_num.resize(finest_level_dim_0);
        Real* ln_avg_global = averaged_grid_level_num.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            ln_avg_local[i]  = Real(0);
            ln_avg_global[i] = Real(0);
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
                 * Get the patch box and grid spacings.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
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
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int idx_lo_2 = index_lo[2];
                    
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
                                 * Compute the linear index and the grid level to add.
                                 */
                                
                                const Real value_to_add = Real(li)*weight/((Real) n_overlapped);
                                
                                for (int ii = 0; ii < ratio_to_finest_level_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratio_to_finest_level_0 + ii;
                                    
                                    ln_avg_local[idx_fine] += value_to_add;
                                }
                            }
                        }
                    }
                }
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        d_mpi.Allreduce(
            ln_avg_local,
            ln_avg_global,
            finest_level_dim_0,
            HAMERS_MPI_REAL,
            MPI_SUM);
        
        std::free(ln_avg_local);
    }
    
    return averaged_grid_level_num;
}


/*
 * Compute averaged grid level number with only y-direction as inhomogeneous direction.
 */
std::vector<Real>
MPIHelperGrid::getAveragedGridLevelNumberWithInhomogeneousYDirection() const
{
    std::vector<Real> averaged_grid_level_num;
    
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
            << ": MPIHelperGrid::getAveragedGridLevelNumberWithInhomogeneousYDirection():\n"
            << "Cannot compute averaged grid level number with only y-direction as inhomogeneous direction for"
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
        
        Real* ln_avg_local = (Real*)std::malloc(finest_level_dim_1*sizeof(Real));
        
        averaged_grid_level_num.resize(finest_level_dim_1);
        Real* ln_avg_global = averaged_grid_level_num.data();
        
        for (int j = 0; j < finest_level_dim_1; j++)
        {
            ln_avg_local[j]  = Real(0);
            ln_avg_global[j] = Real(0);
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
                 * Get the patch box and grid spacings.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
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
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    
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
                             * Compute the linear index and the grid level to add.
                             */
                            
                            const Real value_to_add = Real(li)*weight/((Real) n_overlapped);
                            
                            for (int jj = 0; jj < ratio_to_finest_level_1; jj++)
                            {
                                const int idx_fine = (idx_lo_1 + j)*ratio_to_finest_level_1 + jj;
                                
                                ln_avg_local[idx_fine] += value_to_add;
                            }
                        }
                    }
                }
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        d_mpi.Allreduce(
            ln_avg_local,
            ln_avg_global,
            finest_level_dim_1,
            HAMERS_MPI_REAL,
            MPI_SUM);
        
        std::free(ln_avg_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_1 = d_finest_level_dims[1];
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_x = Real(x_hi[0] - x_lo[0]);
        const Real L_z = Real(x_hi[2] - x_lo[2]);
        
        Real* ln_avg_local = (Real*)std::malloc(finest_level_dim_1*sizeof(Real));
        
        averaged_grid_level_num.resize(finest_level_dim_1);
        Real* ln_avg_global = averaged_grid_level_num.data();
        
        for (int j = 0; j < finest_level_dim_1; j++)
        {
            ln_avg_local[j]  = Real(0);
            ln_avg_global[j] = Real(0);
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
                 * Get the patch box and grid spacings.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
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
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int idx_lo_2 = index_lo[2];
                    
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
                                 * Compute the linear index and the grid level to add.
                                 */
                                
                                const Real value_to_add = Real(li)*weight/((Real) n_overlapped);
                                
                                for (int jj = 0; jj < ratio_to_finest_level_1; jj++)
                                {
                                    const int idx_fine = (idx_lo_1 + j)*ratio_to_finest_level_1 + jj;
                                    
                                    ln_avg_local[idx_fine] += value_to_add;
                                }
                            }
                        }
                    }
                }
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        d_mpi.Allreduce(
            ln_avg_local,
            ln_avg_global,
            finest_level_dim_1,
            HAMERS_MPI_REAL,
            MPI_SUM);
        
        std::free(ln_avg_local);
    }
    
    return averaged_grid_level_num;
}


/*
 * Compute averaged grid level number with only z-direction as inhomogeneous direction.
 */
std::vector<Real>
MPIHelperGrid::getAveragedGridLevelNumberWithInhomogeneousZDirection() const
{
    std::vector<Real> averaged_grid_level_num;
    
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
            << ": MPIHelperAverage::getAveragedGridLevelNumberWithInhomogeneousZDirection():\n"
            << "Cannot compute averaged grid level number with only z-direction as inhomogeneous direction for"
            << " one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": MPIHelperAverage::getAveragedGridLevelNumberWithInhomogeneousZDirection():\n"
            << "Cannot compute averaged grid level number with only z-direction as inhomogeneous direction for"
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
        
        Real* ln_avg_local = (Real*)std::malloc(finest_level_dim_2*sizeof(Real));
        
        averaged_grid_level_num.resize(finest_level_dim_2);
        Real* ln_avg_global = averaged_grid_level_num.data();
        
        for (int k = 0; k < finest_level_dim_2; k++)
        {
            ln_avg_local[k]  = Real(0);
            ln_avg_global[k] = Real(0);
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
                 * Get the patch box and grid spacings.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
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
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int idx_lo_2 = index_lo[2];
                    
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
                                 * Compute the linear index and the grid level to add.
                                 */
                                
                                const Real value_to_add = Real(li)*weight/((Real) n_overlapped);
                                
                                for (int kk = 0; kk < ratio_to_finest_level_2; kk++)
                                {
                                    const int idx_fine = (idx_lo_2 + k)*ratio_to_finest_level_2 + kk;
                                    
                                    ln_avg_local[idx_fine] += value_to_add;
                                }
                            }
                        }
                    }
                }
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        d_mpi.Allreduce(
            ln_avg_local,
            ln_avg_global,
            finest_level_dim_2,
            HAMERS_MPI_REAL,
            MPI_SUM);
        
        std::free(ln_avg_local);
    }
    
    return averaged_grid_level_num;
}
