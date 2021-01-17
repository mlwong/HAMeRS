#include "util/MPI_helpers/MPIHelperCentroid.hpp"

#include "extn/patch_hierarchies/ExtendedFlattenedHierarchy.hpp"

#include "SAMRAI/geom/CartesianPatchGeometry.h"

/*
 * Compute centroid in x-direction.
 */
double
MPIHelperCentroid::getCentroidInXDirection(
    HAMERS_SHARED_PTR<pdat::CellVariable<double> >& variable_quantity,
    const int component_idx,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const
{
    double x_c;
    
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
        double num_local  = double(0);
        double num_global = double(0);
        
        double den_local  = double(0);
        double den_global = double(0);
        
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
                
                const double dx_0 = dx[0];
                
                const double* const x_lo_patch = patch_geom->getXLower();
                
                const double x_lo_patch_0 = x_lo_patch[0];
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_quantity(
                    HAMERS_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch->getPatchData(variable_quantity, data_context)));
                
                double* u = data_quantity->getPointer(component_idx);
                
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
                
                double num_to_add = double(0);
                double den_to_add = double(0);
                
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
                        
                        const double x = (relative_idx_lo_0 + i + 0.5)*dx_0 + x_lo_patch_0;
                        
                        num_to_add += x*u[idx]/((double) n_overlapped);
                        den_to_add += u[idx]/((double) n_overlapped);
                    }
                }
                
                num_to_add = num_to_add*dx[0];
                den_to_add = den_to_add*dx[0];
                
                num_local += num_to_add;
                den_local += den_to_add;
            }
        }
        
        /*
         * Reduction to get the global integrals.
         */
        
        d_mpi.Allreduce(
            &num_local,
            &num_global,
            1,
            MPI_DOUBLE,
            MPI_SUM);
        
        d_mpi.Allreduce(
            &den_local,
            &den_global,
            1,
            MPI_DOUBLE,
            MPI_SUM);
        
        x_c = num_global/den_global;
    }
    else if (d_dim == tbox::Dimension(2))
    {
        double num_local  = double(0);
        double num_global = double(0);
        
        double den_local  = double(0);
        double den_global = double(0);
        
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
                
                const double dx_0 = dx[0];
                
                const double* const x_lo_patch = patch_geom->getXLower();
                
                const double x_lo_patch_0 = x_lo_patch[0];
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_quantity(
                    HAMERS_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch->getPatchData(variable_quantity, data_context)));
                
                double* u = data_quantity->getPointer(component_idx);
                
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
                
                double num_to_add = double(0);
                double den_to_add = double(0);
                
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
                            
                            const double x = (relative_idx_lo_0 + i + 0.5)*dx_0 + x_lo_patch_0;
                            
                            num_to_add += x*u[idx]/((double) n_overlapped);
                            den_to_add += u[idx]/((double) n_overlapped);
                        }
                    }
                }
                
                num_to_add = num_to_add*dx[0]*dx[1];
                den_to_add = den_to_add*dx[0]*dx[1];
                
                num_local += num_to_add;
                den_local += den_to_add;
            }
        }
        
        /*
         * Reduction to get the global integrals.
         */
        
        d_mpi.Allreduce(
            &num_local,
            &num_global,
            1,
            MPI_DOUBLE,
            MPI_SUM);
        
        d_mpi.Allreduce(
            &den_local,
            &den_global,
            1,
            MPI_DOUBLE,
            MPI_SUM);
        
        x_c = num_global/den_global;
    }
    else if (d_dim == tbox::Dimension(3))
    {
        double num_local  = double(0);
        double num_global = double(0);
        
        double den_local  = double(0);
        double den_global = double(0);
        
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
                
                const double dx_0 = dx[0];
                
                const double* const x_lo_patch = patch_geom->getXLower();
                
                const double x_lo_patch_0 = x_lo_patch[0];
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_quantity(
                    HAMERS_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch->getPatchData(variable_quantity, data_context)));
                
                double* u = data_quantity->getPointer(component_idx);
                
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
                
                double num_to_add = double(0);
                double den_to_add = double(0);
                
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
                                
                                const double x = (relative_idx_lo_0 + i + 0.5)*dx_0 + x_lo_patch_0;
                                
                                num_to_add += x*u[idx]/((double) n_overlapped);
                                den_to_add += u[idx]/((double) n_overlapped);
                            }
                        }
                    }
                }
                
                num_to_add = num_to_add*dx[0]*dx[1]*dx[2];
                den_to_add = den_to_add*dx[0]*dx[1]*dx[2];
                
                num_local += num_to_add;
                den_local += den_to_add;
            }
        }
        
        /*
         * Reduction to get the global integrals.
         */
        
        d_mpi.Allreduce(
            &num_local,
            &num_global,
            1,
            MPI_DOUBLE,
            MPI_SUM);
        
        d_mpi.Allreduce(
            &den_local,
            &den_global,
            1,
            MPI_DOUBLE,
            MPI_SUM);
        
        x_c = num_global/den_global;
    }
    
    return x_c;
}

