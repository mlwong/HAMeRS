#include "util/filters/FilterNone.hpp"

#include<cfloat>

/*
 * Apply filter to the given cell data.
 */
void
FilterNone::applyFilter(
    boost::shared_ptr<pdat::CellData<double> >& filtered_cell_data,
    const boost::shared_ptr<pdat::CellData<double> >& cell_data,
    const int depth_filtered_cell_data,
    const int depth_cell_data,
    const hier::Box& domain)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(filtered_cell_data);
    TBOX_ASSERT(cell_data);
    TBOX_ASSERT(depth_filtered_cell_data < filtered_cell_data->getDepth());
    TBOX_ASSERT(depth_cell_data < cell_data->getDepth());
#endif
    
    // Get the dimensions of box that covers interior of patch plus ghost cells.
    const hier::Box ghost_box_cell_data = cell_data->getGhostBox();
    const hier::IntVector ghostcell_dims_cell_data = ghost_box_cell_data.numberCells();
    
    const hier::Box ghost_box_filtered_cell_data = filtered_cell_data->getGhostBox();
    const hier::IntVector ghostcell_dims_filtered_cell_data = ghost_box_filtered_cell_data.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_cell_data(d_dim);
    hier::IntVector offset_filtered_cell_data(d_dim);
    
    if (domain.empty())
    {
        // Get the number of ghost cells of the cell data and filtered cell data.
        const hier::IntVector num_ghosts_cell_data = cell_data->getGhostCellWidth();
        const hier::IntVector num_ghosts_filtered_cell_data = filtered_cell_data->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = cell_data->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(filtered_cell_data->getBox().isSpatiallyEqual(interior_box));
        
        /*
         * Check potential failures.
         */
        
        if (num_ghosts_cell_data < d_num_filter_ghosts)
        {
            TBOX_ERROR(d_object_name
                << ": FilterNone::applyFilter()\n"
                << "The ghost cell width of cell data is smaller than required."
                << std::endl);
        }
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_cell_data - d_num_filter_ghosts;
        num_ghosts_min = hier::IntVector::min(num_ghosts_filtered_cell_data, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_cell_data = num_ghosts_cell_data;
        offset_filtered_cell_data = num_ghosts_filtered_cell_data;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        hier::Box shrinked_ghost_box_cell_data(ghost_box_cell_data);
        shrinked_ghost_box_cell_data.grow(-d_num_filter_ghosts);
        
        TBOX_ASSERT(shrinked_ghost_box_cell_data.contains(domain));
        TBOX_ASSERT(ghost_box_filtered_cell_data.contains(domain));
#endif        
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_cell_data = domain.lower() - ghost_box_cell_data.lower();
        offset_filtered_cell_data = domain.lower() - ghost_box_filtered_cell_data.lower();
    }
    
    // Get the pointers to the depth components of the given cell data.
    double* f_filtered = filtered_cell_data->getPointer(depth_filtered_cell_data);
    double* f = cell_data->getPointer(depth_cell_data);
    
    if (d_direction == DIRECTION::X_DIRECTION)
    {
        if (d_dim == tbox::Dimension(1))
        {
            /*
             * Get the local lower index, numbers of cells in each dimension and offsets.
             */
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_dim_0 = domain_dims[0];
            
            const int offset_0_cell_data = offset_cell_data[0];
            const int offset_0_filtered_cell_data = offset_filtered_cell_data[0];
            
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx           = i + offset_0_filtered_cell_data;
                const int idx_cell_data = i + offset_0_cell_data;
                
                // Filter in the x-direction.
                f_filtered[idx] = f[idx_cell_data];
            }
        }
        else if (d_dim == tbox::Dimension(2))
        {
            /*
             * Get the local lower indices, numbers of cells in each dimension and offsets.
             */
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_lo_1 = domain_lo[1];
            const int domain_dim_0 = domain_dims[0];
            const int domain_dim_1 = domain_dims[1];
            
            const int offset_0_cell_data = offset_cell_data[0];
            const int offset_1_cell_data = offset_cell_data[1];
            const int ghostcell_dim_0_cell_data = ghostcell_dims_cell_data[0];
            
            const int offset_0_filtered_cell_data = offset_filtered_cell_data[0];
            const int offset_1_filtered_cell_data = offset_filtered_cell_data[1];
            const int ghostcell_dim_0_filtered_cell_data = ghostcell_dims_filtered_cell_data[0];
            
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {        
                    // Compute the linear indices.
                    const int idx = (i + offset_0_filtered_cell_data) +
                        (j + offset_1_filtered_cell_data)*ghostcell_dim_0_filtered_cell_data;
                    
                    const int idx_cell_data = (i + offset_0_cell_data) +
                        (j + offset_1_cell_data)*ghostcell_dim_0_cell_data;
                    
                    // Filter in the x-direction.
                    f_filtered[idx] = f[idx_cell_data];
                }
            }
        }
        else if (d_dim == tbox::Dimension(3))
        {
            /*
             * Get the local lower indices, numbers of cells in each dimension and offsets.
             */
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_lo_1 = domain_lo[1];
            const int domain_lo_2 = domain_lo[2];
            const int domain_dim_0 = domain_dims[0];
            const int domain_dim_1 = domain_dims[1];
            const int domain_dim_2 = domain_dims[2];
            
            const int offset_0_cell_data = offset_cell_data[0];
            const int offset_1_cell_data = offset_cell_data[1];
            const int offset_2_cell_data = offset_cell_data[2];
            const int ghostcell_dim_0_cell_data = ghostcell_dims_cell_data[0];
            const int ghostcell_dim_1_cell_data = ghostcell_dims_cell_data[1];
            
            const int offset_0_filtered_cell_data = offset_filtered_cell_data[0];
            const int offset_1_filtered_cell_data = offset_filtered_cell_data[1];
            const int offset_2_filtered_cell_data = offset_filtered_cell_data[2];
            const int ghostcell_dim_0_filtered_cell_data = ghostcell_dims_filtered_cell_data[0];
            const int ghostcell_dim_1_filtered_cell_data = ghostcell_dims_filtered_cell_data[1];
            
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
                        const int idx = (i + offset_0_filtered_cell_data) +
                            (j + offset_1_filtered_cell_data)*ghostcell_dim_0_filtered_cell_data +
                            (k + offset_2_filtered_cell_data)*ghostcell_dim_0_filtered_cell_data*
                                ghostcell_dim_1_filtered_cell_data;
                        
                        const int idx_cell_data = (i + offset_0_cell_data) +
                            (j + offset_1_cell_data)*ghostcell_dim_0_cell_data +
                            (k + offset_2_cell_data)*ghostcell_dim_0_cell_data*
                                ghostcell_dim_1_cell_data;
                        
                        // Filter in the x-direction.
                        f_filtered[idx] = f[idx_cell_data];
                    }
                }
            }
        }
    }
    else if (d_direction == DIRECTION::Y_DIRECTION)
    {
        if (d_dim == tbox::Dimension(2))
        {
            /*
             * Get the local lower indices, numbers of cells in each dimension and offsets.
             */
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_lo_1 = domain_lo[1];
            const int domain_dim_0 = domain_dims[0];
            const int domain_dim_1 = domain_dims[1];
            
            const int offset_0_cell_data = offset_cell_data[0];
            const int offset_1_cell_data = offset_cell_data[1];
            const int ghostcell_dim_0_cell_data = ghostcell_dims_cell_data[0];
            
            const int offset_0_filtered_cell_data = offset_filtered_cell_data[0];
            const int offset_1_filtered_cell_data = offset_filtered_cell_data[1];
            const int ghostcell_dim_0_filtered_cell_data = ghostcell_dims_filtered_cell_data[0];
            
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {        
                    // Compute the linear indices.
                    const int idx = (i + offset_0_filtered_cell_data) +
                        (j + offset_1_filtered_cell_data)*ghostcell_dim_0_filtered_cell_data;
                    
                    const int idx_cell_data = (i + offset_0_cell_data) +
                        (j + offset_1_cell_data)*ghostcell_dim_0_cell_data;
                    
                    // Filter in the y-direction.
                    f_filtered[idx] = f[idx_cell_data];
                }
            }
        }
        else if (d_dim == tbox::Dimension(3))
        {
            /*
             * Get the local lower indices, numbers of cells in each dimension and offsets.
             */
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_lo_1 = domain_lo[1];
            const int domain_lo_2 = domain_lo[2];
            const int domain_dim_0 = domain_dims[0];
            const int domain_dim_1 = domain_dims[1];
            const int domain_dim_2 = domain_dims[2];
            
            const int offset_0_cell_data = offset_cell_data[0];
            const int offset_1_cell_data = offset_cell_data[1];
            const int offset_2_cell_data = offset_cell_data[2];
            const int ghostcell_dim_0_cell_data = ghostcell_dims_cell_data[0];
            const int ghostcell_dim_1_cell_data = ghostcell_dims_cell_data[1];
            
            const int offset_0_filtered_cell_data = offset_filtered_cell_data[0];
            const int offset_1_filtered_cell_data = offset_filtered_cell_data[1];
            const int offset_2_filtered_cell_data = offset_filtered_cell_data[2];
            const int ghostcell_dim_0_filtered_cell_data = ghostcell_dims_filtered_cell_data[0];
            const int ghostcell_dim_1_filtered_cell_data = ghostcell_dims_filtered_cell_data[1];
            
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
                        const int idx = (i + offset_0_filtered_cell_data) +
                            (j + offset_1_filtered_cell_data)*ghostcell_dim_0_filtered_cell_data +
                            (k + offset_2_filtered_cell_data)*ghostcell_dim_0_filtered_cell_data*
                                ghostcell_dim_1_filtered_cell_data;
                        
                        const int idx_cell_data = (i + offset_0_cell_data) +
                            (j + offset_1_cell_data)*ghostcell_dim_0_cell_data +
                            (k + offset_2_cell_data)*ghostcell_dim_0_cell_data*
                                ghostcell_dim_1_cell_data;
                        
                        // Filter in the y-direction.
                        f_filtered[idx] = f[idx_cell_data];
                    }
                }
            }
        }
    }
    else if (d_direction == DIRECTION::Z_DIRECTION)
    {
        /*
         * Get the local lower indices, numbers of cells in each dimension and offsets.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_lo_2 = domain_lo[2];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        const int domain_dim_2 = domain_dims[2];
        
        const int offset_0_cell_data = offset_cell_data[0];
        const int offset_1_cell_data = offset_cell_data[1];
        const int offset_2_cell_data = offset_cell_data[2];
        const int ghostcell_dim_0_cell_data = ghostcell_dims_cell_data[0];
        const int ghostcell_dim_1_cell_data = ghostcell_dims_cell_data[1];
        
        const int offset_0_filtered_cell_data = offset_filtered_cell_data[0];
        const int offset_1_filtered_cell_data = offset_filtered_cell_data[1];
        const int offset_2_filtered_cell_data = offset_filtered_cell_data[2];
        const int ghostcell_dim_0_filtered_cell_data = ghostcell_dims_filtered_cell_data[0];
        const int ghostcell_dim_1_filtered_cell_data = ghostcell_dims_filtered_cell_data[1];
        
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
                    const int idx = (i + offset_0_filtered_cell_data) +
                        (j + offset_1_filtered_cell_data)*ghostcell_dim_0_filtered_cell_data +
                        (k + offset_2_filtered_cell_data)*ghostcell_dim_0_filtered_cell_data*
                            ghostcell_dim_1_filtered_cell_data;
                    
                    const int idx_cell_data = (i + offset_0_cell_data) +
                        (j + offset_1_cell_data)*ghostcell_dim_0_cell_data +
                        (k + offset_2_cell_data)*ghostcell_dim_0_cell_data*
                            ghostcell_dim_1_cell_data;
                    
                    // Filter in the z-direction.
                    f_filtered[idx] = f[idx_cell_data];
                }
            }
        }
    }
}
