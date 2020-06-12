#include "util/filters/FilterTruncatedGaussian.hpp"

#include<cfloat>

FilterTruncatedGaussian::FilterTruncatedGaussian(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const DIRECTION::TYPE& direction):
        Filter(
            object_name,
            dim,
            direction),
        a_G(double(3565)/double( 10368)),
        b_G(double(3091)/double( 12960)),
        c_G(double(1997)/double( 25920)),
        d_G(double( 149)/double( 12960)),
        e_G(double( 107)/double(103680))
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    const double sum_coeffs = a_G + double(2)*(b_G + c_G + d_G + e_G);
    TBOX_ASSERT(fabs(sum_coeffs - double(1)) < DBL_EPSILON);
#endif
    
    switch (direction)
    {
        case (DIRECTION::X_DIRECTION):
        {
            d_num_filter_ghosts[0] = 4;
            
            break;
        }
        case (DIRECTION::Y_DIRECTION):
        {
            d_num_filter_ghosts[1] = 4;
            
            break;
        }
        case (DIRECTION::Z_DIRECTION):
        {
            d_num_filter_ghosts[2] = 4;
            
            break;
        }
    }
}

/*
 * Apply filter to the given cell data.
 */
void
FilterTruncatedGaussian::applyFilter(
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
                << ": FilterTruncatedGaussian::applyFilter()\n"
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
                const int idx = i + offset_0_filtered_cell_data;
                
                const int idx_cell_data        = i     + offset_0_cell_data;
                const int idx_cell_data_x_LLLL = i - 4 + offset_0_cell_data;
                const int idx_cell_data_x_LLL  = i - 3 + offset_0_cell_data;
                const int idx_cell_data_x_LL   = i - 2 + offset_0_cell_data;
                const int idx_cell_data_x_L    = i - 1 + offset_0_cell_data;
                const int idx_cell_data_x_R    = i + 1 + offset_0_cell_data;
                const int idx_cell_data_x_RR   = i + 2 + offset_0_cell_data;
                const int idx_cell_data_x_RRR  = i + 3 + offset_0_cell_data;
                const int idx_cell_data_x_RRRR = i + 4 + offset_0_cell_data;
                
                // Filter in the x-direction.
                f_filtered[idx] = a_G*f[idx_cell_data] + 
                    b_G*(f[idx_cell_data_x_L] + f[idx_cell_data_x_R]) +
                    c_G*(f[idx_cell_data_x_LL] + f[idx_cell_data_x_RR]) +
                    d_G*(f[idx_cell_data_x_LLL] + f[idx_cell_data_x_RRR]) +
                    e_G*(f[idx_cell_data_x_LLLL] + f[idx_cell_data_x_RRRR]);
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
                    
                    const int idx_cell_data_x_LLLL = (i - 4 + offset_0_cell_data) +
                        (j + offset_1_cell_data)*ghostcell_dim_0_cell_data;
                    
                    const int idx_cell_data_x_LLL = (i - 3 + offset_0_cell_data) +
                        (j + offset_1_cell_data)*ghostcell_dim_0_cell_data;
                    
                    const int idx_cell_data_x_LL = (i - 2 + offset_0_cell_data) +
                        (j + offset_1_cell_data)*ghostcell_dim_0_cell_data;
                    
                    const int idx_cell_data_x_L = (i - 1 + offset_0_cell_data) +
                        (j + offset_1_cell_data)*ghostcell_dim_0_cell_data;
                    
                    const int idx_cell_data_x_R = (i + 1 + offset_0_cell_data) +
                        (j + offset_1_cell_data)*ghostcell_dim_0_cell_data;
                    
                    const int idx_cell_data_x_RR = (i + 2 + offset_0_cell_data) +
                        (j + offset_1_cell_data)*ghostcell_dim_0_cell_data;
                    
                    const int idx_cell_data_x_RRR = (i + 3 + offset_0_cell_data) +
                        (j + offset_1_cell_data)*ghostcell_dim_0_cell_data;
                    
                    const int idx_cell_data_x_RRRR = (i + 4 + offset_0_cell_data) +
                        (j + offset_1_cell_data)*ghostcell_dim_0_cell_data;
                    
                    // Filter in the x-direction.
                    f_filtered[idx] = a_G*f[idx_cell_data] + 
                        b_G*(f[idx_cell_data_x_L]    + f[idx_cell_data_x_R]) +
                        c_G*(f[idx_cell_data_x_LL]   + f[idx_cell_data_x_RR]) +
                        d_G*(f[idx_cell_data_x_LLL]  + f[idx_cell_data_x_RRR]) +
                        e_G*(f[idx_cell_data_x_LLLL] + f[idx_cell_data_x_RRRR]);
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
                        
                        const int idx_cell_data_x_LLLL = (i - 4 + offset_0_cell_data) +
                            (j + offset_1_cell_data)*ghostcell_dim_0_cell_data +
                            (k + offset_2_cell_data)*ghostcell_dim_0_cell_data*
                                ghostcell_dim_1_cell_data;
                        
                        const int idx_cell_data_x_LLL = (i - 3 + offset_0_cell_data) +
                            (j + offset_1_cell_data)*ghostcell_dim_0_cell_data +
                            (k + offset_2_cell_data)*ghostcell_dim_0_cell_data*
                                ghostcell_dim_1_cell_data;
                        
                        const int idx_cell_data_x_LL = (i - 2 + offset_0_cell_data) +
                            (j + offset_1_cell_data)*ghostcell_dim_0_cell_data +
                            (k + offset_2_cell_data)*ghostcell_dim_0_cell_data*
                                ghostcell_dim_1_cell_data;
                        
                        const int idx_cell_data_x_L = (i - 1 + offset_0_cell_data) +
                            (j + offset_1_cell_data)*ghostcell_dim_0_cell_data +
                            (k + offset_2_cell_data)*ghostcell_dim_0_cell_data*
                                ghostcell_dim_1_cell_data;
                        
                        const int idx_cell_data_x_R = (i + 1 + offset_0_cell_data) +
                            (j + offset_1_cell_data)*ghostcell_dim_0_cell_data +
                            (k + offset_2_cell_data)*ghostcell_dim_0_cell_data*
                                ghostcell_dim_1_cell_data;
                        
                        const int idx_cell_data_x_RR = (i + 2 + offset_0_cell_data) +
                            (j + offset_1_cell_data)*ghostcell_dim_0_cell_data +
                            (k + offset_2_cell_data)*ghostcell_dim_0_cell_data*
                                ghostcell_dim_1_cell_data;
                        
                        const int idx_cell_data_x_RRR = (i + 3 + offset_0_cell_data) +
                            (j + offset_1_cell_data)*ghostcell_dim_0_cell_data +
                            (k + offset_2_cell_data)*ghostcell_dim_0_cell_data*
                                ghostcell_dim_1_cell_data;
                        
                        const int idx_cell_data_x_RRRR = (i + 4 + offset_0_cell_data) +
                            (j + offset_1_cell_data)*ghostcell_dim_0_cell_data +
                            (k + offset_2_cell_data)*ghostcell_dim_0_cell_data*
                                ghostcell_dim_1_cell_data;
                        
                        // Filter in the x-direction.
                        f_filtered[idx] = a_G*f[idx_cell_data] + 
                            b_G*(f[idx_cell_data_x_L]    + f[idx_cell_data_x_R]) +
                            c_G*(f[idx_cell_data_x_LL]   + f[idx_cell_data_x_RR]) +
                            d_G*(f[idx_cell_data_x_LLL]  + f[idx_cell_data_x_RRR]) +
                            e_G*(f[idx_cell_data_x_LLLL] + f[idx_cell_data_x_RRRR]);
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
                    
                    const int idx_cell_data_y_BBBB = (i + offset_0_cell_data) +
                        (j - 4 + offset_1_cell_data)*ghostcell_dim_0_cell_data;
                    
                    const int idx_cell_data_y_BBB = (i + offset_0_cell_data) +
                        (j - 3 + offset_1_cell_data)*ghostcell_dim_0_cell_data;
                    
                    const int idx_cell_data_y_BB = (i + offset_0_cell_data) +
                        (j - 2 + offset_1_cell_data)*ghostcell_dim_0_cell_data;
                    
                    const int idx_cell_data_y_B = (i + offset_0_cell_data) +
                        (j - 1 + offset_1_cell_data)*ghostcell_dim_0_cell_data;
                    
                    const int idx_cell_data_y_T = (i + offset_0_cell_data) +
                        (j + 1 + offset_1_cell_data)*ghostcell_dim_0_cell_data;
                    
                    const int idx_cell_data_y_TT = (i + offset_0_cell_data) +
                        (j + 2 + offset_1_cell_data)*ghostcell_dim_0_cell_data;
                    
                    const int idx_cell_data_y_TTT = (i + offset_0_cell_data) +
                        (j + 3 + offset_1_cell_data)*ghostcell_dim_0_cell_data;
                    
                    const int idx_cell_data_y_TTTT = (i + offset_0_cell_data) +
                        (j + 4 + offset_1_cell_data)*ghostcell_dim_0_cell_data;
                    
                    // Filter in the y-direction.
                    f_filtered[idx] = a_G*f[idx_cell_data] + 
                        b_G*(f[idx_cell_data_y_B]    + f[idx_cell_data_y_T]) +
                        c_G*(f[idx_cell_data_y_BB]   + f[idx_cell_data_y_TT]) +
                        d_G*(f[idx_cell_data_y_BBB]  + f[idx_cell_data_y_TTT]) +
                        e_G*(f[idx_cell_data_y_BBBB] + f[idx_cell_data_y_TTTT]);
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
                        
                        const int idx_cell_data_y_BBBB = (i + offset_0_cell_data) +
                            (j - 4 + offset_1_cell_data)*ghostcell_dim_0_cell_data +
                            (k + offset_2_cell_data)*ghostcell_dim_0_cell_data*
                                ghostcell_dim_1_cell_data;
                        
                        const int idx_cell_data_y_BBB = (i + offset_0_cell_data) +
                            (j - 3 + offset_1_cell_data)*ghostcell_dim_0_cell_data +
                            (k + offset_2_cell_data)*ghostcell_dim_0_cell_data*
                                ghostcell_dim_1_cell_data;
                        
                        const int idx_cell_data_y_BB = (i + offset_0_cell_data) +
                            (j - 2 + offset_1_cell_data)*ghostcell_dim_0_cell_data +
                            (k + offset_2_cell_data)*ghostcell_dim_0_cell_data*
                                ghostcell_dim_1_cell_data;
                        
                        const int idx_cell_data_y_B = (i + offset_0_cell_data) +
                            (j - 1 + offset_1_cell_data)*ghostcell_dim_0_cell_data +
                            (k + offset_2_cell_data)*ghostcell_dim_0_cell_data*
                                ghostcell_dim_1_cell_data;
                        
                        const int idx_cell_data_y_T = (i + offset_0_cell_data) +
                            (j + 1 + offset_1_cell_data)*ghostcell_dim_0_cell_data +
                            (k + offset_2_cell_data)*ghostcell_dim_0_cell_data*
                                ghostcell_dim_1_cell_data;
                        
                        const int idx_cell_data_y_TT = (i + offset_0_cell_data) +
                            (j + 2 + offset_1_cell_data)*ghostcell_dim_0_cell_data +
                            (k + offset_2_cell_data)*ghostcell_dim_0_cell_data*
                                ghostcell_dim_1_cell_data;
                        
                        const int idx_cell_data_y_TTT = (i + offset_0_cell_data) +
                            (j + 3 + offset_1_cell_data)*ghostcell_dim_0_cell_data +
                            (k + offset_2_cell_data)*ghostcell_dim_0_cell_data*
                                ghostcell_dim_1_cell_data;
                        
                        const int idx_cell_data_y_TTTT = (i + offset_0_cell_data) +
                            (j + 4 + offset_1_cell_data)*ghostcell_dim_0_cell_data +
                            (k + offset_2_cell_data)*ghostcell_dim_0_cell_data*
                                ghostcell_dim_1_cell_data;
                        
                        // Filter in the y-direction.
                        f_filtered[idx] = a_G*f[idx_cell_data] + 
                            b_G*(f[idx_cell_data_y_B]    + f[idx_cell_data_y_T]) +
                            c_G*(f[idx_cell_data_y_BB]   + f[idx_cell_data_y_TT]) +
                            d_G*(f[idx_cell_data_y_BBB]  + f[idx_cell_data_y_TTT]) +
                            e_G*(f[idx_cell_data_y_BBBB] + f[idx_cell_data_y_TTTT]);
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
                    
                    const int idx_cell_data_z_BBBB = (i + offset_0_cell_data) +
                        (j + offset_1_cell_data)*ghostcell_dim_0_cell_data +
                        (k - 4 + offset_2_cell_data)*ghostcell_dim_0_cell_data*
                            ghostcell_dim_1_cell_data;
                    
                    const int idx_cell_data_z_BBB = (i + offset_0_cell_data) +
                        (j + offset_1_cell_data)*ghostcell_dim_0_cell_data +
                        (k - 3 + offset_2_cell_data)*ghostcell_dim_0_cell_data*
                            ghostcell_dim_1_cell_data;
                    
                    const int idx_cell_data_z_BB = (i + offset_0_cell_data) +
                        (j + offset_1_cell_data)*ghostcell_dim_0_cell_data +
                        (k - 2 + offset_2_cell_data)*ghostcell_dim_0_cell_data*
                            ghostcell_dim_1_cell_data;
                    
                    const int idx_cell_data_z_B = (i + offset_0_cell_data) +
                        (j + offset_1_cell_data)*ghostcell_dim_0_cell_data +
                        (k - 1 + offset_2_cell_data)*ghostcell_dim_0_cell_data*
                            ghostcell_dim_1_cell_data;
                    
                    const int idx_cell_data_z_F = (i + offset_0_cell_data) +
                        (j + offset_1_cell_data)*ghostcell_dim_0_cell_data +
                        (k + 1 + offset_2_cell_data)*ghostcell_dim_0_cell_data*
                            ghostcell_dim_1_cell_data;
                    
                    const int idx_cell_data_z_FF = (i + offset_0_cell_data) +
                        (j + offset_1_cell_data)*ghostcell_dim_0_cell_data +
                        (k + 2 + offset_2_cell_data)*ghostcell_dim_0_cell_data*
                            ghostcell_dim_1_cell_data;
                    
                    const int idx_cell_data_z_FFF = (i + offset_0_cell_data) +
                        (j + offset_1_cell_data)*ghostcell_dim_0_cell_data +
                        (k + 3 + offset_2_cell_data)*ghostcell_dim_0_cell_data*
                            ghostcell_dim_1_cell_data;
                    
                    const int idx_cell_data_z_FFFF = (i + offset_0_cell_data) +
                        (j + offset_1_cell_data)*ghostcell_dim_0_cell_data +
                        (k + 4 + offset_2_cell_data)*ghostcell_dim_0_cell_data*
                            ghostcell_dim_1_cell_data;
                    
                    // Filter in the z-direction.
                    f_filtered[idx] = a_G*f[idx_cell_data] + 
                        b_G*(f[idx_cell_data_z_B]    + f[idx_cell_data_z_F]) +
                        c_G*(f[idx_cell_data_z_BB]   + f[idx_cell_data_z_FF]) +
                        d_G*(f[idx_cell_data_z_BBB]  + f[idx_cell_data_z_FFF]) +
                        e_G*(f[idx_cell_data_z_BBBB] + f[idx_cell_data_z_FFFF]);
                }
            }
        }
    }
}
