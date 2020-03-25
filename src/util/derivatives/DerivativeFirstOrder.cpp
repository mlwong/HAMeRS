#include "util/derivatives/DerivativeFirstOrder.hpp"

DerivativeFirstOrder::DerivativeFirstOrder(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const DIRECTION::TYPE& direction,
    const int num_derivative_ghosts):
        Derivative(
            object_name,
            dim,
            direction,
            num_derivative_ghosts)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    if (num_derivative_ghosts < 1)
    {
        TBOX_ERROR(d_object_name
            << ": DerivativeFirstOrder::DerivativeFirstOrder()\n"
            << "At least one ghost cell is needed in any direction to take derivative."
            << std::endl);
    }
#endif
    
    switch (direction)
    {
        case (DIRECTION::X_DIRECTION):
        {
            d_num_derivative_ghosts[0] = num_derivative_ghosts;
            
            break;
        }
        case (DIRECTION::Y_DIRECTION):
        {
            d_num_derivative_ghosts[1] = num_derivative_ghosts;
            
            break;
        }
        case (DIRECTION::Z_DIRECTION):
        {
            d_num_derivative_ghosts[2] = num_derivative_ghosts;
            
            break;
        }
    }
}


/*
 * Compute the derivative with the given cell data.
 */
void
DerivativeFirstOrder::computeDerivative(
    boost::shared_ptr<pdat::CellData<double> >& derivative,
    const boost::shared_ptr<pdat::CellData<double> >& data,
    const double dx,
    const hier::Box& domain,
    const int depth_derivative,
    const int depth_data)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(derivative);
    TBOX_ASSERT(data);
    
    TBOX_ASSERT(depth_derivative < derivative->getDepth());
    TBOX_ASSERT(depth_data < data->getDepth());
#endif
    
    // Get the dimensions of box that covers interior of patch plus ghost cells.
    const hier::Box ghost_box_data = data->getGhostBox();
    const hier::IntVector ghostcell_dims_data = ghost_box_data.numberCells();
    
    const hier::Box ghost_box_derivative = derivative->getGhostBox();
    const hier::IntVector ghostcell_dims_derivative = ghost_box_derivative.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_data(d_dim);
    hier::IntVector offset_derivative(d_dim);
    
    if (domain.empty())
    {
        // Get the number of ghost cells of the cell data and derivative data.
        const hier::IntVector num_ghosts_derivative = derivative->getGhostCellWidth();
        const hier::IntVector num_ghosts_data = data->getGhostCellWidth();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data->getBox();
        
        TBOX_ASSERT(derivative->getBox().isSpatiallyEqual(interior_box));
        
        if (num_ghosts_data - num_ghosts_derivative < d_num_derivative_ghosts)
        {
            TBOX_ERROR(d_object_name
                << ": DerivativeFirstOrder::computeDerivative()\n"
                << "The cell data does not have enough ghost cells for taking derivative."
                << std::endl);
        }
#endif
        
        domain_lo = -num_ghosts_derivative;
        domain_dims = ghost_box_derivative.numberCells();
        
        offset_data = num_ghosts_data;
        offset_derivative = num_ghosts_derivative;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        hier::Box shrinked_ghost_box_data(ghost_box_data);
        shrinked_ghost_box_data.grow(-d_num_derivative_ghosts);
        
        TBOX_ASSERT(shrinked_ghost_box_data.contains(domain));
        TBOX_ASSERT(ghost_box_derivative.contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_data = domain.lower() - ghost_box_data.lower();
        offset_derivative = domain.lower() - ghost_box_derivative.lower();
    }
    
    // Get the pointer to the given cell data.
    double* u = data->getPointer(depth_data);
    
    if (d_direction == DIRECTION::X_DIRECTION)
    {
        // Get the pointer to the derivative.
        double* dudx = derivative->getPointer(depth_derivative);
        
        if (d_dim == tbox::Dimension(1))
        {
            /*
             * Get the local lower index, numbers of cells in each dimension and offsets.
             */
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_dim_0 = domain_dims[0];
            
            const int offset_0_data = offset_data[0];
            const int offset_0_derivative = offset_derivative[0];
            
            if (d_num_derivative_ghosts[0] == 4)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute indices of current and neighboring cells.
                    const int idx_derivative = i + offset_0_derivative;
                    
                    const int idx_x_LLLL = i - 4 + offset_0_data;
                    const int idx_x_LLL  = i - 3 + offset_0_data;
                    const int idx_x_LL   = i - 2 + offset_0_data;
                    const int idx_x_L    = i - 1 + offset_0_data;
                    const int idx_x_R    = i + 1 + offset_0_data;
                    const int idx_x_RR   = i + 2 + offset_0_data;
                    const int idx_x_RRR  = i + 3 + offset_0_data;
                    const int idx_x_RRRR = i + 4 + offset_0_data;
                    
                    dudx[idx_derivative] = (double(4)/double(5)*(u[idx_x_R] - u[idx_x_L]) +
                                            double(-1)/double(5)*(u[idx_x_RR] - u[idx_x_LL]) +
                                            double(4)/double(105)*(u[idx_x_RRR] - u[idx_x_LLL]) +
                                            double(-1)/double(280)*(u[idx_x_RRRR] - u[idx_x_LLLL]))/dx;
                }
            }
            else if (d_num_derivative_ghosts[0] == 3)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute indices of current and neighboring cells.
                    const int idx_derivative = i + offset_0_derivative;
                    
                    const int idx_x_LLL = i - 3 + offset_0_data;
                    const int idx_x_LL  = i - 2 + offset_0_data;
                    const int idx_x_L   = i - 1 + offset_0_data;
                    const int idx_x_R   = i + 1 + offset_0_data;
                    const int idx_x_RR  = i + 2 + offset_0_data;
                    const int idx_x_RRR = i + 3 + offset_0_data;
                    
                    dudx[idx_derivative] = (double(3)/double(4)*(u[idx_x_R] - u[idx_x_L]) +
                                            double(-3)/double(20)*(u[idx_x_RR] - u[idx_x_LL]) +
                                            double(1)/double(60)*(u[idx_x_RRR] - u[idx_x_LLL]))/dx;
                }
            }
            else if (d_num_derivative_ghosts[0] == 2)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute indices of current and neighboring cells.
                    const int idx_derivative = i + offset_0_derivative;
                    
                    const int idx_x_LL = i - 2 + offset_0_data;
                    const int idx_x_L  = i - 1 + offset_0_data;
                    const int idx_x_R  = i + 1 + offset_0_data;
                    const int idx_x_RR = i + 2 + offset_0_data;
                    
                    dudx[idx_derivative] = (double(2)/double(3)*(u[idx_x_R] - u[idx_x_L]) +
                                            double(-1)/double(12)*(u[idx_x_RR] - u[idx_x_LL]))/dx;
                }
            }
            else if (d_num_derivative_ghosts[0] == 1)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute indices of current and neighboring cells.
                    const int idx_derivative = i + offset_0_derivative;
                    
                    const int idx_x_L = i - 1 + offset_0_data;
                    const int idx_x_R = i + 1 + offset_0_data;
                    
                    dudx[idx_derivative] = (double(1)/double(2)*(u[idx_x_R] - u[idx_x_L]))/dx;
                }
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
            
            const int offset_0_data = offset_data[0];
            const int offset_1_data = offset_data[1];
            const int ghostcell_dim_0_data = ghostcell_dims_data[0];
            
            const int offset_0_derivative = offset_derivative[0];
            const int offset_1_derivative = offset_derivative[1];
            const int ghostcell_dim_0_derivative = ghostcell_dims_derivative[0];
            
            if (d_num_derivative_ghosts[0] == 4)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute indices of current and neighboring cells.
                        const int idx_derivative = (i + offset_0_derivative) +
                            (j + offset_1_derivative)*ghostcell_dim_0_derivative;
                        
                        const int idx_x_LLLL = (i - 4 + offset_0_data) +
                            (j + offset_1_data)*ghostcell_dim_0_data;
                        
                        const int idx_x_LLL = (i - 3 + offset_0_data) +
                            (j + offset_1_data)*ghostcell_dim_0_data;
                        
                        const int idx_x_LL = (i - 2 + offset_0_data) +
                            (j + offset_1_data)*ghostcell_dim_0_data;
                        
                        const int idx_x_L = (i - 1 + offset_0_data) +
                            (j + offset_1_data)*ghostcell_dim_0_data;
                        
                        const int idx_x_R = (i + 1 + offset_0_data) +
                            (j + offset_1_data)*ghostcell_dim_0_data;
                        
                        const int idx_x_RR = (i + 2 + offset_0_data) +
                            (j + offset_1_data)*ghostcell_dim_0_data;
                        
                        const int idx_x_RRR = (i + 3 + offset_0_data) +
                            (j + offset_1_data)*ghostcell_dim_0_data;
                        
                        const int idx_x_RRRR = (i + 4 + offset_0_data) +
                            (j + offset_1_data)*ghostcell_dim_0_data;
                        
                        dudx[idx_derivative] = (double(4)/double(5)*(u[idx_x_R] - u[idx_x_L]) +
                                                double(-1)/double(5)*(u[idx_x_RR] - u[idx_x_LL]) +
                                                double(4)/double(105)*(u[idx_x_RRR] - u[idx_x_LLL]) +
                                                double(-1)/double(280)*(u[idx_x_RRRR] - u[idx_x_LLLL]))/dx;
                    }
                }
            }
            else if (d_num_derivative_ghosts[0] == 3)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute indices of current and neighboring cells.
                        const int idx_derivative = (i + offset_0_derivative) +
                            (j + offset_1_derivative)*ghostcell_dim_0_derivative;
                        
                        const int idx_x_LLL = (i - 3 + offset_0_data) +
                            (j + offset_1_data)*ghostcell_dim_0_data;
                        
                        const int idx_x_LL = (i - 2 + offset_0_data) +
                            (j + offset_1_data)*ghostcell_dim_0_data;
                        
                        const int idx_x_L = (i - 1 + offset_0_data) +
                            (j + offset_1_data)*ghostcell_dim_0_data;
                        
                        const int idx_x_R = (i + 1 + offset_0_data) +
                            (j + offset_1_data)*ghostcell_dim_0_data;
                        
                        const int idx_x_RR = (i + 2 + offset_0_data) +
                            (j + offset_1_data)*ghostcell_dim_0_data;
                        
                        const int idx_x_RRR = (i + 3 + offset_0_data) +
                            (j + offset_1_data)*ghostcell_dim_0_data;
                        
                        dudx[idx_derivative] = (double(3)/double(4)*(u[idx_x_R] - u[idx_x_L]) +
                                                double(-3)/double(20)*(u[idx_x_RR] - u[idx_x_LL]) +
                                                double(1)/double(60)*(u[idx_x_RRR] - u[idx_x_LLL]))/dx;
                    }
                }
            }
            else if (d_num_derivative_ghosts[0] == 2)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute indices of current and neighboring cells.
                        const int idx_derivative = (i + offset_0_derivative) +
                            (j + offset_1_derivative)*ghostcell_dim_0_derivative;
                        
                        const int idx_x_LL = (i - 2 + offset_0_data) +
                            (j + offset_1_data)*ghostcell_dim_0_data;
                        
                        const int idx_x_L = (i - 1 + offset_0_data) +
                            (j + offset_1_data)*ghostcell_dim_0_data;
                        
                        const int idx_x_R = (i + 1 + offset_0_data) +
                            (j + offset_1_data)*ghostcell_dim_0_data;
                        
                        const int idx_x_RR = (i + 2 + offset_0_data) +
                            (j + offset_1_data)*ghostcell_dim_0_data;
                        
                        dudx[idx_derivative] = (double(2)/double(3)*(u[idx_x_R] - u[idx_x_L]) +
                                                double(-1)/double(12)*(u[idx_x_RR] - u[idx_x_LL]))/dx;
                    }
                }
            }
            else if (d_num_derivative_ghosts[0] == 1)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute indices of current and neighboring cells.
                        const int idx_derivative = (i + offset_0_derivative) +
                            (j + offset_1_derivative)*ghostcell_dim_0_derivative;
                        
                        const int idx_x_L = (i - 1 + offset_0_data) +
                            (j + offset_1_data)*ghostcell_dim_0_data;
                        
                        const int idx_x_R = (i + 1 + offset_0_data) +
                            (j + offset_1_data)*ghostcell_dim_0_data;
                        
                        dudx[idx_derivative] = (double(1)/double(2)*(u[idx_x_R] - u[idx_x_L]))/dx;
                    }
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
            
            const int offset_0_data = offset_data[0];
            const int offset_1_data = offset_data[1];
            const int offset_2_data = offset_data[2];
            const int ghostcell_dim_0_data = ghostcell_dims_data[0];
            const int ghostcell_dim_1_data = ghostcell_dims_data[1];
            
            const int offset_0_derivative = offset_derivative[0];
            const int offset_1_derivative = offset_derivative[1];
            const int offset_2_derivative = offset_derivative[2];
            const int ghostcell_dim_0_derivative = ghostcell_dims_derivative[0];
            const int ghostcell_dim_1_derivative = ghostcell_dims_derivative[1];
            
            if (d_num_derivative_ghosts[0] == 4)
            {
                for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                {
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute indices of current and neighboring cells.
                            const int idx_derivative = (i + offset_0_derivative) +
                                (j + offset_1_derivative)*ghostcell_dim_0_derivative +
                                (k + offset_2_derivative)*ghostcell_dim_0_derivative*
                                    ghostcell_dim_1_derivative;
                            
                            const int idx_x_LLLL = (i - 4 + offset_0_data) +
                                (j + offset_1_data)*ghostcell_dim_0_data +
                                (k + offset_2_data)*ghostcell_dim_0_data*
                                    ghostcell_dim_1_data;
                            
                            const int idx_x_LLL = (i - 3 + offset_0_data) +
                                (j + offset_1_data)*ghostcell_dim_0_data +
                                (k + offset_2_data)*ghostcell_dim_0_data*
                                    ghostcell_dim_1_data;
                            
                            const int idx_x_LL = (i - 2 + offset_0_data) +
                                (j + offset_1_data)*ghostcell_dim_0_data +
                                (k + offset_2_data)*ghostcell_dim_0_data*
                                    ghostcell_dim_1_data;
                            
                            const int idx_x_L = (i - 1 + offset_0_data) +
                                (j + offset_1_data)*ghostcell_dim_0_data +
                                (k + offset_2_data)*ghostcell_dim_0_data*
                                    ghostcell_dim_1_data;
                            
                            const int idx_x_R = (i + 1 + offset_0_data) +
                                (j + offset_1_data)*ghostcell_dim_0_data +
                                (k + offset_2_data)*ghostcell_dim_0_data*
                                    ghostcell_dim_1_data;
                            
                            const int idx_x_RR = (i + 2 + offset_0_data) +
                                (j + offset_1_data)*ghostcell_dim_0_data +
                                (k + offset_2_data)*ghostcell_dim_0_data*
                                    ghostcell_dim_1_data;
                            
                            const int idx_x_RRR = (i + 3 + offset_0_data) +
                                (j + offset_1_data)*ghostcell_dim_0_data +
                                (k + offset_2_data)*ghostcell_dim_0_data*
                                    ghostcell_dim_1_data;
                            
                            const int idx_x_RRRR = (i + 4 + offset_0_data) +
                                (j + offset_1_data)*ghostcell_dim_0_data +
                                (k + offset_2_data)*ghostcell_dim_0_data*
                                    ghostcell_dim_1_data;
                            
                            dudx[idx_derivative] = (double(4)/double(5)*(u[idx_x_R] - u[idx_x_L]) +
                                                    double(-1)/double(5)*(u[idx_x_RR] - u[idx_x_LL]) +
                                                    double(4)/double(105)*(u[idx_x_RRR] - u[idx_x_LLL]) +
                                                    double(-1)/double(280)*(u[idx_x_RRRR] - u[idx_x_LLLL]))/dx;
                        }
                    }
                }
            }
            else if (d_num_derivative_ghosts[0] == 3)
            {
                for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                {
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute indices of current and neighboring cells.
                            const int idx_derivative = (i + offset_0_derivative) +
                                (j + offset_1_derivative)*ghostcell_dim_0_derivative +
                                (k + offset_2_derivative)*ghostcell_dim_0_derivative*
                                    ghostcell_dim_1_derivative;
                            
                            const int idx_x_LLL = (i - 3 + offset_0_data) +
                                (j + offset_1_data)*ghostcell_dim_0_data +
                                (k + offset_2_data)*ghostcell_dim_0_data*
                                    ghostcell_dim_1_data;
                            
                            const int idx_x_LL = (i - 2 + offset_0_data) +
                                (j + offset_1_data)*ghostcell_dim_0_data +
                                (k + offset_2_data)*ghostcell_dim_0_data*
                                    ghostcell_dim_1_data;
                            
                            const int idx_x_L = (i - 1 + offset_0_data) +
                                (j + offset_1_data)*ghostcell_dim_0_data +
                                (k + offset_2_data)*ghostcell_dim_0_data*
                                    ghostcell_dim_1_data;
                            
                            const int idx_x_R = (i + 1 + offset_0_data) +
                                (j + offset_1_data)*ghostcell_dim_0_data +
                                (k + offset_2_data)*ghostcell_dim_0_data*
                                    ghostcell_dim_1_data;
                            
                            const int idx_x_RR = (i + 2 + offset_0_data) +
                                (j + offset_1_data)*ghostcell_dim_0_data +
                                (k + offset_2_data)*ghostcell_dim_0_data*
                                    ghostcell_dim_1_data;
                            
                            const int idx_x_RRR = (i + 3 + offset_0_data) +
                                (j + offset_1_data)*ghostcell_dim_0_data +
                                (k + offset_2_data)*ghostcell_dim_0_data*
                                    ghostcell_dim_1_data;
                            
                            dudx[idx_derivative] = (double(3)/double(4)*(u[idx_x_R] - u[idx_x_L]) +
                                                    double(-3)/double(20)*(u[idx_x_RR] - u[idx_x_LL]) +
                                                    double(1)/double(60)*(u[idx_x_RRR] - u[idx_x_LLL]))/dx;
                        }
                    }
                }
            }
            else if (d_num_derivative_ghosts[0] == 2)
            {
                for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                {
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute indices of current and neighboring cells.
                            const int idx_derivative = (i + offset_0_derivative) +
                                (j + offset_1_derivative)*ghostcell_dim_0_derivative +
                                (k + offset_2_derivative)*ghostcell_dim_0_derivative*
                                    ghostcell_dim_1_derivative;
                            
                            const int idx_x_LL = (i - 2 + offset_0_data) +
                                (j + offset_1_data)*ghostcell_dim_0_data +
                                (k + offset_2_data)*ghostcell_dim_0_data*
                                    ghostcell_dim_1_data;
                            
                            const int idx_x_L = (i - 1 + offset_0_data) +
                                (j + offset_1_data)*ghostcell_dim_0_data +
                                (k + offset_2_data)*ghostcell_dim_0_data*
                                    ghostcell_dim_1_data;
                            
                            const int idx_x_R = (i + 1 + offset_0_data) +
                                (j + offset_1_data)*ghostcell_dim_0_data +
                                (k + offset_2_data)*ghostcell_dim_0_data*
                                    ghostcell_dim_1_data;
                            
                            const int idx_x_RR = (i + 2 + offset_0_data) +
                                (j + offset_1_data)*ghostcell_dim_0_data +
                                (k + offset_2_data)*ghostcell_dim_0_data*
                                    ghostcell_dim_1_data;
                            
                            dudx[idx_derivative] = (double(2)/double(3)*(u[idx_x_R] - u[idx_x_L]) +
                                                    double(-1)/double(12)*(u[idx_x_RR] - u[idx_x_LL]))/dx;
                        }
                    }
                }
            }
            else if (d_num_derivative_ghosts[0] == 1)
            {
                for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                {
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute indices of current and neighboring cells.
                            const int idx_derivative = (i + offset_0_derivative) +
                                (j + offset_1_derivative)*ghostcell_dim_0_derivative +
                                (k + offset_2_derivative)*ghostcell_dim_0_derivative*
                                    ghostcell_dim_1_derivative;
                            
                            const int idx_x_L = (i - 1 + offset_0_data) +
                                (j + offset_1_data)*ghostcell_dim_0_data +
                                (k + offset_2_data)*ghostcell_dim_0_data*
                                    ghostcell_dim_1_data;
                            
                            const int idx_x_R = (i + 1 + offset_0_data) +
                                (j + offset_1_data)*ghostcell_dim_0_data +
                                (k + offset_2_data)*ghostcell_dim_0_data*
                                    ghostcell_dim_1_data;
                            
                            dudx[idx_derivative] = (double(1)/double(2)*(u[idx_x_R] - u[idx_x_L]))/dx;
                        }
                    }
                }
            }
        }
    }
    else if (d_direction == DIRECTION::Y_DIRECTION)
    {
        // Get the pointer to the derivative.
        double* dudy = derivative->getPointer(depth_derivative);
        
        const double& dy = dx;
        
        if (d_dim == tbox::Dimension(2))
        {
            /*
             * Get the local lower indices, numbers of cells in each dimension and offsets.
             */
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_lo_1 = domain_lo[1];
            const int domain_dim_0 = domain_dims[0];
            const int domain_dim_1 = domain_dims[1];
            
            const int offset_0_data = offset_data[0];
            const int offset_1_data = offset_data[1];
            const int ghostcell_dim_0_data = ghostcell_dims_data[0];
            
            const int offset_0_derivative = offset_derivative[0];
            const int offset_1_derivative = offset_derivative[1];
            const int ghostcell_dim_0_derivative = ghostcell_dims_derivative[0];
            
            if (d_num_derivative_ghosts[1] == 4)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute indices of current and neighboring cells.
                        const int idx_derivative = (i + offset_0_derivative) +
                            (j + offset_1_derivative)*ghostcell_dim_0_derivative;
                        
                        const int idx_y_BBBB = (i + offset_0_data) +
                            (j - 4 + offset_1_data)*ghostcell_dim_0_data;
                        
                        const int idx_y_BBB = (i + offset_0_data) +
                            (j - 3 + offset_1_data)*ghostcell_dim_0_data;
                        
                        const int idx_y_BB = (i + offset_0_data) +
                            (j - 2 + offset_1_data)*ghostcell_dim_0_data;
                        
                        const int idx_y_B = (i + offset_0_data) +
                            (j - 1 + offset_1_data)*ghostcell_dim_0_data;
                        
                        const int idx_y_T = (i + offset_0_data) +
                            (j + 1 + offset_1_data)*ghostcell_dim_0_data;
                        
                        const int idx_y_TT = (i + offset_0_data) +
                            (j + 2 + offset_1_data)*ghostcell_dim_0_data;
                        
                        const int idx_y_TTT = (i + offset_0_data) +
                            (j + 3 + offset_1_data)*ghostcell_dim_0_data;
                        
                        const int idx_y_TTTT = (i + offset_0_data) +
                            (j + 4 + offset_1_data)*ghostcell_dim_0_data;
                        
                        dudy[idx_derivative] = (double(4)/double(5)*(u[idx_y_T] - u[idx_y_B]) +
                                                double(-1)/double(5)*(u[idx_y_TT] - u[idx_y_BB]) +
                                                double(4)/double(105)*(u[idx_y_TTT] - u[idx_y_BBB]) +
                                                double(-1)/double(280)*(u[idx_y_TTTT] - u[idx_y_BBBB]))/dy;
                    }
                }
            }
            else if (d_num_derivative_ghosts[1] == 3)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute indices of current and neighboring cells.
                        const int idx_derivative = (i + offset_0_derivative) +
                            (j + offset_1_derivative)*ghostcell_dim_0_derivative;
                        
                        const int idx_y_BBB = (i + offset_0_data) +
                            (j - 3 + offset_1_data)*ghostcell_dim_0_data;
                        
                        const int idx_y_BB = (i + offset_0_data) +
                            (j - 2 + offset_1_data)*ghostcell_dim_0_data;
                        
                        const int idx_y_B = (i + offset_0_data) +
                            (j - 1 + offset_1_data)*ghostcell_dim_0_data;
                        
                        const int idx_y_T = (i + offset_0_data) +
                            (j + 1 + offset_1_data)*ghostcell_dim_0_data;
                        
                        const int idx_y_TT = (i + offset_0_data) +
                            (j + 2 + offset_1_data)*ghostcell_dim_0_data;
                        
                        const int idx_y_TTT = (i + offset_0_data) +
                            (j + 3 + offset_1_data)*ghostcell_dim_0_data;
                        
                        dudy[idx_derivative] = (double(3)/double(4)*(u[idx_y_T] - u[idx_y_B]) +
                                                double(-3)/double(20)*(u[idx_y_TT] - u[idx_y_BB]) +
                                                double(1)/double(60)*(u[idx_y_TTT] - u[idx_y_BBB]))/dy;
                    }
                }
            }
            else if (d_num_derivative_ghosts[1] == 2)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute indices of current and neighboring cells.
                        const int idx_derivative = (i + offset_0_derivative) +
                            (j + offset_1_derivative)*ghostcell_dim_0_derivative;
                        
                        const int idx_y_BB = (i + offset_0_data) +
                            (j - 2 + offset_1_data)*ghostcell_dim_0_data;
                        
                        const int idx_y_B = (i + offset_0_data) +
                            (j - 1 + offset_1_data)*ghostcell_dim_0_data;
                        
                        const int idx_y_T = (i + offset_0_data) +
                            (j + 1 + offset_1_data)*ghostcell_dim_0_data;
                        
                        const int idx_y_TT = (i + offset_0_data) +
                            (j + 2 + offset_1_data)*ghostcell_dim_0_data;
                        
                        dudy[idx_derivative] = (double(2)/double(3)*(u[idx_y_T] - u[idx_y_B]) +
                                                double(-1)/double(12)*(u[idx_y_TT] - u[idx_y_BB]))/dy;
                    }
                }
            }
            else if (d_num_derivative_ghosts[1] == 1)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute indices of current and neighboring cells.
                        const int idx_derivative = (i + offset_0_derivative) +
                            (j + offset_1_derivative)*ghostcell_dim_0_derivative;
                        
                        const int idx_y_B = (i + offset_0_data) +
                            (j - 1 + offset_1_data)*ghostcell_dim_0_data;
                        
                        const int idx_y_T = (i + offset_0_data) +
                            (j + 1 + offset_1_data)*ghostcell_dim_0_data;
                        
                        dudy[idx_derivative] = (double(1)/double(2)*(u[idx_y_T] - u[idx_y_B]))/dy;
                    }
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
            
            const int offset_0_data = offset_data[0];
            const int offset_1_data = offset_data[1];
            const int offset_2_data = offset_data[2];
            const int ghostcell_dim_0_data = ghostcell_dims_data[0];
            const int ghostcell_dim_1_data = ghostcell_dims_data[1];
            
            const int offset_0_derivative = offset_derivative[0];
            const int offset_1_derivative = offset_derivative[1];
            const int offset_2_derivative = offset_derivative[2];
            const int ghostcell_dim_0_derivative = ghostcell_dims_derivative[0];
            const int ghostcell_dim_1_derivative = ghostcell_dims_derivative[1];
            
            if (d_num_derivative_ghosts[1] == 4)
            {
                for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                {
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute indices of current and neighboring cells.
                            const int idx_derivative = (i + offset_0_derivative) +
                                (j + offset_1_derivative)*ghostcell_dim_0_derivative +
                                (k + offset_2_derivative)*ghostcell_dim_0_derivative*
                                    ghostcell_dim_1_derivative;
                            
                            const int idx_y_BBBB = (i + offset_0_data) +
                                (j - 4 + offset_1_data)*ghostcell_dim_0_data +
                                (k + offset_2_data)*ghostcell_dim_0_data*
                                    ghostcell_dim_1_data;
                            
                            const int idx_y_BBB = (i + offset_0_data) +
                                (j - 3 + offset_1_data)*ghostcell_dim_0_data +
                                (k + offset_2_data)*ghostcell_dim_0_data*
                                    ghostcell_dim_1_data;
                            
                            const int idx_y_BB = (i + offset_0_data) +
                                (j - 2 + offset_1_data)*ghostcell_dim_0_data +
                                (k + offset_2_data)*ghostcell_dim_0_data*
                                    ghostcell_dim_1_data;
                            
                            const int idx_y_B = (i + offset_0_data) +
                                (j - 1 + offset_1_data)*ghostcell_dim_0_data +
                                (k + offset_2_data)*ghostcell_dim_0_data*
                                    ghostcell_dim_1_data;
                            
                            const int idx_y_T = (i + offset_0_data) +
                                (j + 1 + offset_1_data)*ghostcell_dim_0_data +
                                (k + offset_2_data)*ghostcell_dim_0_data*
                                    ghostcell_dim_1_data;
                            
                            const int idx_y_TT = (i + offset_0_data) +
                                (j + 2 + offset_1_data)*ghostcell_dim_0_data +
                                (k + offset_2_data)*ghostcell_dim_0_data*
                                    ghostcell_dim_1_data;
                            
                            const int idx_y_TTT = (i + offset_0_data) +
                                (j + 3 + offset_1_data)*ghostcell_dim_0_data +
                                (k + offset_2_data)*ghostcell_dim_0_data*
                                    ghostcell_dim_1_data;
                            
                            const int idx_y_TTTT = (i + offset_0_data) +
                                (j + 4 + offset_1_data)*ghostcell_dim_0_data +
                                (k + offset_2_data)*ghostcell_dim_0_data*
                                    ghostcell_dim_1_data;
                            
                            dudy[idx_derivative] = (double(4)/double(5)*(u[idx_y_T] - u[idx_y_B]) +
                                                    double(-1)/double(5)*(u[idx_y_TT] - u[idx_y_BB]) +
                                                    double(4)/double(105)*(u[idx_y_TTT] - u[idx_y_BBB]) +
                                                    double(-1)/double(280)*(u[idx_y_TTTT] - u[idx_y_BBBB]))/dy;
                        }
                    }
                }
            }
            else if (d_num_derivative_ghosts[1] == 3)
            {
                for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                {
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute indices of current and neighboring cells.
                            const int idx_derivative = (i + offset_0_derivative) +
                                (j + offset_1_derivative)*ghostcell_dim_0_derivative +
                                (k + offset_2_derivative)*ghostcell_dim_0_derivative*
                                    ghostcell_dim_1_derivative;
                            
                            const int idx_y_BBB = (i + offset_0_data) +
                                (j - 3 + offset_1_data)*ghostcell_dim_0_data +
                                (k + offset_2_data)*ghostcell_dim_0_data*
                                    ghostcell_dim_1_data;
                            
                            const int idx_y_BB = (i + offset_0_data) +
                                (j - 2 + offset_1_data)*ghostcell_dim_0_data +
                                (k + offset_2_data)*ghostcell_dim_0_data*
                                    ghostcell_dim_1_data;
                            
                            const int idx_y_B = (i + offset_0_data) +
                                (j - 1 + offset_1_data)*ghostcell_dim_0_data +
                                (k + offset_2_data)*ghostcell_dim_0_data*
                                    ghostcell_dim_1_data;
                            
                            const int idx_y_T = (i + offset_0_data) +
                                (j + 1 + offset_1_data)*ghostcell_dim_0_data +
                                (k + offset_2_data)*ghostcell_dim_0_data*
                                    ghostcell_dim_1_data;
                            
                            const int idx_y_TT = (i + offset_0_data) +
                                (j + 2 + offset_1_data)*ghostcell_dim_0_data +
                                (k + offset_2_data)*ghostcell_dim_0_data*
                                    ghostcell_dim_1_data;
                            
                            const int idx_y_TTT = (i + offset_0_data) +
                                (j + 3 + offset_1_data)*ghostcell_dim_0_data +
                                (k + offset_2_data)*ghostcell_dim_0_data*
                                    ghostcell_dim_1_data;
                            
                            dudy[idx_derivative] = (double(3)/double(4)*(u[idx_y_T] - u[idx_y_B]) +
                                                    double(-3)/double(20)*(u[idx_y_TT] - u[idx_y_BB]) +
                                                    double(1)/double(60)*(u[idx_y_TTT] - u[idx_y_BBB]))/dy;
                        }
                    }
                }
            }
            else if (d_num_derivative_ghosts[1] == 2)
            {
                for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                {
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute indices of current and neighboring cells.
                            const int idx_derivative = (i + offset_0_derivative) +
                                (j + offset_1_derivative)*ghostcell_dim_0_derivative +
                                (k + offset_2_derivative)*ghostcell_dim_0_derivative*
                                    ghostcell_dim_1_derivative;
                            
                            const int idx_y_BB = (i + offset_0_data) +
                                (j - 2 + offset_1_data)*ghostcell_dim_0_data +
                                (k + offset_2_data)*ghostcell_dim_0_data*
                                    ghostcell_dim_1_data;
                            
                            const int idx_y_B = (i + offset_0_data) +
                                (j - 1 + offset_1_data)*ghostcell_dim_0_data +
                                (k + offset_2_data)*ghostcell_dim_0_data*
                                    ghostcell_dim_1_data;
                            
                            const int idx_y_T = (i + offset_0_data) +
                                (j + 1 + offset_1_data)*ghostcell_dim_0_data +
                                (k + offset_2_data)*ghostcell_dim_0_data*
                                    ghostcell_dim_1_data;
                            
                            const int idx_y_TT = (i + offset_0_data) +
                                (j + 2 + offset_1_data)*ghostcell_dim_0_data +
                                (k + offset_2_data)*ghostcell_dim_0_data*
                                    ghostcell_dim_1_data;
                            
                            dudy[idx_derivative] = (double(2)/double(3)*(u[idx_y_T] - u[idx_y_B]) +
                                                    double(-1)/double(12)*(u[idx_y_TT] - u[idx_y_BB]))/dy;
                        }
                    }
                }
            }
            else if (d_num_derivative_ghosts[1] == 1)
            {
                for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                {
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute indices of current and neighboring cells.
                            const int idx_derivative = (i + offset_0_derivative) +
                                (j + offset_1_derivative)*ghostcell_dim_0_derivative +
                                (k + offset_2_derivative)*ghostcell_dim_0_derivative*
                                    ghostcell_dim_1_derivative;
                            
                            const int idx_y_B = (i + offset_0_data) +
                                (j - 1 + offset_1_data)*ghostcell_dim_0_data +
                                (k + offset_2_data)*ghostcell_dim_0_data*
                                    ghostcell_dim_1_data;
                            
                            const int idx_y_T = (i + offset_0_data) +
                                (j + 1 + offset_1_data)*ghostcell_dim_0_data +
                                (k + offset_2_data)*ghostcell_dim_0_data*
                                    ghostcell_dim_1_data;
                            
                            dudy[idx_derivative] = (double(1)/double(2)*(u[idx_y_T] - u[idx_y_B]))/dy;
                        }
                    }
                }
            }
        }
    }
    else if (d_direction == DIRECTION::Z_DIRECTION)
    {
        // Get the pointer to the derivative.
        double* dudz = derivative->getPointer(depth_derivative);
        
        const double& dz = dx;
        
        /*
         * Get the local lower indices, numbers of cells in each dimension and offsets.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_lo_2 = domain_lo[2];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        const int domain_dim_2 = domain_dims[2];
        
        const int offset_0_data = offset_data[0];
        const int offset_1_data = offset_data[1];
        const int offset_2_data = offset_data[2];
        const int ghostcell_dim_0_data = ghostcell_dims_data[0];
        const int ghostcell_dim_1_data = ghostcell_dims_data[1];
        
        const int offset_0_derivative = offset_derivative[0];
        const int offset_1_derivative = offset_derivative[1];
        const int offset_2_derivative = offset_derivative[2];
        const int ghostcell_dim_0_derivative = ghostcell_dims_derivative[0];
        const int ghostcell_dim_1_derivative = ghostcell_dims_derivative[1];
        
        if (d_num_derivative_ghosts[2] == 4)
        {
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute indices of current and neighboring cells.
                        const int idx_derivative = (i + offset_0_derivative) +
                            (j + offset_1_derivative)*ghostcell_dim_0_derivative +
                            (k + offset_2_derivative)*ghostcell_dim_0_derivative*
                                ghostcell_dim_1_derivative;
                        
                        const int idx_z_BBBB = (i + offset_0_data) +
                            (j + offset_1_data)*ghostcell_dim_0_data +
                            (k - 4 + offset_2_data)*ghostcell_dim_0_data*
                                ghostcell_dim_1_data;
                        
                        const int idx_z_BBB = (i + offset_0_data) +
                            (j + offset_1_data)*ghostcell_dim_0_data +
                            (k - 3 + offset_2_data)*ghostcell_dim_0_data*
                                ghostcell_dim_1_data;
                        
                        const int idx_z_BB = (i + offset_0_data) +
                            (j + offset_1_data)*ghostcell_dim_0_data +
                            (k - 2 + offset_2_data)*ghostcell_dim_0_data*
                                ghostcell_dim_1_data;
                        
                        const int idx_z_B = (i + offset_0_data) +
                            (j + offset_1_data)*ghostcell_dim_0_data +
                            (k - 1 + offset_2_data)*ghostcell_dim_0_data*
                                ghostcell_dim_1_data;
                        
                        const int idx_z_F = (i + offset_0_data) +
                            (j + offset_1_data)*ghostcell_dim_0_data +
                            (k + 1 + offset_2_data)*ghostcell_dim_0_data*
                                ghostcell_dim_1_data;
                        
                        const int idx_z_FF = (i + offset_0_data) +
                            (j + offset_1_data)*ghostcell_dim_0_data +
                            (k + 2 + offset_2_data)*ghostcell_dim_0_data*
                                ghostcell_dim_1_data;
                        
                        const int idx_z_FFF = (i + offset_0_data) +
                            (j + offset_1_data)*ghostcell_dim_0_data +
                            (k + 3 + offset_2_data)*ghostcell_dim_0_data*
                                ghostcell_dim_1_data;
                        
                        const int idx_z_FFFF = (i + offset_0_data) +
                            (j + offset_1_data)*ghostcell_dim_0_data +
                            (k + 4 + offset_2_data)*ghostcell_dim_0_data*
                                ghostcell_dim_1_data;
                        
                        dudz[idx_derivative] = (double(4)/double(5)*(u[idx_z_F] - u[idx_z_B]) +
                                                double(-1)/double(5)*(u[idx_z_FF] - u[idx_z_BB]) +
                                                double(4)/double(105)*(u[idx_z_FFF] - u[idx_z_BBB]) +
                                                double(-1)/double(280)*(u[idx_z_FFFF] - u[idx_z_BBBB]))/dz;
                    }
                }
            }
        }
        else if (d_num_derivative_ghosts[2] == 3)
        {
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute indices of current and neighboring cells.
                        const int idx_derivative = (i + offset_0_derivative) +
                            (j + offset_1_derivative)*ghostcell_dim_0_derivative +
                            (k + offset_2_derivative)*ghostcell_dim_0_derivative*
                                ghostcell_dim_1_derivative;
                        
                        const int idx_z_BBB = (i + offset_0_data) +
                            (j + offset_1_data)*ghostcell_dim_0_data +
                            (k - 3 + offset_2_data)*ghostcell_dim_0_data*
                                ghostcell_dim_1_data;
                        
                        const int idx_z_BB = (i + offset_0_data) +
                            (j + offset_1_data)*ghostcell_dim_0_data +
                            (k - 2 + offset_2_data)*ghostcell_dim_0_data*
                                ghostcell_dim_1_data;
                        
                        const int idx_z_B = (i + offset_0_data) +
                            (j + offset_1_data)*ghostcell_dim_0_data +
                            (k - 1 + offset_2_data)*ghostcell_dim_0_data*
                                ghostcell_dim_1_data;
                        
                        const int idx_z_F = (i + offset_0_data) +
                            (j + offset_1_data)*ghostcell_dim_0_data +
                            (k + 1 + offset_2_data)*ghostcell_dim_0_data*
                                ghostcell_dim_1_data;
                        
                        const int idx_z_FF = (i + offset_0_data) +
                            (j + offset_1_data)*ghostcell_dim_0_data +
                            (k + 2 + offset_2_data)*ghostcell_dim_0_data*
                                ghostcell_dim_1_data;
                        
                        const int idx_z_FFF = (i + offset_0_data) +
                            (j + offset_1_data)*ghostcell_dim_0_data +
                            (k + 3 + offset_2_data)*ghostcell_dim_0_data*
                                ghostcell_dim_1_data;
                        
                        dudz[idx_derivative] = (double(3)/double(4)*(u[idx_z_F] - u[idx_z_B]) +
                                                double(-3)/double(20)*(u[idx_z_FF] - u[idx_z_BB]) +
                                                double(1)/double(60)*(u[idx_z_FFF] - u[idx_z_BBB]))/dz;
                    }
                }
            }
        }
        else if (d_num_derivative_ghosts[2] == 2)
        {
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute indices of current and neighboring cells.
                        const int idx_derivative = (i + offset_0_derivative) +
                            (j + offset_1_derivative)*ghostcell_dim_0_derivative +
                            (k + offset_2_derivative)*ghostcell_dim_0_derivative*
                                ghostcell_dim_1_derivative;
                        
                        const int idx_z_BB = (i + offset_0_data) +
                            (j + offset_1_data)*ghostcell_dim_0_data +
                            (k - 2 + offset_2_data)*ghostcell_dim_0_data*
                                ghostcell_dim_1_data;
                        
                        const int idx_z_B = (i + offset_0_data) +
                            (j + offset_1_data)*ghostcell_dim_0_data +
                            (k - 1 + offset_2_data)*ghostcell_dim_0_data*
                                ghostcell_dim_1_data;
                        
                        const int idx_z_F = (i + offset_0_data) +
                            (j + offset_1_data)*ghostcell_dim_0_data +
                            (k + 1 + offset_2_data)*ghostcell_dim_0_data*
                                ghostcell_dim_1_data;
                        
                        const int idx_z_FF = (i + offset_0_data) +
                            (j + offset_1_data)*ghostcell_dim_0_data +
                            (k + 2 + offset_2_data)*ghostcell_dim_0_data*
                                ghostcell_dim_1_data;
                        
                        dudz[idx_derivative] = (double(2)/double(3)*(u[idx_z_F] - u[idx_z_B]) +
                                                double(-1)/double(12)*(u[idx_z_FF] - u[idx_z_BB]))/dz;
                    }
                }
            }
        }
        else if (d_num_derivative_ghosts[2] == 1)
        {
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute indices of current and neighboring cells.
                        const int idx_derivative = (i + offset_0_derivative) +
                            (j + offset_1_derivative)*ghostcell_dim_0_derivative +
                            (k + offset_2_derivative)*ghostcell_dim_0_derivative*
                                ghostcell_dim_1_derivative;
                        
                        const int idx_z_B = (i + offset_0_data) +
                            (j + offset_1_data)*ghostcell_dim_0_data +
                            (k - 1 + offset_2_data)*ghostcell_dim_0_data*
                                ghostcell_dim_1_data;
                        
                        const int idx_z_F = (i + offset_0_data) +
                            (j + offset_1_data)*ghostcell_dim_0_data +
                            (k + 1 + offset_2_data)*ghostcell_dim_0_data*
                                ghostcell_dim_1_data;
                        
                        dudz[idx_derivative] = (double(1)/double(2)*(u[idx_z_F] - u[idx_z_B]))/dz;
                    }
                }
            }
        }
    }
}
