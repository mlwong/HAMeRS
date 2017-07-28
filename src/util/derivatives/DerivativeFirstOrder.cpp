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
    const boost::shared_ptr<pdat::CellData<double> >& cell_data,
    const double dx,
    const hier::Box& domain,
    const int depth)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(derivative);
    TBOX_ASSERT(cell_data);
    TBOX_ASSERT(depth < cell_data->getDepth());
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = cell_data->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(derivative->getBox().numberCells() == interior_dims);
#endif
    
    // Get the number of ghost cells of the cell data and derivative data.
    const hier::IntVector num_ghosts_cell_data = cell_data->getGhostCellWidth();
    const hier::IntVector num_ghosts_derivative = derivative->getGhostCellWidth();
    
    // Get the dimensions of box that covers interior of patch plus ghost cells.
    const hier::Box ghost_box_cell_data = cell_data->getGhostBox();
    const hier::IntVector ghostcell_dims_cell_data = ghost_box_cell_data.numberCells();
    
    const hier::Box ghost_box_derivative = derivative->getGhostBox();
    const hier::IntVector ghostcell_dims_derivative = ghost_box_derivative.numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        if (num_ghosts_cell_data - num_ghosts_derivative < d_num_derivative_ghosts)
        {
            TBOX_ERROR(d_object_name
                << ": DerivativeFirstOrder::computeDerivative()\n"
                << "The cell data does not have enough ghost cells for taking derivative."
                << std::endl);
        }
#endif
        
        domain_lo = -num_ghosts_derivative;
        domain_dims = ghost_box_derivative.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        hier::Box shrinked_ghost_box_cell_data(ghost_box_cell_data);
        shrinked_ghost_box_cell_data.grow(-d_num_derivative_ghosts);
        
        TBOX_ASSERT(shrinked_ghost_box_cell_data.contains(domain));
        TBOX_ASSERT(ghost_box_derivative.contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
    // Get the pointer to the current depth component of the given cell data.
    double* u = cell_data->getPointer(depth);
    
    if (d_direction == DIRECTION::X_DIRECTION)
    {
        // Get the pointer to the derivative.
        double* dudx = derivative->getPointer(0);
        
        if (d_dim == tbox::Dimension(1))
        {
            /*
             * Get the local lower index, numbers of cells in each dimension and numbers of ghost cells.
             */
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_dim_0 = domain_dims[0];
            
            const int num_ghosts_0_cell_data = num_ghosts_cell_data[0];
            const int num_ghosts_0_derivative = num_ghosts_derivative[0];
            
            if (d_num_derivative_ghosts >= hier::IntVector::getOne(d_dim)*4)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute indices of current and neighboring cells.
                    const int idx_derivative = i + num_ghosts_0_derivative;
                    
                    const int idx_x_LLLL = i - 4 + num_ghosts_0_cell_data;
                    const int idx_x_LLL  = i - 3 + num_ghosts_0_cell_data;
                    const int idx_x_LL   = i - 2 + num_ghosts_0_cell_data;
                    const int idx_x_L    = i - 1 + num_ghosts_0_cell_data;
                    const int idx_x_R    = i + 1 + num_ghosts_0_cell_data;
                    const int idx_x_RR   = i + 2 + num_ghosts_0_cell_data;
                    const int idx_x_RRR  = i + 3 + num_ghosts_0_cell_data;
                    const int idx_x_RRRR = i + 4 + num_ghosts_0_cell_data;
                    
                    dudx[idx_derivative] = (double(-1)/double(280)*u[idx_x_RRRR] + double(4)/double(105)*u[idx_x_RRR] +
                                            double(-1)/double(5)*u[idx_x_RR] + double(4)/double(5)*u[idx_x_R] +
                                            double(-4)/double(5)*u[idx_x_L] + double(1)/double(5)*u[idx_x_LL] +
                                            double(-4)/double(105)*u[idx_x_LLL] + double(1)/double(280)*u[idx_x_LLLL])/dx;
                }
            }
            else if (d_num_derivative_ghosts == hier::IntVector::getOne(d_dim)*3)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute indices of current and neighboring cells.
                    const int idx_derivative = i + num_ghosts_0_derivative;
                    
                    const int idx_x_LLL = i - 3 + num_ghosts_0_cell_data;
                    const int idx_x_LL  = i - 2 + num_ghosts_0_cell_data;
                    const int idx_x_L   = i - 1 + num_ghosts_0_cell_data;
                    const int idx_x_R   = i + 1 + num_ghosts_0_cell_data;
                    const int idx_x_RR  = i + 2 + num_ghosts_0_cell_data;
                    const int idx_x_RRR = i + 3 + num_ghosts_0_cell_data;
                    
                    dudx[idx_derivative] = (double(1)/double(60)*u[idx_x_RRR] + double(-3)/double(20)*u[idx_x_RR] +
                                            double(3)/double(4)*u[idx_x_R] + double(-3)/double(4)*u[idx_x_L] +
                                            double(3)/double(20)*u[idx_x_LL] + double(-1)/double(60)*u[idx_x_LLL])/dx;
                }
            }
            else if (d_num_derivative_ghosts == hier::IntVector::getOne(d_dim)*2)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute indices of current and neighboring cells.
                    const int idx_derivative = i + num_ghosts_0_derivative;
                    
                    const int idx_x_LL = i - 2 + num_ghosts_0_cell_data;
                    const int idx_x_L  = i - 1 + num_ghosts_0_cell_data;
                    const int idx_x_R  = i + 1 + num_ghosts_0_cell_data;
                    const int idx_x_RR = i + 2 + num_ghosts_0_cell_data;
                    
                    dudx[idx_derivative] = (double(-1)/double(12)*u[idx_x_RR] + double(2)/double(3)*u[idx_x_R] +
                                            double(-2)/double(3)*u[idx_x_L] + double(1)/double(12)*u[idx_x_LL])/dx;
                }
            }
            else if (d_num_derivative_ghosts == hier::IntVector::getOne(d_dim))
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute indices of current and neighboring cells.
                    const int idx_derivative = i + num_ghosts_0_derivative;
                    
                    const int idx_x_L = i - 1 + num_ghosts_0_cell_data;
                    const int idx_x_R = i + 1 + num_ghosts_0_cell_data;
                    
                    dudx[idx_derivative] = (double(1)/double(2)*u[idx_x_R] + double(-1)/double(2)*u[idx_x_L])/dx;
                }
            }
        }
        else if (d_dim == tbox::Dimension(2))
        {
            /*
             * Get the local lower indices, numbers of cells in each dimension and numbers of ghost cells.
             */
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_lo_1 = domain_lo[1];
            const int domain_dim_0 = domain_dims[0];
            const int domain_dim_1 = domain_dims[1];
            
            const int num_ghosts_0_cell_data = num_ghosts_cell_data[0];
            const int num_ghosts_1_cell_data = num_ghosts_cell_data[1];
            const int ghostcell_dim_0_cell_data = ghostcell_dims_cell_data[0];
            
            const int num_ghosts_0_derivative = num_ghosts_derivative[0];
            const int num_ghosts_1_derivative = num_ghosts_derivative[1];
            const int ghostcell_dim_0_derivative = ghostcell_dims_derivative[0];
            
            if (d_num_derivative_ghosts >= hier::IntVector::getOne(d_dim)*4)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute indices of current and neighboring cells.
                        const int idx_derivative = (i + num_ghosts_0_derivative) +
                            (j + num_ghosts_1_derivative)*ghostcell_dim_0_derivative;
                        
                        const int idx_x_LLLL = (i - 4 + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        const int idx_x_LLL = (i - 3 + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        const int idx_x_LL = (i - 2 + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        const int idx_x_L = (i - 1 + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        const int idx_x_R = (i + 1 + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        const int idx_x_RR = (i + 2 + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        const int idx_x_RRR = (i + 3 + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        const int idx_x_RRRR = (i + 4 + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        dudx[idx_derivative] = (double(-1)/double(280)*u[idx_x_RRRR] + double(4)/double(105)*u[idx_x_RRR] +
                                                double(-1)/double(5)*u[idx_x_RR] + double(4)/double(5)*u[idx_x_R] +
                                                double(-4)/double(5)*u[idx_x_L] + double(1)/double(5)*u[idx_x_LL] +
                                                double(-4)/double(105)*u[idx_x_LLL] + double(1)/double(280)*u[idx_x_LLLL])/dx;
                    }
                }
            }
            else if (d_num_derivative_ghosts == hier::IntVector::getOne(d_dim)*3)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute indices of current and neighboring cells.
                        const int idx_derivative = (i + num_ghosts_0_derivative) +
                            (j + num_ghosts_1_derivative)*ghostcell_dim_0_derivative;
                        
                        const int idx_x_LLL = (i - 3 + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        const int idx_x_LL = (i - 2 + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        const int idx_x_L = (i - 1 + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        const int idx_x_R = (i + 1 + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        const int idx_x_RR = (i + 2 + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        const int idx_x_RRR = (i + 3 + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        dudx[idx_derivative] = (double(1)/double(60)*u[idx_x_RRR] + double(-3)/double(20)*u[idx_x_RR] +
                                                double(3)/double(4)*u[idx_x_R] + double(-3)/double(4)*u[idx_x_L] +
                                                double(3)/double(20)*u[idx_x_LL] + double(-1)/double(60)*u[idx_x_LLL])/dx;
                    }
                }
            }
            else if (d_num_derivative_ghosts == hier::IntVector::getOne(d_dim)*2)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute indices of current and neighboring cells.
                        const int idx_derivative = (i + num_ghosts_0_derivative) +
                            (j + num_ghosts_1_derivative)*ghostcell_dim_0_derivative;
                        
                        const int idx_x_LL = (i - 2 + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        const int idx_x_L = (i - 1 + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        const int idx_x_R = (i + 1 + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        const int idx_x_RR = (i + 2 + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        dudx[idx_derivative] = (double(-1)/double(12)*u[idx_x_RR] + double(2)/double(3)*u[idx_x_R] +
                                                double(-2)/double(3)*u[idx_x_L] + double(1)/double(12)*u[idx_x_LL])/dx;
                    }
                }
            }
            else if (d_num_derivative_ghosts == hier::IntVector::getOne(d_dim))
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute indices of current and neighboring cells.
                        const int idx_derivative = (i + num_ghosts_0_derivative) +
                            (j + num_ghosts_1_derivative)*ghostcell_dim_0_derivative;
                        
                        const int idx_x_L = (i - 1 + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        const int idx_x_R = (i + 1 + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        dudx[idx_derivative] = (double(1)/double(2)*u[idx_x_R] + double(-1)/double(2)*u[idx_x_L])/dx;
                    }
                }
            }
        }
        else if (d_dim == tbox::Dimension(3))
        {
            /*
             * Get the local lower indices, numbers of cells in each dimension and numbers of ghost cells.
             */
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_lo_1 = domain_lo[1];
            const int domain_lo_2 = domain_lo[2];
            const int domain_dim_0 = domain_dims[0];
            const int domain_dim_1 = domain_dims[1];
            const int domain_dim_2 = domain_dims[2];
            
            const int num_ghosts_0_cell_data = num_ghosts_cell_data[0];
            const int num_ghosts_1_cell_data = num_ghosts_cell_data[1];
            const int num_ghosts_2_cell_data = num_ghosts_cell_data[2];
            const int ghostcell_dim_0_cell_data = ghostcell_dims_cell_data[0];
            const int ghostcell_dim_1_cell_data = ghostcell_dims_cell_data[1];
            
            const int num_ghosts_0_derivative = num_ghosts_derivative[0];
            const int num_ghosts_1_derivative = num_ghosts_derivative[1];
            const int num_ghosts_2_derivative = num_ghosts_derivative[2];
            const int ghostcell_dim_0_derivative = ghostcell_dims_derivative[0];
            const int ghostcell_dim_1_derivative = ghostcell_dims_derivative[1];
            
            if (d_num_derivative_ghosts >= hier::IntVector::getOne(d_dim)*4)
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
                            const int idx_derivative = (i + num_ghosts_0_derivative) +
                                (j + num_ghosts_1_derivative)*ghostcell_dim_0_derivative +
                                (k + num_ghosts_2_derivative)*ghostcell_dim_0_derivative*
                                    ghostcell_dim_1_derivative;
                            
                            const int idx_x_LLLL = (i - 4 + num_ghosts_0_cell_data) +
                                (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_x_LLL = (i - 3 + num_ghosts_0_cell_data) +
                                (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_x_LL = (i - 2 + num_ghosts_0_cell_data) +
                                (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_x_L = (i - 1 + num_ghosts_0_cell_data) +
                                (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_x_R = (i + 1 + num_ghosts_0_cell_data) +
                                (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_x_RR = (i + 2 + num_ghosts_0_cell_data) +
                                (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_x_RRR = (i + 3 + num_ghosts_0_cell_data) +
                                (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_x_RRRR = (i + 4 + num_ghosts_0_cell_data) +
                                (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            dudx[idx_derivative] = (double(-1)/double(280)*u[idx_x_RRRR] + double(4)/double(105)*u[idx_x_RRR] +
                                                    double(-1)/double(5)*u[idx_x_RR] + double(4)/double(5)*u[idx_x_R] +
                                                    double(-4)/double(5)*u[idx_x_L] + double(1)/double(5)*u[idx_x_LL] +
                                                    double(-4)/double(105)*u[idx_x_LLL] + double(1)/double(280)*u[idx_x_LLLL])/dx;
                        }
                    }
                }
            }
            else if (d_num_derivative_ghosts == hier::IntVector::getOne(d_dim)*3)
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
                            const int idx_derivative = (i + num_ghosts_0_derivative) +
                                (j + num_ghosts_1_derivative)*ghostcell_dim_0_derivative +
                                (k + num_ghosts_2_derivative)*ghostcell_dim_0_derivative*
                                    ghostcell_dim_1_derivative;
                            
                            const int idx_x_LLL = (i - 3 + num_ghosts_0_cell_data) +
                                (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_x_LL = (i - 2 + num_ghosts_0_cell_data) +
                                (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_x_L = (i - 1 + num_ghosts_0_cell_data) +
                                (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_x_R = (i + 1 + num_ghosts_0_cell_data) +
                                (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_x_RR = (i + 2 + num_ghosts_0_cell_data) +
                                (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_x_RRR = (i + 3 + num_ghosts_0_cell_data) +
                                (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            dudx[idx_derivative] = (double(1)/double(60)*u[idx_x_RRR] + double(-3)/double(20)*u[idx_x_RR] +
                                                    double(3)/double(4)*u[idx_x_R] + double(-3)/double(4)*u[idx_x_L] +
                                                    double(3)/double(20)*u[idx_x_LL] + double(-1)/double(60)*u[idx_x_LLL])/dx;
                        }
                    }
                }
            }
            else if (d_num_derivative_ghosts == hier::IntVector::getOne(d_dim)*2)
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
                            const int idx_derivative = (i + num_ghosts_0_derivative) +
                                (j + num_ghosts_1_derivative)*ghostcell_dim_0_derivative +
                                (k + num_ghosts_2_derivative)*ghostcell_dim_0_derivative*
                                    ghostcell_dim_1_derivative;
                            
                            const int idx_x_LL = (i - 2 + num_ghosts_0_cell_data) +
                                (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_x_L = (i - 1 + num_ghosts_0_cell_data) +
                                (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_x_R = (i + 1 + num_ghosts_0_cell_data) +
                                (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_x_RR = (i + 2 + num_ghosts_0_cell_data) +
                                (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            dudx[idx_derivative] = (double(-1)/double(12)*u[idx_x_RR] + double(2)/double(3)*u[idx_x_R] +
                                                    double(-2)/double(3)*u[idx_x_L] + double(1)/double(12)*u[idx_x_LL])/dx;
                        }
                    }
                }
            }
            else if (d_num_derivative_ghosts == hier::IntVector::getOne(d_dim))
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
                            const int idx_derivative = (i + num_ghosts_0_derivative) +
                                (j + num_ghosts_1_derivative)*ghostcell_dim_0_derivative +
                                (k + num_ghosts_2_derivative)*ghostcell_dim_0_derivative*
                                    ghostcell_dim_1_derivative;
                            
                            const int idx_x_L = (i - 1 + num_ghosts_0_cell_data) +
                                (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_x_R = (i + 1 + num_ghosts_0_cell_data) +
                                (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            dudx[idx_derivative] = (double(1)/double(2)*u[idx_x_R] + double(-1)/double(2)*u[idx_x_L])/dx;
                        }
                    }
                }
            }
        }
    }
    else if (d_direction == DIRECTION::Y_DIRECTION)
    {
        // Get the pointer to the derivative.
        double* dudy = derivative->getPointer(0);
        
        const double& dy = dx;
        
        if (d_dim == tbox::Dimension(2))
        {
            /*
             * Get the local lower indices, numbers of cells in each dimension and numbers of ghost cells.
             */
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_lo_1 = domain_lo[1];
            const int domain_dim_0 = domain_dims[0];
            const int domain_dim_1 = domain_dims[1];
            
            const int num_ghosts_0_cell_data = num_ghosts_cell_data[0];
            const int num_ghosts_1_cell_data = num_ghosts_cell_data[1];
            const int ghostcell_dim_0_cell_data = ghostcell_dims_cell_data[0];
            
            const int num_ghosts_0_derivative = num_ghosts_derivative[0];
            const int num_ghosts_1_derivative = num_ghosts_derivative[1];
            const int ghostcell_dim_0_derivative = ghostcell_dims_derivative[0];
            
            if (d_num_derivative_ghosts >= hier::IntVector::getOne(d_dim)*4)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute indices of current and neighboring cells.
                        const int idx_derivative = (i + num_ghosts_0_derivative) +
                            (j + num_ghosts_1_derivative)*ghostcell_dim_0_derivative;
                        
                        const int idx_y_BBBB = (i + num_ghosts_0_cell_data) +
                            (j - 4 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        const int idx_y_BBB = (i + num_ghosts_0_cell_data) +
                            (j - 3 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        const int idx_y_BB = (i + num_ghosts_0_cell_data) +
                            (j - 2 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        const int idx_y_B = (i + num_ghosts_0_cell_data) +
                            (j - 1 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        const int idx_y_T = (i + num_ghosts_0_cell_data) +
                            (j + 1 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        const int idx_y_TT = (i + num_ghosts_0_cell_data) +
                            (j + 2 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        const int idx_y_TTT = (i + num_ghosts_0_cell_data) +
                            (j + 3 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        const int idx_y_TTTT = (i + num_ghosts_0_cell_data) +
                            (j + 4 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        dudy[idx_derivative] = (double(-1)/double(280)*u[idx_y_TTTT] + double(4)/double(105)*u[idx_y_TTT] +
                                                double(-1)/double(5)*u[idx_y_TT] + double(4)/double(5)*u[idx_y_T] +
                                                double(-4)/double(5)*u[idx_y_B] + double(1)/double(5)*u[idx_y_BB] +
                                                double(-4)/double(105)*u[idx_y_BBB] + double(1)/double(280)*u[idx_y_BBBB])/dy;
                    }
                }
            }
            else if (d_num_derivative_ghosts == hier::IntVector::getOne(d_dim)*3)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute indices of current and neighboring cells.
                        const int idx_derivative = (i + num_ghosts_0_derivative) +
                            (j + num_ghosts_1_derivative)*ghostcell_dim_0_derivative;
                        
                        const int idx_y_BBB = (i + num_ghosts_0_cell_data) +
                            (j - 3 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        const int idx_y_BB = (i + num_ghosts_0_cell_data) +
                            (j - 2 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        const int idx_y_B = (i + num_ghosts_0_cell_data) +
                            (j - 1 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        const int idx_y_T = (i + num_ghosts_0_cell_data) +
                            (j + 1 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        const int idx_y_TT = (i + num_ghosts_0_cell_data) +
                            (j + 2 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        const int idx_y_TTT = (i + num_ghosts_0_cell_data) +
                            (j + 3 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        dudy[idx_derivative] = (double(1)/double(60)*u[idx_y_TTT] + double(-3)/double(20)*u[idx_y_TT] +
                                                double(3)/double(4)*u[idx_y_T] + double(-3)/double(4)*u[idx_y_B] +
                                                double(3)/double(20)*u[idx_y_BB] + double(-1)/double(60)*u[idx_y_BBB])/dy;
                    }
                }
            }
            else if (d_num_derivative_ghosts == hier::IntVector::getOne(d_dim)*2)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute indices of current and neighboring cells.
                        const int idx_derivative = (i + num_ghosts_0_derivative) +
                            (j + num_ghosts_1_derivative)*ghostcell_dim_0_derivative;
                        
                        const int idx_y_BB = (i + num_ghosts_0_cell_data) +
                            (j - 2 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        const int idx_y_B = (i + num_ghosts_0_cell_data) +
                            (j - 1 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        const int idx_y_T = (i + num_ghosts_0_cell_data) +
                            (j + 1 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        const int idx_y_TT = (i + num_ghosts_0_cell_data) +
                            (j + 2 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        dudy[idx_derivative] = (double(-1)/double(12)*u[idx_y_TT] + double(2)/double(3)*u[idx_y_T] +
                                                double(-2)/double(3)*u[idx_y_B] + double(1)/double(12)*u[idx_y_BB])/dy;
                    }
                }
            }
            else if (d_num_derivative_ghosts == hier::IntVector::getOne(d_dim))
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute indices of current and neighboring cells.
                        const int idx_derivative = (i + num_ghosts_0_derivative) +
                            (j + num_ghosts_1_derivative)*ghostcell_dim_0_derivative;
                        
                        const int idx_y_B = (i + num_ghosts_0_cell_data) +
                            (j - 1 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        const int idx_y_T = (i + num_ghosts_0_cell_data) +
                            (j + 1 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                        
                        dudy[idx_derivative] = (double(1)/double(2)*u[idx_y_T] + double(-1)/double(2)*u[idx_y_B])/dy;
                    }
                }
            }
        }
        else if (d_dim == tbox::Dimension(3))
        {
            /*
             * Get the local lower indices, numbers of cells in each dimension and numbers of ghost cells.
             */
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_lo_1 = domain_lo[1];
            const int domain_lo_2 = domain_lo[2];
            const int domain_dim_0 = domain_dims[0];
            const int domain_dim_1 = domain_dims[1];
            const int domain_dim_2 = domain_dims[2];
            
            const int num_ghosts_0_cell_data = num_ghosts_cell_data[0];
            const int num_ghosts_1_cell_data = num_ghosts_cell_data[1];
            const int num_ghosts_2_cell_data = num_ghosts_cell_data[2];
            const int ghostcell_dim_0_cell_data = ghostcell_dims_cell_data[0];
            const int ghostcell_dim_1_cell_data = ghostcell_dims_cell_data[1];
            
            const int num_ghosts_0_derivative = num_ghosts_derivative[0];
            const int num_ghosts_1_derivative = num_ghosts_derivative[1];
            const int num_ghosts_2_derivative = num_ghosts_derivative[2];
            const int ghostcell_dim_0_derivative = ghostcell_dims_derivative[0];
            const int ghostcell_dim_1_derivative = ghostcell_dims_derivative[1];
            
            if (d_num_derivative_ghosts >= hier::IntVector::getOne(d_dim)*4)
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
                            const int idx_derivative = (i + num_ghosts_0_derivative) +
                                (j + num_ghosts_1_derivative)*ghostcell_dim_0_derivative +
                                (k + num_ghosts_2_derivative)*ghostcell_dim_0_derivative*
                                    ghostcell_dim_1_derivative;
                            
                            const int idx_y_BBBB = (i + num_ghosts_0_cell_data) +
                                (j - 4 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_y_BBB = (i + num_ghosts_0_cell_data) +
                                (j - 3 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_y_BB = (i + num_ghosts_0_cell_data) +
                                (j - 2 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_y_B = (i + num_ghosts_0_cell_data) +
                                (j - 1 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_y_T = (i + num_ghosts_0_cell_data) +
                                (j + 1 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_y_TT = (i + num_ghosts_0_cell_data) +
                                (j + 2 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_y_TTT = (i + num_ghosts_0_cell_data) +
                                (j + 3 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_y_TTTT = (i + num_ghosts_0_cell_data) +
                                (j + 4 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            dudy[idx_derivative] = (double(-1)/double(280)*u[idx_y_TTTT] + double(4)/double(105)*u[idx_y_TTT] +
                                                    double(-1)/double(5)*u[idx_y_TT] + double(4)/double(5)*u[idx_y_T] +
                                                    double(-4)/double(5)*u[idx_y_B] + double(1)/double(5)*u[idx_y_BB] +
                                                    double(-4)/double(105)*u[idx_y_BBB] + double(1)/double(280)*u[idx_y_BBBB])/dy;
                        }
                    }
                }
            }
            else if (d_num_derivative_ghosts == hier::IntVector::getOne(d_dim)*3)
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
                            const int idx_derivative = (i + num_ghosts_0_derivative) +
                                (j + num_ghosts_1_derivative)*ghostcell_dim_0_derivative +
                                (k + num_ghosts_2_derivative)*ghostcell_dim_0_derivative*
                                    ghostcell_dim_1_derivative;
                            
                            const int idx_y_BBB = (i + num_ghosts_0_cell_data) +
                                (j - 3 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_y_BB = (i + num_ghosts_0_cell_data) +
                                (j - 2 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_y_B = (i + num_ghosts_0_cell_data) +
                                (j - 1 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_y_T = (i + num_ghosts_0_cell_data) +
                                (j + 1 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_y_TT = (i + num_ghosts_0_cell_data) +
                                (j + 2 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_y_TTT = (i + num_ghosts_0_cell_data) +
                                (j + 3 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            dudy[idx_derivative] = (double(1)/double(60)*u[idx_y_TTT] + double(-3)/double(20)*u[idx_y_TT] +
                                                    double(3)/double(4)*u[idx_y_T] + double(-3)/double(4)*u[idx_y_B] +
                                                    double(3)/double(20)*u[idx_y_BB] + double(-1)/double(60)*u[idx_y_BBB])/dy;
                        }
                    }
                }
            }
            else if (d_num_derivative_ghosts == hier::IntVector::getOne(d_dim)*2)
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
                            const int idx_derivative = (i + num_ghosts_0_derivative) +
                                (j + num_ghosts_1_derivative)*ghostcell_dim_0_derivative +
                                (k + num_ghosts_2_derivative)*ghostcell_dim_0_derivative*
                                    ghostcell_dim_1_derivative;
                            
                            const int idx_y_BB = (i + num_ghosts_0_cell_data) +
                                (j - 2 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_y_B = (i + num_ghosts_0_cell_data) +
                                (j - 1 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_y_T = (i + num_ghosts_0_cell_data) +
                                (j + 1 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_y_TT = (i + num_ghosts_0_cell_data) +
                                (j + 2 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            dudy[idx_derivative] = (double(-1)/double(12)*u[idx_y_TT] + double(2)/double(3)*u[idx_y_T] +
                                                    double(-2)/double(3)*u[idx_y_B] + double(1)/double(12)*u[idx_y_BB])/dy;
                        }
                    }
                }
            }
            else if (d_num_derivative_ghosts == hier::IntVector::getOne(d_dim))
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
                            const int idx_derivative = (i + num_ghosts_0_derivative) +
                                (j + num_ghosts_1_derivative)*ghostcell_dim_0_derivative +
                                (k + num_ghosts_2_derivative)*ghostcell_dim_0_derivative*
                                    ghostcell_dim_1_derivative;
                            
                            const int idx_y_B = (i + num_ghosts_0_cell_data) +
                                (j - 1 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            const int idx_y_T = (i + num_ghosts_0_cell_data) +
                                (j + 1 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                                (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                    ghostcell_dim_1_cell_data;
                            
                            dudy[idx_derivative] = (double(1)/double(2)*u[idx_y_T] + double(-1)/double(2)*u[idx_y_B])/dy;
                        }
                    }
                }
            }
        }
    }
    else if (d_direction == DIRECTION::Z_DIRECTION)
    {
        // Get the pointer to the derivative.
        double* dudz = derivative->getPointer(0);
        
        const double& dz = dx;
        
        /*
         * Get the local lower indices, numbers of cells in each dimension and numbers of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_lo_2 = domain_lo[2];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        const int domain_dim_2 = domain_dims[2];
        
        const int num_ghosts_0_cell_data = num_ghosts_cell_data[0];
        const int num_ghosts_1_cell_data = num_ghosts_cell_data[1];
        const int num_ghosts_2_cell_data = num_ghosts_cell_data[2];
        const int ghostcell_dim_0_cell_data = ghostcell_dims_cell_data[0];
        const int ghostcell_dim_1_cell_data = ghostcell_dims_cell_data[1];
        
        const int num_ghosts_0_derivative = num_ghosts_derivative[0];
        const int num_ghosts_1_derivative = num_ghosts_derivative[1];
        const int num_ghosts_2_derivative = num_ghosts_derivative[2];
        const int ghostcell_dim_0_derivative = ghostcell_dims_derivative[0];
        const int ghostcell_dim_1_derivative = ghostcell_dims_derivative[1];
        
        if (d_num_derivative_ghosts >= hier::IntVector::getOne(d_dim)*4)
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
                        const int idx_derivative = (i + num_ghosts_0_derivative) +
                            (j + num_ghosts_1_derivative)*ghostcell_dim_0_derivative +
                            (k + num_ghosts_2_derivative)*ghostcell_dim_0_derivative*
                                ghostcell_dim_1_derivative;
                        
                        const int idx_z_BBBB = (i + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                            (k - 4 + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                ghostcell_dim_1_cell_data;
                        
                        const int idx_z_BBB = (i + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                            (k - 3 + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                ghostcell_dim_1_cell_data;
                        
                        const int idx_z_BB = (i + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                            (k - 2 + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                ghostcell_dim_1_cell_data;
                        
                        const int idx_z_B = (i + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                            (k - 1 + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                ghostcell_dim_1_cell_data;
                        
                        const int idx_z_F = (i + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                            (k + 1 + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                ghostcell_dim_1_cell_data;
                        
                        const int idx_z_FF = (i + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                            (k + 2 + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                ghostcell_dim_1_cell_data;
                        
                        const int idx_z_FFF = (i + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                            (k + 3 + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                ghostcell_dim_1_cell_data;
                        
                        const int idx_z_FFFF = (i + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                            (k + 4 + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                ghostcell_dim_1_cell_data;
                        
                        dudz[idx_derivative] = (double(-1)/double(280)*u[idx_z_FFFF] + double(4)/double(105)*u[idx_z_FFF] +
                                                double(-1)/double(5)*u[idx_z_FF] + double(4)/double(5)*u[idx_z_F] +
                                                double(-4)/double(5)*u[idx_z_B] + double(1)/double(5)*u[idx_z_BB] +
                                                double(-4)/double(105)*u[idx_z_BBB] + double(1)/double(280)*u[idx_z_BBBB])/dz;
                    }
                }
            }
        }
        else if (d_num_derivative_ghosts == hier::IntVector::getOne(d_dim)*3)
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
                        const int idx_derivative = (i + num_ghosts_0_derivative) +
                            (j + num_ghosts_1_derivative)*ghostcell_dim_0_derivative +
                            (k + num_ghosts_2_derivative)*ghostcell_dim_0_derivative*
                                ghostcell_dim_1_derivative;
                        
                        const int idx_z_BBB = (i + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                            (k - 3 + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                ghostcell_dim_1_cell_data;
                        
                        const int idx_z_BB = (i + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                            (k - 2 + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                ghostcell_dim_1_cell_data;
                        
                        const int idx_z_B = (i + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                            (k - 1 + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                ghostcell_dim_1_cell_data;
                        
                        const int idx_z_F = (i + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                            (k + 1 + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                ghostcell_dim_1_cell_data;
                        
                        const int idx_z_FF = (i + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                            (k + 2 + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                ghostcell_dim_1_cell_data;
                        
                        const int idx_z_FFF = (i + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                            (k + 3 + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                ghostcell_dim_1_cell_data;
                        
                        dudz[idx_derivative] = (double(1)/double(60)*u[idx_z_FFF] + double(-3)/double(20)*u[idx_z_FF] +
                                                double(3)/double(4)*u[idx_z_F] + double(-3)/double(4)*u[idx_z_B] +
                                                double(3)/double(20)*u[idx_z_BB] + double(-1)/double(60)*u[idx_z_BBB])/dz;
                    }
                }
            }
        }
        else if (d_num_derivative_ghosts == hier::IntVector::getOne(d_dim)*2)
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
                        const int idx_derivative = (i + num_ghosts_0_derivative) +
                            (j + num_ghosts_1_derivative)*ghostcell_dim_0_derivative +
                            (k + num_ghosts_2_derivative)*ghostcell_dim_0_derivative*
                                ghostcell_dim_1_derivative;
                        
                        const int idx_z_BB = (i + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                            (k - 2 + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                ghostcell_dim_1_cell_data;
                        
                        const int idx_z_B = (i + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                            (k - 1 + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                ghostcell_dim_1_cell_data;
                        
                        const int idx_z_F = (i + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                            (k + 1 + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                ghostcell_dim_1_cell_data;
                        
                        const int idx_z_FF = (i + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                            (k + 2 + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                ghostcell_dim_1_cell_data;
                        
                        dudz[idx_derivative] = (double(-1)/double(12)*u[idx_z_FF] + double(2)/double(3)*u[idx_z_F] +
                                                double(-2)/double(3)*u[idx_z_B] + double(1)/double(12)*u[idx_z_BB])/dz;
                    }
                }
            }
        }
        else if (d_num_derivative_ghosts == hier::IntVector::getOne(d_dim))
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
                        const int idx_derivative = (i + num_ghosts_0_derivative) +
                            (j + num_ghosts_1_derivative)*ghostcell_dim_0_derivative +
                            (k + num_ghosts_2_derivative)*ghostcell_dim_0_derivative*
                                ghostcell_dim_1_derivative;
                        
                        const int idx_z_B = (i + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                            (k - 1 + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                ghostcell_dim_1_cell_data;
                        
                        const int idx_z_F = (i + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                            (k + 1 + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                ghostcell_dim_1_cell_data;
                        
                        dudz[idx_derivative] = (double(1)/double(2)*u[idx_z_F] + double(-1)/double(2)*u[idx_z_B])/dz;
                    }
                }
            }
        }
    }
}
