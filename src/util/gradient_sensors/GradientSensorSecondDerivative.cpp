#include "util/gradient_sensors/GradientSensorSecondDerivative.hpp"

GradientSensorSecondDerivative::GradientSensorSecondDerivative(
    const std::string& object_name,
    const tbox::Dimension& dim):
        GradientSensor(
            object_name,
            dim)
{
    d_num_gradient_ghosts = hier::IntVector::getOne(d_dim);
}


/*
 * Compute the gradient with the given cell data.
 */
void
GradientSensorSecondDerivative::computeGradient(
    hier::Patch& patch,
    boost::shared_ptr<pdat::CellData<double> > cell_data,
    boost::shared_ptr<pdat::CellData<double> > gradient,
    int depth)
{
    // Declare a null pointer.
    boost::shared_ptr<pdat::CellData<double> > variable_local_mean;
    
    computeGradientWithVariableLocalMean(
        patch,
        cell_data,
        gradient,
        variable_local_mean,
        depth);
}


/*
 * Compute the gradient and the local mean of the given cell data.
 */
void
GradientSensorSecondDerivative::computeGradientWithVariableLocalMean(
    hier::Patch& patch,
    boost::shared_ptr<pdat::CellData<double> > cell_data,
    boost::shared_ptr<pdat::CellData<double> > gradient,
    boost::shared_ptr<pdat::CellData<double> > variable_local_mean,
    int depth)
{
    if (cell_data->getGhostCellWidth() < d_num_gradient_ghosts)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The ghost cell width is smaller than required."
            << std::endl);
    }
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the number of ghost cells of the cell data and gradient data.
    const hier::IntVector num_ghosts_cell_data = cell_data->getGhostCellWidth();
    const hier::IntVector num_ghosts_gradient = gradient->getGhostCellWidth();
    
    // Get the dimensions of box that covers interior of patch plus ghost cells.
    const hier::Box ghost_box_cell_data = cell_data->getGhostBox();
    const hier::IntVector ghostcell_dims_cell_data = ghost_box_cell_data.numberCells();
    
    const hier::Box ghost_box_gradient = gradient->getGhostBox();
    const hier::IntVector ghostcell_dims_gradient = ghost_box_gradient.numberCells();
    
    // Determine whether local mean is requred to be compted.
    bool compute_variable_local_mean = false;
    if (variable_local_mean)
    {
        compute_variable_local_mean = true;
        
        TBOX_ASSERT(variable_local_mean->getGhostBox().numberCells() == ghostcell_dims_gradient);
    }
    
    // Get the pointer to the current depth component of the given cell data.
    double* f = cell_data->getPointer(depth);
    
    // Get the pointers to the gradient and local mean.
    double* w = gradient->getPointer(0);
    double* f_mean = nullptr;
    
    if (compute_variable_local_mean)
    {
        f_mean = variable_local_mean->getPointer(0);
    }
    
    if (d_dim == tbox::Dimension(1))
    {
        const int interior_dim_0 = interior_dims[0];
        
        const int num_ghosts_0_cell_data = num_ghosts_cell_data[0];
        const int num_ghosts_0_gradient = num_ghosts_gradient[0];
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = 0; i < interior_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx = i + num_ghosts_0_gradient;
            const int idx_cell_data = i + num_ghosts_0_cell_data;
            const int idx_cell_data_x_L = i - 1 + num_ghosts_0_cell_data;
            const int idx_cell_data_x_R = i + 1 + num_ghosts_0_cell_data;
            
            w[idx] = fabs(f[idx_cell_data_x_R] - 2*f[idx_cell_data] + f[idx_cell_data_x_L]);
        }
        
        if (compute_variable_local_mean)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx = i + num_ghosts_0_gradient;
                const int idx_cell_data = i + num_ghosts_0_cell_data;
                const int idx_cell_data_x_L = i - 1 + num_ghosts_0_cell_data;
                const int idx_cell_data_x_R = i + 1 + num_ghosts_0_cell_data;
                
                f_mean[idx] = f[idx_cell_data_x_R] + 2*f[idx_cell_data] + f[idx_cell_data_x_L];
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        const int num_ghosts_0_cell_data = num_ghosts_cell_data[0];
        const int num_ghosts_1_cell_data = num_ghosts_cell_data[1];
        const int ghostcell_dim_0_cell_data = ghostcell_dims_cell_data[0];
        
        const int num_ghosts_0_gradient = num_ghosts_gradient[0];
        const int num_ghosts_1_gradient = num_ghosts_gradient[1];
        const int ghostcell_dim_0_gradient = ghostcell_dims_gradient[0];
        
        for (int j = 0; j < interior_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx = (i + num_ghosts_0_gradient) +
                    (j + num_ghosts_1_gradient)*ghostcell_dim_0_gradient;
                
                const int idx_cell_data = (i + num_ghosts_0_cell_data) +
                    (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                
                const int idx_cell_data_x_L = (i - 1 + num_ghosts_0_cell_data) +
                    (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                
                const int idx_cell_data_x_R = (i + 1 + num_ghosts_0_cell_data) +
                    (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                
                const int idx_cell_data_y_B = (i + num_ghosts_0_cell_data) +
                    (j - 1 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                
                const int idx_cell_data_y_T = (i + num_ghosts_0_cell_data) +
                    (j + 1 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                
                double w_x = f[idx_cell_data_x_R] - 2*f[idx_cell_data] + f[idx_cell_data_x_L];
                double w_y = f[idx_cell_data_y_T] - 2*f[idx_cell_data] + f[idx_cell_data_y_B];
                
                w[idx] = sqrt(w_x*w_x + w_y*w_y);
            }
        }
        
        if (compute_variable_local_mean)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx = (i + num_ghosts_0_gradient) +
                        (j + num_ghosts_1_gradient)*ghostcell_dim_0_gradient;
                    
                    const int idx_cell_data = (i + num_ghosts_0_cell_data) +
                        (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                    
                    const int idx_cell_data_x_L = (i - 1 + num_ghosts_0_cell_data) +
                        (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                    
                    const int idx_cell_data_x_R = (i + 1 + num_ghosts_0_cell_data) +
                        (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                    
                    const int idx_cell_data_y_B = (i + num_ghosts_0_cell_data) +
                        (j - 1 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                    
                    const int idx_cell_data_y_T = (i + num_ghosts_0_cell_data) +
                        (j + 1 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                    
                    double f_mean_x = f[idx_cell_data_x_R] + 2*f[idx_cell_data] + f[idx_cell_data_x_L];
                    double f_mean_y = f[idx_cell_data_y_T] + 2*f[idx_cell_data] + f[idx_cell_data_y_B];
                    
                    f_mean[idx] = sqrt(f_mean_x*f_mean_x + f_mean_y*f_mean_y);
                }
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        const int interior_dim_2 = interior_dims[2];
        
        const int num_ghosts_0_cell_data = num_ghosts_cell_data[0];
        const int num_ghosts_1_cell_data = num_ghosts_cell_data[1];
        const int num_ghosts_2_cell_data = num_ghosts_cell_data[2];
        const int ghostcell_dim_0_cell_data = ghostcell_dims_cell_data[0];
        const int ghostcell_dim_1_cell_data = ghostcell_dims_cell_data[1];
        
        const int num_ghosts_0_gradient = num_ghosts_gradient[0];
        const int num_ghosts_1_gradient = num_ghosts_gradient[1];
        const int num_ghosts_2_gradient = num_ghosts_gradient[2];
        const int ghostcell_dim_0_gradient = ghostcell_dims_gradient[0];
        const int ghostcell_dim_1_gradient = ghostcell_dims_gradient[1];
        
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx = (i + num_ghosts_0_gradient) +
                        (j + num_ghosts_1_gradient)*ghostcell_dim_0_gradient +
                        (k + num_ghosts_2_gradient)*ghostcell_dim_0_gradient*
                            ghostcell_dim_1_gradient;
                    
                    const int idx_cell_data = (i + num_ghosts_0_cell_data) +
                        (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                        (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                            ghostcell_dim_1_cell_data;
                    
                    const int idx_cell_data_x_L = (i - 1 + num_ghosts_0_cell_data) +
                        (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                        (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                            ghostcell_dim_1_cell_data;
                    
                    const int idx_cell_data_x_R = (i + 1 + num_ghosts_0_cell_data) +
                        (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                        (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                            ghostcell_dim_1_cell_data;
                    
                    const int idx_cell_data_y_B = (i + num_ghosts_0_cell_data) +
                        (j - 1 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                        (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                            ghostcell_dim_1_cell_data;
                    
                    const int idx_cell_data_y_T = (i + num_ghosts_0_cell_data) +
                        (j + 1 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                        (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                            ghostcell_dim_1_cell_data;
                    
                    const int idx_cell_data_z_B = (i + num_ghosts_0_cell_data) +
                        (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                        (k - 1 + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                            ghostcell_dim_1_cell_data;
                    
                    const int idx_cell_data_z_F = (i + num_ghosts_0_cell_data) +
                        (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                        (k + 1 + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                            ghostcell_dim_1_cell_data;
                    
                    double w_x = f[idx_cell_data_x_R] - 2*f[idx_cell_data] + f[idx_cell_data_x_L];
                    double w_y = f[idx_cell_data_y_T] - 2*f[idx_cell_data] + f[idx_cell_data_y_B];
                    double w_z = f[idx_cell_data_z_F] - 2*f[idx_cell_data] + f[idx_cell_data_z_B];
                    
                    w[idx] = sqrt(w_x*w_x + w_y*w_y + w_z*w_z);
                }
            }
        }
        
        if (compute_variable_local_mean)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx = (i + num_ghosts_0_gradient) +
                            (j + num_ghosts_1_gradient)*ghostcell_dim_0_gradient +
                            (k + num_ghosts_2_gradient)*ghostcell_dim_0_gradient*
                                ghostcell_dim_1_gradient;
                        
                        const int idx_cell_data = (i + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                            (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                ghostcell_dim_1_cell_data;
                        
                        const int idx_cell_data_x_L = (i - 1 + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                            (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                ghostcell_dim_1_cell_data;
                        
                        const int idx_cell_data_x_R = (i + 1 + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                            (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                ghostcell_dim_1_cell_data;
                        
                        const int idx_cell_data_y_B = (i + num_ghosts_0_cell_data) +
                            (j - 1 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                            (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                ghostcell_dim_1_cell_data;
                        
                        const int idx_cell_data_y_T = (i + num_ghosts_0_cell_data) +
                            (j + 1 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                            (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                ghostcell_dim_1_cell_data;
                        
                        const int idx_cell_data_z_B = (i + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                            (k - 1 + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                ghostcell_dim_1_cell_data;
                        
                        const int idx_cell_data_z_F = (i + num_ghosts_0_cell_data) +
                            (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                            (k + 1 + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                                ghostcell_dim_1_cell_data;
                        
                        double f_mean_x = f[idx_cell_data_x_R] + 2*f[idx_cell_data] + f[idx_cell_data_x_L];
                        double f_mean_y = f[idx_cell_data_y_T] + 2*f[idx_cell_data] + f[idx_cell_data_y_B];
                        double f_mean_z = f[idx_cell_data_z_F] + 2*f[idx_cell_data] + f[idx_cell_data_z_B];
                        
                        f_mean[idx] = sqrt(f_mean_x*f_mean_x + f_mean_y*f_mean_y + f_mean_z*f_mean_z);
                    }
                }
            }
        }
    }
}
