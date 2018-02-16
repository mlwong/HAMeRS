#include "util/gradient_sensors/GradientSensorJameson.hpp"

#include "SAMRAI/geom/CartesianPatchGeometry.h"

#include <cfloat>

#define EPSILON DBL_EPSILON

GradientSensorJameson::GradientSensorJameson(
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
GradientSensorJameson::computeGradient(
    boost::shared_ptr<pdat::CellData<double> >& gradient,
    const boost::shared_ptr<pdat::CellData<double> >& cell_data,
    hier::Patch& patch,
    const int depth)
{
    if (cell_data->getGhostCellWidth() < d_num_gradient_ghosts)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The ghost cell width is smaller than required."
            << std::endl);
    }
    
    // Get the number of ghost cells of the cell data and gradient data.
    const hier::IntVector num_ghosts_cell_data = cell_data->getGhostCellWidth();
    const hier::IntVector num_ghosts_gradient = gradient->getGhostCellWidth();
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus ghost cells.
    const hier::Box ghost_box_cell_data = cell_data->getGhostBox();
    const hier::IntVector ghostcell_dims_cell_data = ghost_box_cell_data.numberCells();
    
    const hier::Box ghost_box_gradient = gradient->getGhostBox();
    const hier::IntVector ghostcell_dims_gradient = ghost_box_gradient.numberCells();
    
    // Get the pointer to the current depth component of the given cell data.
    double* f = cell_data->getPointer(depth);
    
    // Get the pointer to the gradient.
    double* psi = gradient->getPointer(0);
    
    if (d_dim == tbox::Dimension(1))
    {
        const int interior_dim_0 = interior_dims[0];
        
        const int num_ghosts_0_cell_data = num_ghosts_cell_data[0];
        const int num_ghosts_0_gradient = num_ghosts_gradient[0];
        
        // Allocate memory.
        boost::shared_ptr<pdat::CellData<double> > gradient_x(
            new pdat::CellData<double>(interior_box, 1, num_ghosts_gradient));
        boost::shared_ptr<pdat::CellData<double> > local_mean_value_x(
            new pdat::CellData<double>(interior_box, 1, num_ghosts_gradient));
        
        double* psi_x = gradient_x->getPointer(0);
        double* mean_x = local_mean_value_x->getPointer(0);
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = 0; i < interior_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx = i + num_ghosts_0_gradient;
            const int idx_x_L = i - 1 + num_ghosts_0_cell_data;
            const int idx_x   = i + num_ghosts_0_cell_data;
            const int idx_x_R = i + 1 + num_ghosts_0_cell_data;
            
            psi_x[idx] = f[idx_x_R] - 2.0*f[idx_x] + f[idx_x_L];
            mean_x[idx] = f[idx_x_R] + 2.0*f[idx_x] + f[idx_x_L];
        }
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = 0; i < interior_dim_0; i++)
        {
            // Compute the linear index.
            const int idx = i + num_ghosts_0_gradient;
            
            psi[idx] = fabs(psi_x[idx])/(mean_x[idx] + EPSILON);
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
        
        // Allocate memory in different dimensions.
        boost::shared_ptr<pdat::CellData<double> > gradient_x(
            new pdat::CellData<double>(interior_box, 1, num_ghosts_gradient));
        boost::shared_ptr<pdat::CellData<double> > gradient_y(
            new pdat::CellData<double>(interior_box, 1, num_ghosts_gradient));
        boost::shared_ptr<pdat::CellData<double> > local_mean_value_x(
            new pdat::CellData<double>(interior_box, 1, num_ghosts_gradient));
        boost::shared_ptr<pdat::CellData<double> > local_mean_value_y(
            new pdat::CellData<double>(interior_box, 1, num_ghosts_gradient));
        
        double* psi_x = gradient_x->getPointer(0);
        double* psi_y = gradient_y->getPointer(0);
        double* mean_x = local_mean_value_x->getPointer(0);
        double* mean_y = local_mean_value_y->getPointer(0);
        
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
                
                const int idx_x_L = (i - 1 + num_ghosts_0_cell_data) +
                    (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                
                const int idx_x = (i + num_ghosts_0_cell_data) +
                    (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                
                const int idx_x_R = (i + 1 + num_ghosts_0_cell_data) +
                    (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                
                psi_x[idx] = f[idx_x_R] - 2.0*f[idx_x] + f[idx_x_L];
                mean_x[idx] = f[idx_x_R] + 2.0*f[idx_x] + f[idx_x_L];
            }
        }
        
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
                
                const int idx_y_B = (i + num_ghosts_0_cell_data) +
                    (j - 1 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                
                const int idx_y = (i + num_ghosts_0_cell_data) +
                    (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                
                const int idx_y_T = (i + num_ghosts_0_cell_data) +
                    (j + 1 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                
                psi_y[idx] = f[idx_y_T] - 2.0*f[idx_y] + f[idx_y_B];
                mean_y[idx] = f[idx_y_T] + 2.0*f[idx_y] + f[idx_y_B];
            }
        }
        
        for (int j = 0; j < interior_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear index.
                const int idx = (i + num_ghosts_0_gradient) +
                    (j + num_ghosts_1_gradient)*ghostcell_dim_0_gradient;
                
                psi[idx] = sqrt(psi_x[idx]*psi_x[idx] + psi_y[idx]*psi_y[idx])/
                    (sqrt(mean_x[idx]*mean_x[idx] + mean_y[idx]*mean_y[idx]) + EPSILON);
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
        
        // Allocate memory in different dimensions.
        boost::shared_ptr<pdat::CellData<double> > gradient_x(
            new pdat::CellData<double>(interior_box, 1, num_ghosts_gradient));
        boost::shared_ptr<pdat::CellData<double> > gradient_y(
            new pdat::CellData<double>(interior_box, 1, num_ghosts_gradient));
        boost::shared_ptr<pdat::CellData<double> > gradient_z(
            new pdat::CellData<double>(interior_box, 1, num_ghosts_gradient));
        boost::shared_ptr<pdat::CellData<double> > local_mean_value_x(
            new pdat::CellData<double>(interior_box, 1, num_ghosts_gradient));
        boost::shared_ptr<pdat::CellData<double> > local_mean_value_y(
            new pdat::CellData<double>(interior_box, 1, num_ghosts_gradient));
        boost::shared_ptr<pdat::CellData<double> > local_mean_value_z(
            new pdat::CellData<double>(interior_box, 1, num_ghosts_gradient));
        
        double* psi_x = gradient_x->getPointer(0);
        double* psi_y = gradient_y->getPointer(0);
        double* psi_z = gradient_z->getPointer(0);
        double* mean_x = local_mean_value_x->getPointer(0);
        double* mean_y = local_mean_value_y->getPointer(0);
        double* mean_z = local_mean_value_z->getPointer(0);
        
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
                    
                    const int idx_x_L = (i - 1 + num_ghosts_0_cell_data) +
                        (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                        (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                            ghostcell_dim_1_cell_data;
                    
                    const int idx_x = (i + num_ghosts_0_cell_data) +
                        (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                        (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                            ghostcell_dim_1_cell_data;
                    
                    const int idx_x_R = (i + 1 + num_ghosts_0_cell_data) +
                        (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                        (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                            ghostcell_dim_1_cell_data;
                    
                    psi_x[idx] = f[idx_x_R] - 2.0*f[idx_x] + f[idx_x_L];
                    mean_x[idx] = f[idx_x_R] + 2.0*f[idx_x] + f[idx_x_L];
                }
            }
        }
        
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
                    
                    const int idx_y_B = (i + num_ghosts_0_cell_data) +
                        (j - 1 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                        (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                            ghostcell_dim_1_cell_data;
                    
                    const int idx_y = (i + num_ghosts_0_cell_data) +
                        (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                        (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                            ghostcell_dim_1_cell_data;
                    
                    const int idx_y_T = (i + num_ghosts_0_cell_data) +
                        (j + 1 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                        (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                            ghostcell_dim_1_cell_data;
                    
                    psi_y[idx] = f[idx_y_T] - 2.0*f[idx_y] + f[idx_y_B];
                    mean_y[idx] = f[idx_y_T] + 2.0*f[idx_y] + f[idx_y_B];
                }
            }
        }
        
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
                    
                    const int idx_z_B = (i + num_ghosts_0_cell_data) +
                        (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                        (k - 1 + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                            ghostcell_dim_1_cell_data;
                    
                    const int idx_z = (i + num_ghosts_0_cell_data) +
                        (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                        (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                            ghostcell_dim_1_cell_data;
                    
                    const int idx_z_F = (i + num_ghosts_0_cell_data) +
                        (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                        (k + 1 + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                            ghostcell_dim_1_cell_data;
                    
                    psi_z[idx] = f[idx_z_F] - 2.0*f[idx_z] + f[idx_z_B];
                    mean_z[idx] = f[idx_z_F] + 2.0*f[idx_z] + f[idx_z_B];
                }
            }
        }
        
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the index.
                    const int idx = (i + num_ghosts_0_gradient) +
                        (j + num_ghosts_1_gradient)*ghostcell_dim_0_gradient +
                        (k + num_ghosts_2_gradient)*ghostcell_dim_0_gradient*
                            ghostcell_dim_1_gradient;
                    
                    psi[idx] = sqrt(psi_x[idx]*psi_x[idx] + psi_y[idx]*psi_y[idx] + psi_z[idx]*psi_z[idx])/
                        (sqrt(mean_x[idx]*mean_x[idx] + mean_y[idx]*mean_y[idx] + mean_z[idx]*mean_z[idx]) +
                         EPSILON);
                }
            }
        }
    }
}
