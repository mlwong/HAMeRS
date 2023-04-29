#include "util/gradient_sensors/GradientSensorJameson.hpp"

#include "SAMRAI/geom/CartesianPatchGeometry.h"

#define EPSILON HAMERS_EPSILON

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
    HAMERS_SHARED_PTR<pdat::CellData<Real> >& gradient,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& cell_data,
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
    
    // Get the grid spacings.
    
    const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
        HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(patch_geom);
#endif
    
    const double* dx = patch_geom->getDx();
     
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
    Real* f = cell_data->getPointer(depth);
    
    // Get the pointer to the gradient.
    Real* psi = gradient->getPointer(0);
    
    if (d_dim == tbox::Dimension(1))
    {
        const int interior_dim_0 = interior_dims[0];
        
        const int num_ghosts_0_cell_data = num_ghosts_cell_data[0];
        const int num_ghosts_0_gradient = num_ghosts_gradient[0];
        
        HAMERS_PRAGMA_SIMD
        for (int i = 0; i < interior_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_gradient = i + num_ghosts_0_gradient;
            const int idx     = i + num_ghosts_0_cell_data;
            const int idx_x_L = i - 1 + num_ghosts_0_cell_data;
            const int idx_x_R = i + 1 + num_ghosts_0_cell_data;
            
            const Real psi_x  = f[idx_x_R] - Real(2)*f[idx] + f[idx_x_L];
            const Real mean_x = f[idx_x_R] + Real(2)*f[idx] + f[idx_x_L];
            psi[idx_gradient] = std::abs(psi_x)/(mean_x + Real(EPSILON));
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
        
        const Real sq_inv_dx_0 = Real(1)/Real(dx[0]*dx[0]);
        const Real sq_inv_dx_1 = Real(1)/Real(dx[1]*dx[1]);
        
        for (int j = 0; j < interior_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_gradient = (i + num_ghosts_0_gradient) +
                    (j + num_ghosts_1_gradient)*ghostcell_dim_0_gradient;
                
                const int idx = (i + num_ghosts_0_cell_data) +
                    (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                
                const int idx_x_L = (i - 1 + num_ghosts_0_cell_data) +
                    (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                
                const int idx_x_R = (i + 1 + num_ghosts_0_cell_data) +
                    (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                
                const int idx_y_B = (i + num_ghosts_0_cell_data) +
                    (j - 1 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                
                const int idx_y_T = (i + num_ghosts_0_cell_data) +
                    (j + 1 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data;
                
                const Real psi_x  = (f[idx_x_R] - Real(2)*f[idx] + f[idx_x_L])*sq_inv_dx_0;
                const Real mean_x = (f[idx_x_R] + Real(2)*f[idx] + f[idx_x_L])*sq_inv_dx_0;
                
                const Real psi_y  = (f[idx_y_T] - Real(2)*f[idx] + f[idx_y_B])*sq_inv_dx_1;
                const Real mean_y = (f[idx_y_T] + Real(2)*f[idx] + f[idx_y_B])*sq_inv_dx_1;
                
                psi[idx_gradient] = std::sqrt(psi_x*psi_x + psi_y*psi_y)/
                    (std::sqrt(mean_x*mean_x + mean_y*mean_y) + Real(EPSILON));
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
        
        const Real sq_inv_dx_0 = Real(1)/Real(dx[0]*dx[0]);
        const Real sq_inv_dx_1 = Real(1)/Real(dx[1]*dx[1]);
        const Real sq_inv_dx_2 = Real(1)/Real(dx[2]*dx[2]);
        
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_gradient = (i + num_ghosts_0_gradient) +
                        (j + num_ghosts_1_gradient)*ghostcell_dim_0_gradient +
                        (k + num_ghosts_2_gradient)*ghostcell_dim_0_gradient*
                            ghostcell_dim_1_gradient;
                    
                    const int idx = (i + num_ghosts_0_cell_data) +
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
                    
                    const int idx_y_B = (i + num_ghosts_0_cell_data) +
                        (j - 1 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                        (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                            ghostcell_dim_1_cell_data;
                    
                    const int idx_y_T = (i + num_ghosts_0_cell_data) +
                        (j + 1 + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                        (k + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                            ghostcell_dim_1_cell_data;
                    
                    const int idx_z_B = (i + num_ghosts_0_cell_data) +
                        (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                        (k - 1 + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                            ghostcell_dim_1_cell_data;
                    
                    const int idx_z_F = (i + num_ghosts_0_cell_data) +
                        (j + num_ghosts_1_cell_data)*ghostcell_dim_0_cell_data +
                        (k + 1 + num_ghosts_2_cell_data)*ghostcell_dim_0_cell_data*
                            ghostcell_dim_1_cell_data;
                    
                    const Real psi_x  = (f[idx_x_R] - Real(2)*f[idx] + f[idx_x_L])*sq_inv_dx_0;
                    const Real mean_x = (f[idx_x_R] + Real(2)*f[idx] + f[idx_x_L])*sq_inv_dx_0;
                    
                    const Real psi_y  = (f[idx_y_T] - Real(2)*f[idx] + f[idx_y_B])*sq_inv_dx_1;
                    const Real mean_y = (f[idx_y_T] + Real(2)*f[idx] + f[idx_y_B])*sq_inv_dx_1;
                    
                    const Real psi_z  = (f[idx_z_F] - Real(2)*f[idx] + f[idx_z_B])*sq_inv_dx_2;
                    const Real mean_z = (f[idx_z_F] + Real(2)*f[idx] + f[idx_z_B])*sq_inv_dx_2;
                    
                    psi[idx_gradient] = std::sqrt(psi_x*psi_x + psi_y*psi_y + psi_z*psi_z)/
                        (std::sqrt(mean_x*mean_x + mean_y*mean_y + mean_z*mean_z) +
                         Real(EPSILON));
                }
            }
        }
    }
}
