#include "util/gradient_sensors/GradientSensorDucros.hpp"

#include "SAMRAI/geom/CartesianPatchGeometry.h"

#define EPSILON HAMERS_EPSILON

GradientSensorDucros::GradientSensorDucros(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const bool use_dilatation,
    const bool use_strain_rate):
        GradientSensor(
            object_name,
            dim),
        d_use_dilatation(use_dilatation),
        d_use_strain_rate(use_strain_rate)
{
    d_num_gradient_ghosts = hier::IntVector::getOne(d_dim)*3;
    
    if ((!use_dilatation) && (!use_strain_rate))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Cannot have both 'use_dilatation' and 'use_strain_rate' to be false!"
            << std::endl);
    }
}


/*
 * Compute the gradient with the given cell data.
 */
void
GradientSensorDucros::computeGradient(
    HAMERS_SHARED_PTR<pdat::CellData<Real> >& gradient,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& cell_data_velocity,
    hier::Patch& patch,
    const int depth)
{
    NULL_USE(depth);
    
    if (cell_data_velocity->getGhostCellWidth() < d_num_gradient_ghosts)
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
    const hier::IntVector num_ghosts_cell_data_velocity = cell_data_velocity->getGhostCellWidth();
    const hier::IntVector num_ghosts_gradient = gradient->getGhostCellWidth();
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus ghost cells.
    const hier::Box ghost_box_cell_data_velocity = cell_data_velocity->getGhostBox();
    const hier::IntVector ghostcell_dims_cell_data_velocity = ghost_box_cell_data_velocity.numberCells();
    
    const hier::Box ghost_box_gradient = gradient->getGhostBox();
    const hier::IntVector ghostcell_dims_gradient = ghost_box_gradient.numberCells();
    
    // Get the pointer to the gradient.
    Real* psi = gradient->getPointer(0);
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Ducros sensor cannot be used for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        // Get the pointers to the given cell velocity data.
        Real* u = cell_data_velocity->getPointer(0);
        Real* v = cell_data_velocity->getPointer(1);
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        const int num_ghosts_0_cell_data_velocity = num_ghosts_cell_data_velocity[0];
        const int num_ghosts_1_cell_data_velocity = num_ghosts_cell_data_velocity[1];
        const int ghostcell_dim_0_cell_data_velocity = ghostcell_dims_cell_data_velocity[0];
        
        const int num_ghosts_0_gradient = num_ghosts_gradient[0];
        const int num_ghosts_1_gradient = num_ghosts_gradient[1];
        const int ghostcell_dim_0_gradient = ghostcell_dims_gradient[0];
        
        const Real half      = Real(1)/Real(2);
        const Real two       = Real(2);
        const Real one_fifth = Real(1)/(5);
        const Real inv_dx_0 = Real(1)/Real(dx[0]);
        const Real inv_dx_1 = Real(1)/Real(dx[1]);
        
        for (int j = 0; j < interior_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_gradient = (i + num_ghosts_0_gradient) +
                    (j + num_ghosts_1_gradient)*ghostcell_dim_0_gradient;
                
                const int idx = (i + num_ghosts_0_cell_data_velocity) +
                    (j + num_ghosts_1_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity;
                
                const int idx_x_LLL = (i - 3 + num_ghosts_0_cell_data_velocity) +
                    (j + num_ghosts_1_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity;
                
                const int idx_x_LL = (i - 2 + num_ghosts_0_cell_data_velocity) +
                    (j + num_ghosts_1_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity;
                
                const int idx_x_L = (i - 1 + num_ghosts_0_cell_data_velocity) +
                    (j + num_ghosts_1_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity;
                
                const int idx_x_R = (i + 1 + num_ghosts_0_cell_data_velocity) +
                    (j + num_ghosts_1_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity;
                
                const int idx_x_RR = (i + 2 + num_ghosts_0_cell_data_velocity) +
                    (j + num_ghosts_1_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity;
                
                const int idx_x_RRR = (i + 3 + num_ghosts_0_cell_data_velocity) +
                    (j + num_ghosts_1_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity;
                
                const int idx_y_BBB = (i + num_ghosts_0_cell_data_velocity) +
                    (j - 3 + num_ghosts_1_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity;
                
                const int idx_y_BB = (i + num_ghosts_0_cell_data_velocity) +
                    (j - 2 + num_ghosts_1_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity;
                
                const int idx_y_B = (i + num_ghosts_0_cell_data_velocity) +
                    (j - 1 + num_ghosts_1_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity;
                
                const int idx_y_T = (i + num_ghosts_0_cell_data_velocity) +
                    (j + 1 + num_ghosts_1_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity;
                
                const int idx_y_TT = (i + num_ghosts_0_cell_data_velocity) +
                    (j + 2 + num_ghosts_1_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity;
                
                const int idx_y_TTT = (i + num_ghosts_0_cell_data_velocity) +
                    (j + 3 + num_ghosts_1_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity;
                
                // Compute the derivatives in the x-direction.
                const Real du_dx_LL = half*(u[idx_x_L] - u[idx_x_LLL])*inv_dx_0;
                const Real dv_dx_LL = half*(v[idx_x_L] - v[idx_x_LLL])*inv_dx_0;
                
                const Real du_dx_L = half*(u[idx] - u[idx_x_LL])*inv_dx_0;
                const Real dv_dx_L = half*(v[idx] - v[idx_x_LL])*inv_dx_0;
                
                const Real du_dx_M = half*(u[idx_x_R] - u[idx_x_L])*inv_dx_0;
                const Real dv_dx_M = half*(v[idx_x_R] - v[idx_x_L])*inv_dx_0;
                
                const Real du_dx_R = half*(u[idx_x_RR] - u[idx])*inv_dx_0;
                const Real dv_dx_R = half*(v[idx_x_RR] - v[idx])*inv_dx_0;
                
                const Real du_dx_RR = half*(u[idx_x_RRR] - u[idx_x_R])*inv_dx_0;
                const Real dv_dx_RR = half*(v[idx_x_RRR] - v[idx_x_R])*inv_dx_0;
                
                // Compute the derivatives in the y-direction.
                const Real du_dy_BB = half*(u[idx_y_B] - u[idx_y_BBB])*inv_dx_1;
                const Real dv_dy_BB = half*(v[idx_y_B] - v[idx_y_BBB])*inv_dx_1;
                
                const Real du_dy_B = half*(u[idx] - u[idx_y_BB])*inv_dx_1;
                const Real dv_dy_B = half*(v[idx] - v[idx_y_BB])*inv_dx_1;
                
                const Real du_dy_M = half*(u[idx_y_T] - u[idx_y_B])*inv_dx_1;
                const Real dv_dy_M = half*(v[idx_y_T] - v[idx_y_B])*inv_dx_1;
                
                const Real du_dy_T = half*(u[idx_y_TT] - u[idx])*inv_dx_1;
                const Real dv_dy_T = half*(v[idx_y_TT] - v[idx])*inv_dx_1;
                
                const Real du_dy_TT = half*(u[idx_y_TTT] - u[idx_y_T])*inv_dx_1;
                const Real dv_dy_TT = half*(v[idx_y_TTT] - v[idx_y_T])*inv_dx_1;
                
                // Filter the derivatives.
                const Real du_dx = one_fifth*(du_dx_LL + du_dx_L + du_dx_M + du_dx_R + du_dx_RR);
                const Real dv_dx = one_fifth*(dv_dx_LL + dv_dx_L + dv_dx_M + dv_dx_R + dv_dx_RR);
                
                const Real du_dy = one_fifth*(du_dy_BB + du_dy_B + du_dy_M + du_dy_T + du_dy_TT);
                const Real dv_dy = one_fifth*(dv_dy_BB + dv_dy_B + dv_dy_M + dv_dy_T + dv_dy_TT);
                
                const Real S_11 = du_dx;
                const Real S_12 = half*(du_dy + dv_dx);
                const Real S_22 = dv_dy;
                
                const Real S_sq = S_11*S_11 + S_22*S_22 + two*S_12*S_12;
                const Real theta = du_dx + dv_dy;
                const Real omega = dv_dx - du_dy;
                const Real Omega = omega*omega;
                
                psi[idx_gradient] = Real(0);
                if (d_use_dilatation && theta < Real(0))
                {
                    psi[idx_gradient] = std::max(psi[idx_gradient],
                        theta*theta/(theta*theta + Omega + Real(EPSILON)));
                }
                
                if (d_use_strain_rate)
                {
                    psi[idx_gradient] = std::max(psi[idx_gradient],
                        Omega/(S_sq + Omega + Real(EPSILON)));
                }
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        // Get the pointers to the given cell velocity data.
        Real* u = cell_data_velocity->getPointer(0);
        Real* v = cell_data_velocity->getPointer(1);
        Real* w = cell_data_velocity->getPointer(2);
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        const int interior_dim_2 = interior_dims[2];
        
        const int num_ghosts_0_cell_data_velocity = num_ghosts_cell_data_velocity[0];
        const int num_ghosts_1_cell_data_velocity = num_ghosts_cell_data_velocity[1];
        const int num_ghosts_2_cell_data_velocity = num_ghosts_cell_data_velocity[2];
        const int ghostcell_dim_0_cell_data_velocity = ghostcell_dims_cell_data_velocity[0];
        const int ghostcell_dim_1_cell_data_velocity = ghostcell_dims_cell_data_velocity[1];
        
        const int num_ghosts_0_gradient = num_ghosts_gradient[0];
        const int num_ghosts_1_gradient = num_ghosts_gradient[1];
        const int num_ghosts_2_gradient = num_ghosts_gradient[2];
        const int ghostcell_dim_0_gradient = ghostcell_dims_gradient[0];
        const int ghostcell_dim_1_gradient = ghostcell_dims_gradient[1];
        
        const Real half      = Real(1)/Real(2);
        const Real two       = Real(2);
        const Real one_fifth = Real(1)/(5);
        const Real inv_dx_0 = Real(1)/Real(dx[0]);
        const Real inv_dx_1 = Real(1)/Real(dx[1]);
        const Real inv_dx_2 = Real(1)/Real(dx[2]);
        
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
                    
                    const int idx = (i + num_ghosts_0_cell_data_velocity) +
                        (j + num_ghosts_1_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity +
                        (k + num_ghosts_2_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity*
                            ghostcell_dim_1_cell_data_velocity;
                    
                    const int idx_x_LLL = (i - 3 + num_ghosts_0_cell_data_velocity) +
                        (j + num_ghosts_1_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity +
                        (k + num_ghosts_2_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity*
                            ghostcell_dim_1_cell_data_velocity;
                    
                    const int idx_x_LL = (i - 2 + num_ghosts_0_cell_data_velocity) +
                        (j + num_ghosts_1_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity +
                        (k + num_ghosts_2_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity*
                            ghostcell_dim_1_cell_data_velocity;
                    
                    const int idx_x_L = (i - 1 + num_ghosts_0_cell_data_velocity) +
                        (j + num_ghosts_1_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity +
                        (k + num_ghosts_2_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity*
                            ghostcell_dim_1_cell_data_velocity;
                    
                    const int idx_x_R = (i + 1 + num_ghosts_0_cell_data_velocity) +
                        (j + num_ghosts_1_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity +
                        (k + num_ghosts_2_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity*
                            ghostcell_dim_1_cell_data_velocity;
                    
                    const int idx_x_RR = (i + 2 + num_ghosts_0_cell_data_velocity) +
                        (j + num_ghosts_1_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity +
                        (k + num_ghosts_2_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity*
                            ghostcell_dim_1_cell_data_velocity;
                    
                    const int idx_x_RRR = (i + 3 + num_ghosts_0_cell_data_velocity) +
                        (j + num_ghosts_1_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity +
                        (k + num_ghosts_2_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity*
                            ghostcell_dim_1_cell_data_velocity;
                    
                    const int idx_y_BBB = (i + num_ghosts_0_cell_data_velocity) +
                        (j - 3 + num_ghosts_1_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity +
                        (k + num_ghosts_2_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity*
                            ghostcell_dim_1_cell_data_velocity;
                    
                    const int idx_y_BB = (i + num_ghosts_0_cell_data_velocity) +
                        (j - 2 + num_ghosts_1_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity +
                        (k + num_ghosts_2_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity*
                            ghostcell_dim_1_cell_data_velocity;
                    
                    const int idx_y_B = (i + num_ghosts_0_cell_data_velocity) +
                        (j - 1 + num_ghosts_1_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity +
                        (k + num_ghosts_2_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity*
                            ghostcell_dim_1_cell_data_velocity;
                    
                    const int idx_y_T = (i + num_ghosts_0_cell_data_velocity) +
                        (j + 1 + num_ghosts_1_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity +
                        (k + num_ghosts_2_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity*
                            ghostcell_dim_1_cell_data_velocity;
                    
                    const int idx_y_TT = (i + num_ghosts_0_cell_data_velocity) +
                        (j + 2 + num_ghosts_1_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity +
                        (k + num_ghosts_2_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity*
                            ghostcell_dim_1_cell_data_velocity;
                    
                    const int idx_y_TTT = (i + num_ghosts_0_cell_data_velocity) +
                        (j + 3 + num_ghosts_1_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity +
                        (k + num_ghosts_2_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity*
                            ghostcell_dim_1_cell_data_velocity;
                    
                    const int idx_z_BBB = (i + num_ghosts_0_cell_data_velocity) +
                        (j + num_ghosts_1_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity +
                        (k - 3 + num_ghosts_2_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity*
                            ghostcell_dim_1_cell_data_velocity;
                    
                    const int idx_z_BB = (i + num_ghosts_0_cell_data_velocity) +
                        (j + num_ghosts_1_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity +
                        (k - 2 + num_ghosts_2_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity*
                            ghostcell_dim_1_cell_data_velocity;
                    
                    const int idx_z_B = (i + num_ghosts_0_cell_data_velocity) +
                        (j + num_ghosts_1_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity +
                        (k - 1 + num_ghosts_2_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity*
                            ghostcell_dim_1_cell_data_velocity;
                    
                    const int idx_z_F = (i + num_ghosts_0_cell_data_velocity) +
                        (j + num_ghosts_1_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity +
                        (k + 1 + num_ghosts_2_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity*
                            ghostcell_dim_1_cell_data_velocity;
                    
                    const int idx_z_FF = (i + num_ghosts_0_cell_data_velocity) +
                        (j + num_ghosts_1_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity +
                        (k + 2 + num_ghosts_2_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity*
                            ghostcell_dim_1_cell_data_velocity;
                    
                    const int idx_z_FFF = (i + num_ghosts_0_cell_data_velocity) +
                        (j + num_ghosts_1_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity +
                        (k + 3 + num_ghosts_2_cell_data_velocity)*ghostcell_dim_0_cell_data_velocity*
                            ghostcell_dim_1_cell_data_velocity;
                    
                    // Compute the derivatives in the x-direction.
                    const Real du_dx_LL = half*(u[idx_x_L] - u[idx_x_LLL])*inv_dx_0;
                    const Real dv_dx_LL = half*(v[idx_x_L] - v[idx_x_LLL])*inv_dx_0;
                    const Real dw_dx_LL = half*(w[idx_x_L] - w[idx_x_LLL])*inv_dx_0;
                    
                    const Real du_dx_L = half*(u[idx] - u[idx_x_LL])*inv_dx_0;
                    const Real dv_dx_L = half*(v[idx] - v[idx_x_LL])*inv_dx_0;
                    const Real dw_dx_L = half*(w[idx] - w[idx_x_LL])*inv_dx_0;
                    
                    const Real du_dx_M = half*(u[idx_x_R] - u[idx_x_L])*inv_dx_0;
                    const Real dv_dx_M = half*(v[idx_x_R] - v[idx_x_L])*inv_dx_0;
                    const Real dw_dx_M = half*(w[idx_x_R] - w[idx_x_L])*inv_dx_0;
                    
                    const Real du_dx_R = half*(u[idx_x_RR] - u[idx])*inv_dx_0;
                    const Real dv_dx_R = half*(v[idx_x_RR] - v[idx])*inv_dx_0;
                    const Real dw_dx_R = half*(w[idx_x_RR] - w[idx])*inv_dx_0;
                    
                    const Real du_dx_RR = half*(u[idx_x_RRR] - u[idx_x_R])*inv_dx_0;
                    const Real dv_dx_RR = half*(v[idx_x_RRR] - v[idx_x_R])*inv_dx_0;
                    const Real dw_dx_RR = half*(w[idx_x_RRR] - w[idx_x_R])*inv_dx_0;
                    
                    // Compute the derivatives in the y-direction.
                    const Real du_dy_BB = half*(u[idx_y_B] - u[idx_y_BBB])*inv_dx_1;
                    const Real dv_dy_BB = half*(v[idx_y_B] - v[idx_y_BBB])*inv_dx_1;
                    const Real dw_dy_BB = half*(w[idx_y_B] - w[idx_y_BBB])*inv_dx_1;
                    
                    const Real du_dy_B = half*(u[idx] - u[idx_y_BB])*inv_dx_1;
                    const Real dv_dy_B = half*(v[idx] - v[idx_y_BB])*inv_dx_1;
                    const Real dw_dy_B = half*(w[idx] - w[idx_y_BB])*inv_dx_1;
                    
                    const Real du_dy_M = half*(u[idx_y_T] - u[idx_y_B])*inv_dx_1;
                    const Real dv_dy_M = half*(v[idx_y_T] - v[idx_y_B])*inv_dx_1;
                    const Real dw_dy_M = half*(w[idx_y_T] - w[idx_y_B])*inv_dx_1;
                    
                    const Real du_dy_T = half*(u[idx_y_TT] - u[idx])*inv_dx_1;
                    const Real dv_dy_T = half*(v[idx_y_TT] - v[idx])*inv_dx_1;
                    const Real dw_dy_T = half*(w[idx_y_TT] - w[idx])*inv_dx_1;
                    
                    const Real du_dy_TT = half*(u[idx_y_TTT] - u[idx_y_T])*inv_dx_1;
                    const Real dv_dy_TT = half*(v[idx_y_TTT] - v[idx_y_T])*inv_dx_1;
                    const Real dw_dy_TT = half*(w[idx_y_TTT] - w[idx_y_T])*inv_dx_1;
                    
                    // Compute the derivatives in the z-direction.
                    const Real du_dz_BB = half*(u[idx_z_B] - u[idx_z_BBB])*inv_dx_2;
                    const Real dv_dz_BB = half*(v[idx_z_B] - v[idx_z_BBB])*inv_dx_2;
                    const Real dw_dz_BB = half*(w[idx_z_B] - w[idx_z_BBB])*inv_dx_2;
                    
                    const Real du_dz_B = half*(u[idx] - u[idx_z_BB])*inv_dx_2;
                    const Real dv_dz_B = half*(v[idx] - v[idx_z_BB])*inv_dx_2;
                    const Real dw_dz_B = half*(w[idx] - w[idx_z_BB])*inv_dx_2;
                    
                    const Real du_dz_M = half*(u[idx_z_F] - u[idx_z_B])*inv_dx_2;
                    const Real dv_dz_M = half*(v[idx_z_F] - v[idx_z_B])*inv_dx_2;
                    const Real dw_dz_M = half*(w[idx_z_F] - w[idx_z_B])*inv_dx_2;
                    
                    const Real du_dz_F = half*(u[idx_z_FF] - u[idx])*inv_dx_2;
                    const Real dv_dz_F = half*(v[idx_z_FF] - v[idx])*inv_dx_2;
                    const Real dw_dz_F = half*(w[idx_z_FF] - w[idx])*inv_dx_2;
                    
                    const Real du_dz_FF = half*(u[idx_z_FFF] - u[idx_z_F])*inv_dx_2;
                    const Real dv_dz_FF = half*(v[idx_z_FFF] - v[idx_z_F])*inv_dx_2;
                    const Real dw_dz_FF = half*(w[idx_z_FFF] - w[idx_z_F])*inv_dx_2;
                    
                    // Filter the derivatives.
                    const Real du_dx = one_fifth*(du_dx_LL + du_dx_L + du_dx_M + du_dx_R + du_dx_RR);
                    const Real dv_dx = one_fifth*(dv_dx_LL + dv_dx_L + dv_dx_M + dv_dx_R + dv_dx_RR);
                    const Real dw_dx = one_fifth*(dw_dx_LL + dw_dx_L + dw_dx_M + dw_dx_R + dw_dx_RR);
                    
                    const Real du_dy = one_fifth*(du_dy_BB + du_dy_B + du_dy_M + du_dy_T + du_dy_TT);
                    const Real dv_dy = one_fifth*(dv_dy_BB + dv_dy_B + dv_dy_M + dv_dy_T + dv_dy_TT);
                    const Real dw_dy = one_fifth*(dw_dy_BB + dw_dy_B + dw_dy_M + dw_dy_T + dw_dy_TT);
                    
                    const Real du_dz = one_fifth*(du_dz_BB + du_dz_B + du_dz_M + du_dz_F + du_dz_FF);
                    const Real dv_dz = one_fifth*(dv_dz_BB + dv_dz_B + dv_dz_M + dv_dz_F + dv_dz_FF);
                    const Real dw_dz = one_fifth*(dw_dz_BB + dw_dz_B + dw_dz_M + dw_dz_F + dw_dz_FF);
                    
                    const Real S_11 = du_dx;
                    const Real S_12 = half*(du_dy + dv_dx);
                    const Real S_13 = half*(du_dz + dw_dx);
                    const Real S_22 = dv_dy;
                    const Real S_23 = half*(dv_dz + dw_dy);
                    const Real S_33 = dw_dz;
                    
                    const Real S_sq = S_11*S_11 + S_22*S_22 + S_33*S_33 + two*(S_12*S_12 + S_13*S_13 + S_23*S_23);
                    const Real theta = du_dx + dv_dy + dw_dz;
                    const Real omega_x = dw_dy - dv_dz;
                    const Real omega_y = du_dz - dw_dx;
                    const Real omega_z = dv_dx - du_dy;
                    const Real Omega = omega_x*omega_x + omega_y*omega_y + omega_z*omega_z;
                    
                    psi[idx_gradient] = Real(0);
                    if (d_use_dilatation && theta < Real(0))
                    {
                        psi[idx_gradient] = std::max(psi[idx_gradient],
                            theta*theta/(theta*theta + Omega + Real(EPSILON)));
                    }
                    
                    if (d_use_strain_rate)
                    {
                        psi[idx_gradient] = std::max(psi[idx_gradient],
                            Omega/(S_sq + Omega + Real(EPSILON)));
                    }
                }
            }
        }
    }
}
