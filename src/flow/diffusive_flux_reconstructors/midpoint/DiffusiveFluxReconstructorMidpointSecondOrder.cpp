#include "flow/diffusive_flux_reconstructors/midpoint/DiffusiveFluxReconstructorMidpointSecondOrder.hpp"

#include "SAMRAI/geom/CartesianPatchGeometry.h"

#include <map>

DiffusiveFluxReconstructorMidpointSecondOrder::DiffusiveFluxReconstructorMidpointSecondOrder(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_eqn,
    const HAMERS_SHARED_PTR<FlowModel>& flow_model,
    const HAMERS_SHARED_PTR<tbox::Database>& diffusive_flux_reconstructor_db):
        DiffusiveFluxReconstructorMidpoint(
            object_name,
            dim,
            grid_geometry,
            num_eqn,
            flow_model,
            diffusive_flux_reconstructor_db)
{
    d_num_diff_ghosts = hier::IntVector::getOne(d_dim);
    
    // The number of ghost cells used for finite difference and interpolation schemes at midpoints are the same.
    d_num_der_midpoint_ghosts     = hier::IntVector::getOne(d_dim);
    d_num_der_node_ghosts         = hier::IntVector::getOne(d_dim);
    d_num_interp_midpoint_ghosts  = hier::IntVector::getOne(d_dim);
    d_num_flux_reconstruct_ghosts = hier::IntVector::getOne(d_dim);
}


/*
 * Print all characteristics of the diffusive flux reconstruction class.
 */
void
DiffusiveFluxReconstructorMidpointSecondOrder::printClassData(
    std::ostream& os) const
{
    os << "\nPrint DiffusiveFluxReconstructorMidpointSecondOrder object..."
       << std::endl;
    
    os << std::endl;
    
    os << "DiffusiveFluxReconstructorMidpointSecondOrder: this = "
       << (DiffusiveFluxReconstructorMidpointSecondOrder *)this
       << std::endl;
    os << "d_object_name = "
       << d_object_name
       << std::endl;
}


/*
 * Put the characteristics of the diffusive flux reconstruction class
 * into the restart database.
 */
void
DiffusiveFluxReconstructorMidpointSecondOrder::putToRestart(
   const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const
{
    restart_db->putString("d_diffusive_flux_reconstructor", "MIDPOINT_SECOND_ORDER");
}


/*
 * Kernel to compute the derivatives in x-direction at midpoints.
 */
void
DiffusiveFluxReconstructorMidpointSecondOrder::computeFirstDerivativesInXAtMidpointX(
    double* dudx,
    const double* const u,
    const hier::IntVector& num_ghosts_derivative_midpoint,
    const hier::IntVector& num_ghosts_data_node,
    const hier::IntVector& ghostcell_dims_derivative_midpoint,
    const hier::IntVector& ghostcell_dims_data_node,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims,
    const double& dx_0_inv) const
{
    const double a_n = double(1);
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and number of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int num_ghosts_0_derivative_midpoint = num_ghosts_derivative_midpoint[0];
        
        const int num_ghosts_0_data = num_ghosts_data_node[0];
        
        HAMERS_PRAGMA_SIMD
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx = i + num_ghosts_0_derivative_midpoint;
            
            const int idx_data_L = i - 1 + num_ghosts_0_data;
            const int idx_data_R = i + 0 + num_ghosts_0_data;
            
            dudx[idx] = a_n*(u[idx_data_R] - u[idx_data_L])*dx_0_inv;
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and number of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        
        const int num_ghosts_0_derivative_midpoint = num_ghosts_derivative_midpoint[0];
        const int num_ghosts_1_derivative_midpoint = num_ghosts_derivative_midpoint[1];
        const int ghostcell_dim_0_derivative_midpoint = ghostcell_dims_derivative_midpoint[0];
        
        const int num_ghosts_0_data = num_ghosts_data_node[0];
        const int num_ghosts_1_data = num_ghosts_data_node[1];
        const int ghostcell_dim_0_data = ghostcell_dims_data_node[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx = (i + num_ghosts_0_derivative_midpoint) +
                    (j + num_ghosts_1_derivative_midpoint)*(ghostcell_dim_0_derivative_midpoint + 1);
                
                const int idx_data_L = (i - 1 + num_ghosts_0_data) +
                    (j + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                const int idx_data_R = (i + num_ghosts_0_data) +
                    (j + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                dudx[idx] = a_n*(u[idx_data_R] - u[idx_data_L])*dx_0_inv;
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and number of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_lo_2 = domain_lo[2];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        const int domain_dim_2 = domain_dims[2];
        
        const int num_ghosts_0_derivative_midpoint = num_ghosts_derivative_midpoint[0];
        const int num_ghosts_1_derivative_midpoint = num_ghosts_derivative_midpoint[1];
        const int num_ghosts_2_derivative_midpoint = num_ghosts_derivative_midpoint[2];
        const int ghostcell_dim_0_derivative_midpoint = ghostcell_dims_derivative_midpoint[0];
        const int ghostcell_dim_1_derivative_midpoint = ghostcell_dims_derivative_midpoint[1];
        
        const int num_ghosts_0_data = num_ghosts_data_node[0];
        const int num_ghosts_1_data = num_ghosts_data_node[1];
        const int num_ghosts_2_data = num_ghosts_data_node[2];
        const int ghostcell_dim_0_data = ghostcell_dims_data_node[0];
        const int ghostcell_dim_1_data = ghostcell_dims_data_node[1];
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx = (i + num_ghosts_0_derivative_midpoint) +
                        (j + num_ghosts_1_derivative_midpoint)*(ghostcell_dim_0_derivative_midpoint + 1) +
                        (k + num_ghosts_2_derivative_midpoint)*(ghostcell_dim_0_derivative_midpoint + 1)*
                            ghostcell_dim_1_derivative_midpoint;
                    
                    const int idx_data_L = (i - 1 + num_ghosts_0_data) +
                        (j + num_ghosts_1_data)*ghostcell_dim_0_data +
                        (k + num_ghosts_2_data)*ghostcell_dim_0_data*
                            ghostcell_dim_1_data;
                    
                    const int idx_data_R = (i + num_ghosts_0_data) +
                        (j + num_ghosts_1_data)*ghostcell_dim_0_data +
                        (k + num_ghosts_2_data)*ghostcell_dim_0_data*
                            ghostcell_dim_1_data;
                    
                    dudx[idx] = a_n*(u[idx_data_R] - u[idx_data_L])*dx_0_inv;
                }
            }
        }
    }
}


/*
 * Kernel to compute the derivatives in y-direction at midpoints.
 */
void
DiffusiveFluxReconstructorMidpointSecondOrder::computeFirstDerivativesInYAtMidpointY(
    double* dudy,
    const double* const u,
    const hier::IntVector& num_ghosts_derivative_midpoint,
    const hier::IntVector& num_ghosts_data_node,
    const hier::IntVector& ghostcell_dims_derivative_midpoint,
    const hier::IntVector& ghostcell_dims_data_node,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims,
    const double& dx_1_inv) const
{
    const double a_n = double(1);
    
    if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and number of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        
        const int num_ghosts_0_derivative_midpoint = num_ghosts_derivative_midpoint[0];
        const int num_ghosts_1_derivative_midpoint = num_ghosts_derivative_midpoint[1];
        const int ghostcell_dim_0_derivative_midpoint = ghostcell_dims_derivative_midpoint[0];
        
        const int num_ghosts_0_data = num_ghosts_data_node[0];
        const int num_ghosts_1_data = num_ghosts_data_node[1];
        const int ghostcell_dim_0_data = ghostcell_dims_data_node[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx = (i + num_ghosts_0_derivative_midpoint) +
                    (j + num_ghosts_1_derivative_midpoint)*ghostcell_dim_0_derivative_midpoint;
                
                const int idx_data_B = (i + num_ghosts_0_data) +
                    (j - 1 + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                const int idx_data_T = (i + num_ghosts_0_data) +
                    (j + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                dudy[idx] = a_n*(u[idx_data_T] - u[idx_data_B])*dx_1_inv;
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and number of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_lo_2 = domain_lo[2];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        const int domain_dim_2 = domain_dims[2];
        
        const int num_ghosts_0_derivative_midpoint = num_ghosts_derivative_midpoint[0];
        const int num_ghosts_1_derivative_midpoint = num_ghosts_derivative_midpoint[1];
        const int num_ghosts_2_derivative_midpoint = num_ghosts_derivative_midpoint[2];
        const int ghostcell_dim_0_derivative_midpoint = ghostcell_dims_derivative_midpoint[0];
        const int ghostcell_dim_1_derivative_midpoint = ghostcell_dims_derivative_midpoint[1];
        
        const int num_ghosts_0_data = num_ghosts_data_node[0];
        const int num_ghosts_1_data = num_ghosts_data_node[1];
        const int num_ghosts_2_data = num_ghosts_data_node[2];
        const int ghostcell_dim_0_data = ghostcell_dims_data_node[0];
        const int ghostcell_dim_1_data = ghostcell_dims_data_node[1];
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx = (i + num_ghosts_0_derivative_midpoint) +
                        (j + num_ghosts_1_derivative_midpoint)*ghostcell_dim_0_derivative_midpoint +
                        (k + num_ghosts_2_derivative_midpoint)*ghostcell_dim_0_derivative_midpoint*
                            (ghostcell_dim_1_derivative_midpoint + 1);
                    
                    const int idx_data_B = (i + num_ghosts_0_data) +
                        (j - 1 + num_ghosts_1_data)*ghostcell_dim_0_data +
                        (k + num_ghosts_2_data)*ghostcell_dim_0_data*
                            ghostcell_dim_1_data;
                    
                    const int idx_data_T = (i + num_ghosts_0_data) +
                        (j + num_ghosts_1_data)*ghostcell_dim_0_data +
                        (k + num_ghosts_2_data)*ghostcell_dim_0_data*
                            ghostcell_dim_1_data;
                    
                    dudy[idx] = a_n*(u[idx_data_T] - u[idx_data_B])*dx_1_inv;
                }
            }
        }
    }
}


/*
 * Kernel to compute the derivatives in z-direction at midpoints.
 */
void
DiffusiveFluxReconstructorMidpointSecondOrder::computeFirstDerivativesInZAtMidpointZ(
    double* dudz,
    const double* const u,
    const hier::IntVector& num_ghosts_derivative_midpoint,
    const hier::IntVector& num_ghosts_data_node,
    const hier::IntVector& ghostcell_dims_derivative_midpoint,
    const hier::IntVector& ghostcell_dims_data_node,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims,
    const double& dx_2_inv) const
{
    const double a_n = double(1);
    
    /*
     * Get the local lower index, numbers of cells in each dimension and number of ghost cells.
     */
    
    const int domain_lo_0 = domain_lo[0];
    const int domain_lo_1 = domain_lo[1];
    const int domain_lo_2 = domain_lo[2];
    const int domain_dim_0 = domain_dims[0];
    const int domain_dim_1 = domain_dims[1];
    const int domain_dim_2 = domain_dims[2];
    
    const int num_ghosts_0_derivative_midpoint = num_ghosts_derivative_midpoint[0];
    const int num_ghosts_1_derivative_midpoint = num_ghosts_derivative_midpoint[1];
    const int num_ghosts_2_derivative_midpoint = num_ghosts_derivative_midpoint[2];
    const int ghostcell_dim_0_derivative_midpoint = ghostcell_dims_derivative_midpoint[0];
    const int ghostcell_dim_1_derivative_midpoint = ghostcell_dims_derivative_midpoint[1];
    
    const int num_ghosts_0_data = num_ghosts_data_node[0];
    const int num_ghosts_1_data = num_ghosts_data_node[1];
    const int num_ghosts_2_data = num_ghosts_data_node[2];
    const int ghostcell_dim_0_data = ghostcell_dims_data_node[0];
    const int ghostcell_dim_1_data = ghostcell_dims_data_node[1];
    
    for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
    {
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx = (i + num_ghosts_0_derivative_midpoint) +
                    (j + num_ghosts_1_derivative_midpoint)*ghostcell_dim_0_derivative_midpoint +
                    (k + num_ghosts_2_derivative_midpoint)*ghostcell_dim_0_derivative_midpoint*
                        ghostcell_dim_1_derivative_midpoint;
                
                const int idx_data_B = (i + num_ghosts_0_data) +
                    (j + num_ghosts_1_data)*ghostcell_dim_0_data +
                    (k - 1 + num_ghosts_2_data)*ghostcell_dim_0_data*
                        ghostcell_dim_1_data;
                
                const int idx_data_F = (i + num_ghosts_0_data) +
                    (j + num_ghosts_1_data)*ghostcell_dim_0_data +
                    (k + num_ghosts_2_data)*ghostcell_dim_0_data*
                        ghostcell_dim_1_data;
                
                dudz[idx] = a_n*(u[idx_data_F] - u[idx_data_B])*dx_2_inv;
            }
        }
    }
}


/*
 * Kernel to compute the derivatives in x-direction at nodes.
 */
void
DiffusiveFluxReconstructorMidpointSecondOrder::computeFirstDerivativesInXAtNode(
    double* dudx,
    const double* const u,
    const hier::IntVector& num_ghosts_derivative_node,
    const hier::IntVector& num_ghosts_data_node,
    const hier::IntVector& ghostcell_dims_derivative_node,
    const hier::IntVector& ghostcell_dims_data_node,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims,
    const double& dx_0_inv) const
{
    const double a_n = double(1)/double(2);
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and number of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int num_ghosts_0_derivative_node = num_ghosts_derivative_node[0];
        
        const int num_ghosts_0_data = num_ghosts_data_node[0];
        
        
        HAMERS_PRAGMA_SIMD
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx = i + num_ghosts_0_derivative_node;
            
            const int idx_data_L = i - 1 + num_ghosts_0_data;
            const int idx_data_R = i + 1 + num_ghosts_0_data;
            
            dudx[idx] = a_n*(u[idx_data_R] - u[idx_data_L])*dx_0_inv;
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and number of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        
        const int num_ghosts_0_derivative_node = num_ghosts_derivative_node[0];
        const int num_ghosts_1_derivative_node = num_ghosts_derivative_node[1];
        const int ghostcell_dim_0_derivative_node = ghostcell_dims_derivative_node[0];
        
        const int num_ghosts_0_data = num_ghosts_data_node[0];
        const int num_ghosts_1_data = num_ghosts_data_node[1];
        const int ghostcell_dim_0_data = ghostcell_dims_data_node[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx = (i + num_ghosts_0_derivative_node) +
                    (j + num_ghosts_1_derivative_node)*ghostcell_dim_0_derivative_node;
                
                const int idx_data_L = (i - 1 + num_ghosts_0_data) +
                    (j + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                const int idx_data_R = (i + 1 + num_ghosts_0_data) +
                    (j + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                dudx[idx] = a_n*(u[idx_data_R] - u[idx_data_L])*dx_0_inv;
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and number of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_lo_2 = domain_lo[2];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        const int domain_dim_2 = domain_dims[2];
        
        const int num_ghosts_0_derivative_node = num_ghosts_derivative_node[0];
        const int num_ghosts_1_derivative_node = num_ghosts_derivative_node[1];
        const int num_ghosts_2_derivative_node = num_ghosts_derivative_node[2];
        const int ghostcell_dim_0_derivative_node = ghostcell_dims_derivative_node[0];
        const int ghostcell_dim_1_derivative_node = ghostcell_dims_derivative_node[1];
        
        const int num_ghosts_0_data = num_ghosts_data_node[0];
        const int num_ghosts_1_data = num_ghosts_data_node[1];
        const int num_ghosts_2_data = num_ghosts_data_node[2];
        const int ghostcell_dim_0_data = ghostcell_dims_data_node[0];
        const int ghostcell_dim_1_data = ghostcell_dims_data_node[1];
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx = (i + num_ghosts_0_derivative_node) +
                        (j + num_ghosts_1_derivative_node)*ghostcell_dim_0_derivative_node +
                        (k + num_ghosts_2_derivative_node)*ghostcell_dim_0_derivative_node*
                            ghostcell_dim_1_derivative_node;
                    
                    const int idx_data_L = (i - 1 + num_ghosts_0_data) +
                        (j + num_ghosts_1_data)*ghostcell_dim_0_data +
                        (k + num_ghosts_2_data)*ghostcell_dim_0_data*
                            ghostcell_dim_1_data;
                    
                    const int idx_data_R = (i + 1 + num_ghosts_0_data) +
                        (j + num_ghosts_1_data)*ghostcell_dim_0_data +
                        (k + num_ghosts_2_data)*ghostcell_dim_0_data*
                            ghostcell_dim_1_data;
                    
                    dudx[idx] = a_n*(u[idx_data_R] - u[idx_data_L])*dx_0_inv;
                }
            }
        }
    }
}


/*
 * Kernel to compute the derivatives in y-direction at nodes.
 */
void
DiffusiveFluxReconstructorMidpointSecondOrder::computeFirstDerivativesInYAtNode(
    double* dudy,
    const double* const u,
    const hier::IntVector& num_ghosts_derivative_node,
    const hier::IntVector& num_ghosts_data_node,
    const hier::IntVector& ghostcell_dims_derivative_node,
    const hier::IntVector& ghostcell_dims_data_node,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims,
    const double& dx_1_inv) const
{
    const double a_n = double(1)/double(2);
    
    if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and number of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        
        const int num_ghosts_0_derivative_node = num_ghosts_derivative_node[0];
        const int num_ghosts_1_derivative_node = num_ghosts_derivative_node[1];
        const int ghostcell_dim_0_derivative_node = ghostcell_dims_derivative_node[0];
        
        const int num_ghosts_0_data = num_ghosts_data_node[0];
        const int num_ghosts_1_data = num_ghosts_data_node[1];
        const int ghostcell_dim_0_data = ghostcell_dims_data_node[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx = (i + num_ghosts_0_derivative_node) +
                    (j + num_ghosts_1_derivative_node)*ghostcell_dim_0_derivative_node;
                
                const int idx_data_B = (i + num_ghosts_0_data) +
                    (j - 1 + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                const int idx_data_T = (i + num_ghosts_0_data) +
                    (j + 1 + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                dudy[idx] = a_n*(u[idx_data_T] - u[idx_data_B])*dx_1_inv;
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and number of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_lo_2 = domain_lo[2];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        const int domain_dim_2 = domain_dims[2];
        
        const int num_ghosts_0_derivative_node = num_ghosts_derivative_node[0];
        const int num_ghosts_1_derivative_node = num_ghosts_derivative_node[1];
        const int num_ghosts_2_derivative_node = num_ghosts_derivative_node[2];
        const int ghostcell_dim_0_derivative_node = ghostcell_dims_derivative_node[0];
        const int ghostcell_dim_1_derivative_node = ghostcell_dims_derivative_node[1];
        
        const int num_ghosts_0_data = num_ghosts_data_node[0];
        const int num_ghosts_1_data = num_ghosts_data_node[1];
        const int num_ghosts_2_data = num_ghosts_data_node[2];
        const int ghostcell_dim_0_data = ghostcell_dims_data_node[0];
        const int ghostcell_dim_1_data = ghostcell_dims_data_node[1];
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx = (i + num_ghosts_0_derivative_node) +
                        (j + num_ghosts_1_derivative_node)*ghostcell_dim_0_derivative_node +
                        (k + num_ghosts_2_derivative_node)*ghostcell_dim_0_derivative_node*
                            ghostcell_dim_1_derivative_node;
                    
                    const int idx_data_B = (i + num_ghosts_0_data) +
                        (j - 1 + num_ghosts_1_data)*ghostcell_dim_0_data +
                        (k + num_ghosts_2_data)*ghostcell_dim_0_data*
                            ghostcell_dim_1_data;
                    
                    const int idx_data_T = (i + num_ghosts_0_data) +
                        (j + 1 + num_ghosts_1_data)*ghostcell_dim_0_data +
                        (k + num_ghosts_2_data)*ghostcell_dim_0_data*
                            ghostcell_dim_1_data;
                    
                    dudy[idx] = a_n*(u[idx_data_T] - u[idx_data_B])*dx_1_inv;
                }
            }
        }
    }
}


/*
 * Kernel to compute the derivatives in z-direction at nodes.
 */
void
DiffusiveFluxReconstructorMidpointSecondOrder::computeFirstDerivativesInZAtNode(
    double* dudz,
    const double* const u,
    const hier::IntVector& num_ghosts_derivative_node,
    const hier::IntVector& num_ghosts_data_node,
    const hier::IntVector& ghostcell_dims_derivative_node,
    const hier::IntVector& ghostcell_dims_data_node,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims,
    const double& dx_2_inv) const
{
    const double a_n = double(1)/double(2);
    
    /*
     * Get the local lower index, numbers of cells in each dimension and number of ghost cells.
     */
    
    const int domain_lo_0 = domain_lo[0];
    const int domain_lo_1 = domain_lo[1];
    const int domain_lo_2 = domain_lo[2];
    const int domain_dim_0 = domain_dims[0];
    const int domain_dim_1 = domain_dims[1];
    const int domain_dim_2 = domain_dims[2];
    
    const int num_ghosts_0_derivative_node = num_ghosts_derivative_node[0];
    const int num_ghosts_1_derivative_node = num_ghosts_derivative_node[1];
    const int num_ghosts_2_derivative_node = num_ghosts_derivative_node[2];
    const int ghostcell_dim_0_derivative_node = ghostcell_dims_derivative_node[0];
    const int ghostcell_dim_1_derivative_node = ghostcell_dims_derivative_node[1];
    
    const int num_ghosts_0_data = num_ghosts_data_node[0];
    const int num_ghosts_1_data = num_ghosts_data_node[1];
    const int num_ghosts_2_data = num_ghosts_data_node[2];
    const int ghostcell_dim_0_data = ghostcell_dims_data_node[0];
    const int ghostcell_dim_1_data = ghostcell_dims_data_node[1];
    
    for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
    {
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx = (i + num_ghosts_0_derivative_node) +
                    (j + num_ghosts_1_derivative_node)*ghostcell_dim_0_derivative_node +
                    (k + num_ghosts_2_derivative_node)*ghostcell_dim_0_derivative_node*
                        ghostcell_dim_1_derivative_node;
                
                const int idx_data_B = (i + num_ghosts_0_data) +
                    (j + num_ghosts_1_data)*ghostcell_dim_0_data +
                    (k - 1 + num_ghosts_2_data)*ghostcell_dim_0_data*
                        ghostcell_dim_1_data;
                
                const int idx_data_F = (i + num_ghosts_0_data) +
                    (j + num_ghosts_1_data)*ghostcell_dim_0_data +
                    (k + 1 + num_ghosts_2_data)*ghostcell_dim_0_data*
                        ghostcell_dim_1_data;
                
                dudz[idx] = a_n*(u[idx_data_F] - u[idx_data_B])*dx_2_inv;
            }
        }
    }
}


/*
 * Kernel to interpolate the data from nodes to midpoints in x-direction.
 */
void
DiffusiveFluxReconstructorMidpointSecondOrder::interpolateDataFromNodeToMidpointX(
    double* u_midpoint_x,
    const double* const u_node,
    const hier::IntVector& num_ghosts_data_midpoint,
    const hier::IntVector& num_ghosts_data_node,
    const hier::IntVector& ghostcell_dims_data_midpoint,
    const hier::IntVector& ghostcell_dims_data_node,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims) const
{
    const double a_n = double(1)/double(2);
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and number of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int num_ghosts_0_data_midpoint = num_ghosts_data_midpoint[0];
        
        const int num_ghosts_0_data_node = num_ghosts_data_node[0];
        
        HAMERS_PRAGMA_SIMD
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx = i + num_ghosts_0_data_midpoint;
            
            const int idx_node_L = i - 1 + num_ghosts_0_data_node;
            const int idx_node_R = i + 0 + num_ghosts_0_data_node;
            
            u_midpoint_x[idx] = a_n*(u_node[idx_node_R] + u_node[idx_node_L]);
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and number of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        
        const int num_ghosts_0_data_midpoint = num_ghosts_data_midpoint[0];
        const int num_ghosts_1_data_midpoint = num_ghosts_data_midpoint[1];
        const int ghostcell_dim_0_data_midpoint = ghostcell_dims_data_midpoint[0];
        
        const int num_ghosts_0_data_node = num_ghosts_data_node[0];
        const int num_ghosts_1_data_node = num_ghosts_data_node[1];
        const int ghostcell_dim_0_data_node = ghostcell_dims_data_node[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx = (i + num_ghosts_0_data_midpoint) +
                    (j + num_ghosts_1_data_midpoint)*(ghostcell_dim_0_data_midpoint + 1);
                
                const int idx_node_L = (i - 1 + num_ghosts_0_data_node) +
                    (j + num_ghosts_1_data_node)*ghostcell_dim_0_data_node;
                
                const int idx_node_R = (i + num_ghosts_0_data_node) +
                    (j + num_ghosts_1_data_node)*ghostcell_dim_0_data_node;
                
                u_midpoint_x[idx] = a_n*(u_node[idx_node_R] + u_node[idx_node_L]);
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and number of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_lo_2 = domain_lo[2];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        const int domain_dim_2 = domain_dims[2];
        
        const int num_ghosts_0_data_midpoint = num_ghosts_data_midpoint[0];
        const int num_ghosts_1_data_midpoint = num_ghosts_data_midpoint[1];
        const int num_ghosts_2_data_midpoint = num_ghosts_data_midpoint[2];
        const int ghostcell_dim_0_data_midpoint = ghostcell_dims_data_midpoint[0];
        const int ghostcell_dim_1_data_midpoint = ghostcell_dims_data_midpoint[1];
        
        const int num_ghosts_0_data_node = num_ghosts_data_node[0];
        const int num_ghosts_1_data_node = num_ghosts_data_node[1];
        const int num_ghosts_2_data_node = num_ghosts_data_node[2];
        const int ghostcell_dim_0_data_node = ghostcell_dims_data_node[0];
        const int ghostcell_dim_1_data_node = ghostcell_dims_data_node[1];
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx = (i + num_ghosts_0_data_midpoint) +
                        (j + num_ghosts_1_data_midpoint)*(ghostcell_dim_0_data_midpoint + 1) +
                        (k + num_ghosts_2_data_midpoint)*(ghostcell_dim_0_data_midpoint + 1)*
                            ghostcell_dim_1_data_midpoint;
                    
                    const int idx_node_L = (i - 1 + num_ghosts_0_data_node) +
                        (j + num_ghosts_1_data_node)*ghostcell_dim_0_data_node +
                        (k + num_ghosts_2_data_node)*ghostcell_dim_0_data_node*
                            ghostcell_dim_1_data_node;
                    
                    const int idx_node_R = (i + num_ghosts_0_data_node) +
                        (j + num_ghosts_1_data_node)*ghostcell_dim_0_data_node +
                        (k + num_ghosts_2_data_node)*ghostcell_dim_0_data_node*
                            ghostcell_dim_1_data_node;
                    
                    u_midpoint_x[idx] = a_n*(u_node[idx_node_R] + u_node[idx_node_L]);
                }
            }
        }
    }
}


/*
 * Kernel to interpolate the data from nodes to midpoints in y-direction.
 */
void
DiffusiveFluxReconstructorMidpointSecondOrder::interpolateDataFromNodeToMidpointY(
    double* u_midpoint_y,
    const double* const u_node,
    const hier::IntVector& num_ghosts_data_midpoint,
    const hier::IntVector& num_ghosts_data_node,
    const hier::IntVector& ghostcell_dims_data_midpoint,
    const hier::IntVector& ghostcell_dims_data_node,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims) const
{
    const double a_n = double(1)/double(2);
    
    if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and number of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        
        const int num_ghosts_0_data_midpoint = num_ghosts_data_midpoint[0];
        const int num_ghosts_1_data_midpoint = num_ghosts_data_midpoint[1];
        const int ghostcell_dim_0_data_midpoint = ghostcell_dims_data_midpoint[0];
        
        const int num_ghosts_0_data_node = num_ghosts_data_node[0];
        const int num_ghosts_1_data_node = num_ghosts_data_node[1];
        const int ghostcell_dim_0_data_node = ghostcell_dims_data_node[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx = (i + num_ghosts_0_data_midpoint) +
                    (j + num_ghosts_1_data_midpoint)*ghostcell_dim_0_data_midpoint;
                
                const int idx_node_B = (i + num_ghosts_0_data_node) +
                    (j - 1 + num_ghosts_1_data_node)*ghostcell_dim_0_data_node;
                
                const int idx_node_T = (i + num_ghosts_0_data_node) +
                    (j + num_ghosts_1_data_node)*ghostcell_dim_0_data_node;
                
                u_midpoint_y[idx] = a_n*(u_node[idx_node_T] + u_node[idx_node_B]);
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and number of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_lo_2 = domain_lo[2];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        const int domain_dim_2 = domain_dims[2];
        
        const int num_ghosts_0_data_midpoint = num_ghosts_data_midpoint[0];
        const int num_ghosts_1_data_midpoint = num_ghosts_data_midpoint[1];
        const int num_ghosts_2_data_midpoint = num_ghosts_data_midpoint[2];
        const int ghostcell_dim_0_data_midpoint = ghostcell_dims_data_midpoint[0];
        const int ghostcell_dim_1_data_midpoint = ghostcell_dims_data_midpoint[1];
        
        const int num_ghosts_0_data_node = num_ghosts_data_node[0];
        const int num_ghosts_1_data_node = num_ghosts_data_node[1];
        const int num_ghosts_2_data_node = num_ghosts_data_node[2];
        const int ghostcell_dim_0_data_node = ghostcell_dims_data_node[0];
        const int ghostcell_dim_1_data_node = ghostcell_dims_data_node[1];
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx = (i + num_ghosts_0_data_midpoint) +
                        (j + num_ghosts_1_data_midpoint)*ghostcell_dim_0_data_midpoint +
                        (k + num_ghosts_2_data_midpoint)*ghostcell_dim_0_data_midpoint*
                            (ghostcell_dim_1_data_midpoint + 1);
                    
                    const int idx_node_B = (i + num_ghosts_0_data_node) +
                        (j - 1 + num_ghosts_1_data_node)*ghostcell_dim_0_data_node +
                        (k + num_ghosts_2_data_node)*ghostcell_dim_0_data_node*
                            ghostcell_dim_1_data_node;
                    
                    const int idx_node_T = (i + num_ghosts_0_data_node) +
                        (j + num_ghosts_1_data_node)*ghostcell_dim_0_data_node +
                        (k + num_ghosts_2_data_node)*ghostcell_dim_0_data_node*
                            ghostcell_dim_1_data_node;
                    
                    u_midpoint_y[idx] = a_n*(u_node[idx_node_T] + u_node[idx_node_B]);
                }
            }
        }
    }
}


/*
 * Kernel to interpolate the data from nodes to midpoints in z-direction.
 */
void
DiffusiveFluxReconstructorMidpointSecondOrder::interpolateDataFromNodeToMidpointZ(
    double* u_midpoint_z,
    const double* const u_node,
    const hier::IntVector& num_ghosts_data_midpoint,
    const hier::IntVector& num_ghosts_data_node,
    const hier::IntVector& ghostcell_dims_data_midpoint,
    const hier::IntVector& ghostcell_dims_data_node,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims) const
{
    const double a_n = double(1)/double(2);
    
    /*
     * Get the local lower index, numbers of cells in each dimension and number of ghost cells.
     */
    
    const int domain_lo_0 = domain_lo[0];
    const int domain_lo_1 = domain_lo[1];
    const int domain_lo_2 = domain_lo[2];
    const int domain_dim_0 = domain_dims[0];
    const int domain_dim_1 = domain_dims[1];
    const int domain_dim_2 = domain_dims[2];
    
    const int num_ghosts_0_data_midpoint = num_ghosts_data_midpoint[0];
    const int num_ghosts_1_data_midpoint = num_ghosts_data_midpoint[1];
    const int num_ghosts_2_data_midpoint = num_ghosts_data_midpoint[2];
    const int ghostcell_dim_0_data_midpoint = ghostcell_dims_data_midpoint[0];
    const int ghostcell_dim_1_data_midpoint = ghostcell_dims_data_midpoint[1];
    
    const int num_ghosts_0_data_node = num_ghosts_data_node[0];
    const int num_ghosts_1_data_node = num_ghosts_data_node[1];
    const int num_ghosts_2_data_node = num_ghosts_data_node[2];
    const int ghostcell_dim_0_data_node = ghostcell_dims_data_node[0];
    const int ghostcell_dim_1_data_node = ghostcell_dims_data_node[1];
    
    for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
    {
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx = (i + num_ghosts_0_data_midpoint) +
                    (j + num_ghosts_1_data_midpoint)*ghostcell_dim_0_data_midpoint +
                    (k + num_ghosts_2_data_midpoint)*ghostcell_dim_0_data_midpoint*
                        ghostcell_dim_1_data_midpoint;
                
                const int idx_node_B = (i + num_ghosts_0_data_node) +
                    (j + num_ghosts_1_data_node)*ghostcell_dim_0_data_node +
                    (k - 1 + num_ghosts_2_data_node)*ghostcell_dim_0_data_node*
                        ghostcell_dim_1_data_node;
                
                const int idx_node_F = (i + num_ghosts_0_data_node) +
                    (j + num_ghosts_1_data_node)*ghostcell_dim_0_data_node +
                    (k + num_ghosts_2_data_node)*ghostcell_dim_0_data_node*
                        ghostcell_dim_1_data_node;
                
                u_midpoint_z[idx] = a_n*(u_node[idx_node_F] + u_node[idx_node_B]);
            }
        }
    }
}


/*
 * Kernel to reconstruct the flux using flux at midpoints in x-direction.
 */
void
DiffusiveFluxReconstructorMidpointSecondOrder::reconstructFluxX(
    double* F_face_x,
    const double* const F_midpoint_x,
    const hier::IntVector& num_ghosts_flux_midpoint,
    const hier::IntVector& ghostcell_dims_flux_midpoint,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims,
    const hier::IntVector& interior_dims,
    const double& dt) const
{
    const double a_m = double(1);
    
    const double a_r = a_m;
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and number of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int num_ghosts_0_flux_midpoint = num_ghosts_flux_midpoint[0];
        
        HAMERS_PRAGMA_SIMD
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_face_x = i;
            
            const int idx_midpoint = i + 0 + num_ghosts_0_flux_midpoint;
            
            F_face_x[idx_face_x] += dt*a_r*F_midpoint_x[idx_midpoint];
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and number of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        
        const int interior_dim_0 = interior_dims[0];
        
        const int num_ghosts_0_flux_midpoint = num_ghosts_flux_midpoint[0];
        const int num_ghosts_1_flux_midpoint = num_ghosts_flux_midpoint[1];
        const int ghostcell_dim_0_flux_midpoint = ghostcell_dims_flux_midpoint[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_face_x = i +
                    j*(interior_dim_0 + 1);
                
                const int idx_midpoint = (i + num_ghosts_0_flux_midpoint) +
                    (j + num_ghosts_1_flux_midpoint)*(ghostcell_dim_0_flux_midpoint + 1);
                
                F_face_x[idx_face_x] += dt*a_r*F_midpoint_x[idx_midpoint];
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and number of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_lo_2 = domain_lo[2];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        const int domain_dim_2 = domain_dims[2];
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        const int num_ghosts_0_flux_midpoint = num_ghosts_flux_midpoint[0];
        const int num_ghosts_1_flux_midpoint = num_ghosts_flux_midpoint[1];
        const int num_ghosts_2_flux_midpoint = num_ghosts_flux_midpoint[2];
        const int ghostcell_dim_0_flux_midpoint = ghostcell_dims_flux_midpoint[0];
        const int ghostcell_dim_1_flux_midpoint = ghostcell_dims_flux_midpoint[1];
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_face_x = i +
                        j*(interior_dim_0 + 1) +
                        k*(interior_dim_0 + 1)*interior_dim_1;
                    
                    const int idx_midpoint = (i + num_ghosts_0_flux_midpoint) +
                        (j + num_ghosts_1_flux_midpoint)*(ghostcell_dim_0_flux_midpoint + 1) +
                        (k + num_ghosts_2_flux_midpoint)*(ghostcell_dim_0_flux_midpoint + 1)*
                            ghostcell_dim_1_flux_midpoint;
                    
                    F_face_x[idx_face_x] += dt*a_r*F_midpoint_x[idx_midpoint];
                }
            }
        }
    }
}


/*
 * Kernel to reconstruct the flux using flux at midpoints in y-direction.
 */
void
DiffusiveFluxReconstructorMidpointSecondOrder::reconstructFluxY(
    double* F_face_y,
    const double* const F_midpoint_y,
    const hier::IntVector& num_ghosts_flux_midpoint,
    const hier::IntVector& ghostcell_dims_flux_midpoint,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims,
    const hier::IntVector& interior_dims,
    const double& dt) const
{
    const double a_m = double(1);
    
    const double a_r = a_m;
    
    if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and number of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        
        const int interior_dim_0 = interior_dims[0];
        
        const int num_ghosts_0_flux_midpoint = num_ghosts_flux_midpoint[0];
        const int num_ghosts_1_flux_midpoint = num_ghosts_flux_midpoint[1];
        const int ghostcell_dim_0_flux_midpoint = ghostcell_dims_flux_midpoint[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_face_y = i +
                    j*interior_dim_0;
                
                const int idx_midpoint = (i + num_ghosts_0_flux_midpoint) +
                    (j + num_ghosts_1_flux_midpoint)*ghostcell_dim_0_flux_midpoint;
                
                F_face_y[idx_face_y] += dt*a_r*F_midpoint_y[idx_midpoint];
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and number of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_lo_2 = domain_lo[2];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        const int domain_dim_2 = domain_dims[2];
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        const int num_ghosts_0_flux_midpoint = num_ghosts_flux_midpoint[0];
        const int num_ghosts_1_flux_midpoint = num_ghosts_flux_midpoint[1];
        const int num_ghosts_2_flux_midpoint = num_ghosts_flux_midpoint[2];
        const int ghostcell_dim_0_flux_midpoint = ghostcell_dims_flux_midpoint[0];
        const int ghostcell_dim_1_flux_midpoint = ghostcell_dims_flux_midpoint[1];
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_face_y = i +
                        j*interior_dim_0 +
                        k*interior_dim_0*(interior_dim_1 + 1);
                    
                    const int idx_midpoint = (i + num_ghosts_0_flux_midpoint) +
                        (j + num_ghosts_1_flux_midpoint)*ghostcell_dim_0_flux_midpoint +
                        (k + num_ghosts_2_flux_midpoint)*ghostcell_dim_0_flux_midpoint*
                            (ghostcell_dim_1_flux_midpoint + 1);
                    
                    F_face_y[idx_face_y] += dt*a_r*F_midpoint_y[idx_midpoint];
                }
            }
        }
    }
}


/*
 * Kernel to reconstruct the flux using flux at midpoints in z-direction.
 */
void
DiffusiveFluxReconstructorMidpointSecondOrder::reconstructFluxZ(
    double* F_face_z,
    const double* const F_midpoint_z,
    const hier::IntVector& num_ghosts_flux_midpoint,
    const hier::IntVector& ghostcell_dims_flux_midpoint,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims,
    const hier::IntVector& interior_dims,
    const double& dt) const
{
    const double a_m = double(1);
    
    const double a_r = a_m;
    
    /*
     * Get the local lower index, numbers of cells in each dimension and number of ghost cells.
     */
    
    const int domain_lo_0 = domain_lo[0];
    const int domain_lo_1 = domain_lo[1];
    const int domain_lo_2 = domain_lo[2];
    const int domain_dim_0 = domain_dims[0];
    const int domain_dim_1 = domain_dims[1];
    const int domain_dim_2 = domain_dims[2];
    
    const int interior_dim_0 = interior_dims[0];
    const int interior_dim_1 = interior_dims[1];
    
    const int num_ghosts_0_flux_midpoint = num_ghosts_flux_midpoint[0];
    const int num_ghosts_1_flux_midpoint = num_ghosts_flux_midpoint[1];
    const int num_ghosts_2_flux_midpoint = num_ghosts_flux_midpoint[2];
    const int ghostcell_dim_0_flux_midpoint = ghostcell_dims_flux_midpoint[0];
    const int ghostcell_dim_1_flux_midpoint = ghostcell_dims_flux_midpoint[1];
    
    for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
    {
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_face_z = i +
                    j*interior_dim_0 +
                    k*interior_dim_0*interior_dim_1;
                
                const int idx_midpoint = (i + num_ghosts_0_flux_midpoint) +
                    (j + num_ghosts_1_flux_midpoint)*ghostcell_dim_0_flux_midpoint +
                    (k + num_ghosts_2_flux_midpoint)*ghostcell_dim_0_flux_midpoint*
                        ghostcell_dim_1_flux_midpoint;
                
                F_face_z[idx_face_z] += dt*a_r*F_midpoint_z[idx_midpoint];
            }
        }
    }
}
