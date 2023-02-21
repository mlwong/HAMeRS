#include "flow/diffusive_flux_reconstructors/node/DiffusiveFluxReconstructorNodeFourthOrder.hpp"

#include "SAMRAI/geom/CartesianPatchGeometry.h"

#include <map>

DiffusiveFluxReconstructorNodeFourthOrder::DiffusiveFluxReconstructorNodeFourthOrder(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_eqn,
    const HAMERS_SHARED_PTR<FlowModel>& flow_model,
    const HAMERS_SHARED_PTR<tbox::Database>& diffusive_flux_reconstructor_db):
        DiffusiveFluxReconstructorNode(
            object_name,
            dim,
            grid_geometry,
            num_eqn,
            flow_model,
            diffusive_flux_reconstructor_db)
{
    d_num_diff_ghosts = hier::IntVector::getOne(d_dim)*4;
    
    d_num_der_node_ghosts         = hier::IntVector::getOne(d_dim)*2;
    d_num_flux_reconstruct_ghosts = hier::IntVector::getOne(d_dim)*2;
}


/*
 * Print all characteristics of the diffusive flux reconstruction class.
 */
void
DiffusiveFluxReconstructorNodeFourthOrder::printClassData(
    std::ostream& os) const
{
    os << "\nPrint DiffusiveFluxReconstructorNodeFourthOrder object..."
       << std::endl;
    
    os << std::endl;
    
    os << "DiffusiveFluxReconstructorNodeFourthOrder: this = "
       << (DiffusiveFluxReconstructorNodeFourthOrder *)this
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
DiffusiveFluxReconstructorNodeFourthOrder::putToRestart(
   const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const
{
    restart_db->putString("d_diffusive_flux_reconstructor", "FOURTH_ORDER");
}


/*
 * Kernel to compute the first derivatives in the x-direction.
 */
void
DiffusiveFluxReconstructorNodeFourthOrder::computeFirstDerivativesInX(
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
    const double a_n =  double(2)/double(3);
    const double b_n = -double(1)/double(12);
    
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
            
            const int idx_data_LL = i - 2 + num_ghosts_0_data;
            const int idx_data_L  = i - 1 + num_ghosts_0_data;
            const int idx_data_R  = i + 1 + num_ghosts_0_data;
            const int idx_data_RR = i + 2 + num_ghosts_0_data;
            
            dudx[idx] = (a_n*(u[idx_data_R]  - u[idx_data_L]) +
                         b_n*(u[idx_data_RR] - u[idx_data_LL])
                        )*dx_0_inv;
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
                
                const int idx_data_LL = (i - 2 + num_ghosts_0_data) +
                    (j + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                const int idx_data_L = (i - 1 + num_ghosts_0_data) +
                    (j + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                const int idx_data_R = (i + 1 + num_ghosts_0_data) +
                    (j + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                const int idx_data_RR = (i + 2 + num_ghosts_0_data) +
                    (j + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                dudx[idx] = (a_n*(u[idx_data_R]  - u[idx_data_L]) +
                             b_n*(u[idx_data_RR] - u[idx_data_LL])
                            )*dx_0_inv;
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
                    
                    const int idx_data_LL = (i - 2 + num_ghosts_0_data) +
                        (j + num_ghosts_1_data)*ghostcell_dim_0_data +
                        (k + num_ghosts_2_data)*ghostcell_dim_0_data*
                            ghostcell_dim_1_data;
                    
                    const int idx_data_L = (i - 1 + num_ghosts_0_data) +
                        (j + num_ghosts_1_data)*ghostcell_dim_0_data +
                        (k + num_ghosts_2_data)*ghostcell_dim_0_data*
                            ghostcell_dim_1_data;
                    
                    const int idx_data_R = (i + 1 + num_ghosts_0_data) +
                        (j + num_ghosts_1_data)*ghostcell_dim_0_data +
                        (k + num_ghosts_2_data)*ghostcell_dim_0_data*
                            ghostcell_dim_1_data;
                    
                    const int idx_data_RR = (i + 2 + num_ghosts_0_data) +
                        (j + num_ghosts_1_data)*ghostcell_dim_0_data +
                        (k + num_ghosts_2_data)*ghostcell_dim_0_data*
                            ghostcell_dim_1_data;
                    
                    dudx[idx] = (a_n*(u[idx_data_R]  - u[idx_data_L]) +
                                 b_n*(u[idx_data_RR] - u[idx_data_LL])
                                )*dx_0_inv;
                }
            }
        }
    }
}


/*
 * Kernel to compute the first derivatives in the y-direction.
 */
void
DiffusiveFluxReconstructorNodeFourthOrder::computeFirstDerivativesInY(
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
    const double a_n =  double(2)/double(3);
    const double b_n = -double(1)/double(12);
    
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
                
                const int idx_data_BB = (i + num_ghosts_0_data) +
                    (j - 2 + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                const int idx_data_B = (i + num_ghosts_0_data) +
                    (j - 1 + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                const int idx_data_T = (i + num_ghosts_0_data) +
                    (j + 1 + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                const int idx_data_TT = (i + num_ghosts_0_data) +
                    (j + 2 + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                dudy[idx] = (a_n*(u[idx_data_T]  - u[idx_data_B]) +
                             b_n*(u[idx_data_TT] - u[idx_data_BB])
                            )*dx_1_inv;
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
                    
                    const int idx_data_BB = (i + num_ghosts_0_data) +
                        (j - 2 + num_ghosts_1_data)*ghostcell_dim_0_data +
                        (k + num_ghosts_2_data)*ghostcell_dim_0_data*
                            ghostcell_dim_1_data;
                    
                    const int idx_data_B = (i + num_ghosts_0_data) +
                        (j - 1 + num_ghosts_1_data)*ghostcell_dim_0_data +
                        (k + num_ghosts_2_data)*ghostcell_dim_0_data*
                            ghostcell_dim_1_data;
                    
                    const int idx_data_T = (i + num_ghosts_0_data) +
                        (j + 1 + num_ghosts_1_data)*ghostcell_dim_0_data +
                        (k + num_ghosts_2_data)*ghostcell_dim_0_data*
                            ghostcell_dim_1_data;
                    
                    const int idx_data_TT = (i + num_ghosts_0_data) +
                        (j + 2 + num_ghosts_1_data)*ghostcell_dim_0_data +
                        (k + num_ghosts_2_data)*ghostcell_dim_0_data*
                            ghostcell_dim_1_data;
                    
                    dudy[idx] = (a_n*(u[idx_data_T]  - u[idx_data_B]) +
                                 b_n*(u[idx_data_TT] - u[idx_data_BB])
                                )*dx_1_inv;
                }
            }
        }
    }
}


/*
 * Kernel to compute the first derivatives in the z-direction.
 */
void
DiffusiveFluxReconstructorNodeFourthOrder::computeFirstDerivativesInZ(
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
    const double a_n =  double(2)/double(3);
    const double b_n = -double(1)/double(12);
    
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
                
                const int idx_data_BB = (i + num_ghosts_0_data) +
                    (j + num_ghosts_1_data)*ghostcell_dim_0_data +
                    (k - 2 + num_ghosts_2_data)*ghostcell_dim_0_data*
                        ghostcell_dim_1_data;
                
                const int idx_data_B = (i + num_ghosts_0_data) +
                    (j + num_ghosts_1_data)*ghostcell_dim_0_data +
                    (k - 1 + num_ghosts_2_data)*ghostcell_dim_0_data*
                        ghostcell_dim_1_data;
                
                const int idx_data_F = (i + num_ghosts_0_data) +
                    (j + num_ghosts_1_data)*ghostcell_dim_0_data +
                    (k + 1 + num_ghosts_2_data)*ghostcell_dim_0_data*
                        ghostcell_dim_1_data;
                
                const int idx_data_FF = (i + num_ghosts_0_data) +
                    (j + num_ghosts_1_data)*ghostcell_dim_0_data +
                    (k + 2 + num_ghosts_2_data)*ghostcell_dim_0_data*
                        ghostcell_dim_1_data;
                
                dudz[idx] = (a_n*(u[idx_data_F]  - u[idx_data_B]) +
                             b_n*(u[idx_data_FF] - u[idx_data_BB])
                            )*dx_2_inv;
            }
        }
    }
}


/*
 * Kernel to reconstruct the flux using flux at nodes in x-direction.
 */
void
DiffusiveFluxReconstructorNodeFourthOrder::reconstructFluxX(
    double* F_face_x,
    const double* const F_node_x,
    const hier::IntVector& num_ghosts_flux_node,
    const hier::IntVector& ghostcell_dims_flux_node,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims,
    const hier::IntVector& interior_dims,
    const double& dt) const
{
    const double a_n =  double(2)/double(3);
    const double b_n = -double(1)/double(12);
    
    const double a_r = a_n + b_n;
    const double b_r = b_n;
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and number of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int num_ghosts_0_flux_node = num_ghosts_flux_node[0];
        
        HAMERS_PRAGMA_SIMD
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_face_x = i;
            
            const int idx_node_LL = i - 2 + num_ghosts_0_flux_node;
            const int idx_node_L  = i - 1 + num_ghosts_0_flux_node;
            const int idx_node_R  = i + 0 + num_ghosts_0_flux_node;
            const int idx_node_RR = i + 1 + num_ghosts_0_flux_node;
            
            F_face_x[idx_face_x] += dt*(
                a_r*(F_node_x[idx_node_L]  + F_node_x[idx_node_R]) +
                b_r*(F_node_x[idx_node_LL] + F_node_x[idx_node_RR])
                );
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
        
        const int num_ghosts_0_flux_node = num_ghosts_flux_node[0];
        const int num_ghosts_1_flux_node = num_ghosts_flux_node[1];
        const int ghostcell_dim_0_flux_node = ghostcell_dims_flux_node[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_face_x = i +
                    j*(interior_dim_0 + 1);
                
                const int idx_node_LL = (i - 2 + num_ghosts_0_flux_node) +
                    (j + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node;
                
                const int idx_node_L = (i - 1 + num_ghosts_0_flux_node) +
                    (j + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node;
                
                const int idx_node_R = (i + num_ghosts_0_flux_node) +
                    (j + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node;
                
                const int idx_node_RR = (i + 1 + num_ghosts_0_flux_node) +
                    (j + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node;
                
                F_face_x[idx_face_x] += dt*(
                    a_r*(F_node_x[idx_node_L]  + F_node_x[idx_node_R]) +
                    b_r*(F_node_x[idx_node_LL] + F_node_x[idx_node_RR])
                    );
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
        
        const int num_ghosts_0_flux_node = num_ghosts_flux_node[0];
        const int num_ghosts_1_flux_node = num_ghosts_flux_node[1];
        const int num_ghosts_2_flux_node = num_ghosts_flux_node[2];
        const int ghostcell_dim_0_flux_node = ghostcell_dims_flux_node[0];
        const int ghostcell_dim_1_flux_node = ghostcell_dims_flux_node[1];
        
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
                    
                    const int idx_node_LL = (i - 2 + num_ghosts_0_flux_node) +
                        (j + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node +
                        (k + num_ghosts_2_flux_node)*ghostcell_dim_0_flux_node*
                            ghostcell_dim_1_flux_node;
                    
                    const int idx_node_L = (i - 1 + num_ghosts_0_flux_node) +
                        (j + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node +
                        (k + num_ghosts_2_flux_node)*ghostcell_dim_0_flux_node*
                            ghostcell_dim_1_flux_node;
                    
                    const int idx_node_R = (i + num_ghosts_0_flux_node) +
                        (j + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node +
                        (k + num_ghosts_2_flux_node)*ghostcell_dim_0_flux_node*
                            ghostcell_dim_1_flux_node;
                    
                    const int idx_node_RR = (i + 1 + num_ghosts_0_flux_node) +
                        (j + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node +
                        (k + num_ghosts_2_flux_node)*ghostcell_dim_0_flux_node*
                            ghostcell_dim_1_flux_node;
                    
                    F_face_x[idx_face_x] += dt*(
                        a_r*(F_node_x[idx_node_L]  + F_node_x[idx_node_R]) +
                        b_r*(F_node_x[idx_node_LL] + F_node_x[idx_node_RR])
                        );
                }
            }
        }
    }
}


/*
 * Kernel to reconstruct the flux using flux at nodes in y-direction.
 */
void
DiffusiveFluxReconstructorNodeFourthOrder::reconstructFluxY(
    double* F_face_y,
    const double* const F_node_y,
    const hier::IntVector& num_ghosts_flux_node,
    const hier::IntVector& ghostcell_dims_flux_node,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims,
    const hier::IntVector& interior_dims,
    const double& dt) const
{
    const double a_n =  double(2)/double(3);
    const double b_n = -double(1)/double(12);
    
    const double a_r = a_n + b_n;
    const double b_r = b_n;
    
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
        
        const int num_ghosts_0_flux_node = num_ghosts_flux_node[0];
        const int num_ghosts_1_flux_node = num_ghosts_flux_node[1];
        const int ghostcell_dim_0_flux_node = ghostcell_dims_flux_node[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_face_y = i +
                    j*interior_dim_0;
                
                const int idx_node_BB = (i + num_ghosts_0_flux_node) +
                    (j - 2 + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node;
                
                const int idx_node_B = (i + num_ghosts_0_flux_node) +
                    (j - 1 + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node;
                
                const int idx_node_T = (i + num_ghosts_0_flux_node) +
                    (j + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node;
                
                const int idx_node_TT = (i + num_ghosts_0_flux_node) +
                    (j + 1 + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node;
                
                F_face_y[idx_face_y] += dt*(
                    a_r*(F_node_y[idx_node_B]  + F_node_y[idx_node_T]) +
                    b_r*(F_node_y[idx_node_BB] + F_node_y[idx_node_TT])
                    );
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
        
        const int num_ghosts_0_flux_node = num_ghosts_flux_node[0];
        const int num_ghosts_1_flux_node = num_ghosts_flux_node[1];
        const int num_ghosts_2_flux_node = num_ghosts_flux_node[2];
        const int ghostcell_dim_0_flux_node = ghostcell_dims_flux_node[0];
        const int ghostcell_dim_1_flux_node = ghostcell_dims_flux_node[1];
        
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
                    
                    const int idx_node_BB = (i + num_ghosts_0_flux_node) +
                        (j - 2 + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node +
                        (k + num_ghosts_2_flux_node)*ghostcell_dim_0_flux_node*
                            ghostcell_dim_1_flux_node;
                    
                    const int idx_node_B = (i + num_ghosts_0_flux_node) +
                        (j - 1 + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node +
                        (k + num_ghosts_2_flux_node)*ghostcell_dim_0_flux_node*
                            ghostcell_dim_1_flux_node;
                    
                    const int idx_node_T = (i + num_ghosts_0_flux_node) +
                        (j + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node +
                        (k + num_ghosts_2_flux_node)*ghostcell_dim_0_flux_node*
                            ghostcell_dim_1_flux_node;
                    
                    const int idx_node_TT = (i + num_ghosts_0_flux_node) +
                        (j + 1 + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node +
                        (k + num_ghosts_2_flux_node)*ghostcell_dim_0_flux_node*
                            ghostcell_dim_1_flux_node;
                    
                    F_face_y[idx_face_y] += dt*(
                        a_r*(F_node_y[idx_node_B]  + F_node_y[idx_node_T]) +
                        b_r*(F_node_y[idx_node_BB] + F_node_y[idx_node_TT])
                        );
                }
            }
        }
    }
}


/*
 * Kernel to reconstruct the flux using flux at nodes in z-direction.
 */
void
DiffusiveFluxReconstructorNodeFourthOrder::reconstructFluxZ(
    double* F_face_z,
    const double* const F_node_z,
    const hier::IntVector& num_ghosts_flux_node,
    const hier::IntVector& ghostcell_dims_flux_node,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims,
    const hier::IntVector& interior_dims,
    const double& dt) const
{
    const double a_n =  double(2)/double(3);
    const double b_n = -double(1)/double(12);
    
    const double a_r = a_n + b_n;
    const double b_r = b_n;
    
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
    
    const int num_ghosts_0_flux_node = num_ghosts_flux_node[0];
    const int num_ghosts_1_flux_node = num_ghosts_flux_node[1];
    const int num_ghosts_2_flux_node = num_ghosts_flux_node[2];
    const int ghostcell_dim_0_flux_node = ghostcell_dims_flux_node[0];
    const int ghostcell_dim_1_flux_node = ghostcell_dims_flux_node[1];
    
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
                
                const int idx_node_BB = (i + num_ghosts_0_flux_node) +
                    (j + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node +
                    (k - 2 + num_ghosts_2_flux_node)*ghostcell_dim_0_flux_node*
                        ghostcell_dim_1_flux_node;
                
                const int idx_node_B = (i + num_ghosts_0_flux_node) +
                    (j + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node +
                    (k - 1 + num_ghosts_2_flux_node)*ghostcell_dim_0_flux_node*
                        ghostcell_dim_1_flux_node;
                
                const int idx_node_F = (i + num_ghosts_0_flux_node) +
                    (j + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node +
                    (k + num_ghosts_2_flux_node)*ghostcell_dim_0_flux_node*
                        ghostcell_dim_1_flux_node;
                
                const int idx_node_FF = (i + num_ghosts_0_flux_node) +
                    (j + num_ghosts_1_flux_node)*ghostcell_dim_0_flux_node +
                    (k + 1 + num_ghosts_2_flux_node)*ghostcell_dim_0_flux_node*
                        ghostcell_dim_1_flux_node;
                
                F_face_z[idx_face_z] += dt*(
                    a_r*(F_node_z[idx_node_B]  + F_node_z[idx_node_F]) +
                    b_r*(F_node_z[idx_node_BB] + F_node_z[idx_node_FF])
                    );
            }
        }
    }
}
