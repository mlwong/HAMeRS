#include "flow/diffusive_flux_reconstructors/midpoint/DiffusiveFluxReconstructorMidpointSixthOrder.hpp"

#include "SAMRAI/geom/CartesianPatchGeometry.h"

#include <map>

DiffusiveFluxReconstructorMidpointSixthOrder::DiffusiveFluxReconstructorMidpointSixthOrder(
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
    d_num_diff_ghosts = hier::IntVector::getOne(d_dim)*5;
    
    // The number of ghost cells used for finite difference and interpolation schemes at midpoints are the same.
    d_num_der_midpoint_ghosts     = hier::IntVector::getOne(d_dim)*3;
    d_num_der_node_ghosts         = hier::IntVector::getOne(d_dim)*3;
    d_num_interp_midpoint_ghosts  = hier::IntVector::getOne(d_dim)*3;
    d_num_flux_reconstruct_ghosts = hier::IntVector::getOne(d_dim)*3;
}


/*
 * Print all characteristics of the diffusive flux reconstruction class.
 */
void
DiffusiveFluxReconstructorMidpointSixthOrder::printClassData(
    std::ostream& os) const
{
    os << "\nPrint DiffusiveFluxReconstructorMidpointSixthOrder object..."
       << std::endl;
    
    os << std::endl;
    
    os << "DiffusiveFluxReconstructorMidpointSixthOrder: this = "
       << (DiffusiveFluxReconstructorMidpointSixthOrder *)this
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
DiffusiveFluxReconstructorMidpointSixthOrder::putToRestart(
   const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const
{
    restart_db->putString("d_diffusive_flux_reconstructor", "MIDPOINT_SIXTH_ORDER");
}


/*
 * Kernel to compute the derivatives in x-direction at midpoints.
 */
void
DiffusiveFluxReconstructorMidpointSixthOrder::computeFirstDerivativesInXAtMidpointX(
    Real* dudx,
    const Real* const u,
    const hier::IntVector& num_ghosts_derivative_midpoint,
    const hier::IntVector& num_ghosts_data_node,
    const hier::IntVector& ghostcell_dims_derivative_midpoint,
    const hier::IntVector& ghostcell_dims_data_node,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims,
    const Real& dx_0_inv) const
{
    const Real a_n =  Real(75)/Real(64);
    const Real b_n = -Real(25)/Real(384);
    const Real c_n =  Real(3)/Real(640);
    
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
            
            const int idx_data_LLL = i - 3 + num_ghosts_0_data;
            const int idx_data_LL  = i - 2 + num_ghosts_0_data;
            const int idx_data_L   = i - 1 + num_ghosts_0_data;
            const int idx_data_R   = i + 0 + num_ghosts_0_data;
            const int idx_data_RR  = i + 1 + num_ghosts_0_data;
            const int idx_data_RRR = i + 2 + num_ghosts_0_data;
            
            dudx[idx] = (a_n*(u[idx_data_R]   - u[idx_data_L]) +
                         b_n*(u[idx_data_RR]  - u[idx_data_LL]) +
                         c_n*(u[idx_data_RRR] - u[idx_data_LLL])
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
                
                const int idx_data_LLL = (i - 3 + num_ghosts_0_data) +
                    (j + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                const int idx_data_LL = (i - 2 + num_ghosts_0_data) +
                    (j + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                const int idx_data_L = (i - 1 + num_ghosts_0_data) +
                    (j + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                const int idx_data_R = (i + num_ghosts_0_data) +
                    (j + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                const int idx_data_RR = (i + 1 + num_ghosts_0_data) +
                    (j + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                const int idx_data_RRR = (i + 2 + num_ghosts_0_data) +
                    (j + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                dudx[idx] = (a_n*(u[idx_data_R]   - u[idx_data_L]) +
                             b_n*(u[idx_data_RR]  - u[idx_data_LL]) +
                             c_n*(u[idx_data_RRR] - u[idx_data_LLL])
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
                    
                    const int idx_data_LLL = (i - 3 + num_ghosts_0_data) +
                        (j + num_ghosts_1_data)*ghostcell_dim_0_data +
                        (k + num_ghosts_2_data)*ghostcell_dim_0_data*
                            ghostcell_dim_1_data;
                    
                    const int idx_data_LL = (i - 2 + num_ghosts_0_data) +
                        (j + num_ghosts_1_data)*ghostcell_dim_0_data +
                        (k + num_ghosts_2_data)*ghostcell_dim_0_data*
                            ghostcell_dim_1_data;
                    
                    const int idx_data_L = (i - 1 + num_ghosts_0_data) +
                        (j + num_ghosts_1_data)*ghostcell_dim_0_data +
                        (k + num_ghosts_2_data)*ghostcell_dim_0_data*
                            ghostcell_dim_1_data;
                    
                    const int idx_data_R = (i + num_ghosts_0_data) +
                        (j + num_ghosts_1_data)*ghostcell_dim_0_data +
                        (k + num_ghosts_2_data)*ghostcell_dim_0_data*
                            ghostcell_dim_1_data;
                    
                    const int idx_data_RR = (i + 1 + num_ghosts_0_data) +
                        (j + num_ghosts_1_data)*ghostcell_dim_0_data +
                        (k + num_ghosts_2_data)*ghostcell_dim_0_data*
                            ghostcell_dim_1_data;
                    
                    const int idx_data_RRR = (i + 2 + num_ghosts_0_data) +
                        (j + num_ghosts_1_data)*ghostcell_dim_0_data +
                        (k + num_ghosts_2_data)*ghostcell_dim_0_data*
                            ghostcell_dim_1_data;
                    
                    dudx[idx] = (a_n*(u[idx_data_R]   - u[idx_data_L]) +
                                 b_n*(u[idx_data_RR]  - u[idx_data_LL]) +
                                 c_n*(u[idx_data_RRR] - u[idx_data_LLL])
                                )*dx_0_inv;
                }
            }
        }
    }
}


/*
 * Kernel to compute the derivatives in y-direction at midpoints.
 */
void
DiffusiveFluxReconstructorMidpointSixthOrder::computeFirstDerivativesInYAtMidpointY(
    Real* dudy,
    const Real* const u,
    const hier::IntVector& num_ghosts_derivative_midpoint,
    const hier::IntVector& num_ghosts_data_node,
    const hier::IntVector& ghostcell_dims_derivative_midpoint,
    const hier::IntVector& ghostcell_dims_data_node,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims,
    const Real& dx_1_inv) const
{
    const Real a_n =  Real(75)/Real(64);
    const Real b_n = -Real(25)/Real(384);
    const Real c_n =  Real(3)/Real(640);
    
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
                
                const int idx_data_BBB = (i + num_ghosts_0_data) +
                    (j - 3 + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                const int idx_data_BB = (i + num_ghosts_0_data) +
                    (j - 2 + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                const int idx_data_B = (i + num_ghosts_0_data) +
                    (j - 1 + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                const int idx_data_T = (i + num_ghosts_0_data) +
                    (j + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                const int idx_data_TT = (i + num_ghosts_0_data) +
                    (j + 1 + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                const int idx_data_TTT = (i + num_ghosts_0_data) +
                    (j + 2 + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                dudy[idx] = (a_n*(u[idx_data_T]   - u[idx_data_B]) +
                             b_n*(u[idx_data_TT]  - u[idx_data_BB]) +
                             c_n*(u[idx_data_TTT] - u[idx_data_BBB])
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
                    
                    const int idx_data_BBB = (i + num_ghosts_0_data) +
                        (j - 3 + num_ghosts_1_data)*ghostcell_dim_0_data +
                        (k + num_ghosts_2_data)*ghostcell_dim_0_data*
                            ghostcell_dim_1_data;
                    
                    const int idx_data_BB = (i + num_ghosts_0_data) +
                        (j - 2 + num_ghosts_1_data)*ghostcell_dim_0_data +
                        (k + num_ghosts_2_data)*ghostcell_dim_0_data*
                            ghostcell_dim_1_data;
                    
                    const int idx_data_B = (i + num_ghosts_0_data) +
                        (j - 1 + num_ghosts_1_data)*ghostcell_dim_0_data +
                        (k + num_ghosts_2_data)*ghostcell_dim_0_data*
                            ghostcell_dim_1_data;
                    
                    const int idx_data_T = (i + num_ghosts_0_data) +
                        (j + num_ghosts_1_data)*ghostcell_dim_0_data +
                        (k + num_ghosts_2_data)*ghostcell_dim_0_data*
                            ghostcell_dim_1_data;
                    
                    const int idx_data_TT = (i + num_ghosts_0_data) +
                        (j + 1 + num_ghosts_1_data)*ghostcell_dim_0_data +
                        (k + num_ghosts_2_data)*ghostcell_dim_0_data*
                            ghostcell_dim_1_data;
                    
                    const int idx_data_TTT = (i + num_ghosts_0_data) +
                        (j + 2 + num_ghosts_1_data)*ghostcell_dim_0_data +
                        (k + num_ghosts_2_data)*ghostcell_dim_0_data*
                            ghostcell_dim_1_data;
                    
                    dudy[idx] = (a_n*(u[idx_data_T]   - u[idx_data_B]) +
                                 b_n*(u[idx_data_TT]  - u[idx_data_BB]) +
                                 c_n*(u[idx_data_TTT] - u[idx_data_BBB])
                                )*dx_1_inv;
                }
            }
        }
    }
}


/*
 * Kernel to compute the derivatives in z-direction at midpoints.
 */
void
DiffusiveFluxReconstructorMidpointSixthOrder::computeFirstDerivativesInZAtMidpointZ(
    Real* dudz,
    const Real* const u,
    const hier::IntVector& num_ghosts_derivative_midpoint,
    const hier::IntVector& num_ghosts_data_node,
    const hier::IntVector& ghostcell_dims_derivative_midpoint,
    const hier::IntVector& ghostcell_dims_data_node,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims,
    const Real& dx_2_inv) const
{
    const Real a_n =  Real(75)/Real(64);
    const Real b_n = -Real(25)/Real(384);
    const Real c_n =  Real(3)/Real(640);
    
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
                
                const int idx_data_BBB = (i + num_ghosts_0_data) +
                    (j + num_ghosts_1_data)*ghostcell_dim_0_data +
                    (k - 3 + num_ghosts_2_data)*ghostcell_dim_0_data*
                        ghostcell_dim_1_data;
                
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
                    (k + num_ghosts_2_data)*ghostcell_dim_0_data*
                        ghostcell_dim_1_data;
                
                const int idx_data_FF = (i + num_ghosts_0_data) +
                    (j + num_ghosts_1_data)*ghostcell_dim_0_data +
                    (k + 1 + num_ghosts_2_data)*ghostcell_dim_0_data*
                        ghostcell_dim_1_data;
                
                const int idx_data_FFF = (i + num_ghosts_0_data) +
                    (j + num_ghosts_1_data)*ghostcell_dim_0_data +
                    (k + 2 + num_ghosts_2_data)*ghostcell_dim_0_data*
                        ghostcell_dim_1_data;
                
                dudz[idx] = (a_n*(u[idx_data_F]   - u[idx_data_B]) +
                             b_n*(u[idx_data_FF]  - u[idx_data_BB]) +
                             c_n*(u[idx_data_FFF] - u[idx_data_BBB])
                            )*dx_2_inv;
            }
        }
    }
}


/*
 * Kernel to compute the derivatives in x-direction at nodes.
 */
void
DiffusiveFluxReconstructorMidpointSixthOrder::computeFirstDerivativesInXAtNode(
    Real* dudx,
    const Real* const u,
    const hier::IntVector& num_ghosts_derivative_node,
    const hier::IntVector& num_ghosts_data_node,
    const hier::IntVector& ghostcell_dims_derivative_node,
    const hier::IntVector& ghostcell_dims_data_node,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims,
    const Real& dx_0_inv) const
{
    const Real a_n =  Real(3)/Real(4);
    const Real b_n = -Real(3)/Real(20);
    const Real c_n =  Real(1)/Real(60);
    
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
            
            const int idx_data_LLL = i - 3 + num_ghosts_0_data;
            const int idx_data_LL  = i - 2 + num_ghosts_0_data;
            const int idx_data_L   = i - 1 + num_ghosts_0_data;
            const int idx_data_R   = i + 1 + num_ghosts_0_data;
            const int idx_data_RR  = i + 2 + num_ghosts_0_data;
            const int idx_data_RRR = i + 3 + num_ghosts_0_data;
            
            dudx[idx] = (a_n*(u[idx_data_R]   - u[idx_data_L]) +
                         b_n*(u[idx_data_RR]  - u[idx_data_LL]) +
                         c_n*(u[idx_data_RRR] - u[idx_data_LLL])
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
                
                const int idx_data_LLL = (i - 3 + num_ghosts_0_data) +
                    (j + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                const int idx_data_LL = (i - 2 + num_ghosts_0_data) +
                    (j + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                const int idx_data_L = (i - 1 + num_ghosts_0_data) +
                    (j + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                const int idx_data_R = (i + 1 + num_ghosts_0_data) +
                    (j + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                const int idx_data_RR = (i + 2 + num_ghosts_0_data) +
                    (j + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                const int idx_data_RRR = (i + 3 + num_ghosts_0_data) +
                    (j + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                dudx[idx] = (a_n*(u[idx_data_R]   - u[idx_data_L]) +
                             b_n*(u[idx_data_RR]  - u[idx_data_LL]) +
                             c_n*(u[idx_data_RRR] - u[idx_data_LLL])
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
                    
                    const int idx_data_LLL = (i - 3 + num_ghosts_0_data) +
                        (j + num_ghosts_1_data)*ghostcell_dim_0_data +
                        (k + num_ghosts_2_data)*ghostcell_dim_0_data*
                            ghostcell_dim_1_data;
                    
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
                    
                    const int idx_data_RRR = (i + 3 + num_ghosts_0_data) +
                        (j + num_ghosts_1_data)*ghostcell_dim_0_data +
                        (k + num_ghosts_2_data)*ghostcell_dim_0_data*
                            ghostcell_dim_1_data;
                    
                    dudx[idx] = (a_n*(u[idx_data_R]   - u[idx_data_L]) +
                                 b_n*(u[idx_data_RR]  - u[idx_data_LL]) +
                                 c_n*(u[idx_data_RRR] - u[idx_data_LLL])
                                )*dx_0_inv;
                }
            }
        }
    }
}


/*
 * Kernel to compute the derivatives in y-direction at nodes.
 */
void
DiffusiveFluxReconstructorMidpointSixthOrder::computeFirstDerivativesInYAtNode(
    Real* dudy,
    const Real* const u,
    const hier::IntVector& num_ghosts_derivative_node,
    const hier::IntVector& num_ghosts_data_node,
    const hier::IntVector& ghostcell_dims_derivative_node,
    const hier::IntVector& ghostcell_dims_data_node,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims,
    const Real& dx_1_inv) const
{
    const Real a_n =  Real(3)/Real(4);
    const Real b_n = -Real(3)/Real(20);
    const Real c_n =  Real(1)/Real(60);
    
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
                
                const int idx_data_BBB = (i + num_ghosts_0_data) +
                    (j - 3 + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                const int idx_data_BB = (i + num_ghosts_0_data) +
                    (j - 2 + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                const int idx_data_B = (i + num_ghosts_0_data) +
                    (j - 1 + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                const int idx_data_T = (i + num_ghosts_0_data) +
                    (j + 1 + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                const int idx_data_TT = (i + num_ghosts_0_data) +
                    (j + 2 + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                const int idx_data_TTT = (i + num_ghosts_0_data) +
                    (j + 3 + num_ghosts_1_data)*ghostcell_dim_0_data;
                
                dudy[idx] = (a_n*(u[idx_data_T]   - u[idx_data_B]) +
                             b_n*(u[idx_data_TT]  - u[idx_data_BB]) +
                             c_n*(u[idx_data_TTT] - u[idx_data_BBB])
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
                    
                    const int idx_data_BBB = (i + num_ghosts_0_data) +
                        (j - 3 + num_ghosts_1_data)*ghostcell_dim_0_data +
                        (k + num_ghosts_2_data)*ghostcell_dim_0_data*
                            ghostcell_dim_1_data;
                    
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
                    
                    const int idx_data_TTT = (i + num_ghosts_0_data) +
                        (j + 3 + num_ghosts_1_data)*ghostcell_dim_0_data +
                        (k + num_ghosts_2_data)*ghostcell_dim_0_data*
                            ghostcell_dim_1_data;
                    
                    dudy[idx] = (a_n*(u[idx_data_T]   - u[idx_data_B]) +
                                 b_n*(u[idx_data_TT]  - u[idx_data_BB]) +
                                 c_n*(u[idx_data_TTT] - u[idx_data_BBB])
                                )*dx_1_inv;
                }
            }
        }
    }
}


/*
 * Kernel to compute the derivatives in z-direction at nodes.
 */
void
DiffusiveFluxReconstructorMidpointSixthOrder::computeFirstDerivativesInZAtNode(
    Real* dudz,
    const Real* const u,
    const hier::IntVector& num_ghosts_derivative_node,
    const hier::IntVector& num_ghosts_data_node,
    const hier::IntVector& ghostcell_dims_derivative_node,
    const hier::IntVector& ghostcell_dims_data_node,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims,
    const Real& dx_2_inv) const
{
    const Real a_n =  Real(3)/Real(4);
    const Real b_n = -Real(3)/Real(20);
    const Real c_n =  Real(1)/Real(60);
    
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
                
                const int idx_data_BBB = (i + num_ghosts_0_data) +
                    (j + num_ghosts_1_data)*ghostcell_dim_0_data +
                    (k - 3 + num_ghosts_2_data)*ghostcell_dim_0_data*
                        ghostcell_dim_1_data;
                
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
                
                const int idx_data_FFF = (i + num_ghosts_0_data) +
                    (j + num_ghosts_1_data)*ghostcell_dim_0_data +
                    (k + 3 + num_ghosts_2_data)*ghostcell_dim_0_data*
                        ghostcell_dim_1_data;
                
                dudz[idx] = (a_n*(u[idx_data_F]   - u[idx_data_B]) +
                             b_n*(u[idx_data_FF]  - u[idx_data_BB]) +
                             c_n*(u[idx_data_FFF] - u[idx_data_BBB])
                            )*dx_2_inv;
            }
        }
    }
}


/*
 * Kernel to interpolate the data from nodes to midpoints in x-direction.
 */
void
DiffusiveFluxReconstructorMidpointSixthOrder::interpolateDataFromNodeToMidpointX(
    Real* u_midpoint_x,
    const Real* const u_node,
    const hier::IntVector& num_ghosts_data_midpoint,
    const hier::IntVector& num_ghosts_data_node,
    const hier::IntVector& ghostcell_dims_data_midpoint,
    const hier::IntVector& ghostcell_dims_data_node,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims) const
{
    const Real a_n =  Real(75)/Real(128);
    const Real b_n = -Real(25)/Real(256);
    const Real c_n =  Real(3)/Real(256);
    
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
            
            const int idx_node_LLL = i - 3 + num_ghosts_0_data_node;
            const int idx_node_LL  = i - 2 + num_ghosts_0_data_node;
            const int idx_node_L   = i - 1 + num_ghosts_0_data_node;
            const int idx_node_R   = i + 0 + num_ghosts_0_data_node;
            const int idx_node_RR  = i + 1 + num_ghosts_0_data_node;
            const int idx_node_RRR = i + 2 + num_ghosts_0_data_node;
            
            u_midpoint_x[idx] = (a_n*(u_node[idx_node_R]   + u_node[idx_node_L]) +
                                 b_n*(u_node[idx_node_RR]  + u_node[idx_node_LL]) +
                                 c_n*(u_node[idx_node_RRR] + u_node[idx_node_LLL])
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
                
                const int idx_node_LLL = (i - 3 + num_ghosts_0_data_node) +
                    (j + num_ghosts_1_data_node)*ghostcell_dim_0_data_node;
                
                const int idx_node_LL = (i - 2 + num_ghosts_0_data_node) +
                    (j + num_ghosts_1_data_node)*ghostcell_dim_0_data_node;
                
                const int idx_node_L = (i - 1 + num_ghosts_0_data_node) +
                    (j + num_ghosts_1_data_node)*ghostcell_dim_0_data_node;
                
                const int idx_node_R = (i + num_ghosts_0_data_node) +
                    (j + num_ghosts_1_data_node)*ghostcell_dim_0_data_node;
                
                const int idx_node_RR = (i + 1 + num_ghosts_0_data_node) +
                    (j + num_ghosts_1_data_node)*ghostcell_dim_0_data_node;
                
                const int idx_node_RRR = (i + 2 + num_ghosts_0_data_node) +
                    (j + num_ghosts_1_data_node)*ghostcell_dim_0_data_node;
                
                u_midpoint_x[idx] = (a_n*(u_node[idx_node_R]   + u_node[idx_node_L]) +
                                     b_n*(u_node[idx_node_RR]  + u_node[idx_node_LL]) +
                                     c_n*(u_node[idx_node_RRR] + u_node[idx_node_LLL])
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
                    
                    const int idx_node_LLL = (i - 3 + num_ghosts_0_data_node) +
                        (j + num_ghosts_1_data_node)*ghostcell_dim_0_data_node +
                        (k + num_ghosts_2_data_node)*ghostcell_dim_0_data_node*
                            ghostcell_dim_1_data_node;
                    
                    const int idx_node_LL = (i - 2 + num_ghosts_0_data_node) +
                        (j + num_ghosts_1_data_node)*ghostcell_dim_0_data_node +
                        (k + num_ghosts_2_data_node)*ghostcell_dim_0_data_node*
                            ghostcell_dim_1_data_node;
                    
                    const int idx_node_L = (i - 1 + num_ghosts_0_data_node) +
                        (j + num_ghosts_1_data_node)*ghostcell_dim_0_data_node +
                        (k + num_ghosts_2_data_node)*ghostcell_dim_0_data_node*
                            ghostcell_dim_1_data_node;
                    
                    const int idx_node_R = (i + num_ghosts_0_data_node) +
                        (j + num_ghosts_1_data_node)*ghostcell_dim_0_data_node +
                        (k + num_ghosts_2_data_node)*ghostcell_dim_0_data_node*
                            ghostcell_dim_1_data_node;
                    
                    const int idx_node_RR = (i + 1 + num_ghosts_0_data_node) +
                        (j + num_ghosts_1_data_node)*ghostcell_dim_0_data_node +
                        (k + num_ghosts_2_data_node)*ghostcell_dim_0_data_node*
                            ghostcell_dim_1_data_node;
                    
                    const int idx_node_RRR = (i + 2 + num_ghosts_0_data_node) +
                        (j + num_ghosts_1_data_node)*ghostcell_dim_0_data_node +
                        (k + num_ghosts_2_data_node)*ghostcell_dim_0_data_node*
                            ghostcell_dim_1_data_node;
                    
                    u_midpoint_x[idx] = (a_n*(u_node[idx_node_R]   + u_node[idx_node_L]) +
                                         b_n*(u_node[idx_node_RR]  + u_node[idx_node_LL]) +
                                         c_n*(u_node[idx_node_RRR] + u_node[idx_node_LLL])
                                        );
                }
            }
        }
    }
}


/*
 * Kernel to interpolate the data from nodes to midpoints in y-direction.
 */
void
DiffusiveFluxReconstructorMidpointSixthOrder::interpolateDataFromNodeToMidpointY(
    Real* u_midpoint_y,
    const Real* const u_node,
    const hier::IntVector& num_ghosts_data_midpoint,
    const hier::IntVector& num_ghosts_data_node,
    const hier::IntVector& ghostcell_dims_data_midpoint,
    const hier::IntVector& ghostcell_dims_data_node,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims) const
{
    const Real a_n =  Real(75)/Real(128);
    const Real b_n = -Real(25)/Real(256);
    const Real c_n =  Real(3)/Real(256);
    
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
                
                const int idx_node_BBB = (i + num_ghosts_0_data_node) +
                    (j - 3 + num_ghosts_1_data_node)*ghostcell_dim_0_data_node;
                
                const int idx_node_BB = (i + num_ghosts_0_data_node) +
                    (j - 2 + num_ghosts_1_data_node)*ghostcell_dim_0_data_node;
                
                const int idx_node_B = (i + num_ghosts_0_data_node) +
                    (j - 1 + num_ghosts_1_data_node)*ghostcell_dim_0_data_node;
                
                const int idx_node_T = (i + num_ghosts_0_data_node) +
                    (j + num_ghosts_1_data_node)*ghostcell_dim_0_data_node;
                
                const int idx_node_TT = (i + num_ghosts_0_data_node) +
                    (j + 1 + num_ghosts_1_data_node)*ghostcell_dim_0_data_node;
                
                const int idx_node_TTT = (i + num_ghosts_0_data_node) +
                    (j + 2 + num_ghosts_1_data_node)*ghostcell_dim_0_data_node;
                
                u_midpoint_y[idx] = (a_n*(u_node[idx_node_T]   + u_node[idx_node_B]) +
                                     b_n*(u_node[idx_node_TT]  + u_node[idx_node_BB]) +
                                     c_n*(u_node[idx_node_TTT] + u_node[idx_node_BBB])
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
                    
                    const int idx_node_BBB = (i + num_ghosts_0_data_node) +
                        (j - 3 + num_ghosts_1_data_node)*ghostcell_dim_0_data_node +
                        (k + num_ghosts_2_data_node)*ghostcell_dim_0_data_node*
                            ghostcell_dim_1_data_node;
                    
                    const int idx_node_BB = (i + num_ghosts_0_data_node) +
                        (j - 2 + num_ghosts_1_data_node)*ghostcell_dim_0_data_node +
                        (k + num_ghosts_2_data_node)*ghostcell_dim_0_data_node*
                            ghostcell_dim_1_data_node;
                    
                    const int idx_node_B = (i + num_ghosts_0_data_node) +
                        (j - 1 + num_ghosts_1_data_node)*ghostcell_dim_0_data_node +
                        (k + num_ghosts_2_data_node)*ghostcell_dim_0_data_node*
                            ghostcell_dim_1_data_node;
                    
                    const int idx_node_T = (i + num_ghosts_0_data_node) +
                        (j + num_ghosts_1_data_node)*ghostcell_dim_0_data_node +
                        (k + num_ghosts_2_data_node)*ghostcell_dim_0_data_node*
                            ghostcell_dim_1_data_node;
                    
                    const int idx_node_TT = (i + num_ghosts_0_data_node) +
                        (j + 1 + num_ghosts_1_data_node)*ghostcell_dim_0_data_node +
                        (k + num_ghosts_2_data_node)*ghostcell_dim_0_data_node*
                            ghostcell_dim_1_data_node;
                    
                    const int idx_node_TTT = (i + num_ghosts_0_data_node) +
                        (j + 2 + num_ghosts_1_data_node)*ghostcell_dim_0_data_node +
                        (k + num_ghosts_2_data_node)*ghostcell_dim_0_data_node*
                            ghostcell_dim_1_data_node;
                    
                    u_midpoint_y[idx] = (a_n*(u_node[idx_node_T]   + u_node[idx_node_B]) +
                                         b_n*(u_node[idx_node_TT]  + u_node[idx_node_BB]) +
                                         c_n*(u_node[idx_node_TTT] + u_node[idx_node_BBB])
                                        );
                }
            }
        }
    }
}


/*
 * Kernel to interpolate the data from nodes to midpoints in z-direction.
 */
void
DiffusiveFluxReconstructorMidpointSixthOrder::interpolateDataFromNodeToMidpointZ(
    Real* u_midpoint_z,
    const Real* const u_node,
    const hier::IntVector& num_ghosts_data_midpoint,
    const hier::IntVector& num_ghosts_data_node,
    const hier::IntVector& ghostcell_dims_data_midpoint,
    const hier::IntVector& ghostcell_dims_data_node,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims) const
{
    const Real a_n =  Real(75)/Real(128);
    const Real b_n = -Real(25)/Real(256);
    const Real c_n =  Real(3)/Real(256);
    
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
                
                const int idx_node_BBB = (i + num_ghosts_0_data_node) +
                    (j + num_ghosts_1_data_node)*ghostcell_dim_0_data_node +
                    (k - 3 + num_ghosts_2_data_node)*ghostcell_dim_0_data_node*
                        ghostcell_dim_1_data_node;
                
                const int idx_node_BB = (i + num_ghosts_0_data_node) +
                    (j + num_ghosts_1_data_node)*ghostcell_dim_0_data_node +
                    (k - 2 + num_ghosts_2_data_node)*ghostcell_dim_0_data_node*
                        ghostcell_dim_1_data_node;
                
                const int idx_node_B = (i + num_ghosts_0_data_node) +
                    (j + num_ghosts_1_data_node)*ghostcell_dim_0_data_node +
                    (k - 1 + num_ghosts_2_data_node)*ghostcell_dim_0_data_node*
                        ghostcell_dim_1_data_node;
                
                const int idx_node_F = (i + num_ghosts_0_data_node) +
                    (j + num_ghosts_1_data_node)*ghostcell_dim_0_data_node +
                    (k + num_ghosts_2_data_node)*ghostcell_dim_0_data_node*
                        ghostcell_dim_1_data_node;
                
                const int idx_node_FF = (i + num_ghosts_0_data_node) +
                    (j + num_ghosts_1_data_node)*ghostcell_dim_0_data_node +
                    (k + 1 + num_ghosts_2_data_node)*ghostcell_dim_0_data_node*
                        ghostcell_dim_1_data_node;
                
                const int idx_node_FFF = (i + num_ghosts_0_data_node) +
                    (j + num_ghosts_1_data_node)*ghostcell_dim_0_data_node +
                    (k + 2 + num_ghosts_2_data_node)*ghostcell_dim_0_data_node*
                        ghostcell_dim_1_data_node;
                
                u_midpoint_z[idx] = (a_n*(u_node[idx_node_F]   + u_node[idx_node_B]) +
                                     b_n*(u_node[idx_node_FF]  + u_node[idx_node_BB]) +
                                     c_n*(u_node[idx_node_FFF] + u_node[idx_node_BBB])
                                    );
            }
        }
    }
}


/*
 * Kernel to reconstruct the flux using flux at midpoints in x-direction.
 */
void
DiffusiveFluxReconstructorMidpointSixthOrder::reconstructFluxX(
    Real* F_face_x,
    const Real* const F_midpoint_x,
    const hier::IntVector& num_ghosts_flux_midpoint,
    const hier::IntVector& ghostcell_dims_flux_midpoint,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims,
    const hier::IntVector& interior_dims,
    const double& dt) const
{
    const Real a_m =  Real(75)/Real(64);
    const Real b_m = -Real(25)/Real(384);
    const Real c_m =  Real(3)/Real(640);
    
    const Real a_r = a_m + b_m + c_m;
    const Real b_r = b_m + c_m;
    const Real c_r = c_m;
    
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
            
            const int idx_midpoint_LL = i - 2 + num_ghosts_0_flux_midpoint;
            const int idx_midpoint_L  = i - 1 + num_ghosts_0_flux_midpoint;
            const int idx_midpoint    = i + 0 + num_ghosts_0_flux_midpoint;
            const int idx_midpoint_R  = i + 1 + num_ghosts_0_flux_midpoint;
            const int idx_midpoint_RR = i + 2 + num_ghosts_0_flux_midpoint;
            
            F_face_x[idx_face_x] += Real(dt)*(
                a_r*(F_midpoint_x[idx_midpoint]) +
                b_r*(F_midpoint_x[idx_midpoint_L]  + F_midpoint_x[idx_midpoint_R]) +
                c_r*(F_midpoint_x[idx_midpoint_LL] + F_midpoint_x[idx_midpoint_RR])
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
                
                const int idx_midpoint_LL = (i - 2 + num_ghosts_0_flux_midpoint) +
                    (j + num_ghosts_1_flux_midpoint)*(ghostcell_dim_0_flux_midpoint + 1);
                
                const int idx_midpoint_L = (i - 1 + num_ghosts_0_flux_midpoint) +
                    (j + num_ghosts_1_flux_midpoint)*(ghostcell_dim_0_flux_midpoint + 1);
                
                const int idx_midpoint = (i + num_ghosts_0_flux_midpoint) +
                    (j + num_ghosts_1_flux_midpoint)*(ghostcell_dim_0_flux_midpoint + 1);
                
                const int idx_midpoint_R = (i + 1 + num_ghosts_0_flux_midpoint) +
                    (j + num_ghosts_1_flux_midpoint)*(ghostcell_dim_0_flux_midpoint + 1);
                
                const int idx_midpoint_RR = (i + 2 + num_ghosts_0_flux_midpoint) +
                    (j + num_ghosts_1_flux_midpoint)*(ghostcell_dim_0_flux_midpoint + 1);
                
                F_face_x[idx_face_x] += Real(dt)*(
                    a_r*(F_midpoint_x[idx_midpoint]) +
                    b_r*(F_midpoint_x[idx_midpoint_L]  + F_midpoint_x[idx_midpoint_R]) +
                    c_r*(F_midpoint_x[idx_midpoint_LL] + F_midpoint_x[idx_midpoint_RR])
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
                    
                    const int idx_midpoint_LL = (i - 2 + num_ghosts_0_flux_midpoint) +
                        (j + num_ghosts_1_flux_midpoint)*(ghostcell_dim_0_flux_midpoint + 1) +
                        (k + num_ghosts_2_flux_midpoint)*(ghostcell_dim_0_flux_midpoint + 1)*
                            ghostcell_dim_1_flux_midpoint;
                    
                    const int idx_midpoint_L = (i - 1 + num_ghosts_0_flux_midpoint) +
                        (j + num_ghosts_1_flux_midpoint)*(ghostcell_dim_0_flux_midpoint + 1) +
                        (k + num_ghosts_2_flux_midpoint)*(ghostcell_dim_0_flux_midpoint + 1)*
                            ghostcell_dim_1_flux_midpoint;
                    
                    const int idx_midpoint = (i + num_ghosts_0_flux_midpoint) +
                        (j + num_ghosts_1_flux_midpoint)*(ghostcell_dim_0_flux_midpoint + 1) +
                        (k + num_ghosts_2_flux_midpoint)*(ghostcell_dim_0_flux_midpoint + 1)*
                            ghostcell_dim_1_flux_midpoint;
                    
                    const int idx_midpoint_R = (i + 1 + num_ghosts_0_flux_midpoint) +
                        (j + num_ghosts_1_flux_midpoint)*(ghostcell_dim_0_flux_midpoint + 1) +
                        (k + num_ghosts_2_flux_midpoint)*(ghostcell_dim_0_flux_midpoint + 1)*
                            ghostcell_dim_1_flux_midpoint;
                    
                    const int idx_midpoint_RR = (i + 2 + num_ghosts_0_flux_midpoint) +
                        (j + num_ghosts_1_flux_midpoint)*(ghostcell_dim_0_flux_midpoint + 1) +
                        (k + num_ghosts_2_flux_midpoint)*(ghostcell_dim_0_flux_midpoint + 1)*
                            ghostcell_dim_1_flux_midpoint;
                    
                    F_face_x[idx_face_x] += Real(dt)*(
                        a_r*(F_midpoint_x[idx_midpoint]) +
                        b_r*(F_midpoint_x[idx_midpoint_L]  + F_midpoint_x[idx_midpoint_R]) +
                        c_r*(F_midpoint_x[idx_midpoint_LL] + F_midpoint_x[idx_midpoint_RR])
                        );
                }
            }
        }
    }
}


/*
 * Kernel to reconstruct the flux using flux at midpoints in y-direction.
 */
void
DiffusiveFluxReconstructorMidpointSixthOrder::reconstructFluxY(
    Real* F_face_y,
    const Real* const F_midpoint_y,
    const hier::IntVector& num_ghosts_flux_midpoint,
    const hier::IntVector& ghostcell_dims_flux_midpoint,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims,
    const hier::IntVector& interior_dims,
    const double& dt) const
{
    const Real a_m =  Real(75)/Real(64);
    const Real b_m = -Real(25)/Real(384);
    const Real c_m =  Real(3)/Real(640);
    
    const Real a_r = a_m + b_m + c_m;
    const Real b_r = b_m + c_m;
    const Real c_r = c_m;
    
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
                
                const int idx_midpoint_BB = (i + num_ghosts_0_flux_midpoint) +
                    (j - 2 + num_ghosts_1_flux_midpoint)*ghostcell_dim_0_flux_midpoint;
                
                const int idx_midpoint_B = (i + num_ghosts_0_flux_midpoint) +
                    (j - 1 + num_ghosts_1_flux_midpoint)*ghostcell_dim_0_flux_midpoint;
                
                const int idx_midpoint = (i + num_ghosts_0_flux_midpoint) +
                    (j + num_ghosts_1_flux_midpoint)*ghostcell_dim_0_flux_midpoint;
                
                const int idx_midpoint_T = (i + num_ghosts_0_flux_midpoint) +
                    (j + 1 + num_ghosts_1_flux_midpoint)*ghostcell_dim_0_flux_midpoint;
                
                const int idx_midpoint_TT = (i + num_ghosts_0_flux_midpoint) +
                    (j + 2 + num_ghosts_1_flux_midpoint)*ghostcell_dim_0_flux_midpoint;
                
                F_face_y[idx_face_y] += Real(dt)*(
                    a_r*(F_midpoint_y[idx_midpoint]) +
                    b_r*(F_midpoint_y[idx_midpoint_B]  + F_midpoint_y[idx_midpoint_T]) +
                    c_r*(F_midpoint_y[idx_midpoint_BB] + F_midpoint_y[idx_midpoint_TT])
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
                    
                    const int idx_midpoint_BB = (i + num_ghosts_0_flux_midpoint) +
                        (j - 2 + num_ghosts_1_flux_midpoint)*ghostcell_dim_0_flux_midpoint +
                        (k + num_ghosts_2_flux_midpoint)*ghostcell_dim_0_flux_midpoint*
                            (ghostcell_dim_1_flux_midpoint + 1);
                    
                    const int idx_midpoint_B = (i + num_ghosts_0_flux_midpoint) +
                        (j - 1 + num_ghosts_1_flux_midpoint)*ghostcell_dim_0_flux_midpoint +
                        (k + num_ghosts_2_flux_midpoint)*ghostcell_dim_0_flux_midpoint*
                            (ghostcell_dim_1_flux_midpoint + 1);
                    
                    const int idx_midpoint = (i + num_ghosts_0_flux_midpoint) +
                        (j + num_ghosts_1_flux_midpoint)*ghostcell_dim_0_flux_midpoint +
                        (k + num_ghosts_2_flux_midpoint)*ghostcell_dim_0_flux_midpoint*
                            (ghostcell_dim_1_flux_midpoint + 1);
                    
                    const int idx_midpoint_T = (i + num_ghosts_0_flux_midpoint) +
                        (j + 1 + num_ghosts_1_flux_midpoint)*ghostcell_dim_0_flux_midpoint +
                        (k + num_ghosts_2_flux_midpoint)*ghostcell_dim_0_flux_midpoint*
                            (ghostcell_dim_1_flux_midpoint + 1);
                    
                    const int idx_midpoint_TT = (i + num_ghosts_0_flux_midpoint) +
                        (j + 2 + num_ghosts_1_flux_midpoint)*ghostcell_dim_0_flux_midpoint +
                        (k + num_ghosts_2_flux_midpoint)*ghostcell_dim_0_flux_midpoint*
                            (ghostcell_dim_1_flux_midpoint + 1);
                    
                    F_face_y[idx_face_y] += Real(dt)*(
                        a_r*(F_midpoint_y[idx_midpoint]) +
                        b_r*(F_midpoint_y[idx_midpoint_B]  + F_midpoint_y[idx_midpoint_T]) +
                        c_r*(F_midpoint_y[idx_midpoint_BB] + F_midpoint_y[idx_midpoint_TT])
                        );
                }
            }
        }
    }
}


/*
 * Kernel to reconstruct the flux using flux at midpoints in z-direction.
 */
void
DiffusiveFluxReconstructorMidpointSixthOrder::reconstructFluxZ(
    Real* F_face_z,
    const Real* const F_midpoint_z,
    const hier::IntVector& num_ghosts_flux_midpoint,
    const hier::IntVector& ghostcell_dims_flux_midpoint,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims,
    const hier::IntVector& interior_dims,
    const double& dt) const
{
    const Real a_m =  Real(75)/Real(64);
    const Real b_m = -Real(25)/Real(384);
    const Real c_m =  Real(3)/Real(640);
    
    const Real a_r = a_m + b_m + c_m;
    const Real b_r = b_m + c_m;
    const Real c_r = c_m;
    
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
                
                const int idx_midpoint_BB = (i + num_ghosts_0_flux_midpoint) +
                    (j + num_ghosts_1_flux_midpoint)*ghostcell_dim_0_flux_midpoint +
                    (k - 2 + num_ghosts_2_flux_midpoint)*ghostcell_dim_0_flux_midpoint*
                        ghostcell_dim_1_flux_midpoint;
                
                const int idx_midpoint_B = (i + num_ghosts_0_flux_midpoint) +
                    (j + num_ghosts_1_flux_midpoint)*ghostcell_dim_0_flux_midpoint +
                    (k - 1 + num_ghosts_2_flux_midpoint)*ghostcell_dim_0_flux_midpoint*
                        ghostcell_dim_1_flux_midpoint;
                
                const int idx_midpoint = (i + num_ghosts_0_flux_midpoint) +
                    (j + num_ghosts_1_flux_midpoint)*ghostcell_dim_0_flux_midpoint +
                    (k + num_ghosts_2_flux_midpoint)*ghostcell_dim_0_flux_midpoint*
                        ghostcell_dim_1_flux_midpoint;
                
                const int idx_midpoint_F = (i + num_ghosts_0_flux_midpoint) +
                    (j + num_ghosts_1_flux_midpoint)*ghostcell_dim_0_flux_midpoint +
                    (k + 1 + num_ghosts_2_flux_midpoint)*ghostcell_dim_0_flux_midpoint*
                        ghostcell_dim_1_flux_midpoint;
                
                const int idx_midpoint_FF = (i + num_ghosts_0_flux_midpoint) +
                    (j + num_ghosts_1_flux_midpoint)*ghostcell_dim_0_flux_midpoint +
                    (k + 2 + num_ghosts_2_flux_midpoint)*ghostcell_dim_0_flux_midpoint*
                        ghostcell_dim_1_flux_midpoint;
                
                F_face_z[idx_face_z] += Real(dt)*(
                    a_r*(F_midpoint_z[idx_midpoint]) +
                    b_r*(F_midpoint_z[idx_midpoint_B]  + F_midpoint_z[idx_midpoint_F]) +
                    c_r*(F_midpoint_z[idx_midpoint_BB] + F_midpoint_z[idx_midpoint_FF])
                    );
            }
        }
    }
}
