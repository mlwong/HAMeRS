#include "flow/convective_flux_reconstructors/WCNS56/ConvectiveFluxReconstructorWCNS6-CU-M2-HLLC-HLL.hpp"

#include <cfloat>

#include "SAMRAI/geom/CartesianPatchGeometry.h"


/*
 * Interger based power function.
 */
static inline __attribute__((always_inline)) Real ipow(Real base, const int& exp)
{
    Real result = base;
    for (int i = 1; i < exp; i++)
    {
        result *= base;
    }

    return result;
}


/*
 * Compute local beta's.
 */
static inline __attribute__((always_inline)) void computeLocalBeta(
    Real& beta_0,
    Real& beta_1,
    Real& beta_2,
    Real& beta_3,
    Real** U_array,
    const int& idx_side)
{
    beta_0 = Real(1)/Real(3)*(U_array[0][idx_side]*(Real(4)*U_array[0][idx_side] -
         Real(19)*U_array[1][idx_side] + Real(11)*U_array[2][idx_side]) +
         U_array[1][idx_side]*(Real(25)*U_array[1][idx_side] - Real(31)*U_array[2][idx_side]) +
         Real(10)*U_array[2][idx_side]*U_array[2][idx_side]);
    
    beta_1 = Real(1)/Real(3)*(U_array[1][idx_side]*(Real(4)*U_array[1][idx_side] -
         Real(13)*U_array[2][idx_side] + Real(5)*U_array[3][idx_side]) +
         Real(13)*U_array[2][idx_side]*(U_array[2][idx_side] - U_array[3][idx_side]) +
         Real(4)*U_array[3][idx_side]*U_array[3][idx_side]);
    
    beta_2 = Real(1)/Real(3)*(U_array[2][idx_side]*(Real(10)*U_array[2][idx_side] -
         Real(31)*U_array[3][idx_side] + Real(11)*U_array[4][idx_side]) +
         U_array[3][idx_side]*(Real(25)*U_array[3][idx_side] - Real(19)*U_array[4][idx_side]) +
         Real(4)*U_array[4][idx_side]*U_array[4][idx_side]);
    
    beta_3 = Real(1)/Real(232243200)*(U_array[0][idx_side]*(Real(525910327)*U_array[0][idx_side] -
         Real(4562164630)*U_array[1][idx_side] + Real(7799501420)*U_array[2][idx_side] -
         Real(6610694540)*U_array[3][idx_side] + Real(2794296070)*U_array[4][idx_side] -
         Real(472758974)*U_array[5][idx_side]) + Real(5)*U_array[1][idx_side]*
        (Real(2146987907)*U_array[1][idx_side] - Real(7722406988)*U_array[2][idx_side] +
         Real(6763559276)*U_array[3][idx_side] - Real(2926461814)*U_array[4][idx_side] +
         Real(503766638)*U_array[5][idx_side]) + Real(20)*U_array[2][idx_side]*
        (Real(1833221603)*U_array[2][idx_side] - Real(3358664662)*U_array[3][idx_side] +
         Real(1495974539)*U_array[4][idx_side] - Real(263126407)*U_array[5][idx_side]) +
        Real(20)*U_array[3][idx_side]*(Real(1607794163)*U_array[3][idx_side] -
         Real(1486026707)*U_array[4][idx_side] + Real(268747951)*U_array[5][idx_side]) +
        Real(5)*U_array[4][idx_side]*(Real(1432381427)*U_array[4][idx_side] -
         Real(536951582)*U_array[5][idx_side]) +
        Real(263126407)*U_array[5][idx_side]*U_array[5][idx_side]);
}


/*
 * Compute local beta_tilde's.
 */
static inline __attribute__((always_inline)) void computeLocalBetaTilde(
    Real& beta_tilde_0,
    Real& beta_tilde_1,
    Real& beta_tilde_2,
    Real& beta_tilde_3,
    Real** U_array,
    const int& idx_side)
{
    beta_tilde_0 = Real(1)/Real(3)*(U_array[5][idx_side]*(Real(4)*U_array[5][idx_side] -
         Real(19)*U_array[4][idx_side] + Real(11)*U_array[3][idx_side]) +
         U_array[4][idx_side]*(Real(25)*U_array[4][idx_side] - Real(31)*U_array[3][idx_side]) +
         Real(10)*U_array[3][idx_side]*U_array[3][idx_side]);
    
    beta_tilde_1 = Real(1)/Real(3)*(U_array[4][idx_side]*(Real(4)*U_array[4][idx_side] -
         Real(13)*U_array[3][idx_side] + Real(5)*U_array[2][idx_side]) +
         Real(13)*U_array[3][idx_side]*(U_array[3][idx_side] - U_array[2][idx_side]) +
         Real(4)*U_array[2][idx_side]*U_array[2][idx_side]);
    
    beta_tilde_2 = Real(1)/Real(3)*(U_array[3][idx_side]*(Real(10)*U_array[3][idx_side] -
         Real(31)*U_array[2][idx_side] + Real(11)*U_array[1][idx_side]) +
         U_array[2][idx_side]*(Real(25)*U_array[2][idx_side] - Real(19)*U_array[1][idx_side]) +
         Real(4)*U_array[1][idx_side]*U_array[1][idx_side]);
    
    beta_tilde_3 = Real(1)/Real(232243200)*(U_array[5][idx_side]*(Real(525910327)*U_array[5][idx_side] -
         Real(4562164630)*U_array[4][idx_side] + Real(7799501420)*U_array[3][idx_side] -
         Real(6610694540)*U_array[2][idx_side] + Real(2794296070)*U_array[1][idx_side] -
         Real(472758974)*U_array[0][idx_side]) + Real(5)*U_array[4][idx_side]*
        (Real(2146987907)*U_array[4][idx_side] - Real(7722406988)*U_array[3][idx_side] +
         Real(6763559276)*U_array[2][idx_side] - Real(2926461814)*U_array[1][idx_side] +
         Real(503766638)*U_array[0][idx_side]) + Real(20)*U_array[3][idx_side]*
        (Real(1833221603)*U_array[3][idx_side] - Real(3358664662)*U_array[2][idx_side] +
         Real(1495974539)*U_array[1][idx_side] - Real(263126407)*U_array[0][idx_side]) +
        Real(20)*U_array[2][idx_side]*(Real(1607794163)*U_array[2][idx_side] -
         Real(1486026707)*U_array[1][idx_side] + Real(268747951)*U_array[0][idx_side]) +
        Real(5)*U_array[1][idx_side]*(Real(1432381427)*U_array[1][idx_side] -
         Real(536951582)*U_array[0][idx_side]) +
        Real(263126407)*U_array[0][idx_side]*U_array[0][idx_side]);
}


/*
 * Perform local WENO interpolation of U_minus.
 */
static inline __attribute__((always_inline)) void performLocalWENOInterpolationMinus(
    Real* U_minus,
    Real** U_array,
    const int& idx_side,
    const int& q,
    const Real& C,
    const Real& Chi,
    const Real& epsilon,
    const Real& dx)
{
    /*
     * Compute beta's.
     */
    
    Real beta_0, beta_1, beta_2, beta_3;
    
    computeLocalBeta(beta_0, beta_1, beta_2, beta_3, U_array, idx_side);
    
    /*
     * Compute the weights omega.
     */
    
    Real omega_0, omega_1, omega_2, omega_3;
    
    Real beta_avg = Real(1)/Real(8)*(beta_0 + beta_2 + 6*beta_1);
    Real tau_6 = std::abs(beta_3 - beta_avg);
    
    omega_0 = Real(1)/Real(32)*ipow(C + tau_6/(beta_0 + epsilon*dx*dx)*
        (beta_avg + Chi*dx*dx)/(beta_0 + Chi*dx*dx), q);
    omega_1 = Real(15)/Real(32)*ipow(C + tau_6/(beta_1 + epsilon*dx*dx)*
        (beta_avg + Chi*dx*dx)/(beta_1 + Chi*dx*dx), q);
    omega_2 = Real(15)/Real(32)*ipow(C + tau_6/(beta_2 + epsilon*dx*dx)*
        (beta_avg + Chi*dx*dx)/(beta_2 + Chi*dx*dx), q);
    omega_3 = Real(1)/Real(32)*ipow(C + tau_6/(beta_3 + epsilon*dx*dx)*
        (beta_avg + Chi*dx*dx)/(beta_3 + Chi*dx*dx), q);
    
    Real omega_sum = omega_0 + omega_1 + omega_2 + omega_3;
    
    omega_0 = omega_0/omega_sum;
    omega_1 = omega_1/omega_sum;
    omega_2 = omega_2/omega_sum;
    omega_3 = omega_3/omega_sum;
    
    /*
     * Compute U_minus.
     */
    
    U_minus[idx_side] = Real(3)/Real(8)*omega_0*U_array[0][idx_side] +
        (-Real(10)/Real(8)*omega_0 - Real(1)/Real(8)*omega_1)*U_array[1][idx_side] +
        (Real(15)/Real(8)*omega_0 + Real(6)/Real(8)*omega_1 + Real(3)/Real(8)*omega_2)*
            U_array[2][idx_side] +
        (Real(3)/Real(8)*omega_1 + Real(6)/Real(8)*omega_2 + Real(15)/Real(8)*omega_3)*
            U_array[3][idx_side] +
        (-Real(1)/Real(8)*omega_2 - Real(10)/Real(8)*omega_3)*U_array[4][idx_side] +
        Real(3)/Real(8)*omega_3*U_array[5][idx_side];
}


/*
 * Perform local WENO interpolation of U_plus.
 */
static inline __attribute__((always_inline)) void performLocalWENOInterpolationPlus(
    Real* U_plus,
    Real** U_array,
    const int& idx_side,
    const int& q,
    const Real& C,
    const Real& Chi,
    const Real& epsilon,
    const Real& dx)
{
    /*
     * Compute beta_tilde's.
     */
    
    Real beta_tilde_0, beta_tilde_1, beta_tilde_2, beta_tilde_3;
    
    computeLocalBetaTilde(beta_tilde_0, beta_tilde_1, beta_tilde_2, beta_tilde_3, U_array, idx_side);
    
    /*
     * Compute the weights omega_tilde.
     */
    
    Real omega_tilde_0, omega_tilde_1, omega_tilde_2, omega_tilde_3;
    
    Real beta_avg_tilde = Real(1)/Real(8)*(beta_tilde_0 + beta_tilde_2 + 6*beta_tilde_1);
    Real tau_6_tilde = std::abs(beta_tilde_3 - beta_avg_tilde);
    
    omega_tilde_0 = Real(1)/Real(32)*ipow(C + tau_6_tilde/(beta_tilde_0 + epsilon*dx*dx)*
        (beta_avg_tilde + Chi*dx*dx)/(beta_tilde_0 + Chi*dx*dx), q);
    omega_tilde_1 = Real(15)/Real(32)*ipow(C + tau_6_tilde/(beta_tilde_1 + epsilon*dx*dx)*
        (beta_avg_tilde + Chi*dx*dx)/(beta_tilde_1 + Chi*dx*dx), q);
    omega_tilde_2 = Real(15)/Real(32)*ipow(C + tau_6_tilde/(beta_tilde_2 + epsilon*dx*dx)*
        (beta_avg_tilde + Chi*dx*dx)/(beta_tilde_2 + Chi*dx*dx), q);
    omega_tilde_3 = Real(1)/Real(32)*ipow(C + tau_6_tilde/(beta_tilde_3 + epsilon*dx*dx)*
        (beta_avg_tilde + Chi*dx*dx)/(beta_tilde_3 + Chi*dx*dx), q);
    
    Real omega_tilde_sum = omega_tilde_0 + omega_tilde_1 + omega_tilde_2 + omega_tilde_3;
    
    omega_tilde_0 = omega_tilde_0/omega_tilde_sum;
    omega_tilde_1 = omega_tilde_1/omega_tilde_sum;
    omega_tilde_2 = omega_tilde_2/omega_tilde_sum;
    omega_tilde_3 = omega_tilde_3/omega_tilde_sum;
    
    /*
     * Compute U_plus.
     */
    
    U_plus[idx_side] = Real(3)/Real(8)*omega_tilde_0*U_array[5][idx_side] +
        (-Real(10)/Real(8)*omega_tilde_0 - Real(1)/Real(8)*omega_tilde_1)*U_array[4][idx_side] +
        (Real(15)/Real(8)*omega_tilde_0 + Real(6)/Real(8)*omega_tilde_1 + Real(3)/Real(8)*omega_tilde_2)*
            U_array[3][idx_side] +
        (Real(3)/Real(8)*omega_tilde_1 + Real(6)/Real(8)*omega_tilde_2 + Real(15)/Real(8)*omega_tilde_3)*
            U_array[2][idx_side] +
        (-Real(1)/Real(8)*omega_tilde_2 - Real(10)/Real(8)*omega_tilde_3)*U_array[1][idx_side] +
        Real(3)/Real(8)*omega_tilde_3*U_array[0][idx_side];
}


ConvectiveFluxReconstructorWCNS6_CU_M2_HLLC_HLL::ConvectiveFluxReconstructorWCNS6_CU_M2_HLLC_HLL(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_eqn,
    const FLOW_MODEL::TYPE& flow_model_type,
    const HAMERS_SHARED_PTR<FlowModel>& flow_model,
    const HAMERS_SHARED_PTR<tbox::Database>& convective_flux_reconstructor_db):
        ConvectiveFluxReconstructorWCNS56(
            object_name,
            dim,
            grid_geometry,
            num_eqn,
            flow_model_type,
            flow_model,
            convective_flux_reconstructor_db)
{
    /*
     * Set the constants that are used in the scheme.
     */
    
    d_constant_q = convective_flux_reconstructor_db->
        getIntegerWithDefault("constant_q", 4);
    d_constant_q = convective_flux_reconstructor_db->
        getIntegerWithDefault("d_constant_q", d_constant_q);
    
    d_constant_C = convective_flux_reconstructor_db->
        getRealWithDefault("constant_C", Real(1000));
    d_constant_C = convective_flux_reconstructor_db->
        getRealWithDefault("d_constant_C", d_constant_C);
    
    d_constant_Chi = convective_flux_reconstructor_db->
        getRealWithDefault("constant_Chi", Real(1.0e8));
    d_constant_Chi = convective_flux_reconstructor_db->
        getRealWithDefault("d_constant_Chi", d_constant_Chi);
    
    d_constant_epsilon = convective_flux_reconstructor_db->
        getRealWithDefault("constant_epsilon", Real(1.0e-8));
    d_constant_epsilon = convective_flux_reconstructor_db->
        getRealWithDefault("d_constant_epsilon", d_constant_epsilon);
}


/*
 * Print all characteristics of the convective flux reconstruction class.
 */
void
ConvectiveFluxReconstructorWCNS6_CU_M2_HLLC_HLL::printClassData(
    std::ostream& os) const
{
    os << "\nPrint ConvectiveFluxReconstructorWCNS6_CU_M2_HLLC_HLL object..."
       << std::endl;
    
    os << std::endl;
    
    os << "ConvectiveFluxReconstructorWCNS6_CU_M2_HLLC_HLL: this = "
       << (ConvectiveFluxReconstructorWCNS6_CU_M2_HLLC_HLL *)this
       << std::endl;
    os << "d_object_name = "
       << d_object_name
       << std::endl;
    os << "d_constant_q = "
       << d_constant_q
       << std::endl;
    os << "d_constant_C = "
       << d_constant_C
       << std::endl;
    os << "d_constant_Chi = "
       << d_constant_Chi
       << std::endl;
    os << "d_constant_epsilon = "
       << d_constant_epsilon
       << std::endl;
}


/*
 * Put the characteristics of the convective flux reconstruction class
 * into the restart database.
 */
void
ConvectiveFluxReconstructorWCNS6_CU_M2_HLLC_HLL::putToRestart(
   const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const
{
    restart_db->putInteger("d_constant_q", d_constant_q);
    restart_db->putReal("d_constant_C", d_constant_C);
    restart_db->putReal("d_constant_Chi", d_constant_Chi);
    restart_db->putReal("d_constant_epsilon", d_constant_epsilon);
}


/*
 * Perform WENO interpolation.
 */
void
ConvectiveFluxReconstructorWCNS6_CU_M2_HLLC_HLL::performWENOInterpolation(
    std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& variables_minus,
    std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& variables_plus,
    const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > >& variables)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(variables_minus.size()) == d_num_eqn);
    TBOX_ASSERT(static_cast<int>(variables_plus.size()) == d_num_eqn);
    
    TBOX_ASSERT(static_cast<int>(variables.size()) == 6);
#endif
    
    /*
     * Get the interior dimensions.
     */
    
    const hier::IntVector interior_dims = variables_minus[0]->getBox().numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        TBOX_ASSERT(variables_minus[ei]->getBox().numberCells() == interior_dims);
        TBOX_ASSERT(variables_plus[ei]->getBox().numberCells() == interior_dims);
        
        TBOX_ASSERT(variables_minus[ei]->getGhostCellWidth() == hier::IntVector::getOne(d_dim));
        TBOX_ASSERT(variables_plus[ei]->getGhostCellWidth() == hier::IntVector::getOne(d_dim));
    }
    
    TBOX_ASSERT(static_cast<int>(variables.size()) == 6);
    
    for (int m = 0; m < 6; m++)
    {
        TBOX_ASSERT(static_cast<int>(variables[m].size()) == d_num_eqn);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ASSERT(variables[m][ei]->getBox().numberCells() == interior_dims);
            TBOX_ASSERT(variables[m][ei]->getGhostCellWidth() == hier::IntVector::getOne(d_dim));
        }
    }
#endif
    
    const double* const grid_spacing = d_grid_geometry->getDx();
    Real dx = 0.0;
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the dimension.
         */
        
        const int interior_dim_0 = interior_dims[0];
        
        /*
         * Peform WENO interpolation in the x-direction.
         */
        
        dx = Real(grid_spacing[0]);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<Real*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(0));
            }
            
            Real* U_L = variables_minus[ei]->getPointer(0);
            
            HAMERS_PRAGMA_SIMD
            for (int i = -1; i < interior_dim_0 + 2; i++)
            {
                // Compute the linear index of the mid-point.
                const int idx_midpoint_x = i + 1;
                
                performLocalWENOInterpolationMinus(
                    U_L,
                    U_array.data(),
                    idx_midpoint_x,
                    d_constant_q,
                    d_constant_C,
                    d_constant_Chi,
                    d_constant_epsilon,
                    dx);
            }
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<Real*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(0));
            }
            
            Real* U_R = variables_plus[ei]->getPointer(0);
            
            HAMERS_PRAGMA_SIMD
            for (int i = -1; i < interior_dim_0 + 2; i++)
            {
                // Compute the linear index of the mid-point.
                const int idx_midpoint_x = i + 1;
                
                performLocalWENOInterpolationPlus(
                    U_R,
                    U_array.data(),
                    idx_midpoint_x,
                    d_constant_q,
                    d_constant_C,
                    d_constant_Chi,
                    d_constant_epsilon,
                    dx);
            }
        }
        
    } // if (d_dim == tbox::Dimension(1))
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the dimensions.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        /*
         * Peform WENO interpolation in the x-direction.
         */
        
        dx = Real(grid_spacing[0]);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<Real*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(0));
            }
            
            Real* U_L = variables_minus[ei]->getPointer(0);
            
            for (int j = 0; j < interior_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = -1; i < interior_dim_0 + 2; i++)
                {
                    // Compute the linear index of the mid-point.
                    const int idx_midpoint_x = (i + 1) +
                        (j + 1)*(interior_dim_0 + 3);
                    
                    performLocalWENOInterpolationMinus(
                        U_L,
                        U_array.data(),
                        idx_midpoint_x,
                        d_constant_q,
                        d_constant_C,
                        d_constant_Chi,
                        d_constant_epsilon,
                        dx);
                }
            }
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<Real*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(0));
            }
            
            Real* U_R = variables_plus[ei]->getPointer(0);
            
            for (int j = 0; j < interior_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = -1; i < interior_dim_0 + 2; i++)
                {
                    // Compute the linear index of the mid-point.
                    const int idx_midpoint_x = (i + 1) +
                        (j + 1)*(interior_dim_0 + 3);
                    
                    performLocalWENOInterpolationPlus(
                        U_R,
                        U_array.data(),
                        idx_midpoint_x,
                        d_constant_q,
                        d_constant_C,
                        d_constant_Chi,
                        d_constant_epsilon,
                        dx);
                }
            }
        }
        
        /*
         * Peform WENO interpolation in the y-direction.
         */
        
        dx = Real(grid_spacing[1]);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<Real*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(1));
            }
            
            Real* U_B = variables_minus[ei]->getPointer(1);
            
            for (int j = -1; j < interior_dim_1 + 2; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index of the mid-point.
                    const int idx_midpoint_y = (i + 1) +
                        (j + 1)*(interior_dim_0 + 2);
                    
                    performLocalWENOInterpolationMinus(
                        U_B,
                        U_array.data(),
                        idx_midpoint_y,
                        d_constant_q,
                        d_constant_C,
                        d_constant_Chi,
                        d_constant_epsilon,
                        dx);
                }
            }
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<Real*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(1));
            }
            
            Real* U_T = variables_plus[ei]->getPointer(1);
            
            for (int j = -1; j < interior_dim_1 + 2; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index of the mid-point.
                    const int idx_midpoint_y = (i + 1) +
                        (j + 1)*(interior_dim_0 + 2);
                    
                    performLocalWENOInterpolationPlus(
                        U_T,
                        U_array.data(),
                        idx_midpoint_y,
                        d_constant_q,
                        d_constant_C,
                        d_constant_Chi,
                        d_constant_epsilon,
                        dx);
                }
            }
        }
        
    } // if (d_dim == tbox::Dimension(2))
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the dimensions.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        const int interior_dim_2 = interior_dims[2];
        
        /*
         * Peform WENO interpolation in the x-direction.
         */
        
        dx = Real(grid_spacing[0]);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<Real*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(0));
            }
            
            Real* U_L = variables_minus[ei]->getPointer(0);
            
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = -1; i < interior_dim_0 + 2; i++)
                    {
                        // Compute the linear index of the mid-point.
                        const int idx_midpoint_x = (i + 1) +
                            (j + 1)*(interior_dim_0 + 3) +
                            (k + 1)*(interior_dim_0 + 3)*
                                (interior_dim_1 + 2);
                        
                        performLocalWENOInterpolationMinus(
                            U_L,
                            U_array.data(),
                            idx_midpoint_x,
                            d_constant_q,
                            d_constant_C,
                            d_constant_Chi,
                            d_constant_epsilon,
                            dx);
                    }
                }
            }
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<Real*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(0));
            }
            
            Real* U_R = variables_plus[ei]->getPointer(0);
            
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = -1; i < interior_dim_0 + 2; i++)
                    {
                        // Compute the linear index of the mid-point.
                        const int idx_midpoint_x = (i + 1) +
                            (j + 1)*(interior_dim_0 + 3) +
                            (k + 1)*(interior_dim_0 + 3)*
                                (interior_dim_1 + 2);
                        
                        performLocalWENOInterpolationPlus(
                            U_R,
                            U_array.data(),
                            idx_midpoint_x,
                            d_constant_q,
                            d_constant_C,
                            d_constant_Chi,
                            d_constant_epsilon,
                            dx);
                    }
                }
            }
        }
        
        /*
         * Peform WENO interpolation in the y-direction.
         */
        
        dx = Real(grid_spacing[1]);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<Real*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(1));
            }
            
            Real* U_B = variables_minus[ei]->getPointer(1);
            
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = -1; j < interior_dim_1 + 2; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index of the mid-point.
                        const int idx_midpoint_y = (i + 1) +
                            (j + 1)*(interior_dim_0 + 2) +
                            (k + 1)*(interior_dim_0 + 2)*
                                (interior_dim_1 + 3);
                        
                        performLocalWENOInterpolationMinus(
                            U_B,
                            U_array.data(),
                            idx_midpoint_y,
                            d_constant_q,
                            d_constant_C,
                            d_constant_Chi,
                            d_constant_epsilon,
                            dx);
                    }
                }
            }
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<Real*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(1));
            }
            
            Real* U_T = variables_plus[ei]->getPointer(1);
            
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = -1; j < interior_dim_1 + 2; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index of the mid-point.
                        const int idx_midpoint_y = (i + 1) +
                            (j + 1)*(interior_dim_0 + 2) +
                            (k + 1)*(interior_dim_0 + 2)*
                                (interior_dim_1 + 3);
                        
                        performLocalWENOInterpolationPlus(
                            U_T,
                            U_array.data(),
                            idx_midpoint_y,
                            d_constant_q,
                            d_constant_C,
                            d_constant_Chi,
                            d_constant_epsilon,
                            dx);
                    }
                }
            }
        }
        
        /*
         * Peform WENO interpolation in the z-direction.
         */
        
        dx = Real(grid_spacing[2]);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<Real*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(2));
            }
            
            Real* U_B = variables_minus[ei]->getPointer(2);
            
            for (int k = -1; k < interior_dim_2 + 2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index of the mid-point.
                        const int idx_midpoint_z = (i + 1) +
                            (j + 1)*(interior_dim_0 + 2) +
                            (k + 1)*(interior_dim_0 + 2)*
                                (interior_dim_1 + 2);
                        
                        performLocalWENOInterpolationMinus(
                            U_B,
                            U_array.data(),
                            idx_midpoint_z,
                            d_constant_q,
                            d_constant_C,
                            d_constant_Chi,
                            d_constant_epsilon,
                            dx);
                    }
                }
            }
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<Real*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(2));
            }
            
            Real* U_F = variables_plus[ei]->getPointer(2);
            
            for (int k = -1; k < interior_dim_2 + 2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index of the mid-point.
                        const int idx_midpoint_z = (i + 1) +
                            (j + 1)*(interior_dim_0 + 2) +
                            (k + 1)*(interior_dim_0 + 2)*
                                (interior_dim_1 + 2);
                        
                        performLocalWENOInterpolationPlus(
                            U_F,
                            U_array.data(),
                            idx_midpoint_z,
                            d_constant_q,
                            d_constant_C,
                            d_constant_Chi,
                            d_constant_epsilon,
                            dx);
                    }
                }
            }
        }
        
    } // if (d_dim == tbox::Dimension(3))
}
