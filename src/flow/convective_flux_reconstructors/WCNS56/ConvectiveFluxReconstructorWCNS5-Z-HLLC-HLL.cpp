#include "flow/convective_flux_reconstructors/WCNS56/ConvectiveFluxReconstructorWCNS5-Z-HLLC-HLL.hpp"

#define EPSILON HAMERS_EPSILON


/*
 * Interger based power function.
 */
static inline __attribute__((always_inline)) double ipow(double base, const int& exp)
{
    double result = base;
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
    double& beta_0,
    double& beta_1,
    double& beta_2,
    double** U_array,
    const int& idx_side)
{
    beta_0 = double(1)/double(3)*(U_array[0][idx_side]*(double(4)*U_array[0][idx_side] -
         double(19)*U_array[1][idx_side] + double(11)*U_array[2][idx_side]) +
         U_array[1][idx_side]*(double(25)*U_array[1][idx_side] - double(31)*U_array[2][idx_side]) +
         double(10)*U_array[2][idx_side]*U_array[2][idx_side]);
    
    beta_1 = double(1)/double(3)*(U_array[1][idx_side]*(double(4)*U_array[1][idx_side] -
         double(13)*U_array[2][idx_side] + double(5)*U_array[3][idx_side]) +
         double(13)*U_array[2][idx_side]*(U_array[2][idx_side] - U_array[3][idx_side]) +
         double(4)*U_array[3][idx_side]*U_array[3][idx_side]);
    
    beta_2 = double(1)/double(3)*(U_array[2][idx_side]*(double(10)*U_array[2][idx_side] -
         double(31)*U_array[3][idx_side] + double(11)*U_array[4][idx_side]) +
         U_array[3][idx_side]*(double(25)*U_array[3][idx_side] - double(19)*U_array[4][idx_side]) +
         double(4)*U_array[4][idx_side]*U_array[4][idx_side]);
}


/*
 * Compute local beta_tilde's.
 */
static inline __attribute__((always_inline)) void computeLocalBetaTilde(
    double& beta_tilde_0,
    double& beta_tilde_1,
    double& beta_tilde_2,
    double** U_array,
    const int& idx_side)
{
    beta_tilde_0 = double(1)/double(3)*(U_array[5][idx_side]*(double(4)*U_array[5][idx_side] -
         double(19)*U_array[4][idx_side] + double(11)*U_array[3][idx_side]) +
         U_array[4][idx_side]*(double(25)*U_array[4][idx_side] - double(31)*U_array[3][idx_side]) +
         double(10)*U_array[3][idx_side]*U_array[3][idx_side]);
    
    beta_tilde_1 = double(1)/double(3)*(U_array[4][idx_side]*(double(4)*U_array[4][idx_side] -
         double(13)*U_array[3][idx_side] + double(5)*U_array[2][idx_side]) +
         double(13)*U_array[3][idx_side]*(U_array[3][idx_side] - U_array[2][idx_side]) +
         double(4)*U_array[2][idx_side]*U_array[2][idx_side]);
    
    beta_tilde_2 = double(1)/double(3)*(U_array[3][idx_side]*(double(10)*U_array[3][idx_side] -
         double(31)*U_array[2][idx_side] + double(11)*U_array[1][idx_side]) +
         U_array[2][idx_side]*(double(25)*U_array[2][idx_side] - double(19)*U_array[1][idx_side]) +
         double(4)*U_array[1][idx_side]*U_array[1][idx_side]);
}


/*
 * Perform local WENO interpolation of U_minus.
 */
static inline __attribute__((always_inline)) void performLocalWENOInterpolationMinus(
    double* U_minus,
    double** U_array,
    const int& idx_side,
    const int& p)
{
    /*
     * Compute beta's.
     */
    
    double beta_0, beta_1, beta_2;
    
    computeLocalBeta(beta_0, beta_1, beta_2, U_array, idx_side);
    
    /*
     * Compute the weights omega.
     */
    
    double omega_0, omega_1, omega_2;
    
    double tau_5 = fabs(beta_0 - beta_2);
    
    omega_0 = double(1)/double(16)*(double(1) + ipow(tau_5/(beta_0 + EPSILON), p));
    omega_1 = double(5)/double(8)*(double(1) + ipow(tau_5/(beta_1 + EPSILON), p));
    omega_2 = double(5)/double(16)*(double(1) + ipow(tau_5/(beta_2 + EPSILON), p));
    
    double omega_sum = omega_0 + omega_1 + omega_2;
    
    omega_0 = omega_0/omega_sum;
    omega_1 = omega_1/omega_sum;
    omega_2 = omega_2/omega_sum;
    
    /*
     * Compute U_minus.
     */
    
    U_minus[idx_side] = double(3)/double(8)*omega_0*U_array[0][idx_side] +
        (-double(10)/double(8)*omega_0 - double(1)/double(8)*omega_1)*U_array[1][idx_side] +
        (double(15)/double(8)*omega_0 + double(6)/double(8)*omega_1 +
        double(3)/double(8)*omega_2)*U_array[2][idx_side] +
        (double(3)/double(8)*omega_1 + double(6)/double(8)*omega_2)*U_array[3][idx_side] -
        double(1)/double(8)*omega_2*U_array[4][idx_side];
}


/*
 * Perform local WENO interpolation of U_plus.
 */
static inline __attribute__((always_inline)) void performLocalWENOInterpolationPlus(
    double* U_plus,
    double** U_array,
    const int& idx_side,
    const int& p)
{
    /*
     * Compute beta_tilde's.
     */
    
    double beta_tilde_0, beta_tilde_1, beta_tilde_2;
    
    computeLocalBetaTilde(beta_tilde_0, beta_tilde_1, beta_tilde_2, U_array, idx_side);
    
    /*
     * Compute the weights omega_upwind_tilde.
     */
    
    double omega_tilde_0, omega_tilde_1, omega_tilde_2;
    
    double tau_5_tilde = fabs(beta_tilde_0 - beta_tilde_2);
    
    omega_tilde_0 = double(1)/double(16)*(double(1) + ipow(tau_5_tilde/(beta_tilde_0 + EPSILON), p));
    omega_tilde_1 = double(5)/double(8)*(double(1) + ipow(tau_5_tilde/(beta_tilde_1 + EPSILON), p));
    omega_tilde_2 = double(5)/double(16)*(double(1) + ipow(tau_5_tilde/(beta_tilde_2 + EPSILON), p));
    
    double omega_tilde_sum = omega_tilde_0 + omega_tilde_1 + omega_tilde_2;
    
    omega_tilde_0 = omega_tilde_0/omega_tilde_sum;
    omega_tilde_1 = omega_tilde_1/omega_tilde_sum;
    omega_tilde_2 = omega_tilde_2/omega_tilde_sum;
    
    /*
     * Compute U_plus.
     */
    
    U_plus[idx_side] = double(3)/double(8)*omega_tilde_0*U_array[5][idx_side] +
        (-double(10)/double(8)*omega_tilde_0 - double(1)/double(8)*omega_tilde_1)*U_array[4][idx_side] +
        (double(15)/double(8)*omega_tilde_0 + double(6)/double(8)*omega_tilde_1 +
        double(3)/double(8)*omega_tilde_2)*U_array[3][idx_side] +
        (double(3)/double(8)*omega_tilde_1 + double(6)/double(8)*omega_tilde_2)*U_array[2][idx_side] -
        double(1)/double(8)*omega_tilde_2*U_array[1][idx_side];
}


ConvectiveFluxReconstructorWCNS5_Z_HLLC_HLL::ConvectiveFluxReconstructorWCNS5_Z_HLLC_HLL(
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
     * Set the constant that is used in the scheme.
     */
    
    d_constant_p = d_convective_flux_reconstructor_db->
        getIntegerWithDefault("constant_p", 2);
    d_constant_p = d_convective_flux_reconstructor_db->
        getIntegerWithDefault("d_constant_p", d_constant_p);
}


/*
 * Print all characteristics of the convective flux reconstruction class.
 */
void
ConvectiveFluxReconstructorWCNS5_Z_HLLC_HLL::printClassData(
    std::ostream& os) const
{
    os << "\nPrint ConvectiveFluxReconstructorWCNS5_Z_HLLC_HLL object..."
       << std::endl;
    
    os << std::endl;
    
    os << "ConvectiveFluxReconstructorWCNS5_Z_HLLC_HLL: this = "
       << (ConvectiveFluxReconstructorWCNS5_Z_HLLC_HLL *)this
       << std::endl;
    os << "d_object_name = "
       << d_object_name
       << std::endl;
    os << "d_constant_p = "
       << d_constant_p
       << std::endl;
}


/*
 * Put the characteristics of the convective flux reconstruction class
 * into the restart database.
 */
void
ConvectiveFluxReconstructorWCNS5_Z_HLLC_HLL::putToRestart(
   const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const
{
    restart_db->putInteger("d_constant_p", d_constant_p);
}


/*
 * Perform WENO interpolation.
 */
void
ConvectiveFluxReconstructorWCNS5_Z_HLLC_HLL::performWENOInterpolation(
    std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& variables_minus,
    std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& variables_plus,
    const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > >& variables)
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
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the dimension.
         */
        
        const int interior_dim_0 = interior_dims[0];
        
        /*
         * Peform WENO interpolation in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<double*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(0));
            }
            
            double* U_L = variables_minus[ei]->getPointer(0);
            
            HAMERS_PRAGMA_SIMD
            for (int i = -1; i < interior_dim_0 + 2; i++)
            {
                // Compute the linear index of the mid-point.
                const int idx_midpoint_x = i + 1;
                
                performLocalWENOInterpolationMinus(
                    U_L,
                    U_array.data(),
                    idx_midpoint_x,
                    d_constant_p);
            }
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<double*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(0));
            }
            
            double* U_R = variables_plus[ei]->getPointer(0);
            
            HAMERS_PRAGMA_SIMD
            for (int i = -1; i < interior_dim_0 + 2; i++)
            {
                // Compute the linear index of the mid-point.
                const int idx_midpoint_x = i + 1;
                
                performLocalWENOInterpolationPlus(
                    U_R,
                    U_array.data(),
                    idx_midpoint_x,
                    d_constant_p);
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
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<double*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(0));
            }
            
            double* U_L = variables_minus[ei]->getPointer(0);
            
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
                        d_constant_p);
                }
            }
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<double*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(0));
            }
            
            double* U_R = variables_plus[ei]->getPointer(0);
            
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
                        d_constant_p);
                }
            }
        }
        
        /*
         * Peform WENO interpolation in the y-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<double*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(1));
            }
            
            double* U_B = variables_minus[ei]->getPointer(1);
            
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
                        d_constant_p);
                }
            }
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<double*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(1));
            }
            
            double* U_T = variables_plus[ei]->getPointer(1);
            
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
                        d_constant_p);
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
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<double*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(0));
            }
            
            double* U_L = variables_minus[ei]->getPointer(0);
            
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
                            d_constant_p);
                    }
                }
            }
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<double*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(0));
            }
            
            double* U_R = variables_plus[ei]->getPointer(0);
            
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
                            d_constant_p);
                    }
                }
            }
        }
        
        /*
         * Peform WENO interpolation in the y-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<double*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(1));
            }
            
            double* U_B = variables_minus[ei]->getPointer(1);
            
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
                            d_constant_p);
                    }
                }
            }
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<double*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(1));
            }
            
            double* U_T = variables_plus[ei]->getPointer(1);
            
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
                            d_constant_p);
                    }
                }
            }
        }
        
        /*
         * Peform WENO interpolation in the z-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<double*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(2));
            }
            
            double* U_B = variables_minus[ei]->getPointer(2);
            
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
                            d_constant_p);
                    }
                }
            }
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<double*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(2));
            }
            
            double* U_F = variables_plus[ei]->getPointer(2);
            
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
                            d_constant_p);
                    }
                }
            }
        }
        
    } // if (d_dim == tbox::Dimension(3))
}
