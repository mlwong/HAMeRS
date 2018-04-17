#include "flow/convective_flux_reconstructors/WCNS56/ConvectiveFluxReconstructorWCNS6-CU-M2-HLLC-HLL.hpp"

#include <cfloat>

#include "SAMRAI/geom/CartesianPatchGeometry.h"


/*
 * Interger based power function.
 */
static inline __attribute__((always_inline)) double ipow(double base, int exp)
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
    double* beta_0,
    double* beta_1,
    double* beta_2,
    double* beta_3,
    double** U_array,
    int idx_side)
{
    *beta_0 = double(1)/double(3)*(U_array[0][idx_side]*(double(4)*U_array[0][idx_side] -
         double(19)*U_array[1][idx_side] + double(11)*U_array[2][idx_side]) +
         U_array[1][idx_side]*(double(25)*U_array[1][idx_side] - double(31)*U_array[2][idx_side]) +
         double(10)*U_array[2][idx_side]*U_array[2][idx_side]);
    
    *beta_1 = double(1)/double(3)*(U_array[1][idx_side]*(double(4)*U_array[1][idx_side] -
         double(13)*U_array[2][idx_side] + double(5)*U_array[3][idx_side]) +
         double(13)*U_array[2][idx_side]*(U_array[2][idx_side] - U_array[3][idx_side]) +
         double(4)*U_array[3][idx_side]*U_array[3][idx_side]);
    
    *beta_2 = double(1)/double(3)*(U_array[2][idx_side]*(double(10)*U_array[2][idx_side] -
         double(31)*U_array[3][idx_side] + double(11)*U_array[4][idx_side]) +
         U_array[3][idx_side]*(double(25)*U_array[3][idx_side] - double(19)*U_array[4][idx_side]) +
         double(4)*U_array[4][idx_side]*U_array[4][idx_side]);
    
    *beta_3 = double(1)/double(232243200)*(U_array[0][idx_side]*(double(525910327)*U_array[0][idx_side] -
         double(4562164630)*U_array[1][idx_side] + double(7799501420)*U_array[2][idx_side] -
         double(6610694540)*U_array[3][idx_side] + double(2794296070)*U_array[4][idx_side] -
         double(472758974)*U_array[5][idx_side]) + double(5)*U_array[1][idx_side]*
        (double(2146987907)*U_array[1][idx_side] - double(7722406988)*U_array[2][idx_side] +
         double(6763559276)*U_array[3][idx_side] - double(2926461814)*U_array[4][idx_side] +
         double(503766638)*U_array[5][idx_side]) + double(20)*U_array[2][idx_side]*
        (double(1833221603)*U_array[2][idx_side] - double(3358664662)*U_array[3][idx_side] +
         double(1495974539)*U_array[4][idx_side] - double(263126407)*U_array[5][idx_side]) +
        double(20)*U_array[3][idx_side]*(double(1607794163)*U_array[3][idx_side] -
         double(1486026707)*U_array[4][idx_side] + double(268747951)*U_array[5][idx_side]) +
        double(5)*U_array[4][idx_side]*(double(1432381427)*U_array[4][idx_side] -
         double(536951582)*U_array[5][idx_side]) +
        double(263126407)*U_array[5][idx_side]*U_array[5][idx_side]);
}


/*
 * Compute local beta_tilde's.
 */
static inline __attribute__((always_inline)) void computeLocalBetaTilde(
    double* beta_tilde_0,
    double* beta_tilde_1,
    double* beta_tilde_2,
    double* beta_tilde_3,
    double** U_array,
    int idx_side)
{
    *beta_tilde_0 = double(1)/double(3)*(U_array[5][idx_side]*(double(4)*U_array[5][idx_side] -
         double(19)*U_array[4][idx_side] + double(11)*U_array[3][idx_side]) +
         U_array[4][idx_side]*(double(25)*U_array[4][idx_side] - double(31)*U_array[3][idx_side]) +
         double(10)*U_array[3][idx_side]*U_array[3][idx_side]);
    
    *beta_tilde_1 = double(1)/double(3)*(U_array[4][idx_side]*(double(4)*U_array[4][idx_side] -
         double(13)*U_array[3][idx_side] + double(5)*U_array[2][idx_side]) +
         double(13)*U_array[3][idx_side]*(U_array[3][idx_side] - U_array[2][idx_side]) +
         double(4)*U_array[2][idx_side]*U_array[2][idx_side]);
    
    *beta_tilde_2 = double(1)/double(3)*(U_array[3][idx_side]*(double(10)*U_array[3][idx_side] -
         double(31)*U_array[2][idx_side] + double(11)*U_array[1][idx_side]) +
         U_array[2][idx_side]*(double(25)*U_array[2][idx_side] - double(19)*U_array[1][idx_side]) +
         double(4)*U_array[1][idx_side]*U_array[1][idx_side]);
    
    *beta_tilde_3 = double(1)/double(232243200)*(U_array[5][idx_side]*(double(525910327)*U_array[5][idx_side] -
         double(4562164630)*U_array[4][idx_side] + double(7799501420)*U_array[3][idx_side] -
         double(6610694540)*U_array[2][idx_side] + double(2794296070)*U_array[1][idx_side] -
         double(472758974)*U_array[0][idx_side]) + double(5)*U_array[4][idx_side]*
        (double(2146987907)*U_array[4][idx_side] - double(7722406988)*U_array[3][idx_side] +
         double(6763559276)*U_array[2][idx_side] - double(2926461814)*U_array[1][idx_side] +
         double(503766638)*U_array[0][idx_side]) + double(20)*U_array[3][idx_side]*
        (double(1833221603)*U_array[3][idx_side] - double(3358664662)*U_array[2][idx_side] +
         double(1495974539)*U_array[1][idx_side] - double(263126407)*U_array[0][idx_side]) +
        double(20)*U_array[2][idx_side]*(double(1607794163)*U_array[2][idx_side] -
         double(1486026707)*U_array[1][idx_side] + double(268747951)*U_array[0][idx_side]) +
        double(5)*U_array[1][idx_side]*(double(1432381427)*U_array[1][idx_side] -
         double(536951582)*U_array[0][idx_side]) +
        double(263126407)*U_array[0][idx_side]*U_array[0][idx_side]);
}


/*
 * Perform local WENO interpolation of U_minus.
 */
static inline __attribute__((always_inline)) void performLocalWENOInterpolationMinus(
   double* U_minus,
   double** U_array,
   int idx_side,
   int q,
   double C,
   double Chi,
   double epsilon,
   double dx)
{
    /*
     * Compute beta's.
     */
    
    double beta_0, beta_1, beta_2, beta_3;
    
    computeLocalBeta(&beta_0, &beta_1, &beta_2, &beta_3, U_array, idx_side);
    
    /*
     * Compute the weights omega.
     */
    
    double omega_0, omega_1, omega_2, omega_3;
    
    double beta_avg = double(1)/double(8)*(beta_0 + beta_2 + 6*beta_1);
    double tau_6 = fabs(beta_3 - beta_avg);
    
    omega_0 = double(1)/double(32)*ipow(C + tau_6/(beta_0 + epsilon*dx*dx)*
        (beta_avg + Chi*dx*dx)/(beta_0 + Chi*dx*dx), q);
    omega_1 = double(15)/double(32)*ipow(C + tau_6/(beta_1 + epsilon*dx*dx)*
        (beta_avg + Chi*dx*dx)/(beta_1 + Chi*dx*dx), q);
    omega_2 = double(15)/double(32)*ipow(C + tau_6/(beta_2 + epsilon*dx*dx)*
        (beta_avg + Chi*dx*dx)/(beta_2 + Chi*dx*dx), q);
    omega_3 = double(1)/double(32)*ipow(C + tau_6/(beta_3 + epsilon*dx*dx)*
        (beta_avg + Chi*dx*dx)/(beta_3 + Chi*dx*dx), q);
    
    double omega_sum = omega_0 + omega_1 + omega_2 + omega_3;
    
    omega_0 = omega_0/omega_sum;
    omega_1 = omega_1/omega_sum;
    omega_2 = omega_2/omega_sum;
    omega_3 = omega_3/omega_sum;
    
    /*
     * Compute U_minus.
     */
    
    U_minus[idx_side] = double(3)/double(8)*omega_0*U_array[0][idx_side] +
        (-double(10)/double(8)*omega_0 - double(1)/double(8)*omega_1)*U_array[1][idx_side] +
        (double(15)/double(8)*omega_0 + double(6)/double(8)*omega_1 + double(3)/double(8)*omega_2)*
            U_array[2][idx_side] +
        (double(3)/double(8)*omega_1 + double(6)/double(8)*omega_2 + double(15)/double(8)*omega_3)*
            U_array[3][idx_side] +
        (-double(1)/double(8)*omega_2 - double(10)/double(8)*omega_3)*U_array[4][idx_side] +
        double(3)/double(8)*omega_3*U_array[5][idx_side];
}


/*
 * Perform local WENO interpolation of U_plus.
 */
static inline __attribute__((always_inline)) void performLocalWENOInterpolationPlus(
   double* U_plus,
   double** U_array,
   int idx_side,
   int q,
   double C,
   double Chi,
   double epsilon,
   double dx)
{
    /*
     * Compute beta_tilde's.
     */
    
    double beta_tilde_0, beta_tilde_1, beta_tilde_2, beta_tilde_3;
    
    computeLocalBetaTilde(&beta_tilde_0, &beta_tilde_1, &beta_tilde_2, &beta_tilde_3, U_array, idx_side);
    
    /*
     * Compute the weights omega_tilde.
     */
    
    double omega_tilde_0, omega_tilde_1, omega_tilde_2, omega_tilde_3;
    
    double beta_avg_tilde = double(1)/double(8)*(beta_tilde_0 + beta_tilde_2 + 6*beta_tilde_1);
    double tau_6_tilde = fabs(beta_tilde_3 - beta_avg_tilde);
    
    omega_tilde_0 = double(1)/double(32)*ipow(C + tau_6_tilde/(beta_tilde_0 + epsilon*dx*dx)*
        (beta_avg_tilde + Chi*dx*dx)/(beta_tilde_0 + Chi*dx*dx), q);
    omega_tilde_1 = double(15)/double(32)*ipow(C + tau_6_tilde/(beta_tilde_1 + epsilon*dx*dx)*
        (beta_avg_tilde + Chi*dx*dx)/(beta_tilde_1 + Chi*dx*dx), q);
    omega_tilde_2 = double(15)/double(32)*ipow(C + tau_6_tilde/(beta_tilde_2 + epsilon*dx*dx)*
        (beta_avg_tilde + Chi*dx*dx)/(beta_tilde_2 + Chi*dx*dx), q);
    omega_tilde_3 = double(1)/double(32)*ipow(C + tau_6_tilde/(beta_tilde_3 + epsilon*dx*dx)*
        (beta_avg_tilde + Chi*dx*dx)/(beta_tilde_3 + Chi*dx*dx), q);
    
    double omega_tilde_sum = omega_tilde_0 + omega_tilde_1 + omega_tilde_2 + omega_tilde_3;
    
    omega_tilde_0 = omega_tilde_0/omega_tilde_sum;
    omega_tilde_1 = omega_tilde_1/omega_tilde_sum;
    omega_tilde_2 = omega_tilde_2/omega_tilde_sum;
    omega_tilde_3 = omega_tilde_3/omega_tilde_sum;
    
    /*
     * Compute U_plus.
     */
    
    U_plus[idx_side] = double(3)/double(8)*omega_tilde_0*U_array[5][idx_side] +
        (-double(10)/double(8)*omega_tilde_0 - double(1)/double(8)*omega_tilde_1)*U_array[4][idx_side] +
        (double(15)/double(8)*omega_tilde_0 + double(6)/double(8)*omega_tilde_1 + double(3)/double(8)*omega_tilde_2)*
            U_array[3][idx_side] +
        (double(3)/double(8)*omega_tilde_1 + double(6)/double(8)*omega_tilde_2 + double(15)/double(8)*omega_tilde_3)*
            U_array[2][idx_side] +
        (-double(1)/double(8)*omega_tilde_2 - double(10)/double(8)*omega_tilde_3)*U_array[1][idx_side] +
        double(3)/double(8)*omega_tilde_3*U_array[0][idx_side];
}


ConvectiveFluxReconstructorWCNS6_CU_M2_HLLC_HLL::ConvectiveFluxReconstructorWCNS6_CU_M2_HLLC_HLL(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_eqn,
    const int& num_species,
    const boost::shared_ptr<FlowModel>& flow_model,
    const boost::shared_ptr<tbox::Database>& convective_flux_reconstructor_db):
        ConvectiveFluxReconstructorWCNS56(
            object_name,
            dim,
            grid_geometry,
            num_eqn,
            num_species,
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
        getDoubleWithDefault("constant_C", double(1000));
    d_constant_C = convective_flux_reconstructor_db->
        getDoubleWithDefault("d_constant_C", d_constant_C);
    
    d_constant_Chi = convective_flux_reconstructor_db->
        getDoubleWithDefault("constant_Chi", double(1.0e8));
    d_constant_Chi = convective_flux_reconstructor_db->
        getDoubleWithDefault("d_constant_Chi", d_constant_Chi);
    
    d_constant_epsilon = convective_flux_reconstructor_db->
        getDoubleWithDefault("constant_epsilon", double(1.0e-8));
    d_constant_epsilon = convective_flux_reconstructor_db->
        getDoubleWithDefault("d_constant_epsilon", d_constant_epsilon);
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
   const boost::shared_ptr<tbox::Database>& restart_db) const
{
    restart_db->putInteger("d_constant_q", d_constant_q);
    restart_db->putDouble("d_constant_C", d_constant_C);
    restart_db->putDouble("d_constant_Chi", d_constant_Chi);
    restart_db->putDouble("d_constant_epsilon", d_constant_epsilon);
}


/*
 * Perform WENO interpolation.
 */
void
ConvectiveFluxReconstructorWCNS6_CU_M2_HLLC_HLL::performWENOInterpolation(
    std::vector<boost::shared_ptr<pdat::SideData<double> > >& variables_minus,
    std::vector<boost::shared_ptr<pdat::SideData<double> > >& variables_plus,
    const std::vector<std::vector<boost::shared_ptr<pdat::SideData<double> > > >& variables)
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
    double dx = 0.0;
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the dimension.
         */
        
        const int interior_dim_0 = interior_dims[0];
        
        /*
         * Peform WENO interpolation in the x-direction.
         */
        
        dx = grid_spacing[0];
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<double*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(0));
            }
            
            double* U_L = variables_minus[ei]->getPointer(0);
            
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
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
            std::vector<double*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(0));
            }
            
            double* U_R = variables_plus[ei]->getPointer(0);
            
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
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
        
        dx = grid_spacing[0];
        
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
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
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
            std::vector<double*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(0));
            }
            
            double* U_R = variables_plus[ei]->getPointer(0);
            
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
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
        
        dx = grid_spacing[1];
        
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
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
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
            std::vector<double*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(1));
            }
            
            double* U_T = variables_plus[ei]->getPointer(1);
            
            for (int j = -1; j < interior_dim_1 + 2; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
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
        
        dx = grid_spacing[0];
        
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
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
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
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
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
        
        dx = grid_spacing[1];
        
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
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
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
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
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
        
        dx = grid_spacing[2];
        
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
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
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
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
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
