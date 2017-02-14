#include "flow/convective_flux_reconstructors/WCNS56/ConvectiveFluxReconstructorWCNS6-LD-HLLC-HLL.hpp"

#define EPSILON 1e-40

/*
 * Compute local sigma.
 */
static inline __attribute__((always_inline)) void computeLocalSigma(
    double* sigma,
    double** U_array,
    int idx_side)
{
    /*
     * Compute the sigma.
     */
    
    const double alpha_1 = U_array[2][idx_side] - U_array[1][idx_side];
    const double alpha_2 = U_array[3][idx_side] - U_array[2][idx_side];
    const double alpha_3 = U_array[4][idx_side] - U_array[3][idx_side];
    
    const double theta_1 = fabs(alpha_1 - alpha_2)/(fabs(alpha_1) + fabs(alpha_2) + EPSILON);
    const double theta_2 = fabs(alpha_2 - alpha_3)/(fabs(alpha_2) + fabs(alpha_3) + EPSILON);
    
    *sigma = fmax(theta_1, theta_2);
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
    *beta_0 = 1.0/3.0*(U_array[0][idx_side]*(4.0*U_array[0][idx_side] - 19.0*U_array[1][idx_side] +
         11.0*U_array[2][idx_side]) + U_array[1][idx_side]*(25.0*U_array[1][idx_side] -
         31.0*U_array[2][idx_side]) + 10.0*U_array[2][idx_side]*U_array[2][idx_side]);
    
    *beta_1 = 1.0/3.0*(U_array[1][idx_side]*(4.0*U_array[1][idx_side] - 13.0*U_array[2][idx_side] +
         5.0*U_array[3][idx_side]) + 13.0*U_array[2][idx_side]*(U_array[2][idx_side] -
         U_array[3][idx_side]) + 4.0*U_array[3][idx_side]*U_array[3][idx_side]);
    
    *beta_2 = 1.0/3.0*(U_array[2][idx_side]*(10.0*U_array[2][idx_side] - 31.0*U_array[3][idx_side] +
         11.0*U_array[4][idx_side]) + U_array[3][idx_side]*(25.0*U_array[3][idx_side] -
         19.0*U_array[4][idx_side]) + 4.0*U_array[4][idx_side]*U_array[4][idx_side]);
    
    *beta_3 = 1.0/232243200.0*(U_array[0][idx_side]*(525910327.0*U_array[0][idx_side] -
         4562164630.0*U_array[1][idx_side] + 7799501420.0*U_array[2][idx_side] -
         6610694540.0*U_array[3][idx_side] + 2794296070.0*U_array[4][idx_side] -
         472758974.0*U_array[5][idx_side]) + 5.0*U_array[1][idx_side]*
        (2146987907.0*U_array[1][idx_side] - 7722406988.0*U_array[2][idx_side] +
         6763559276.0*U_array[3][idx_side] - 2926461814.0*U_array[4][idx_side] +
         503766638.0*U_array[5][idx_side]) + 20.0*U_array[2][idx_side]*
        (1833221603.0*U_array[2][idx_side] - 3358664662.0*U_array[3][idx_side] +
         1495974539.0*U_array[4][idx_side] - 263126407.0*U_array[5][idx_side]) +
        20.0*U_array[3][idx_side]*(1607794163.0*U_array[3][idx_side] -
         1486026707.0*U_array[4][idx_side] + 268747951.0*U_array[5][idx_side]) +
        5.0*U_array[4][idx_side]*(1432381427.0*U_array[4][idx_side] -
         536951582.0*U_array[5][idx_side]) +
        263126407.0*U_array[5][idx_side]*U_array[5][idx_side]);
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
    *beta_tilde_0 = 1.0/3.0*(U_array[5][idx_side]*(4.0*U_array[5][idx_side] - 19.0*U_array[4][idx_side] +
         11.0*U_array[3][idx_side]) + U_array[4][idx_side]*(25.0*U_array[4][idx_side] -
         31.0*U_array[3][idx_side]) + 10.0*U_array[3][idx_side]*U_array[3][idx_side]);
    
    *beta_tilde_1 = 1.0/3.0*(U_array[4][idx_side]*(4.0*U_array[4][idx_side] - 13.0*U_array[3][idx_side] +
         5.0*U_array[2][idx_side]) + 13.0*U_array[3][idx_side]*(U_array[3][idx_side] -
         U_array[2][idx_side]) + 4.0*U_array[2][idx_side]*U_array[2][idx_side]);
    
    *beta_tilde_2 = 1.0/3.0*(U_array[3][idx_side]*(10.0*U_array[3][idx_side] - 31.0*U_array[2][idx_side] +
         11.0*U_array[1][idx_side]) + U_array[2][idx_side]*(25.0*U_array[2][idx_side] -
         19.0*U_array[1][idx_side]) + 4.0*U_array[1][idx_side]*U_array[1][idx_side]);
    
    *beta_tilde_3 = 1.0/232243200.0*(U_array[5][idx_side]*(525910327.0*U_array[5][idx_side] -
         4562164630.0*U_array[4][idx_side] + 7799501420.0*U_array[3][idx_side] -
         6610694540.0*U_array[2][idx_side] + 2794296070.0*U_array[1][idx_side] -
         472758974.0*U_array[0][idx_side]) + 5.0*U_array[4][idx_side]*
        (2146987907.0*U_array[4][idx_side] - 7722406988.0*U_array[3][idx_side] +
         6763559276.0*U_array[2][idx_side] - 2926461814.0*U_array[1][idx_side] +
         503766638.0*U_array[0][idx_side]) + 20.0*U_array[3][idx_side]*
        (1833221603.0*U_array[3][idx_side] - 3358664662.0*U_array[2][idx_side] +
         1495974539.0*U_array[1][idx_side] - 263126407.0*U_array[0][idx_side]) +
        20.0*U_array[2][idx_side]*(1607794163.0*U_array[2][idx_side] -
         1486026707.0*U_array[1][idx_side] + 268747951.0*U_array[0][idx_side]) +
        5.0*U_array[1][idx_side]*(1432381427.0*U_array[1][idx_side] -
         536951582.0*U_array[0][idx_side]) +
        263126407.0*U_array[0][idx_side]*U_array[0][idx_side]);
}


/*
 * Perform local WENO interpolation.
 */
static inline __attribute__((always_inline)) void performLocalWENOInterpolation(
   double* U_minus,
   double* U_plus,
   double** U_array,
   int idx_side,
   int p,
   int q,
   double C,
   double alpha_tau)
{
    /*
     * Compute sigma.
     */
    
    double sigma;
    
    computeLocalSigma(&sigma, U_array, idx_side);
    
    /*
     * Compute beta's.
     */
    
    double beta_0, beta_1, beta_2, beta_3;
    
    computeLocalBeta(&beta_0, &beta_1, &beta_2, &beta_3, U_array, idx_side);
    
    /*
     * Compute the weights omega_upwind.
     */
    
    double omega_upwind_0, omega_upwind_1, omega_upwind_2;
    
    double tau_5 = fabs(beta_0 - beta_2);
    
    omega_upwind_0 = 1.0/16.0*(1.0 + pow(tau_5/(beta_0 + EPSILON), p));
    omega_upwind_1 = 5.0/8.0*(1.0 + pow(tau_5/(beta_1 + EPSILON), p));
    omega_upwind_2 = 5.0/16.0*(1.0 + pow(tau_5/(beta_2 + EPSILON), p));
    
    double omega_upwind_sum = omega_upwind_0 + omega_upwind_1 + omega_upwind_2;
    
    omega_upwind_0 = omega_upwind_0/omega_upwind_sum;
    omega_upwind_1 = omega_upwind_1/omega_upwind_sum;
    omega_upwind_2 = omega_upwind_2/omega_upwind_sum;
    
    /*
     * Compute the weights omega_central (store in omega first).
     */
    
    double omega_0, omega_1, omega_2, omega_3;
    
    double beta_avg = 1.0/8*(beta_0 + beta_2 + 6*beta_1);
    double tau_6 = fabs(beta_3 - beta_avg);
    
    omega_0 = 1.0/32.0*(C + pow(tau_6/(beta_0 + EPSILON), q));
    omega_1 = 15.0/32.0*(C + pow(tau_6/(beta_1 + EPSILON), q));
    omega_2 = 15.0/32.0*(C + pow(tau_6/(beta_2 + EPSILON), q));
    omega_3 = 1.0/32.0*(C + pow(tau_6/(beta_3 + EPSILON), q));
    
    double omega_sum = omega_0 + omega_1 + omega_2 + omega_3;
    
    omega_0 = omega_0/omega_sum;
    omega_1 = omega_1/omega_sum;
    omega_2 = omega_2/omega_sum;
    omega_3 = omega_3/omega_sum;
    
    /*
     * Compute the weights omega.
     */
    
    double R_tau = fabs(tau_6/(beta_avg + EPSILON));
    
    if (R_tau > alpha_tau)
    {
        omega_0 = sigma*omega_upwind_0 + (1.0 - sigma)*omega_0;
        omega_1 = sigma*omega_upwind_1 + (1.0 - sigma)*omega_1;
        omega_2 = sigma*omega_upwind_2 + (1.0 - sigma)*omega_2;
        omega_3 = (1.0 - sigma)*omega_3;
    }
    
    /*
     * Compute U_minus.
     */
    
    U_minus[idx_side] = 3.0/8.0*omega_0*U_array[0][idx_side] +
        (-10.0/8.0*omega_0 - 1.0/8.0*omega_1)*U_array[1][idx_side] +
        (15.0/8.0*omega_0 + 6.0/8.0*omega_1 + 3.0/8.0*omega_2)*U_array[2][idx_side] +
        (3.0/8.0*omega_1 + 6.0/8.0*omega_2 + 15.0/8.0*omega_3)*U_array[3][idx_side] +
        (-1.0/8.0*omega_2 - 10.0/8.0*omega_3)*U_array[4][idx_side] +
        3.0/8.0*omega_3*U_array[5][idx_side];
    
    /*
     * Compute beta_tilde's.
     */
    
    double beta_tilde_0, beta_tilde_1, beta_tilde_2, beta_tilde_3;
    
    computeLocalBetaTilde(&beta_tilde_0, &beta_tilde_1, &beta_tilde_2, &beta_tilde_3, U_array, idx_side);
    
    /*
     * Compute the weights omega_upwind_tilde.
     */
    
    double omega_upwind_tilde_0, omega_upwind_tilde_1, omega_upwind_tilde_2;
    
    double tau_5_tilde = fabs(beta_tilde_0 - beta_tilde_2);
    
    omega_upwind_tilde_0 = 1.0/16.0*(1.0 + pow(tau_5_tilde/(beta_tilde_0 + EPSILON), p));
    omega_upwind_tilde_1 = 5.0/8.0*(1.0 + pow(tau_5_tilde/(beta_tilde_1 + EPSILON), p));
    omega_upwind_tilde_2 = 5.0/16.0*(1.0 + pow(tau_5_tilde/(beta_tilde_2 + EPSILON), p));
    
    double omega_upwind_tilde_sum = omega_upwind_tilde_0 + omega_upwind_tilde_1 + omega_upwind_tilde_2;
    
    omega_upwind_tilde_0 = omega_upwind_tilde_0/omega_upwind_tilde_sum;
    omega_upwind_tilde_1 = omega_upwind_tilde_1/omega_upwind_tilde_sum;
    omega_upwind_tilde_2 = omega_upwind_tilde_2/omega_upwind_tilde_sum;
    
    /*
     * Compute the weights omega_central_tilde (store in omega_tilde first).
     */
    
    double omega_tilde_0, omega_tilde_1, omega_tilde_2, omega_tilde_3;
    
    double beta_avg_tilde = 1.0/8*(beta_tilde_0 + beta_tilde_2 + 6*beta_tilde_1);
    double tau_6_tilde = fabs(beta_tilde_3 - beta_avg_tilde);
    
    omega_tilde_0 = 1.0/32.0*(C + pow(tau_6_tilde/(beta_tilde_0 + EPSILON), q));
    omega_tilde_1 = 15.0/32.0*(C + pow(tau_6_tilde/(beta_tilde_1 + EPSILON), q));
    omega_tilde_2 = 15.0/32.0*(C + pow(tau_6_tilde/(beta_tilde_2 + EPSILON), q));
    omega_tilde_3 = 1.0/32.0*(C + pow(tau_6_tilde/(beta_tilde_3 + EPSILON), q));
    
    double omega_tilde_sum = omega_tilde_0 + omega_tilde_1 + omega_tilde_2 + omega_tilde_3;
    
    omega_tilde_0 = omega_tilde_0/omega_tilde_sum;
    omega_tilde_1 = omega_tilde_1/omega_tilde_sum;
    omega_tilde_2 = omega_tilde_2/omega_tilde_sum;
    omega_tilde_3 = omega_tilde_3/omega_tilde_sum;
    
    /*
     * Compute the weights omega_tilde.
     */
    
    double R_tau_tilde = fabs(tau_6_tilde/(beta_avg_tilde + EPSILON));
    
    if (R_tau_tilde > alpha_tau)
    {
        omega_tilde_0 = sigma*omega_upwind_tilde_0 + (1.0 - sigma)*omega_tilde_0;
        omega_tilde_1 = sigma*omega_upwind_tilde_1 + (1.0 - sigma)*omega_tilde_1;
        omega_tilde_2 = sigma*omega_upwind_tilde_2 + (1.0 - sigma)*omega_tilde_2;
        omega_tilde_3 = (1.0 - sigma)*omega_tilde_3;
    }
    
    U_plus[idx_side] = 3.0/8.0*omega_tilde_0*U_array[5][idx_side] +
        (-10.0/8.0*omega_tilde_0 - 1.0/8.0*omega_tilde_1)*U_array[4][idx_side] +
        (15.0/8.0*omega_tilde_0 + 6.0/8.0*omega_tilde_1 + 3.0/8.0*omega_tilde_2)*U_array[3][idx_side] +
        (3.0/8.0*omega_tilde_1 + 6.0/8.0*omega_tilde_2 + 15.0/8.0*omega_tilde_3)*U_array[2][idx_side] +
        (-1.0/8.0*omega_tilde_2 - 10.0/8.0*omega_tilde_3)*U_array[1][idx_side] +
        3.0/8.0*omega_tilde_3*U_array[0][idx_side];
}


ConvectiveFluxReconstructorWCNS6_LD_HLLC_HLL::ConvectiveFluxReconstructorWCNS6_LD_HLLC_HLL(
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
    
    d_constant_p = d_convective_flux_reconstructor_db->
        getIntegerWithDefault("constant_p", 2);
    d_constant_p = d_convective_flux_reconstructor_db->
        getIntegerWithDefault("d_constant_p", d_constant_p);
    
    d_constant_q = d_convective_flux_reconstructor_db->
        getIntegerWithDefault("constant_q", 4);
    d_constant_q = d_convective_flux_reconstructor_db->
        getIntegerWithDefault("d_constant_q", d_constant_q);
    
    d_constant_C = d_convective_flux_reconstructor_db->
        getDoubleWithDefault("constant_C", 1.0e9);
    d_constant_C = d_convective_flux_reconstructor_db->
        getDoubleWithDefault("d_constant_C", d_constant_C);
    
    d_constant_alpha_tau = d_convective_flux_reconstructor_db->
        getDoubleWithDefault("constant_alpha_tau", 35.0);
    d_constant_alpha_tau = d_convective_flux_reconstructor_db->
        getDoubleWithDefault("d_constant_alpha_tau", d_constant_alpha_tau);
}


/*
 * Print all characteristics of the convective flux reconstruction class.
 */
void
ConvectiveFluxReconstructorWCNS6_LD_HLLC_HLL::printClassData(
    std::ostream& os) const
{
    os << "\nPrint ConvectiveFluxReconstructorWCNS6_LD_HLLC_HLL object..."
       << std::endl;
    
    os << std::endl;
    
    os << "ConvectiveFluxReconstructorWCNS6_LD_HLLC_HLL: this = "
       << (ConvectiveFluxReconstructorWCNS6_LD_HLLC_HLL *)this
       << std::endl;
    os << "d_object_name = "
       << d_object_name
       << std::endl;
    os << "d_constant_p = "
       << d_constant_p
       << std::endl;
    os << "d_constant_q = "
       << d_constant_q
       << std::endl;
    os << "d_constant_C = "
       << d_constant_C
       << std::endl;
    os << "d_constant_alpha_tau = "
       << d_constant_alpha_tau
       << std::endl;
}


/*
 * Put the characteristics of the convective flux reconstruction class
 * into the restart database.
 */
void
ConvectiveFluxReconstructorWCNS6_LD_HLLC_HLL::putToRestart(
   const boost::shared_ptr<tbox::Database>& restart_db) const
{
    restart_db->putInteger("d_constant_p", d_constant_p);
    restart_db->putInteger("d_constant_q", d_constant_q);
    restart_db->putDouble("d_constant_C", d_constant_C);
    restart_db->putDouble("d_constant_alpha_tau", d_constant_alpha_tau);
}


/*
 * Perform WENO interpolation.
 */
void
ConvectiveFluxReconstructorWCNS6_LD_HLLC_HLL::performWENOInterpolation(
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
            double* U_R = variables_plus[ei]->getPointer(0);
            
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
            #pragma ivdep
#endif
            for (int i = -1; i < interior_dim_0 + 2; i++)
            {
                // Compute the linear index of the mid-point.
                const int idx_midpoint_x = i + 1;
                
                performLocalWENOInterpolation(
                    U_L,
                    U_R,
                    U_array.data(),
                    idx_midpoint_x,
                    d_constant_p,
                    d_constant_q,
                    d_constant_C,
                    d_constant_alpha_tau);
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
            double* U_R = variables_plus[ei]->getPointer(0);
            
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd safelen(1024)
#endif
                for (int i = -1; i < interior_dim_0 + 2; i++)
                {
                    // Compute the linear index of the mid-point.
                    const int idx_midpoint_x = (i + 1) +
                        (j + 1)*(interior_dim_0 + 3);
                    
                    performLocalWENOInterpolation(
                        U_L,
                        U_R,
                        U_array.data(),
                        idx_midpoint_x,
                        d_constant_p,
                        d_constant_q,
                        d_constant_C,
                        d_constant_alpha_tau);
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
            double* U_T = variables_plus[ei]->getPointer(1);
            
            for (int j = -1; j < interior_dim_1 + 2; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd safelen(1024)
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index of the mid-point.
                    const int idx_midpoint_y = (i + 1) +
                        (j + 1)*(interior_dim_0 + 2);
                    
                    performLocalWENOInterpolation(
                        U_B,
                        U_T,
                        U_array.data(),
                        idx_midpoint_y,
                        d_constant_p,
                        d_constant_q,
                        d_constant_C,
                        d_constant_alpha_tau);
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
            double* U_R = variables_plus[ei]->getPointer(0);
            
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd safelen(1024)
#endif
                    for (int i = -1; i < interior_dim_0 + 2; i++)
                    {
                        // Compute the linear index of the mid-point.
                        const int idx_midpoint_x = (i + 1) +
                            (j + 1)*(interior_dim_0 + 3) +
                            (k + 1)*(interior_dim_0 + 3)*
                                (interior_dim_1 + 2);
                        
                        performLocalWENOInterpolation(
                            U_L,
                            U_R,
                            U_array.data(),
                            idx_midpoint_x,
                            d_constant_p,
                            d_constant_q,
                            d_constant_C,
                            d_constant_alpha_tau);
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
            double* U_T = variables_plus[ei]->getPointer(1);
            
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = -1; j < interior_dim_1 + 2; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd safelen(1024)
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index of the mid-point.
                        const int idx_midpoint_y = (i + 1) +
                            (j + 1)*(interior_dim_0 + 2) +
                            (k + 1)*(interior_dim_0 + 2)*
                                (interior_dim_1 + 3);
                        
                        performLocalWENOInterpolation(
                            U_B,
                            U_T,
                            U_array.data(),
                            idx_midpoint_y,
                            d_constant_p,
                            d_constant_q,
                            d_constant_C,
                            d_constant_alpha_tau);
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
            double* U_F = variables_plus[ei]->getPointer(2);
            
            for (int k = -1; k < interior_dim_2 + 2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd safelen(1024)
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index of the mid-point.
                        const int idx_midpoint_z = (i + 1) +
                            (j + 1)*(interior_dim_0 + 2) +
                            (k + 1)*(interior_dim_0 + 2)*
                                (interior_dim_1 + 2);
                        
                        performLocalWENOInterpolation(
                            U_B,
                            U_F,
                            U_array.data(),
                            idx_midpoint_z,
                            d_constant_p,
                            d_constant_q,
                            d_constant_C,
                            d_constant_alpha_tau);
                    }
                }
            }
        }
        
    } // if (d_dim == tbox::Dimension(3))
}