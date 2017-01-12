#ifndef CONVECTIVE_FLUX_RECONSTRUCTOR_WCNS6_TEST_HPP
#define CONVECTIVE_FLUX_RECONSTRUCTOR_WCNS6_TEST_HPP

#include "SAMRAI/pdat/SideVariable.h"

#include "flow/convective_flux_reconstructors/ConvectiveFluxReconstructor.hpp"
#include "util/Directions.hpp"

#include "boost/multi_array.hpp"

#define EPSILON 1e-40

/*
 * Compute local sigma's.
 */
#ifdef HAMERS_ENABLE_SIMD
#pragma omp declare simd
#endif
inline void computeLocalSigma(
    double* sigma,
    const std::vector<double*>& U_array,
    const int& idx_side)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(U_array.size()) == 6);
#endif
    
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
#ifdef HAMERS_ENABLE_SIMD
#pragma omp declare simd
#endif
inline void computeLocalBeta(
    double* beta_0,
    double* beta_1,
    double* beta_2,
    double* beta_3,
    const std::vector<double*>& U_array,
    const int& idx_side)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(U_array.size()) == 6);
#endif
    
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
#ifdef HAMERS_ENABLE_SIMD
#pragma omp declare simd
#endif
inline void computeLocalBetaTilde(
    double* beta_tilde_0,
    double* beta_tilde_1,
    double* beta_tilde_2,
    double* beta_tilde_3,
    const std::vector<double*>& U_array,
    const int& idx_side)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(U_array.size()) == 6);
#endif
    
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


class ConvectiveFluxReconstructorWCNS6_Test: public ConvectiveFluxReconstructor
{
    public:
        ConvectiveFluxReconstructorWCNS6_Test(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_eqn,
            const int& num_species,
            const boost::shared_ptr<FlowModel>& flow_model,
            const boost::shared_ptr<tbox::Database>& convective_flux_reconstructor_db);
        
        ~ConvectiveFluxReconstructorWCNS6_Test();
        
        /*
         * Print all characteristics of the convective flux reconstruction class.
         */
        void
        printClassData(std::ostream& os) const;
        
        /*
         * Put the characteristics of the convective flux reconstruction class
         * into the restart database.
         */
        void
        putToRestart(
            const boost::shared_ptr<tbox::Database>& restart_db) const;
        
        /*
         * Compute the convective fluxes and sources due to hyperbolization
         * of the equations.
         */
        void
        computeConvectiveFluxesAndSources(
            hier::Patch& patch,
            const double time,
            const double dt,
            const int RK_step_number,
            const boost::shared_ptr<pdat::FaceVariable<double> >& variable_convective_flux,
            const boost::shared_ptr<pdat::CellVariable<double> >& variable_source,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
    private:
        /*
         * Perform WENO interpolation.
         */
        void
        performWENOInterpolation(
            std::vector<boost::shared_ptr<pdat::SideData<double> > >& variables_minus,
            std::vector<boost::shared_ptr<pdat::SideData<double> > >& variables_plus,
            const std::vector<std::vector<boost::shared_ptr<pdat::SideData<double> > > >& variables);
        
        /*
         * Constants used by the scheme.
         */
        double d_constant_C;
        int    d_constant_p;
        int    d_constant_q;
        double d_constant_alpha_tau;
        
        /*
         * Weights used in WENO interpolations.
         */
        boost::multi_array<double, 2> d_weights_c;
        
        /*
         * Projection matrix and its inversed used in WENO interpolation.
         */
        boost::multi_array<double, 2> R_inv_intercell;
        boost::multi_array<double, 2> R_intercell;
        
        /*
         * Container to store the characteristic variables in WENO interpolation.
         */
        boost::multi_array<double, 2> W_array;
        
        /*
         * Containers to store the interpolated characteristic variables in WENO interpolation.
         */
        std::vector<double> W_minus;
        std::vector<double> W_plus;
        
        /*
         * Containers used in WENO interpolation.
         */
        std::vector<double> beta;
        std::vector<double> beta_tilde;
        
        /*
         * Timers interspersed throughout the class.
         */
        
        static boost::shared_ptr<tbox::Timer> t_characteristic_decomposition;
        static boost::shared_ptr<tbox::Timer> t_WENO_interpolation;
        static boost::shared_ptr<tbox::Timer> t_Riemann_solver;
        static boost::shared_ptr<tbox::Timer> t_reconstruct_flux;
        static boost::shared_ptr<tbox::Timer> t_compute_source;
        
};


/*
 * Perform local WENO interpolation.
 */
#pragma omp declare simd
inline void performLocalWENOInterpolation(
   double* U_minus,
   double* U_plus,
   const std::vector<double*>& U_array,
   const int& idx_side)
{
    #ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(U_array.size()) == 6);
    #endif
   
    /*
     * Perform the WENO interpolation.
     */
    
    /*
    const double& C          = d_constant_C;
    const int& p             = d_constant_p;
    const int& q             = d_constant_q;
    const double& alpha_tau  = d_constant_alpha_tau;
    */
    
    const double C          = 1.0e9;
    const int p             = 2;
    const int q             = 4;
    const double alpha_tau  = 35.0;
    
    /*
     * Compute sigma.
     */
    
    const double alpha_1 = U_array[2][idx_side] - U_array[1][idx_side];
    const double alpha_2 = U_array[3][idx_side] - U_array[2][idx_side];
    const double alpha_3 = U_array[4][idx_side] - U_array[3][idx_side];
    
    const double theta_1 = fabs(alpha_1 - alpha_2)/(fabs(alpha_1) + fabs(alpha_2) + EPSILON);
    const double theta_2 = fabs(alpha_2 - alpha_3)/(fabs(alpha_2) + fabs(alpha_3) + EPSILON);
    
    double sigma = fmax(theta_1, theta_2);
    
    /*
     * Compute beta's.
     */
    
    double beta_0, beta_1, beta_2, beta_3;
    
    beta_0 = 1.0/3.0*(U_array[0][idx_side]*(4.0*U_array[0][idx_side] - 19.0*U_array[1][idx_side] +
         11.0*U_array[2][idx_side]) + U_array[1][idx_side]*(25.0*U_array[1][idx_side] -
         31.0*U_array[2][idx_side]) + 10.0*U_array[2][idx_side]*U_array[2][idx_side]);
    
    beta_1 = 1.0/3.0*(U_array[1][idx_side]*(4.0*U_array[1][idx_side] - 13.0*U_array[2][idx_side] +
         5.0*U_array[3][idx_side]) + 13.0*U_array[2][idx_side]*(U_array[2][idx_side] -
         U_array[3][idx_side]) + 4.0*U_array[3][idx_side]*U_array[3][idx_side]);
    
    beta_2 = 1.0/3.0*(U_array[2][idx_side]*(10.0*U_array[2][idx_side] - 31.0*U_array[3][idx_side] +
         11.0*U_array[4][idx_side]) + U_array[3][idx_side]*(25.0*U_array[3][idx_side] -
         19.0*U_array[4][idx_side]) + 4.0*U_array[4][idx_side]*U_array[4][idx_side]);
    
    beta_3 = 1.0/232243200.0*(U_array[0][idx_side]*(525910327.0*U_array[0][idx_side] -
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
    
    beta_tilde_0 = 1.0/3.0*(U_array[5][idx_side]*(4.0*U_array[5][idx_side] - 19.0*U_array[4][idx_side] +
         11.0*U_array[3][idx_side]) + U_array[4][idx_side]*(25.0*U_array[4][idx_side] -
         31.0*U_array[3][idx_side]) + 10.0*U_array[3][idx_side]*U_array[3][idx_side]);
    
    beta_tilde_1 = 1.0/3.0*(U_array[4][idx_side]*(4.0*U_array[4][idx_side] - 13.0*U_array[3][idx_side] +
         5.0*U_array[2][idx_side]) + 13.0*U_array[3][idx_side]*(U_array[3][idx_side] -
         U_array[2][idx_side]) + 4.0*U_array[2][idx_side]*U_array[2][idx_side]);
    
    beta_tilde_2 = 1.0/3.0*(U_array[3][idx_side]*(10.0*U_array[3][idx_side] - 31.0*U_array[2][idx_side] +
         11.0*U_array[1][idx_side]) + U_array[2][idx_side]*(25.0*U_array[2][idx_side] -
         19.0*U_array[1][idx_side]) + 4.0*U_array[1][idx_side]*U_array[1][idx_side]);
    
    beta_tilde_3 = 1.0/232243200.0*(U_array[5][idx_side]*(525910327.0*U_array[5][idx_side] -
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
     * Compute the weights omega_central_tilde  (store in omega_tilde first).
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

#endif /* CONVECTIVE_FLUX_RECONSTRUCTOR_WCNS6_TEST_HPP */
