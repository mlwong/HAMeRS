#include "flow/convective_flux_reconstructors/ConvectiveFluxReconstructor.hpp"

#define EPSILON HAMERS_EPSILON

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
 * Compute local sigma.
 */
static inline __attribute__((always_inline)) void computeLocalSigma(
    Real& sigma,
    Real** U_array,
    const int& idx_side)
{
    /*
     * Compute the sigma.
     */
    
    const Real alpha_1 = U_array[2][idx_side] - U_array[1][idx_side];
    const Real alpha_2 = U_array[3][idx_side] - U_array[2][idx_side];
    const Real alpha_3 = U_array[4][idx_side] - U_array[3][idx_side];
    
    const Real theta_1 = std::abs(alpha_1 - alpha_2)/(std::abs(alpha_1) + std::abs(alpha_2) + EPSILON);
    const Real theta_2 = std::abs(alpha_2 - alpha_3)/(std::abs(alpha_2) + std::abs(alpha_3) + EPSILON);
    
    sigma = std::max(theta_1, theta_2);
}


/*
 * Compute local beta's.
 */
static inline __attribute__((always_inline)) void computeLocalBeta6(
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
static inline __attribute__((always_inline)) void computeLocalBetaTilde6(
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
static inline __attribute__((always_inline)) void performLocalWENOInterpolationMinusLD(
    Real* U_minus,
    Real** U_array,
    const int& idx_side,
    const int& p,
    const int& q,
    const Real& C,
    const Real& alpha_tau)
{
    /*
     * Compute sigma.
     */
    
    Real sigma;
    
    computeLocalSigma(sigma, U_array, idx_side);
    
    /*
     * Compute beta's.
     */
    
    Real beta_0, beta_1, beta_2, beta_3;
    
    computeLocalBeta6(beta_0, beta_1, beta_2, beta_3, U_array, idx_side);
    
    /*
     * Compute the weights omega_upwind.
     */
    
    Real omega_upwind_0, omega_upwind_1, omega_upwind_2;
    
    Real tau_5 = std::abs(beta_0 - beta_2);
    
    omega_upwind_0 = Real(1)/Real(16)*(Real(1) + ipow(tau_5/(beta_0 + EPSILON), p));
    omega_upwind_1 = Real(5)/Real(8)*(Real(1) + ipow(tau_5/(beta_1 + EPSILON), p));
    omega_upwind_2 = Real(5)/Real(16)*(Real(1) + ipow(tau_5/(beta_2 + EPSILON), p));
    
    Real omega_upwind_sum = omega_upwind_0 + omega_upwind_1 + omega_upwind_2;
    
    omega_upwind_0 = omega_upwind_0/omega_upwind_sum;
    omega_upwind_1 = omega_upwind_1/omega_upwind_sum;
    omega_upwind_2 = omega_upwind_2/omega_upwind_sum;
    
    /*
     * Compute the weights omega_central (store in omega first).
     */
    
    Real omega_0, omega_1, omega_2, omega_3;
    
    Real beta_avg = Real(1)/Real(8)*(beta_0 + beta_2 + Real(6)*beta_1);
    Real tau_6 = std::abs(beta_3 - beta_avg);
    
    omega_0 = Real(1)/Real(32)*(C + ipow(tau_6/(beta_0 + EPSILON), q));
    omega_1 = Real(15)/Real(32)*(C + ipow(tau_6/(beta_1 + EPSILON), q));
    omega_2 = Real(15)/Real(32)*(C + ipow(tau_6/(beta_2 + EPSILON), q));
    omega_3 = Real(1)/Real(32)*(C + ipow(tau_6/(beta_3 + EPSILON), q));
    
    Real omega_sum = omega_0 + omega_1 + omega_2 + omega_3;
    
    omega_0 = omega_0/omega_sum;
    omega_1 = omega_1/omega_sum;
    omega_2 = omega_2/omega_sum;
    omega_3 = omega_3/omega_sum;
    
    /*
     * Compute the weights omega.
     */
    
    Real R_tau = tau_6/(beta_avg + EPSILON);
    
    if (R_tau > alpha_tau)
    {
        omega_0 = sigma*omega_upwind_0 + (Real(1) - sigma)*omega_0;
        omega_1 = sigma*omega_upwind_1 + (Real(1) - sigma)*omega_1;
        omega_2 = sigma*omega_upwind_2 + (Real(1) - sigma)*omega_2;
        omega_3 = (Real(1) - sigma)*omega_3;
    }
    
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
static inline __attribute__((always_inline)) void performLocalWENOInterpolationPlusLD(
    Real* U_plus,
    Real** U_array,
    const int& idx_side,
    const int& p,
    const int& q,
    const Real& C,
    const Real& alpha_tau)
{
    /*
     * Compute sigma.
     */
    
    Real sigma;
    
    computeLocalSigma(sigma, U_array, idx_side);
    
    /*
     * Compute beta_tilde's.
     */
    
    Real beta_tilde_0, beta_tilde_1, beta_tilde_2, beta_tilde_3;
    
    computeLocalBetaTilde6(beta_tilde_0, beta_tilde_1, beta_tilde_2, beta_tilde_3, U_array, idx_side);
    
    /*
     * Compute the weights omega_upwind_tilde.
     */
    
    Real omega_upwind_tilde_0, omega_upwind_tilde_1, omega_upwind_tilde_2;
    
    Real tau_5_tilde = std::abs(beta_tilde_0 - beta_tilde_2);
    
    omega_upwind_tilde_0 = Real(1)/Real(16)*(Real(1) + ipow(tau_5_tilde/(beta_tilde_0 + EPSILON), p));
    omega_upwind_tilde_1 = Real(5)/Real(8)*(Real(1) + ipow(tau_5_tilde/(beta_tilde_1 + EPSILON), p));
    omega_upwind_tilde_2 = Real(5)/Real(16)*(Real(1) + ipow(tau_5_tilde/(beta_tilde_2 + EPSILON), p));
    
    Real omega_upwind_tilde_sum = omega_upwind_tilde_0 + omega_upwind_tilde_1 + omega_upwind_tilde_2;
    
    omega_upwind_tilde_0 = omega_upwind_tilde_0/omega_upwind_tilde_sum;
    omega_upwind_tilde_1 = omega_upwind_tilde_1/omega_upwind_tilde_sum;
    omega_upwind_tilde_2 = omega_upwind_tilde_2/omega_upwind_tilde_sum;
    
    /*
     * Compute the weights omega_central_tilde (store in omega_tilde first).
     */
    
    Real omega_tilde_0, omega_tilde_1, omega_tilde_2, omega_tilde_3;
    
    Real beta_avg_tilde = Real(1)/Real(8)*(beta_tilde_0 + beta_tilde_2 + Real(6)*beta_tilde_1);
    Real tau_6_tilde = std::abs(beta_tilde_3 - beta_avg_tilde);
    
    omega_tilde_0 = Real(1)/Real(32)*(C + ipow(tau_6_tilde/(beta_tilde_0 + EPSILON), q));
    omega_tilde_1 = Real(15)/Real(32)*(C + ipow(tau_6_tilde/(beta_tilde_1 + EPSILON), q));
    omega_tilde_2 = Real(15)/Real(32)*(C + ipow(tau_6_tilde/(beta_tilde_2 + EPSILON), q));
    omega_tilde_3 = Real(1)/Real(32)*(C + ipow(tau_6_tilde/(beta_tilde_3 + EPSILON), q));
    
    Real omega_tilde_sum = omega_tilde_0 + omega_tilde_1 + omega_tilde_2 + omega_tilde_3;
    
    omega_tilde_0 = omega_tilde_0/omega_tilde_sum;
    omega_tilde_1 = omega_tilde_1/omega_tilde_sum;
    omega_tilde_2 = omega_tilde_2/omega_tilde_sum;
    omega_tilde_3 = omega_tilde_3/omega_tilde_sum;
    
    /*
     * Compute the weights omega_tilde.
     */
    
    Real R_tau_tilde = tau_6_tilde/(beta_avg_tilde + EPSILON);
    
    if (R_tau_tilde > alpha_tau)
    {
        omega_tilde_0 = sigma*omega_upwind_tilde_0 + (Real(1) - sigma)*omega_tilde_0;
        omega_tilde_1 = sigma*omega_upwind_tilde_1 + (Real(1) - sigma)*omega_tilde_1;
        omega_tilde_2 = sigma*omega_upwind_tilde_2 + (Real(1) - sigma)*omega_tilde_2;
        omega_tilde_3 = (Real(1) - sigma)*omega_tilde_3;
    }
    
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


/*
 * Compute local beta's.
 */
static inline __attribute__((always_inline)) void computeLocalBeta5(
    Real& beta_0,
    Real& beta_1,
    Real& beta_2,
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
}


/*
 * Compute local beta_tilde's.
 */
static inline __attribute__((always_inline)) void computeLocalBetaTilde5(
    Real& beta_tilde_0,
    Real& beta_tilde_1,
    Real& beta_tilde_2,
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
}


/*
 * Perform local WENO interpolation of U_minus.
 */
static inline __attribute__((always_inline)) void performLocalWENOInterpolationMinusZ(
    Real* U_minus,
    Real** U_array,
    const int& idx_side,
    const int& p)
{
    /*
     * Compute beta's.
     */
    
    Real beta_0, beta_1, beta_2;
    
    computeLocalBeta5(beta_0, beta_1, beta_2, U_array, idx_side);
    
    /*
     * Compute the weights omega.
     */
    
    Real omega_0, omega_1, omega_2;
    
    Real tau_5 = std::abs(beta_0 - beta_2);
    
    omega_0 = Real(1)/Real(16)*(Real(1) + ipow(tau_5/(beta_0 + EPSILON), p));
    omega_1 = Real(5)/Real(8)*(Real(1) + ipow(tau_5/(beta_1 + EPSILON), p));
    omega_2 = Real(5)/Real(16)*(Real(1) + ipow(tau_5/(beta_2 + EPSILON), p));
    
    Real omega_sum = omega_0 + omega_1 + omega_2;
    
    omega_0 = omega_0/omega_sum;
    omega_1 = omega_1/omega_sum;
    omega_2 = omega_2/omega_sum;
    
    /*
     * Compute U_minus.
     */
    
    U_minus[idx_side] = Real(3)/Real(8)*omega_0*U_array[0][idx_side] +
        (-Real(10)/Real(8)*omega_0 - Real(1)/Real(8)*omega_1)*U_array[1][idx_side] +
        (Real(15)/Real(8)*omega_0 + Real(6)/Real(8)*omega_1 +
        Real(3)/Real(8)*omega_2)*U_array[2][idx_side] +
        (Real(3)/Real(8)*omega_1 + Real(6)/Real(8)*omega_2)*U_array[3][idx_side] -
        Real(1)/Real(8)*omega_2*U_array[4][idx_side];
}


/*
 * Perform local WENO interpolation of U_plus.
 */
static inline __attribute__((always_inline)) void performLocalWENOInterpolationPlusZ(
    Real* U_plus,
    Real** U_array,
    const int& idx_side,
    const int& p)
{
    /*
     * Compute beta_tilde's.
     */
    
    Real beta_tilde_0, beta_tilde_1, beta_tilde_2;
    
    computeLocalBetaTilde5(beta_tilde_0, beta_tilde_1, beta_tilde_2, U_array, idx_side);
    
    /*
     * Compute the weights omega_upwind_tilde.
     */
    
    Real omega_tilde_0, omega_tilde_1, omega_tilde_2;
    
    Real tau_5_tilde = std::abs(beta_tilde_0 - beta_tilde_2);
    
    omega_tilde_0 = Real(1)/Real(16)*(Real(1) + ipow(tau_5_tilde/(beta_tilde_0 + EPSILON), p));
    omega_tilde_1 = Real(5)/Real(8)*(Real(1) + ipow(tau_5_tilde/(beta_tilde_1 + EPSILON), p));
    omega_tilde_2 = Real(5)/Real(16)*(Real(1) + ipow(tau_5_tilde/(beta_tilde_2 + EPSILON), p));
    
    Real omega_tilde_sum = omega_tilde_0 + omega_tilde_1 + omega_tilde_2;
    
    omega_tilde_0 = omega_tilde_0/omega_tilde_sum;
    omega_tilde_1 = omega_tilde_1/omega_tilde_sum;
    omega_tilde_2 = omega_tilde_2/omega_tilde_sum;
    
    /*
     * Compute U_plus.
     */
    
    U_plus[idx_side] = Real(3)/Real(8)*omega_tilde_0*U_array[5][idx_side] +
        (-Real(10)/Real(8)*omega_tilde_0 - Real(1)/Real(8)*omega_tilde_1)*U_array[4][idx_side] +
        (Real(15)/Real(8)*omega_tilde_0 + Real(6)/Real(8)*omega_tilde_1 +
        Real(3)/Real(8)*omega_tilde_2)*U_array[3][idx_side] +
        (Real(3)/Real(8)*omega_tilde_1 + Real(6)/Real(8)*omega_tilde_2)*U_array[2][idx_side] -
        Real(1)/Real(8)*omega_tilde_2*U_array[1][idx_side];
}


ConvectiveFluxReconstructor::ConvectiveFluxReconstructor(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_eqn,
    const FLOW_MODEL::TYPE& flow_model_type,
    const HAMERS_SHARED_PTR<FlowModel>& flow_model,
    const HAMERS_SHARED_PTR<tbox::Database>& convective_flux_reconstructor_db):
        d_object_name(object_name),
        d_dim(dim),
        d_grid_geometry(grid_geometry),
        d_num_conv_ghosts(hier::IntVector::getZero(d_dim)),
        d_num_eqn(num_eqn),
        d_flow_model_type(flow_model_type),
        d_flow_model(flow_model),
        d_convective_flux_reconstructor_db(convective_flux_reconstructor_db),
        d_num_ghosts_shock_interface_capturing(3)
{
    d_threshold_sensor_shock = Real(3)/Real(10);
    d_threshold_sensor_interface = Real(5)/Real(10000);
    
    d_threshold_sensor_shock = d_convective_flux_reconstructor_db->getRealWithDefault("threshold_sensor_shock", d_threshold_sensor_shock);
    d_threshold_sensor_shock = d_convective_flux_reconstructor_db->getRealWithDefault("d_threshold_sensor_shock", d_threshold_sensor_shock);
    
    d_threshold_sensor_interface = d_convective_flux_reconstructor_db->getRealWithDefault("threshold_sensor_interface", d_threshold_sensor_interface);
    d_threshold_sensor_interface = d_convective_flux_reconstructor_db->getRealWithDefault("d_threshold_sensor_interface", d_threshold_sensor_interface);
    
    
    std::string weno_interp_str = "WENO5Z";
    
    weno_interp_str = d_convective_flux_reconstructor_db->getStringWithDefault("weno_interp", weno_interp_str);
    weno_interp_str = d_convective_flux_reconstructor_db->getStringWithDefault("d_weno_interp", weno_interp_str);
    
    if (weno_interp_str == "WENO5Z")
    {
        d_weno_interp = WENO_INTERP::WENO5Z;
    }
    else if (weno_interp_str == "WENO6LD")
    {
        d_weno_interp = WENO_INTERP::WENO6LD;
    }
    else
    {
        TBOX_ERROR(d_object_name << ": "
            << "Unknown WENO interpolation method specified in input database"
            << std::endl);
    }
    
    d_use_MND_finite_differencing = d_convective_flux_reconstructor_db->getBoolWithDefault(
        "use_MND_finite_differencing", false);
    d_use_MND_finite_differencing = d_convective_flux_reconstructor_db->getBoolWithDefault(
        "d_use_MND_finite_differencing", d_use_MND_finite_differencing);
    
    if (d_use_MND_finite_differencing)
    {
        d_num_ghosts_shock_interface_capturing = 4;
    }
    else
    {
        d_num_ghosts_shock_interface_capturing = 3;
    }
    
    d_eqn_form = d_flow_model->getEquationsForm();
    d_has_advective_eqn_form = false;
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        if (d_eqn_form[ei] == EQN_FORM::ADVECTIVE)
        {
            d_has_advective_eqn_form = true;
        }
    }
}


/*
 * Put the characteristics of the base convective flux reconstruction class into the restart database.
 */
void
ConvectiveFluxReconstructor::putToRestartBase(
    const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const
{
    restart_db->putReal("d_threshold_sensor_shock", d_threshold_sensor_shock);
    restart_db->putReal("d_threshold_sensor_interface", d_threshold_sensor_interface);
    
    if (d_weno_interp == WENO_INTERP::WENO5Z) {
        restart_db->putString("d_weno_interp", "WENO5Z");
    }
    else if (d_weno_interp == WENO_INTERP::WENO6LD) {
        restart_db->putString("d_weno_interp", "WENO6LD");
    }
    
    restart_db->putBool("d_use_MND_finite_differencing", d_use_MND_finite_differencing);
}


/*
 * (Old) Compute the convective flux and source due to splitting using shock-capturing scheme.
 */
void
ConvectiveFluxReconstructor::computeConvectiveFluxAndSourceOnPatchShockCapturingOld(
    hier::Patch& patch,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& convective_flux,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& source_scratch,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const hier::Box& domain,
    const double dt,
    const bool use_shock_capturing,
    const bool use_interface_capturing) const
{
    if (!use_shock_capturing && !use_interface_capturing)
    {
        return;
    }
    
    d_flow_model->setupRiemannSolver();
    d_flow_model->setupBasicUtilities();
    
    HAMERS_SHARED_PTR<FlowModelRiemannSolver> riemann_solver = d_flow_model->getFlowModelRiemannSolver();
    HAMERS_SHARED_PTR<FlowModelBasicUtilities> basic_utilities = d_flow_model->getFlowModelBasicUtilities();
    
    // Get the grid spacing.
    const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
        HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
    const double* const dx = patch_geom->getDx();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    domain_lo = domain.lower() - interior_box.lower();
    domain_dims = domain.numberCells();
    
    // Create domains in different directions.
    std::vector<hier::Box> domains(d_dim.getValue(), domain);
    hier::Index lower = domain.lower();
    hier::Index upper = domain.upper();
    lower[0] = lower[0] - 1;
    upper[0] = upper[0] + 1;
    domains[0].setLower(lower);
    domains[0].setUpper(upper);
    if (d_dim > tbox::Dimension(1))
    {
        lower = domain.lower();
        upper = domain.upper();
        lower[1] = lower[1] - 1;
        upper[1] = upper[1] + 1;
        domains[1].setLower(lower);
        domains[1].setUpper(upper);
    }
    if (d_dim > tbox::Dimension(2))
    {
        lower = domain.lower();
        upper = domain.upper();
        lower[2] = lower[2] - 1;
        upper[2] = upper[2] + 1;
        domains[2].setLower(lower);
        domains[2].setUpper(upper);
    }
    
    // Allocate temporary patch data.
    HAMERS_SHARED_PTR<pdat::SideData<Real> > velocity_midpoint;
    
    if (d_has_advective_eqn_form)
    {
        velocity_midpoint.reset(new pdat::SideData<Real>(
            interior_box, d_dim.getValue(), hier::IntVector::getOne(d_dim)));
    }
    
    HAMERS_SHARED_PTR<pdat::SideData<Real> > convective_flux_midpoint(
        new pdat::SideData<Real>(interior_box, d_num_eqn, hier::IntVector::getOne(d_dim)));
    
    HAMERS_SHARED_PTR<pdat::SideData<Real> > discontinuity_sensor_side(
        new pdat::SideData<Real>(interior_box, 1, hier::IntVector::getZero(d_dim)));
    
    HAMERS_SHARED_PTR<pdat::CellData<Real> > discontinuity_sensor_cell(
        new pdat::CellData<Real>(interior_box, 1, hier::IntVector::getZero(d_dim)));
    
    HAMERS_SHARED_PTR<pdat::CellData<Real> > velocity_derivatives;
    HAMERS_SHARED_PTR<pdat::CellData<Real> > dilatation;
    HAMERS_SHARED_PTR<pdat::CellData<Real> > enstrophy;
    std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > density_sensors;
    
    bool perform_WCNS = false;
    
    if (d_dim > tbox::Dimension(1))
    {
        velocity_derivatives.reset(new pdat::CellData<Real>(
            interior_box, d_dim.getValue()*d_dim.getValue(), hier::IntVector::getOne(d_dim)));
        
        dilatation.reset(new pdat::CellData<Real>(
            interior_box, 1, hier::IntVector::getOne(d_dim)));
        
        enstrophy.reset(new pdat::CellData<Real>(
            interior_box, 1, hier::IntVector::getOne(d_dim)));
        
        density_sensors.resize(d_dim.getValue());
        for (int di = 0; di < d_dim.getValue(); di++)
        {
            density_sensors[di].reset(new pdat::CellData<Real>(
                interior_box, 1, hier::IntVector::getOne(d_dim)));
        }
    }
    
    /*
     * Register the patch and derived cell variables in the flow model and compute the corresponding cell data.
     */
    
    d_flow_model->registerPatchWithDataContext(patch, data_context);
    
    std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
    
    num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("DENSITY", d_num_conv_ghosts));
    num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("VELOCITY", d_num_conv_ghosts));
    num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("CONVECTIVE_FLUX_X", d_num_conv_ghosts));
    if (d_dim > tbox::Dimension(1))
    {
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("CONVECTIVE_FLUX_Y", d_num_conv_ghosts));
    }
    if (d_dim > tbox::Dimension(2))
    {
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("CONVECTIVE_FLUX_Z", d_num_conv_ghosts));
    }
    num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("PRIMITIVE_VARIABLES", d_num_conv_ghosts));
    
    d_flow_model->registerDerivedVariables(num_subghosts_of_data);
    
    basic_utilities->registerDerivedVariablesForCharacteristicProjectionOfPrimitiveVariables(
        d_num_conv_ghosts,
        AVERAGING::SIMPLE);
    
    d_flow_model->allocateMemoryForDerivedCellData();
    
    d_flow_model->computeDerivedCellData();
    
    HAMERS_SHARED_PTR<pdat::CellData<Real> > density = d_flow_model->getCellData("DENSITY");
    HAMERS_SHARED_PTR<pdat::CellData<Real> > velocity = d_flow_model->getCellData("VELOCITY");
    
    std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > convective_flux_node(d_dim.getValue());
    convective_flux_node[0] = d_flow_model->getCellData("CONVECTIVE_FLUX_X");
    if (d_dim > tbox::Dimension(1))
    {
        convective_flux_node[1] = d_flow_model->getCellData("CONVECTIVE_FLUX_Y");
    }
    if (d_dim > tbox::Dimension(2))
    {
        convective_flux_node[2] = d_flow_model->getCellData("CONVECTIVE_FLUX_Z");
    }
    
    /*
     * Get the pointers to the conservative variables and primitive variables.
     * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
     */
    
    std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > conservative_variables =
        d_flow_model->getCellDataOfConservativeVariables();
    
    std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > primitive_variables =
        d_flow_model->getCellDataOfPrimitiveVariables();
    
    std::vector<hier::IntVector> num_subghosts_conservative_var;
    num_subghosts_conservative_var.reserve(d_num_eqn);
    
    std::vector<hier::IntVector> num_subghosts_primitive_var;
        num_subghosts_primitive_var.reserve(d_num_eqn);
    
    std::vector<hier::IntVector> subghostcell_dims_conservative_var;
    subghostcell_dims_conservative_var.reserve(d_num_eqn);
    
    std::vector<hier::IntVector> subghostcell_dims_primitive_var;
        subghostcell_dims_primitive_var.reserve(d_num_eqn);
    
    std::vector<Real*> Q;
    Q.reserve(d_num_eqn);
    
    std::vector<Real*> V;
    V.reserve(d_num_eqn);
    
    int count_eqn = 0;
    
    for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
    {
        int depth = conservative_variables[vi]->getDepth();
        
        for (int di = 0; di < depth; di++)
        {
            // If the last element of the conservative variable vector is not in the system of equations,
            // ignore it.
            if (count_eqn >= d_num_eqn)
                break;
            
            Q.push_back(conservative_variables[vi]->getPointer(di));
            num_subghosts_conservative_var.push_back(conservative_variables[vi]->getGhostCellWidth());
            subghostcell_dims_conservative_var.push_back(
                conservative_variables[vi]->getGhostBox().numberCells());
            
            count_eqn++;
        }
    }
    
    count_eqn = 0;
    
    for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
    {
        int depth = primitive_variables[vi]->getDepth();
        
        for (int di = 0; di < depth; di++)
        {
            // If the last element of the primitive variable vector is not in the system of equations,
            // ignore it.
            if (count_eqn >= d_num_eqn)
                break;
            
            V.push_back(primitive_variables[vi]->getPointer(di));
            num_subghosts_primitive_var.push_back(primitive_variables[vi]->getGhostCellWidth());
            subghostcell_dims_primitive_var.push_back(
                primitive_variables[vi]->getGhostBox().numberCells());
            
            count_eqn++;
        }
    }
    
    /*
     * Pointers to shock sensors.
     */
    
    Real* s   = discontinuity_sensor_cell->getPointer(0);
    Real* s_x = discontinuity_sensor_side->getPointer(0);
    Real* s_y = nullptr;
    Real* s_z = nullptr;
    
    if (d_dim == tbox::Dimension(1))
    {
        discontinuity_sensor_cell->fillAll(Real(1));
        discontinuity_sensor_side->fillAll(Real(1));
        perform_WCNS = true;
    }
    else
    {
        /*
         * Get the numbers of ghost cells of the variables.
         */
        
        const hier::IntVector num_ghosts_density = density->getGhostCellWidth();
        
        /*
         * Get the ghost cell dimensions of of the variables.
         */
        
        const hier::IntVector ghostcell_dims_density = density->getGhostBox().numberCells();
        
        // Get the pointers to the data.
        Real* rho   = density->getPointer(0);
        Real* theta = dilatation->getPointer(0);
        Real* Omega = enstrophy->getPointer(0);
        Real* rho_s_x = density_sensors[0]->getPointer(0);
        Real* rho_s_y = density_sensors[1]->getPointer(0);
        
        s_y = discontinuity_sensor_side->getPointer(1);
        
        const hier::Box empty_box(d_dim);
        
        if (d_dim == tbox::Dimension(2))
        {
            /*
             * Get the local lower indices and the number of cells in each dimension.
             */
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_lo_1 = domain_lo[1];
            const int domain_dim_0 = domain_dims[0];
            const int domain_dim_1 = domain_dims[1];
            
            const int num_ghosts_0_density = num_ghosts_density[0];
            const int num_ghosts_1_density = num_ghosts_density[1];
            const int ghostcell_dim_0_density = ghostcell_dims_density[0];
            
            /*
             * Get the interior dimensions.
             */
            
            const int interior_dim_0 = interior_dims[0];
            const int interior_dim_1 = interior_dims[1];
            
            /*
             * Compute the derivatives of velocity, dilatation and vorticity magnitude.
             */
            
            HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_x(
                new DerivativeFirstOrder("first order derivative in x-direction", d_dim, DIRECTION::X_DIRECTION, 1));
            
            HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_y(
                new DerivativeFirstOrder("first order derivative in y-direction", d_dim, DIRECTION::Y_DIRECTION, 1));
            
            // Compute dudx.
            derivative_first_order_x->computeDerivative(
                velocity_derivatives,
                velocity,
                Real(dx[0]),
                empty_box,
                0,
                0);
            
            // Compute dudy.
            derivative_first_order_y->computeDerivative(
                velocity_derivatives,
                velocity,
                Real(dx[1]),
                empty_box,
                1,
                0);
            
            // Compute dvdx.
            derivative_first_order_x->computeDerivative(
                velocity_derivatives,
                velocity,
                Real(dx[0]),
                empty_box,
                2,
                1);
            
            // Compute dvdy.
            derivative_first_order_y->computeDerivative(
                velocity_derivatives,
                velocity,
                Real(dx[1]),
                empty_box,
                3,
                1);
            
            // Get the pointers to the cell data of velocity derivatives.
            Real* dudx = velocity_derivatives->getPointer(0);
            Real* dudy = velocity_derivatives->getPointer(1);
            Real* dvdx = velocity_derivatives->getPointer(2);
            Real* dvdy = velocity_derivatives->getPointer(3);
            
            // Compute the dilatation and the enstrophy.
            for (int j = domain_lo_1 - 1; j < domain_lo_1 + domain_dim_1 + 1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0 - 1; i < domain_lo_0 + domain_dim_0 + 1; i++)
                {
                    // Compute the linear indices.
                    const int idx = (i + 1) +
                        (j + 1)*(interior_dim_0 + 2);
                    
                    const int idx_rho = (i + num_ghosts_0_density) +
                        (j + num_ghosts_1_density)*ghostcell_dim_0_density;
                    
                    const int idx_rho_x_L = (i - 1 + num_ghosts_0_density) +
                        (j + num_ghosts_1_density)*ghostcell_dim_0_density;
                    
                    const int idx_rho_x_R = (i + 1 + num_ghosts_0_density) +
                        (j + num_ghosts_1_density)*ghostcell_dim_0_density;
                    
                    const int idx_rho_y_B = (i + num_ghosts_0_density) +
                        (j - 1 + num_ghosts_1_density)*ghostcell_dim_0_density;
                    
                    const int idx_rho_y_T = (i + num_ghosts_0_density) +
                        (j + 1 + num_ghosts_1_density)*ghostcell_dim_0_density;
                    
                    theta[idx] = dudx[idx] + dvdy[idx];
                    Omega[idx] = (dvdx[idx] - dudy[idx])*(dvdx[idx] - dudy[idx]);
                    
                    rho_s_x[idx] = std::abs(rho[idx_rho_x_R] - Real(2)*rho[idx_rho] + rho[idx_rho_x_L])/
                        (rho[idx_rho_x_R] + Real(2)*rho[idx_rho] + rho[idx_rho_x_L] + EPSILON);
                    
                    rho_s_y[idx] = std::abs(rho[idx_rho_y_T] - Real(2)*rho[idx_rho] + rho[idx_rho_y_B])/
                        (rho[idx_rho_y_T] + Real(2)*rho[idx_rho] + rho[idx_rho_y_B] + EPSILON);
                }
            }
            
            int count_valid_sensor = 0;
            
            const Real half = Real(1)/Real(2);
            
            // Compute the shock sensor on the side in the x-direction.
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                {
                    // Compute the linear indices.
                    const int idx_sensor = i +
                        j*(interior_dim_0 + 1);
                    
                    const int idx_L = (i - 1 + 1) +
                        (j + 1)*(interior_dim_0 + 2);
                    
                    const int idx_R = (i + 0 + 1) +
                        (j + 1)*(interior_dim_0 + 2);
                    
                    const Real Ducros_value_L = (-theta[idx_L]*std::abs(theta[idx_L]))/
                        (theta[idx_L]*theta[idx_L] + Omega[idx_L] + EPSILON);
                    
                    const Real Ducros_value_R = (-theta[idx_R]*std::abs(theta[idx_R]))/
                        (theta[idx_R]*theta[idx_R] + Omega[idx_R] + EPSILON);
                    
                    const Real Ducros_value_midpoint = half*(Ducros_value_L + Ducros_value_R);
                    
                    s_x[idx_sensor] = Real(0);
                    if (Ducros_value_midpoint > d_threshold_sensor_shock && use_shock_capturing)
                    {
                        s_x[idx_sensor] = Real(1);
                    }
                    
                    const Real rho_s_x_midpoint = half*(rho_s_x[idx_L] + rho_s_x[idx_R]);
                    if (rho_s_x_midpoint > d_threshold_sensor_interface && use_interface_capturing)
                    {
                        s_x[idx_sensor] = Real(1);
                    }
                    
                    if (s_x[idx_sensor] > Real(0))
                    {
                        count_valid_sensor++;
                    }
                }
            }
            
            // Compute the shock sensor on the side in the y-direction.
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_sensor = i +
                        j*interior_dim_0;
                    
                    const int idx_B = (i + 1) +
                        (j - 1 + 1)*(interior_dim_0 + 2);
                    
                    const int idx_T = (i + 1) +
                        (j + 0 + 1)*(interior_dim_0 + 2);
                    
                    const Real Ducros_value_B = (-theta[idx_B]*std::abs(theta[idx_B]))/
                        (theta[idx_B]*theta[idx_B] + Omega[idx_B] + EPSILON);
                    
                    const Real Ducros_value_T = (-theta[idx_T]*std::abs(theta[idx_T]))/
                        (theta[idx_T]*theta[idx_T] + Omega[idx_T] + EPSILON);
                    
                    const Real Ducros_value_midpoint = half*(Ducros_value_B + Ducros_value_T);
                    
                    s_y[idx_sensor] = Real(0);
                    if (Ducros_value_midpoint > d_threshold_sensor_shock && use_shock_capturing)
                    {
                        s_y[idx_sensor] = Real(1);
                    }
                    
                    const Real rho_s_y_midpoint = half*(rho_s_y[idx_B] + rho_s_y[idx_T]);
                    if (rho_s_y_midpoint > d_threshold_sensor_interface && use_interface_capturing)
                    {
                        s_y[idx_sensor] = Real(1);
                    }
                    
                    if (s_y[idx_sensor] > Real(0))
                    {
                        count_valid_sensor++;
                    }
                }
            }
            
            if (d_has_advective_eqn_form)
            {
                // Compute the cell-based Ducros-like shock sensor.
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_sensor = i +
                            j*interior_dim_0;
                        
                        const int idx = (i + 1) +
                            (j + 1)*(interior_dim_0 + 2);
                        
                        const Real Ducros_value = (-theta[idx]*std::abs(theta[idx]))/
                            (theta[idx]*theta[idx] + Omega[idx] + EPSILON);
                        
                        s[idx_sensor] = Real(0);
                        if (Ducros_value > d_threshold_sensor_shock && use_shock_capturing)
                        {
                            s[idx_sensor] = Real(1);
                        }
                        {
                            s[idx_sensor] = Real(1);
                        }
                        if (rho_s_x[idx] > d_threshold_sensor_interface && use_interface_capturing)
                        {
                            s[idx_sensor] = Real(1);
                        }
                        if (rho_s_y[idx] > d_threshold_sensor_interface && use_interface_capturing)
                        {
                            s[idx_sensor] = Real(1);
                        }
                        
                        if (s[idx_sensor] > Real(0))
                        {
                            count_valid_sensor++;
                        }
                    }
                }
            }
            
            if (count_valid_sensor > 0)
            {
                perform_WCNS = true;
            }
        }
        else if (d_dim == tbox::Dimension(3))
        {
            Real* rho_s_z = density_sensors[2]->getPointer(0);
            Real* s_z = discontinuity_sensor_side->getPointer(2);
            
            /*
             * Get the local lower indices and the number of cells in each dimension.
             */
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_lo_1 = domain_lo[1];
            const int domain_lo_2 = domain_lo[2];
            const int domain_dim_0 = domain_dims[0];
            const int domain_dim_1 = domain_dims[1];
            const int domain_dim_2 = domain_dims[2];
            
            const int num_ghosts_0_density = num_ghosts_density[0];
            const int num_ghosts_1_density = num_ghosts_density[1];
            const int num_ghosts_2_density = num_ghosts_density[2];
            const int ghostcell_dim_0_density = ghostcell_dims_density[0];
            const int ghostcell_dim_1_density = ghostcell_dims_density[1];
            
            /*
             * Get the interior dimensions.
             */
            
            const int interior_dim_0 = interior_dims[0];
            const int interior_dim_1 = interior_dims[1];
            const int interior_dim_2 = interior_dims[2];
            
            /*
             * Compute the derivatives of velocity, dilatation and vorticity magnitude.
             */
            
            HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_x(
                new DerivativeFirstOrder("first order derivative in x-direction", d_dim, DIRECTION::X_DIRECTION, 1));
            
            HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_y(
                new DerivativeFirstOrder("first order derivative in y-direction", d_dim, DIRECTION::Y_DIRECTION, 1));
            
            HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_z(
                new DerivativeFirstOrder("first order derivative in z-direction", d_dim, DIRECTION::Z_DIRECTION, 1));
            
            // Compute dudx.
            derivative_first_order_x->computeDerivative(
                velocity_derivatives,
                velocity,
                Real(dx[0]),
                empty_box,
                0,
                0);
            
            // Compute dudy.
            derivative_first_order_y->computeDerivative(
                velocity_derivatives,
                velocity,
                Real(dx[1]),
                empty_box,
                1,
                0);
            
            // Compute dudz.
            derivative_first_order_z->computeDerivative(
                velocity_derivatives,
                velocity,
                Real(dx[2]),
                empty_box,
                2,
                0);
            
            // Compute dvdx.
            derivative_first_order_x->computeDerivative(
                velocity_derivatives,
                velocity,
                Real(dx[0]),
                empty_box,
                3,
                1);
            
            // Compute dvdy.
            derivative_first_order_y->computeDerivative(
                velocity_derivatives,
                velocity,
                Real(dx[1]),
                empty_box,
                4,
                1);
            
            // Compute dvdz.
            derivative_first_order_z->computeDerivative(
                velocity_derivatives,
                velocity,
                Real(dx[2]),
                empty_box,
                5,
                1);
            
            // Compute dwdx.
            derivative_first_order_x->computeDerivative(
                velocity_derivatives,
                velocity,
                Real(dx[0]),
                empty_box,
                6,
                2);
            
            // Compute dwdy.
            derivative_first_order_y->computeDerivative(
                velocity_derivatives,
                velocity,
                Real(dx[1]),
                empty_box,
                7,
                2);
            
            // Compute dwdz.
            derivative_first_order_z->computeDerivative(
                velocity_derivatives,
                velocity,
                Real(dx[2]),
                empty_box,
                8,
                2);
            
            // Get the pointers to the cell data of velocity derivatives.
            Real* dudx = velocity_derivatives->getPointer(0);
            Real* dudy = velocity_derivatives->getPointer(1);
            Real* dudz = velocity_derivatives->getPointer(2);
            Real* dvdx = velocity_derivatives->getPointer(3);
            Real* dvdy = velocity_derivatives->getPointer(4);
            Real* dvdz = velocity_derivatives->getPointer(5);
            Real* dwdx = velocity_derivatives->getPointer(6);
            Real* dwdy = velocity_derivatives->getPointer(7);
            Real* dwdz = velocity_derivatives->getPointer(8);
            
            // Compute the dilatation and the enstrophy.
            for (int k = domain_lo_2 - 1; k < domain_lo_2 + domain_dim_2 + 1; k++)
            {
                for (int j = domain_lo_1 - 1; j < domain_lo_1 + domain_dim_1 + 1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0 - 1; i < domain_lo_0 + domain_dim_0 + 1; i++)
                    {
                        // Compute the linear index.
                        const int idx = (i + 1) +
                            (j + 1)*(interior_dim_0 + 2) +
                            (k + 1)*(interior_dim_0 + 2)*
                                (interior_dim_1 + 2);
                        
                        theta[idx] = dudx[idx] + dvdy[idx] + dwdz[idx];
                        
                        const Real omega_x = dwdy[idx] - dvdz[idx];
                        const Real omega_y = dudz[idx] - dwdx[idx];
                        const Real omega_z = dvdx[idx] - dudy[idx];
                        
                        Omega[idx] = omega_x*omega_x + omega_y*omega_y + omega_z*omega_z;
                    }
                }
            }
            
            int count_valid_sensor = 0;
            
            const Real half = Real(1)/Real(2);
            
            // Compute the shock sensor on the side in the x-direction.
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                    {
                        // Compute the linear indices.
                        const int idx_sensor = i +
                            j*(interior_dim_0 + 1) +
                            k*(interior_dim_0 + 1)*
                                interior_dim_1;
                        
                        const int idx_L = (i - 1 + 1) +
                            (j + 1)*(interior_dim_0 + 2) +
                            (k + 1)*(interior_dim_0 + 2)*
                                (interior_dim_1 + 2);
                        
                        const int idx_R = (i + 0 + 1) +
                            (j + 1)*(interior_dim_0 + 2) +
                            (k + 1)*(interior_dim_0 + 2)*
                                (interior_dim_1 + 2);
                        
                        const Real Ducros_value_L = (-theta[idx_L]*std::abs(theta[idx_L]))/
                            (theta[idx_L]*theta[idx_L] + Omega[idx_L] + EPSILON);
                        
                        const Real Ducros_value_R = (-theta[idx_R]*std::abs(theta[idx_R]))/
                            (theta[idx_R]*theta[idx_R] + Omega[idx_R] + EPSILON);
                        
                        const Real Ducros_value_midpoint = half*(Ducros_value_L + Ducros_value_R);
                        
                        s_x[idx_sensor] = Real(0);
                        if (Ducros_value_midpoint > d_threshold_sensor_shock && use_shock_capturing)
                        {
                            s_x[idx_sensor] = Real(1);
                        }
                        
                        const Real rho_s_x_midpoint = half*(rho_s_x[idx_L] + rho_s_x[idx_R]);
                        if (rho_s_x_midpoint > d_threshold_sensor_interface && use_interface_capturing)
                        {
                            s_x[idx_sensor] = Real(1);
                        }
                        
                        if (s_x[idx_sensor] > Real(0))
                        {
                            count_valid_sensor++;
                        }
                    }
                }
            }
            
            // Compute the shock sensor on the side in the y-direction.
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_sensor = i +
                            j*interior_dim_0 +
                            k*interior_dim_0*
                                (interior_dim_1 + 1);
                        
                        const int idx_B = (i + 1) +
                            (j - 1 + 1)*(interior_dim_0 + 2) +
                            (k + 1)*(interior_dim_0 + 2)*
                                (interior_dim_1 + 2);
                        
                        const int idx_T = (i + 1) +
                            (j + 0 + 1)*(interior_dim_0 + 2) +
                            (k + 1)*(interior_dim_0 + 2)*
                                (interior_dim_1 + 2);
                        
                        const Real Ducros_value_B = (-theta[idx_B]*std::abs(theta[idx_B]))/
                            (theta[idx_B]*theta[idx_B] + Omega[idx_B] + EPSILON);
                        
                        const Real Ducros_value_T = (-theta[idx_T]*std::abs(theta[idx_T]))/
                            (theta[idx_T]*theta[idx_T] + Omega[idx_T] + EPSILON);
                        
                        const Real Ducros_value_midpoint = half*(Ducros_value_B + Ducros_value_T);
                        
                        s_y[idx_sensor] = Real(0);
                        if (Ducros_value_midpoint > d_threshold_sensor_shock && use_shock_capturing)
                        {
                            s_y[idx_sensor] = Real(1);
                        }
                        
                        const Real rho_s_y_midpoint = half*(rho_s_y[idx_B] + rho_s_y[idx_T]);
                        if (rho_s_y_midpoint > d_threshold_sensor_interface && use_interface_capturing)
                        {
                            s_y[idx_sensor] = Real(1);
                        }
                        
                        if (s_y[idx_sensor] > Real(0))
                        {
                            count_valid_sensor++;
                        }
                    }
                }
            }
            
            // Compute the shock sensor on the side in the z-direction.
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2 + 1; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_sensor = i +
                            j*interior_dim_0 +
                            k*interior_dim_0*
                                interior_dim_1;
                        
                        const int idx_B = (i + 1) +
                            (j + 1)*(interior_dim_0 + 2) +
                            (k - 1 + 1)*(interior_dim_0 + 2)*
                                (interior_dim_1 + 2);
                        
                        const int idx_T = (i + 1) +
                            (j + 1)*(interior_dim_0 + 2) +
                            (k + 0 + 1)*(interior_dim_0 + 2)*
                                (interior_dim_1 + 2);
                        
                        const Real Ducros_value_B = (-theta[idx_B]*std::abs(theta[idx_B]))/
                            (theta[idx_B]*theta[idx_B] + Omega[idx_B] + EPSILON);
                        
                        const Real Ducros_value_T = (-theta[idx_T]*std::abs(theta[idx_T]))/
                            (theta[idx_T]*theta[idx_T] + Omega[idx_T] + EPSILON);
                        
                        const Real Ducros_value_midpoint = half*(Ducros_value_B + Ducros_value_T);
                        
                        s_z[idx_sensor] = Real(0);
                        if (Ducros_value_midpoint > d_threshold_sensor_shock && use_shock_capturing)
                        {
                            s_z[idx_sensor] = Real(1);
                        }
                        
                        const Real rho_s_z_midpoint = half*(rho_s_z[idx_B] + rho_s_z[idx_T]);
                        if (rho_s_z_midpoint > d_threshold_sensor_interface && use_interface_capturing)
                        {
                            s_z[idx_sensor] = Real(1);
                        }
                        
                        if (s_z[idx_sensor] > Real(0))
                        {
                            count_valid_sensor++;
                        }
                    }
                }
            }
            
            if (d_has_advective_eqn_form)
            {
                // Compute the cell-based Ducros-like shock sensor.
                for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                {
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_sensor = i +
                                j*interior_dim_0 +
                                k*interior_dim_0*
                                    interior_dim_1;
                            
                            const int idx = (i + 1) +
                                (j + 1)*(interior_dim_0 + 2) +
                                (k + 1)*(interior_dim_0 + 2)*
                                    (interior_dim_1 + 2);
                            
                            const Real Ducros_value = (-theta[idx]*std::abs(theta[idx]))/
                                (theta[idx]*theta[idx] + Omega[idx] + EPSILON);
                            
                            s[idx_sensor] = Real(0);
                            if (Ducros_value > d_threshold_sensor_shock && use_shock_capturing)
                            {
                                s[idx_sensor] = Real(1);
                            }
                            
                            if (rho_s_x[idx] > d_threshold_sensor_interface && use_interface_capturing)
                            {
                                s[idx_sensor] = Real(1);
                            }
                            if (rho_s_y[idx] > d_threshold_sensor_interface && use_interface_capturing)
                            {
                                s[idx_sensor] = Real(1);
                            }
                            if (rho_s_z[idx] > d_threshold_sensor_interface && use_interface_capturing)
                            {
                                s[idx_sensor] = Real(1);
                            }
                            
                            if (s[idx_sensor] > Real(0))
                            {
                                count_valid_sensor++;
                            }
                        }
                    }
                }
            }
            
            if (count_valid_sensor > 0)
            {
                perform_WCNS = true;
            }
        } // if (d_dim == tbox::Dimension(3))
    } // if (d_dim == tbox::Dimension(1))
    
    if (perform_WCNS)
    {
        /*
         * Declare temporary data containers for WENO interpolation.
         */
        
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > projection_variables;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > > characteristic_variables;
        
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > characteristic_variables_minus;
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > characteristic_variables_plus;
        
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > primitive_variables_minus;
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > primitive_variables_plus;
        
        HAMERS_SHARED_PTR<pdat::SideData<int> > bounded_flag_minus;
        HAMERS_SHARED_PTR<pdat::SideData<int> > bounded_flag_plus;
        
        /*
         * Initialize temporary data containers for WENO interpolation.
         */
        
        const int num_projection_var = basic_utilities->getNumberOfProjectionVariablesForPrimitiveVariables();
        projection_variables.reserve(num_projection_var);
        
        for (int vi = 0; vi < num_projection_var; vi++)
        {
            projection_variables.push_back(HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
                interior_box, 1, hier::IntVector::getOne(d_dim)));
        }
        
        characteristic_variables.resize(6);
        
        for (int m = 0; m < 6; m++)
        {
            characteristic_variables[m].reserve(d_num_eqn);
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                characteristic_variables[m].push_back(HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
                    interior_box, 1, hier::IntVector::getOne(d_dim)));
            }
        }
        
        characteristic_variables_minus.reserve(d_num_eqn);
        characteristic_variables_plus.reserve(d_num_eqn);
        primitive_variables_minus.reserve(d_num_eqn);
        primitive_variables_plus.reserve(d_num_eqn);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            characteristic_variables_minus.push_back(HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
                interior_box, 1, hier::IntVector::getOne(d_dim)));
            
            characteristic_variables_plus.push_back(HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
                interior_box, 1, hier::IntVector::getOne(d_dim)));
            
            primitive_variables_minus.push_back(HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
                interior_box, 1, hier::IntVector::getOne(d_dim)));
            
            primitive_variables_plus.push_back(HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
                interior_box, 1, hier::IntVector::getOne(d_dim)));
        }
        
        bounded_flag_minus.reset(
            new pdat::SideData<int>(interior_box, 1, hier::IntVector::getOne(d_dim)));
        
        bounded_flag_plus.reset(
            new pdat::SideData<int>(interior_box, 1, hier::IntVector::getOne(d_dim)));
        
        /*
         * Compute the side data of the projection variables for transformation between primitive variables and
         * characteristic variables.
         */
        
        basic_utilities->computeSideDataOfProjectionVariablesForPrimitiveVariables(
            projection_variables,
            domains);
        
        /*
         * Transform primitive variables to characteristic variables.
         */
        
        for (int m = 0; m < 6; m++)
        {
            basic_utilities->computeSideDataOfCharacteristicVariablesFromPrimitiveVariables(
                characteristic_variables[m],
                primitive_variables,
                projection_variables,
                m - 3,
                domains);
        }
        
        /*
         * Peform WENO interpolation.
         */
        
        performWENOInterpolation(
            characteristic_variables_minus,
            characteristic_variables_plus,
            characteristic_variables,
            domains);
        
        /*
         * Transform characteristic variables back to primitive variables.
         */
        
        basic_utilities->computeSideDataOfPrimitiveVariablesFromCharacteristicVariables(
            primitive_variables_minus,
            characteristic_variables_minus,
            projection_variables,
            domains);
        
        basic_utilities->computeSideDataOfPrimitiveVariablesFromCharacteristicVariables(
            primitive_variables_plus,
            characteristic_variables_plus,
            projection_variables,
            domains);
        
        /*
         * Check whether the interpolated side primitive variables are within the bounds.
         */
        
        basic_utilities->checkSideDataOfPrimitiveVariablesBounded(
            bounded_flag_minus,
            primitive_variables_minus,
            domains);
        
        basic_utilities->checkSideDataOfPrimitiveVariablesBounded(
            bounded_flag_plus,
            primitive_variables_plus,
            domains);
        
        /*
         * Declare containers to store pointers for computing mid-point fluxes.
         */
        
        std::vector<Real*> V_minus;
        std::vector<Real*> V_plus;
        V_minus.resize(d_num_eqn);
        V_plus.resize(d_num_eqn);
        
        int* flag_minus = nullptr;
        int* flag_plus = nullptr;
        
        /*
         * Coefficients for finite differencing.
         */
        
        const Real a_midpoint_r = Real(23)/Real(15);
        const Real b_midpoint_r = Real(1)/Real(30);
        
        const Real a_node_r = -Real(3)/Real(10);
        
        const Real a_midpoint = Real(3)/Real(2);
        const Real b_midpoint = Real(1)/Real(30);
        
        const Real a_node =  -Real(3)/Real(10);
        
        /*
         * Compute the convective flux and source using shock-capturing scheme.
         */
        if (d_dim == tbox::Dimension(1))
        {
            /*
             * Get the local lower index and the number of cells in each dimension.
             */
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_dim_0 = domain_dims[0];
            
            /*
             * Get the pointers to the velocity and convective flux cell data inside the flow model.
             * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
             */
            
            HAMERS_SHARED_PTR<pdat::CellData<Real> > velocity = d_flow_model->getCellData("VELOCITY");
            
            std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > convective_flux_node(1);
            convective_flux_node[0] = d_flow_model->getCellData("CONVECTIVE_FLUX_X");
            
            hier::IntVector num_subghosts_velocity = velocity->getGhostCellWidth();
            hier::IntVector subghostcell_dims_velocity = velocity->getGhostBox().numberCells();
            
            hier::IntVector num_subghosts_convective_flux_x = convective_flux_node[0]->getGhostCellWidth();
            hier::IntVector subghostcell_dims_convective_flux_x = convective_flux_node[0]->getGhostBox().numberCells();
            
            const int num_subghosts_0_convective_flux_x = num_subghosts_convective_flux_x[0];
            const int subghostcell_dim_0_convective_flux_x = subghostcell_dims_convective_flux_x[0];
            
            const int num_subghosts_0_velocity = num_subghosts_velocity[0];
            const int subghostcell_dim_0_velocity = subghostcell_dims_velocity[0];
            
            Real* u = velocity->getPointer(0);
            
            std::vector<Real*> F_node_x;
            F_node_x.reserve(d_num_eqn);
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_node_x.push_back(convective_flux_node[0]->getPointer(ei));
            }
            
            std::vector<Real*> F_midpoint_x;
            F_midpoint_x.reserve(d_num_eqn);
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_midpoint_x.push_back(convective_flux_midpoint->getPointer(0, ei));
            }
            
            /*
             * Use first order interpolation if interpolated side primitive variables in x-direction
             * are out of bounds.
             */
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                V_minus[ei] = primitive_variables_minus[ei]->getPointer(0);
                V_plus[ei] = primitive_variables_plus[ei]->getPointer(0);
            }
            
            flag_minus = bounded_flag_minus->getPointer(0);
            flag_plus = bounded_flag_plus->getPointer(0);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                const int num_subghosts_0_primitive_var = num_subghosts_primitive_var[ei][0];
                const int num_subghosts_1_primitive_var = num_subghosts_primitive_var[ei][1];
                const int subghostcell_dim_0_primitive_var = subghostcell_dims_primitive_var[ei][0];
                
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0 - 1; i < domain_lo_0 + domain_dim_0 + 2; i++)
                {
                    // Compute the linear indices.
                    const int idx_midpoint_x = i + 1;
                    
                    const int idx_cell_L = i - 1 + num_subghosts_0_primitive_var;
                    const int idx_cell_R = i     + num_subghosts_0_primitive_var;
                    
                    if (flag_minus[idx_midpoint_x] == 0 || flag_plus[idx_midpoint_x] == 0)
                    {
                        V_minus[ei][idx_midpoint_x] = V[ei][idx_cell_L];
                        V_plus[ei][idx_midpoint_x]  = V[ei][idx_cell_R];
                    }
                }
            }
            
            /*
             * Compute mid-point flux in the x-direction.
             */
            
            if (d_has_advective_eqn_form)
            {
                riemann_solver->computeConvectiveFluxAndVelocityFromPrimitiveVariables(
                    convective_flux_midpoint,
                    velocity_midpoint,
                    primitive_variables_minus,
                    primitive_variables_plus,
                    DIRECTION::X_DIRECTION,
                    RIEMANN_SOLVER::HLLC,
                    domains[0]);
            }
            else
            {
                riemann_solver->computeConvectiveFluxFromPrimitiveVariables(
                    convective_flux_midpoint,
                    primitive_variables_minus,
                    primitive_variables_plus,
                    DIRECTION::X_DIRECTION,
                    RIEMANN_SOLVER::HLLC,
                    domains[0]);
            }
            
            /*
             * Reconstruct the flux in the x-direction.
             */
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                Real* F_face_x = convective_flux->getPointer(0, ei);
                
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                {
                    // Compute the linear indices.
                    const int idx_face_x = i;
                    
                    if (s_x[idx_face_x] > Real(0))
                    {
                        const int idx_midpoint_x   = i + 1;
                        const int idx_midpoint_x_L = i;
                        const int idx_midpoint_x_R = i + 2;
                        
                        const int idx_node_L = i - 1 + num_subghosts_0_convective_flux_x;
                        const int idx_node_R = i     + num_subghosts_0_convective_flux_x;
                        
                        F_face_x[idx_face_x] = Real(dt)*(
                            a_midpoint_r*F_midpoint_x[ei][idx_midpoint_x] +
                            b_midpoint_r*(F_midpoint_x[ei][idx_midpoint_x_L] + F_midpoint_x[ei][idx_midpoint_x_R]) +
                            a_node_r*(F_node_x[ei][idx_node_L] + F_node_x[ei][idx_node_R])
                            );
                    }
                }
            }
            
            /*
             * Compute the source.
             */
            
            if (d_has_advective_eqn_form)
            {
                Real* u_midpoint_x = velocity_midpoint->getPointer(0, 0);
                
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    if (d_eqn_form[ei] == EQN_FORM::ADVECTIVE)
                    {
                        Real* S = source_scratch->getPointer(ei);
                        
                        const int num_subghosts_0_conservative_var = num_subghosts_conservative_var[ei][0];
                        
                        HAMERS_PRAGMA_SIMD
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_cell_nghost = i;
                            
                            if (s[idx_cell_nghost] > Real(0))
                            {
                                const int idx_cell_wghost = i + num_subghosts_0_conservative_var;
                                
                                const int idx_cell_wghost_x_L = i - 1 + num_subghosts_0_velocity;
                                const int idx_cell_wghost_x_R = i + 1 + num_subghosts_0_velocity;
                                
                                const int idx_midpoint_x_LL = i;
                                const int idx_midpoint_x_L  = i + 1;
                                const int idx_midpoint_x_R  = i + 2;
                                const int idx_midpoint_x_RR = i + 3;
                                
                                S[idx_cell_nghost] = Real(dt)*Q[ei][idx_cell_wghost]*((
                                    a_midpoint*(u_midpoint_x[idx_midpoint_x_R]  - u_midpoint_x[idx_midpoint_x_L]) +
                                    b_midpoint*(u_midpoint_x[idx_midpoint_x_RR] - u_midpoint_x[idx_midpoint_x_LL]) +
                                    a_node*(u[idx_cell_wghost_x_R] - u[idx_cell_wghost_x_L]))/Real(dx[0]));
                            }
                        }
                    }
                }
            }
        }
        else if (d_dim == tbox::Dimension(2))
        {
            /*
             * Get the local lower indices and the number of cells in each dimension.
             */
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_lo_1 = domain_lo[1];
            const int domain_dim_0 = domain_dims[0];
            const int domain_dim_1 = domain_dims[1];
            
            /*
             * Get the interior dimension.
             */
            
            const int interior_dim_0 = interior_dims[0];
            
            /*
             * Get the pointers to the velocity and convective flux cell data inside the flow model.
             * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
             */
            
            HAMERS_SHARED_PTR<pdat::CellData<Real> > velocity = d_flow_model->getCellData("VELOCITY");
            
            std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > convective_flux_node(2);
            convective_flux_node[0] = d_flow_model->getCellData("CONVECTIVE_FLUX_X");
            convective_flux_node[1] = d_flow_model->getCellData("CONVECTIVE_FLUX_Y");
            
            hier::IntVector num_subghosts_velocity = velocity->getGhostCellWidth();
            hier::IntVector subghostcell_dims_velocity = velocity->getGhostBox().numberCells();
            
            hier::IntVector num_subghosts_convective_flux_x = convective_flux_node[0]->getGhostCellWidth();
            hier::IntVector subghostcell_dims_convective_flux_x = convective_flux_node[0]->getGhostBox().numberCells();
            
            hier::IntVector num_subghosts_convective_flux_y = convective_flux_node[1]->getGhostCellWidth();
            hier::IntVector subghostcell_dims_convective_flux_y = convective_flux_node[1]->getGhostBox().numberCells();
            
            const int num_subghosts_0_velocity = num_subghosts_velocity[0];
            const int num_subghosts_1_velocity = num_subghosts_velocity[1];
            const int subghostcell_dim_0_velocity = subghostcell_dims_velocity[0];
            
            const int num_subghosts_0_convective_flux_x = num_subghosts_convective_flux_x[0];
            const int num_subghosts_1_convective_flux_x = num_subghosts_convective_flux_x[1];
            const int subghostcell_dim_0_convective_flux_x = subghostcell_dims_convective_flux_x[0];
            
            const int num_subghosts_0_convective_flux_y = num_subghosts_convective_flux_y[0];
            const int num_subghosts_1_convective_flux_y = num_subghosts_convective_flux_y[1];
            const int subghostcell_dim_0_convective_flux_y = subghostcell_dims_convective_flux_y[0];
            
            Real* u = velocity->getPointer(0);
            Real* v = velocity->getPointer(1);
            
            std::vector<Real*> F_node_x;
            std::vector<Real*> F_node_y;
            F_node_x.reserve(d_num_eqn);
            F_node_y.reserve(d_num_eqn);
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_node_x.push_back(convective_flux_node[0]->getPointer(ei));
                F_node_y.push_back(convective_flux_node[1]->getPointer(ei));
            }
            
            std::vector<Real*> F_midpoint_x;
            std::vector<Real*> F_midpoint_y;
            F_midpoint_x.reserve(d_num_eqn);
            F_midpoint_y.reserve(d_num_eqn);
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_midpoint_x.push_back(convective_flux_midpoint->getPointer(0, ei));
                F_midpoint_y.push_back(convective_flux_midpoint->getPointer(1, ei));
            }
            
            /*
             * Use first order interpolation if interpolated side primitive variables in x-direction
             * are out of bounds.
             */
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                V_minus[ei] = primitive_variables_minus[ei]->getPointer(0);
                V_plus[ei] = primitive_variables_plus[ei]->getPointer(0);
            }
            
            flag_minus = bounded_flag_minus->getPointer(0);
            flag_plus = bounded_flag_plus->getPointer(0);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                const int num_subghosts_0_primitive_var = num_subghosts_primitive_var[ei][0];
                const int num_subghosts_1_primitive_var = num_subghosts_primitive_var[ei][1];
                const int subghostcell_dim_0_primitive_var = subghostcell_dims_primitive_var[ei][0];
                
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0 - 1; i < domain_lo_0 + domain_dim_0 + 2; i++)
                    {
                        // Compute the linear indices.
                        const int idx_midpoint_x = (i + 1) +
                            (j + 1)*(interior_dim_0 + 3);
                        
                        const int idx_cell_L = (i - 1 + num_subghosts_0_primitive_var) +
                            (j + num_subghosts_1_primitive_var)*subghostcell_dim_0_primitive_var;
                        
                        const int idx_cell_R = (i + num_subghosts_0_primitive_var) +
                            (j + num_subghosts_1_primitive_var)*subghostcell_dim_0_primitive_var;
                        
                        if (flag_minus[idx_midpoint_x] == 0 || flag_plus[idx_midpoint_x] == 0)
                        {
                            V_minus[ei][idx_midpoint_x] = V[ei][idx_cell_L];
                            V_plus[ei][idx_midpoint_x] = V[ei][idx_cell_R];
                        }
                    }
                }
            }
            
            /*
             * Use first order interpolation if interpolated side primitive variables in y-direction
             * are out of bounds.
             */
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                V_minus[ei] = primitive_variables_minus[ei]->getPointer(1);
                V_plus[ei] = primitive_variables_plus[ei]->getPointer(1);
            }
            
            flag_minus = bounded_flag_minus->getPointer(1);
            flag_plus = bounded_flag_plus->getPointer(1);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                const int num_subghosts_0_primitive_var = num_subghosts_primitive_var[ei][0];
                const int num_subghosts_1_primitive_var = num_subghosts_primitive_var[ei][1];
                const int subghostcell_dim_0_primitive_var = subghostcell_dims_primitive_var[ei][0];
                
                for (int j = domain_lo_1 - 1; j < domain_lo_1 + domain_dim_1 + 2; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_midpoint_y = (i + 1) +
                            (j + 1)*(interior_dim_0 + 2);
                        
                        const int idx_cell_B = (i + num_subghosts_0_primitive_var) +
                            (j - 1 + num_subghosts_1_primitive_var)*subghostcell_dim_0_primitive_var;
                        
                        const int idx_cell_T = (i + num_subghosts_0_primitive_var) +
                            (j + num_subghosts_1_primitive_var)*subghostcell_dim_0_primitive_var;
                        
                        if (flag_minus[idx_midpoint_y] == 0 || flag_plus[idx_midpoint_y] == 0)
                        {
                            V_minus[ei][idx_midpoint_y] = V[ei][idx_cell_B];
                            V_plus[ei][idx_midpoint_y] = V[ei][idx_cell_T];
                        }
                    }
                }
            }
            
            /*
             * Compute mid-point flux in the x-direction.
             */
            
            if (d_has_advective_eqn_form)
            {
                riemann_solver->computeConvectiveFluxAndVelocityFromPrimitiveVariables(
                    convective_flux_midpoint,
                    velocity_midpoint,
                    primitive_variables_minus,
                    primitive_variables_plus,
                    DIRECTION::X_DIRECTION,
                    RIEMANN_SOLVER::HLLC,
                    domains[0]);
            }
            else
            {
                riemann_solver->computeConvectiveFluxFromPrimitiveVariables(
                    convective_flux_midpoint,
                    primitive_variables_minus,
                    primitive_variables_plus,
                    DIRECTION::X_DIRECTION,
                    RIEMANN_SOLVER::HLLC,
                    domains[0]);
            }
            
            /*
             * Compute mid-point flux in the y-direction.
             */
            
            if (d_has_advective_eqn_form)
            {
                riemann_solver->computeConvectiveFluxAndVelocityFromPrimitiveVariables(
                    convective_flux_midpoint,
                    velocity_midpoint,
                    primitive_variables_minus,
                    primitive_variables_plus,
                    DIRECTION::Y_DIRECTION,
                    RIEMANN_SOLVER::HLLC,
                    domains[1]);
            }
            else
            {
                riemann_solver->computeConvectiveFluxFromPrimitiveVariables(
                    convective_flux_midpoint,
                    primitive_variables_minus,
                    primitive_variables_plus,
                    DIRECTION::Y_DIRECTION,
                    RIEMANN_SOLVER::HLLC,
                    domains[1]);
            }
            
            /*
             * Reconstruct the flux in the x-direction.
             */
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                Real* F_face_x = convective_flux->getPointer(0, ei);
                
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                    {
                        // Compute the linear indices.
                        const int idx_face_x = i +
                            j*(interior_dim_0 + 1);
                        
                        if (s_x[idx_face_x] > Real(0))
                        {
                            const int idx_midpoint_x = (i + 1) +
                                (j + 1)*(interior_dim_0 + 3);
                            
                            const int idx_midpoint_x_L = i +
                                (j + 1)*(interior_dim_0 + 3);
                            
                            const int idx_midpoint_x_R = (i + 2) +
                                (j + 1)*(interior_dim_0 + 3);
                            
                            const int idx_node_L = (i - 1 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                            
                            const int idx_node_R = (i + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                            
                            F_face_x[idx_face_x] = Real(dt)*(
                                a_midpoint_r*F_midpoint_x[ei][idx_midpoint_x] +
                                b_midpoint_r*(F_midpoint_x[ei][idx_midpoint_x_L] + F_midpoint_x[ei][idx_midpoint_x_R]) +
                                a_node_r*(F_node_x[ei][idx_node_L] + F_node_x[ei][idx_node_R])
                                );
                        }
                    }
                }
            }
            
            /*
             * Reconstruct the flux in the y-direction.
             */
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                Real* F_face_y = convective_flux->getPointer(1, ei);
                
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_face_y = i +
                            j*interior_dim_0;
                        
                        if (s_y[idx_face_y] > Real(0))
                        {
                            const int idx_midpoint_y = (i + 1) +
                                (j + 1)*(interior_dim_0 + 2);
                            
                            const int idx_midpoint_y_B = (i + 1) +
                                j*(interior_dim_0 + 2);
                            
                            const int idx_midpoint_y_T = (i + 1) +
                                (j + 2)*(interior_dim_0 + 2);
                            
                            const int idx_node_B = (i + num_subghosts_0_convective_flux_y) +
                                (j - 1 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                            
                            const int idx_node_T = (i + num_subghosts_0_convective_flux_y) +
                                (j + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                            
                            
                            F_face_y[idx_face_y] = Real(dt)*(
                                a_midpoint_r*F_midpoint_y[ei][idx_midpoint_y] +
                                b_midpoint_r*(F_midpoint_y[ei][idx_midpoint_y_B] + F_midpoint_y[ei][idx_midpoint_y_T]) +
                                a_node_r*(F_node_y[ei][idx_node_B] + F_node_x[ei][idx_node_T])
                                );
                        }
                    }
                }
            }
            
            /*
             * Compute the source.
             */
            
            if (d_has_advective_eqn_form)
            {
                Real* u_midpoint_x = velocity_midpoint->getPointer(0, 0);
                Real* v_midpoint_y = velocity_midpoint->getPointer(1, 1);
                
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    if (d_eqn_form[ei] == EQN_FORM::ADVECTIVE)
                    {
                        Real* S = source_scratch->getPointer(ei);
                        
                        const int num_subghosts_0_conservative_var = num_subghosts_conservative_var[ei][0];
                        const int num_subghosts_1_conservative_var = num_subghosts_conservative_var[ei][1];
                        const int subghostcell_dim_0_conservative_var = subghostcell_dims_conservative_var[ei][0];
                        
                        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                        {
                            HAMERS_PRAGMA_SIMD
                            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                            {
                                // Compute the linear indices.
                                const int idx_cell_nghost = i + j*interior_dim_0;
                                
                                if (s[idx_cell_nghost] > Real(0))
                                {
                                    const int idx_cell_wghost = (i + num_subghosts_0_conservative_var) +
                                        (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                                    
                                    const int idx_cell_wghost_x_L = (i - 1 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                    
                                    const int idx_cell_wghost_x_R = (i + 1 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                    
                                    const int idx_cell_wghost_y_B = (i + num_subghosts_0_velocity) +
                                        (j - 1 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                    
                                    const int idx_cell_wghost_y_T = (i + num_subghosts_0_velocity) +
                                        (j + 1 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                    
                                    const int idx_midpoint_x_LL = i + 
                                        (j + 1)*(interior_dim_0 + 3);
                                    
                                    const int idx_midpoint_x_L = (i + 1) +
                                        (j + 1)*(interior_dim_0 + 3);
                                    
                                    const int idx_midpoint_x_R = (i + 2) +
                                        (j + 1)*(interior_dim_0 + 3);
                                    
                                    const int idx_midpoint_x_RR = (i + 3) +
                                        (j + 1)*(interior_dim_0 + 3);
                                    
                                    const int idx_midpoint_y_BB = (i + 1) +
                                        j*(interior_dim_0 + 2);
                                    
                                    const int idx_midpoint_y_B = (i + 1) +
                                        (j + 1)*(interior_dim_0 + 2);
                                    
                                    const int idx_midpoint_y_T = (i + 1) +
                                        (j + 2)*(interior_dim_0 + 2);
                                    
                                    const int idx_midpoint_y_TT = (i + 1) +
                                        (j + 3)*(interior_dim_0 + 2);
                                    
                                    S[idx_cell_nghost] = Real(dt)*Q[ei][idx_cell_wghost]*((
                                        a_midpoint*(u_midpoint_x[idx_midpoint_x_R]  - u_midpoint_x[idx_midpoint_x_L]) +
                                        b_midpoint*(u_midpoint_x[idx_midpoint_x_RR] - u_midpoint_x[idx_midpoint_x_LL]) +
                                        a_node*(u[idx_cell_wghost_x_R] - u[idx_cell_wghost_x_L])
                                        )/Real(dx[0]) + (
                                        a_midpoint*(v_midpoint_y[idx_midpoint_y_T]  - v_midpoint_y[idx_midpoint_y_B]) +
                                        b_midpoint*(v_midpoint_y[idx_midpoint_y_TT] - v_midpoint_y[idx_midpoint_y_BB]) +
                                        a_node*(v[idx_cell_wghost_y_T] - v[idx_cell_wghost_y_B])
                                        )/Real(dx[1]));
                                }
                            }
                        }
                    }
                }
            }
        }
        else if (d_dim == tbox::Dimension(3))
        {
            /*
             * Get the local lower indices and the number of cells in each dimension.
             */
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_lo_1 = domain_lo[1];
            const int domain_lo_2 = domain_lo[2];
            const int domain_dim_0 = domain_dims[0];
            const int domain_dim_1 = domain_dims[1];
            const int domain_dim_2 = domain_dims[2];
            
            /*
             * Get the interior dimensions.
             */
            
            const int interior_dim_0 = interior_dims[0];
            const int interior_dim_1 = interior_dims[1];
            
            /*
             * Get the pointers to the velocity and convective flux cell data inside the flow model.
             * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
             */
            
            HAMERS_SHARED_PTR<pdat::CellData<Real> > velocity = d_flow_model->getCellData("VELOCITY");
            
            std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > convective_flux_node(3);
            convective_flux_node[0] = d_flow_model->getCellData("CONVECTIVE_FLUX_X");
            convective_flux_node[1] = d_flow_model->getCellData("CONVECTIVE_FLUX_Y");
            convective_flux_node[2] = d_flow_model->getCellData("CONVECTIVE_FLUX_Z");
            
            hier::IntVector num_subghosts_velocity = velocity->getGhostCellWidth();
            hier::IntVector subghostcell_dims_velocity = velocity->getGhostBox().numberCells();
            
            hier::IntVector num_subghosts_convective_flux_x = convective_flux_node[0]->getGhostCellWidth();
            hier::IntVector subghostcell_dims_convective_flux_x = convective_flux_node[0]->getGhostBox().numberCells();
            
            hier::IntVector num_subghosts_convective_flux_y = convective_flux_node[1]->getGhostCellWidth();
            hier::IntVector subghostcell_dims_convective_flux_y = convective_flux_node[1]->getGhostBox().numberCells();
            
            hier::IntVector num_subghosts_convective_flux_z = convective_flux_node[2]->getGhostCellWidth();
            hier::IntVector subghostcell_dims_convective_flux_z = convective_flux_node[2]->getGhostBox().numberCells();
            
            const int num_subghosts_0_velocity = num_subghosts_velocity[0];
            const int num_subghosts_1_velocity = num_subghosts_velocity[1];
            const int num_subghosts_2_velocity = num_subghosts_velocity[2];
            const int subghostcell_dim_0_velocity = subghostcell_dims_velocity[0];
            const int subghostcell_dim_1_velocity = subghostcell_dims_velocity[1];
            
            const int num_subghosts_0_convective_flux_x = num_subghosts_convective_flux_x[0];
            const int num_subghosts_1_convective_flux_x = num_subghosts_convective_flux_x[1];
            const int num_subghosts_2_convective_flux_x = num_subghosts_convective_flux_x[2];
            const int subghostcell_dim_0_convective_flux_x = subghostcell_dims_convective_flux_x[0];
            const int subghostcell_dim_1_convective_flux_x = subghostcell_dims_convective_flux_x[1];
            
            const int num_subghosts_0_convective_flux_y = num_subghosts_convective_flux_y[0];
            const int num_subghosts_1_convective_flux_y = num_subghosts_convective_flux_y[1];
            const int num_subghosts_2_convective_flux_y = num_subghosts_convective_flux_y[2];
            const int subghostcell_dim_0_convective_flux_y = subghostcell_dims_convective_flux_y[0];
            const int subghostcell_dim_1_convective_flux_y = subghostcell_dims_convective_flux_y[1];
            
            const int num_subghosts_0_convective_flux_z = num_subghosts_convective_flux_z[0];
            const int num_subghosts_1_convective_flux_z = num_subghosts_convective_flux_z[1];
            const int num_subghosts_2_convective_flux_z = num_subghosts_convective_flux_z[2];
            const int subghostcell_dim_0_convective_flux_z = subghostcell_dims_convective_flux_z[0];
            const int subghostcell_dim_1_convective_flux_z = subghostcell_dims_convective_flux_z[1];
            
            Real* u = velocity->getPointer(0);
            Real* v = velocity->getPointer(1);
            Real* w = velocity->getPointer(2);
            
            std::vector<Real*> F_node_x;
            std::vector<Real*> F_node_y;
            std::vector<Real*> F_node_z;
            F_node_x.reserve(d_num_eqn);
            F_node_y.reserve(d_num_eqn);
            F_node_z.reserve(d_num_eqn);
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_node_x.push_back(convective_flux_node[0]->getPointer(ei));
                F_node_y.push_back(convective_flux_node[1]->getPointer(ei));
                F_node_z.push_back(convective_flux_node[2]->getPointer(ei));
            }
            
            std::vector<Real*> F_midpoint_x;
            std::vector<Real*> F_midpoint_y;
            std::vector<Real*> F_midpoint_z;
            F_midpoint_x.reserve(d_num_eqn);
            F_midpoint_y.reserve(d_num_eqn);
            F_midpoint_z.reserve(d_num_eqn);
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_midpoint_x.push_back(convective_flux_midpoint->getPointer(0, ei));
                F_midpoint_y.push_back(convective_flux_midpoint->getPointer(1, ei));
                F_midpoint_z.push_back(convective_flux_midpoint->getPointer(2, ei));
            }
            
            /*
             * Use first order interpolation if interpolated side primitive variables in x-direction
             * are out of bounds.
             */
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                V_minus[ei] = primitive_variables_minus[ei]->getPointer(0);
                V_plus[ei] = primitive_variables_plus[ei]->getPointer(0);
            }
            
            flag_minus = bounded_flag_minus->getPointer(0);
            flag_plus = bounded_flag_plus->getPointer(0);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                const int num_subghosts_0_primitive_var = num_subghosts_primitive_var[ei][0];
                const int num_subghosts_1_primitive_var = num_subghosts_primitive_var[ei][1];
                const int num_subghosts_2_primitive_var = num_subghosts_primitive_var[ei][2];
                const int subghostcell_dim_0_primitive_var = subghostcell_dims_primitive_var[ei][0];
                const int subghostcell_dim_1_primitive_var = subghostcell_dims_primitive_var[ei][1];
                
                for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                {
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = domain_lo_0 - 1; i < domain_lo_0 + domain_dim_0 + 2; i++)
                        {
                            // Compute the linear indices.
                            const int idx_midpoint_x = (i + 1) +
                                (j + 1)*(interior_dim_0 + 3) +
                                (k + 1)*(interior_dim_0 + 3)*
                                    (interior_dim_1 + 2);
                            
                            const int idx_cell_L = (i - 1 + num_subghosts_0_primitive_var) +
                                (j + num_subghosts_1_primitive_var)*subghostcell_dim_0_primitive_var +
                                (k + num_subghosts_2_primitive_var)*subghostcell_dim_0_primitive_var*
                                    subghostcell_dim_1_primitive_var;
                            
                            const int idx_cell_R = (i + num_subghosts_0_primitive_var) +
                                (j + num_subghosts_1_primitive_var)*subghostcell_dim_0_primitive_var +
                                (k + num_subghosts_2_primitive_var)*subghostcell_dim_0_primitive_var*
                                    subghostcell_dim_1_primitive_var;
                            
                            if (flag_minus[idx_midpoint_x] == 0 || flag_plus[idx_midpoint_x] == 0)
                            {
                                V_minus[ei][idx_midpoint_x] = V[ei][idx_cell_L];
                                V_plus[ei][idx_midpoint_x] = V[ei][idx_cell_R];
                            }
                        }
                    }
                }
            }
            
            /*
             * Use first order interpolation if interpolated side primitive variables in y-direction
             * are out of bounds.
             */
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                V_minus[ei] = primitive_variables_minus[ei]->getPointer(1);
                V_plus[ei] = primitive_variables_plus[ei]->getPointer(1);
            }
            
            flag_minus = bounded_flag_minus->getPointer(1);
            flag_plus = bounded_flag_plus->getPointer(1);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                const int num_subghosts_0_primitive_var = num_subghosts_primitive_var[ei][0];
                const int num_subghosts_1_primitive_var = num_subghosts_primitive_var[ei][1];
                const int num_subghosts_2_primitive_var = num_subghosts_primitive_var[ei][2];
                const int subghostcell_dim_0_primitive_var = subghostcell_dims_primitive_var[ei][0];
                const int subghostcell_dim_1_primitive_var = subghostcell_dims_primitive_var[ei][1];
                
                for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                {
                    for (int j = domain_lo_1 - 1; j < domain_lo_1 + domain_dim_1 + 2; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_midpoint_y = (i + 1) +
                                (j + 1)*(interior_dim_0 + 2) +
                                (k + 1)*(interior_dim_0 + 2)*
                                    (interior_dim_1 + 3);
                            
                            const int idx_cell_B = (i + num_subghosts_0_primitive_var) +
                                (j - 1 + num_subghosts_1_primitive_var)*subghostcell_dim_0_primitive_var +
                                (k + num_subghosts_2_primitive_var)*subghostcell_dim_0_primitive_var*
                                    subghostcell_dim_1_primitive_var;
                            
                            const int idx_cell_T = (i + num_subghosts_0_primitive_var) +
                                (j + num_subghosts_1_primitive_var)*subghostcell_dim_0_primitive_var +
                                (k + num_subghosts_2_primitive_var)*subghostcell_dim_0_primitive_var*
                                    subghostcell_dim_1_primitive_var;
                            
                            if (flag_minus[idx_midpoint_y] == 0 || flag_plus[idx_midpoint_y] == 0)
                            {
                                V_minus[ei][idx_midpoint_y] = V[ei][idx_cell_B];
                                V_plus[ei][idx_midpoint_y] = V[ei][idx_cell_T];
                            }
                        }
                    }
                }
            }
            
            /*
             * Use first order interpolation if interpolated side primitive variables in z-direction
             * are out of bounds.
             */
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                V_minus[ei] = primitive_variables_minus[ei]->getPointer(2);
                V_plus[ei] = primitive_variables_plus[ei]->getPointer(2);
            }
            
            flag_minus = bounded_flag_minus->getPointer(2);
            flag_plus = bounded_flag_plus->getPointer(2);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                const int num_subghosts_0_primitive_var = num_subghosts_primitive_var[ei][0];
                const int num_subghosts_1_primitive_var = num_subghosts_primitive_var[ei][1];
                const int num_subghosts_2_primitive_var = num_subghosts_primitive_var[ei][2];
                const int subghostcell_dim_0_primitive_var = subghostcell_dims_primitive_var[ei][0];
                const int subghostcell_dim_1_primitive_var = subghostcell_dims_primitive_var[ei][1];
                
                for (int k = domain_lo_2 - 1; k < domain_lo_2 + domain_dim_2 + 2; k++)
                {
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_midpoint_z = (i + 1) +
                                (j + 1)*(interior_dim_0 + 2) +
                                (k + 1)*(interior_dim_0 + 2)*
                                    (interior_dim_1 + 2);
                            
                            const int idx_cell_B = (i + num_subghosts_0_primitive_var) +
                                (j + num_subghosts_1_primitive_var)*subghostcell_dim_0_primitive_var +
                                (k - 1 + num_subghosts_2_primitive_var)*subghostcell_dim_0_primitive_var*
                                    subghostcell_dim_1_primitive_var;
                            
                            const int idx_cell_F = (i + num_subghosts_0_primitive_var) +
                                (j + num_subghosts_1_primitive_var)*subghostcell_dim_0_primitive_var +
                                (k + num_subghosts_2_primitive_var)*subghostcell_dim_0_primitive_var*
                                    subghostcell_dim_1_primitive_var;
                            
                            if (flag_minus[idx_midpoint_z] == 0 || flag_plus[idx_midpoint_z] == 0)
                            {
                                V_minus[ei][idx_midpoint_z] = V[ei][idx_cell_B];
                                V_plus[ei][idx_midpoint_z] = V[ei][idx_cell_F];
                            }
                        }
                    }
                }
            }
            
            /*
             * Compute mid-point flux in the x-direction.
             */
            
            if (d_has_advective_eqn_form)
            {
                riemann_solver->computeConvectiveFluxAndVelocityFromPrimitiveVariables(
                    convective_flux_midpoint,
                    velocity_midpoint,
                    primitive_variables_minus,
                    primitive_variables_plus,
                    DIRECTION::X_DIRECTION,
                    RIEMANN_SOLVER::HLLC);
            }
            else
            {
                riemann_solver->computeConvectiveFluxFromPrimitiveVariables(
                    convective_flux_midpoint,
                    primitive_variables_minus,
                    primitive_variables_plus,
                    DIRECTION::X_DIRECTION,
                    RIEMANN_SOLVER::HLLC);
            }
            
            /*
             * Compute mid-point flux in the y-direction.
             */
            
            if (d_has_advective_eqn_form)
            {
                riemann_solver->computeConvectiveFluxAndVelocityFromPrimitiveVariables(
                    convective_flux_midpoint,
                    velocity_midpoint,
                    primitive_variables_minus,
                    primitive_variables_plus,
                    DIRECTION::Y_DIRECTION,
                    RIEMANN_SOLVER::HLLC);
            }
            else
            {
                riemann_solver->computeConvectiveFluxFromPrimitiveVariables(
                    convective_flux_midpoint,
                    primitive_variables_minus,
                    primitive_variables_plus,
                    DIRECTION::Y_DIRECTION,
                    RIEMANN_SOLVER::HLLC);
            }
            
            /*
             * Compute mid-point flux in the z-direction.
             */
            
            if (d_has_advective_eqn_form)
            {
                riemann_solver->computeConvectiveFluxAndVelocityFromPrimitiveVariables(
                    convective_flux_midpoint,
                    velocity_midpoint,
                    primitive_variables_minus,
                    primitive_variables_plus,
                    DIRECTION::Z_DIRECTION,
                    RIEMANN_SOLVER::HLLC);
            }
            else
            {
                riemann_solver->computeConvectiveFluxFromPrimitiveVariables(
                    convective_flux_midpoint,
                    primitive_variables_minus,
                    primitive_variables_plus,
                    DIRECTION::Z_DIRECTION,
                    RIEMANN_SOLVER::HLLC);
            }
            
            /*
             * Reconstruct the flux in the x-direction.
             */
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                Real* F_face_x = convective_flux->getPointer(0, ei);
                
                for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                {
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                        {
                            // Compute the linear indices.
                            const int idx_face_x = i +
                                j*(interior_dim_0 + 1) +
                                k*(interior_dim_0 + 1)*
                                    interior_dim_1;
                            
                            if (s_x[idx_face_x] > Real(0))
                            {
                                const int idx_midpoint_x = (i + 1) +
                                    (j + 1)*(interior_dim_0 + 3) +
                                    (k + 1)*(interior_dim_0 + 3)*
                                        (interior_dim_1 + 2);
                                
                                const int idx_midpoint_x_L = i +
                                    (j + 1)*(interior_dim_0 + 3) +
                                    (k + 1)*(interior_dim_0 + 3)*
                                        (interior_dim_1 + 2);
                                
                                const int idx_midpoint_x_R = (i + 2) +
                                    (j + 1)*(interior_dim_0 + 3) +
                                    (k + 1)*(interior_dim_0 + 3)*
                                        (interior_dim_1 + 2);
                                
                                const int idx_node_L = (i - 1 + num_subghosts_0_convective_flux_x) +
                                    (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                    (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                        subghostcell_dim_1_convective_flux_x;
                                
                                const int idx_node_R = (i + num_subghosts_0_convective_flux_x) +
                                    (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                    (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                        subghostcell_dim_1_convective_flux_x;
                            
                                F_face_x[idx_face_x] = Real(dt)*(
                                    a_midpoint_r*F_midpoint_x[ei][idx_midpoint_x] +
                                    b_midpoint_r*(F_midpoint_x[ei][idx_midpoint_x_L] + F_midpoint_x[ei][idx_midpoint_x_R]) +
                                    a_node_r*(F_node_x[ei][idx_node_L] + F_node_x[ei][idx_node_R])
                                    );
                            }
                        }
                    }
                }
            }
            
            /*
             * Reconstruct the flux in the y-direction.
             */
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                Real* F_face_y = convective_flux->getPointer(1, ei);
                
                for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                {
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_face_y = i +
                                j*interior_dim_0 +
                                k*interior_dim_0*(interior_dim_1 + 1);
                            
                            if (s_y[idx_face_y] > Real(0))
                            {
                                const int idx_midpoint_y = (i + 1) +
                                    (j + 1)*(interior_dim_0 + 2) +
                                    (k + 1)*(interior_dim_0 + 2)*(interior_dim_1 + 3);
                                
                                const int idx_midpoint_y_B = (i + 1) +
                                    j*(interior_dim_0 + 2) +
                                    (k + 1)*(interior_dim_0 + 2)*(interior_dim_1 + 3);
                                
                                const int idx_midpoint_y_T = (i + 1) +
                                    (j + 2)*(interior_dim_0 + 2) +
                                    (k + 1)*(interior_dim_0 + 2)*(interior_dim_1 + 3);
                                
                                const int idx_node_B = (i + num_subghosts_0_convective_flux_y) +
                                    (j - 1 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                    (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                        subghostcell_dim_1_convective_flux_y;
                                
                                const int idx_node_T = (i + num_subghosts_0_convective_flux_y) +
                                    (j + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                    (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                        subghostcell_dim_1_convective_flux_y;
                                
                                F_face_y[idx_face_y] = Real(dt)*(
                                    a_midpoint_r*F_midpoint_y[ei][idx_midpoint_y] +
                                    b_midpoint_r*(F_midpoint_y[ei][idx_midpoint_y_B] + F_midpoint_y[ei][idx_midpoint_y_T]) +
                                    a_node_r*(F_node_y[ei][idx_node_B] + F_node_x[ei][idx_node_T])
                                    );
                                }
                        }
                    }
                }
            }
            
            /*
             * Reconstruct the flux in the z-direction.
             */
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                Real* F_face_z = convective_flux->getPointer(2, ei);
                
                for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2 + 1; k++)
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
                            
                            if (s_y[idx_face_z] > Real(0))
                            {
                                const int idx_midpoint_z = (i + 1) +
                                    (j + 1)*(interior_dim_0 + 2) +
                                    (k + 1)*(interior_dim_0 + 2)*(interior_dim_1 + 2);
                                
                                const int idx_midpoint_z_B = (i + 1) +
                                    (j + 1)*(interior_dim_0 + 2) +
                                    k*(interior_dim_0 + 2)*(interior_dim_1 + 2);
                                
                                const int idx_midpoint_z_F = (i + 1) +
                                    (j + 1)*(interior_dim_0 + 2) +
                                    (k + 2)*(interior_dim_0 + 2)*(interior_dim_1 + 2);
                                
                                const int idx_node_B = (i + num_subghosts_0_convective_flux_z) +
                                    (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                    (k - 1 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                        subghostcell_dim_1_convective_flux_z;
                                
                                const int idx_node_F = (i + num_subghosts_0_convective_flux_z) +
                                    (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                    (k + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                        subghostcell_dim_1_convective_flux_z;
                                
                                F_face_z[idx_face_z] = Real(dt)*(
                                    a_midpoint_r*F_midpoint_z[ei][idx_midpoint_z] +
                                    b_midpoint_r*(F_midpoint_z[ei][idx_midpoint_z_B] + F_midpoint_z[ei][idx_midpoint_z_F]) +
                                    a_node_r*(F_node_z[ei][idx_node_B] + F_node_z[ei][idx_node_F])
                                    );
                            }
                        }
                    }
                }
            }
            
            /*
             * Compute the source.
             */
            
            if (d_has_advective_eqn_form)
            {
                Real* u_midpoint_x = velocity_midpoint->getPointer(0, 0);
                Real* v_midpoint_y = velocity_midpoint->getPointer(1, 1);
                Real* w_midpoint_z = velocity_midpoint->getPointer(2, 2);
                
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    if (d_eqn_form[ei] == EQN_FORM::ADVECTIVE)
                    {
                        Real* S = source_scratch->getPointer(ei);
                        
                        const int num_subghosts_0_conservative_var = num_subghosts_conservative_var[ei][0];
                        const int num_subghosts_1_conservative_var = num_subghosts_conservative_var[ei][1];
                        const int num_subghosts_2_conservative_var = num_subghosts_conservative_var[ei][2];
                        const int subghostcell_dim_0_conservative_var = subghostcell_dims_conservative_var[ei][0];
                        const int subghostcell_dim_1_conservative_var = subghostcell_dims_conservative_var[ei][1];
                        
                        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                        {
                            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                            {
                                HAMERS_PRAGMA_SIMD
                                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                                {
                                    // Compute the linear indices.
                                    const int idx_cell_nghost = i +
                                        j*interior_dim_0 +
                                        k*interior_dim_0*
                                            interior_dim_1;
                                    
                                    if (s[idx_cell_nghost] > Real(0))
                                    {
                                        const int idx_cell_wghost = (i + num_subghosts_0_conservative_var) +
                                            (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                            (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                                subghostcell_dim_1_conservative_var;
                                        
                                        const int idx_cell_wghost_x_L = (i - 1 + num_subghosts_0_velocity) +
                                            (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                            (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                                subghostcell_dim_1_velocity;
                                        
                                        const int idx_cell_wghost_x_R = (i + 1 + num_subghosts_0_velocity) +
                                            (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                            (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                                subghostcell_dim_1_velocity;
                                        
                                        const int idx_cell_wghost_y_B = (i + num_subghosts_0_velocity) +
                                            (j - 1 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                            (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                                subghostcell_dim_1_velocity;
                                        
                                        const int idx_cell_wghost_y_T = (i + num_subghosts_0_velocity) +
                                            (j + 1 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                            (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                                subghostcell_dim_1_velocity;
                                        
                                        const int idx_cell_wghost_z_B = (i + num_subghosts_0_velocity) +
                                            (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                            (k - 1 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                                subghostcell_dim_1_velocity;
                                        
                                        const int idx_cell_wghost_z_F = (i + num_subghosts_0_velocity) +
                                            (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                            (k + 1 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                                subghostcell_dim_1_velocity;
                                        
                                        const int idx_midpoint_x_LL = i +
                                            (j + 1)*(interior_dim_0 + 3) +
                                            (k + 1)*(interior_dim_0 + 3)*
                                                (interior_dim_1 + 2);
                                        
                                        const int idx_midpoint_x_L = (i + 1) +
                                            (j + 1)*(interior_dim_0 + 3) +
                                            (k + 1)*(interior_dim_0 + 3)*
                                                (interior_dim_1 + 2);
                                        
                                        const int idx_midpoint_x_R = (i + 2) +
                                            (j + 1)*(interior_dim_0 + 3) +
                                            (k + 1)*(interior_dim_0 + 3)*
                                                (interior_dim_1 + 2);
                                        
                                        const int idx_midpoint_x_RR = (i + 3) +
                                            (j + 1)*(interior_dim_0 + 3) +
                                            (k + 1)*(interior_dim_0 + 3)*
                                                (interior_dim_1 + 2);
                                        
                                        const int idx_midpoint_y_BB = (i + 1) +
                                            j*(interior_dim_0 + 2) +
                                            (k + 1)*(interior_dim_0 + 2)*
                                                (interior_dim_1 + 3);
                                        
                                        const int idx_midpoint_y_B = (i + 1) +
                                            (j + 1)*(interior_dim_0 + 2) +
                                            (k + 1)*(interior_dim_0 + 2)*
                                                (interior_dim_1 + 3);
                                        
                                        const int idx_midpoint_y_T = (i + 1) +
                                            (j + 2)*(interior_dim_0 + 2) +
                                            (k + 1)*(interior_dim_0 + 2)*
                                                (interior_dim_1 + 3);
                                        
                                        const int idx_midpoint_y_TT = (i + 1) +
                                            (j + 3)*(interior_dim_0 + 2) +
                                            (k + 1)*(interior_dim_0 + 2)*
                                                (interior_dim_1 + 3);
                                        
                                        const int idx_midpoint_z_BB = (i + 1) +
                                            (j + 1)*(interior_dim_0 + 2) +
                                            k*(interior_dim_0 + 2)*
                                                (interior_dim_1 + 2);
                                        
                                        const int idx_midpoint_z_B = (i + 1) +
                                            (j + 1)*(interior_dim_0 + 2) +
                                            (k + 1)*(interior_dim_0 + 2)*
                                                (interior_dim_1 + 2);
                                        
                                        const int idx_midpoint_z_F = (i + 1) +
                                            (j + 1)*(interior_dim_0 + 2) +
                                            (k + 2)*(interior_dim_0 + 2)*
                                                (interior_dim_1 + 2);
                                        
                                        const int idx_midpoint_z_FF = (i + 1) +
                                            (j + 1)*(interior_dim_0 + 2) +
                                            (k + 3)*(interior_dim_0 + 2)*
                                                (interior_dim_1 + 2);
                                        
                                        S[idx_cell_nghost] = Real(dt)*Q[ei][idx_cell_wghost]*((
                                            a_midpoint*(u_midpoint_x[idx_midpoint_x_R]  - u_midpoint_x[idx_midpoint_x_L]) +
                                            b_midpoint*(u_midpoint_x[idx_midpoint_x_RR] - u_midpoint_x[idx_midpoint_x_LL]) +
                                            a_node*(u[idx_cell_wghost_x_R] - u[idx_cell_wghost_x_L])
                                            )/Real(dx[0]) + (
                                            a_midpoint*(v_midpoint_y[idx_midpoint_y_T]  - v_midpoint_y[idx_midpoint_y_B]) +
                                            b_midpoint*(v_midpoint_y[idx_midpoint_y_TT] - v_midpoint_y[idx_midpoint_y_BB]) +
                                            a_node*(v[idx_cell_wghost_y_T] - v[idx_cell_wghost_y_B])
                                            )/Real(dx[1]) + (
                                            a_midpoint*(w_midpoint_z[idx_midpoint_z_F]  - w_midpoint_z[idx_midpoint_z_B]) +
                                            b_midpoint*(w_midpoint_z[idx_midpoint_z_FF] - w_midpoint_z[idx_midpoint_z_BB]) +
                                            a_node*(w[idx_cell_wghost_z_F] - w[idx_cell_wghost_z_B])
                                            )/Real(dx[2]));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    /*
     * Unregister the patch and data of all registered derived cell variables in the flow model.
     */
    
    d_flow_model->unregisterPatch();
}


/*
 * Compute the convective flux and source due to splitting using shock-capturing scheme.
 */
void
ConvectiveFluxReconstructor::computeConvectiveFluxAndSourceOnPatchShockCapturing(
    hier::Patch& patch,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> >& convective_flux,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& source_scratch,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const hier::Box& domain,
    const double dt,
    const bool use_shock_capturing,
    const bool use_interface_capturing) const
{
    if (!use_shock_capturing && !use_interface_capturing)
    {
        return;
    }
    
    d_flow_model->setupRiemannSolver();
    d_flow_model->setupBasicUtilities();
    
    HAMERS_SHARED_PTR<FlowModelRiemannSolver> riemann_solver = d_flow_model->getFlowModelRiemannSolver();
    HAMERS_SHARED_PTR<FlowModelBasicUtilities> basic_utilities = d_flow_model->getFlowModelBasicUtilities();
    
    // Get the grid spacing.
    const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
        HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
    const double* const dx = patch_geom->getDx();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    domain_lo = domain.lower() - interior_box.lower();
    domain_dims = domain.numberCells();
    
    // Create domains in different directions.
    std::vector<hier::Box> domains(d_dim.getValue(), domain);
    
    // Allocate temporary patch data.
    HAMERS_SHARED_PTR<pdat::SideData<Real> > velocity_midpoint;
    
    if (d_has_advective_eqn_form)
    {
        velocity_midpoint.reset(new pdat::SideData<Real>(
            interior_box, d_dim.getValue(), hier::IntVector::getZero(d_dim))); // No ghost midpoints.
    }
    
    HAMERS_SHARED_PTR<pdat::SideData<Real> > convective_flux_midpoint(
        new pdat::SideData<Real>(interior_box, d_num_eqn, hier::IntVector::getZero(d_dim))); // No ghost midpoints.
    
    HAMERS_SHARED_PTR<pdat::SideData<Real> > discontinuity_sensor_side(
        new pdat::SideData<Real>(interior_box, 1, hier::IntVector::getZero(d_dim))); // No ghost midpoints.
    
    HAMERS_SHARED_PTR<pdat::CellData<Real> > discontinuity_sensor_cell(
        new pdat::CellData<Real>(interior_box, 1, hier::IntVector::getZero(d_dim)));
    
    HAMERS_SHARED_PTR<pdat::CellData<Real> > velocity_derivatives;
    HAMERS_SHARED_PTR<pdat::CellData<Real> > dilatation;
    HAMERS_SHARED_PTR<pdat::CellData<Real> > enstrophy;
    std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > density_sensors;
    
    bool perform_WCNS = false;
    
    if (d_dim > tbox::Dimension(1))
    {
        velocity_derivatives.reset(new pdat::CellData<Real>(
            interior_box, d_dim.getValue()*d_dim.getValue(), hier::IntVector::getOne(d_dim)));
        
        dilatation.reset(new pdat::CellData<Real>(
            interior_box, 1, hier::IntVector::getOne(d_dim)));
        
        enstrophy.reset(new pdat::CellData<Real>(
            interior_box, 1, hier::IntVector::getOne(d_dim)));
        
        density_sensors.resize(d_dim.getValue());
        for (int di = 0; di < d_dim.getValue(); di++)
        {
            density_sensors[di].reset(new pdat::CellData<Real>(
                interior_box, 1, hier::IntVector::getOne(d_dim)));
        }
    }
    
    /*
     * Register the patch and derived cell variables in the flow model and compute the corresponding cell data.
     */
    
    d_flow_model->registerPatchWithDataContext(patch, data_context);
    
    std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
    
    num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("DENSITY", d_num_conv_ghosts));
    num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("VELOCITY", d_num_conv_ghosts));
    num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("CONVECTIVE_FLUX_X", d_num_conv_ghosts));
    if (d_dim > tbox::Dimension(1))
    {
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("CONVECTIVE_FLUX_Y", d_num_conv_ghosts));
    }
    if (d_dim > tbox::Dimension(2))
    {
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("CONVECTIVE_FLUX_Z", d_num_conv_ghosts));
    }
    num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("PRIMITIVE_VARIABLES", d_num_conv_ghosts));
    
    d_flow_model->registerDerivedVariables(num_subghosts_of_data);
    
    basic_utilities->registerDerivedVariablesForCharacteristicProjectionOfPrimitiveVariables(
        d_num_conv_ghosts,
        AVERAGING::SIMPLE);
    
    d_flow_model->allocateMemoryForDerivedCellData();
    
    d_flow_model->computeDerivedCellData();
    
    HAMERS_SHARED_PTR<pdat::CellData<Real> > density = d_flow_model->getCellData("DENSITY");
    HAMERS_SHARED_PTR<pdat::CellData<Real> > velocity = d_flow_model->getCellData("VELOCITY");
    
    std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > convective_flux_node(d_dim.getValue());
    convective_flux_node[0] = d_flow_model->getCellData("CONVECTIVE_FLUX_X");
    if (d_dim > tbox::Dimension(1))
    {
        convective_flux_node[1] = d_flow_model->getCellData("CONVECTIVE_FLUX_Y");
    }
    if (d_dim > tbox::Dimension(2))
    {
        convective_flux_node[2] = d_flow_model->getCellData("CONVECTIVE_FLUX_Z");
    }
    
    /*
     * Get the pointers to the conservative variables and primitive variables.
     * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
     */
    
    std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > conservative_variables =
        d_flow_model->getCellDataOfConservativeVariables();
    
    std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > primitive_variables =
        d_flow_model->getCellDataOfPrimitiveVariables();
    
    std::vector<hier::IntVector> num_subghosts_conservative_var;
    num_subghosts_conservative_var.reserve(d_num_eqn);
    
    std::vector<hier::IntVector> num_subghosts_primitive_var;
        num_subghosts_primitive_var.reserve(d_num_eqn);
    
    std::vector<hier::IntVector> subghostcell_dims_conservative_var;
    subghostcell_dims_conservative_var.reserve(d_num_eqn);
    
    std::vector<hier::IntVector> subghostcell_dims_primitive_var;
        subghostcell_dims_primitive_var.reserve(d_num_eqn);
    
    std::vector<Real*> Q;
    Q.reserve(d_num_eqn);
    
    std::vector<Real*> V;
    V.reserve(d_num_eqn);
    
    int count_eqn = 0;
    
    for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
    {
        int depth = conservative_variables[vi]->getDepth();
        
        for (int di = 0; di < depth; di++)
        {
            // If the last element of the conservative variable vector is not in the system of equations,
            // ignore it.
            if (count_eqn >= d_num_eqn)
                break;
            
            Q.push_back(conservative_variables[vi]->getPointer(di));
            num_subghosts_conservative_var.push_back(conservative_variables[vi]->getGhostCellWidth());
            subghostcell_dims_conservative_var.push_back(
                conservative_variables[vi]->getGhostBox().numberCells());
            
            count_eqn++;
        }
    }
    
    count_eqn = 0;
    
    for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
    {
        int depth = primitive_variables[vi]->getDepth();
        
        for (int di = 0; di < depth; di++)
        {
            // If the last element of the primitive variable vector is not in the system of equations,
            // ignore it.
            if (count_eqn >= d_num_eqn)
                break;
            
            V.push_back(primitive_variables[vi]->getPointer(di));
            num_subghosts_primitive_var.push_back(primitive_variables[vi]->getGhostCellWidth());
            subghostcell_dims_primitive_var.push_back(
                primitive_variables[vi]->getGhostBox().numberCells());
            
            count_eqn++;
        }
    }
    
    /*
     * Pointers to shock sensors.
     */
    
    Real* s   = discontinuity_sensor_cell->getPointer(0);
    Real* s_x = discontinuity_sensor_side->getPointer(0);
    Real* s_y = nullptr;
    Real* s_z = nullptr;
    
    if (d_dim == tbox::Dimension(1))
    {
        discontinuity_sensor_cell->fillAll(Real(1));
        discontinuity_sensor_side->fillAll(Real(1));
        perform_WCNS = true;
    }
    else
    {
        /*
         * Get the numbers of ghost cells of the variables.
         */
        
        const hier::IntVector num_ghosts_density = density->getGhostCellWidth();
        
        /*
         * Get the ghost cell dimensions of of the variables.
         */
        
        const hier::IntVector ghostcell_dims_density = density->getGhostBox().numberCells();
        
        // Get the pointers to the data.
        Real* rho   = density->getPointer(0);
        Real* theta = dilatation->getPointer(0);
        Real* Omega = enstrophy->getPointer(0);
        Real* rho_s_x = density_sensors[0]->getPointer(0);
        Real* rho_s_y = density_sensors[1]->getPointer(0);
        
        s_y = discontinuity_sensor_side->getPointer(1);
        
        const hier::Box empty_box(d_dim);
        
        if (d_dim == tbox::Dimension(2))
        {
            /*
             * Get the local lower indices and the number of cells in each dimension.
             */
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_lo_1 = domain_lo[1];
            const int domain_dim_0 = domain_dims[0];
            const int domain_dim_1 = domain_dims[1];
            
            const int num_ghosts_0_density = num_ghosts_density[0];
            const int num_ghosts_1_density = num_ghosts_density[1];
            const int ghostcell_dim_0_density = ghostcell_dims_density[0];
            
            /*
             * Get the interior dimensions.
             */
            
            const int interior_dim_0 = interior_dims[0];
            const int interior_dim_1 = interior_dims[1];
            
            /*
             * Compute the derivatives of velocity, dilatation and vorticity magnitude.
             */
            
            HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_x(
                new DerivativeFirstOrder("first order derivative in x-direction", d_dim, DIRECTION::X_DIRECTION, 1));
            
            HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_y(
                new DerivativeFirstOrder("first order derivative in y-direction", d_dim, DIRECTION::Y_DIRECTION, 1));
            
            // Compute dudx.
            derivative_first_order_x->computeDerivative(
                velocity_derivatives,
                velocity,
                Real(dx[0]),
                empty_box,
                0,
                0);
            
            // Compute dudy.
            derivative_first_order_y->computeDerivative(
                velocity_derivatives,
                velocity,
                Real(dx[1]),
                empty_box,
                1,
                0);
            
            // Compute dvdx.
            derivative_first_order_x->computeDerivative(
                velocity_derivatives,
                velocity,
                Real(dx[0]),
                empty_box,
                2,
                1);
            
            // Compute dvdy.
            derivative_first_order_y->computeDerivative(
                velocity_derivatives,
                velocity,
                Real(dx[1]),
                empty_box,
                3,
                1);
            
            // Get the pointers to the cell data of velocity derivatives.
            Real* dudx = velocity_derivatives->getPointer(0);
            Real* dudy = velocity_derivatives->getPointer(1);
            Real* dvdx = velocity_derivatives->getPointer(2);
            Real* dvdy = velocity_derivatives->getPointer(3);
            
            // Compute the dilatation and the enstrophy.
            for (int j = domain_lo_1 - 1; j < domain_lo_1 + domain_dim_1 + 1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0 - 1; i < domain_lo_0 + domain_dim_0 + 1; i++)
                {
                    // Compute the linear indices.
                    const int idx = (i + 1) +
                        (j + 1)*(interior_dim_0 + 2);
                    
                    const int idx_rho = (i + num_ghosts_0_density) +
                        (j + num_ghosts_1_density)*ghostcell_dim_0_density;
                    
                    const int idx_rho_x_L = (i - 1 + num_ghosts_0_density) +
                        (j + num_ghosts_1_density)*ghostcell_dim_0_density;
                    
                    const int idx_rho_x_R = (i + 1 + num_ghosts_0_density) +
                        (j + num_ghosts_1_density)*ghostcell_dim_0_density;
                    
                    const int idx_rho_y_B = (i + num_ghosts_0_density) +
                        (j - 1 + num_ghosts_1_density)*ghostcell_dim_0_density;
                    
                    const int idx_rho_y_T = (i + num_ghosts_0_density) +
                        (j + 1 + num_ghosts_1_density)*ghostcell_dim_0_density;
                    
                    theta[idx] = dudx[idx] + dvdy[idx];
                    Omega[idx] = (dvdx[idx] - dudy[idx])*(dvdx[idx] - dudy[idx]);
                    
                    rho_s_x[idx] = std::abs(rho[idx_rho_x_R] - Real(2)*rho[idx_rho] + rho[idx_rho_x_L])/
                        (rho[idx_rho_x_R] + Real(2)*rho[idx_rho] + rho[idx_rho_x_L] + EPSILON);
                    
                    rho_s_y[idx] = std::abs(rho[idx_rho_y_T] - Real(2)*rho[idx_rho] + rho[idx_rho_y_B])/
                        (rho[idx_rho_y_T] + Real(2)*rho[idx_rho] + rho[idx_rho_y_B] + EPSILON);
                }
            }
            
            int count_valid_sensor = 0;
            
            const Real half = Real(1)/Real(2);
            
            // Compute the shock sensor on the side in the x-direction.
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                {
                    // Compute the linear indices.
                    const int idx_sensor = i +
                        j*(interior_dim_0 + 1);
                    
                    const int idx_L = (i - 1 + 1) +
                        (j + 1)*(interior_dim_0 + 2);
                    
                    const int idx_R = (i + 0 + 1) +
                        (j + 1)*(interior_dim_0 + 2);
                    
                    const Real Ducros_value_L = (-theta[idx_L]*std::abs(theta[idx_L]))/
                        (theta[idx_L]*theta[idx_L] + Omega[idx_L] + EPSILON);
                    
                    const Real Ducros_value_R = (-theta[idx_R]*std::abs(theta[idx_R]))/
                        (theta[idx_R]*theta[idx_R] + Omega[idx_R] + EPSILON);
                    
                    const Real Ducros_value_midpoint = half*(Ducros_value_L + Ducros_value_R);
                    
                    s_x[idx_sensor] = Real(0);
                    if (Ducros_value_midpoint > d_threshold_sensor_shock && use_shock_capturing)
                    {
                        s_x[idx_sensor] = Real(1);
                    }
                    
                    const Real rho_s_x_midpoint = half*(rho_s_x[idx_L] + rho_s_x[idx_R]);
                    if (rho_s_x_midpoint > d_threshold_sensor_interface && use_interface_capturing)
                    {
                        s_x[idx_sensor] = Real(1);
                    }
                    
                    if (s_x[idx_sensor] > Real(0))
                    {
                        count_valid_sensor++;
                    }
                }
            }
            
            // Compute the shock sensor on the side in the y-direction.
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_sensor = i +
                        j*interior_dim_0;
                    
                    const int idx_B = (i + 1) +
                        (j - 1 + 1)*(interior_dim_0 + 2);
                    
                    const int idx_T = (i + 1) +
                        (j + 0 + 1)*(interior_dim_0 + 2);
                    
                    const Real Ducros_value_B = (-theta[idx_B]*std::abs(theta[idx_B]))/
                        (theta[idx_B]*theta[idx_B] + Omega[idx_B] + EPSILON);
                    
                    const Real Ducros_value_T = (-theta[idx_T]*std::abs(theta[idx_T]))/
                        (theta[idx_T]*theta[idx_T] + Omega[idx_T] + EPSILON);
                    
                    const Real Ducros_value_midpoint = half*(Ducros_value_B + Ducros_value_T);
                    
                    s_y[idx_sensor] = Real(0);
                    if (Ducros_value_midpoint > d_threshold_sensor_shock && use_shock_capturing)
                    {
                        s_y[idx_sensor] = Real(1);
                    }
                    
                    const Real rho_s_y_midpoint = half*(rho_s_y[idx_B] + rho_s_y[idx_T]);
                    if (rho_s_y_midpoint > d_threshold_sensor_interface && use_interface_capturing)
                    {
                        s_y[idx_sensor] = Real(1);
                    }
                    
                    if (s_y[idx_sensor] > Real(0))
                    {
                        count_valid_sensor++;
                    }
                }
            }
            
            if (d_has_advective_eqn_form)
            {
                // Compute the cell-based Ducros-like shock sensor.
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_sensor = i +
                            j*interior_dim_0;
                        
                        const int idx = (i + 1) +
                            (j + 1)*(interior_dim_0 + 2);
                        
                        const Real Ducros_value = (-theta[idx]*std::abs(theta[idx]))/
                            (theta[idx]*theta[idx] + Omega[idx] + EPSILON);
                        
                        s[idx_sensor] = Real(0);
                        if (Ducros_value > d_threshold_sensor_shock && use_shock_capturing)
                        {
                            s[idx_sensor] = Real(1);
                        }
                        {
                            s[idx_sensor] = Real(1);
                        }
                        if (rho_s_x[idx] > d_threshold_sensor_interface && use_interface_capturing)
                        {
                            s[idx_sensor] = Real(1);
                        }
                        if (rho_s_y[idx] > d_threshold_sensor_interface && use_interface_capturing)
                        {
                            s[idx_sensor] = Real(1);
                        }
                        
                        if (s[idx_sensor] > Real(0))
                        {
                            count_valid_sensor++;
                        }
                    }
                }
            }
            
            if (count_valid_sensor > 0)
            {
                perform_WCNS = true;
            }
        }
        else if (d_dim == tbox::Dimension(3))
        {
            Real* rho_s_z = density_sensors[2]->getPointer(0);
            Real* s_z = discontinuity_sensor_side->getPointer(2);
            
            /*
             * Get the local lower indices and the number of cells in each dimension.
             */
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_lo_1 = domain_lo[1];
            const int domain_lo_2 = domain_lo[2];
            const int domain_dim_0 = domain_dims[0];
            const int domain_dim_1 = domain_dims[1];
            const int domain_dim_2 = domain_dims[2];
            
            const int num_ghosts_0_density = num_ghosts_density[0];
            const int num_ghosts_1_density = num_ghosts_density[1];
            const int num_ghosts_2_density = num_ghosts_density[2];
            const int ghostcell_dim_0_density = ghostcell_dims_density[0];
            const int ghostcell_dim_1_density = ghostcell_dims_density[1];
            
            /*
             * Get the interior dimensions.
             */
            
            const int interior_dim_0 = interior_dims[0];
            const int interior_dim_1 = interior_dims[1];
            const int interior_dim_2 = interior_dims[2];
            
            /*
             * Compute the derivatives of velocity, dilatation and vorticity magnitude.
             */
            
            HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_x(
                new DerivativeFirstOrder("first order derivative in x-direction", d_dim, DIRECTION::X_DIRECTION, 1));
            
            HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_y(
                new DerivativeFirstOrder("first order derivative in y-direction", d_dim, DIRECTION::Y_DIRECTION, 1));
            
            HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_z(
                new DerivativeFirstOrder("first order derivative in z-direction", d_dim, DIRECTION::Z_DIRECTION, 1));
            
            // Compute dudx.
            derivative_first_order_x->computeDerivative(
                velocity_derivatives,
                velocity,
                Real(dx[0]),
                empty_box,
                0,
                0);
            
            // Compute dudy.
            derivative_first_order_y->computeDerivative(
                velocity_derivatives,
                velocity,
                Real(dx[1]),
                empty_box,
                1,
                0);
            
            // Compute dudz.
            derivative_first_order_z->computeDerivative(
                velocity_derivatives,
                velocity,
                Real(dx[2]),
                empty_box,
                2,
                0);
            
            // Compute dvdx.
            derivative_first_order_x->computeDerivative(
                velocity_derivatives,
                velocity,
                Real(dx[0]),
                empty_box,
                3,
                1);
            
            // Compute dvdy.
            derivative_first_order_y->computeDerivative(
                velocity_derivatives,
                velocity,
                Real(dx[1]),
                empty_box,
                4,
                1);
            
            // Compute dvdz.
            derivative_first_order_z->computeDerivative(
                velocity_derivatives,
                velocity,
                Real(dx[2]),
                empty_box,
                5,
                1);
            
            // Compute dwdx.
            derivative_first_order_x->computeDerivative(
                velocity_derivatives,
                velocity,
                Real(dx[0]),
                empty_box,
                6,
                2);
            
            // Compute dwdy.
            derivative_first_order_y->computeDerivative(
                velocity_derivatives,
                velocity,
                Real(dx[1]),
                empty_box,
                7,
                2);
            
            // Compute dwdz.
            derivative_first_order_z->computeDerivative(
                velocity_derivatives,
                velocity,
                Real(dx[2]),
                empty_box,
                8,
                2);
            
            // Get the pointers to the cell data of velocity derivatives.
            Real* dudx = velocity_derivatives->getPointer(0);
            Real* dudy = velocity_derivatives->getPointer(1);
            Real* dudz = velocity_derivatives->getPointer(2);
            Real* dvdx = velocity_derivatives->getPointer(3);
            Real* dvdy = velocity_derivatives->getPointer(4);
            Real* dvdz = velocity_derivatives->getPointer(5);
            Real* dwdx = velocity_derivatives->getPointer(6);
            Real* dwdy = velocity_derivatives->getPointer(7);
            Real* dwdz = velocity_derivatives->getPointer(8);
            
            // Compute the dilatation and the enstrophy.
            for (int k = domain_lo_2 - 1; k < domain_lo_2 + domain_dim_2 + 1; k++)
            {
                for (int j = domain_lo_1 - 1; j < domain_lo_1 + domain_dim_1 + 1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0 - 1; i < domain_lo_0 + domain_dim_0 + 1; i++)
                    {
                        // Compute the linear index.
                        const int idx = (i + 1) +
                            (j + 1)*(interior_dim_0 + 2) +
                            (k + 1)*(interior_dim_0 + 2)*
                                (interior_dim_1 + 2);
                        
                        theta[idx] = dudx[idx] + dvdy[idx] + dwdz[idx];
                        
                        const Real omega_x = dwdy[idx] - dvdz[idx];
                        const Real omega_y = dudz[idx] - dwdx[idx];
                        const Real omega_z = dvdx[idx] - dudy[idx];
                        
                        Omega[idx] = omega_x*omega_x + omega_y*omega_y + omega_z*omega_z;
                    }
                }
            }
            
            int count_valid_sensor = 0;
            
            const Real half = Real(1)/Real(2);
            
            // Compute the shock sensor on the side in the x-direction.
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                    {
                        // Compute the linear indices.
                        const int idx_sensor = i +
                            j*(interior_dim_0 + 1) +
                            k*(interior_dim_0 + 1)*
                                interior_dim_1;
                        
                        const int idx_L = (i - 1 + 1) +
                            (j + 1)*(interior_dim_0 + 2) +
                            (k + 1)*(interior_dim_0 + 2)*
                                (interior_dim_1 + 2);
                        
                        const int idx_R = (i + 0 + 1) +
                            (j + 1)*(interior_dim_0 + 2) +
                            (k + 1)*(interior_dim_0 + 2)*
                                (interior_dim_1 + 2);
                        
                        const Real Ducros_value_L = (-theta[idx_L]*std::abs(theta[idx_L]))/
                            (theta[idx_L]*theta[idx_L] + Omega[idx_L] + EPSILON);
                        
                        const Real Ducros_value_R = (-theta[idx_R]*std::abs(theta[idx_R]))/
                            (theta[idx_R]*theta[idx_R] + Omega[idx_R] + EPSILON);
                        
                        const Real Ducros_value_midpoint = half*(Ducros_value_L + Ducros_value_R);
                        
                        s_x[idx_sensor] = Real(0);
                        if (Ducros_value_midpoint > d_threshold_sensor_shock && use_shock_capturing)
                        {
                            s_x[idx_sensor] = Real(1);
                        }
                        
                        const Real rho_s_x_midpoint = half*(rho_s_x[idx_L] + rho_s_x[idx_R]);
                        if (rho_s_x_midpoint > d_threshold_sensor_interface && use_interface_capturing)
                        {
                            s_x[idx_sensor] = Real(1);
                        }
                        
                        if (s_x[idx_sensor] > Real(0))
                        {
                            count_valid_sensor++;
                        }
                    }
                }
            }
            
            // Compute the shock sensor on the side in the y-direction.
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_sensor = i +
                            j*interior_dim_0 +
                            k*interior_dim_0*
                                (interior_dim_1 + 1);
                        
                        const int idx_B = (i + 1) +
                            (j - 1 + 1)*(interior_dim_0 + 2) +
                            (k + 1)*(interior_dim_0 + 2)*
                                (interior_dim_1 + 2);
                        
                        const int idx_T = (i + 1) +
                            (j + 0 + 1)*(interior_dim_0 + 2) +
                            (k + 1)*(interior_dim_0 + 2)*
                                (interior_dim_1 + 2);
                        
                        const Real Ducros_value_B = (-theta[idx_B]*std::abs(theta[idx_B]))/
                            (theta[idx_B]*theta[idx_B] + Omega[idx_B] + EPSILON);
                        
                        const Real Ducros_value_T = (-theta[idx_T]*std::abs(theta[idx_T]))/
                            (theta[idx_T]*theta[idx_T] + Omega[idx_T] + EPSILON);
                        
                        const Real Ducros_value_midpoint = half*(Ducros_value_B + Ducros_value_T);
                        
                        s_y[idx_sensor] = Real(0);
                        if (Ducros_value_midpoint > d_threshold_sensor_shock && use_shock_capturing)
                        {
                            s_y[idx_sensor] = Real(1);
                        }
                        
                        const Real rho_s_y_midpoint = half*(rho_s_y[idx_B] + rho_s_y[idx_T]);
                        if (rho_s_y_midpoint > d_threshold_sensor_interface && use_interface_capturing)
                        {
                            s_y[idx_sensor] = Real(1);
                        }
                        
                        if (s_y[idx_sensor] > Real(0))
                        {
                            count_valid_sensor++;
                        }
                    }
                }
            }
            
            // Compute the shock sensor on the side in the z-direction.
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2 + 1; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_sensor = i +
                            j*interior_dim_0 +
                            k*interior_dim_0*
                                interior_dim_1;
                        
                        const int idx_B = (i + 1) +
                            (j + 1)*(interior_dim_0 + 2) +
                            (k - 1 + 1)*(interior_dim_0 + 2)*
                                (interior_dim_1 + 2);
                        
                        const int idx_T = (i + 1) +
                            (j + 1)*(interior_dim_0 + 2) +
                            (k + 0 + 1)*(interior_dim_0 + 2)*
                                (interior_dim_1 + 2);
                        
                        const Real Ducros_value_B = (-theta[idx_B]*std::abs(theta[idx_B]))/
                            (theta[idx_B]*theta[idx_B] + Omega[idx_B] + EPSILON);
                        
                        const Real Ducros_value_T = (-theta[idx_T]*std::abs(theta[idx_T]))/
                            (theta[idx_T]*theta[idx_T] + Omega[idx_T] + EPSILON);
                        
                        const Real Ducros_value_midpoint = half*(Ducros_value_B + Ducros_value_T);
                        
                        s_z[idx_sensor] = Real(0);
                        if (Ducros_value_midpoint > d_threshold_sensor_shock && use_shock_capturing)
                        {
                            s_z[idx_sensor] = Real(1);
                        }
                        
                        const Real rho_s_z_midpoint = half*(rho_s_z[idx_B] + rho_s_z[idx_T]);
                        if (rho_s_z_midpoint > d_threshold_sensor_interface && use_interface_capturing)
                        {
                            s_z[idx_sensor] = Real(1);
                        }
                        
                        if (s_z[idx_sensor] > Real(0))
                        {
                            count_valid_sensor++;
                        }
                    }
                }
            }
            
            if (d_has_advective_eqn_form)
            {
                // Compute the cell-based Ducros-like shock sensor.
                for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                {
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_sensor = i +
                                j*interior_dim_0 +
                                k*interior_dim_0*
                                    interior_dim_1;
                            
                            const int idx = (i + 1) +
                                (j + 1)*(interior_dim_0 + 2) +
                                (k + 1)*(interior_dim_0 + 2)*
                                    (interior_dim_1 + 2);
                            
                            const Real Ducros_value = (-theta[idx]*std::abs(theta[idx]))/
                                (theta[idx]*theta[idx] + Omega[idx] + EPSILON);
                            
                            s[idx_sensor] = Real(0);
                            if (Ducros_value > d_threshold_sensor_shock && use_shock_capturing)
                            {
                                s[idx_sensor] = Real(1);
                            }
                            
                            if (rho_s_x[idx] > d_threshold_sensor_interface && use_interface_capturing)
                            {
                                s[idx_sensor] = Real(1);
                            }
                            if (rho_s_y[idx] > d_threshold_sensor_interface && use_interface_capturing)
                            {
                                s[idx_sensor] = Real(1);
                            }
                            if (rho_s_z[idx] > d_threshold_sensor_interface && use_interface_capturing)
                            {
                                s[idx_sensor] = Real(1);
                            }
                            
                            if (s[idx_sensor] > Real(0))
                            {
                                count_valid_sensor++;
                            }
                        }
                    }
                }
            }
            
            if (count_valid_sensor > 0)
            {
                perform_WCNS = true;
            }
        } // if (d_dim == tbox::Dimension(3))
    } // if (d_dim == tbox::Dimension(1))
    
    if (perform_WCNS)
    {
        /*
         * Declare temporary data containers for WENO interpolation.
         */
        
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > projection_variables;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > > characteristic_variables;
        
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > characteristic_variables_minus;
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > characteristic_variables_plus;
        
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > primitive_variables_minus;
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > primitive_variables_plus;
        
        HAMERS_SHARED_PTR<pdat::SideData<int> > bounded_flag_minus;
        HAMERS_SHARED_PTR<pdat::SideData<int> > bounded_flag_plus;
        
        /*
         * Initialize temporary data containers for WENO interpolation.
         */
        
        const int num_projection_var = basic_utilities->getNumberOfProjectionVariablesForPrimitiveVariables();
        projection_variables.reserve(num_projection_var);
        
        for (int vi = 0; vi < num_projection_var; vi++)
        {
            projection_variables.push_back(HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
                interior_box, 1, hier::IntVector::getZero(d_dim))); // No ghost midpoints.
        }
        
        characteristic_variables.resize(6);
        
        for (int m = 0; m < 6; m++)
        {
            characteristic_variables[m].reserve(d_num_eqn);
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                characteristic_variables[m].push_back(HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
                    interior_box, 1, hier::IntVector::getZero(d_dim))); // No ghost midpoints.
            }
        }
        
        characteristic_variables_minus.reserve(d_num_eqn);
        characteristic_variables_plus.reserve(d_num_eqn);
        primitive_variables_minus.reserve(d_num_eqn);
        primitive_variables_plus.reserve(d_num_eqn);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            characteristic_variables_minus.push_back(HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
                interior_box, 1, hier::IntVector::getZero(d_dim))); // No ghost midpoints.
            
            characteristic_variables_plus.push_back(HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
                interior_box, 1, hier::IntVector::getZero(d_dim))); // No ghost midpoints.
            
            primitive_variables_minus.push_back(HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
                interior_box, 1, hier::IntVector::getZero(d_dim))); // No ghost midpoints.
            
            primitive_variables_plus.push_back(HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
                interior_box, 1, hier::IntVector::getZero(d_dim))); // No ghost midpoints.
        }
        
        bounded_flag_minus.reset(
            new pdat::SideData<int>(interior_box, 1, hier::IntVector::getZero(d_dim)));
        
        bounded_flag_plus.reset(
            new pdat::SideData<int>(interior_box, 1, hier::IntVector::getZero(d_dim)));
        
        /*
         * Compute the side data of the projection variables for transformation between primitive variables and
         * characteristic variables.
         */
        
        basic_utilities->computeSideDataOfProjectionVariablesForPrimitiveVariables(
            projection_variables,
            domains);
        
        /*
         * Transform primitive variables to characteristic variables.
         */
        
        for (int m = 0; m < 6; m++)
        {
            basic_utilities->computeSideDataOfCharacteristicVariablesFromPrimitiveVariables(
                characteristic_variables[m],
                primitive_variables,
                projection_variables,
                m - 3,
                domains);
        }
        
        /*
         * Peform WENO interpolation.
         */
        
        performWENOInterpolation(
            characteristic_variables_minus,
            characteristic_variables_plus,
            characteristic_variables,
            domains);
        
        /*
         * Transform characteristic variables back to primitive variables.
         */
        
        basic_utilities->computeSideDataOfPrimitiveVariablesFromCharacteristicVariables(
            primitive_variables_minus,
            characteristic_variables_minus,
            projection_variables,
            domains);
        
        basic_utilities->computeSideDataOfPrimitiveVariablesFromCharacteristicVariables(
            primitive_variables_plus,
            characteristic_variables_plus,
            projection_variables,
            domains);
        
        /*
         * Check whether the interpolated side primitive variables are within the bounds.
         */
        
        basic_utilities->checkSideDataOfPrimitiveVariablesBounded(
            bounded_flag_minus,
            primitive_variables_minus,
            domains);
        
        basic_utilities->checkSideDataOfPrimitiveVariablesBounded(
            bounded_flag_plus,
            primitive_variables_plus,
            domains);
        
        /*
         * Declare containers to store pointers for computing mid-point fluxes.
         */
        
        std::vector<Real*> V_minus;
        std::vector<Real*> V_plus;
        V_minus.resize(d_num_eqn);
        V_plus.resize(d_num_eqn);
        
        int* flag_minus = nullptr;
        int* flag_plus = nullptr;
        
        /*
         * Coefficients for finite differencing.
         */
        
        const Real phi = Real(256)/Real(175);
        
        const Real a_midpoint_r = phi;
        
        const Real a_node_r = -(Real(75)/Real(128)*phi - Real(37)/Real(60));
        const Real b_node_r = Real(25)/Real(256)*phi - Real(2)/Real(15);
        const Real c_node_r = -(Real(3)/Real(256)*phi - Real(1)/Real(60));
        
        const Real a_midpoint = phi;
        
        const Real a_node = -(Real(175)*phi - Real(192))/Real(256);
        const Real b_node = (Real(35)*phi - Real(48))/Real(320);
        const Real c_node = -(Real(45)*phi - Real(64))/Real(3840);
        
        /*
         * Compute the convective flux and source using shock-capturing scheme.
         */
        if (d_dim == tbox::Dimension(1))
        {
            /*
             * Get the local lower index and the number of cells in each dimension.
             */
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_dim_0 = domain_dims[0];
            
            /*
             * Get the pointers to the velocity and convective flux cell data inside the flow model.
             * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
             */
            
            HAMERS_SHARED_PTR<pdat::CellData<Real> > velocity = d_flow_model->getCellData("VELOCITY");
            
            std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > convective_flux_node(1);
            convective_flux_node[0] = d_flow_model->getCellData("CONVECTIVE_FLUX_X");
            
            hier::IntVector num_subghosts_velocity = velocity->getGhostCellWidth();
            hier::IntVector subghostcell_dims_velocity = velocity->getGhostBox().numberCells();
            
            hier::IntVector num_subghosts_convective_flux_x = convective_flux_node[0]->getGhostCellWidth();
            hier::IntVector subghostcell_dims_convective_flux_x = convective_flux_node[0]->getGhostBox().numberCells();
            
            const int num_subghosts_0_convective_flux_x = num_subghosts_convective_flux_x[0];
            const int subghostcell_dim_0_convective_flux_x = subghostcell_dims_convective_flux_x[0];
            
            const int num_subghosts_0_velocity = num_subghosts_velocity[0];
            const int subghostcell_dim_0_velocity = subghostcell_dims_velocity[0];
            
            Real* u = velocity->getPointer(0);
            
            std::vector<Real*> F_node_x;
            F_node_x.reserve(d_num_eqn);
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_node_x.push_back(convective_flux_node[0]->getPointer(ei));
            }
            
            std::vector<Real*> F_midpoint_x;
            F_midpoint_x.reserve(d_num_eqn);
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_midpoint_x.push_back(convective_flux_midpoint->getPointer(0, ei));
            }
            
            /*
             * Use first order interpolation if interpolated side primitive variables in x-direction
             * are out of bounds.
             */
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                V_minus[ei] = primitive_variables_minus[ei]->getPointer(0);
                V_plus[ei] = primitive_variables_plus[ei]->getPointer(0);
            }
            
            flag_minus = bounded_flag_minus->getPointer(0);
            flag_plus = bounded_flag_plus->getPointer(0);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                const int num_subghosts_0_primitive_var = num_subghosts_primitive_var[ei][0];
                const int num_subghosts_1_primitive_var = num_subghosts_primitive_var[ei][1];
                const int subghostcell_dim_0_primitive_var = subghostcell_dims_primitive_var[ei][0];
                
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                {
                    // Compute the linear indices.
                    const int idx_midpoint_x = i + 1;
                    
                    const int idx_cell_L = i - 1 + num_subghosts_0_primitive_var;
                    const int idx_cell_R = i     + num_subghosts_0_primitive_var;
                    
                    if (flag_minus[idx_midpoint_x] == 0 || flag_plus[idx_midpoint_x] == 0)
                    {
                        V_minus[ei][idx_midpoint_x] = V[ei][idx_cell_L];
                        V_plus[ei][idx_midpoint_x]  = V[ei][idx_cell_R];
                    }
                }
            }
            
            /*
             * Compute mid-point flux in the x-direction.
             */
            
            if (d_has_advective_eqn_form)
            {
                riemann_solver->computeConvectiveFluxAndVelocityFromPrimitiveVariables(
                    convective_flux_midpoint,
                    velocity_midpoint,
                    primitive_variables_minus,
                    primitive_variables_plus,
                    DIRECTION::X_DIRECTION,
                    RIEMANN_SOLVER::HLLC,
                    domains[0]);
            }
            else
            {
                riemann_solver->computeConvectiveFluxFromPrimitiveVariables(
                    convective_flux_midpoint,
                    primitive_variables_minus,
                    primitive_variables_plus,
                    DIRECTION::X_DIRECTION,
                    RIEMANN_SOLVER::HLLC,
                    domains[0]);
            }
            
            /*
             * Reconstruct the flux in the x-direction.
             */
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                Real* F_face_x = convective_flux->getPointer(0, ei);
                
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                {
                    // Compute the linear indices.
                    const int idx_face_x = i;
                    
                    if (s_x[idx_face_x] > Real(0))
                    {
                        const int idx_midpoint_x = i;
                        
                        const int idx_node_LLL = i - 3 + num_subghosts_0_convective_flux_x;
                        const int idx_node_LL  = i - 2 + num_subghosts_0_convective_flux_x;
                        const int idx_node_L   = i - 1 + num_subghosts_0_convective_flux_x;
                        const int idx_node_R   = i     + num_subghosts_0_convective_flux_x;
                        const int idx_node_RR  = i + 1 + num_subghosts_0_convective_flux_x;
                        const int idx_node_RRR = i + 2 + num_subghosts_0_convective_flux_x;
                        
                        F_face_x[idx_face_x] = Real(dt)*(
                            a_midpoint_r*F_midpoint_x[ei][idx_midpoint_x] +
                            a_node_r*(F_node_x[ei][idx_node_L]   + F_node_x[ei][idx_node_R]) +
                            b_node_r*(F_node_x[ei][idx_node_LL]  + F_node_x[ei][idx_node_RR]) +
                            c_node_r*(F_node_x[ei][idx_node_LLL] + F_node_x[ei][idx_node_RRR])
                            );
                    }
                }
            }
            
            /*
             * Compute the source.
             */
            
            if (d_has_advective_eqn_form)
            {
                Real* u_midpoint_x = velocity_midpoint->getPointer(0, 0);
                
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    if (d_eqn_form[ei] == EQN_FORM::ADVECTIVE)
                    {
                        Real* S = source_scratch->getPointer(ei);
                        
                        const int num_subghosts_0_conservative_var = num_subghosts_conservative_var[ei][0];
                        
                        HAMERS_PRAGMA_SIMD
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_cell_nghost = i;
                            
                            if (s[idx_cell_nghost] > Real(0))
                            {
                                const int idx_cell_wghost = i + num_subghosts_0_conservative_var;
                                
                                const int idx_cell_wghost_x_LLL = i - 3 + num_subghosts_0_velocity;
                                const int idx_cell_wghost_x_LL  = i - 2 + num_subghosts_0_velocity;
                                const int idx_cell_wghost_x_L   = i - 1 + num_subghosts_0_velocity;
                                const int idx_cell_wghost_x_R   = i + 1 + num_subghosts_0_velocity;
                                const int idx_cell_wghost_x_RR  = i + 2 + num_subghosts_0_velocity;
                                const int idx_cell_wghost_x_RRR = i + 3 + num_subghosts_0_velocity;
                                
                                const int idx_midpoint_x_L = i;
                                const int idx_midpoint_x_R = i + 1;
                                
                                S[idx_cell_nghost] = Real(dt)*Q[ei][idx_cell_wghost]*((
                                    a_midpoint*(u_midpoint_x[idx_midpoint_x_R] - u_midpoint_x[idx_midpoint_x_L]) +
                                    a_node*(u[idx_cell_wghost_x_R]   - u[idx_cell_wghost_x_L]) +
                                    b_node*(u[idx_cell_wghost_x_RR]  - u[idx_cell_wghost_x_LL]) +
                                    c_node*(u[idx_cell_wghost_x_RRR] - u[idx_cell_wghost_x_LLL])
                                    )/Real(dx[0]));
                            }
                        }
                    }
                }
            }
        }
        else if (d_dim == tbox::Dimension(2))
        {
            /*
             * Get the local lower indices and the number of cells in each dimension.
             */
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_lo_1 = domain_lo[1];
            const int domain_dim_0 = domain_dims[0];
            const int domain_dim_1 = domain_dims[1];
            
            /*
             * Get the interior dimension.
             */
            
            const int interior_dim_0 = interior_dims[0];
            
            /*
             * Get the pointers to the velocity and convective flux cell data inside the flow model.
             * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
             */
            
            HAMERS_SHARED_PTR<pdat::CellData<Real> > velocity = d_flow_model->getCellData("VELOCITY");
            
            std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > convective_flux_node(2);
            convective_flux_node[0] = d_flow_model->getCellData("CONVECTIVE_FLUX_X");
            convective_flux_node[1] = d_flow_model->getCellData("CONVECTIVE_FLUX_Y");
            
            hier::IntVector num_subghosts_velocity = velocity->getGhostCellWidth();
            hier::IntVector subghostcell_dims_velocity = velocity->getGhostBox().numberCells();
            
            hier::IntVector num_subghosts_convective_flux_x = convective_flux_node[0]->getGhostCellWidth();
            hier::IntVector subghostcell_dims_convective_flux_x = convective_flux_node[0]->getGhostBox().numberCells();
            
            hier::IntVector num_subghosts_convective_flux_y = convective_flux_node[1]->getGhostCellWidth();
            hier::IntVector subghostcell_dims_convective_flux_y = convective_flux_node[1]->getGhostBox().numberCells();
            
            const int num_subghosts_0_velocity = num_subghosts_velocity[0];
            const int num_subghosts_1_velocity = num_subghosts_velocity[1];
            const int subghostcell_dim_0_velocity = subghostcell_dims_velocity[0];
            
            const int num_subghosts_0_convective_flux_x = num_subghosts_convective_flux_x[0];
            const int num_subghosts_1_convective_flux_x = num_subghosts_convective_flux_x[1];
            const int subghostcell_dim_0_convective_flux_x = subghostcell_dims_convective_flux_x[0];
            
            const int num_subghosts_0_convective_flux_y = num_subghosts_convective_flux_y[0];
            const int num_subghosts_1_convective_flux_y = num_subghosts_convective_flux_y[1];
            const int subghostcell_dim_0_convective_flux_y = subghostcell_dims_convective_flux_y[0];
            
            Real* u = velocity->getPointer(0);
            Real* v = velocity->getPointer(1);
            
            std::vector<Real*> F_node_x;
            std::vector<Real*> F_node_y;
            F_node_x.reserve(d_num_eqn);
            F_node_y.reserve(d_num_eqn);
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_node_x.push_back(convective_flux_node[0]->getPointer(ei));
                F_node_y.push_back(convective_flux_node[1]->getPointer(ei));
            }
            
            std::vector<Real*> F_midpoint_x;
            std::vector<Real*> F_midpoint_y;
            F_midpoint_x.reserve(d_num_eqn);
            F_midpoint_y.reserve(d_num_eqn);
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_midpoint_x.push_back(convective_flux_midpoint->getPointer(0, ei));
                F_midpoint_y.push_back(convective_flux_midpoint->getPointer(1, ei));
            }
            
            /*
             * Use first order interpolation if interpolated side primitive variables in x-direction
             * are out of bounds.
             */
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                V_minus[ei] = primitive_variables_minus[ei]->getPointer(0);
                V_plus[ei] = primitive_variables_plus[ei]->getPointer(0);
            }
            
            flag_minus = bounded_flag_minus->getPointer(0);
            flag_plus = bounded_flag_plus->getPointer(0);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                const int num_subghosts_0_primitive_var = num_subghosts_primitive_var[ei][0];
                const int num_subghosts_1_primitive_var = num_subghosts_primitive_var[ei][1];
                const int subghostcell_dim_0_primitive_var = subghostcell_dims_primitive_var[ei][0];
                
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                    {
                        // Compute the linear indices.
                        const int idx_midpoint_x = i +
                            j*(interior_dim_0 + 1);
                        
                        const int idx_cell_L = (i - 1 + num_subghosts_0_primitive_var) +
                            (j + num_subghosts_1_primitive_var)*subghostcell_dim_0_primitive_var;
                        
                        const int idx_cell_R = (i + num_subghosts_0_primitive_var) +
                            (j + num_subghosts_1_primitive_var)*subghostcell_dim_0_primitive_var;
                        
                        if (flag_minus[idx_midpoint_x] == 0 || flag_plus[idx_midpoint_x] == 0)
                        {
                            V_minus[ei][idx_midpoint_x] = V[ei][idx_cell_L];
                            V_plus[ei][idx_midpoint_x] = V[ei][idx_cell_R];
                        }
                    }
                }
            }
            
            /*
             * Use first order interpolation if interpolated side primitive variables in y-direction
             * are out of bounds.
             */
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                V_minus[ei] = primitive_variables_minus[ei]->getPointer(1);
                V_plus[ei] = primitive_variables_plus[ei]->getPointer(1);
            }
            
            flag_minus = bounded_flag_minus->getPointer(1);
            flag_plus = bounded_flag_plus->getPointer(1);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                const int num_subghosts_0_primitive_var = num_subghosts_primitive_var[ei][0];
                const int num_subghosts_1_primitive_var = num_subghosts_primitive_var[ei][1];
                const int subghostcell_dim_0_primitive_var = subghostcell_dims_primitive_var[ei][0];
                
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_midpoint_y = i +
                            j*interior_dim_0;
                        
                        const int idx_cell_B = (i + num_subghosts_0_primitive_var) +
                            (j - 1 + num_subghosts_1_primitive_var)*subghostcell_dim_0_primitive_var;
                        
                        const int idx_cell_T = (i + num_subghosts_0_primitive_var) +
                            (j + num_subghosts_1_primitive_var)*subghostcell_dim_0_primitive_var;
                        
                        if (flag_minus[idx_midpoint_y] == 0 || flag_plus[idx_midpoint_y] == 0)
                        {
                            V_minus[ei][idx_midpoint_y] = V[ei][idx_cell_B];
                            V_plus[ei][idx_midpoint_y] = V[ei][idx_cell_T];
                        }
                    }
                }
            }
            
            /*
             * Compute mid-point flux in the x-direction.
             */
            
            if (d_has_advective_eqn_form)
            {
                riemann_solver->computeConvectiveFluxAndVelocityFromPrimitiveVariables(
                    convective_flux_midpoint,
                    velocity_midpoint,
                    primitive_variables_minus,
                    primitive_variables_plus,
                    DIRECTION::X_DIRECTION,
                    RIEMANN_SOLVER::HLLC,
                    domains[0]);
            }
            else
            {
                riemann_solver->computeConvectiveFluxFromPrimitiveVariables(
                    convective_flux_midpoint,
                    primitive_variables_minus,
                    primitive_variables_plus,
                    DIRECTION::X_DIRECTION,
                    RIEMANN_SOLVER::HLLC,
                    domains[0]);
            }
            
            /*
             * Compute mid-point flux in the y-direction.
             */
            
            if (d_has_advective_eqn_form)
            {
                riemann_solver->computeConvectiveFluxAndVelocityFromPrimitiveVariables(
                    convective_flux_midpoint,
                    velocity_midpoint,
                    primitive_variables_minus,
                    primitive_variables_plus,
                    DIRECTION::Y_DIRECTION,
                    RIEMANN_SOLVER::HLLC,
                    domains[1]);
            }
            else
            {
                riemann_solver->computeConvectiveFluxFromPrimitiveVariables(
                    convective_flux_midpoint,
                    primitive_variables_minus,
                    primitive_variables_plus,
                    DIRECTION::Y_DIRECTION,
                    RIEMANN_SOLVER::HLLC,
                    domains[1]);
            }
            
            /*
             * Reconstruct the flux in the x-direction.
             */
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                Real* F_face_x = convective_flux->getPointer(0, ei);
                
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                    {
                        // Compute the linear indices.
                        const int idx_face_x = i +
                            j*(interior_dim_0 + 1);
                        
                        if (s_x[idx_face_x] > Real(0))
                        {
                            const int idx_midpoint_x = i +
                                j*(interior_dim_0 + 1);
                            
                            const int idx_node_LLL = (i - 3 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                            
                            const int idx_node_LL = (i - 2 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                            
                            const int idx_node_L = (i - 1 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                            
                            const int idx_node_R = (i + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                            
                            const int idx_node_RR = (i + 1 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                            
                            const int idx_node_RRR = (i + 2 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                            
                            F_face_x[idx_face_x] = Real(dt)*(
                                a_midpoint_r*F_midpoint_x[ei][idx_midpoint_x] +
                                a_node_r*(F_node_x[ei][idx_node_L]   + F_node_x[ei][idx_node_R]) +
                                b_node_r*(F_node_x[ei][idx_node_LL]  + F_node_x[ei][idx_node_RR]) +
                                c_node_r*(F_node_x[ei][idx_node_LLL] + F_node_x[ei][idx_node_RRR])
                                );
                        }
                    }
                }
            }
            
            /*
             * Reconstruct the flux in the y-direction.
             */
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                Real* F_face_y = convective_flux->getPointer(1, ei);
                
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_face_y = i +
                            j*interior_dim_0;
                        
                        if (s_y[idx_face_y] > Real(0))
                        {
                            const int idx_midpoint_y = i +
                                j*interior_dim_0;
                            
                            const int idx_node_BBB = (i + num_subghosts_0_convective_flux_y) +
                                (j - 3 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                            
                            const int idx_node_BB = (i + num_subghosts_0_convective_flux_y) +
                                (j - 2 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                            
                            const int idx_node_B = (i + num_subghosts_0_convective_flux_y) +
                                (j - 1 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                            
                            const int idx_node_T = (i + num_subghosts_0_convective_flux_y) +
                                (j + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                            
                            const int idx_node_TT = (i + num_subghosts_0_convective_flux_y) +
                                (j + 1 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                            
                            const int idx_node_TTT = (i + num_subghosts_0_convective_flux_y) +
                                (j + 2 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                            
                            F_face_y[idx_face_y] = Real(dt)*(
                                a_midpoint_r*F_midpoint_y[ei][idx_midpoint_y] +
                                a_node_r*(F_node_y[ei][idx_node_B]   + F_node_y[ei][idx_node_T]) +
                                b_node_r*(F_node_y[ei][idx_node_BB]  + F_node_y[ei][idx_node_TT]) +
                                c_node_r*(F_node_y[ei][idx_node_BBB] + F_node_y[ei][idx_node_TTT])
                                );
                        }
                    }
                }
            }
            
            /*
             * Compute the source.
             */
            
            if (d_has_advective_eqn_form)
            {
                Real* u_midpoint_x = velocity_midpoint->getPointer(0, 0);
                Real* v_midpoint_y = velocity_midpoint->getPointer(1, 1);
                
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    if (d_eqn_form[ei] == EQN_FORM::ADVECTIVE)
                    {
                        Real* S = source_scratch->getPointer(ei);
                        
                        const int num_subghosts_0_conservative_var = num_subghosts_conservative_var[ei][0];
                        const int num_subghosts_1_conservative_var = num_subghosts_conservative_var[ei][1];
                        const int subghostcell_dim_0_conservative_var = subghostcell_dims_conservative_var[ei][0];
                        
                        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                        {
                            HAMERS_PRAGMA_SIMD
                            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                            {
                                // Compute the linear indices.
                                const int idx_cell_nghost = i + j*interior_dim_0;
                                
                                if (s[idx_cell_nghost] > Real(0))
                                {
                                    const int idx_cell_wghost = (i + num_subghosts_0_conservative_var) +
                                        (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                                    
                                    const int idx_cell_wghost_x_LLL = (i - 3 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                    
                                    const int idx_cell_wghost_x_LL = (i - 2 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                    
                                    const int idx_cell_wghost_x_L = (i - 1 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                    
                                    const int idx_cell_wghost_x_R = (i + 1 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                    
                                    const int idx_cell_wghost_x_RR = (i + 2 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                    
                                    const int idx_cell_wghost_x_RRR = (i + 3 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                    
                                    const int idx_cell_wghost_y_BBB = (i + num_subghosts_0_velocity) +
                                        (j - 3 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                    
                                    const int idx_cell_wghost_y_BB = (i + num_subghosts_0_velocity) +
                                        (j - 2 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                    
                                    const int idx_cell_wghost_y_B = (i + num_subghosts_0_velocity) +
                                        (j - 1 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                    
                                    const int idx_cell_wghost_y_T = (i + num_subghosts_0_velocity) +
                                        (j + 1 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                    
                                    const int idx_cell_wghost_y_TT = (i + num_subghosts_0_velocity) +
                                        (j + 2 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                    
                                    const int idx_cell_wghost_y_TTT = (i + num_subghosts_0_velocity) +
                                        (j + 3 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                    
                                    const int idx_midpoint_x_L = i +
                                        j*(interior_dim_0 + 1);
                                    
                                    const int idx_midpoint_x_R = (i + 1) +
                                        j*(interior_dim_0 + 1);
                                    
                                    const int idx_midpoint_y_B = i +
                                        j*interior_dim_0;
                                    
                                    const int idx_midpoint_y_T = i +
                                        (j + 1)*interior_dim_0;
                                    
                                    S[idx_cell_nghost] = Real(dt)*Q[ei][idx_cell_wghost]*((
                                        a_midpoint*(u_midpoint_x[idx_midpoint_x_R] - u_midpoint_x[idx_midpoint_x_L]) +
                                        a_node*(u[idx_cell_wghost_x_R]   - u[idx_cell_wghost_x_L]) +
                                        b_node*(u[idx_cell_wghost_x_RR]  - u[idx_cell_wghost_x_LL]) +
                                        c_node*(u[idx_cell_wghost_x_RRR] - u[idx_cell_wghost_x_LLL])
                                        )/Real(dx[0]) + (
                                        a_midpoint*(v_midpoint_y[idx_midpoint_y_T] - v_midpoint_y[idx_midpoint_y_B]) +
                                        a_node*(v[idx_cell_wghost_y_T]   - v[idx_cell_wghost_y_B]) +
                                        b_node*(v[idx_cell_wghost_y_TT]  - v[idx_cell_wghost_y_BB]) +
                                        c_node*(v[idx_cell_wghost_y_TTT] - v[idx_cell_wghost_y_BBB])
                                        )/Real(dx[1]));
                                }
                            }
                        }
                    }
                }
            }
        }
        else if (d_dim == tbox::Dimension(3))
        {
            /*
             * Get the local lower indices and the number of cells in each dimension.
             */
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_lo_1 = domain_lo[1];
            const int domain_lo_2 = domain_lo[2];
            const int domain_dim_0 = domain_dims[0];
            const int domain_dim_1 = domain_dims[1];
            const int domain_dim_2 = domain_dims[2];
            
            /*
             * Get the interior dimensions.
             */
            
            const int interior_dim_0 = interior_dims[0];
            const int interior_dim_1 = interior_dims[1];
            
            /*
             * Get the pointers to the velocity and convective flux cell data inside the flow model.
             * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
             */
            
            HAMERS_SHARED_PTR<pdat::CellData<Real> > velocity = d_flow_model->getCellData("VELOCITY");
            
            std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > convective_flux_node(3);
            convective_flux_node[0] = d_flow_model->getCellData("CONVECTIVE_FLUX_X");
            convective_flux_node[1] = d_flow_model->getCellData("CONVECTIVE_FLUX_Y");
            convective_flux_node[2] = d_flow_model->getCellData("CONVECTIVE_FLUX_Z");
            
            hier::IntVector num_subghosts_velocity = velocity->getGhostCellWidth();
            hier::IntVector subghostcell_dims_velocity = velocity->getGhostBox().numberCells();
            
            hier::IntVector num_subghosts_convective_flux_x = convective_flux_node[0]->getGhostCellWidth();
            hier::IntVector subghostcell_dims_convective_flux_x = convective_flux_node[0]->getGhostBox().numberCells();
            
            hier::IntVector num_subghosts_convective_flux_y = convective_flux_node[1]->getGhostCellWidth();
            hier::IntVector subghostcell_dims_convective_flux_y = convective_flux_node[1]->getGhostBox().numberCells();
            
            hier::IntVector num_subghosts_convective_flux_z = convective_flux_node[2]->getGhostCellWidth();
            hier::IntVector subghostcell_dims_convective_flux_z = convective_flux_node[2]->getGhostBox().numberCells();
            
            const int num_subghosts_0_velocity = num_subghosts_velocity[0];
            const int num_subghosts_1_velocity = num_subghosts_velocity[1];
            const int num_subghosts_2_velocity = num_subghosts_velocity[2];
            const int subghostcell_dim_0_velocity = subghostcell_dims_velocity[0];
            const int subghostcell_dim_1_velocity = subghostcell_dims_velocity[1];
            
            const int num_subghosts_0_convective_flux_x = num_subghosts_convective_flux_x[0];
            const int num_subghosts_1_convective_flux_x = num_subghosts_convective_flux_x[1];
            const int num_subghosts_2_convective_flux_x = num_subghosts_convective_flux_x[2];
            const int subghostcell_dim_0_convective_flux_x = subghostcell_dims_convective_flux_x[0];
            const int subghostcell_dim_1_convective_flux_x = subghostcell_dims_convective_flux_x[1];
            
            const int num_subghosts_0_convective_flux_y = num_subghosts_convective_flux_y[0];
            const int num_subghosts_1_convective_flux_y = num_subghosts_convective_flux_y[1];
            const int num_subghosts_2_convective_flux_y = num_subghosts_convective_flux_y[2];
            const int subghostcell_dim_0_convective_flux_y = subghostcell_dims_convective_flux_y[0];
            const int subghostcell_dim_1_convective_flux_y = subghostcell_dims_convective_flux_y[1];
            
            const int num_subghosts_0_convective_flux_z = num_subghosts_convective_flux_z[0];
            const int num_subghosts_1_convective_flux_z = num_subghosts_convective_flux_z[1];
            const int num_subghosts_2_convective_flux_z = num_subghosts_convective_flux_z[2];
            const int subghostcell_dim_0_convective_flux_z = subghostcell_dims_convective_flux_z[0];
            const int subghostcell_dim_1_convective_flux_z = subghostcell_dims_convective_flux_z[1];
            
            Real* u = velocity->getPointer(0);
            Real* v = velocity->getPointer(1);
            Real* w = velocity->getPointer(2);
            
            std::vector<Real*> F_node_x;
            std::vector<Real*> F_node_y;
            std::vector<Real*> F_node_z;
            F_node_x.reserve(d_num_eqn);
            F_node_y.reserve(d_num_eqn);
            F_node_z.reserve(d_num_eqn);
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_node_x.push_back(convective_flux_node[0]->getPointer(ei));
                F_node_y.push_back(convective_flux_node[1]->getPointer(ei));
                F_node_z.push_back(convective_flux_node[2]->getPointer(ei));
            }
            
            std::vector<Real*> F_midpoint_x;
            std::vector<Real*> F_midpoint_y;
            std::vector<Real*> F_midpoint_z;
            F_midpoint_x.reserve(d_num_eqn);
            F_midpoint_y.reserve(d_num_eqn);
            F_midpoint_z.reserve(d_num_eqn);
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_midpoint_x.push_back(convective_flux_midpoint->getPointer(0, ei));
                F_midpoint_y.push_back(convective_flux_midpoint->getPointer(1, ei));
                F_midpoint_z.push_back(convective_flux_midpoint->getPointer(2, ei));
            }
            
            /*
             * Use first order interpolation if interpolated side primitive variables in x-direction
             * are out of bounds.
             */
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                V_minus[ei] = primitive_variables_minus[ei]->getPointer(0);
                V_plus[ei] = primitive_variables_plus[ei]->getPointer(0);
            }
            
            flag_minus = bounded_flag_minus->getPointer(0);
            flag_plus = bounded_flag_plus->getPointer(0);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                const int num_subghosts_0_primitive_var = num_subghosts_primitive_var[ei][0];
                const int num_subghosts_1_primitive_var = num_subghosts_primitive_var[ei][1];
                const int num_subghosts_2_primitive_var = num_subghosts_primitive_var[ei][2];
                const int subghostcell_dim_0_primitive_var = subghostcell_dims_primitive_var[ei][0];
                const int subghostcell_dim_1_primitive_var = subghostcell_dims_primitive_var[ei][1];
                
                for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                {
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                        {
                            // Compute the linear indices.
                            const int idx_midpoint_x = i +
                                j*(interior_dim_0 + 1) +
                                k*(interior_dim_0 + 1)*
                                    interior_dim_1;
                            
                            const int idx_cell_L = (i - 1 + num_subghosts_0_primitive_var) +
                                (j + num_subghosts_1_primitive_var)*subghostcell_dim_0_primitive_var +
                                (k + num_subghosts_2_primitive_var)*subghostcell_dim_0_primitive_var*
                                    subghostcell_dim_1_primitive_var;
                            
                            const int idx_cell_R = (i + num_subghosts_0_primitive_var) +
                                (j + num_subghosts_1_primitive_var)*subghostcell_dim_0_primitive_var +
                                (k + num_subghosts_2_primitive_var)*subghostcell_dim_0_primitive_var*
                                    subghostcell_dim_1_primitive_var;
                            
                            if (flag_minus[idx_midpoint_x] == 0 || flag_plus[idx_midpoint_x] == 0)
                            {
                                V_minus[ei][idx_midpoint_x] = V[ei][idx_cell_L];
                                V_plus[ei][idx_midpoint_x] = V[ei][idx_cell_R];
                            }
                        }
                    }
                }
            }
            
            /*
             * Use first order interpolation if interpolated side primitive variables in y-direction
             * are out of bounds.
             */
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                V_minus[ei] = primitive_variables_minus[ei]->getPointer(1);
                V_plus[ei] = primitive_variables_plus[ei]->getPointer(1);
            }
            
            flag_minus = bounded_flag_minus->getPointer(1);
            flag_plus = bounded_flag_plus->getPointer(1);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                const int num_subghosts_0_primitive_var = num_subghosts_primitive_var[ei][0];
                const int num_subghosts_1_primitive_var = num_subghosts_primitive_var[ei][1];
                const int num_subghosts_2_primitive_var = num_subghosts_primitive_var[ei][2];
                const int subghostcell_dim_0_primitive_var = subghostcell_dims_primitive_var[ei][0];
                const int subghostcell_dim_1_primitive_var = subghostcell_dims_primitive_var[ei][1];
                
                for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                {
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_midpoint_y = i +
                                j*interior_dim_0 +
                                k*interior_dim_0*
                                    (interior_dim_1 + 1);
                            
                            const int idx_cell_B = (i + num_subghosts_0_primitive_var) +
                                (j - 1 + num_subghosts_1_primitive_var)*subghostcell_dim_0_primitive_var +
                                (k + num_subghosts_2_primitive_var)*subghostcell_dim_0_primitive_var*
                                    subghostcell_dim_1_primitive_var;
                            
                            const int idx_cell_T = (i + num_subghosts_0_primitive_var) +
                                (j + num_subghosts_1_primitive_var)*subghostcell_dim_0_primitive_var +
                                (k + num_subghosts_2_primitive_var)*subghostcell_dim_0_primitive_var*
                                    subghostcell_dim_1_primitive_var;
                            
                            if (flag_minus[idx_midpoint_y] == 0 || flag_plus[idx_midpoint_y] == 0)
                            {
                                V_minus[ei][idx_midpoint_y] = V[ei][idx_cell_B];
                                V_plus[ei][idx_midpoint_y] = V[ei][idx_cell_T];
                            }
                        }
                    }
                }
            }
            
            /*
             * Use first order interpolation if interpolated side primitive variables in z-direction
             * are out of bounds.
             */
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                V_minus[ei] = primitive_variables_minus[ei]->getPointer(2);
                V_plus[ei] = primitive_variables_plus[ei]->getPointer(2);
            }
            
            flag_minus = bounded_flag_minus->getPointer(2);
            flag_plus = bounded_flag_plus->getPointer(2);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                const int num_subghosts_0_primitive_var = num_subghosts_primitive_var[ei][0];
                const int num_subghosts_1_primitive_var = num_subghosts_primitive_var[ei][1];
                const int num_subghosts_2_primitive_var = num_subghosts_primitive_var[ei][2];
                const int subghostcell_dim_0_primitive_var = subghostcell_dims_primitive_var[ei][0];
                const int subghostcell_dim_1_primitive_var = subghostcell_dims_primitive_var[ei][1];
                
                for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2 + 1; k++)
                {
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_midpoint_z = i +
                                j*interior_dim_0 +
                                k*interior_dim_0*
                                    interior_dim_1;
                            
                            const int idx_cell_B = (i + num_subghosts_0_primitive_var) +
                                (j + num_subghosts_1_primitive_var)*subghostcell_dim_0_primitive_var +
                                (k - 1 + num_subghosts_2_primitive_var)*subghostcell_dim_0_primitive_var*
                                    subghostcell_dim_1_primitive_var;
                            
                            const int idx_cell_F = (i + num_subghosts_0_primitive_var) +
                                (j + num_subghosts_1_primitive_var)*subghostcell_dim_0_primitive_var +
                                (k + num_subghosts_2_primitive_var)*subghostcell_dim_0_primitive_var*
                                    subghostcell_dim_1_primitive_var;
                            
                            if (flag_minus[idx_midpoint_z] == 0 || flag_plus[idx_midpoint_z] == 0)
                            {
                                V_minus[ei][idx_midpoint_z] = V[ei][idx_cell_B];
                                V_plus[ei][idx_midpoint_z] = V[ei][idx_cell_F];
                            }
                        }
                    }
                }
            }
            
            /*
             * Compute mid-point flux in the x-direction.
             */
            
            if (d_has_advective_eqn_form)
            {
                riemann_solver->computeConvectiveFluxAndVelocityFromPrimitiveVariables(
                    convective_flux_midpoint,
                    velocity_midpoint,
                    primitive_variables_minus,
                    primitive_variables_plus,
                    DIRECTION::X_DIRECTION,
                    RIEMANN_SOLVER::HLLC);
            }
            else
            {
                riemann_solver->computeConvectiveFluxFromPrimitiveVariables(
                    convective_flux_midpoint,
                    primitive_variables_minus,
                    primitive_variables_plus,
                    DIRECTION::X_DIRECTION,
                    RIEMANN_SOLVER::HLLC);
            }
            
            /*
             * Compute mid-point flux in the y-direction.
             */
            
            if (d_has_advective_eqn_form)
            {
                riemann_solver->computeConvectiveFluxAndVelocityFromPrimitiveVariables(
                    convective_flux_midpoint,
                    velocity_midpoint,
                    primitive_variables_minus,
                    primitive_variables_plus,
                    DIRECTION::Y_DIRECTION,
                    RIEMANN_SOLVER::HLLC);
            }
            else
            {
                riemann_solver->computeConvectiveFluxFromPrimitiveVariables(
                    convective_flux_midpoint,
                    primitive_variables_minus,
                    primitive_variables_plus,
                    DIRECTION::Y_DIRECTION,
                    RIEMANN_SOLVER::HLLC);
            }
            
            /*
             * Compute mid-point flux in the z-direction.
             */
            
            if (d_has_advective_eqn_form)
            {
                riemann_solver->computeConvectiveFluxAndVelocityFromPrimitiveVariables(
                    convective_flux_midpoint,
                    velocity_midpoint,
                    primitive_variables_minus,
                    primitive_variables_plus,
                    DIRECTION::Z_DIRECTION,
                    RIEMANN_SOLVER::HLLC);
            }
            else
            {
                riemann_solver->computeConvectiveFluxFromPrimitiveVariables(
                    convective_flux_midpoint,
                    primitive_variables_minus,
                    primitive_variables_plus,
                    DIRECTION::Z_DIRECTION,
                    RIEMANN_SOLVER::HLLC);
            }
            
            /*
             * Reconstruct the flux in the x-direction.
             */
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                Real* F_face_x = convective_flux->getPointer(0, ei);
                
                for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                {
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                        {
                            // Compute the linear indices.
                            const int idx_face_x = i +
                                j*(interior_dim_0 + 1) +
                                k*(interior_dim_0 + 1)*
                                    interior_dim_1;
                            
                            if (s_x[idx_face_x] > Real(0))
                            {
                                const int idx_midpoint_x = i +
                                    j*(interior_dim_0 + 1) +
                                    k*(interior_dim_0 + 1)*
                                        interior_dim_1;
                                
                                const int idx_node_LLL = (i - 3 + num_subghosts_0_convective_flux_x) +
                                    (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                    (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                        subghostcell_dim_1_convective_flux_x;
                                
                                const int idx_node_LL = (i - 2 + num_subghosts_0_convective_flux_x) +
                                    (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                    (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                        subghostcell_dim_1_convective_flux_x;
                                
                                const int idx_node_L = (i - 1 + num_subghosts_0_convective_flux_x) +
                                    (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                    (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                        subghostcell_dim_1_convective_flux_x;
                                
                                const int idx_node_R = (i + num_subghosts_0_convective_flux_x) +
                                    (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                    (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                        subghostcell_dim_1_convective_flux_x;
                                
                                const int idx_node_RR = (i + 1 + num_subghosts_0_convective_flux_x) +
                                    (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                    (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                        subghostcell_dim_1_convective_flux_x;
                                
                                const int idx_node_RRR = (i + 2 + num_subghosts_0_convective_flux_x) +
                                    (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                    (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                        subghostcell_dim_1_convective_flux_x;
                                
                                F_face_x[idx_face_x] = Real(dt)*(
                                    a_midpoint_r*F_midpoint_x[ei][idx_midpoint_x] +
                                    a_node_r*(F_node_x[ei][idx_node_L]   + F_node_x[ei][idx_node_R]) +
                                    b_node_r*(F_node_x[ei][idx_node_LL]  + F_node_x[ei][idx_node_RR]) +
                                    c_node_r*(F_node_x[ei][idx_node_LLL] + F_node_x[ei][idx_node_RRR])
                                    );
                            }
                        }
                    }
                }
            }
            
            /*
             * Reconstruct the flux in the y-direction.
             */
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                Real* F_face_y = convective_flux->getPointer(1, ei);
                
                for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                {
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_face_y = i +
                                j*interior_dim_0 +
                                k*interior_dim_0*(interior_dim_1 + 1);
                            
                            if (s_y[idx_face_y] > Real(0))
                            {
                                const int idx_midpoint_y = i +
                                    j*interior_dim_0 +
                                    k*interior_dim_0*(interior_dim_1 + 1);
                                
                                const int idx_node_BBB = (i + num_subghosts_0_convective_flux_y) +
                                    (j - 3 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                    (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                        subghostcell_dim_1_convective_flux_y;
                                
                                const int idx_node_BB = (i + num_subghosts_0_convective_flux_y) +
                                    (j - 2 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                    (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                        subghostcell_dim_1_convective_flux_y;
                                
                                const int idx_node_B = (i + num_subghosts_0_convective_flux_y) +
                                    (j - 1 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                    (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                        subghostcell_dim_1_convective_flux_y;
                                
                                const int idx_node_T = (i + num_subghosts_0_convective_flux_y) +
                                    (j + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                    (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                        subghostcell_dim_1_convective_flux_y;
                                
                                const int idx_node_TT = (i + num_subghosts_0_convective_flux_y) +
                                    (j + 1 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                    (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                        subghostcell_dim_1_convective_flux_y;
                                
                                const int idx_node_TTT = (i + num_subghosts_0_convective_flux_y) +
                                    (j + 2 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                    (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                        subghostcell_dim_1_convective_flux_y;
                                
                                F_face_y[idx_face_y] = Real(dt)*(
                                    a_midpoint_r*F_midpoint_y[ei][idx_midpoint_y] +
                                    a_node_r*(F_node_y[ei][idx_node_B]   + F_node_y[ei][idx_node_T]) +
                                    b_node_r*(F_node_y[ei][idx_node_BB]  + F_node_y[ei][idx_node_TT]) +
                                    c_node_r*(F_node_y[ei][idx_node_BBB] + F_node_y[ei][idx_node_TTT])
                                    );
                            }
                        }
                    }
                }
            }
            
            
            /*
             * Reconstruct the flux in the z-direction.
             */
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                Real* F_face_z = convective_flux->getPointer(2, ei);
                
                for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2 + 1; k++)
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
                            
                            if (s_y[idx_face_z] > Real(0))
                            {
                                const int idx_midpoint_z =  i +
                                    j*interior_dim_0 +
                                    k*interior_dim_0*interior_dim_1;
                                
                                const int idx_node_BBB = (i + num_subghosts_0_convective_flux_z) +
                                    (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                    (k - 3 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                        subghostcell_dim_1_convective_flux_z;
                                
                                const int idx_node_BB = (i + num_subghosts_0_convective_flux_z) +
                                    (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                    (k - 2 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                        subghostcell_dim_1_convective_flux_z;
                                
                                const int idx_node_B = (i + num_subghosts_0_convective_flux_z) +
                                    (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                    (k - 1 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                        subghostcell_dim_1_convective_flux_z;
                                
                                const int idx_node_F = (i + num_subghosts_0_convective_flux_z) +
                                    (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                    (k + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                        subghostcell_dim_1_convective_flux_z;
                                
                                const int idx_node_FF = (i + num_subghosts_0_convective_flux_z) +
                                    (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                    (k + 1 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                        subghostcell_dim_1_convective_flux_z;
                                
                                const int idx_node_FFF = (i + num_subghosts_0_convective_flux_z) +
                                    (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                    (k + 2 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                        subghostcell_dim_1_convective_flux_z;
                                
                                F_face_z[idx_face_z] = Real(dt)*(
                                    a_midpoint_r*F_midpoint_z[ei][idx_midpoint_z] +
                                    a_node_r*(F_node_z[ei][idx_node_B]   + F_node_z[ei][idx_node_F]) +
                                    b_node_r*(F_node_z[ei][idx_node_BB]  + F_node_z[ei][idx_node_FF]) +
                                    c_node_r*(F_node_z[ei][idx_node_BBB] + F_node_z[ei][idx_node_FFF])
                                    );
                            }
                        }
                    }
                }
            }
            
            /*
             * Compute the source.
             */
            
            if (d_has_advective_eqn_form)
            {
                Real* u_midpoint_x = velocity_midpoint->getPointer(0, 0);
                Real* v_midpoint_y = velocity_midpoint->getPointer(1, 1);
                Real* w_midpoint_z = velocity_midpoint->getPointer(2, 2);
                
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    if (d_eqn_form[ei] == EQN_FORM::ADVECTIVE)
                    {
                        Real* S = source_scratch->getPointer(ei);
                        
                        const int num_subghosts_0_conservative_var = num_subghosts_conservative_var[ei][0];
                        const int num_subghosts_1_conservative_var = num_subghosts_conservative_var[ei][1];
                        const int num_subghosts_2_conservative_var = num_subghosts_conservative_var[ei][2];
                        const int subghostcell_dim_0_conservative_var = subghostcell_dims_conservative_var[ei][0];
                        const int subghostcell_dim_1_conservative_var = subghostcell_dims_conservative_var[ei][1];
                        
                        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                        {
                            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                            {
                                HAMERS_PRAGMA_SIMD
                                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                                {
                                    // Compute the linear indices.
                                    const int idx_cell_nghost = i +
                                        j*interior_dim_0 +
                                        k*interior_dim_0*
                                            interior_dim_1;
                                    
                                    if (s[idx_cell_nghost] > Real(0))
                                    {
                                        const int idx_cell_wghost = (i + num_subghosts_0_conservative_var) +
                                            (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                            (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                                subghostcell_dim_1_conservative_var;
                                        
                                        const int idx_cell_wghost_x_LLL = (i - 3 + num_subghosts_0_velocity) +
                                            (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                            (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                                subghostcell_dim_1_velocity;
                                        
                                        const int idx_cell_wghost_x_LL = (i - 2 + num_subghosts_0_velocity) +
                                            (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                            (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                                subghostcell_dim_1_velocity;
                                        
                                        const int idx_cell_wghost_x_L = (i - 1 + num_subghosts_0_velocity) +
                                            (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                            (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                                subghostcell_dim_1_velocity;
                                        
                                        const int idx_cell_wghost_x_R = (i + 1 + num_subghosts_0_velocity) +
                                            (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                            (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                                subghostcell_dim_1_velocity;
                                        
                                        const int idx_cell_wghost_x_RR = (i + 2 + num_subghosts_0_velocity) +
                                            (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                            (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                                subghostcell_dim_1_velocity;
                                        
                                        const int idx_cell_wghost_x_RRR = (i + 3 + num_subghosts_0_velocity) +
                                            (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                            (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                                subghostcell_dim_1_velocity;
                                        
                                        const int idx_cell_wghost_y_BBB = (i + num_subghosts_0_velocity) +
                                            (j - 3 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                            (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                                subghostcell_dim_1_velocity;
                                        
                                        const int idx_cell_wghost_y_BB = (i + num_subghosts_0_velocity) +
                                            (j - 2 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                            (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                                subghostcell_dim_1_velocity;
                                        
                                        const int idx_cell_wghost_y_B = (i + num_subghosts_0_velocity) +
                                            (j - 1 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                            (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                                subghostcell_dim_1_velocity;
                                        
                                        const int idx_cell_wghost_y_T = (i + num_subghosts_0_velocity) +
                                            (j + 1 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                            (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                                subghostcell_dim_1_velocity;
                                        
                                        const int idx_cell_wghost_y_TT = (i + num_subghosts_0_velocity) +
                                            (j + 2 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                            (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                                subghostcell_dim_1_velocity;
                                        
                                        const int idx_cell_wghost_y_TTT = (i + num_subghosts_0_velocity) +
                                            (j + 3 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                            (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                                subghostcell_dim_1_velocity;
                                        
                                        const int idx_cell_wghost_z_BBB = (i + num_subghosts_0_velocity) +
                                            (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                            (k - 3 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                                subghostcell_dim_1_velocity;
                                        
                                        const int idx_cell_wghost_z_BB = (i + num_subghosts_0_velocity) +
                                            (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                            (k - 2 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                                subghostcell_dim_1_velocity;
                                        
                                        const int idx_cell_wghost_z_B = (i + num_subghosts_0_velocity) +
                                            (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                            (k - 1 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                                subghostcell_dim_1_velocity;
                                        
                                        const int idx_cell_wghost_z_F = (i + num_subghosts_0_velocity) +
                                            (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                            (k + 1 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                                subghostcell_dim_1_velocity;
                                        
                                        const int idx_cell_wghost_z_FF = (i + num_subghosts_0_velocity) +
                                            (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                            (k + 2 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                                subghostcell_dim_1_velocity;
                                        
                                        const int idx_cell_wghost_z_FFF = (i + num_subghosts_0_velocity) +
                                            (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                            (k + 3 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                                subghostcell_dim_1_velocity;
                                        
                                        const int idx_midpoint_x_L = i +
                                            j*(interior_dim_0 + 1) +
                                            k*(interior_dim_0 + 1)*
                                                interior_dim_1;
                                        
                                        const int idx_midpoint_x_R = (i + 1) +
                                            j*(interior_dim_0 + 1) +
                                            k*(interior_dim_0 + 1)*
                                                interior_dim_1;
                                        
                                        const int idx_midpoint_y_B = i +
                                            j*interior_dim_0 +
                                            k*interior_dim_0*
                                                (interior_dim_1 + 1);
                                        
                                        const int idx_midpoint_y_T = i +
                                            (j + 1)*interior_dim_0 +
                                            k*interior_dim_0*
                                                (interior_dim_1 + 1);
                                        
                                        const int idx_midpoint_z_B = i +
                                            j*interior_dim_0 +
                                            k*interior_dim_0*
                                                interior_dim_1;
                                        
                                        const int idx_midpoint_z_F = i +
                                            j*interior_dim_0 +
                                            (k + 1)*interior_dim_0*
                                                interior_dim_1;
                                        
                                        S[idx_cell_nghost] = Real(dt)*Q[ei][idx_cell_wghost]*((
                                            a_midpoint*(u_midpoint_x[idx_midpoint_x_R] - u_midpoint_x[idx_midpoint_x_L]) +
                                            a_node*(u[idx_cell_wghost_x_R]   - u[idx_cell_wghost_x_L]) +
                                            b_node*(u[idx_cell_wghost_x_RR]  - u[idx_cell_wghost_x_LL]) +
                                            c_node*(u[idx_cell_wghost_x_RRR] - u[idx_cell_wghost_x_LLL])
                                            )/Real(dx[0]) + (
                                            a_midpoint*(v_midpoint_y[idx_midpoint_y_T] - v_midpoint_y[idx_midpoint_y_B]) +
                                            a_node*(v[idx_cell_wghost_y_T]   - v[idx_cell_wghost_y_B]) +
                                            b_node*(v[idx_cell_wghost_y_TT]  - v[idx_cell_wghost_y_BB]) +
                                            c_node*(v[idx_cell_wghost_y_TTT] - v[idx_cell_wghost_y_BBB])
                                            )/Real(dx[1]) + (
                                            a_midpoint*(w_midpoint_z[idx_midpoint_z_F] - w_midpoint_z[idx_midpoint_z_B]) +
                                            a_node*(w[idx_cell_wghost_z_F]   - w[idx_cell_wghost_z_B]) +
                                            b_node*(w[idx_cell_wghost_z_FF]  - w[idx_cell_wghost_z_BB]) +
                                            c_node*(w[idx_cell_wghost_z_FFF] - w[idx_cell_wghost_z_BBB])
                                            )/Real(dx[2]));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    /*
     * Unregister the patch and data of all registered derived cell variables in the flow model.
     */
    
    d_flow_model->unregisterPatch();
}


/*
 * Perform WENO interpolation.
 */
void
ConvectiveFluxReconstructor::performWENOInterpolation(
    std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& variables_minus,
    std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& variables_plus,
    const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > >& variables,
    const std::vector<hier::Box>& domains) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(variables_minus.size()) == d_num_eqn);
    TBOX_ASSERT(static_cast<int>(variables_plus.size()) == d_num_eqn);
    
    TBOX_ASSERT(static_cast<int>(variables.size()) == 6);
#endif
    
    // Default scheme constants.
    const int constant_p = 2;
    const int constant_q = 4;
    const Real constant_C = Real(1.0e9);
    const Real constant_alpha_tau = Real(35);
    
    hier::Box interior_box = variables_minus[0]->getBox();
    
    const hier::IntVector num_ghosts = variables_minus[0]->getGhostCellWidth();
    const hier::IntVector ghostcell_dims = variables_minus[0]->getGhostBox().numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    const hier::IntVector interior_dims = variables_minus[0]->getBox().numberCells();
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        TBOX_ASSERT(variables_minus[ei]->getBox().numberCells() == interior_dims);
        TBOX_ASSERT(variables_plus[ei]->getBox().numberCells() == interior_dims);
        
        TBOX_ASSERT(variables_minus[ei]->getGhostCellWidth() == num_ghosts);
        TBOX_ASSERT(variables_plus[ei]->getGhostCellWidth() == num_ghosts);
    }
    
    TBOX_ASSERT(static_cast<int>(variables.size()) == 6);
    
    for (int m = 0; m < 6; m++)
    {
        TBOX_ASSERT(static_cast<int>(variables[m].size()) == d_num_eqn);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ASSERT(variables[m][ei]->getBox().numberCells() == interior_dims);
            TBOX_ASSERT(variables[m][ei]->getGhostCellWidth() == num_ghosts);
        }
    }
#endif
    
    /*
     * Get the domain dimensions.
     */
    
    hier::IntVector domain_x_lo(d_dim);
    hier::IntVector domain_y_lo(d_dim);
    hier::IntVector domain_z_lo(d_dim);
    hier::IntVector domain_x_dims(d_dim);
    hier::IntVector domain_y_dims(d_dim);
    hier::IntVector domain_z_dims(d_dim);
    
    domain_x_lo = domains[0].lower() - interior_box.lower();
    domain_x_dims = domains[0].numberCells();
    if (d_dim > tbox::Dimension(1))
    {
        domain_y_lo = domains[1].lower() - interior_box.lower();
        domain_y_dims = domains[1].numberCells();
    }
    if (d_dim > tbox::Dimension(2))
    {
        domain_z_lo = domains[2].lower() - interior_box.lower();
        domain_z_dims = domains[2].numberCells();
    }
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower index and the number of cells in each dimension.
         */
        
        const int domain_x_lo_0  = domain_x_lo[0];
        const int domain_x_dim_0 = domain_x_dims[0];
        
        const int num_ghosts_0 = num_ghosts[0];
        
        /*
         * Peform WENO interpolation in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<Real*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(0));
            }
            
            Real* U_L = variables_minus[ei]->getPointer(0);
            
            switch(d_weno_interp)
            {
                case WENO_INTERP::TYPE::WENO5Z:
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_x_lo_0; i < domain_x_lo_0 + domain_x_dim_0 + 1; i++)
                    {
                        // Compute the linear index of the mid-point.
                        const int idx_midpoint_x = i + num_ghosts_0;
                        
                        performLocalWENOInterpolationMinusZ(
                            U_L,
                            U_array.data(),
                            idx_midpoint_x,
                            constant_p);
                    }
                    break;
                }
                case WENO_INTERP::TYPE::WENO6LD:
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_x_lo_0; i < domain_x_lo_0 + domain_x_dim_0 + 1; i++)
                    {
                        // Compute the linear index of the mid-point.
                        const int idx_midpoint_x = i + num_ghosts_0;
                        
                        performLocalWENOInterpolationMinusLD(
                            U_L,
                            U_array.data(),
                            idx_midpoint_x,
                            constant_p,
                            constant_q,
                            constant_C,
                            constant_alpha_tau);
                    }
                    break;
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
            
            switch(d_weno_interp)
            {
                case WENO_INTERP::TYPE::WENO5Z:
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_x_lo_0; i < domain_x_lo_0 + domain_x_dim_0 + 1; i++)
                    {
                        // Compute the linear index of the mid-point.
                        const int idx_midpoint_x = i + num_ghosts_0;
                        
                        performLocalWENOInterpolationPlusZ(
                            U_R,
                            U_array.data(),
                            idx_midpoint_x,
                            constant_p);
                    }
                    break;
                }
                case WENO_INTERP::TYPE::WENO6LD:
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_x_lo_0; i < domain_x_lo_0 + domain_x_dim_0 + 1; i++)
                    {
                        // Compute the linear index of the mid-point.
                        const int idx_midpoint_x = i + num_ghosts_0;
                        
                        performLocalWENOInterpolationPlusLD(
                            U_R,
                            U_array.data(),
                            idx_midpoint_x,
                            constant_p,
                            constant_q,
                            constant_C,
                            constant_alpha_tau);
                    }
                    break;
                }
            }
        }
    } // if (d_dim == tbox::Dimension(1))
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the local lower indices and the number of cells in each dimension.
         */
        
        const int domain_x_lo_0  = domain_x_lo[0];
        const int domain_x_lo_1  = domain_x_lo[1];
        const int domain_x_dim_0 = domain_x_dims[0];
        const int domain_x_dim_1 = domain_x_dims[1];
        
        const int domain_y_lo_0  = domain_y_lo[0];
        const int domain_y_lo_1  = domain_y_lo[1];
        const int domain_y_dim_0 = domain_y_dims[0];
        const int domain_y_dim_1 = domain_y_dims[1];
        
        const int num_ghosts_0 = num_ghosts[0];
        const int num_ghosts_1 = num_ghosts[1];
        const int ghostcell_dim_0 = ghostcell_dims[0];
        
        /*
         * Peform WENO interpolation in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<Real*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(0));
            }
            
            Real* U_L = variables_minus[ei]->getPointer(0);
            
            switch(d_weno_interp)
            {
                case WENO_INTERP::TYPE::WENO5Z:
                {
                    for (int j = domain_x_lo_1; j < domain_x_lo_1 + domain_x_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = domain_x_lo_0; i < domain_x_lo_0 + domain_x_dim_0 + 1; i++)
                        {
                            // Compute the linear index of the mid-point.
                            const int idx_midpoint_x = (i + num_ghosts_0) +
                                (j + num_ghosts_1)*(ghostcell_dim_0 + 1);
                            
                            performLocalWENOInterpolationMinusZ(
                                U_L,
                                U_array.data(),
                                idx_midpoint_x,
                                constant_p);
                        }
                    }
                    break;
                }
                case WENO_INTERP::TYPE::WENO6LD:
                {
                    for (int j = domain_x_lo_1; j < domain_x_lo_1 + domain_x_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = domain_x_lo_0; i < domain_x_lo_0 + domain_x_dim_0 + 1; i++)
                        {
                            // Compute the linear index of the mid-point.
                            const int idx_midpoint_x = (i + num_ghosts_0) +
                                (j + num_ghosts_1)*(ghostcell_dim_0 + 1);
                            
                            performLocalWENOInterpolationMinusLD(
                                U_L,
                                U_array.data(),
                                idx_midpoint_x,
                                constant_p,
                                constant_q,
                                constant_C,
                                constant_alpha_tau);
                        }
                    }
                    break;
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
            
            switch(d_weno_interp)
            {
                case WENO_INTERP::TYPE::WENO5Z:
                {
                    for (int j = domain_x_lo_1; j < domain_x_lo_1 + domain_x_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = domain_x_lo_0; i < domain_x_lo_0 + domain_x_dim_0 + 1; i++)
                        {
                            // Compute the linear index of the mid-point.
                            const int idx_midpoint_x = (i + num_ghosts_0) +
                                (j + num_ghosts_1)*(ghostcell_dim_0 + 1);
                            
                            performLocalWENOInterpolationPlusZ(
                                U_R,
                                U_array.data(),
                                idx_midpoint_x,
                                constant_p);
                        }
                    }
                    break;
                }
                case WENO_INTERP::TYPE::WENO6LD:
                {
                    for (int j = domain_x_lo_1; j < domain_x_lo_1 + domain_x_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = domain_x_lo_0; i < domain_x_lo_0 + domain_x_dim_0 + 1; i++)
                        {
                            // Compute the linear index of the mid-point.
                            const int idx_midpoint_x = (i + num_ghosts_0) +
                                (j + num_ghosts_1)*(ghostcell_dim_0 + 1);
                            
                            performLocalWENOInterpolationPlusLD(
                                U_R,
                                U_array.data(),
                                idx_midpoint_x,
                                constant_p,
                                constant_q,
                                constant_C,
                                constant_alpha_tau);
                        }
                    }
                    break;
                }
            }
        }
        
        /*
         * Peform WENO interpolation in the y-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<Real*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(1));
            }
            
            Real* U_B = variables_minus[ei]->getPointer(1);
            
            switch(d_weno_interp)
            {
                case WENO_INTERP::TYPE::WENO5Z:
                {
                    for (int j = domain_y_lo_1; j < domain_y_lo_1 + domain_y_dim_1 + 1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = domain_y_lo_0; i < domain_y_lo_0 + domain_y_dim_0; i++)
                        {
                            // Compute the linear index of the mid-point.
                            const int idx_midpoint_y = (i + num_ghosts_0) +
                                (j + num_ghosts_1)*ghostcell_dim_0;
                            
                            performLocalWENOInterpolationMinusZ(
                                U_B,
                                U_array.data(),
                                idx_midpoint_y,
                                constant_p);
                        }
                    }
                    break;
                }
                case WENO_INTERP::TYPE::WENO6LD:
                {
                    for (int j = domain_y_lo_1; j < domain_y_lo_1 + domain_y_dim_1 + 1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = domain_y_lo_0; i < domain_y_lo_0 + domain_y_dim_0; i++)
                        {
                            // Compute the linear index of the mid-point.
                            const int idx_midpoint_y = (i + num_ghosts_0) +
                                (j + num_ghosts_1)*ghostcell_dim_0;
                            
                            performLocalWENOInterpolationMinusLD(
                                U_B,
                                U_array.data(),
                                idx_midpoint_y,
                                constant_p,
                                constant_q,
                                constant_C,
                                constant_alpha_tau);
                        }
                    }
                    break;
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
            
            switch(d_weno_interp)
            {
                case WENO_INTERP::TYPE::WENO5Z:
                {
                    for (int j = domain_y_lo_1; j < domain_y_lo_1 + domain_y_dim_1 + 1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = domain_y_lo_0; i < domain_y_lo_0 + domain_y_dim_0; i++)
                        {
                            // Compute the linear index of the mid-point.
                            const int idx_midpoint_y = (i + num_ghosts_0) +
                                (j + num_ghosts_1)*ghostcell_dim_0;
                            
                            performLocalWENOInterpolationPlusZ(
                                U_T,
                                U_array.data(),
                                idx_midpoint_y,
                                constant_p);
                        }
                    }
                    break;
                }
                case WENO_INTERP::TYPE::WENO6LD:
                {
                    for (int j = domain_y_lo_1; j < domain_y_lo_1 + domain_y_dim_1 + 1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = domain_y_lo_0; i < domain_y_lo_0 + domain_y_dim_0; i++)
                        {
                            // Compute the linear index of the mid-point.
                            const int idx_midpoint_y = (i + num_ghosts_0) +
                                (j + num_ghosts_1)*ghostcell_dim_0;
                            
                            performLocalWENOInterpolationPlusLD(
                                U_T,
                                U_array.data(),
                                idx_midpoint_y,
                                constant_p,
                                constant_q,
                                constant_C,
                                constant_alpha_tau);
                        }
                    }
                    break;
                }
            }
        }
    } // if (d_dim == tbox::Dimension(2))
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the local lower indices and the number of cells in each dimension.
         */
        
        const int domain_x_lo_0  = domain_x_lo[0];
        const int domain_x_lo_1  = domain_x_lo[1];
        const int domain_x_lo_2  = domain_x_lo[2];
        const int domain_x_dim_0 = domain_x_dims[0];
        const int domain_x_dim_1 = domain_x_dims[1];
        const int domain_x_dim_2 = domain_x_dims[2];
        
        const int domain_y_lo_0  = domain_y_lo[0];
        const int domain_y_lo_1  = domain_y_lo[1];
        const int domain_y_lo_2  = domain_y_lo[2];
        const int domain_y_dim_0 = domain_y_dims[0];
        const int domain_y_dim_1 = domain_y_dims[1];
        const int domain_y_dim_2 = domain_y_dims[2];
        
        const int domain_z_lo_0  = domain_z_lo[0];
        const int domain_z_lo_1  = domain_z_lo[1];
        const int domain_z_lo_2  = domain_z_lo[2];
        const int domain_z_dim_0 = domain_z_dims[0];
        const int domain_z_dim_1 = domain_z_dims[1];
        const int domain_z_dim_2 = domain_z_dims[2];
        
        const int num_ghosts_0 = num_ghosts[0];
        const int num_ghosts_1 = num_ghosts[1];
        const int num_ghosts_2 = num_ghosts[2];
        const int ghostcell_dim_0 = ghostcell_dims[0];
        const int ghostcell_dim_1 = ghostcell_dims[1];
        
        /*
         * Peform WENO interpolation in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<Real*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(0));
            }
            
            Real* U_L = variables_minus[ei]->getPointer(0);
            
            switch(d_weno_interp)
            {
                case WENO_INTERP::TYPE::WENO5Z:
                {
                    for (int k = domain_x_lo_2; k < domain_x_lo_2 + domain_x_dim_2; k++)
                    {
                        for (int j = domain_x_lo_1; j < domain_x_lo_1 + domain_x_dim_1; j++)
                        {
                            HAMERS_PRAGMA_SIMD
                            for (int i = domain_x_lo_0; i < domain_x_lo_0 + domain_x_dim_0 + 1; i++)
                            {
                                // Compute the linear index of the mid-point.
                                const int idx_midpoint_x = (i + num_ghosts_0) +
                                    (j + num_ghosts_1)*(ghostcell_dim_0 + 1) +
                                    (k + num_ghosts_2)*(ghostcell_dim_0 + 1)*
                                        ghostcell_dim_1;
                                
                                performLocalWENOInterpolationMinusZ(
                                    U_L,
                                    U_array.data(),
                                    idx_midpoint_x,
                                    constant_p);
                            }
                        }
                    }
                    break;
                }
                case WENO_INTERP::TYPE::WENO6LD:
                {
                    for (int k = domain_x_lo_2; k < domain_x_lo_2 + domain_x_dim_2; k++)
                    {
                        for (int j = domain_x_lo_1; j < domain_x_lo_1 + domain_x_dim_1; j++)
                        {
                            HAMERS_PRAGMA_SIMD
                            for (int i = domain_x_lo_0; i < domain_x_lo_0 + domain_x_dim_0 + 1; i++)
                            {
                                // Compute the linear index of the mid-point.
                                const int idx_midpoint_x = (i + num_ghosts_0) +
                                    (j + num_ghosts_1)*(ghostcell_dim_0 + 1) +
                                    (k + num_ghosts_2)*(ghostcell_dim_0 + 1)*
                                        ghostcell_dim_1;
                                
                                performLocalWENOInterpolationMinusLD(
                                    U_L,
                                    U_array.data(),
                                    idx_midpoint_x,
                                    constant_p,
                                    constant_q,
                                    constant_C,
                                    constant_alpha_tau);
                            }
                        }
                    }
                    break;
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
            
            switch(d_weno_interp)
            {
                case WENO_INTERP::TYPE::WENO5Z:
                {
                    for (int k = domain_x_lo_2; k < domain_x_lo_2 + domain_x_dim_2; k++)
                    {
                        for (int j = domain_x_lo_1; j < domain_x_lo_1 + domain_x_dim_1; j++)
                        {
                            HAMERS_PRAGMA_SIMD
                            for (int i = domain_x_lo_0; i < domain_x_lo_0 + domain_x_dim_0 + 1; i++)
                            {
                                // Compute the linear index of the mid-point.
                                const int idx_midpoint_x = (i + num_ghosts_0) +
                                    (j + num_ghosts_1)*(ghostcell_dim_0 + 1) +
                                    (k + num_ghosts_2)*(ghostcell_dim_0 + 1)*
                                        ghostcell_dim_1;
                                
                                performLocalWENOInterpolationPlusZ(
                                    U_R,
                                    U_array.data(),
                                    idx_midpoint_x,
                                    constant_p);
                            }
                        }
                    }
                    break;
                }
                case WENO_INTERP::TYPE::WENO6LD:
                {
                    for (int k = domain_x_lo_2; k < domain_x_lo_2 + domain_x_dim_2; k++)
                    {
                        for (int j = domain_x_lo_1; j < domain_x_lo_1 + domain_x_dim_1; j++)
                        {
                            HAMERS_PRAGMA_SIMD
                            for (int i = domain_x_lo_0; i < domain_x_lo_0 + domain_x_dim_0 + 1; i++)
                            {
                                // Compute the linear index of the mid-point.
                                const int idx_midpoint_x = (i + num_ghosts_0) +
                                    (j + num_ghosts_1)*(ghostcell_dim_0 + 1) +
                                    (k + num_ghosts_2)*(ghostcell_dim_0 + 1)*
                                        ghostcell_dim_1;
                                
                                performLocalWENOInterpolationPlusLD(
                                    U_R,
                                    U_array.data(),
                                    idx_midpoint_x,
                                    constant_p,
                                    constant_q,
                                    constant_C,
                                    constant_alpha_tau);
                            }
                        }
                    }
                    break;
                }
            }
        }
        
        /*
         * Peform WENO interpolation in the y-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<Real*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(1));
            }
            
            Real* U_B = variables_minus[ei]->getPointer(1);
            
            switch(d_weno_interp)
            {
                case WENO_INTERP::TYPE::WENO5Z:
                {
                    for (int k = domain_y_lo_2; k < domain_y_lo_2 + domain_y_dim_2; k++)
                    {
                        for (int j = domain_y_lo_1; j < domain_y_lo_1 + domain_y_dim_1 + 1; j++)
                        {
                            HAMERS_PRAGMA_SIMD
                            for (int i = domain_y_lo_0; i < domain_y_lo_0 + domain_y_dim_0; i++)
                            {
                                // Compute the linear index of the mid-point.
                                const int idx_midpoint_y = (i + num_ghosts_0) +
                                    (j + num_ghosts_1)*ghostcell_dim_0 +
                                    (k + num_ghosts_2)*ghostcell_dim_0*
                                        (ghostcell_dim_1 + 1);
                                
                                performLocalWENOInterpolationMinusZ(
                                    U_B,
                                    U_array.data(),
                                    idx_midpoint_y,
                                    constant_p);
                            }
                        }
                    }
                    break;
                }
                case WENO_INTERP::TYPE::WENO6LD:
                {
                    for (int k = domain_y_lo_2; k < domain_y_lo_2 + domain_y_dim_2; k++)
                    {
                        for (int j = domain_y_lo_1; j < domain_y_lo_1 + domain_y_dim_1 + 1; j++)
                        {
                            HAMERS_PRAGMA_SIMD
                            for (int i = domain_y_lo_0; i < domain_y_lo_0 + domain_y_dim_0; i++)
                            {
                                // Compute the linear index of the mid-point.
                                const int idx_midpoint_y = (i + num_ghosts_0) +
                                    (j + num_ghosts_1)*ghostcell_dim_0 +
                                    (k + num_ghosts_2)*ghostcell_dim_0*
                                        (ghostcell_dim_1 + 1);
                                
                                performLocalWENOInterpolationMinusLD(
                                    U_B,
                                    U_array.data(),
                                    idx_midpoint_y,
                                    constant_p,
                                    constant_q,
                                    constant_C,
                                    constant_alpha_tau);
                            }
                        }
                    }
                    break;
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
            
            switch(d_weno_interp)
            {
                case WENO_INTERP::TYPE::WENO5Z:
                {
                    for (int k = domain_y_lo_2; k < domain_y_lo_2 + domain_y_dim_2; k++)
                    {
                        for (int j = domain_y_lo_1; j < domain_y_lo_1 + domain_y_dim_1 + 1; j++)
                        {
                            HAMERS_PRAGMA_SIMD
                            for (int i = domain_y_lo_0; i < domain_y_lo_0 + domain_y_dim_0; i++)
                            {
                                // Compute the linear index of the mid-point.
                                const int idx_midpoint_y = (i + num_ghosts_0) +
                                    (j + num_ghosts_1)*ghostcell_dim_0 +
                                    (k + num_ghosts_2)*ghostcell_dim_0*
                                        (ghostcell_dim_1 + 1);
                                
                                performLocalWENOInterpolationPlusZ(
                                    U_T,
                                    U_array.data(),
                                    idx_midpoint_y,
                                    constant_p);
                            }
                        }
                    }
                    break;
                }
                case WENO_INTERP::TYPE::WENO6LD:
                {
                    for (int k = domain_y_lo_2; k < domain_y_lo_2 + domain_y_dim_2; k++)
                    {
                        for (int j = domain_y_lo_1; j < domain_y_lo_1 + domain_y_dim_1 + 1; j++)
                        {
                            HAMERS_PRAGMA_SIMD
                            for (int i = domain_y_lo_0; i < domain_y_lo_0 + domain_y_dim_0; i++)
                            {
                                // Compute the linear index of the mid-point.
                                const int idx_midpoint_y = (i + num_ghosts_0) +
                                    (j + num_ghosts_1)*ghostcell_dim_0 +
                                    (k + num_ghosts_2)*ghostcell_dim_0*
                                        (ghostcell_dim_1 + 1);
                                
                                performLocalWENOInterpolationPlusLD(
                                    U_T,
                                    U_array.data(),
                                    idx_midpoint_y,
                                    constant_p,
                                    constant_q,
                                    constant_C,
                                    constant_alpha_tau);
                            }
                        }
                    }
                    break;
                }
            }
        }
        
        /*
         * Peform WENO interpolation in the z-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<Real*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(2));
            }
            
            Real* U_B = variables_minus[ei]->getPointer(2);
            
            switch(d_weno_interp)
            {
                case WENO_INTERP::TYPE::WENO5Z:
                {
                    for (int k = domain_z_lo_2; k < domain_z_lo_2 + domain_z_dim_2 + 1; k++)
                    {
                        for (int j = domain_z_lo_1; j < domain_z_lo_1 + domain_z_dim_1; j++)
                        {
                            HAMERS_PRAGMA_SIMD
                            for (int i = domain_z_lo_0; i < domain_z_lo_0 + domain_z_dim_0; i++)
                            {
                                // Compute the linear index of the mid-point.
                                const int idx_midpoint_z = (i + num_ghosts_0) +
                                    (j + num_ghosts_1)*ghostcell_dim_0 +
                                    (k + num_ghosts_2)*ghostcell_dim_0*
                                        ghostcell_dim_1;
                                
                                performLocalWENOInterpolationMinusZ(
                                    U_B,
                                    U_array.data(),
                                    idx_midpoint_z,
                                    constant_p);
                            }
                        }
                    }
                    break;
                }
                case WENO_INTERP::TYPE::WENO6LD:
                {
                    for (int k = domain_z_lo_2; k < domain_z_lo_2 + domain_z_dim_2 + 1; k++)
                    {
                        for (int j = domain_z_lo_1; j < domain_z_lo_1 + domain_z_dim_1; j++)
                        {
                            HAMERS_PRAGMA_SIMD
                            for (int i = domain_z_lo_0; i < domain_z_lo_0 + domain_z_dim_0; i++)
                            {
                                // Compute the linear index of the mid-point.
                                const int idx_midpoint_z = (i + num_ghosts_0) +
                                    (j + num_ghosts_1)*ghostcell_dim_0 +
                                    (k + num_ghosts_2)*ghostcell_dim_0*
                                        ghostcell_dim_1;
                                
                                performLocalWENOInterpolationMinusLD(
                                    U_B,
                                    U_array.data(),
                                    idx_midpoint_z,
                                    constant_p,
                                    constant_q,
                                    constant_C,
                                    constant_alpha_tau);
                            }
                        }
                    }
                    break;
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
            
            switch(d_weno_interp)
            {
                case WENO_INTERP::TYPE::WENO5Z:
                {
                    for (int k = domain_z_lo_2; k < domain_z_lo_2 + domain_z_dim_2 + 1; k++)
                    {
                        for (int j = domain_z_lo_1; j < domain_z_lo_1 + domain_z_dim_1; j++)
                        {
                            HAMERS_PRAGMA_SIMD
                            for (int i = domain_z_lo_0; i < domain_z_lo_0 + domain_z_dim_0; i++)
                            {
                                // Compute the linear index of the mid-point.
                                const int idx_midpoint_z = (i + num_ghosts_0) +
                                    (j + num_ghosts_1)*ghostcell_dim_0 +
                                    (k + num_ghosts_2)*ghostcell_dim_0*
                                        ghostcell_dim_1;
                                
                                performLocalWENOInterpolationPlusZ(
                                    U_F,
                                    U_array.data(),
                                    idx_midpoint_z,
                                    constant_p);
                            }
                        }
                    }
                    break;
                }
                case WENO_INTERP::TYPE::WENO6LD:
                {
                    for (int k = domain_z_lo_2; k < domain_z_lo_2 + domain_z_dim_2 + 1; k++)
                    {
                        for (int j = domain_z_lo_1; j < domain_z_lo_1 + domain_z_dim_1; j++)
                        {
                            HAMERS_PRAGMA_SIMD
                            for (int i = domain_z_lo_0; i < domain_z_lo_0 + domain_z_dim_0; i++)
                            {
                                // Compute the linear index of the mid-point.
                                const int idx_midpoint_z = (i + num_ghosts_0) +
                                    (j + num_ghosts_1)*ghostcell_dim_0 +
                                    (k + num_ghosts_2)*ghostcell_dim_0*
                                        ghostcell_dim_1;
                                
                                performLocalWENOInterpolationPlusLD(
                                    U_F,
                                    U_array.data(),
                                    idx_midpoint_z,
                                    constant_p,
                                    constant_q,
                                    constant_C,
                                    constant_alpha_tau);
                            }
                        }
                    }
                    break;
                }
            }
        }
    } // if (d_dim == tbox::Dimension(3))
}
