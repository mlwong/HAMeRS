#include "flow/flow_models/single-species/FlowModelRiemannSolverSingleSpecies.hpp"

#define EPSILON HAMERS_EPSILON


/*
 * Compute the local convective flux in the x-direction from conservative variables with
 * 1D HLLC-HLL Riemann solver.
 */
static inline __attribute__((always_inline)) void
computeLocalConvectiveFluxInXDirectionFromConservativeVariablesHLLC_HLL1D(
    Real** F_x,
    Real** Q_x_L,
    Real** Q_x_R,
    Real* p_x_L,
    Real* p_x_R,
    Real* c_x_L,
    Real* c_x_R,
    Real& u_x_L,
    Real& u_x_R,
    Real& s_x_minus,
    Real& s_x_plus,
    Real& s_x_star,
    Real& Chi_x_star_LR,
    const int& idx_flux,
    const int& idx)
{
    u_x_L = Q_x_L[1][idx]/Q_x_L[0][idx];
    u_x_R = Q_x_R[1][idx]/Q_x_R[0][idx];
    
    const Real u_x_average = Real(1)/Real(2)*(u_x_L + u_x_R);
    const Real c_x_average = Real(1)/Real(2)*(c_x_L[idx] + c_x_R[idx]);
    
    const Real s_x_L = std::min(u_x_average - c_x_average, u_x_L - c_x_L[idx]);
    const Real s_x_R = std::max(u_x_average + c_x_average, u_x_R + c_x_R[idx]);
    
    s_x_minus = std::min(Real(0), s_x_L);
    s_x_plus  = std::max(Real(0), s_x_R);
    
    s_x_star = (p_x_R[idx] - p_x_L[idx] +
        Q_x_L[1][idx]*(s_x_L - u_x_L) - Q_x_R[1][idx]*(s_x_R - u_x_R))/
        (Q_x_L[0][idx]*(s_x_L - u_x_L) - Q_x_R[0][idx]*(s_x_R - u_x_R));
    
    Real Q_x_star_LR[3];
    Real F_x_LR[3];
    
    if (s_x_star > Real(0))
    {
        Chi_x_star_LR = (s_x_L - u_x_L)/(s_x_L - s_x_star);
        
        Q_x_star_LR[0] = Chi_x_star_LR*Q_x_L[0][idx];
        Q_x_star_LR[1] = Chi_x_star_LR*Q_x_L[0][idx]*s_x_star;
        Q_x_star_LR[2] = Chi_x_star_LR*(Q_x_L[2][idx] + (s_x_star - u_x_L)*(Q_x_L[0][idx]*s_x_star +
            p_x_L[idx]/(s_x_L - u_x_L)));
        
        F_x_LR[0] = Q_x_L[1][idx];
        F_x_LR[1] = u_x_L*Q_x_L[1][idx] + p_x_L[idx];
        F_x_LR[2] = u_x_L*(Q_x_L[2][idx] + p_x_L[idx]);
        
        for (int ei = 0; ei < 3; ei++)
        {
            F_x[ei][idx_flux] = F_x_LR[ei] + s_x_minus*(Q_x_star_LR[ei] - Q_x_L[ei][idx]);
        }
    }
    else
    {
        Chi_x_star_LR = (s_x_R - u_x_R)/(s_x_R - s_x_star);
        
        Q_x_star_LR[0] = Chi_x_star_LR*Q_x_R[0][idx];
        Q_x_star_LR[1] = Chi_x_star_LR*Q_x_R[0][idx]*s_x_star;
        Q_x_star_LR[2] = Chi_x_star_LR*(Q_x_R[2][idx] + (s_x_star - u_x_R)*(Q_x_R[0][idx]*s_x_star +
            p_x_R[idx]/(s_x_R - u_x_R)));
        
        F_x_LR[0] = Q_x_R[1][idx];
        F_x_LR[1] = u_x_R*Q_x_R[1][idx] + p_x_R[idx];
        F_x_LR[2] = u_x_R*(Q_x_R[2][idx] + p_x_R[idx]);
        
        for (int ei = 0; ei < 3; ei++)
        {
            F_x[ei][idx_flux] = F_x_LR[ei] + s_x_plus*(Q_x_star_LR[ei] - Q_x_R[ei][idx]);
        }
    }
}


/*
 * Compute the local convective flux in the x-direction from conservative variables with
 * 2D HLLC-HLL Riemann solver.
 */
static inline __attribute__((always_inline)) void
computeLocalConvectiveFluxInXDirectionFromConservativeVariablesHLLC_HLL2D(
    Real** F_x,
    Real** Q_x_L,
    Real** Q_x_R,
    Real* p_x_L,
    Real* p_x_R,
    Real* c_x_L,
    Real* c_x_R,
    Real& u_x_L,
    Real& u_x_R,
    Real& s_x_minus,
    Real& s_x_plus,
    Real& s_x_star,
    Real& Chi_x_star_LR,
    const int& idx_flux,
    const int& idx)
{
    u_x_L = Q_x_L[1][idx]/Q_x_L[0][idx];
    u_x_R = Q_x_R[1][idx]/Q_x_R[0][idx];
    
    const Real v_x_L = Q_x_L[2][idx]/Q_x_L[0][idx];
    const Real v_x_R = Q_x_R[2][idx]/Q_x_R[0][idx];
    
    const Real u_x_average = Real(1)/Real(2)*(u_x_L + u_x_R);
    const Real c_x_average = Real(1)/Real(2)*(c_x_L[idx] + c_x_R[idx]);
    
    const Real s_x_L = std::min(u_x_average - c_x_average, u_x_L - c_x_L[idx]);
    const Real s_x_R = std::max(u_x_average + c_x_average, u_x_R + c_x_R[idx]);
    
    s_x_minus = std::min(Real(0), s_x_L);
    s_x_plus  = std::max(Real(0), s_x_R);
    
    s_x_star = (p_x_R[idx] - p_x_L[idx] +
        Q_x_L[1][idx]*(s_x_L - u_x_L) - Q_x_R[1][idx]*(s_x_R - u_x_R))/
        (Q_x_L[0][idx]*(s_x_L - u_x_L) - Q_x_R[0][idx]*(s_x_R - u_x_R));
    
    Real F_x_L[4];
    Real F_x_R[4];
    Real F_x_HLL[2];
    Real F_x_HLLC[4];
    Real Q_x_star_LR[4];
    
    F_x_L[0] = Q_x_L[1][idx];
    F_x_L[1] = u_x_L*Q_x_L[1][idx] + p_x_L[idx];
    F_x_L[2] = u_x_L*Q_x_L[2][idx];
    F_x_L[3] = u_x_L*(Q_x_L[3][idx] + p_x_L[idx]);
    
    F_x_R[0] = Q_x_R[1][idx];
    F_x_R[1] = u_x_R*Q_x_R[1][idx] + p_x_R[idx];
    F_x_R[2] = u_x_R*Q_x_R[2][idx];
    F_x_R[3] = u_x_R*(Q_x_R[3][idx] + p_x_R[idx]);
    
    F_x_HLL[0] = (s_x_R*F_x_L[0] - s_x_L*F_x_R[0] + s_x_R*s_x_L*(Q_x_R[0][idx] - Q_x_L[0][idx]))/
        (s_x_R - s_x_L);
    F_x_HLL[1] = (s_x_R*F_x_L[2] - s_x_L*F_x_R[2] + s_x_R*s_x_L*(Q_x_R[2][idx] - Q_x_L[2][idx]))/
        (s_x_R - s_x_L);
    
    if (s_x_L > Real(0))
    {
        F_x_HLL[0] = F_x_L[0];
        F_x_HLL[1] = F_x_L[2];
    }
    
    if (s_x_R < Real(0))
    {
        F_x_HLL[0] = F_x_R[0];
        F_x_HLL[1] = F_x_R[2];
    }
    
    if (s_x_star > Real(0))
    {
        Chi_x_star_LR = (s_x_L - u_x_L)/(s_x_L - s_x_star);
        
        Q_x_star_LR[0] = Chi_x_star_LR*Q_x_L[0][idx];
        Q_x_star_LR[1] = Chi_x_star_LR*Q_x_L[0][idx]*s_x_star;
        Q_x_star_LR[2] = Chi_x_star_LR*Q_x_L[2][idx];
        Q_x_star_LR[3] = Chi_x_star_LR*(Q_x_L[3][idx] + (s_x_star - u_x_L)*(Q_x_L[0][idx]*s_x_star +
            p_x_L[idx]/(s_x_L - u_x_L)));
        
        for (int ei = 0; ei < 4; ei++)
        {
            F_x_HLLC[ei] = F_x_L[ei] + s_x_minus*(Q_x_star_LR[ei] - Q_x_L[ei][idx]);
        }
    }
    else
    {
        Chi_x_star_LR = (s_x_R - u_x_R)/(s_x_R - s_x_star);
        
        Q_x_star_LR[0] = Chi_x_star_LR*Q_x_R[0][idx];
        Q_x_star_LR[1] = Chi_x_star_LR*Q_x_R[0][idx]*s_x_star;
        Q_x_star_LR[2] = Chi_x_star_LR*Q_x_R[2][idx];
        Q_x_star_LR[3] = Chi_x_star_LR*(Q_x_R[3][idx] + (s_x_star - u_x_R)*(Q_x_R[0][idx]*s_x_star +
            p_x_R[idx]/(s_x_R - u_x_R)));
        
        for (int ei = 0; ei < 4; ei++)
        {
            F_x_HLLC[ei] = F_x_R[ei] + s_x_plus*(Q_x_star_LR[ei] - Q_x_R[ei][idx]);
        }
    }
    
    /*
     * Calulate the weights beta for hybridization.
     */
    
    const Real u_x_diff = u_x_R - u_x_L;
    const Real v_x_diff = v_x_R - v_x_L;
    const Real vel_mag = std::sqrt(u_x_diff*u_x_diff + v_x_diff*v_x_diff);
    
    Real alpha_1, alpha_2;
    if (vel_mag < EPSILON)
    {
        alpha_1 = Real(1);
        alpha_2 = Real(0);
    }
    else
    {
        alpha_1 = std::abs(u_x_diff)/vel_mag;
        alpha_2 = std::sqrt(Real(1) - alpha_1*alpha_1);
    }
    
    const Real beta_1 = Real(1)/Real(2)*(Real(1) + alpha_1/(alpha_1 + alpha_2));
    const Real beta_2 = Real(1) - beta_1;
    
    F_x[0][idx_flux] = beta_1*F_x_HLLC[0] + beta_2*F_x_HLL[0];
    F_x[1][idx_flux] = F_x_HLLC[1];
    F_x[2][idx_flux] = beta_1*F_x_HLLC[2] + beta_2*F_x_HLL[1];
    F_x[3][idx_flux] = F_x_HLLC[3];
}


/*
 * Compute the local convective flux in the x-direction from conservative variables with
 * 3D HLLC-HLL Riemann solver.
 */
static inline __attribute__((always_inline)) void
computeLocalConvectiveFluxInXDirectionFromConservativeVariablesHLLC_HLL3D(
    Real** F_x,
    Real** Q_x_L,
    Real** Q_x_R,
    Real* p_x_L,
    Real* p_x_R,
    Real* c_x_L,
    Real* c_x_R,
    Real& u_x_L,
    Real& u_x_R,
    Real& s_x_minus,
    Real& s_x_plus,
    Real& s_x_star,
    Real& Chi_x_star_LR,
    const int& idx_flux,
    const int& idx)
{
    u_x_L = Q_x_L[1][idx]/Q_x_L[0][idx];
    u_x_R = Q_x_R[1][idx]/Q_x_R[0][idx];
    
    const Real v_x_L = Q_x_L[2][idx]/Q_x_L[0][idx];
    const Real v_x_R = Q_x_R[2][idx]/Q_x_R[0][idx];
    
    const Real w_x_L = Q_x_L[3][idx]/Q_x_L[0][idx];
    const Real w_x_R = Q_x_R[3][idx]/Q_x_R[0][idx];
    
    const Real u_x_average = Real(1)/Real(2)*(u_x_L + u_x_R);
    const Real c_x_average = Real(1)/Real(2)*(c_x_L[idx] + c_x_R[idx]);
    
    const Real s_x_L = std::min(u_x_average - c_x_average, u_x_L - c_x_L[idx]);
    const Real s_x_R = std::max(u_x_average + c_x_average, u_x_R + c_x_R[idx]);
    
    s_x_minus = std::min(Real(0), s_x_L);
    s_x_plus  = std::max(Real(0), s_x_R);
    
    s_x_star = (p_x_R[idx] - p_x_L[idx] +
        Q_x_L[1][idx]*(s_x_L - u_x_L) - Q_x_R[1][idx]*(s_x_R - u_x_R))/
        (Q_x_L[0][idx]*(s_x_L - u_x_L) - Q_x_R[0][idx]*(s_x_R - u_x_R));
    
    Real F_x_L[5];
    Real F_x_R[5];
    Real F_x_HLL[3];
    Real F_x_HLLC[5];
    Real Q_x_star_LR[5];
    
    F_x_L[0] = Q_x_L[1][idx];
    F_x_L[1] = u_x_L*Q_x_L[1][idx] + p_x_L[idx];
    F_x_L[2] = u_x_L*Q_x_L[2][idx];
    F_x_L[3] = u_x_L*Q_x_L[3][idx];
    F_x_L[4] = u_x_L*(Q_x_L[4][idx] + p_x_L[idx]);
    
    F_x_R[0] = Q_x_R[1][idx];
    F_x_R[1] = u_x_R*Q_x_R[1][idx] + p_x_R[idx];
    F_x_R[2] = u_x_R*Q_x_R[2][idx];
    F_x_R[3] = u_x_R*Q_x_R[3][idx];
    F_x_R[4] = u_x_R*(Q_x_R[4][idx] + p_x_R[idx]);
    
    F_x_HLL[0] = (s_x_R*F_x_L[0] - s_x_L*F_x_R[0] + s_x_R*s_x_L*(Q_x_R[0][idx] - Q_x_L[0][idx]))/
        (s_x_R - s_x_L);
    F_x_HLL[1] = (s_x_R*F_x_L[2] - s_x_L*F_x_R[2] + s_x_R*s_x_L*(Q_x_R[2][idx] - Q_x_L[2][idx]))/
        (s_x_R - s_x_L);
    F_x_HLL[2] = (s_x_R*F_x_L[3] - s_x_L*F_x_R[3] + s_x_R*s_x_L*(Q_x_R[3][idx] - Q_x_L[3][idx]))/
        (s_x_R - s_x_L);
    
    if (s_x_L > Real(0))
    {
        F_x_HLL[0] = F_x_L[0];
        F_x_HLL[1] = F_x_L[2];
        F_x_HLL[2] = F_x_L[3];
    }
    
    if (s_x_R < Real(0))
    {
        F_x_HLL[0] = F_x_R[0];
        F_x_HLL[1] = F_x_R[2];
        F_x_HLL[2] = F_x_R[3];
    }
    
    if (s_x_star > Real(0))
    {
        Chi_x_star_LR = (s_x_L - u_x_L)/(s_x_L - s_x_star);
        
        Q_x_star_LR[0] = Chi_x_star_LR*Q_x_L[0][idx];
        Q_x_star_LR[1] = Chi_x_star_LR*Q_x_L[0][idx]*s_x_star;
        Q_x_star_LR[2] = Chi_x_star_LR*Q_x_L[2][idx];
        Q_x_star_LR[3] = Chi_x_star_LR*Q_x_L[3][idx];
        Q_x_star_LR[4] = Chi_x_star_LR*(Q_x_L[4][idx] + (s_x_star - u_x_L)*(Q_x_L[0][idx]*s_x_star +
            p_x_L[idx]/(s_x_L - u_x_L)));
        
        for (int ei = 0; ei < 5; ei++)
        {
            F_x_HLLC[ei] = F_x_L[ei] + s_x_minus*(Q_x_star_LR[ei] - Q_x_L[ei][idx]);
        }
    }
    else
    {
        Chi_x_star_LR = (s_x_R - u_x_R)/(s_x_R - s_x_star);
        
        Q_x_star_LR[0] = Chi_x_star_LR*Q_x_R[0][idx];
        Q_x_star_LR[1] = Chi_x_star_LR*Q_x_R[0][idx]*s_x_star;
        Q_x_star_LR[2] = Chi_x_star_LR*Q_x_R[2][idx];
        Q_x_star_LR[3] = Chi_x_star_LR*Q_x_R[3][idx];
        Q_x_star_LR[4] = Chi_x_star_LR*(Q_x_R[4][idx] + (s_x_star - u_x_R)*(Q_x_R[0][idx]*s_x_star +
            p_x_R[idx]/(s_x_R - u_x_R)));
        
        for (int ei = 0; ei < 5; ei++)
        {
            F_x_HLLC[ei] = F_x_R[ei] + s_x_plus*(Q_x_star_LR[ei] - Q_x_R[ei][idx]);
        }
    }

    /*
     * Calulate the weights beta for hybridization.
     */
    
    const Real u_x_diff = u_x_R - u_x_L;
    const Real v_x_diff = v_x_R - v_x_L;
    const Real w_x_diff = w_x_R - w_x_L;
    const Real vel_mag = std::sqrt(u_x_diff*u_x_diff + v_x_diff*v_x_diff + w_x_diff*w_x_diff);
    
    Real alpha_1, alpha_2;
    if (vel_mag < EPSILON)
    {
        alpha_1 = Real(1);
        alpha_2 = Real(0);
    }
    else
    {
        alpha_1 = std::abs(u_x_diff)/vel_mag;
        alpha_2 = std::sqrt(Real(1) - alpha_1*alpha_1);
    }
    
    const Real beta_1 = Real(1)/Real(2)*(Real(1) + alpha_1/(alpha_1 + alpha_2));
    const Real beta_2 = Real(1) - beta_1;
    
    F_x[0][idx_flux] = beta_1*F_x_HLLC[0] + beta_2*F_x_HLL[0];
    F_x[1][idx_flux] = F_x_HLLC[1];
    F_x[2][idx_flux] = beta_1*F_x_HLLC[2] + beta_2*F_x_HLL[1];
    F_x[3][idx_flux] = beta_1*F_x_HLLC[3] + beta_2*F_x_HLL[2];
    F_x[4][idx_flux] = F_x_HLLC[4];
}


/*
 * Compute the local convective flux in the y-direction from conservative variables with
 * 2D HLLC-HLL Riemann solver.
 */
static inline __attribute__((always_inline)) void
computeLocalConvectiveFluxInYDirectionFromConservativeVariablesHLLC_HLL2D(
    Real** F_y,
    Real** Q_y_B,
    Real** Q_y_T,
    Real* p_y_B,
    Real* p_y_T,
    Real* c_y_B,
    Real* c_y_T,
    Real& v_y_B,
    Real& v_y_T,
    Real& s_y_minus,
    Real& s_y_plus,
    Real& s_y_star,
    Real& Chi_y_star_BT,
    const int& idx_flux,
    const int& idx)
{
    v_y_B = Q_y_B[2][idx]/Q_y_B[0][idx];
    v_y_T = Q_y_T[2][idx]/Q_y_T[0][idx];
    
    const Real u_y_B = Q_y_B[1][idx]/Q_y_B[0][idx];
    const Real u_y_T = Q_y_T[1][idx]/Q_y_T[0][idx];
    
    const Real v_y_average = Real(1)/Real(2)*(v_y_B + v_y_T);
    const Real c_y_average = Real(1)/Real(2)*(c_y_B[idx] + c_y_T[idx]);
    
    const Real s_y_B = std::min(v_y_average - c_y_average, v_y_B - c_y_B[idx]);
    const Real s_y_T = std::max(v_y_average + c_y_average, v_y_T + c_y_T[idx]);
    
    s_y_minus = std::min(Real(0), s_y_B);
    s_y_plus  = std::max(Real(0), s_y_T);
    
    s_y_star = (p_y_T[idx] - p_y_B[idx] +
        Q_y_B[2][idx]*(s_y_B - v_y_B) - Q_y_T[2][idx]*(s_y_T - v_y_T))/
        (Q_y_B[0][idx]*(s_y_B - v_y_B) - Q_y_T[0][idx]*(s_y_T - v_y_T));
    
    Real F_y_B[4];
    Real F_y_T[4];
    Real F_y_HLL[2];
    Real F_y_HLLC[4];
    Real Q_y_star_BT[4];
    
    F_y_B[0] = Q_y_B[2][idx];
    F_y_B[1] = v_y_B*Q_y_B[1][idx];
    F_y_B[2] = v_y_B*Q_y_B[2][idx] + p_y_B[idx];
    F_y_B[3] = v_y_B*(Q_y_B[3][idx] + p_y_B[idx]);
    
    F_y_T[0] = Q_y_T[2][idx];
    F_y_T[1] = v_y_T*Q_y_T[1][idx];
    F_y_T[2] = v_y_T*Q_y_T[2][idx] + p_y_T[idx];
    F_y_T[3] = v_y_T*(Q_y_T[3][idx] + p_y_T[idx]);
    
    F_y_HLL[0] = (s_y_T*F_y_B[0] - s_y_B*F_y_T[0] + s_y_T*s_y_B*(Q_y_T[0][idx] - Q_y_B[0][idx]))/
        (s_y_T - s_y_B);
    F_y_HLL[1] = (s_y_T*F_y_B[1] - s_y_B*F_y_T[1] + s_y_T*s_y_B*(Q_y_T[1][idx] - Q_y_B[1][idx]))/
        (s_y_T - s_y_B);
    
    if (s_y_B > Real(0))
    {
        F_y_HLL[0] = F_y_B[0];
        F_y_HLL[1] = F_y_B[1];
    }
    
    if (s_y_T < Real(0))
    {
        F_y_HLL[0] = F_y_T[0];
        F_y_HLL[1] = F_y_T[1];
    }
    
    if (s_y_star > Real(0))
    {
        Chi_y_star_BT = (s_y_B - v_y_B)/(s_y_B - s_y_star);
        
        Q_y_star_BT[0] = Chi_y_star_BT*Q_y_B[0][idx];
        Q_y_star_BT[1] = Chi_y_star_BT*Q_y_B[1][idx];
        Q_y_star_BT[2] = Chi_y_star_BT*Q_y_B[0][idx]*s_y_star;
        Q_y_star_BT[3] = Chi_y_star_BT*(Q_y_B[3][idx] + (s_y_star - v_y_B)*(Q_y_B[0][idx]*s_y_star +
            p_y_B[idx]/(s_y_B - v_y_B)));
        
        for (int ei = 0; ei < 4; ei++)
        {
            F_y_HLLC[ei] = F_y_B[ei] + s_y_minus*(Q_y_star_BT[ei] - Q_y_B[ei][idx]);
        }
    }
    else
    {
        Chi_y_star_BT = (s_y_T - v_y_T)/(s_y_T - s_y_star);
        
        Q_y_star_BT[0] = Chi_y_star_BT*Q_y_T[0][idx];
        Q_y_star_BT[1] = Chi_y_star_BT*Q_y_T[1][idx];
        Q_y_star_BT[2] = Chi_y_star_BT*Q_y_T[0][idx]*s_y_star;
        Q_y_star_BT[3] = Chi_y_star_BT*(Q_y_T[3][idx] + (s_y_star - v_y_T)*(Q_y_T[0][idx]*s_y_star +
            p_y_T[idx]/(s_y_T - v_y_T)));
        
        for (int ei = 0; ei < 4; ei++)
        {
            F_y_HLLC[ei] = F_y_T[ei] + s_y_plus*(Q_y_star_BT[ei] - Q_y_T[ei][idx]);
        }
    }
    
    /*
     * Calulate the weights beta for hybridization.
     */
    
    const Real u_y_diff = u_y_T - u_y_B;
    const Real v_y_diff = v_y_T - v_y_B;
    const Real vel_mag = std::sqrt(u_y_diff*u_y_diff + v_y_diff*v_y_diff);
    
    Real alpha_1, alpha_2;
    if (vel_mag < EPSILON)
    {
        alpha_1 = Real(1);
        alpha_2 = Real(0);
    }
    else
    {
        alpha_1 = std::abs(v_y_diff)/vel_mag;
        alpha_2 = std::sqrt(Real(1) - alpha_1*alpha_1);
    }
    
    const Real beta_1 = Real(1)/Real(2)*(Real(1) + alpha_1/(alpha_1 + alpha_2));
    const Real beta_2 = Real(1) - beta_1;
    
    F_y[0][idx_flux] = beta_1*F_y_HLLC[0] + beta_2*F_y_HLL[0];
    F_y[1][idx_flux] = beta_1*F_y_HLLC[1] + beta_2*F_y_HLL[1];
    F_y[2][idx_flux] = F_y_HLLC[2];
    F_y[3][idx_flux] = F_y_HLLC[3];
}


/*
 * Compute the local convective flux in the y-direction from conservative variables with
 * 3D HLLC-HLL Riemann solver.
 */
static inline __attribute__((always_inline)) void
computeLocalConvectiveFluxInYDirectionFromConservativeVariablesHLLC_HLL3D(
    Real** F_y,
    Real** Q_y_B,
    Real** Q_y_T,
    Real* p_y_B,
    Real* p_y_T,
    Real* c_y_B,
    Real* c_y_T,
    Real& v_y_B,
    Real& v_y_T,
    Real& s_y_minus,
    Real& s_y_plus,
    Real& s_y_star,
    Real& Chi_y_star_BT,
    const int& idx_flux,
    const int& idx)
{
    v_y_B = Q_y_B[2][idx]/Q_y_B[0][idx];
    v_y_T = Q_y_T[2][idx]/Q_y_T[0][idx];
    
    const Real u_y_B = Q_y_B[1][idx]/Q_y_B[0][idx];
    const Real u_y_T = Q_y_T[1][idx]/Q_y_T[0][idx];
    
    const Real w_y_B = Q_y_B[3][idx]/Q_y_B[0][idx];
    const Real w_y_T = Q_y_T[3][idx]/Q_y_T[0][idx];
    
    const Real v_y_average = Real(1)/Real(2)*(v_y_B + v_y_T);
    const Real c_y_average = Real(1)/Real(2)*(c_y_B[idx] + c_y_T[idx]);
    
    const Real s_y_B = std::min(v_y_average - c_y_average, v_y_B - c_y_B[idx]);
    const Real s_y_T = std::max(v_y_average + c_y_average, v_y_T + c_y_T[idx]);
    
    s_y_minus = std::min(Real(0), s_y_B);
    s_y_plus  = std::max(Real(0), s_y_T);
    
    s_y_star = (p_y_T[idx] - p_y_B[idx] +
        Q_y_B[2][idx]*(s_y_B - v_y_B) - Q_y_T[2][idx]*(s_y_T - v_y_T))/
        (Q_y_B[0][idx]*(s_y_B - v_y_B) - Q_y_T[0][idx]*(s_y_T - v_y_T));
    
    Real F_y_B[5];
    Real F_y_T[5];
    Real F_y_HLL[3];
    Real F_y_HLLC[5];
    Real Q_y_star_BT[5];
    
    F_y_B[0] = Q_y_B[2][idx];
    F_y_B[1] = v_y_B*Q_y_B[1][idx];
    F_y_B[2] = v_y_B*Q_y_B[2][idx] + p_y_B[idx];
    F_y_B[3] = v_y_B*Q_y_B[3][idx];
    F_y_B[4] = v_y_B*(Q_y_B[4][idx] + p_y_B[idx]);
    
    F_y_T[0] = Q_y_T[2][idx];
    F_y_T[1] = v_y_T*Q_y_T[1][idx];
    F_y_T[2] = v_y_T*Q_y_T[2][idx] + p_y_T[idx];
    F_y_T[3] = v_y_T*Q_y_T[3][idx];
    F_y_T[4] = v_y_T*(Q_y_T[4][idx] + p_y_T[idx]);
    
    F_y_HLL[0] = (s_y_T*F_y_B[0] - s_y_B*F_y_T[0] + s_y_T*s_y_B*(Q_y_T[0][idx] - Q_y_B[0][idx]))/
        (s_y_T - s_y_B);
    F_y_HLL[1] = (s_y_T*F_y_B[1] - s_y_B*F_y_T[1] + s_y_T*s_y_B*(Q_y_T[1][idx] - Q_y_B[1][idx]))/
        (s_y_T - s_y_B);
    F_y_HLL[2] = (s_y_T*F_y_B[3] - s_y_B*F_y_T[3] + s_y_T*s_y_B*(Q_y_T[3][idx] - Q_y_B[3][idx]))/
        (s_y_T - s_y_B);
    
    if (s_y_B > Real(0))
    {
        F_y_HLL[0] = F_y_B[0];
        F_y_HLL[1] = F_y_B[1];
        F_y_HLL[2] = F_y_B[3];
    }
    
    if (s_y_T < Real(0))
    {
        F_y_HLL[0] = F_y_T[0];
        F_y_HLL[1] = F_y_T[1];
        F_y_HLL[2] = F_y_T[3];
    }
    
    if (s_y_star > Real(0))
    {
        Chi_y_star_BT = (s_y_B - v_y_B)/(s_y_B - s_y_star);
        
        Q_y_star_BT[0] = Chi_y_star_BT*Q_y_B[0][idx];
        Q_y_star_BT[1] = Chi_y_star_BT*Q_y_B[1][idx];
        Q_y_star_BT[2] = Chi_y_star_BT*Q_y_B[0][idx]*s_y_star;
        Q_y_star_BT[3] = Chi_y_star_BT*Q_y_B[3][idx];
        Q_y_star_BT[4] = Chi_y_star_BT*(Q_y_B[4][idx] + (s_y_star - v_y_B)*(Q_y_B[0][idx]*s_y_star +
            p_y_B[idx]/(s_y_B - v_y_B)));
        
        for (int ei = 0; ei < 5; ei++)
        {
            F_y_HLLC[ei] = F_y_B[ei] + s_y_minus*(Q_y_star_BT[ei] - Q_y_B[ei][idx]);
        }
    }
    else
    {
        Chi_y_star_BT = (s_y_T - v_y_T)/(s_y_T - s_y_star);
        
        Q_y_star_BT[0] = Chi_y_star_BT*Q_y_T[0][idx];
        Q_y_star_BT[1] = Chi_y_star_BT*Q_y_T[1][idx];
        Q_y_star_BT[2] = Chi_y_star_BT*Q_y_T[0][idx]*s_y_star;
        Q_y_star_BT[3] = Chi_y_star_BT*Q_y_T[3][idx];
        Q_y_star_BT[4] = Chi_y_star_BT*(Q_y_T[4][idx] + (s_y_star - v_y_T)*(Q_y_T[0][idx]*s_y_star +
            p_y_T[idx]/(s_y_T - v_y_T)));
        
        for (int ei = 0; ei < 5; ei++)
        {
            F_y_HLLC[ei] = F_y_T[ei] + s_y_plus*(Q_y_star_BT[ei] - Q_y_T[ei][idx]);
        }
    }
    
    /*
     * Calulate the weights beta for hybridization.
     */
    
    const Real u_y_diff = u_y_T - u_y_B;
    const Real v_y_diff = v_y_T - v_y_B;
    const Real w_y_diff = w_y_T - w_y_B;
    const Real vel_mag = std::sqrt(u_y_diff*u_y_diff + v_y_diff*v_y_diff + w_y_diff*w_y_diff);
    
    Real alpha_1, alpha_2;
    if (vel_mag < EPSILON)
    {
        alpha_1 = Real(1);
        alpha_2 = Real(0);
    }
    else
    {
        alpha_1 = std::abs(v_y_diff)/vel_mag;
        alpha_2 = std::sqrt(Real(1) - alpha_1*alpha_1);
    }
    
    const Real beta_1 = Real(1)/Real(2)*(Real(1) + alpha_1/(alpha_1 + alpha_2));
    const Real beta_2 = Real(1) - beta_1;
    
    F_y[0][idx_flux] = beta_1*F_y_HLLC[0] + beta_2*F_y_HLL[0];
    F_y[1][idx_flux] = beta_1*F_y_HLLC[1] + beta_2*F_y_HLL[1];
    F_y[2][idx_flux] = F_y_HLLC[2];
    F_y[3][idx_flux] = beta_1*F_y_HLLC[3] + beta_2*F_y_HLL[2];
    F_y[4][idx_flux] = F_y_HLLC[4];
}


/*
 * Compute the local convective flux in the z-direction from conservative variables with
 * 3D HLLC-HLL Riemann solver.
 */
static inline __attribute__((always_inline)) void
computeLocalConvectiveFluxInZDirectionFromConservativeVariablesHLLC_HLL3D(
    Real** F_z,
    Real** Q_z_B,
    Real** Q_z_F,
    Real* p_z_B,
    Real* p_z_F,
    Real* c_z_B,
    Real* c_z_F,
    Real& w_z_B,
    Real& w_z_F,
    Real& s_z_minus,
    Real& s_z_plus,
    Real& s_z_star,
    Real& Chi_z_star_BF,
    const int& idx_flux,
    const int& idx)
{
    w_z_B = Q_z_B[3][idx]/Q_z_B[0][idx];
    w_z_F = Q_z_F[3][idx]/Q_z_F[0][idx];
    
    const Real u_z_B = Q_z_B[1][idx]/Q_z_B[0][idx];
    const Real u_z_F = Q_z_F[1][idx]/Q_z_F[0][idx];
    
    const Real v_z_B = Q_z_B[2][idx]/Q_z_B[0][idx];
    const Real v_z_F = Q_z_F[2][idx]/Q_z_F[0][idx];
    
    const Real w_z_average = Real(1)/Real(2)*(w_z_B + w_z_F);
    const Real c_z_average = Real(1)/Real(2)*(c_z_B[idx] + c_z_F[idx]);
    
    const Real s_z_B = std::min(w_z_average - c_z_average, w_z_B - c_z_B[idx]);
    const Real s_z_F = std::max(w_z_average + c_z_average, w_z_F + c_z_F[idx]);
    
    s_z_minus = std::min(Real(0), s_z_B);
    s_z_plus  = std::max(Real(0), s_z_F);
    
    s_z_star = (p_z_F[idx] - p_z_B[idx] +
        Q_z_B[3][idx]*(s_z_B - w_z_B) - Q_z_F[3][idx]*(s_z_F - w_z_F))/
        (Q_z_B[0][idx]*(s_z_B - w_z_B) - Q_z_F[0][idx]*(s_z_F - w_z_F));
    
    Real F_z_B[5];
    Real F_z_F[5];
    Real F_z_HLL[3];
    Real F_z_HLLC[5];
    Real Q_z_star_BF[5];
    
    F_z_B[0] = Q_z_B[3][idx];
    F_z_B[1] = w_z_B*Q_z_B[1][idx];
    F_z_B[2] = w_z_B*Q_z_B[2][idx];
    F_z_B[3] = w_z_B*Q_z_B[3][idx] + p_z_B[idx];
    F_z_B[4] = w_z_B*(Q_z_B[4][idx] + p_z_B[idx]);
    
    F_z_F[0] = Q_z_F[3][idx];
    F_z_F[1] = w_z_F*Q_z_F[1][idx];
    F_z_F[2] = w_z_F*Q_z_F[2][idx];
    F_z_F[3] = w_z_F*Q_z_F[3][idx] + p_z_F[idx];
    F_z_F[4] = w_z_F*(Q_z_F[4][idx] + p_z_F[idx]);
    
    F_z_HLL[0] = (s_z_F*F_z_B[0] - s_z_B*F_z_F[0] + s_z_F*s_z_B*(Q_z_F[0][idx] - Q_z_B[0][idx]))/
        (s_z_F - s_z_B);
    F_z_HLL[1] = (s_z_F*F_z_B[1] - s_z_B*F_z_F[1] + s_z_F*s_z_B*(Q_z_F[1][idx] - Q_z_B[1][idx]))/
        (s_z_F - s_z_B);
    F_z_HLL[2] = (s_z_F*F_z_B[2] - s_z_B*F_z_F[2] + s_z_F*s_z_B*(Q_z_F[2][idx] - Q_z_B[2][idx]))/
        (s_z_F - s_z_B);
    
    if (s_z_B > Real(0))
    {
        F_z_HLL[0] = F_z_B[0];
        F_z_HLL[1] = F_z_B[1];
        F_z_HLL[2] = F_z_B[2];
    }
    
    if (s_z_F < Real(0))
    {
        F_z_HLL[0] = F_z_F[0];
        F_z_HLL[1] = F_z_F[1];
        F_z_HLL[2] = F_z_F[2];
    }
    
    if (s_z_star > Real(0))
    {
        Chi_z_star_BF = (s_z_B - w_z_B)/(s_z_B - s_z_star);
        
        Q_z_star_BF[0] = Chi_z_star_BF*Q_z_B[0][idx];
        Q_z_star_BF[1] = Chi_z_star_BF*Q_z_B[1][idx];
        Q_z_star_BF[2] = Chi_z_star_BF*Q_z_B[2][idx];
        Q_z_star_BF[3] = Chi_z_star_BF*Q_z_B[0][idx]*s_z_star;
        Q_z_star_BF[4] = Chi_z_star_BF*(Q_z_B[4][idx] + (s_z_star - w_z_B)*(Q_z_B[0][idx]*s_z_star +
            p_z_B[idx]/(s_z_B - w_z_B)));
        
        for (int ei = 0; ei < 5; ei++)
        {
            F_z_HLLC[ei] = F_z_B[ei] + s_z_minus*(Q_z_star_BF[ei] - Q_z_B[ei][idx]);
        }
    }
    else
    {
        Chi_z_star_BF = (s_z_F - w_z_F)/(s_z_F - s_z_star);
        
        Q_z_star_BF[0] = Chi_z_star_BF*Q_z_F[0][idx];
        Q_z_star_BF[1] = Chi_z_star_BF*Q_z_F[1][idx];
        Q_z_star_BF[2] = Chi_z_star_BF*Q_z_F[2][idx];
        Q_z_star_BF[3] = Chi_z_star_BF*Q_z_F[0][idx]*s_z_star;
        Q_z_star_BF[4] = Chi_z_star_BF*(Q_z_F[4][idx] + (s_z_star - w_z_F)*(Q_z_F[0][idx]*s_z_star +
            p_z_F[idx]/(s_z_F - w_z_F)));
        
        for (int ei = 0; ei < 5; ei++)
        {
            F_z_HLLC[ei] = F_z_F[ei] + s_z_plus*(Q_z_star_BF[ei] - Q_z_F[ei][idx]);
        }
    }
    
    /*
     * Calulate the weights beta for hybridization.
     */
    
    const Real u_z_diff = u_z_F - u_z_B;
    const Real v_z_diff = v_z_F - v_z_B;
    const Real w_z_diff = w_z_F - w_z_B;
    const Real vel_mag = std::sqrt(u_z_diff*u_z_diff + v_z_diff*v_z_diff + w_z_diff*w_z_diff);
    
    Real alpha_1, alpha_2;
    if (vel_mag < EPSILON)
    {
        alpha_1 = Real(1);
        alpha_2 = Real(0);
    }
    else
    {
        alpha_1 = std::abs(w_z_diff)/vel_mag;
        alpha_2 = std::sqrt(Real(1) - alpha_1*alpha_1);
    }
    
    const Real beta_1 = Real(1)/Real(2)*(Real(1) + alpha_1/(alpha_1 + alpha_2));
    const Real beta_2 = Real(1) - beta_1;
    
    F_z[0][idx_flux] = beta_1*F_z_HLLC[0] + beta_2*F_z_HLL[0];
    F_z[1][idx_flux] = beta_1*F_z_HLLC[1] + beta_2*F_z_HLL[1];
    F_z[2][idx_flux] = beta_1*F_z_HLLC[2] + beta_2*F_z_HLL[2];
    F_z[3][idx_flux] = F_z_HLLC[3];
    F_z[4][idx_flux] = F_z_HLLC[4];
}


/*
 * Compute the local convective flux in the x-direction from primitive variables with
 * 1D HLLC-HLL Riemann solver.
 */
static inline __attribute__((always_inline)) void
computeLocalConvectiveFluxInXDirectionFromPrimitiveVariablesHLLC_HLL1D(
    Real** F_x,
    Real** V_x_L,
    Real** V_x_R,
    Real* c_x_L,
    Real* c_x_R,
    Real* epsilon_x_L,
    Real* epsilon_x_R,
    Real& s_x_minus,
    Real& s_x_plus,
    Real& s_x_star,
    Real& Chi_x_star_LR,
    const int& idx_flux,
    const int& idx)
{
    const Real u_x_average = Real(1)/Real(2)*(V_x_L[1][idx] + V_x_R[1][idx]);
    const Real c_x_average = Real(1)/Real(2)*(c_x_L[idx] + c_x_R[idx]);
    
    const Real s_x_L = std::min(u_x_average - c_x_average, V_x_L[1][idx] - c_x_L[idx]);
    const Real s_x_R = std::max(u_x_average + c_x_average, V_x_R[1][idx] + c_x_R[idx]);
    
    s_x_minus = std::min(Real(0), s_x_L);
    s_x_plus  = std::max(Real(0), s_x_R);
    
    s_x_star = (V_x_R[2][idx] - V_x_L[2][idx] + V_x_L[0][idx]*V_x_L[1][idx]*(s_x_L - V_x_L[1][idx]) -
        V_x_R[0][idx]*V_x_R[1][idx]*(s_x_R - V_x_R[1][idx]))/
        (V_x_L[0][idx]*(s_x_L - V_x_L[1][idx]) - V_x_R[0][idx]*(s_x_R - V_x_R[1][idx]));
    
    Real Q_x_LR[3];
    Real Q_x_star_LR[3];
    Real F_x_LR[3];
    
    if (s_x_star > Real(0))
    {
        Chi_x_star_LR = (s_x_L - V_x_L[1][idx])/(s_x_L - s_x_star);
        
        Q_x_LR[0] = V_x_L[0][idx];
        Q_x_LR[1] = V_x_L[0][idx]*V_x_L[1][idx];
        Q_x_LR[2] = V_x_L[0][idx]*(epsilon_x_L[idx] + Real(1)/Real(2)*V_x_L[1][idx]*V_x_L[1][idx]);
        
        Q_x_star_LR[0] = Chi_x_star_LR*V_x_L[0][idx];
        Q_x_star_LR[1] = Chi_x_star_LR*V_x_L[0][idx]*s_x_star;
        Q_x_star_LR[2] = Chi_x_star_LR*(Q_x_LR[2] + (s_x_star - V_x_L[1][idx])*(V_x_L[0][idx]*s_x_star +
            V_x_L[2][idx]/(s_x_L - V_x_L[1][idx])));
        
        F_x_LR[0] = Q_x_LR[1];
        F_x_LR[1] = Q_x_LR[1]*V_x_L[1][idx] + V_x_L[2][idx];
        F_x_LR[2] = V_x_L[1][idx]*(Q_x_LR[2] + V_x_L[2][idx]);
        
        for (int ei = 0; ei < 3; ei++)
        {
            F_x[ei][idx_flux] = F_x_LR[ei] + s_x_minus*(Q_x_star_LR[ei] - Q_x_LR[ei]);
        }
    }
    else
    {
        Chi_x_star_LR = (s_x_R - V_x_R[1][idx])/(s_x_R - s_x_star);
        
        Q_x_LR[0] = V_x_R[0][idx];
        Q_x_LR[1] = V_x_R[0][idx]*V_x_R[1][idx];
        Q_x_LR[2] = V_x_R[0][idx]*(epsilon_x_R[idx] + Real(1)/Real(2)*V_x_R[1][idx]*V_x_R[1][idx]);
        
        Q_x_star_LR[0] = Chi_x_star_LR*V_x_R[0][idx];
        Q_x_star_LR[1] = Chi_x_star_LR*V_x_R[0][idx]*s_x_star;
        Q_x_star_LR[2] = Chi_x_star_LR*(Q_x_LR[2] + (s_x_star - V_x_R[1][idx])*(V_x_R[0][idx]*s_x_star +
            V_x_R[2][idx]/(s_x_R - V_x_R[1][idx])));
        
        F_x_LR[0] = Q_x_LR[1];
        F_x_LR[1] = Q_x_LR[1]*V_x_R[1][idx] + V_x_R[2][idx];
        F_x_LR[2] = V_x_R[1][idx]*(Q_x_LR[2] + V_x_R[2][idx]);
        
        for (int ei = 0; ei < 3; ei++)
        {
            F_x[ei][idx_flux] = F_x_LR[ei] + s_x_plus*(Q_x_star_LR[ei] - Q_x_LR[ei]);
        }
    }
}


/*
 * Compute the local convective flux in the x-direction from primitive variables with
 * 2D HLLC-HLL Riemann solver.
 */
static inline __attribute__((always_inline)) void
computeLocalConvectiveFluxInXDirectionFromPrimitiveVariablesHLLC_HLL2D(
    Real** F_x,
    Real** V_x_L,
    Real** V_x_R,
    Real* c_x_L,
    Real* c_x_R,
    Real* epsilon_x_L,
    Real* epsilon_x_R,
    Real& s_x_minus,
    Real& s_x_plus,
    Real& s_x_star,
    Real& Chi_x_star_LR,
    const int& idx_flux,
    const int& idx)
{
    const Real u_x_average = Real(1)/Real(2)*(V_x_L[1][idx] + V_x_R[1][idx]);
    const Real c_x_average = Real(1)/Real(2)*(c_x_L[idx] + c_x_R[idx]);
    
    const Real s_x_L = std::min(u_x_average - c_x_average, V_x_L[1][idx] - c_x_L[idx]);
    const Real s_x_R = std::max(u_x_average + c_x_average, V_x_R[1][idx] + c_x_R[idx]);
    
    s_x_minus = std::min(Real(0), s_x_L);
    s_x_plus  = std::max(Real(0), s_x_R);
    
    s_x_star = (V_x_R[3][idx] - V_x_L[3][idx] + V_x_L[0][idx]*V_x_L[1][idx]*(s_x_L - V_x_L[1][idx]) -
        V_x_R[0][idx]*V_x_R[1][idx]*(s_x_R - V_x_R[1][idx]))/
        (V_x_L[0][idx]*(s_x_L - V_x_L[1][idx]) - V_x_R[0][idx]*(s_x_R - V_x_R[1][idx]));
    
    Real Q_x_L[4];
    Real Q_x_R[4];
    Real F_x_L[4];
    Real F_x_R[4];
    Real F_x_HLL[2];
    Real F_x_HLLC[4];
    Real Q_x_star_LR[4];
    
    Q_x_L[0] = V_x_L[0][idx];
    Q_x_L[1] = V_x_L[0][idx]*V_x_L[1][idx];
    Q_x_L[2] = V_x_L[0][idx]*V_x_L[2][idx];
    Q_x_L[3] = V_x_L[0][idx]*(epsilon_x_L[idx] + Real(1)/Real(2)*(V_x_L[1][idx]*V_x_L[1][idx]
        + V_x_L[2][idx]*V_x_L[2][idx]));
    
    Q_x_R[0] = V_x_R[0][idx];
    Q_x_R[1] = V_x_R[0][idx]*V_x_R[1][idx];
    Q_x_R[2] = V_x_R[0][idx]*V_x_R[2][idx];
    Q_x_R[3] = V_x_R[0][idx]*(epsilon_x_R[idx] + Real(1)/Real(2)*(V_x_R[1][idx]*V_x_R[1][idx]
        + V_x_R[2][idx]*V_x_R[2][idx]));
    
    F_x_L[0] = Q_x_L[1];
    F_x_L[1] = Q_x_L[1]*V_x_L[1][idx] + V_x_L[3][idx];
    F_x_L[2] = Q_x_L[1]*V_x_L[2][idx];
    F_x_L[3] = V_x_L[1][idx]*(Q_x_L[3] + V_x_L[3][idx]);
    
    F_x_R[0] = Q_x_R[1];
    F_x_R[1] = Q_x_R[1]*V_x_R[1][idx] + V_x_R[3][idx];
    F_x_R[2] = Q_x_R[1]*V_x_R[2][idx];
    F_x_R[3] = V_x_R[1][idx]*(Q_x_R[3] + V_x_R[3][idx]);
    
    F_x_HLL[0] = (s_x_R*F_x_L[0] - s_x_L*F_x_R[0] + s_x_R*s_x_L*(Q_x_R[0] - Q_x_L[0]))/
        (s_x_R - s_x_L);
    F_x_HLL[1] = (s_x_R*F_x_L[2] - s_x_L*F_x_R[2] + s_x_R*s_x_L*(Q_x_R[2] - Q_x_L[2]))/
        (s_x_R - s_x_L);
    
    if (s_x_L > Real(0))
    {
        F_x_HLL[0] = F_x_L[0];
        F_x_HLL[1] = F_x_L[2];
    }
    
    if (s_x_R < Real(0))
    {
        F_x_HLL[0] = F_x_R[0];
        F_x_HLL[1] = F_x_R[2];
    }
    
    if (s_x_star > Real(0))
    {
        Chi_x_star_LR = (s_x_L - V_x_L[1][idx])/(s_x_L - s_x_star);
        
        Q_x_star_LR[0] = Chi_x_star_LR*V_x_L[0][idx];
        Q_x_star_LR[1] = Chi_x_star_LR*V_x_L[0][idx]*s_x_star;
        Q_x_star_LR[2] = Chi_x_star_LR*Q_x_L[2];
        Q_x_star_LR[3] = Chi_x_star_LR*(Q_x_L[3] + (s_x_star - V_x_L[1][idx])*(V_x_L[0][idx]*s_x_star +
            V_x_L[3][idx]/(s_x_L - V_x_L[1][idx])));
        
        for (int ei = 0; ei < 4; ei++)
        {
            F_x_HLLC[ei] = F_x_L[ei] + s_x_minus*(Q_x_star_LR[ei] - Q_x_L[ei]);
        }
    }
    else
    {
        Chi_x_star_LR = (s_x_R - V_x_R[1][idx])/(s_x_R - s_x_star);
        
        Q_x_star_LR[0] = Chi_x_star_LR*V_x_R[0][idx];
        Q_x_star_LR[1] = Chi_x_star_LR*V_x_R[0][idx]*s_x_star;
        Q_x_star_LR[2] = Chi_x_star_LR*Q_x_R[2];
        Q_x_star_LR[3] = Chi_x_star_LR*(Q_x_R[3] + (s_x_star - V_x_R[1][idx])*(V_x_R[0][idx]*s_x_star +
            V_x_R[3][idx]/(s_x_R - V_x_R[1][idx])));
        
        for (int ei = 0; ei < 4; ei++)
        {
            F_x_HLLC[ei] = F_x_R[ei] + s_x_plus*(Q_x_star_LR[ei] - Q_x_R[ei]);
        }
    }
    
    /*
     * Calulate the weights beta for hybridization.
     */
    
    const Real u_x_diff = V_x_R[1][idx] - V_x_L[1][idx];
    const Real v_x_diff = V_x_R[2][idx] - V_x_L[2][idx];
    const Real vel_mag = std::sqrt(u_x_diff*u_x_diff + v_x_diff*v_x_diff);
    
    Real alpha_1, alpha_2;
    if (vel_mag < EPSILON)
    {
        alpha_1 = Real(1);
        alpha_2 = Real(0);
    }
    else
    {
        alpha_1 = std::abs(u_x_diff)/vel_mag;
        alpha_2 = std::sqrt(Real(1) - alpha_1*alpha_1);
    }
    
    const Real beta_1 = Real(1)/Real(2)*(Real(1) + alpha_1/(alpha_1 + alpha_2));
    const Real beta_2 = Real(1) - beta_1;
    
    F_x[0][idx_flux] = beta_1*F_x_HLLC[0] + beta_2*F_x_HLL[0];
    F_x[1][idx_flux] = F_x_HLLC[1];
    F_x[2][idx_flux] = beta_1*F_x_HLLC[2] + beta_2*F_x_HLL[1];
    F_x[3][idx_flux] = F_x_HLLC[3];
}


/*
 * Compute the local convective flux in the x-direction from primitive variables with
 * 3D HLLC-HLL Riemann solver.
 */
static inline __attribute__((always_inline)) void
computeLocalConvectiveFluxInXDirectionFromPrimitiveVariablesHLLC_HLL3D(
    Real** F_x,
    Real** V_x_L,
    Real** V_x_R,
    Real* c_x_L,
    Real* c_x_R,
    Real* epsilon_x_L,
    Real* epsilon_x_R,
    Real& s_x_minus,
    Real& s_x_plus,
    Real& s_x_star,
    Real& Chi_x_star_LR,
    const int& idx_flux,
    const int& idx)
{
    const Real u_x_average = Real(1)/Real(2)*(V_x_L[1][idx] + V_x_R[1][idx]);
    const Real c_x_average = Real(1)/Real(2)*(c_x_L[idx] + c_x_R[idx]);
    
    const Real s_x_L = std::min(u_x_average - c_x_average, V_x_L[1][idx] - c_x_L[idx]);
    const Real s_x_R = std::max(u_x_average + c_x_average, V_x_R[1][idx] + c_x_R[idx]);
    
    s_x_minus = std::min(Real(0), s_x_L);
    s_x_plus  = std::max(Real(0), s_x_R);
    
    s_x_star = (V_x_R[4][idx] - V_x_L[4][idx] + V_x_L[0][idx]*V_x_L[1][idx]*(s_x_L - V_x_L[1][idx]) -
        V_x_R[0][idx]*V_x_R[1][idx]*(s_x_R - V_x_R[1][idx]))/
        (V_x_L[0][idx]*(s_x_L - V_x_L[1][idx]) - V_x_R[0][idx]*(s_x_R - V_x_R[1][idx]));
    
    Real Q_x_L[5];
    Real Q_x_R[5];
    Real F_x_L[5];
    Real F_x_R[5];
    Real F_x_HLL[3];
    Real F_x_HLLC[5];
    Real Q_x_star_LR[5];
    
    Q_x_L[0] = V_x_L[0][idx];
    Q_x_L[1] = V_x_L[0][idx]*V_x_L[1][idx];
    Q_x_L[2] = V_x_L[0][idx]*V_x_L[2][idx];
    Q_x_L[3] = V_x_L[0][idx]*V_x_L[3][idx];
    Q_x_L[4] = V_x_L[0][idx]*(epsilon_x_L[idx] + Real(1)/Real(2)*(V_x_L[1][idx]*V_x_L[1][idx]
        + V_x_L[2][idx]*V_x_L[2][idx] + V_x_L[3][idx]*V_x_L[3][idx]));
    
    Q_x_R[0] = V_x_R[0][idx];
    Q_x_R[1] = V_x_R[0][idx]*V_x_R[1][idx];
    Q_x_R[2] = V_x_R[0][idx]*V_x_R[2][idx];
    Q_x_R[3] = V_x_R[0][idx]*V_x_R[3][idx];
    Q_x_R[4] = V_x_R[0][idx]*(epsilon_x_R[idx] + Real(1)/Real(2)*(V_x_R[1][idx]*V_x_R[1][idx]
        + V_x_R[2][idx]*V_x_R[2][idx] + V_x_R[3][idx]*V_x_R[3][idx]));
    
    F_x_L[0] = Q_x_L[1];
    F_x_L[1] = Q_x_L[1]*V_x_L[1][idx] + V_x_L[4][idx];
    F_x_L[2] = Q_x_L[1]*V_x_L[2][idx];
    F_x_L[3] = Q_x_L[1]*V_x_L[3][idx];
    F_x_L[4] = V_x_L[1][idx]*(Q_x_L[4] + V_x_L[4][idx]);
    
    F_x_R[0] = Q_x_R[1];
    F_x_R[1] = Q_x_R[1]*V_x_R[1][idx] + V_x_R[4][idx];
    F_x_R[2] = Q_x_R[1]*V_x_R[2][idx];
    F_x_R[3] = Q_x_R[1]*V_x_R[3][idx];
    F_x_R[4] = V_x_R[1][idx]*(Q_x_R[4] + V_x_R[4][idx]);
    
    F_x_HLL[0] = (s_x_R*F_x_L[0] - s_x_L*F_x_R[0] + s_x_R*s_x_L*(Q_x_R[0] - Q_x_L[0]))/
        (s_x_R - s_x_L);
    F_x_HLL[1] = (s_x_R*F_x_L[2] - s_x_L*F_x_R[2] + s_x_R*s_x_L*(Q_x_R[2] - Q_x_L[2]))/
        (s_x_R - s_x_L);
    F_x_HLL[2] = (s_x_R*F_x_L[3] - s_x_L*F_x_R[3] + s_x_R*s_x_L*(Q_x_R[3] - Q_x_L[3]))/
        (s_x_R - s_x_L);
    
    if (s_x_L > Real(0))
    {
        F_x_HLL[0] = F_x_L[0];
        F_x_HLL[1] = F_x_L[2];
        F_x_HLL[2] = F_x_L[3];
    }
    
    if (s_x_R < Real(0))
    {
        F_x_HLL[0] = F_x_R[0];
        F_x_HLL[1] = F_x_R[2];
        F_x_HLL[2] = F_x_R[3];
    }
    
    if (s_x_star > Real(0))
    {
        Chi_x_star_LR = (s_x_L - V_x_L[1][idx])/(s_x_L - s_x_star);
        
        Q_x_star_LR[0] = Chi_x_star_LR*V_x_L[0][idx];
        Q_x_star_LR[1] = Chi_x_star_LR*V_x_L[0][idx]*s_x_star;
        Q_x_star_LR[2] = Chi_x_star_LR*Q_x_L[2];
        Q_x_star_LR[3] = Chi_x_star_LR*Q_x_L[3];
        Q_x_star_LR[4] = Chi_x_star_LR*(Q_x_L[4] + (s_x_star - V_x_L[1][idx])*(V_x_L[0][idx]*s_x_star +
            V_x_L[4][idx]/(s_x_L - V_x_L[1][idx])));
        
        for (int ei = 0; ei < 5; ei++)
        {
            F_x_HLLC[ei] = F_x_L[ei] + s_x_minus*(Q_x_star_LR[ei] - Q_x_L[ei]);
        }
    }
    else
    {
        Chi_x_star_LR = (s_x_R - V_x_R[1][idx])/(s_x_R - s_x_star);
        
        Q_x_star_LR[0] = Chi_x_star_LR*V_x_R[0][idx];
        Q_x_star_LR[1] = Chi_x_star_LR*V_x_R[0][idx]*s_x_star;
        Q_x_star_LR[2] = Chi_x_star_LR*Q_x_R[2];
        Q_x_star_LR[3] = Chi_x_star_LR*Q_x_R[3];
        Q_x_star_LR[4] = Chi_x_star_LR*(Q_x_R[4] + (s_x_star - V_x_R[1][idx])*(V_x_R[0][idx]*s_x_star +
            V_x_R[4][idx]/(s_x_R - V_x_R[1][idx])));
        
        for (int ei = 0; ei < 5; ei++)
        {
            F_x_HLLC[ei] = F_x_R[ei] + s_x_plus*(Q_x_star_LR[ei] - Q_x_R[ei]);
        }
    }
    
    /*
     * Calulate the weights beta for hybridization.
     */
    
    const Real u_x_diff = V_x_R[1][idx] - V_x_L[1][idx];
    const Real v_x_diff = V_x_R[2][idx] - V_x_L[2][idx];
    const Real w_x_diff = V_x_R[3][idx] - V_x_L[3][idx];
    const Real vel_mag = std::sqrt(u_x_diff*u_x_diff + v_x_diff*v_x_diff + w_x_diff*w_x_diff);
    
    Real alpha_1, alpha_2;
    if (vel_mag < EPSILON)
    {
        alpha_1 = Real(1);
        alpha_2 = Real(0);
    }
    else
    {
        alpha_1 = std::abs(u_x_diff)/vel_mag;
        alpha_2 = std::sqrt(Real(1) - alpha_1*alpha_1);
    }
    
    const Real beta_1 = Real(1)/Real(2)*(Real(1) + alpha_1/(alpha_1 + alpha_2));
    const Real beta_2 = Real(1) - beta_1;
    
    F_x[0][idx_flux] = beta_1*F_x_HLLC[0] + beta_2*F_x_HLL[0];
    F_x[1][idx_flux] = F_x_HLLC[1];
    F_x[2][idx_flux] = beta_1*F_x_HLLC[2] + beta_2*F_x_HLL[1];
    F_x[3][idx_flux] = beta_1*F_x_HLLC[3] + beta_2*F_x_HLL[2];
    F_x[4][idx_flux] = F_x_HLLC[4];
}


/*
 * Compute the local convective flux in the y-direction from primitive variables with
 * 2D HLLC-HLL Riemann solver.
 */
static inline __attribute__((always_inline)) void
computeLocalConvectiveFluxInYDirectionFromPrimitiveVariablesHLLC_HLL2D(
    Real** F_y,
    Real** V_y_B,
    Real** V_y_T,
    Real* c_y_B,
    Real* c_y_T,
    Real* epsilon_y_B,
    Real* epsilon_y_T,
    Real& s_y_minus,
    Real& s_y_plus,
    Real& s_y_star,
    Real& Chi_y_star_BT,
    const int& idx_flux,
    const int& idx)
{
    const Real v_y_average = Real(1)/Real(2)*(V_y_B[2][idx] + V_y_T[2][idx]);
    const Real c_y_average = Real(1)/Real(2)*(c_y_B[idx] + c_y_T[idx]);
    
    const Real s_y_B = std::min(v_y_average - c_y_average, V_y_B[2][idx] - c_y_B[idx]);
    const Real s_y_T = std::max(v_y_average + c_y_average, V_y_T[2][idx] + c_y_T[idx]);
    
    s_y_minus = std::min(Real(0), s_y_B);
    s_y_plus  = std::max(Real(0), s_y_T);
    
    s_y_star = (V_y_T[3][idx] - V_y_B[3][idx] + V_y_B[0][idx]*V_y_B[2][idx]*(s_y_B - V_y_B[2][idx]) -
        V_y_T[0][idx]*V_y_T[2][idx]*(s_y_T - V_y_T[2][idx]))/
        (V_y_B[0][idx]*(s_y_B - V_y_B[2][idx]) - V_y_T[0][idx]*(s_y_T - V_y_T[2][idx]));
    
    Real Q_y_B[4];
    Real Q_y_T[4];
    Real F_y_B[4];
    Real F_y_T[4];
    Real F_y_HLL[2];
    Real F_y_HLLC[4];
    Real Q_y_star_BT[4];
    
    Q_y_B[0] = V_y_B[0][idx];
    Q_y_B[1] = V_y_B[0][idx]*V_y_B[1][idx];
    Q_y_B[2] = V_y_B[0][idx]*V_y_B[2][idx];
    Q_y_B[3] = V_y_B[0][idx]*(epsilon_y_B[idx] + Real(1)/Real(2)*(V_y_B[1][idx]*V_y_B[1][idx]
        + V_y_B[2][idx]*V_y_B[2][idx]));
    
    Q_y_T[0] = V_y_T[0][idx];
    Q_y_T[1] = V_y_T[0][idx]*V_y_T[1][idx];
    Q_y_T[2] = V_y_T[0][idx]*V_y_T[2][idx];
    Q_y_T[3] = V_y_T[0][idx]*(epsilon_y_T[idx] + Real(1)/Real(2)*(V_y_T[1][idx]*V_y_T[1][idx]
        + V_y_T[2][idx]*V_y_T[2][idx]));
    
    F_y_B[0] = Q_y_B[2];
    F_y_B[1] = Q_y_B[2]*V_y_B[1][idx];
    F_y_B[2] = Q_y_B[2]*V_y_B[2][idx] + V_y_B[3][idx];
    F_y_B[3] = V_y_B[2][idx]*(Q_y_B[3] + V_y_B[3][idx]);
    
    F_y_T[0] = Q_y_T[2];
    F_y_T[1] = Q_y_T[2]*V_y_T[1][idx];
    F_y_T[2] = Q_y_T[2]*V_y_T[2][idx] + V_y_T[3][idx];
    F_y_T[3] = V_y_T[2][idx]*(Q_y_T[3] + V_y_T[3][idx]);
    
    F_y_HLL[0] = (s_y_T*F_y_B[0] - s_y_B*F_y_T[0] + s_y_T*s_y_B*(Q_y_T[0] - Q_y_B[0]))/
        (s_y_T - s_y_B);
    F_y_HLL[1] = (s_y_T*F_y_B[1] - s_y_B*F_y_T[1] + s_y_T*s_y_B*(Q_y_T[1] - Q_y_B[1]))/
        (s_y_T - s_y_B);
    
    if (s_y_B > Real(0))
    {
        F_y_HLL[0] = F_y_B[0];
        F_y_HLL[1] = F_y_B[1];
    }
    
    if (s_y_T < Real(0))
    {
        F_y_HLL[0] = F_y_T[0];
        F_y_HLL[1] = F_y_T[1];
    }
    
    if (s_y_star > Real(0))
    {
        Chi_y_star_BT = (s_y_B - V_y_B[2][idx])/(s_y_B - s_y_star);
        
        Q_y_star_BT[0] = Chi_y_star_BT*V_y_B[0][idx];
        Q_y_star_BT[1] = Chi_y_star_BT*Q_y_B[1];
        Q_y_star_BT[2] = Chi_y_star_BT*V_y_B[0][idx]*s_y_star;
        Q_y_star_BT[3] = Chi_y_star_BT*(Q_y_B[3] + (s_y_star - V_y_B[2][idx])*(V_y_B[0][idx]*s_y_star +
            V_y_B[3][idx]/(s_y_B - V_y_B[2][idx])));
        
        for (int ei = 0; ei < 4; ei++)
        {
            F_y_HLLC[ei] = F_y_B[ei] + s_y_minus*(Q_y_star_BT[ei] - Q_y_B[ei]);
        }
    }
    else
    {
        Chi_y_star_BT = (s_y_T - V_y_T[2][idx])/(s_y_T - s_y_star);
        
        Q_y_star_BT[0] = Chi_y_star_BT*V_y_T[0][idx];
        Q_y_star_BT[1] = Chi_y_star_BT*Q_y_T[1];
        Q_y_star_BT[2] = Chi_y_star_BT*V_y_T[0][idx]*s_y_star;
        Q_y_star_BT[3] = Chi_y_star_BT*(Q_y_T[3] + (s_y_star - V_y_T[2][idx])*(V_y_T[0][idx]*s_y_star +
            V_y_T[3][idx]/(s_y_T - V_y_T[2][idx])));
        
        for (int ei = 0; ei < 4; ei++)
        {
            F_y_HLLC[ei] = F_y_T[ei] + s_y_plus*(Q_y_star_BT[ei] - Q_y_T[ei]);
        }
    }
    
    /*
     * Calulate the weights beta for hybridization.
     */
    
    const Real u_y_diff = V_y_T[1][idx] - V_y_B[1][idx];
    const Real v_y_diff = V_y_T[2][idx] - V_y_B[2][idx];
    const Real vel_mag = std::sqrt(u_y_diff*u_y_diff + v_y_diff*v_y_diff);
    
    Real alpha_1, alpha_2;
    if (vel_mag < EPSILON)
    {
        alpha_1 = Real(1);
        alpha_2 = Real(0);
    }
    else
    {
        alpha_1 = std::abs(v_y_diff)/vel_mag;
        alpha_2 = std::sqrt(Real(1) - alpha_1*alpha_1);
    }
    
    const Real beta_1 = Real(1)/Real(2)*(Real(1) + alpha_1/(alpha_1 + alpha_2));
    const Real beta_2 = Real(1) - beta_1;
    
    F_y[0][idx_flux] = beta_1*F_y_HLLC[0] + beta_2*F_y_HLL[0];
    F_y[1][idx_flux] = beta_1*F_y_HLLC[1] + beta_2*F_y_HLL[1];
    F_y[2][idx_flux] = F_y_HLLC[2];
    F_y[3][idx_flux] = F_y_HLLC[3];
}


/*
 * Compute the local convective flux in the y-direction from primitive variables with
 * 3D HLLC-HLL Riemann solver.
 */
static inline __attribute__((always_inline)) void
computeLocalConvectiveFluxInYDirectionFromPrimitiveVariablesHLLC_HLL3D(
    Real** F_y,
    Real** V_y_B,
    Real** V_y_T,
    Real* c_y_B,
    Real* c_y_T,
    Real* epsilon_y_B,
    Real* epsilon_y_T,
    Real& s_y_minus,
    Real& s_y_plus,
    Real& s_y_star,
    Real& Chi_y_star_BT,
    const int& idx_flux,
    const int& idx)
{
    const Real v_y_average = Real(1)/Real(2)*(V_y_B[2][idx] + V_y_T[2][idx]);
    const Real c_y_average = Real(1)/Real(2)*(c_y_B[idx] + c_y_T[idx]);
    
    const Real s_y_B = std::min(v_y_average - c_y_average, V_y_B[2][idx] - c_y_B[idx]);
    const Real s_y_T = std::max(v_y_average + c_y_average, V_y_T[2][idx] + c_y_T[idx]);
    
    s_y_minus = std::min(Real(0), s_y_B);
    s_y_plus  = std::max(Real(0), s_y_T);
    
    s_y_star = (V_y_T[4][idx] - V_y_B[4][idx] + V_y_B[0][idx]*V_y_B[2][idx]*(s_y_B - V_y_B[2][idx]) -
        V_y_T[0][idx]*V_y_T[2][idx]*(s_y_T - V_y_T[2][idx]))/
        (V_y_B[0][idx]*(s_y_B - V_y_B[2][idx]) - V_y_T[0][idx]*(s_y_T - V_y_T[2][idx]));
    
    Real Q_y_B[5];
    Real Q_y_T[5];
    Real F_y_B[5];
    Real F_y_T[5];
    Real F_y_HLL[3];
    Real F_y_HLLC[5];
    Real Q_y_star_BT[5];
    
    Q_y_B[0] = V_y_B[0][idx];
    Q_y_B[1] = V_y_B[0][idx]*V_y_B[1][idx];
    Q_y_B[2] = V_y_B[0][idx]*V_y_B[2][idx];
    Q_y_B[3] = V_y_B[0][idx]*V_y_B[3][idx];
    Q_y_B[4] = V_y_B[0][idx]*(epsilon_y_B[idx] + Real(1)/Real(2)*(V_y_B[1][idx]*V_y_B[1][idx]
        + V_y_B[2][idx]*V_y_B[2][idx] + V_y_B[3][idx]*V_y_B[3][idx]));
    
    Q_y_T[0] = V_y_T[0][idx];
    Q_y_T[1] = V_y_T[0][idx]*V_y_T[1][idx];
    Q_y_T[2] = V_y_T[0][idx]*V_y_T[2][idx];
    Q_y_T[3] = V_y_T[0][idx]*V_y_T[3][idx];
    Q_y_T[4] = V_y_T[0][idx]*(epsilon_y_T[idx] + Real(1)/Real(2)*(V_y_T[1][idx]*V_y_T[1][idx]
        + V_y_T[2][idx]*V_y_T[2][idx] + V_y_T[3][idx]*V_y_T[3][idx]));
    
    F_y_B[0] = Q_y_B[2];
    F_y_B[1] = Q_y_B[2]*V_y_B[1][idx];
    F_y_B[2] = Q_y_B[2]*V_y_B[2][idx] + V_y_B[4][idx];
    F_y_B[3] = Q_y_B[2]*V_y_B[3][idx];
    F_y_B[4] = V_y_B[2][idx]*(Q_y_B[4] + V_y_B[4][idx]);
    
    F_y_T[0] = Q_y_T[2];
    F_y_T[1] = Q_y_T[2]*V_y_T[1][idx];
    F_y_T[2] = Q_y_T[2]*V_y_T[2][idx] + V_y_T[4][idx];
    F_y_T[3] = Q_y_T[2]*V_y_T[3][idx];
    F_y_T[4] = V_y_T[2][idx]*(Q_y_T[4] + V_y_T[4][idx]);
    
    F_y_HLL[0] = (s_y_T*F_y_B[0] - s_y_B*F_y_T[0] + s_y_T*s_y_B*(Q_y_T[0] - Q_y_B[0]))/
        (s_y_T - s_y_B);
    F_y_HLL[1] = (s_y_T*F_y_B[1] - s_y_B*F_y_T[1] + s_y_T*s_y_B*(Q_y_T[1] - Q_y_B[1]))/
        (s_y_T - s_y_B);
    F_y_HLL[2] = (s_y_T*F_y_B[3] - s_y_B*F_y_T[3] + s_y_T*s_y_B*(Q_y_T[3] - Q_y_B[3]))/
        (s_y_T - s_y_B);
    
    if (s_y_B > Real(0))
    {
        F_y_HLL[0] = F_y_B[0];
        F_y_HLL[1] = F_y_B[1];
        F_y_HLL[2] = F_y_B[3];
    }
    
    if (s_y_T < Real(0))
    {
        F_y_HLL[0] = F_y_T[0];
        F_y_HLL[1] = F_y_T[1];
        F_y_HLL[2] = F_y_T[3];
    }
    
    if (s_y_star > Real(0))
    {
        Chi_y_star_BT = (s_y_B - V_y_B[2][idx])/(s_y_B - s_y_star);
        
        Q_y_star_BT[0] = Chi_y_star_BT*V_y_B[0][idx];
        Q_y_star_BT[1] = Chi_y_star_BT*Q_y_B[1];
        Q_y_star_BT[2] = Chi_y_star_BT*V_y_B[0][idx]*s_y_star;
        Q_y_star_BT[3] = Chi_y_star_BT*Q_y_B[3];
        Q_y_star_BT[4] = Chi_y_star_BT*(Q_y_B[4] + (s_y_star - V_y_B[2][idx])*(V_y_B[0][idx]*s_y_star +
            V_y_B[4][idx]/(s_y_B - V_y_B[2][idx])));
        
        for (int ei = 0; ei < 5; ei++)
        {
            F_y_HLLC[ei] = F_y_B[ei] + s_y_minus*(Q_y_star_BT[ei] - Q_y_B[ei]);
        }
    }
    else
    {
        Chi_y_star_BT = (s_y_T - V_y_T[2][idx])/(s_y_T - s_y_star);
        
        Q_y_star_BT[0] = Chi_y_star_BT*V_y_T[0][idx];
        Q_y_star_BT[1] = Chi_y_star_BT*Q_y_T[1];
        Q_y_star_BT[2] = Chi_y_star_BT*V_y_T[0][idx]*s_y_star;
        Q_y_star_BT[3] = Chi_y_star_BT*Q_y_T[3];
        Q_y_star_BT[4] = Chi_y_star_BT*(Q_y_T[4] + (s_y_star - V_y_T[2][idx])*(V_y_T[0][idx]*s_y_star +
            V_y_T[4][idx]/(s_y_T - V_y_T[2][idx])));
        
        for (int ei = 0; ei < 5; ei++)
        {
            F_y_HLLC[ei] = F_y_T[ei] + s_y_plus*(Q_y_star_BT[ei] - Q_y_T[ei]);
        }
    }
    
    /*
     * Calulate the weights beta for hybridization.
     */
    
    const Real u_y_diff = V_y_T[1][idx] - V_y_B[1][idx];
    const Real v_y_diff = V_y_T[2][idx] - V_y_B[2][idx];
    const Real w_y_diff = V_y_T[3][idx] - V_y_B[3][idx];
    const Real vel_mag = std::sqrt(u_y_diff*u_y_diff + v_y_diff*v_y_diff + w_y_diff*w_y_diff);
    
    Real alpha_1, alpha_2;
    if (vel_mag < EPSILON)
    {
        alpha_1 = Real(1);
        alpha_2 = Real(0);
    }
    else
    {
        alpha_1 = std::abs(v_y_diff)/vel_mag;
        alpha_2 = std::sqrt(Real(1) - alpha_1*alpha_1);
    }
    
    const Real beta_1 = Real(1)/Real(2)*(Real(1) + alpha_1/(alpha_1 + alpha_2));
    const Real beta_2 = Real(1) - beta_1;
    
    F_y[0][idx_flux] = beta_1*F_y_HLLC[0] + beta_2*F_y_HLL[0];
    F_y[1][idx_flux] = beta_1*F_y_HLLC[1] + beta_2*F_y_HLL[1];
    F_y[2][idx_flux] = F_y_HLLC[2];
    F_y[3][idx_flux] = beta_1*F_y_HLLC[3] + beta_2*F_y_HLL[2];
    F_y[4][idx_flux] = F_y_HLLC[4];
}


/*
 * Compute the local convective flux in the z-direction from primitive variables with
 * 3D HLLC-HLL Riemann solver.
 */
static inline __attribute__((always_inline)) void
computeLocalConvectiveFluxInZDirectionFromPrimitiveVariablesHLLC_HLL3D(
    Real** F_z,
    Real** V_z_B,
    Real** V_z_F,
    Real* c_z_B,
    Real* c_z_F,
    Real* epsilon_z_B,
    Real* epsilon_z_F,
    Real& s_z_minus,
    Real& s_z_plus,
    Real& s_z_star,
    Real& Chi_z_star_BF,
    const int& idx_flux,
    const int& idx)
{
    const Real w_z_average = Real(1)/Real(2)*(V_z_B[3][idx] + V_z_F[3][idx]);
    const Real c_z_average = Real(1)/Real(2)*(c_z_B[idx] + c_z_F[idx]);
    
    const Real s_z_B = std::min(w_z_average - c_z_average, V_z_B[3][idx] - c_z_B[idx]);
    const Real s_z_F = std::max(w_z_average + c_z_average, V_z_F[3][idx] + c_z_F[idx]);
    
    s_z_minus = std::min(Real(0), s_z_B);
    s_z_plus  = std::max(Real(0), s_z_F);
    
    s_z_star = (V_z_F[4][idx] - V_z_B[4][idx] + V_z_B[0][idx]*V_z_B[3][idx]*(s_z_B - V_z_B[3][idx]) -
        V_z_F[0][idx]*V_z_F[3][idx]*(s_z_F - V_z_F[3][idx]))/
        (V_z_B[0][idx]*(s_z_B - V_z_B[3][idx]) - V_z_F[0][idx]*(s_z_F - V_z_F[3][idx]));
    
    Real Q_z_B[5];
    Real Q_z_F[5];
    Real F_z_B[5];
    Real F_z_F[5];
    Real F_z_HLL[3];
    Real F_z_HLLC[5];
    Real Q_z_star_BF[5];
    
    Q_z_B[0] = V_z_B[0][idx];
    Q_z_B[1] = V_z_B[0][idx]*V_z_B[1][idx];
    Q_z_B[2] = V_z_B[0][idx]*V_z_B[2][idx];
    Q_z_B[3] = V_z_B[0][idx]*V_z_B[3][idx];
    Q_z_B[4] = V_z_B[0][idx]*(epsilon_z_B[idx] + Real(1)/Real(2)*(V_z_B[1][idx]*V_z_B[1][idx]
        + V_z_B[2][idx]*V_z_B[2][idx] + V_z_B[3][idx]*V_z_B[3][idx]));
    
    Q_z_F[0] = V_z_F[0][idx];
    Q_z_F[1] = V_z_F[0][idx]*V_z_F[1][idx];
    Q_z_F[2] = V_z_F[0][idx]*V_z_F[2][idx];
    Q_z_F[3] = V_z_F[0][idx]*V_z_F[3][idx];
    Q_z_F[4] = V_z_F[0][idx]*(epsilon_z_F[idx] + Real(1)/Real(2)*(V_z_F[1][idx]*V_z_F[1][idx]
        + V_z_F[2][idx]*V_z_F[2][idx] + V_z_F[3][idx]*V_z_F[3][idx]));
    
    F_z_B[0] = Q_z_B[3];
    F_z_B[1] = Q_z_B[3]*V_z_B[1][idx];
    F_z_B[2] = Q_z_B[3]*V_z_B[2][idx];
    F_z_B[3] = Q_z_B[3]*V_z_B[3][idx] + V_z_B[4][idx];
    F_z_B[4] = V_z_B[3][idx]*(Q_z_B[4] + V_z_B[4][idx]);
    
    F_z_F[0] = Q_z_F[3];
    F_z_F[1] = Q_z_F[3]*V_z_F[1][idx];
    F_z_F[2] = Q_z_F[3]*V_z_F[2][idx];
    F_z_F[3] = Q_z_F[3]*V_z_F[3][idx] + V_z_F[4][idx];
    F_z_F[4] = V_z_F[3][idx]*(Q_z_F[4] + V_z_F[4][idx]);
    
    F_z_HLL[0] = (s_z_F*F_z_B[0] - s_z_B*F_z_F[0] + s_z_F*s_z_B*(Q_z_F[0] - Q_z_B[0]))/
        (s_z_F - s_z_B);
    F_z_HLL[1] = (s_z_F*F_z_B[1] - s_z_B*F_z_F[1] + s_z_F*s_z_B*(Q_z_F[1] - Q_z_B[1]))/
        (s_z_F - s_z_B);
    F_z_HLL[2] = (s_z_F*F_z_B[2] - s_z_B*F_z_F[2] + s_z_F*s_z_B*(Q_z_F[2] - Q_z_B[2]))/
        (s_z_F - s_z_B);
    
    if (s_z_B > Real(0))
    {
        F_z_HLL[0] = F_z_B[0];
        F_z_HLL[1] = F_z_B[1];
        F_z_HLL[2] = F_z_B[2];
    }
    
    if (s_z_F < Real(0))
    {
        F_z_HLL[0] = F_z_F[0];
        F_z_HLL[1] = F_z_F[1];
        F_z_HLL[2] = F_z_F[2];
    }
    
    if (s_z_star > Real(0))
    {
        Chi_z_star_BF = (s_z_B - V_z_B[3][idx])/(s_z_B - s_z_star);
        
        Q_z_star_BF[0] = Chi_z_star_BF*V_z_B[0][idx];
        Q_z_star_BF[1] = Chi_z_star_BF*Q_z_B[1];
        Q_z_star_BF[2] = Chi_z_star_BF*Q_z_B[2];
        Q_z_star_BF[3] = Chi_z_star_BF*V_z_B[0][idx]*s_z_star;
        Q_z_star_BF[4] = Chi_z_star_BF*(Q_z_B[4] + (s_z_star - V_z_B[3][idx])*(V_z_B[0][idx]*s_z_star +
            V_z_B[4][idx]/(s_z_B - V_z_B[3][idx])));
        
        for (int ei = 0; ei < 5; ei++)
        {
            F_z_HLLC[ei] = F_z_B[ei] + s_z_minus*(Q_z_star_BF[ei] - Q_z_B[ei]);
        }
    }
    else
    {
        Chi_z_star_BF = (s_z_F - V_z_F[3][idx])/(s_z_F - s_z_star);
        
        Q_z_star_BF[0] = Chi_z_star_BF*V_z_F[0][idx];
        Q_z_star_BF[1] = Chi_z_star_BF*Q_z_F[1];
        Q_z_star_BF[2] = Chi_z_star_BF*Q_z_F[2];
        Q_z_star_BF[3] = Chi_z_star_BF*V_z_F[0][idx]*s_z_star;
        Q_z_star_BF[4] = Chi_z_star_BF*(Q_z_F[4] + (s_z_star - V_z_F[3][idx])*(V_z_F[0][idx]*s_z_star +
            V_z_F[4][idx]/(s_z_F - V_z_F[3][idx])));
        
        for (int ei = 0; ei < 5; ei++)
        {
            F_z_HLLC[ei] = F_z_F[ei] + s_z_plus*(Q_z_star_BF[ei] - Q_z_F[ei]);
        }
    }
    
    /*
     * Calulate the weights beta for hybridization.
     */
    
    const Real u_z_diff = V_z_F[1][idx] - V_z_B[1][idx];
    const Real v_z_diff = V_z_F[2][idx] - V_z_B[2][idx];
    const Real w_z_diff = V_z_F[3][idx] - V_z_B[3][idx];
    const Real vel_mag = std::sqrt(u_z_diff*u_z_diff + v_z_diff*v_z_diff + w_z_diff*w_z_diff);
    
    Real alpha_1, alpha_2;
    if (vel_mag < EPSILON)
    {
        alpha_1 = Real(1);
        alpha_2 = Real(0);
    }
    else
    {
        alpha_1 = std::abs(w_z_diff)/vel_mag;
        alpha_2 = std::sqrt(Real(1) - alpha_1*alpha_1);
    }
    
    const Real beta_1 = Real(1)/Real(2)*(Real(1) + alpha_1/(alpha_1 + alpha_2));
    const Real beta_2 = Real(1) - beta_1;
    
    F_z[0][idx_flux] = beta_1*F_z_HLLC[0] + beta_2*F_z_HLL[0];
    F_z[1][idx_flux] = beta_1*F_z_HLLC[1] + beta_2*F_z_HLL[1];
    F_z[2][idx_flux] = beta_1*F_z_HLLC[2] + beta_2*F_z_HLL[2];
    F_z[3][idx_flux] = F_z_HLLC[3];
    F_z[4][idx_flux] = F_z_HLLC[4];
}


/*
 * Compute the convective flux and velocity in the x-direction from conservative variables with
 * HLLC-HLL Riemann solver.
 */
void
FlowModelRiemannSolverSingleSpecies::computeConvectiveFluxAndVelocityInXDirectionFromConservativeVariablesHLLC_HLL(
    HAMERS_SHARED_PTR<pdat::SideData<Real> > convective_flux,
    HAMERS_SHARED_PTR<pdat::SideData<Real> > velocity,
    const std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& conservative_variables_L,
    const std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& conservative_variables_R,
    const hier::Box& domain,
    bool compute_velocity) const
{
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    const int num_eqn = flow_model_tmp->getNumberOfEquations();
    
    // Get the box that covers the interior of patch.
    const hier::Box interior_box = convective_flux->getBox();
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_convective_flux = convective_flux->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_convective_flux =
        convective_flux->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_conservative_variables =
        conservative_variables_L[0]->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_conservative_variables =
        conservative_variables_L[0]->getGhostBox().numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        const hier::IntVector num_ghosts_min =
            hier::IntVector::min(num_ghosts_convective_flux, num_ghosts_conservative_variables);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(convective_flux->getGhostBox().contains(domain));
        TBOX_ASSERT(conservative_variables_L[0]->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
    /*
     * Get the equation of state mixing rules and the thermodynamic properties of the species.
     */
    
    const HAMERS_SHARED_PTR<EquationOfStateMixingRules> equation_of_state_mixing_rules =
        flow_model_tmp->getEquationOfStateMixingRules();
    
    const int num_thermo_properties = equation_of_state_mixing_rules->
        getNumberOfSpeciesThermodynamicProperties();
    
    std::vector<Real> thermo_properties;
    std::vector<Real*> thermo_properties_ptr;
    std::vector<const Real*> thermo_properties_const_ptr;
    
    thermo_properties.resize(num_thermo_properties);
    thermo_properties_ptr.reserve(num_thermo_properties);
    thermo_properties_const_ptr.reserve(num_thermo_properties);
    
    for (int ti = 0; ti < num_thermo_properties; ti++)
    {
        thermo_properties_ptr.push_back(&thermo_properties[ti]);
        thermo_properties_const_ptr.push_back(&thermo_properties[ti]);
    }
    
    equation_of_state_mixing_rules->getSpeciesThermodynamicProperties(
        thermo_properties_ptr,
        0);
    
    /*
     * Get the pointers to the side data of convective flux and conservative variables.
     */
    
    std::vector<Real*> F_x;
    F_x.reserve(num_eqn);
    for (int ei = 0; ei < num_eqn; ei++)
    {
        F_x.push_back(convective_flux->getPointer(0, ei));
    }
    
    std::vector<Real*> Q_x_L;
    std::vector<Real*> Q_x_R;
    Q_x_L.reserve(num_eqn);
    Q_x_R.reserve(num_eqn);
    for (int ei = 0; ei < num_eqn; ei++)
    {
        Q_x_L.push_back(conservative_variables_L[ei]->getPointer(0, 0));
        Q_x_R.push_back(conservative_variables_R[ei]->getPointer(0, 0));
    }
    
    /*
     * Allocate temporary data.
     */
    
    hier::IntVector direction_x = hier::IntVector::getZero(d_dim);
    direction_x[0] = 1;
    
    HAMERS_SHARED_PTR<pdat::SideData<Real> > internal_energy_x_L(
        new pdat::SideData<Real>(interior_box, 1, num_ghosts_conservative_variables,
            direction_x));

    HAMERS_SHARED_PTR<pdat::SideData<Real> > internal_energy_x_R(
        new pdat::SideData<Real>(interior_box, 1, num_ghosts_conservative_variables,
            direction_x));
    
    HAMERS_SHARED_PTR<pdat::SideData<Real> > pressure_x_L(
        new pdat::SideData<Real>(interior_box, 1, num_ghosts_conservative_variables,
            direction_x));
    
    HAMERS_SHARED_PTR<pdat::SideData<Real> > pressure_x_R(
        new pdat::SideData<Real>(interior_box, 1, num_ghosts_conservative_variables,
            direction_x));
    
    HAMERS_SHARED_PTR<pdat::SideData<Real> > sound_speed_x_L(
        new pdat::SideData<Real>(interior_box, 1, num_ghosts_conservative_variables,
            direction_x));
    
    HAMERS_SHARED_PTR<pdat::SideData<Real> > sound_speed_x_R(
        new pdat::SideData<Real>(interior_box, 1, num_ghosts_conservative_variables,
            direction_x));
    
    /*
     * Get the pointers to the temporary data.
     */
    
    Real* epsilon_x_L = internal_energy_x_L->getPointer(0, 0);
    Real* epsilon_x_R = internal_energy_x_R->getPointer(0, 0);
    
    Real* p_x_L = pressure_x_L->getPointer(0, 0);
    Real* p_x_R = pressure_x_R->getPointer(0, 0);
    
    Real* c_x_L = sound_speed_x_L->getPointer(0, 0);
    Real* c_x_R = sound_speed_x_R->getPointer(0, 0);
    
    Real u_x_L = Real(0);
    Real u_x_R = Real(0);
    
    Real s_x_minus = Real(0);
    Real s_x_plus  = Real(0);
    Real s_x_star  = Real(0);
    
    Real Chi_x_star_LR = Real(0);
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and numbers of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int num_ghosts_0_convective_flux = num_ghosts_convective_flux[0];
        const int num_ghosts_0_conservative_variables = num_ghosts_conservative_variables[0];
        
        /*
         * Compute the internal energy field.
         */
        
        HAMERS_PRAGMA_SIMD
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
        {
            // Compute the linear index.
            const int idx = i + num_ghosts_0_conservative_variables;
            
            epsilon_x_L[idx] = (Q_x_L[2][idx] -
                Real(1)/Real(2)*Q_x_L[1][idx]*Q_x_L[1][idx]/Q_x_L[0][idx])/Q_x_L[0][idx];
        }
        
        HAMERS_PRAGMA_SIMD
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
        {
            // Compute the linear index.
            const int idx = i + num_ghosts_0_conservative_variables;
            
            epsilon_x_R[idx] = (Q_x_R[2][idx] -
                Real(1)/Real(2)*Q_x_R[1][idx]*Q_x_R[1][idx]/Q_x_R[0][idx])/Q_x_R[0][idx];
        }
        
        flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computePressure(
                pressure_x_L,
                conservative_variables_L[0],
                internal_energy_x_L,
                thermo_properties_const_ptr,
                0,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computePressure(
                pressure_x_R,
                conservative_variables_R[0],
                internal_energy_x_R,
                thermo_properties_const_ptr,
                0,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computeSoundSpeed(
                sound_speed_x_L,
                conservative_variables_L[0],
                pressure_x_L,
                thermo_properties_const_ptr,
                0,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computeSoundSpeed(
                sound_speed_x_R,
                conservative_variables_R[0],
                pressure_x_R,
                thermo_properties_const_ptr,
                0,
                domain);
        
        if (compute_velocity)
        {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(velocity);
#endif
            
            /*
             * Get the number of ghost cells of velocity and the pointer to velocity.
             */
            
            const hier::IntVector num_ghosts_velocity = velocity->getGhostCellWidth();
            
            const int num_ghosts_0_velocity = num_ghosts_velocity[0];
            
            Real* u = velocity->getPointer(0, 0);
            
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
            {
                // Compute the linear indices.
                const int idx_flux = i + num_ghosts_0_convective_flux;
                const int idx_velocity = i + num_ghosts_0_velocity;
                const int idx = i + num_ghosts_0_conservative_variables;
                
                computeLocalConvectiveFluxInXDirectionFromConservativeVariablesHLLC_HLL1D(
                    F_x.data(),
                    Q_x_L.data(),
                    Q_x_R.data(),
                    p_x_L,
                    p_x_R,
                    c_x_L,
                    c_x_R,
                    u_x_L,
                    u_x_R,
                    s_x_minus,
                    s_x_plus,
                    s_x_star,
                    Chi_x_star_LR,
                    idx_flux,
                    idx);
                
                if (s_x_star > Real(0))
                {
                    u[idx_velocity] = u_x_L + s_x_minus*(Chi_x_star_LR - Real(1));
                }
                else
                {
                    u[idx_velocity] = u_x_R + s_x_plus*(Chi_x_star_LR - Real(1));
                }
            }
        }
        else
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
            {
                // Compute the linear indices.
                const int idx_flux = i + num_ghosts_0_convective_flux;
                const int idx = i + num_ghosts_0_conservative_variables;
                
                computeLocalConvectiveFluxInXDirectionFromConservativeVariablesHLLC_HLL1D(
                    F_x.data(),
                    Q_x_L.data(),
                    Q_x_R.data(),
                    p_x_L,
                    p_x_R,
                    c_x_L,
                    c_x_R,
                    u_x_L,
                    u_x_R,
                    s_x_minus,
                    s_x_plus,
                    s_x_star,
                    Chi_x_star_LR,
                    idx_flux,
                    idx);
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the local lower indices, numbers of cells in each dimension and numbers of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        
        const int num_ghosts_0_convective_flux = num_ghosts_convective_flux[0];
        const int num_ghosts_1_convective_flux = num_ghosts_convective_flux[1];
        const int ghostcell_dim_0_convective_flux = ghostcell_dims_convective_flux[0] + 1;
        
        const int num_ghosts_0_conservative_variables = num_ghosts_conservative_variables[0];
        const int num_ghosts_1_conservative_variables = num_ghosts_conservative_variables[1];
        const int ghostcell_dim_0_conservative_variables = ghostcell_dims_conservative_variables[0] + 1;
        
        /*
         * Compute the internal energy field.
         */
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
            {
                // Compute the linear index.
                const int idx = (i + num_ghosts_0_conservative_variables) +
                    (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables;
                
                epsilon_x_L[idx] = (Q_x_L[3][idx] -
                    Real(1)/Real(2)*(Q_x_L[1][idx]*Q_x_L[1][idx] + Q_x_L[2][idx]*Q_x_L[2][idx])/
                    Q_x_L[0][idx])/Q_x_L[0][idx];
            }
        }
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
            {
                // Compute the linear index.
                const int idx = (i + num_ghosts_0_conservative_variables) +
                    (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables;
                
                epsilon_x_R[idx] = (Q_x_R[3][idx] -
                    Real(1)/Real(2)*(Q_x_R[1][idx]*Q_x_R[1][idx] + Q_x_R[2][idx]*Q_x_R[2][idx])/
                    Q_x_R[0][idx])/Q_x_R[0][idx];
            }
        }
        
        flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computePressure(
                pressure_x_L,
                conservative_variables_L[0],
                internal_energy_x_L,
                thermo_properties_const_ptr,
                0,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computePressure(
                pressure_x_R,
                conservative_variables_R[0],
                internal_energy_x_R,
                thermo_properties_const_ptr,
                0,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computeSoundSpeed(
                sound_speed_x_L,
                conservative_variables_L[0],
                pressure_x_L,
                thermo_properties_const_ptr,
                0,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computeSoundSpeed(
                sound_speed_x_R,
                conservative_variables_R[0],
                pressure_x_R,
                thermo_properties_const_ptr,
                0,
                domain);
        
        if (compute_velocity)
        {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(velocity);
#endif
            
            /*
             * Get the numbers of cells in each dimension, number of ghost cells of velocity and
             * the pointer to velocity.
             */
            
            const hier::IntVector num_ghosts_velocity = velocity->getGhostCellWidth();
            const hier::IntVector ghostcell_dims_velocity = velocity->getGhostBox().numberCells();
            
            const int num_ghosts_0_velocity = num_ghosts_velocity[0];
            const int num_ghosts_1_velocity = num_ghosts_velocity[1];
            const int ghostcell_dim_0_velocity = ghostcell_dims_velocity[0] + 1;
            
            Real* u = velocity->getPointer(0, 0);
            
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                {
                    // Compute the linear indices.
                    const int idx_flux = (i + num_ghosts_0_convective_flux) +
                        (j + num_ghosts_1_convective_flux)*ghostcell_dim_0_convective_flux;
                    
                    const int idx_velocity = (i + num_ghosts_0_velocity) +
                        (j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                    
                    const int idx = (i + num_ghosts_0_conservative_variables) +
                        (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables;
                    
                    computeLocalConvectiveFluxInXDirectionFromConservativeVariablesHLLC_HLL2D(
                        F_x.data(),
                        Q_x_L.data(),
                        Q_x_R.data(),
                        p_x_L,
                        p_x_R,
                        c_x_L,
                        c_x_R,
                        u_x_L,
                        u_x_R,
                        s_x_minus,
                        s_x_plus,
                        s_x_star,
                        Chi_x_star_LR,
                        idx_flux,
                        idx);
                    
                    if (s_x_star > Real(0))
                    {
                        u[idx_velocity] = u_x_L + s_x_minus*(Chi_x_star_LR - Real(1));
                    }
                    else
                    {
                        u[idx_velocity] = u_x_R + s_x_plus*(Chi_x_star_LR - Real(1));
                    }
                }
            }
        }
        else
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                {
                    // Compute the linear indices.
                    const int idx_flux = (i + num_ghosts_0_convective_flux) +
                        (j + num_ghosts_1_convective_flux)*ghostcell_dim_0_convective_flux;
                    
                    const int idx = (i + num_ghosts_0_conservative_variables) +
                        (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables;
                    
                    computeLocalConvectiveFluxInXDirectionFromConservativeVariablesHLLC_HLL2D(
                        F_x.data(),
                        Q_x_L.data(),
                        Q_x_R.data(),
                        p_x_L,
                        p_x_R,
                        c_x_L,
                        c_x_R,
                        u_x_L,
                        u_x_R,
                        s_x_minus,
                        s_x_plus,
                        s_x_star,
                        Chi_x_star_LR,
                        idx_flux,
                        idx);
                }
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the local lower indices, numbers of cells in each dimension and numbers of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_lo_2 = domain_lo[2];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        const int domain_dim_2 = domain_dims[2];
        
        const int num_ghosts_0_convective_flux = num_ghosts_convective_flux[0];
        const int num_ghosts_1_convective_flux = num_ghosts_convective_flux[1];
        const int num_ghosts_2_convective_flux = num_ghosts_convective_flux[2];
        const int ghostcell_dim_0_convective_flux = ghostcell_dims_convective_flux[0] + 1;
        const int ghostcell_dim_1_convective_flux = ghostcell_dims_convective_flux[1];
        
        const int num_ghosts_0_conservative_variables = num_ghosts_conservative_variables[0];
        const int num_ghosts_1_conservative_variables = num_ghosts_conservative_variables[1];
        const int num_ghosts_2_conservative_variables = num_ghosts_conservative_variables[2];
        const int ghostcell_dim_0_conservative_variables = ghostcell_dims_conservative_variables[0] + 1;
        const int ghostcell_dim_1_conservative_variables = ghostcell_dims_conservative_variables[1];
        
        /*
         * Compute the internal energy field.
         */
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_conservative_variables) +
                        (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables +
                        (k + num_ghosts_2_conservative_variables)*ghostcell_dim_0_conservative_variables*
                            ghostcell_dim_1_conservative_variables;
                    
                    epsilon_x_L[idx] = (Q_x_L[4][idx] -
                        Real(1)/Real(2)*(Q_x_L[1][idx]*Q_x_L[1][idx] + Q_x_L[2][idx]*Q_x_L[2][idx] +
                        Q_x_L[3][idx]*Q_x_L[3][idx])/Q_x_L[0][idx])/Q_x_L[0][idx];
                }
            }
        }
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_conservative_variables) +
                        (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables +
                        (k + num_ghosts_2_conservative_variables)*ghostcell_dim_0_conservative_variables*
                            ghostcell_dim_1_conservative_variables;
                    
                    epsilon_x_R[idx] = (Q_x_R[4][idx] -
                        Real(1)/Real(2)*(Q_x_R[1][idx]*Q_x_R[1][idx] + Q_x_R[2][idx]*Q_x_R[2][idx] +
                        Q_x_R[3][idx]*Q_x_R[3][idx])/Q_x_R[0][idx])/Q_x_R[0][idx];
                }
            }
        }
        
        flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computePressure(
                pressure_x_L,
                conservative_variables_L[0],
                internal_energy_x_L,
                thermo_properties_const_ptr,
                0,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computePressure(
                pressure_x_R,
                conservative_variables_R[0],
                internal_energy_x_R,
                thermo_properties_const_ptr,
                0,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computeSoundSpeed(
                sound_speed_x_L,
                conservative_variables_L[0],
                pressure_x_L,
                thermo_properties_const_ptr,
                0,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computeSoundSpeed(
                sound_speed_x_R,
                conservative_variables_R[0],
                pressure_x_R,
                thermo_properties_const_ptr,
                0,
                domain);
        
        if (compute_velocity)
        {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(velocity);
#endif
            
            /*
             * Get the numbers of cells in each dimension, number of ghost cells of velocity and
             * the pointer to velocity.
             */
            
            const hier::IntVector num_ghosts_velocity = velocity->getGhostCellWidth();
            const hier::IntVector ghostcell_dims_velocity = velocity->getGhostBox().numberCells();
            
            const int num_ghosts_0_velocity = num_ghosts_velocity[0];
            const int num_ghosts_1_velocity = num_ghosts_velocity[1];
            const int num_ghosts_2_velocity = num_ghosts_velocity[2];
            const int ghostcell_dim_0_velocity = ghostcell_dims_velocity[0] + 1;
            const int ghostcell_dim_1_velocity = ghostcell_dims_velocity[1];
            
            Real* u = velocity->getPointer(0, 0);
            
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                    {
                        // Compute the linear indices.
                        const int idx_flux = (i + num_ghosts_0_convective_flux) +
                            (j + num_ghosts_1_convective_flux)*ghostcell_dim_0_convective_flux +
                            (k + num_ghosts_2_convective_flux)*ghostcell_dim_0_convective_flux*
                                ghostcell_dim_1_convective_flux;
                        
                        const int idx_velocity = (i + num_ghosts_0_velocity) +
                            (j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                            (k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                ghostcell_dim_1_velocity;
                        
                        const int idx = (i + num_ghosts_0_conservative_variables) +
                            (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables +
                            (k + num_ghosts_2_conservative_variables)*ghostcell_dim_0_conservative_variables*
                                ghostcell_dim_1_conservative_variables;
                        
                        computeLocalConvectiveFluxInXDirectionFromConservativeVariablesHLLC_HLL3D(
                            F_x.data(),
                            Q_x_L.data(),
                            Q_x_R.data(),
                            p_x_L,
                            p_x_R,
                            c_x_L,
                            c_x_R,
                            u_x_L,
                            u_x_R,
                            s_x_minus,
                            s_x_plus,
                            s_x_star,
                            Chi_x_star_LR,
                            idx_flux,
                            idx);
                        
                        if (s_x_star > Real(0))
                        {
                            u[idx_velocity] = u_x_L + s_x_minus*(Chi_x_star_LR - Real(1));
                        }
                        else
                        {
                            u[idx_velocity] = u_x_R + s_x_plus*(Chi_x_star_LR - Real(1));
                        }
                    }
                }
            }
        }
        else
        {
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                    {
                        // Compute the linear indices.
                        const int idx_flux = (i + num_ghosts_0_convective_flux) +
                            (j + num_ghosts_1_convective_flux)*ghostcell_dim_0_convective_flux +
                            (k + num_ghosts_2_convective_flux)*ghostcell_dim_0_convective_flux*
                                ghostcell_dim_1_convective_flux;
                        
                        const int idx = (i + num_ghosts_0_conservative_variables) +
                            (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables +
                            (k + num_ghosts_2_conservative_variables)*ghostcell_dim_0_conservative_variables*
                                ghostcell_dim_1_conservative_variables;
                        
                        computeLocalConvectiveFluxInXDirectionFromConservativeVariablesHLLC_HLL3D(
                            F_x.data(),
                            Q_x_L.data(),
                            Q_x_R.data(),
                            p_x_L,
                            p_x_R,
                            c_x_L,
                            c_x_R,
                            u_x_L,
                            u_x_R,
                            s_x_minus,
                            s_x_plus,
                            s_x_star,
                            Chi_x_star_LR,
                            idx_flux,
                            idx);
                    }
                }
            }
        }
    }
}


/*
 * Compute the convective flux and velocity in the y-direction from conservative variables with
 * HLLC-HLL Riemann solver.
 */
void
FlowModelRiemannSolverSingleSpecies::computeConvectiveFluxAndVelocityInYDirectionFromConservativeVariablesHLLC_HLL(
    HAMERS_SHARED_PTR<pdat::SideData<Real> > convective_flux,
    HAMERS_SHARED_PTR<pdat::SideData<Real> > velocity,
    const std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& conservative_variables_B,
    const std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& conservative_variables_T,
    const hier::Box& domain,
    bool compute_velocity) const
{
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    const int num_eqn = flow_model_tmp->getNumberOfEquations();
    
    // Get the box that covers the interior of patch.
    const hier::Box interior_box = convective_flux->getBox();
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_convective_flux = convective_flux->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_convective_flux =
        convective_flux->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_conservative_variables =
        conservative_variables_B[0]->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_conservative_variables =
        conservative_variables_B[0]->getGhostBox().numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        const hier::IntVector num_ghosts_min =
            hier::IntVector::min(num_ghosts_convective_flux, num_ghosts_conservative_variables);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(convective_flux->getGhostBox().contains(domain));
        TBOX_ASSERT(conservative_variables_B[0]->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
    /*
     * Get the equation of state mixing rules and the thermodynamic properties of the species.
     */
    
    const HAMERS_SHARED_PTR<EquationOfStateMixingRules> equation_of_state_mixing_rules =
        flow_model_tmp->getEquationOfStateMixingRules();
    
    const int num_thermo_properties = equation_of_state_mixing_rules->
        getNumberOfSpeciesThermodynamicProperties();
    
    std::vector<Real> thermo_properties;
    std::vector<Real*> thermo_properties_ptr;
    std::vector<const Real*> thermo_properties_const_ptr;
    
    thermo_properties.resize(num_thermo_properties);
    thermo_properties_ptr.reserve(num_thermo_properties);
    thermo_properties_const_ptr.reserve(num_thermo_properties);
    
    for (int ti = 0; ti < num_thermo_properties; ti++)
    {
        thermo_properties_ptr.push_back(&thermo_properties[ti]);
        thermo_properties_const_ptr.push_back(&thermo_properties[ti]);
    }
    
    equation_of_state_mixing_rules->getSpeciesThermodynamicProperties(
        thermo_properties_ptr,
        0);
    
    /*
     * Get the pointers to the side data of convective flux and conservative variables.
     */
    
    std::vector<Real*> F_y;
    F_y.reserve(num_eqn);
    for (int ei = 0; ei < num_eqn; ei++)
    {
        F_y.push_back(convective_flux->getPointer(1, ei));
    }
    
    std::vector<Real*> Q_y_B;
    std::vector<Real*> Q_y_T;
    Q_y_B.reserve(num_eqn);
    Q_y_T.reserve(num_eqn);
    for (int ei = 0; ei < num_eqn; ei++)
    {
        Q_y_B.push_back(conservative_variables_B[ei]->getPointer(1, 0));
        Q_y_T.push_back(conservative_variables_T[ei]->getPointer(1, 0));
    }
    
    /*
     * Allocate temporary data.
     */
    
    hier::IntVector direction_y = hier::IntVector::getZero(d_dim);
    direction_y[1] = 1;
    
    HAMERS_SHARED_PTR<pdat::SideData<Real> > internal_energy_y_B(
        new pdat::SideData<Real>(interior_box, 1, num_ghosts_conservative_variables,
            direction_y));

    HAMERS_SHARED_PTR<pdat::SideData<Real> > internal_energy_y_T(
        new pdat::SideData<Real>(interior_box, 1, num_ghosts_conservative_variables,
            direction_y));
    
    HAMERS_SHARED_PTR<pdat::SideData<Real> > pressure_y_B(
        new pdat::SideData<Real>(interior_box, 1, num_ghosts_conservative_variables,
            direction_y));
    
    HAMERS_SHARED_PTR<pdat::SideData<Real> > pressure_y_T(
        new pdat::SideData<Real>(interior_box, 1, num_ghosts_conservative_variables,
            direction_y));
    
    HAMERS_SHARED_PTR<pdat::SideData<Real> > sound_speed_y_B(
        new pdat::SideData<Real>(interior_box, 1, num_ghosts_conservative_variables,
            direction_y));
    
    HAMERS_SHARED_PTR<pdat::SideData<Real> > sound_speed_y_T(
        new pdat::SideData<Real>(interior_box, 1, num_ghosts_conservative_variables,
            direction_y));
    
    /*
     * Get the pointers to the temporary data.
     */
    
    Real* epsilon_y_B = internal_energy_y_B->getPointer(1, 0);
    Real* epsilon_y_T = internal_energy_y_T->getPointer(1, 0);
    
    Real* p_y_B = pressure_y_B->getPointer(1, 0);
    Real* p_y_T = pressure_y_T->getPointer(1, 0);
    
    Real* c_y_B = sound_speed_y_B->getPointer(1, 0);
    Real* c_y_T = sound_speed_y_T->getPointer(1, 0);
    
    Real v_y_B = Real(0);
    Real v_y_T = Real(0);
    
    Real s_y_minus = Real(0);
    Real s_y_plus  = Real(0);
    Real s_y_star  = Real(0);
    
    Real Chi_y_star_BT = Real(0);
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelRiemannSolverSingleSpecies::"
            << "computeConvectiveFluxAndVelocityInYDirectionFromConservativeVariablesHLLC_HLL()\n"
            << "There is no y direction for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the local lower indices, numbers of cells in each dimension and numbers of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        
        const int num_ghosts_0_convective_flux = num_ghosts_convective_flux[0];
        const int num_ghosts_1_convective_flux = num_ghosts_convective_flux[1];
        const int ghostcell_dim_0_convective_flux = ghostcell_dims_convective_flux[0];
        
        const int num_ghosts_0_conservative_variables = num_ghosts_conservative_variables[0];
        const int num_ghosts_1_conservative_variables = num_ghosts_conservative_variables[1];
        const int ghostcell_dim_0_conservative_variables = ghostcell_dims_conservative_variables[0];
        
        /*
         * Compute the internal energy field.
         */
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear index.
                const int idx = (i + num_ghosts_0_conservative_variables) +
                    (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables;
                
                epsilon_y_B[idx] = (Q_y_B[3][idx] -
                    Real(1)/Real(2)*(Q_y_B[1][idx]*Q_y_B[1][idx] + Q_y_B[2][idx]*Q_y_B[2][idx])/
                    Q_y_B[0][idx])/Q_y_B[0][idx];
            }
        }
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear index.
                const int idx = (i + num_ghosts_0_conservative_variables) +
                    (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables;
                
                epsilon_y_T[idx] = (Q_y_T[3][idx] -
                    Real(1)/Real(2)*(Q_y_T[1][idx]*Q_y_T[1][idx] + Q_y_T[2][idx]*Q_y_T[2][idx])/
                    Q_y_T[0][idx])/Q_y_T[0][idx];
            }
        }
        
        flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computePressure(
                pressure_y_B,
                conservative_variables_B[0],
                internal_energy_y_B,
                thermo_properties_const_ptr,
                1,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computePressure(
                pressure_y_T,
                conservative_variables_T[0],
                internal_energy_y_T,
                thermo_properties_const_ptr,
                1,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computeSoundSpeed(
                sound_speed_y_B,
                conservative_variables_B[0],
                pressure_y_B,
                thermo_properties_const_ptr,
                1,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computeSoundSpeed(
                sound_speed_y_T,
                conservative_variables_T[0],
                pressure_y_T,
                thermo_properties_const_ptr,
                1,
                domain);

        if (compute_velocity)
        {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(velocity);
#endif
            
            /*
             * Get the numbers of cells in each dimension, number of ghost cells of velocity and
             * the pointer to velocity.
             */
            
            const hier::IntVector num_ghosts_velocity = velocity->getGhostCellWidth();
            const hier::IntVector ghostcell_dims_velocity = velocity->getGhostBox().numberCells();
            
            const int num_ghosts_0_velocity = num_ghosts_velocity[0];
            const int num_ghosts_1_velocity = num_ghosts_velocity[1];
            const int ghostcell_dim_0_velocity = ghostcell_dims_velocity[0];
            
            Real* v = velocity->getPointer(1, 1);
            
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_flux = (i + num_ghosts_0_convective_flux) +
                        (j + num_ghosts_1_convective_flux)*ghostcell_dim_0_convective_flux;
                    
                    const int idx_velocity = (i + num_ghosts_0_velocity) +
                        (j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                    
                    const int idx = (i + num_ghosts_0_conservative_variables) +
                        (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables;
                    
                    computeLocalConvectiveFluxInYDirectionFromConservativeVariablesHLLC_HLL2D(
                        F_y.data(),
                        Q_y_B.data(),
                        Q_y_T.data(),
                        p_y_B,
                        p_y_T,
                        c_y_B,
                        c_y_T,
                        v_y_B,
                        v_y_T,
                        s_y_minus,
                        s_y_plus,
                        s_y_star,
                        Chi_y_star_BT,
                        idx_flux,
                        idx);
                    
                    if (s_y_star > Real(0))
                    {
                        v[idx_velocity] = v_y_B + s_y_minus*(Chi_y_star_BT - Real(1));
                    }
                    else
                    {
                        v[idx_velocity] = v_y_T + s_y_plus*(Chi_y_star_BT - Real(1));
                    }
                }
            }
        }
        else
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_flux = (i + num_ghosts_0_convective_flux) +
                        (j + num_ghosts_1_convective_flux)*ghostcell_dim_0_convective_flux;
                    
                    const int idx = (i + num_ghosts_0_conservative_variables) +
                        (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables;
                    
                    computeLocalConvectiveFluxInYDirectionFromConservativeVariablesHLLC_HLL2D(
                        F_y.data(),
                        Q_y_B.data(),
                        Q_y_T.data(),
                        p_y_B,
                        p_y_T,
                        c_y_B,
                        c_y_T,
                        v_y_B,
                        v_y_T,
                        s_y_minus,
                        s_y_plus,
                        s_y_star,
                        Chi_y_star_BT,
                        idx_flux,
                        idx);
                }
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the local lower indices, numbers of cells in each dimension and numbers of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_lo_2 = domain_lo[2];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        const int domain_dim_2 = domain_dims[2];
        
        const int num_ghosts_0_convective_flux = num_ghosts_convective_flux[0];
        const int num_ghosts_1_convective_flux = num_ghosts_convective_flux[1];
        const int num_ghosts_2_convective_flux = num_ghosts_convective_flux[2];
        const int ghostcell_dim_0_convective_flux = ghostcell_dims_convective_flux[0];
        const int ghostcell_dim_1_convective_flux = ghostcell_dims_convective_flux[1] + 1;
        
        const int num_ghosts_0_conservative_variables = num_ghosts_conservative_variables[0];
        const int num_ghosts_1_conservative_variables = num_ghosts_conservative_variables[1];
        const int num_ghosts_2_conservative_variables = num_ghosts_conservative_variables[2];
        const int ghostcell_dim_0_conservative_variables = ghostcell_dims_conservative_variables[0];
        const int ghostcell_dim_1_conservative_variables = ghostcell_dims_conservative_variables[1] + 1;
        
        /*
         * Compute the internal energy field.
         */
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_conservative_variables) +
                        (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables +
                        (k + num_ghosts_2_conservative_variables)*ghostcell_dim_0_conservative_variables*
                            ghostcell_dim_1_conservative_variables;
                    
                    epsilon_y_B[idx] = (Q_y_B[4][idx] -
                        Real(1)/Real(2)*(Q_y_B[1][idx]*Q_y_B[1][idx] + Q_y_B[2][idx]*Q_y_B[2][idx] +
                        Q_y_B[3][idx]*Q_y_B[3][idx])/Q_y_B[0][idx])/Q_y_B[0][idx];
                }
            }
        }
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_conservative_variables) +
                        (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables +
                        (k + num_ghosts_2_conservative_variables)*ghostcell_dim_0_conservative_variables*
                            ghostcell_dim_1_conservative_variables;
                    
                    epsilon_y_T[idx] = (Q_y_T[4][idx] -
                        Real(1)/Real(2)*(Q_y_T[1][idx]*Q_y_T[1][idx] + Q_y_T[2][idx]*Q_y_T[2][idx] +
                        Q_y_T[3][idx]*Q_y_T[3][idx])/Q_y_T[0][idx])/Q_y_T[0][idx];
                }
            }
        }
        
        flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computePressure(
                pressure_y_B,
                conservative_variables_B[0],
                internal_energy_y_B,
                thermo_properties_const_ptr,
                1,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computePressure(
                pressure_y_T,
                conservative_variables_T[0],
                internal_energy_y_T,
                thermo_properties_const_ptr,
                1,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computeSoundSpeed(
                sound_speed_y_B,
                conservative_variables_B[0],
                pressure_y_B,
                thermo_properties_const_ptr,
                1,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computeSoundSpeed(
                sound_speed_y_T,
                conservative_variables_T[0],
                pressure_y_T,
                thermo_properties_const_ptr,
                1,
                domain);
        
        if (compute_velocity)
        {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(velocity);
#endif
            
            /*
             * Get the numbers of cells in each dimension, number of ghost cells of velocity and
             * the pointer to velocity.
             */
            
            const hier::IntVector num_ghosts_velocity = velocity->getGhostCellWidth();
            const hier::IntVector ghostcell_dims_velocity = velocity->getGhostBox().numberCells();
            
            const int num_ghosts_0_velocity = num_ghosts_velocity[0];
            const int num_ghosts_1_velocity = num_ghosts_velocity[1];
            const int num_ghosts_2_velocity = num_ghosts_velocity[2];
            const int ghostcell_dim_0_velocity = ghostcell_dims_velocity[0];
            const int ghostcell_dim_1_velocity = ghostcell_dims_velocity[1] + 1;
            
            Real* v = velocity->getPointer(1, 1);
            
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_flux = (i + num_ghosts_0_convective_flux) +
                            (j + num_ghosts_1_convective_flux)*ghostcell_dim_0_convective_flux +
                            (k + num_ghosts_2_convective_flux)*ghostcell_dim_0_convective_flux*
                                ghostcell_dim_1_convective_flux;
                        
                        const int idx_velocity = (i + num_ghosts_0_velocity) +
                            (j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                            (k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                ghostcell_dim_1_velocity;
                        
                        const int idx = (i + num_ghosts_0_conservative_variables) +
                            (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables +
                            (k + num_ghosts_2_conservative_variables)*ghostcell_dim_0_conservative_variables*
                                ghostcell_dim_1_conservative_variables;
                        
                        computeLocalConvectiveFluxInYDirectionFromConservativeVariablesHLLC_HLL3D(
                            F_y.data(),
                            Q_y_B.data(),
                            Q_y_T.data(),
                            p_y_B,
                            p_y_T,
                            c_y_B,
                            c_y_T,
                            v_y_B,
                            v_y_T,
                            s_y_minus,
                            s_y_plus,
                            s_y_star,
                            Chi_y_star_BT,
                            idx_flux,
                            idx);
                        
                        if (s_y_star > Real(0))
                        {
                            v[idx_velocity] = v_y_B + s_y_minus*(Chi_y_star_BT - Real(1));
                        }
                        else
                        {
                            v[idx_velocity] = v_y_T + s_y_plus*(Chi_y_star_BT - Real(1));
                        }
                    }
                }
            }
        }
        else
        {
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_flux = (i + num_ghosts_0_convective_flux) +
                            (j + num_ghosts_1_convective_flux)*ghostcell_dim_0_convective_flux +
                            (k + num_ghosts_2_convective_flux)*ghostcell_dim_0_convective_flux*
                                ghostcell_dim_1_convective_flux;
                        
                        const int idx = (i + num_ghosts_0_conservative_variables) +
                            (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables +
                            (k + num_ghosts_2_conservative_variables)*ghostcell_dim_0_conservative_variables*
                                ghostcell_dim_1_conservative_variables;
                        
                        computeLocalConvectiveFluxInYDirectionFromConservativeVariablesHLLC_HLL3D(
                            F_y.data(),
                            Q_y_B.data(),
                            Q_y_T.data(),
                            p_y_B,
                            p_y_T,
                            c_y_B,
                            c_y_T,
                            v_y_B,
                            v_y_T,
                            s_y_minus,
                            s_y_plus,
                            s_y_star,
                            Chi_y_star_BT,
                            idx_flux,
                            idx);
                    }
                }
            }
        }
    }
}


/*
 * Compute the convective flux and velocity in the z-direction from conservative variables with
 * HLLC-HLL Riemann solver.
 */
void
FlowModelRiemannSolverSingleSpecies::computeConvectiveFluxAndVelocityInZDirectionFromConservativeVariablesHLLC_HLL(
    HAMERS_SHARED_PTR<pdat::SideData<Real> > convective_flux,
    HAMERS_SHARED_PTR<pdat::SideData<Real> > velocity,
    const std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& conservative_variables_B,
    const std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& conservative_variables_F,
    const hier::Box& domain,
    bool compute_velocity) const
{
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    const int num_eqn = flow_model_tmp->getNumberOfEquations();
    
    // Get the box that covers the interior of patch.
    const hier::Box interior_box = convective_flux->getBox();
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_convective_flux = convective_flux->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_convective_flux =
        convective_flux->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_conservative_variables =
        conservative_variables_B[0]->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_conservative_variables =
        conservative_variables_B[0]->getGhostBox().numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        const hier::IntVector num_ghosts_min =
            hier::IntVector::min(num_ghosts_convective_flux, num_ghosts_conservative_variables);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(convective_flux->getGhostBox().contains(domain));
        TBOX_ASSERT(conservative_variables_B[0]->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
    /*
     * Get the equation of state mixing rules and the thermodynamic properties of the species.
     */
    
    const HAMERS_SHARED_PTR<EquationOfStateMixingRules> equation_of_state_mixing_rules =
        flow_model_tmp->getEquationOfStateMixingRules();
    
    const int num_thermo_properties = equation_of_state_mixing_rules->
        getNumberOfSpeciesThermodynamicProperties();
    
    std::vector<Real> thermo_properties;
    std::vector<Real*> thermo_properties_ptr;
    std::vector<const Real*> thermo_properties_const_ptr;
    
    thermo_properties.resize(num_thermo_properties);
    thermo_properties_ptr.reserve(num_thermo_properties);
    thermo_properties_const_ptr.reserve(num_thermo_properties);
    
    for (int ti = 0; ti < num_thermo_properties; ti++)
    {
        thermo_properties_ptr.push_back(&thermo_properties[ti]);
        thermo_properties_const_ptr.push_back(&thermo_properties[ti]);
    }
    
    equation_of_state_mixing_rules->getSpeciesThermodynamicProperties(
        thermo_properties_ptr,
        0);
    
    /*
     * Get the pointers to the side data of convective flux and conservative variables.
     */
    
    std::vector<Real*> F_z;
    F_z.reserve(num_eqn);
    for (int ei = 0; ei < num_eqn; ei++)
    {
        F_z.push_back(convective_flux->getPointer(2, ei));
    }
    
    std::vector<Real*> Q_z_B;
    std::vector<Real*> Q_z_F;
    Q_z_B.reserve(num_eqn);
    Q_z_F.reserve(num_eqn);
    for (int ei = 0; ei < num_eqn; ei++)
    {
        Q_z_B.push_back(conservative_variables_B[ei]->getPointer(2, 0));
        Q_z_F.push_back(conservative_variables_F[ei]->getPointer(2, 0));
    }
    
    /*
     * Allocate temporary data.
     */
    
    hier::IntVector direction_z = hier::IntVector::getZero(d_dim);
    direction_z[2] = 1;
    
    HAMERS_SHARED_PTR<pdat::SideData<Real> > internal_energy_z_B(
        new pdat::SideData<Real>(interior_box, 1, num_ghosts_conservative_variables,
            direction_z));

    HAMERS_SHARED_PTR<pdat::SideData<Real> > internal_energy_z_F(
        new pdat::SideData<Real>(interior_box, 1, num_ghosts_conservative_variables,
            direction_z));
    
    HAMERS_SHARED_PTR<pdat::SideData<Real> > pressure_z_B(
        new pdat::SideData<Real>(interior_box, 1, num_ghosts_conservative_variables,
            direction_z));
    
    HAMERS_SHARED_PTR<pdat::SideData<Real> > pressure_z_F(
        new pdat::SideData<Real>(interior_box, 1, num_ghosts_conservative_variables,
            direction_z));
    
    HAMERS_SHARED_PTR<pdat::SideData<Real> > sound_speed_z_B(
        new pdat::SideData<Real>(interior_box, 1, num_ghosts_conservative_variables,
            direction_z));
    
    HAMERS_SHARED_PTR<pdat::SideData<Real> > sound_speed_z_F(
        new pdat::SideData<Real>(interior_box, 1, num_ghosts_conservative_variables,
            direction_z));
    
    /*
     * Get the pointers to the temporary data.
     */
    
    Real* epsilon_z_B = internal_energy_z_B->getPointer(2, 0);
    Real* epsilon_z_F = internal_energy_z_F->getPointer(2, 0);
    
    Real* p_z_B = pressure_z_B->getPointer(2, 0);
    Real* p_z_F = pressure_z_F->getPointer(2, 0);
    
    Real* c_z_B = sound_speed_z_B->getPointer(2, 0);
    Real* c_z_F = sound_speed_z_F->getPointer(2, 0);
    
    Real w_z_B = Real(0);
    Real w_z_F = Real(0);
    
    Real s_z_minus = Real(0);
    Real s_z_plus  = Real(0);
    Real s_z_star  = Real(0);
    
    Real Chi_z_star_BF = Real(0);
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelRiemannSolverSingleSpecies::"
            << "computeConvectiveFluxAndVelocityInZDirectionFromConservativeVariablesHLLC_HLL()\n"
            << "There is no z direction for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelRiemannSolverSingleSpecies::"
            << "computeConvectiveFluxAndVelocityInZDirectionFromConservativeVariablesHLLC_HLL()\n"
            << "There is no z direction for two-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the local lower indices, numbers of cells in each dimension and numbers of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_lo_2 = domain_lo[2];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        const int domain_dim_2 = domain_dims[2];
        
        const int num_ghosts_0_convective_flux = num_ghosts_convective_flux[0];
        const int num_ghosts_1_convective_flux = num_ghosts_convective_flux[1];
        const int num_ghosts_2_convective_flux = num_ghosts_convective_flux[2];
        const int ghostcell_dim_0_convective_flux = ghostcell_dims_convective_flux[0];
        const int ghostcell_dim_1_convective_flux = ghostcell_dims_convective_flux[1];
        
        const int num_ghosts_0_conservative_variables = num_ghosts_conservative_variables[0];
        const int num_ghosts_1_conservative_variables = num_ghosts_conservative_variables[1];
        const int num_ghosts_2_conservative_variables = num_ghosts_conservative_variables[2];
        const int ghostcell_dim_0_conservative_variables = ghostcell_dims_conservative_variables[0];
        const int ghostcell_dim_1_conservative_variables = ghostcell_dims_conservative_variables[1];
        
        /*
         * Compute the internal energy field.
         */
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2 + 1; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_conservative_variables) +
                        (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables +
                        (k + num_ghosts_2_conservative_variables)*ghostcell_dim_0_conservative_variables*
                            ghostcell_dim_1_conservative_variables;
                    
                    epsilon_z_B[idx] = (Q_z_B[4][idx] -
                        Real(1)/Real(2)*(Q_z_B[1][idx]*Q_z_B[1][idx] + Q_z_B[2][idx]*Q_z_B[2][idx] +
                        Q_z_B[3][idx]*Q_z_B[3][idx])/Q_z_B[0][idx])/Q_z_B[0][idx];
                }
            }
        }
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2 + 1; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_conservative_variables) +
                        (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables +
                        (k + num_ghosts_2_conservative_variables)*ghostcell_dim_0_conservative_variables*
                            ghostcell_dim_1_conservative_variables;
                    
                    epsilon_z_F[idx] = (Q_z_F[4][idx] -
                        Real(1)/Real(2)*(Q_z_F[1][idx]*Q_z_F[1][idx] + Q_z_F[2][idx]*Q_z_F[2][idx] +
                        Q_z_F[3][idx]*Q_z_F[3][idx])/Q_z_F[0][idx])/Q_z_F[0][idx];
                }
            }
        }
        
        flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computePressure(
                pressure_z_B,
                conservative_variables_B[0],
                internal_energy_z_B,
                thermo_properties_const_ptr,
                2,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computePressure(
                pressure_z_F,
                conservative_variables_F[0],
                internal_energy_z_F,
                thermo_properties_const_ptr,
                2,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computeSoundSpeed(
                sound_speed_z_B,
                conservative_variables_B[0],
                pressure_z_B,
                thermo_properties_const_ptr,
                2,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computeSoundSpeed(
                sound_speed_z_F,
                conservative_variables_F[0],
                pressure_z_F,
                thermo_properties_const_ptr,
                2,
                domain);
        
        if (compute_velocity)
        {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(velocity);
#endif
            
            /*
             * Get the numbers of cells in each dimension, number of ghost cells of velocity and
             * the pointer to velocity.
             */
            
            const hier::IntVector num_ghosts_velocity = velocity->getGhostCellWidth();
            const hier::IntVector ghostcell_dims_velocity = velocity->getGhostBox().numberCells();
            
            const int num_ghosts_0_velocity = num_ghosts_velocity[0];
            const int num_ghosts_1_velocity = num_ghosts_velocity[1];
            const int num_ghosts_2_velocity = num_ghosts_velocity[2];
            const int ghostcell_dim_0_velocity = ghostcell_dims_velocity[0];
            const int ghostcell_dim_1_velocity = ghostcell_dims_velocity[1];
            
            Real* w = velocity->getPointer(2, 2);
            
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2 + 1; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_flux = (i + num_ghosts_0_convective_flux) +
                            (j + num_ghosts_1_convective_flux)*ghostcell_dim_0_convective_flux +
                            (k + num_ghosts_2_convective_flux)*ghostcell_dim_0_convective_flux*
                                ghostcell_dim_1_convective_flux;
                        
                        const int idx_velocity = (i + num_ghosts_0_velocity) +
                            (j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                            (k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                ghostcell_dim_1_velocity;
                        
                        const int idx = (i + num_ghosts_0_conservative_variables) +
                            (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables +
                            (k + num_ghosts_2_conservative_variables)*ghostcell_dim_0_conservative_variables*
                                ghostcell_dim_1_conservative_variables;
                        
                        computeLocalConvectiveFluxInZDirectionFromConservativeVariablesHLLC_HLL3D(
                            F_z.data(),
                            Q_z_B.data(),
                            Q_z_F.data(),
                            p_z_B,
                            p_z_F,
                            c_z_B,
                            c_z_F,
                            w_z_B,
                            w_z_F,
                            s_z_minus,
                            s_z_plus,
                            s_z_star,
                            Chi_z_star_BF,
                            idx_flux,
                            idx);
                        
                        if (s_z_star > Real(0))
                        {
                            w[idx_velocity] = w_z_B + s_z_minus*(Chi_z_star_BF - Real(1));
                        }
                        else
                        {
                            w[idx_velocity] = w_z_F + s_z_plus*(Chi_z_star_BF - Real(1));
                        }
                    }
                }
            }
        }
        else
        {
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2 + 1; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_flux = (i + num_ghosts_0_convective_flux) +
                            (j + num_ghosts_1_convective_flux)*ghostcell_dim_0_convective_flux +
                            (k + num_ghosts_2_convective_flux)*ghostcell_dim_0_convective_flux*
                                ghostcell_dim_1_convective_flux;
                        
                        const int idx = (i + num_ghosts_0_conservative_variables) +
                            (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables +
                            (k + num_ghosts_2_conservative_variables)*ghostcell_dim_0_conservative_variables*
                                ghostcell_dim_1_conservative_variables;
                        
                        computeLocalConvectiveFluxInZDirectionFromConservativeVariablesHLLC_HLL3D(
                            F_z.data(),
                            Q_z_B.data(),
                            Q_z_F.data(),
                            p_z_B,
                            p_z_F,
                            c_z_B,
                            c_z_F,
                            w_z_B,
                            w_z_F,
                            s_z_minus,
                            s_z_plus,
                            s_z_star,
                            Chi_z_star_BF,
                            idx_flux,
                            idx);
                    }
                }
            }
        }
    }
}


/*
 * Compute the convective flux and velocity in the x-direction from primitive variables with
 * HLLC-HLL Riemann solver.
 */
void
FlowModelRiemannSolverSingleSpecies::computeConvectiveFluxAndVelocityInXDirectionFromPrimitiveVariablesHLLC_HLL(
    HAMERS_SHARED_PTR<pdat::SideData<Real> > convective_flux,
    HAMERS_SHARED_PTR<pdat::SideData<Real> > velocity,
    const std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& primitive_variables_L,
    const std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& primitive_variables_R,
    const hier::Box& domain,
    bool compute_velocity) const
{
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    const int num_eqn = flow_model_tmp->getNumberOfEquations();
    
    // Get the box that covers the interior of patch.
    const hier::Box interior_box = convective_flux->getBox();
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_convective_flux = convective_flux->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_convective_flux =
        convective_flux->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_primitive_variables =
        primitive_variables_L[0]->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_primitive_variables =
        primitive_variables_L[0]->getGhostBox().numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        const hier::IntVector num_ghosts_min =
            hier::IntVector::min(num_ghosts_convective_flux, num_ghosts_primitive_variables);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(convective_flux->getGhostBox().contains(domain));
        TBOX_ASSERT(primitive_variables_L[0]->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
    /*
     * Get the equation of state mixing rules and the thermodynamic properties of the species.
     */
    
    const HAMERS_SHARED_PTR<EquationOfStateMixingRules> equation_of_state_mixing_rules =
        flow_model_tmp->getEquationOfStateMixingRules();
    
    const int num_thermo_properties = equation_of_state_mixing_rules->
        getNumberOfSpeciesThermodynamicProperties();
    
    std::vector<Real> thermo_properties;
    std::vector<Real*> thermo_properties_ptr;
    std::vector<const Real*> thermo_properties_const_ptr;
    
    thermo_properties.resize(num_thermo_properties);
    thermo_properties_ptr.reserve(num_thermo_properties);
    thermo_properties_const_ptr.reserve(num_thermo_properties);
    
    for (int ti = 0; ti < num_thermo_properties; ti++)
    {
        thermo_properties_ptr.push_back(&thermo_properties[ti]);
        thermo_properties_const_ptr.push_back(&thermo_properties[ti]);
    }
    
    equation_of_state_mixing_rules->getSpeciesThermodynamicProperties(
        thermo_properties_ptr,
        0);
    
    /*
     * Get the pointers to the side data of convective flux and primitive variables.
     */
    
    std::vector<Real*> F_x;
    F_x.reserve(num_eqn);
    for (int ei = 0; ei < num_eqn; ei++)
    {
        F_x.push_back(convective_flux->getPointer(0, ei));
    }
    
    std::vector<Real*> V_x_L;
    std::vector<Real*> V_x_R;
    V_x_L.reserve(num_eqn);
    V_x_R.reserve(num_eqn);
    for (int ei = 0; ei < num_eqn; ei++)
    {
        V_x_L.push_back(primitive_variables_L[ei]->getPointer(0, 0));
        V_x_R.push_back(primitive_variables_R[ei]->getPointer(0, 0));
    }
    
    /*
     * Allocate temporary data.
     */
    
    hier::IntVector direction_x = hier::IntVector::getZero(d_dim);
    direction_x[0] = 1;
    
    HAMERS_SHARED_PTR<pdat::SideData<Real> > sound_speed_x_L(
        new pdat::SideData<Real>(interior_box, 1, num_ghosts_primitive_variables,
            direction_x));
    
    HAMERS_SHARED_PTR<pdat::SideData<Real> > sound_speed_x_R(
        new pdat::SideData<Real>(interior_box, 1, num_ghosts_primitive_variables,
            direction_x));
    
    HAMERS_SHARED_PTR<pdat::SideData<Real> > internal_energy_x_L(
        new pdat::SideData<Real>(interior_box, 1, num_ghosts_primitive_variables,
            direction_x));

    HAMERS_SHARED_PTR<pdat::SideData<Real> > internal_energy_x_R(
        new pdat::SideData<Real>(interior_box, 1, num_ghosts_primitive_variables,
            direction_x));
    
    /*
     * Get the pointers to the temporary data.
     */
    
    Real* c_x_L = sound_speed_x_L->getPointer(0, 0);
    Real* c_x_R = sound_speed_x_R->getPointer(0, 0);
    
    Real* epsilon_x_L = internal_energy_x_L->getPointer(0, 0);
    Real* epsilon_x_R = internal_energy_x_R->getPointer(0, 0);
    
    Real s_x_minus = Real(0);
    Real s_x_plus  = Real(0);
    Real s_x_star  = Real(0);
    
    Real Chi_x_star_LR = Real(0);
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and numbers of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int num_ghosts_0_convective_flux = num_ghosts_convective_flux[0];
        const int num_ghosts_0_primitive_variables = num_ghosts_primitive_variables[0];
        
        flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computeSoundSpeed(
                sound_speed_x_L,
                primitive_variables_L[0],
                primitive_variables_L[2],
                thermo_properties_const_ptr,
                0,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computeSoundSpeed(
                sound_speed_x_R,
                primitive_variables_R[0],
                primitive_variables_R[2],
                thermo_properties_const_ptr,
                0,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computeInternalEnergy(
                internal_energy_x_L,
                primitive_variables_L[0],
                primitive_variables_L[2],
                thermo_properties_const_ptr,
                0,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computeInternalEnergy(
                internal_energy_x_R,
                primitive_variables_R[0],
                primitive_variables_R[2],
                thermo_properties_const_ptr,
                0,
                domain);
        
        if (compute_velocity)
        {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(velocity);
#endif
            
            /*
             * Get the number of ghost cells of velocity and the pointer to velocity.
             */
            
            const hier::IntVector num_ghosts_velocity = velocity->getGhostCellWidth();
            
            const int num_ghosts_0_velocity = num_ghosts_velocity[0];
            
            Real* u = velocity->getPointer(0, 0);
            
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
            {
                // Compute the linear indices.
                const int idx_flux = i + num_ghosts_0_convective_flux;
                const int idx_velocity = i + num_ghosts_0_velocity;
                const int idx = i + num_ghosts_0_primitive_variables;
                
                computeLocalConvectiveFluxInXDirectionFromPrimitiveVariablesHLLC_HLL1D(
                    F_x.data(),
                    V_x_L.data(),
                    V_x_R.data(),
                    c_x_L,
                    c_x_R,
                    epsilon_x_L,
                    epsilon_x_R,
                    s_x_minus,
                    s_x_plus,
                    s_x_star,
                    Chi_x_star_LR,
                    idx_flux,
                    idx);
                
                if (s_x_star > Real(0))
                {
                    u[idx_velocity] = V_x_L[1][idx] + s_x_minus*(Chi_x_star_LR - Real(1));
                }
                else
                {
                    u[idx_velocity] = V_x_R[1][idx] + s_x_plus*(Chi_x_star_LR - Real(1));
                }
            }
        }
        else
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
            {
                // Compute the linear indices.
                const int idx_flux = i + num_ghosts_0_convective_flux;
                const int idx = i + num_ghosts_0_primitive_variables;
                
                computeLocalConvectiveFluxInXDirectionFromPrimitiveVariablesHLLC_HLL1D(
                    F_x.data(),
                    V_x_L.data(),
                    V_x_R.data(),
                    c_x_L,
                    c_x_R,
                    epsilon_x_L,
                    epsilon_x_R,
                    s_x_minus,
                    s_x_plus,
                    s_x_star,
                    Chi_x_star_LR,
                    idx_flux,
                    idx);
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the local lower indices, numbers of cells in each dimension and numbers of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        
        const int num_ghosts_0_convective_flux = num_ghosts_convective_flux[0];
        const int num_ghosts_1_convective_flux = num_ghosts_convective_flux[1];
        const int ghostcell_dim_0_convective_flux = ghostcell_dims_convective_flux[0] + 1;
        
        const int num_ghosts_0_primitive_variables = num_ghosts_primitive_variables[0];
        const int num_ghosts_1_primitive_variables = num_ghosts_primitive_variables[1];
        const int ghostcell_dim_0_primitive_variables = ghostcell_dims_primitive_variables[0] + 1;
        
        flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computeSoundSpeed(
                sound_speed_x_L,
                primitive_variables_L[0],
                primitive_variables_L[3],
                thermo_properties_const_ptr,
                0,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computeSoundSpeed(
                sound_speed_x_R,
                primitive_variables_R[0],
                primitive_variables_R[3],
                thermo_properties_const_ptr,
                0,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computeInternalEnergy(
                internal_energy_x_L,
                primitive_variables_L[0],
                primitive_variables_L[3],
                thermo_properties_const_ptr,
                0,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computeInternalEnergy(
                internal_energy_x_R,
                primitive_variables_R[0],
                primitive_variables_R[3],
                thermo_properties_const_ptr,
                0,
                domain);
        
        if (compute_velocity)
        {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(velocity);
#endif
            
            /*
             * Get the numbers of cells in each dimension, number of ghost cells of velocity and
             * the pointer to velocity.
             */
            
            const hier::IntVector num_ghosts_velocity = velocity->getGhostCellWidth();
            const hier::IntVector ghostcell_dims_velocity = velocity->getGhostBox().numberCells();
            
            const int num_ghosts_0_velocity = num_ghosts_velocity[0];
            const int num_ghosts_1_velocity = num_ghosts_velocity[1];
            const int ghostcell_dim_0_velocity = ghostcell_dims_velocity[0] + 1;
            
            Real* u = velocity->getPointer(0, 0);
            
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                {
                    // Compute the linear indices.
                    const int idx_flux = (i + num_ghosts_0_convective_flux) +
                        (j + num_ghosts_1_convective_flux)*ghostcell_dim_0_convective_flux;
                    
                    const int idx_velocity = (i + num_ghosts_0_velocity) +
                        (j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                    
                    const int idx = (i + num_ghosts_0_primitive_variables) +
                        (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables;
                    
                    computeLocalConvectiveFluxInXDirectionFromPrimitiveVariablesHLLC_HLL2D(
                        F_x.data(),
                        V_x_L.data(),
                        V_x_R.data(),
                        c_x_L,
                        c_x_R,
                        epsilon_x_L,
                        epsilon_x_R,
                        s_x_minus,
                        s_x_plus,
                        s_x_star,
                        Chi_x_star_LR,
                        idx_flux,
                        idx);
                    
                    if (s_x_star > Real(0))
                    {
                        u[idx_velocity] = V_x_L[1][idx] + s_x_minus*(Chi_x_star_LR - Real(1));
                    }
                    else
                    {
                        u[idx_velocity] = V_x_R[1][idx] + s_x_plus*(Chi_x_star_LR - Real(1));
                    }
                }
            }
        }
        else
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                {
                    // Compute the linear indices.
                    const int idx_flux = (i + num_ghosts_0_convective_flux) +
                        (j + num_ghosts_1_convective_flux)*ghostcell_dim_0_convective_flux;
                    
                    const int idx = (i + num_ghosts_0_primitive_variables) +
                        (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables;
                    
                    computeLocalConvectiveFluxInXDirectionFromPrimitiveVariablesHLLC_HLL2D(
                        F_x.data(),
                        V_x_L.data(),
                        V_x_R.data(),
                        c_x_L,
                        c_x_R,
                        epsilon_x_L,
                        epsilon_x_R,
                        s_x_minus,
                        s_x_plus,
                        s_x_star,
                        Chi_x_star_LR,
                        idx_flux,
                        idx);
                }
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the local lower indices, numbers of cells in each dimension and numbers of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_lo_2 = domain_lo[2];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        const int domain_dim_2 = domain_dims[2];
        
        const int num_ghosts_0_convective_flux = num_ghosts_convective_flux[0];
        const int num_ghosts_1_convective_flux = num_ghosts_convective_flux[1];
        const int num_ghosts_2_convective_flux = num_ghosts_convective_flux[2];
        const int ghostcell_dim_0_convective_flux = ghostcell_dims_convective_flux[0] + 1;
        const int ghostcell_dim_1_convective_flux = ghostcell_dims_convective_flux[1];
        
        const int num_ghosts_0_primitive_variables = num_ghosts_primitive_variables[0];
        const int num_ghosts_1_primitive_variables = num_ghosts_primitive_variables[1];
        const int num_ghosts_2_primitive_variables = num_ghosts_primitive_variables[2];
        const int ghostcell_dim_0_primitive_variables = ghostcell_dims_primitive_variables[0] + 1;
        const int ghostcell_dim_1_primitive_variables = ghostcell_dims_primitive_variables[1];
        
        flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computeSoundSpeed(
                sound_speed_x_L,
                primitive_variables_L[0],
                primitive_variables_L[4],
                thermo_properties_const_ptr,
                0,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computeSoundSpeed(
                sound_speed_x_R,
                primitive_variables_R[0],
                primitive_variables_R[4],
                thermo_properties_const_ptr,
                0,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computeInternalEnergy(
                internal_energy_x_L,
                primitive_variables_L[0],
                primitive_variables_L[4],
                thermo_properties_const_ptr,
                0,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computeInternalEnergy(
                internal_energy_x_R,
                primitive_variables_R[0],
                primitive_variables_R[4],
                thermo_properties_const_ptr,
                0,
                domain);
        
        if (compute_velocity)
        {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(velocity);
#endif
            
            /*
             * Get the numbers of cells in each dimension, number of ghost cells of velocity and
             * the pointer to velocity.
             */
            
            const hier::IntVector num_ghosts_velocity = velocity->getGhostCellWidth();
            const hier::IntVector ghostcell_dims_velocity = velocity->getGhostBox().numberCells();
            
            const int num_ghosts_0_velocity = num_ghosts_velocity[0];
            const int num_ghosts_1_velocity = num_ghosts_velocity[1];
            const int num_ghosts_2_velocity = num_ghosts_velocity[2];
            const int ghostcell_dim_0_velocity = ghostcell_dims_velocity[0] + 1;
            const int ghostcell_dim_1_velocity = ghostcell_dims_velocity[1];
            
            Real* u = velocity->getPointer(0, 0);
            
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                    {
                        // Compute the linear indices.
                        const int idx_flux = (i + num_ghosts_0_convective_flux) +
                            (j + num_ghosts_1_convective_flux)*ghostcell_dim_0_convective_flux +
                            (k + num_ghosts_2_convective_flux)*ghostcell_dim_0_convective_flux*
                                ghostcell_dim_1_convective_flux;
                        
                        const int idx_velocity = (i + num_ghosts_0_velocity) +
                            (j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                            (k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                ghostcell_dim_1_velocity;
                        
                        const int idx = (i + num_ghosts_0_primitive_variables) +
                            (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables +
                            (k + num_ghosts_2_primitive_variables)*ghostcell_dim_0_primitive_variables*
                                ghostcell_dim_1_primitive_variables;
                        
                        computeLocalConvectiveFluxInXDirectionFromPrimitiveVariablesHLLC_HLL3D(
                            F_x.data(),
                            V_x_L.data(),
                            V_x_R.data(),
                            c_x_L,
                            c_x_R,
                            epsilon_x_L,
                            epsilon_x_R,
                            s_x_minus,
                            s_x_plus,
                            s_x_star,
                            Chi_x_star_LR,
                            idx_flux,
                            idx);
                        
                        if (s_x_star > Real(0))
                        {
                            u[idx_velocity] = V_x_L[1][idx] + s_x_minus*(Chi_x_star_LR - Real(1));
                        }
                        else
                        {
                            u[idx_velocity] = V_x_R[1][idx] + s_x_plus*(Chi_x_star_LR - Real(1));
                        }
                    }
                }
            }
        }
        else
        {
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                    {
                        // Compute the linear indices.
                        const int idx_flux = (i + num_ghosts_0_convective_flux) +
                            (j + num_ghosts_1_convective_flux)*ghostcell_dim_0_convective_flux +
                            (k + num_ghosts_2_convective_flux)*ghostcell_dim_0_convective_flux*
                                ghostcell_dim_1_convective_flux;
                        
                        const int idx = (i + num_ghosts_0_primitive_variables) +
                            (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables +
                            (k + num_ghosts_2_primitive_variables)*ghostcell_dim_0_primitive_variables*
                                ghostcell_dim_1_primitive_variables;
                        
                        computeLocalConvectiveFluxInXDirectionFromPrimitiveVariablesHLLC_HLL3D(
                            F_x.data(),
                            V_x_L.data(),
                            V_x_R.data(),
                            c_x_L,
                            c_x_R,
                            epsilon_x_L,
                            epsilon_x_R,
                            s_x_minus,
                            s_x_plus,
                            s_x_star,
                            Chi_x_star_LR,
                            idx_flux,
                            idx);
                    }
                }
            }
        }
    }
}


/*
 * Compute the convective flux and velocity in the y-direction from primitive variables with
 * HLLC-HLL Riemann solver.
 */
void
FlowModelRiemannSolverSingleSpecies::computeConvectiveFluxAndVelocityInYDirectionFromPrimitiveVariablesHLLC_HLL(
    HAMERS_SHARED_PTR<pdat::SideData<Real> > convective_flux,
    HAMERS_SHARED_PTR<pdat::SideData<Real> > velocity,
    const std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& primitive_variables_B,
    const std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& primitive_variables_T,
    const hier::Box& domain,
    bool compute_velocity) const
{
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    const int num_eqn = flow_model_tmp->getNumberOfEquations();
    
    // Get the box that covers the interior of patch.
    const hier::Box interior_box = convective_flux->getBox();
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_convective_flux = convective_flux->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_convective_flux =
        convective_flux->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_primitive_variables =
        primitive_variables_B[0]->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_primitive_variables =
        primitive_variables_B[0]->getGhostBox().numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        const hier::IntVector num_ghosts_min =
            hier::IntVector::min(num_ghosts_convective_flux, num_ghosts_primitive_variables);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(convective_flux->getGhostBox().contains(domain));
        TBOX_ASSERT(primitive_variables_B[0]->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
    /*
     * Get the equation of state mixing rules and the thermodynamic properties of the species.
     */
    
    const HAMERS_SHARED_PTR<EquationOfStateMixingRules> equation_of_state_mixing_rules =
        flow_model_tmp->getEquationOfStateMixingRules();
    
    const int num_thermo_properties = equation_of_state_mixing_rules->
        getNumberOfSpeciesThermodynamicProperties();
    
    std::vector<Real> thermo_properties;
    std::vector<Real*> thermo_properties_ptr;
    std::vector<const Real*> thermo_properties_const_ptr;
    
    thermo_properties.resize(num_thermo_properties);
    thermo_properties_ptr.reserve(num_thermo_properties);
    thermo_properties_const_ptr.reserve(num_thermo_properties);
    
    for (int ti = 0; ti < num_thermo_properties; ti++)
    {
        thermo_properties_ptr.push_back(&thermo_properties[ti]);
        thermo_properties_const_ptr.push_back(&thermo_properties[ti]);
    }
    
    equation_of_state_mixing_rules->getSpeciesThermodynamicProperties(
        thermo_properties_ptr,
        0);
    
    /*
     * Get the pointers to the side data of convective flux and primitive variables.
     */
    
    std::vector<Real*> F_y;
    F_y.reserve(num_eqn);
    for (int ei = 0; ei < num_eqn; ei++)
    {
        F_y.push_back(convective_flux->getPointer(1, ei));
    }
    
    std::vector<Real*> V_y_B;
    std::vector<Real*> V_y_T;
    V_y_B.reserve(num_eqn);
    V_y_T.reserve(num_eqn);
    for (int ei = 0; ei < num_eqn; ei++)
    {
        V_y_B.push_back(primitive_variables_B[ei]->getPointer(1, 0));
        V_y_T.push_back(primitive_variables_T[ei]->getPointer(1, 0));
    }
    
    /*
     * Allocate temporary data.
     */
    
    hier::IntVector direction_y = hier::IntVector::getZero(d_dim);
    direction_y[1] = 1;
    
    HAMERS_SHARED_PTR<pdat::SideData<Real> > sound_speed_y_B(
        new pdat::SideData<Real>(interior_box, 1, num_ghosts_primitive_variables,
            direction_y));
    
    HAMERS_SHARED_PTR<pdat::SideData<Real> > sound_speed_y_T(
        new pdat::SideData<Real>(interior_box, 1, num_ghosts_primitive_variables,
            direction_y));
    
    HAMERS_SHARED_PTR<pdat::SideData<Real> > internal_energy_y_B(
        new pdat::SideData<Real>(interior_box, 1, num_ghosts_primitive_variables,
            direction_y));

    HAMERS_SHARED_PTR<pdat::SideData<Real> > internal_energy_y_T(
        new pdat::SideData<Real>(interior_box, 1, num_ghosts_primitive_variables,
            direction_y));
    
    /*
     * Get the pointers to the temporary data.
     */
    
    Real* c_y_B = sound_speed_y_B->getPointer(1, 0);
    Real* c_y_T = sound_speed_y_T->getPointer(1, 0);
    
    Real* epsilon_y_B = internal_energy_y_B->getPointer(1, 0);
    Real* epsilon_y_T = internal_energy_y_T->getPointer(1, 0);
    
    Real s_y_minus = Real(0);
    Real s_y_plus  = Real(0);
    Real s_y_star  = Real(0);
    
    Real Chi_y_star_BT = Real(0);
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelRiemannSolverSingleSpecies::"
            << "computeConvectiveFluxAndVelocityInYDirectionFromPrimitiveVariablesHLLC_HLL()\n"
            << "There is no y direction for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the local lower indices, numbers of cells in each dimension and numbers of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        
        const int num_ghosts_0_convective_flux = num_ghosts_convective_flux[0];
        const int num_ghosts_1_convective_flux = num_ghosts_convective_flux[1];
        const int ghostcell_dim_0_convective_flux = ghostcell_dims_convective_flux[0];
        
        const int num_ghosts_0_primitive_variables = num_ghosts_primitive_variables[0];
        const int num_ghosts_1_primitive_variables = num_ghosts_primitive_variables[1];
        const int ghostcell_dim_0_primitive_variables = ghostcell_dims_primitive_variables[0];
        
        flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computeSoundSpeed(
                sound_speed_y_B,
                primitive_variables_B[0],
                primitive_variables_B[3],
                thermo_properties_const_ptr,
                1,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computeSoundSpeed(
                sound_speed_y_T,
                primitive_variables_T[0],
                primitive_variables_T[3],
                thermo_properties_const_ptr,
                1,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computeInternalEnergy(
                internal_energy_y_B,
                primitive_variables_B[0],
                primitive_variables_B[3],
                thermo_properties_const_ptr,
                1,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computeInternalEnergy(
                internal_energy_y_T,
                primitive_variables_T[0],
                primitive_variables_T[3],
                thermo_properties_const_ptr,
                1,
                domain);
        
        if (compute_velocity)
        {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(velocity);
#endif
            
            /*
             * Get the numbers of cells in each dimension, number of ghost cells of velocity and
             * the pointer to velocity.
             */
            
            const hier::IntVector num_ghosts_velocity = velocity->getGhostCellWidth();
            const hier::IntVector ghostcell_dims_velocity = velocity->getGhostBox().numberCells();
            
            const int num_ghosts_0_velocity = num_ghosts_velocity[0];
            const int num_ghosts_1_velocity = num_ghosts_velocity[1];
            const int ghostcell_dim_0_velocity = ghostcell_dims_velocity[0];
            
            Real* v = velocity->getPointer(1, 1);
            
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_flux = (i + num_ghosts_0_convective_flux) +
                        (j + num_ghosts_1_convective_flux)*ghostcell_dim_0_convective_flux;
                    
                    const int idx_velocity = (i + num_ghosts_0_velocity) +
                        (j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                    
                    const int idx = (i + num_ghosts_0_primitive_variables) +
                        (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables;
                    
                    computeLocalConvectiveFluxInYDirectionFromPrimitiveVariablesHLLC_HLL2D(
                        F_y.data(),
                        V_y_B.data(),
                        V_y_T.data(),
                        c_y_B,
                        c_y_T,
                        epsilon_y_B,
                        epsilon_y_T,
                        s_y_minus,
                        s_y_plus,
                        s_y_star,
                        Chi_y_star_BT,
                        idx_flux,
                        idx);
                    
                    if (s_y_star > Real(0))
                    {
                        v[idx_velocity] = V_y_B[2][idx] + s_y_minus*(Chi_y_star_BT - Real(1));
                    }
                    else
                    {
                        v[idx_velocity] = V_y_T[2][idx] + s_y_plus*(Chi_y_star_BT - Real(1));
                    }
                }
            }
        }
        else
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_flux = (i + num_ghosts_0_convective_flux) +
                        (j + num_ghosts_1_convective_flux)*ghostcell_dim_0_convective_flux;
                    
                    const int idx = (i + num_ghosts_0_primitive_variables) +
                        (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables;
                    
                    computeLocalConvectiveFluxInYDirectionFromPrimitiveVariablesHLLC_HLL2D(
                        F_y.data(),
                        V_y_B.data(),
                        V_y_T.data(),
                        c_y_B,
                        c_y_T,
                        epsilon_y_B,
                        epsilon_y_T,
                        s_y_minus,
                        s_y_plus,
                        s_y_star,
                        Chi_y_star_BT,
                        idx_flux,
                        idx);
                }
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the local lower indices, numbers of cells in each dimension and numbers of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_lo_2 = domain_lo[2];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        const int domain_dim_2 = domain_dims[2];
        
        const int num_ghosts_0_convective_flux = num_ghosts_convective_flux[0];
        const int num_ghosts_1_convective_flux = num_ghosts_convective_flux[1];
        const int num_ghosts_2_convective_flux = num_ghosts_convective_flux[2];
        const int ghostcell_dim_0_convective_flux = ghostcell_dims_convective_flux[0];
        const int ghostcell_dim_1_convective_flux = ghostcell_dims_convective_flux[1] + 1;
        
        const int num_ghosts_0_primitive_variables = num_ghosts_primitive_variables[0];
        const int num_ghosts_1_primitive_variables = num_ghosts_primitive_variables[1];
        const int num_ghosts_2_primitive_variables = num_ghosts_primitive_variables[2];
        const int ghostcell_dim_0_primitive_variables = ghostcell_dims_primitive_variables[0];
        const int ghostcell_dim_1_primitive_variables = ghostcell_dims_primitive_variables[1] + 1;
        
        flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computeSoundSpeed(
                sound_speed_y_B,
                primitive_variables_B[0],
                primitive_variables_B[4],
                thermo_properties_const_ptr,
                1,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computeSoundSpeed(
                sound_speed_y_T,
                primitive_variables_T[0],
                primitive_variables_T[4],
                thermo_properties_const_ptr,
                1,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computeInternalEnergy(
                internal_energy_y_B,
                primitive_variables_B[0],
                primitive_variables_B[4],
                thermo_properties_const_ptr,
                1,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computeInternalEnergy(
                internal_energy_y_T,
                primitive_variables_T[0],
                primitive_variables_T[4],
                thermo_properties_const_ptr,
                1,
                domain);
        
        if (compute_velocity)
        {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(velocity);
#endif
            
            /*
             * Get the numbers of cells in each dimension, number of ghost cells of velocity and
             * the pointer to velocity.
             */
            
            const hier::IntVector num_ghosts_velocity = velocity->getGhostCellWidth();
            const hier::IntVector ghostcell_dims_velocity = velocity->getGhostBox().numberCells();
            
            const int num_ghosts_0_velocity = num_ghosts_velocity[0];
            const int num_ghosts_1_velocity = num_ghosts_velocity[1];
            const int num_ghosts_2_velocity = num_ghosts_velocity[2];
            const int ghostcell_dim_0_velocity = ghostcell_dims_velocity[0];
            const int ghostcell_dim_1_velocity = ghostcell_dims_velocity[1] + 1;
            
            Real* v = velocity->getPointer(1, 1);
            
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_flux = (i + num_ghosts_0_convective_flux) +
                            (j + num_ghosts_1_convective_flux)*ghostcell_dim_0_convective_flux +
                            (k + num_ghosts_2_convective_flux)*ghostcell_dim_0_convective_flux*
                                ghostcell_dim_1_convective_flux;
                        
                        const int idx_velocity = (i + num_ghosts_0_velocity) +
                            (j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                            (k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                ghostcell_dim_1_velocity;
                        
                        const int idx = (i + num_ghosts_0_primitive_variables) +
                            (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables +
                            (k + num_ghosts_2_primitive_variables)*ghostcell_dim_0_primitive_variables*
                                ghostcell_dim_1_primitive_variables;
                        
                        computeLocalConvectiveFluxInYDirectionFromPrimitiveVariablesHLLC_HLL3D(
                            F_y.data(),
                            V_y_B.data(),
                            V_y_T.data(),
                            c_y_B,
                            c_y_T,
                            epsilon_y_B,
                            epsilon_y_T,
                            s_y_minus,
                            s_y_plus,
                            s_y_star,
                            Chi_y_star_BT,
                            idx_flux,
                            idx);
                        
                        if (s_y_star > Real(0))
                        {
                            v[idx_velocity] = V_y_B[2][idx] + s_y_minus*(Chi_y_star_BT - Real(1));
                        }
                        else
                        {
                            v[idx_velocity] = V_y_T[2][idx] + s_y_plus*(Chi_y_star_BT - Real(1));
                        }
                    }
                }
            }
        }
        else
        {
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_flux = (i + num_ghosts_0_convective_flux) +
                            (j + num_ghosts_1_convective_flux)*ghostcell_dim_0_convective_flux +
                            (k + num_ghosts_2_convective_flux)*ghostcell_dim_0_convective_flux*
                                ghostcell_dim_1_convective_flux;
                        
                        const int idx = (i + num_ghosts_0_primitive_variables) +
                            (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables +
                            (k + num_ghosts_2_primitive_variables)*ghostcell_dim_0_primitive_variables*
                                ghostcell_dim_1_primitive_variables;
                        
                        computeLocalConvectiveFluxInYDirectionFromPrimitiveVariablesHLLC_HLL3D(
                            F_y.data(),
                            V_y_B.data(),
                            V_y_T.data(),
                            c_y_B,
                            c_y_T,
                            epsilon_y_B,
                            epsilon_y_T,
                            s_y_minus,
                            s_y_plus,
                            s_y_star,
                            Chi_y_star_BT,
                            idx_flux,
                            idx);
                    }
                }
            }
        }
    }
}


/*
 * Compute the convective flux and velocity in the z-direction from primitive variables with
 * HLLC-HLL Riemann solver.
 */
void
FlowModelRiemannSolverSingleSpecies::computeConvectiveFluxAndVelocityInZDirectionFromPrimitiveVariablesHLLC_HLL(
    HAMERS_SHARED_PTR<pdat::SideData<Real> > convective_flux,
    HAMERS_SHARED_PTR<pdat::SideData<Real> > velocity,
    const std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& primitive_variables_B,
    const std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& primitive_variables_F,
    const hier::Box& domain,
    bool compute_velocity) const
{
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    const int num_eqn = flow_model_tmp->getNumberOfEquations();
    
    // Get the box that covers the interior of patch.
    const hier::Box interior_box = convective_flux->getBox();
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_convective_flux = convective_flux->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_convective_flux =
        convective_flux->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_primitive_variables =
        primitive_variables_B[0]->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_primitive_variables =
        primitive_variables_B[0]->getGhostBox().numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        const hier::IntVector num_ghosts_min =
            hier::IntVector::min(num_ghosts_convective_flux, num_ghosts_primitive_variables);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(convective_flux->getGhostBox().contains(domain));
        TBOX_ASSERT(primitive_variables_B[0]->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
    /*
     * Get the equation of state mixing rules and the thermodynamic properties of the species.
     */
    
    const HAMERS_SHARED_PTR<EquationOfStateMixingRules> equation_of_state_mixing_rules =
        flow_model_tmp->getEquationOfStateMixingRules();
    
    const int num_thermo_properties = equation_of_state_mixing_rules->
        getNumberOfSpeciesThermodynamicProperties();
    
    std::vector<Real> thermo_properties;
    std::vector<Real*> thermo_properties_ptr;
    std::vector<const Real*> thermo_properties_const_ptr;
    
    thermo_properties.resize(num_thermo_properties);
    thermo_properties_ptr.reserve(num_thermo_properties);
    thermo_properties_const_ptr.reserve(num_thermo_properties);
    
    for (int ti = 0; ti < num_thermo_properties; ti++)
    {
        thermo_properties_ptr.push_back(&thermo_properties[ti]);
        thermo_properties_const_ptr.push_back(&thermo_properties[ti]);
    }
    
    equation_of_state_mixing_rules->getSpeciesThermodynamicProperties(
        thermo_properties_ptr,
        0);
    
    /*
     * Get the pointers to the side data of convective flux and primitive variables.
     */
    
    std::vector<Real*> F_z;
    F_z.reserve(num_eqn);
    for (int ei = 0; ei < num_eqn; ei++)
    {
        F_z.push_back(convective_flux->getPointer(2, ei));
    }
    
    std::vector<Real*> V_z_B;
    std::vector<Real*> V_z_F;
    V_z_B.reserve(num_eqn);
    V_z_F.reserve(num_eqn);
    for (int ei = 0; ei < num_eqn; ei++)
    {
        V_z_B.push_back(primitive_variables_B[ei]->getPointer(2, 0));
        V_z_F.push_back(primitive_variables_F[ei]->getPointer(2, 0));
    }
    
    /*
     * Allocate temporary data.
     */
    
    hier::IntVector direction_z = hier::IntVector::getZero(d_dim);
    direction_z[2] = 1;
    
    HAMERS_SHARED_PTR<pdat::SideData<Real> > sound_speed_z_B(
        new pdat::SideData<Real>(interior_box, 1, num_ghosts_primitive_variables,
            direction_z));
    
    HAMERS_SHARED_PTR<pdat::SideData<Real> > sound_speed_z_F(
        new pdat::SideData<Real>(interior_box, 1, num_ghosts_primitive_variables,
            direction_z));
    
    HAMERS_SHARED_PTR<pdat::SideData<Real> > internal_energy_z_B(
        new pdat::SideData<Real>(interior_box, 1, num_ghosts_primitive_variables,
            direction_z));

    HAMERS_SHARED_PTR<pdat::SideData<Real> > internal_energy_z_F(
        new pdat::SideData<Real>(interior_box, 1, num_ghosts_primitive_variables,
            direction_z));
    
    /*
     * Get the pointers to the temporary data.
     */
    
    Real* c_z_B = sound_speed_z_B->getPointer(2, 0);
    Real* c_z_F = sound_speed_z_F->getPointer(2, 0);
    
    Real* epsilon_z_B = internal_energy_z_B->getPointer(2, 0);
    Real* epsilon_z_F = internal_energy_z_F->getPointer(2, 0);
    
    Real s_z_minus = Real(0);
    Real s_z_plus  = Real(0);
    Real s_z_star  = Real(0);
    
    Real Chi_z_star_BF = Real(0);
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelRiemannSolverSingleSpecies::"
            << "computeConvectiveFluxAndVelocityInZDirectionFromPrimitiveVariablesHLLC_HLL()\n"
            << "There is no z direction for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelRiemannSolverSingleSpecies::"
            << "computeConvectiveFluxAndVelocityInZDirectionFromPrimitiveVariablesHLLC_HLL()\n"
            << "There is no z direction for two-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the local lower indices, numbers of cells in each dimension and numbers of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_lo_2 = domain_lo[2];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        const int domain_dim_2 = domain_dims[2];
        
        const int num_ghosts_0_convective_flux = num_ghosts_convective_flux[0];
        const int num_ghosts_1_convective_flux = num_ghosts_convective_flux[1];
        const int num_ghosts_2_convective_flux = num_ghosts_convective_flux[2];
        const int ghostcell_dim_0_convective_flux = ghostcell_dims_convective_flux[0];
        const int ghostcell_dim_1_convective_flux = ghostcell_dims_convective_flux[1];
        
        const int num_ghosts_0_primitive_variables = num_ghosts_primitive_variables[0];
        const int num_ghosts_1_primitive_variables = num_ghosts_primitive_variables[1];
        const int num_ghosts_2_primitive_variables = num_ghosts_primitive_variables[2];
        const int ghostcell_dim_0_primitive_variables = ghostcell_dims_primitive_variables[0];
        const int ghostcell_dim_1_primitive_variables = ghostcell_dims_primitive_variables[1];
        
        flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computeSoundSpeed(
                sound_speed_z_B,
                primitive_variables_B[0],
                primitive_variables_B[4],
                thermo_properties_const_ptr,
                2,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computeSoundSpeed(
                sound_speed_z_F,
                primitive_variables_F[0],
                primitive_variables_F[4],
                thermo_properties_const_ptr,
                2,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computeInternalEnergy(
                internal_energy_z_B,
                primitive_variables_B[0],
                primitive_variables_B[4],
                thermo_properties_const_ptr,
                2,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computeInternalEnergy(
                internal_energy_z_F,
                primitive_variables_F[0],
                primitive_variables_F[4],
                thermo_properties_const_ptr,
                2,
                domain);
        
        if (compute_velocity)
        {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(velocity);
#endif
            
            /*
             * Get the numbers of cells in each dimension, number of ghost cells of velocity and
             * the pointer to velocity.
             */
            
            const hier::IntVector num_ghosts_velocity = velocity->getGhostCellWidth();
            const hier::IntVector ghostcell_dims_velocity = velocity->getGhostBox().numberCells();
            
            const int num_ghosts_0_velocity = num_ghosts_velocity[0];
            const int num_ghosts_1_velocity = num_ghosts_velocity[1];
            const int num_ghosts_2_velocity = num_ghosts_velocity[2];
            const int ghostcell_dim_0_velocity = ghostcell_dims_velocity[0];
            const int ghostcell_dim_1_velocity = ghostcell_dims_velocity[1];
            
            Real* w = velocity->getPointer(2, 2);
            
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2 + 1; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_flux = (i + num_ghosts_0_convective_flux) +
                            (j + num_ghosts_1_convective_flux)*ghostcell_dim_0_convective_flux +
                            (k + num_ghosts_2_convective_flux)*ghostcell_dim_0_convective_flux*
                                ghostcell_dim_1_convective_flux;
                        
                        const int idx_velocity = (i + num_ghosts_0_velocity) +
                            (j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                            (k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                ghostcell_dim_1_velocity;
                        
                        const int idx = (i + num_ghosts_0_primitive_variables) +
                            (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables +
                            (k + num_ghosts_2_primitive_variables)*ghostcell_dim_0_primitive_variables*
                                ghostcell_dim_1_primitive_variables;
                        
                        computeLocalConvectiveFluxInZDirectionFromPrimitiveVariablesHLLC_HLL3D(
                            F_z.data(),
                            V_z_B.data(),
                            V_z_F.data(),
                            c_z_B,
                            c_z_F,
                            epsilon_z_B,
                            epsilon_z_F,
                            s_z_minus,
                            s_z_plus,
                            s_z_star,
                            Chi_z_star_BF,
                            idx_flux,
                            idx);
                        
                        if (s_z_star > Real(0))
                        {
                            w[idx_velocity] = V_z_B[3][idx] + s_z_minus*(Chi_z_star_BF - Real(1));
                        }
                        else
                        {
                            w[idx_velocity] = V_z_F[3][idx] + s_z_plus*(Chi_z_star_BF - Real(1));
                        }
                    }
                }
            }
        }
        else
        {
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2 + 1; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_flux = (i + num_ghosts_0_convective_flux) +
                            (j + num_ghosts_1_convective_flux)*ghostcell_dim_0_convective_flux +
                            (k + num_ghosts_2_convective_flux)*ghostcell_dim_0_convective_flux*
                                ghostcell_dim_1_convective_flux;
                        
                        const int idx = (i + num_ghosts_0_primitive_variables) +
                            (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables +
                            (k + num_ghosts_2_primitive_variables)*ghostcell_dim_0_primitive_variables*
                                ghostcell_dim_1_primitive_variables;
                        
                        computeLocalConvectiveFluxInZDirectionFromPrimitiveVariablesHLLC_HLL3D(
                            F_z.data(),
                            V_z_B.data(),
                            V_z_F.data(),
                            c_z_B,
                            c_z_F,
                            epsilon_z_B,
                            epsilon_z_F,
                            s_z_minus,
                            s_z_plus,
                            s_z_star,
                            Chi_z_star_BF,
                            idx_flux,
                            idx);
                    }
                }
            }
        }
    }
}
