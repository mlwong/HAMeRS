#include "flow/flow_models/five-eqn_Allaire/FlowModelRiemannSolverFiveEqnAllaire.hpp"

#define EPSILON HAMERS_EPSILON


/*
 * Compute the local convective flux in the x-direction from conservative variables with
 * 1D HLLC-HLL Riemann solver.
 */
static inline __attribute__((always_inline)) void
computeLocalConvectiveFluxInXDirectionFromConservativeVariablesHLLC_HLL1D(
    double** F_x,
    double** Q_x_L,
    double** Q_x_R,
    double* rho_x_L,
    double* rho_x_R,
    double* p_x_L,
    double* p_x_R,
    double* c_x_L,
    double* c_x_R,
    double& u_x_L,
    double& u_x_R,
    double& s_x_minus,
    double& s_x_plus,
    double& s_x_star,
    double& Chi_x_star_LR,
    const int& idx_flux,
    const int& idx,
    const int& num_species,
    const int& num_eqn)
{
    u_x_L = Q_x_L[num_species][idx]/rho_x_L[idx];
    u_x_R = Q_x_R[num_species][idx]/rho_x_R[idx];
    
    const double u_x_average = double(1)/double(2)*(u_x_L + u_x_R);
    const double c_x_average = double(1)/double(2)*(c_x_L[idx] + c_x_R[idx]);
    
    const double s_x_L = fmin(u_x_average - c_x_average, u_x_L - c_x_L[idx]);
    const double s_x_R = fmax(u_x_average + c_x_average, u_x_R + c_x_R[idx]);
    
    s_x_minus = fmin(double(0), s_x_L);
    s_x_plus  = fmax(double(0), s_x_R);
    
    s_x_star = (p_x_R[idx] - p_x_L[idx] +
        Q_x_L[num_species][idx]*(s_x_L - u_x_L) - Q_x_R[num_species][idx]*(s_x_R - u_x_R))/
        (rho_x_L[idx]*(s_x_L - u_x_L) - rho_x_R[idx]*(s_x_R - u_x_R));
    
    double Q_x_star_LR[num_eqn];
    double F_x_LR[num_eqn];
    
    if (s_x_star > double(0))
    {
        Chi_x_star_LR = (s_x_L - u_x_L)/(s_x_L - s_x_star);
        
        for (int si = 0; si < num_species; si++)
        {
            Q_x_star_LR[si] = Chi_x_star_LR*Q_x_L[si][idx];
        }
        Q_x_star_LR[num_species] = Chi_x_star_LR*rho_x_L[idx]*s_x_star;
        Q_x_star_LR[num_species + 1] = Chi_x_star_LR*(Q_x_L[num_species + 1][idx] +
            (s_x_star - u_x_L)*(rho_x_L[idx]*s_x_star + p_x_L[idx]/(s_x_L - u_x_L)));
        for (int si = 0; si < num_species - 1; si++)
        {
            Q_x_star_LR[num_species + 2 + si] = Chi_x_star_LR*Q_x_L[num_species + 2 + si][idx];
        }
        
        for (int si = 0; si < num_species; si++)
        {
            F_x_LR[si] = u_x_L*Q_x_L[si][idx];
        }
        F_x_LR[num_species] = u_x_L*Q_x_L[num_species][idx] + p_x_L[idx];
        F_x_LR[num_species + 1] = u_x_L*(Q_x_L[num_species + 1][idx] + p_x_L[idx]);
        for (int si = 0; si < num_species - 1; si++)
        {
            F_x_LR[num_species + 2 + si] = u_x_L*(Q_x_L[num_species + 2 + si][idx]);
        }
        
        for (int ei = 0; ei < num_eqn; ei++)
        {
            F_x[ei][idx_flux] = F_x_LR[ei] + s_x_minus*(Q_x_star_LR[ei] - Q_x_L[ei][idx]);
        }
    }
    else
    {
        Chi_x_star_LR = (s_x_R - u_x_R)/(s_x_R - s_x_star);
        
        for (int si = 0; si < num_species; si++)
        {
            Q_x_star_LR[si] = Chi_x_star_LR*Q_x_R[si][idx];
        }
        Q_x_star_LR[num_species] = Chi_x_star_LR*rho_x_R[idx]*s_x_star;
        Q_x_star_LR[num_species + 1] = Chi_x_star_LR*(Q_x_R[num_species + 1][idx] +
            (s_x_star - u_x_R)*(rho_x_R[idx]*s_x_star + p_x_R[idx]/(s_x_R - u_x_R)));
        for (int si = 0; si < num_species - 1; si++)
        {
            Q_x_star_LR[num_species + 2 + si] = Chi_x_star_LR*Q_x_R[num_species + 2 + si][idx];
        }
        
        for (int si = 0; si < num_species; si++)
        {
            F_x_LR[si] = u_x_R*Q_x_R[si][idx];
        }
        F_x_LR[num_species] = u_x_R*Q_x_R[num_species][idx] + p_x_R[idx];
        F_x_LR[num_species + 1] = u_x_R*(Q_x_R[num_species + 1][idx] + p_x_R[idx]);
        for (int si = 0; si < num_species - 1; si++)
        {
            F_x_LR[num_species + 2 + si] = u_x_R*(Q_x_R[num_species + 2 + si][idx]);
        }
        
        for (int ei = 0; ei < num_eqn; ei++)
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
    double** F_x,
    double** Q_x_L,
    double** Q_x_R,
    double* rho_x_L,
    double* rho_x_R,
    double* p_x_L,
    double* p_x_R,
    double* c_x_L,
    double* c_x_R,
    double& u_x_L,
    double& u_x_R,
    double& s_x_minus,
    double& s_x_plus,
    double& s_x_star,
    double& Chi_x_star_LR,
    const int& idx_flux,
    const int& idx,
    const int& num_species,
    const int& num_eqn)
{
    u_x_L = Q_x_L[num_species][idx]/rho_x_L[idx];
    u_x_R = Q_x_R[num_species][idx]/rho_x_R[idx];
    
    const double v_x_L = Q_x_L[num_species + 1][idx]/rho_x_L[idx];
    const double v_x_R = Q_x_R[num_species + 1][idx]/rho_x_R[idx];
    
    const double u_x_average = double(1)/double(2)*(u_x_L + u_x_R);
    const double c_x_average = double(1)/double(2)*(c_x_L[idx] + c_x_R[idx]);
    
    const double s_x_L = fmin(u_x_average - c_x_average, u_x_L - c_x_L[idx]);
    const double s_x_R = fmax(u_x_average + c_x_average, u_x_R + c_x_R[idx]);
    
    s_x_minus = fmin(double(0), s_x_L);
    s_x_plus  = fmax(double(0), s_x_R);
    
    s_x_star = (p_x_R[idx] - p_x_L[idx] +
        Q_x_L[num_species][idx]*(s_x_L - u_x_L) - Q_x_R[num_species][idx]*(s_x_R - u_x_R))/
        (rho_x_L[idx]*(s_x_L - u_x_L) - rho_x_R[idx]*(s_x_R - u_x_R));
    
    double F_x_L[num_eqn];
    double F_x_R[num_eqn];
    double F_x_HLL[2*num_species];
    double F_x_HLLC[num_eqn];
    double Q_x_star_LR[num_eqn];
    
    for (int si = 0; si < num_species; si++)
    {
        F_x_L[si] = u_x_L*Q_x_L[si][idx];
    }
    F_x_L[num_species] = u_x_L*Q_x_L[num_species][idx] + p_x_L[idx];
    F_x_L[num_species + 1] = u_x_L*Q_x_L[num_species + 1][idx];
    F_x_L[num_species + 2] = u_x_L*(Q_x_L[num_species + 2][idx] + p_x_L[idx]);
    for (int si = 0; si < num_species - 1; si++)
    {
        F_x_L[num_species + 3 + si] = u_x_L*(Q_x_L[num_species + 3 + si][idx]);
    }
    
    for (int si = 0; si < num_species; si++)
    {
        F_x_R[si] = u_x_R*Q_x_R[si][idx];
    }
    F_x_R[num_species] = u_x_R*Q_x_R[num_species][idx] + p_x_R[idx];
    F_x_R[num_species + 1] = u_x_R*Q_x_R[num_species + 1][idx];
    F_x_R[num_species + 2] = u_x_R*(Q_x_R[num_species + 2][idx] + p_x_R[idx]);
    for (int si = 0; si < num_species - 1; si++)
    {
        F_x_R[num_species + 3 + si] = u_x_R*(Q_x_R[num_species + 3 + si][idx]);
    }
    
    for (int si = 0; si < num_species; si++)
    {
        F_x_HLL[si] = (s_x_R*F_x_L[si] - s_x_L*F_x_R[si] + s_x_R*s_x_L*(Q_x_R[si][idx] - Q_x_L[si][idx]))/
            (s_x_R - s_x_L);
    }
    F_x_HLL[num_species] = (s_x_R*F_x_L[num_species + 1] - s_x_L*F_x_R[num_species + 1] + s_x_R*s_x_L*
        (Q_x_R[num_species + 1][idx] - Q_x_L[num_species + 1][idx]))/(s_x_R - s_x_L);
    for (int si = 0; si < num_species - 1; si++)
    {
        F_x_HLL[num_species + 1 + si] = (s_x_R*F_x_L[num_species + 3 + si] - s_x_L*F_x_R[num_species + 3 + si] +
            s_x_R*s_x_L*(Q_x_R[num_species + 3 + si][idx] - Q_x_L[num_species + 3 + si][idx]))/(s_x_R - s_x_L);
    }
    
    if (s_x_L > double(0))
    {
        for (int si = 0; si < num_species; si++)
        {
            F_x_HLL[si] = F_x_L[si];
        }
        F_x_HLL[num_species] = F_x_L[num_species + 1];
        for (int si = 0; si < num_species - 1; si++)
        {
            F_x_HLL[num_species + 1 + si] = F_x_L[num_species + 3 + si];
        }
    }
    
    if (s_x_R < double(0))
    {
        for (int si = 0; si < num_species; si++)
        {
            F_x_HLL[si] = F_x_R[si];
        }
        F_x_HLL[num_species] = F_x_R[num_species + 1];
        for (int si = 0; si < num_species - 1; si++)
        {
            F_x_HLL[num_species + 1 + si] = F_x_R[num_species + 3 + si];
        }
    }
    
    if (s_x_star > double(0))
    {
        Chi_x_star_LR = (s_x_L - u_x_L)/(s_x_L - s_x_star);
        
        for (int si = 0; si < num_species; si++)
        {
            Q_x_star_LR[si] = Chi_x_star_LR*Q_x_L[si][idx];
        }
        Q_x_star_LR[num_species] = Chi_x_star_LR*rho_x_L[idx]*s_x_star;
        Q_x_star_LR[num_species + 1] = Chi_x_star_LR*Q_x_L[num_species + 1][idx];
        Q_x_star_LR[num_species + 2] = Chi_x_star_LR*(Q_x_L[num_species + 2][idx] +
            (s_x_star - u_x_L)*(rho_x_L[idx]*s_x_star + p_x_L[idx]/(s_x_L - u_x_L)));
        for (int si = 0; si < num_species - 1; si++)
        {
            Q_x_star_LR[num_species + 3 + si] = Chi_x_star_LR*Q_x_L[num_species + 3 + si][idx];
        }
        
        for (int ei = 0; ei < num_eqn; ei++)
        {
            F_x_HLLC[ei] = F_x_L[ei] + s_x_minus*(Q_x_star_LR[ei] - Q_x_L[ei][idx]);
        }
    }
    else
    {
        Chi_x_star_LR = (s_x_R - u_x_R)/(s_x_R - s_x_star);
        
        for (int si = 0; si < num_species; si++)
        {
            Q_x_star_LR[si] = Chi_x_star_LR*Q_x_R[si][idx];
        }
        Q_x_star_LR[num_species] = Chi_x_star_LR*rho_x_R[idx]*s_x_star;
        Q_x_star_LR[num_species + 1] = Chi_x_star_LR*Q_x_R[num_species + 1][idx];
        Q_x_star_LR[num_species + 2] = Chi_x_star_LR*(Q_x_R[num_species + 2][idx] +
            (s_x_star - u_x_R)*(rho_x_R[idx]*s_x_star + p_x_R[idx]/(s_x_R - u_x_R)));
        for (int si = 0; si < num_species - 1; si++)
        {
            Q_x_star_LR[num_species + 3 + si] = Chi_x_star_LR*Q_x_R[num_species + 3 + si][idx];
        }
        
        for (int ei = 0; ei < num_eqn; ei++)
        {
            F_x_HLLC[ei] = F_x_R[ei] + s_x_plus*(Q_x_star_LR[ei] - Q_x_R[ei][idx]);
        }
    }
    
    /*
     * Calulate the weights beta for hybridization.
     */
    
    const double u_x_diff = u_x_R - u_x_L;
    const double v_x_diff = v_x_R - v_x_L;
    const double vel_mag = sqrt(u_x_diff*u_x_diff + v_x_diff*v_x_diff);
    
    double alpha_1, alpha_2;
    if (vel_mag < EPSILON)
    {
        alpha_1 = double(1);
        alpha_2 = double(0);
    }
    else
    {
        alpha_1 = fabs(u_x_diff)/vel_mag;
        alpha_2 = sqrt(double(1) - alpha_1*alpha_1);
    }
    
    const double beta_1 = double(1)/double(2)*(double(1) + alpha_1/(alpha_1 + alpha_2));
    const double beta_2 = double(1) - beta_1;
    
    for (int si = 0; si < num_species; si++)
    {
        F_x[si][idx_flux] = beta_1*F_x_HLLC[si] + beta_2*F_x_HLL[si];
    }
    F_x[num_species][idx_flux] = F_x_HLLC[num_species];
    F_x[num_species + 1][idx_flux] = beta_1*F_x_HLLC[num_species + 1] + beta_2*F_x_HLL[num_species];
    F_x[num_species + 2][idx_flux] = F_x_HLLC[num_species + 2];
    for (int si = 0; si < num_species - 1; si++)
    {
        F_x[num_species + 3 + si][idx_flux] = beta_1*F_x_HLLC[num_species + 3 + si] +
            beta_2*F_x_HLL[num_species + 1 + si];
    }
}


/*
 * Compute the local convective flux in the x-direction from conservative variables with
 * 3D HLLC-HLL Riemann solver.
 */
static inline __attribute__((always_inline)) void
computeLocalConvectiveFluxInXDirectionFromConservativeVariablesHLLC_HLL3D(
    double** F_x,
    double** Q_x_L,
    double** Q_x_R,
    double* rho_x_L,
    double* rho_x_R,
    double* p_x_L,
    double* p_x_R,
    double* c_x_L,
    double* c_x_R,
    double& u_x_L,
    double& u_x_R,
    double& s_x_minus,
    double& s_x_plus,
    double& s_x_star,
    double& Chi_x_star_LR,
    const int& idx_flux,
    const int& idx,
    const int& num_species,
    const int& num_eqn)
{
    u_x_L = Q_x_L[num_species][idx]/rho_x_L[idx];
    u_x_R = Q_x_R[num_species][idx]/rho_x_R[idx];
    
    const double v_x_L = Q_x_L[num_species + 1][idx]/rho_x_L[idx];
    const double v_x_R = Q_x_R[num_species + 1][idx]/rho_x_R[idx];
    
    const double w_x_L = Q_x_L[num_species + 2][idx]/rho_x_L[idx];
    const double w_x_R = Q_x_R[num_species + 2][idx]/rho_x_R[idx];
    
    const double u_x_average = double(1)/double(2)*(u_x_L + u_x_R);
    const double c_x_average = double(1)/double(2)*(c_x_L[idx] + c_x_R[idx]);
    
    const double s_x_L = fmin(u_x_average - c_x_average, u_x_L - c_x_L[idx]);
    const double s_x_R = fmax(u_x_average + c_x_average, u_x_R + c_x_R[idx]);
    
    s_x_minus = fmin(double(0), s_x_L);
    s_x_plus  = fmax(double(0), s_x_R);
    
    s_x_star = (p_x_R[idx] - p_x_L[idx] +
        Q_x_L[num_species][idx]*(s_x_L - u_x_L) - Q_x_R[num_species][idx]*(s_x_R - u_x_R))/
        (rho_x_L[idx]*(s_x_L - u_x_L) - rho_x_R[idx]*(s_x_R - u_x_R));
    
    double F_x_L[num_eqn];
    double F_x_R[num_eqn];
    double F_x_HLL[2*num_species + 1];
    double F_x_HLLC[num_eqn];
    double Q_x_star_LR[num_eqn];
    
    for (int si = 0; si < num_species; si++)
    {
        F_x_L[si] = u_x_L*Q_x_L[si][idx];
    }
    F_x_L[num_species] = u_x_L*Q_x_L[num_species][idx] + p_x_L[idx];
    F_x_L[num_species + 1] = u_x_L*Q_x_L[num_species + 1][idx];
    F_x_L[num_species + 2] = u_x_L*Q_x_L[num_species + 2][idx];
    F_x_L[num_species + 3] = u_x_L*(Q_x_L[num_species + 3][idx] + p_x_L[idx]);
    for (int si = 0; si < num_species - 1; si++)
    {
        F_x_L[num_species + 4 + si] = u_x_L*(Q_x_L[num_species + 4 + si][idx]);
    }
    
    for (int si = 0; si < num_species; si++)
    {
        F_x_R[si] = u_x_R*Q_x_R[si][idx];
    }
    F_x_R[num_species] = u_x_R*Q_x_R[num_species][idx] + p_x_R[idx];
    F_x_R[num_species + 1] = u_x_R*Q_x_R[num_species + 1][idx];
    F_x_R[num_species + 2] = u_x_R*Q_x_R[num_species + 2][idx];
    F_x_R[num_species + 3] = u_x_R*(Q_x_R[num_species + 3][idx] + p_x_R[idx]);
    for (int si = 0; si < num_species - 1; si++)
    {
        F_x_R[num_species + 4 + si] = u_x_R*(Q_x_R[num_species + 4 + si][idx]);
    }
    
    for (int si = 0; si < num_species; si++)
    {
        F_x_HLL[si] = (s_x_R*F_x_L[si] - s_x_L*F_x_R[si] + s_x_R*s_x_L*(Q_x_R[si][idx] - Q_x_L[si][idx]))/
            (s_x_R - s_x_L);
    }
    F_x_HLL[num_species] = (s_x_R*F_x_L[num_species + 1] - s_x_L*F_x_R[num_species + 1] + s_x_R*s_x_L*
        (Q_x_R[num_species + 1][idx] - Q_x_L[num_species + 1][idx]))/(s_x_R - s_x_L);
    F_x_HLL[num_species + 1] = (s_x_R*F_x_L[num_species + 2] - s_x_L*F_x_R[num_species + 2] + s_x_R*s_x_L*
        (Q_x_R[num_species + 2][idx] - Q_x_L[num_species + 2][idx]))/(s_x_R - s_x_L);
    for (int si = 0; si < num_species - 1; si++)
    {
        F_x_HLL[num_species + 2 + si] = (s_x_R*F_x_L[num_species + 4 + si] - s_x_L*F_x_R[num_species + 4 + si] +
            s_x_R*s_x_L*(Q_x_R[num_species + 4 + si][idx] - Q_x_L[num_species + 4 + si][idx]))/(s_x_R - s_x_L);
    }
    
    if (s_x_L > double(0))
    {
        for (int si = 0; si < num_species; si++)
        {
            F_x_HLL[si] = F_x_L[si];
        }
        F_x_HLL[num_species] = F_x_L[num_species + 1];
        F_x_HLL[num_species + 1] = F_x_L[num_species + 2];
        for (int si = 0; si < num_species - 1; si++)
        {
            F_x_HLL[num_species + 2 + si] = F_x_L[num_species + 4 + si];
        }
    }
    
    if (s_x_R < double(0))
    {
        for (int si = 0; si < num_species; si++)
        {
            F_x_HLL[si] = F_x_R[si];
        }
        F_x_HLL[num_species] = F_x_R[num_species + 1];
        F_x_HLL[num_species + 1] = F_x_R[num_species + 2];
        for (int si = 0; si < num_species - 1; si++)
        {
            F_x_HLL[num_species + 2 + si] = F_x_R[num_species + 4 + si];
        }
    }
    
    if (s_x_star > double(0))
    {
        Chi_x_star_LR = (s_x_L - u_x_L)/(s_x_L - s_x_star);
        
        for (int si = 0; si < num_species; si++)
        {
            Q_x_star_LR[si] = Chi_x_star_LR*Q_x_L[si][idx];
        }
        Q_x_star_LR[num_species] = Chi_x_star_LR*rho_x_L[idx]*s_x_star;
        Q_x_star_LR[num_species + 1] = Chi_x_star_LR*Q_x_L[num_species + 1][idx];
        Q_x_star_LR[num_species + 2] = Chi_x_star_LR*Q_x_L[num_species + 2][idx];
        Q_x_star_LR[num_species + 3] = Chi_x_star_LR*(Q_x_L[num_species + 3][idx] +
            (s_x_star - u_x_L)*(rho_x_L[idx]*s_x_star + p_x_L[idx]/(s_x_L - u_x_L)));
        for (int si = 0; si < num_species - 1; si++)
        {
            Q_x_star_LR[num_species + 4 + si] = Chi_x_star_LR*Q_x_L[num_species + 4 + si][idx];
        }
        
        for (int ei = 0; ei < num_eqn; ei++)
        {
            F_x_HLLC[ei] = F_x_L[ei] + s_x_minus*(Q_x_star_LR[ei] - Q_x_L[ei][idx]);
        }
    }
    else
    {
        Chi_x_star_LR = (s_x_R - u_x_R)/(s_x_R - s_x_star);
        
        for (int si = 0; si < num_species; si++)
        {
            Q_x_star_LR[si] = Chi_x_star_LR*Q_x_R[si][idx];
        }
        Q_x_star_LR[num_species] = Chi_x_star_LR*rho_x_R[idx]*s_x_star;
        Q_x_star_LR[num_species + 1] = Chi_x_star_LR*Q_x_R[num_species + 1][idx];
        Q_x_star_LR[num_species + 2] = Chi_x_star_LR*Q_x_R[num_species + 2][idx];
        Q_x_star_LR[num_species + 3] = Chi_x_star_LR*(Q_x_R[num_species + 3][idx] +
            (s_x_star - u_x_R)*(rho_x_R[idx]*s_x_star + p_x_R[idx]/(s_x_R - u_x_R)));
        for (int si = 0; si < num_species - 1; si++)
        {
            Q_x_star_LR[num_species + 4 + si] = Chi_x_star_LR*Q_x_R[num_species + 4 + si][idx];
        }
        
        for (int ei = 0; ei < num_eqn; ei++)
        {
            F_x_HLLC[ei] = F_x_R[ei] + s_x_plus*(Q_x_star_LR[ei] - Q_x_R[ei][idx]);
        }
    }
    
    /*
     * Calulate the weights beta for hybridization.
     */
    
    const double u_x_diff = u_x_R - u_x_L;
    const double v_x_diff = v_x_R - v_x_L;
    const double w_x_diff = w_x_R - w_x_L;
    const double vel_mag = sqrt(u_x_diff*u_x_diff + v_x_diff*v_x_diff + w_x_diff*w_x_diff);
    
    double alpha_1, alpha_2;
    if (vel_mag < EPSILON)
    {
        alpha_1 = double(1);
        alpha_2 = double(0);
    }
    else
    {
        alpha_1 = fabs(u_x_diff)/vel_mag;
        alpha_2 = sqrt(double(1) - alpha_1*alpha_1);
    }
    
    const double beta_1 = double(1)/double(2)*(double(1) + alpha_1/(alpha_1 + alpha_2));
    const double beta_2 = double(1) - beta_1;
    
    for (int si = 0; si < num_species; si++)
    {
        F_x[si][idx_flux] = beta_1*F_x_HLLC[si] + beta_2*F_x_HLL[si];
    }
    F_x[num_species][idx_flux] = F_x_HLLC[num_species];
    F_x[num_species + 1][idx_flux] = beta_1*F_x_HLLC[num_species + 1] + beta_2*F_x_HLL[num_species];
    F_x[num_species + 2][idx_flux] = beta_1*F_x_HLLC[num_species + 2] + beta_2*F_x_HLL[num_species + 1];
    F_x[num_species + 3][idx_flux] = F_x_HLLC[num_species + 3];
    for (int si = 0; si < num_species - 1; si++)
    {
        F_x[num_species + 4 + si][idx_flux] = beta_1*F_x_HLLC[num_species + 4 + si] +
            beta_2*F_x_HLL[num_species + 2 + si];
    }
}


/*
 * Compute the local convective flux in the y-direction from conservative variables with
 * 2D HLLC-HLL Riemann solver.
 */
static inline __attribute__((always_inline)) void
computeLocalConvectiveFluxInYDirectionFromConservativeVariablesHLLC_HLL2D(
    double** F_y,
    double** Q_y_B,
    double** Q_y_T,
    double* rho_y_B,
    double* rho_y_T,
    double* p_y_B,
    double* p_y_T,
    double* c_y_B,
    double* c_y_T,
    double& v_y_B,
    double& v_y_T,
    double& s_y_minus,
    double& s_y_plus,
    double& s_y_star,
    double& Chi_y_star_BT,
    const int& idx_flux,
    const int& idx,
    const int& num_species,
    const int& num_eqn)
{
    v_y_B = Q_y_B[num_species + 1][idx]/rho_y_B[idx];
    v_y_T = Q_y_T[num_species + 1][idx]/rho_y_T[idx];
    
    const double u_y_B = Q_y_B[num_species][idx]/rho_y_B[idx];
    const double u_y_T = Q_y_T[num_species][idx]/rho_y_T[idx];
    
    const double v_y_average = double(1)/double(2)*(v_y_B + v_y_T);
    const double c_y_average = double(1)/double(2)*(c_y_B[idx] + c_y_T[idx]);
    
    const double s_y_B = fmin(v_y_average - c_y_average, v_y_B - c_y_B[idx]);
    const double s_y_T = fmax(v_y_average + c_y_average, v_y_T + c_y_T[idx]);
    
    s_y_minus = fmin(double(0), s_y_B);
    s_y_plus  = fmax(double(0), s_y_T);
    
    s_y_star = (p_y_T[idx] - p_y_B[idx] +
        Q_y_B[num_species + 1][idx]*(s_y_B - v_y_B) - Q_y_T[num_species + 1][idx]*(s_y_T - v_y_T))/
        (rho_y_B[idx]*(s_y_B - v_y_B) - rho_y_T[idx]*(s_y_T - v_y_T));
    
    double F_y_B[num_eqn];
    double F_y_T[num_eqn];
    double F_y_HLL[2*num_species];
    double F_y_HLLC[num_eqn];
    double Q_y_star_BT[num_eqn];
    
    for (int si = 0; si < num_species; si++)
    {
        F_y_B[si] = v_y_B*Q_y_B[si][idx];
    }
    F_y_B[num_species] = v_y_B*Q_y_B[num_species][idx];
    F_y_B[num_species + 1] = v_y_B*Q_y_B[num_species + 1][idx] + p_y_B[idx];
    F_y_B[num_species + 2] = v_y_B*(Q_y_B[num_species + 2][idx] + p_y_B[idx]);
    for (int si = 0; si < num_species - 1; si++)
    {
        F_y_B[num_species + 3 + si] = v_y_B*(Q_y_B[num_species + 3 + si][idx]);
    }
    
    for (int si = 0; si < num_species; si++)
    {
        F_y_T[si] = v_y_T*Q_y_T[si][idx];
    }
    F_y_T[num_species] = v_y_T*Q_y_T[num_species][idx];
    F_y_T[num_species + 1] = v_y_T*Q_y_T[num_species + 1][idx] + p_y_T[idx];
    F_y_T[num_species + 2] = v_y_T*(Q_y_T[num_species + 2][idx] + p_y_T[idx]);
    for (int si = 0; si < num_species - 1; si++)
    {
        F_y_T[num_species + 3 + si] = v_y_T*(Q_y_T[num_species + 3 + si][idx]);
    }
    
    for (int si = 0; si < num_species; si++)
    {
        F_y_HLL[si] = (s_y_T*F_y_B[si] - s_y_B*F_y_T[si] + s_y_T*s_y_B*(Q_y_T[si][idx] - Q_y_B[si][idx]))/
            (s_y_T - s_y_B);
    }
    F_y_HLL[num_species] = (s_y_T*F_y_B[num_species] - s_y_B*F_y_T[num_species] + s_y_T*s_y_B*
        (Q_y_T[num_species][idx] - Q_y_B[num_species][idx]))/(s_y_T - s_y_B);
    for (int si = 0; si < num_species - 1; si++)
    {
        F_y_HLL[num_species + 1 + si] = (s_y_T*F_y_B[num_species + 3 + si] - s_y_B*F_y_T[num_species + 3 + si] +
            s_y_T*s_y_B*(Q_y_T[num_species + 3 + si][idx] - Q_y_B[num_species + 3 + si][idx]))/(s_y_T - s_y_B);
    }
    
    if (s_y_B > double(0))
    {
        for (int si = 0; si < num_species; si++)
        {
            F_y_HLL[si] = F_y_B[si];
        }
        F_y_HLL[num_species] = F_y_B[num_species];
        for (int si = 0; si < num_species - 1; si++)
        {
            F_y_HLL[num_species + 1 + si] = F_y_B[num_species + 3 + si];
        }
    }
    
    if (s_y_T < double(0))
    {
        for (int si = 0; si < num_species; si++)
        {
            F_y_HLL[si] = F_y_T[si];
        }
        F_y_HLL[num_species] = F_y_T[num_species];
        for (int si = 0; si < num_species - 1; si++)
        {
            F_y_HLL[num_species + 1 + si] = F_y_T[num_species + 3 + si];
        }
    }
    
    if (s_y_star > double(0))
    {
        Chi_y_star_BT = (s_y_B - v_y_B)/(s_y_B - s_y_star);
        
        for (int si = 0; si < num_species; si++)
        {
            Q_y_star_BT[si] = Chi_y_star_BT*Q_y_B[si][idx];
        }
        Q_y_star_BT[num_species] = Chi_y_star_BT*Q_y_B[num_species][idx];
        Q_y_star_BT[num_species + 1] = Chi_y_star_BT*rho_y_B[idx]*s_y_star;
        Q_y_star_BT[num_species + 2] = Chi_y_star_BT*(Q_y_B[num_species + 2][idx] +
            (s_y_star - v_y_B)*(rho_y_B[idx]*s_y_star + p_y_B[idx]/(s_y_B - v_y_B)));
        for (int si = 0; si < num_species - 1; si++)
        {
            Q_y_star_BT[num_species + 3 + si] = Chi_y_star_BT*Q_y_B[num_species + 3 + si][idx];
        }
        
        for (int ei = 0; ei < num_eqn; ei++)
        {
            F_y_HLLC[ei] = F_y_B[ei] + s_y_minus*(Q_y_star_BT[ei] - Q_y_B[ei][idx]);
        }
    }
    else
    {
        Chi_y_star_BT = (s_y_T - v_y_T)/(s_y_T - s_y_star);
        
        for (int si = 0; si < num_species; si++)
        {
            Q_y_star_BT[si] = Chi_y_star_BT*Q_y_T[si][idx];
        }
        Q_y_star_BT[num_species] = Chi_y_star_BT*Q_y_T[num_species][idx];
        Q_y_star_BT[num_species + 1] = Chi_y_star_BT*rho_y_T[idx]*s_y_star;
        Q_y_star_BT[num_species + 2] = Chi_y_star_BT*(Q_y_T[num_species + 2][idx] +
            (s_y_star - v_y_T)*(rho_y_T[idx]*s_y_star + p_y_T[idx]/(s_y_T - v_y_T)));
        for (int si = 0; si < num_species - 1; si++)
        {
            Q_y_star_BT[num_species + 3 + si] = Chi_y_star_BT*Q_y_T[num_species + 3 + si][idx];
        }
        
        for (int ei = 0; ei < num_eqn; ei++)
        {
            F_y_HLLC[ei] = F_y_T[ei] + s_y_plus*(Q_y_star_BT[ei] - Q_y_T[ei][idx]);
        }
    }
    
    /*
     * Calulate the weights beta for hybridization.
     */
    
    const double u_y_diff = u_y_T - u_y_B;
    const double v_y_diff = v_y_T - v_y_B;
    const double vel_mag = sqrt(u_y_diff*u_y_diff + v_y_diff*v_y_diff);
    
    double alpha_1, alpha_2;
    if (vel_mag < EPSILON)
    {
        alpha_1 = double(1);
        alpha_2 = double(0);
    }
    else
    {
        alpha_1 = fabs(v_y_diff)/vel_mag;
        alpha_2 = sqrt(double(1) - alpha_1*alpha_1);
    }
    
    const double beta_1 = double(1)/double(2)*(double(1) + alpha_1/(alpha_1 + alpha_2));
    const double beta_2 = double(1) - beta_1;
    
    for (int si = 0; si < num_species; si++)
    {
        F_y[si][idx_flux] = beta_1*F_y_HLLC[si] + beta_2*F_y_HLL[si];
    }
    F_y[num_species][idx_flux] = beta_1*F_y_HLLC[num_species] + beta_2*F_y_HLL[num_species];
    F_y[num_species + 1][idx_flux] = F_y_HLLC[num_species + 1];
    F_y[num_species + 2][idx_flux] = F_y_HLLC[num_species + 2];
    for (int si = 0; si < num_species - 1; si++)
    {
        F_y[num_species + 3 + si][idx_flux] = beta_1*F_y_HLLC[num_species + 3 + si] +
            beta_2*F_y_HLL[num_species + 1 + si];
    }
}


/*
 * Compute the local convective flux in the y-direction from conservative variables with
 * 3D HLLC-HLL Riemann solver.
 */
static inline __attribute__((always_inline)) void
computeLocalConvectiveFluxInYDirectionFromConservativeVariablesHLLC_HLL3D(
    double** F_y,
    double** Q_y_B,
    double** Q_y_T,
    double* rho_y_B,
    double* rho_y_T,
    double* p_y_B,
    double* p_y_T,
    double* c_y_B,
    double* c_y_T,
    double& v_y_B,
    double& v_y_T,
    double& s_y_minus,
    double& s_y_plus,
    double& s_y_star,
    double& Chi_y_star_BT,
    const int& idx_flux,
    const int& idx,
    const int& num_species,
    const int& num_eqn)
{
    v_y_B = Q_y_B[num_species + 1][idx]/rho_y_B[idx];
    v_y_T = Q_y_T[num_species + 1][idx]/rho_y_T[idx];
    
    const double u_y_B = Q_y_B[num_species][idx]/rho_y_B[idx];
    const double u_y_T = Q_y_T[num_species][idx]/rho_y_T[idx];
    
    const double w_y_B = Q_y_B[num_species + 2][idx]/rho_y_B[idx];
    const double w_y_T = Q_y_T[num_species + 2][idx]/rho_y_T[idx];
    
    const double v_y_average = double(1)/double(2)*(v_y_B + v_y_T);
    const double c_y_average = double(1)/double(2)*(c_y_B[idx] + c_y_T[idx]);
    
    const double s_y_B = fmin(v_y_average - c_y_average, v_y_B - c_y_B[idx]);
    const double s_y_T = fmax(v_y_average + c_y_average, v_y_T + c_y_T[idx]);
    
    s_y_minus = fmin(double(0), s_y_B);
    s_y_plus  = fmax(double(0), s_y_T);
    
    s_y_star = (p_y_T[idx] - p_y_B[idx] +
        Q_y_B[num_species + 1][idx]*(s_y_B - v_y_B) - Q_y_T[num_species + 1][idx]*(s_y_T - v_y_T))/
        (rho_y_B[idx]*(s_y_B - v_y_B) - rho_y_T[idx]*(s_y_T - v_y_T));
    
    double F_y_B[num_eqn];
    double F_y_T[num_eqn];
    double F_y_HLL[2*num_species + 1];
    double F_y_HLLC[num_eqn];
    double Q_y_star_BT[num_eqn];
    
    for (int si = 0; si < num_species; si++)
    {
        F_y_B[si] = v_y_B*Q_y_B[si][idx];
    }
    F_y_B[num_species] = v_y_B*Q_y_B[num_species][idx];
    F_y_B[num_species + 1] = v_y_B*Q_y_B[num_species + 1][idx] + p_y_B[idx];
    F_y_B[num_species + 2] = v_y_B*Q_y_B[num_species + 2][idx];
    F_y_B[num_species + 3] = v_y_B*(Q_y_B[num_species + 3][idx] + p_y_B[idx]);
    for (int si = 0; si < num_species - 1; si++)
    {
        F_y_B[num_species + 4 + si] = v_y_B*(Q_y_B[num_species + 4 + si][idx]);
    }
    
    for (int si = 0; si < num_species; si++)
    {
        F_y_T[si] = v_y_T*Q_y_T[si][idx];
    }
    F_y_T[num_species] = v_y_T*Q_y_T[num_species][idx];
    F_y_T[num_species + 1] = v_y_T*Q_y_T[num_species + 1][idx] + p_y_T[idx];
    F_y_T[num_species + 2] = v_y_T*Q_y_T[num_species + 2][idx];
    F_y_T[num_species + 3] = v_y_T*(Q_y_T[num_species + 3][idx] + p_y_T[idx]);
    for (int si = 0; si < num_species - 1; si++)
    {
        F_y_T[num_species + 4 + si] = v_y_T*(Q_y_T[num_species + 4 + si][idx]);
    }
    
    for (int si = 0; si < num_species; si++)
    {
        F_y_HLL[si] = (s_y_T*F_y_B[si] - s_y_B*F_y_T[si] + s_y_T*s_y_B*(Q_y_T[si][idx] - Q_y_B[si][idx]))/
            (s_y_T - s_y_B);
    }
    F_y_HLL[num_species] = (s_y_T*F_y_B[num_species] - s_y_B*F_y_T[num_species] + s_y_T*s_y_B*
        (Q_y_T[num_species][idx] - Q_y_B[num_species][idx]))/(s_y_T - s_y_B);
    F_y_HLL[num_species + 1] = (s_y_T*F_y_B[num_species + 2] - s_y_B*F_y_T[num_species + 2] + s_y_T*s_y_B*
        (Q_y_T[num_species + 2][idx] - Q_y_B[num_species + 2][idx]))/(s_y_T - s_y_B);
    for (int si = 0; si < num_species - 1; si++)
    {
        F_y_HLL[num_species + 2 + si] = (s_y_T*F_y_B[num_species + 4 + si] - s_y_B*F_y_T[num_species + 4 + si] +
            s_y_T*s_y_B*(Q_y_T[num_species + 4 + si][idx] - Q_y_B[num_species + 4 + si][idx]))/(s_y_T - s_y_B);
    }
    
    if (s_y_B > double(0))
    {
        for (int si = 0; si < num_species; si++)
        {
            F_y_HLL[si] = F_y_B[si];
        }
        F_y_HLL[num_species] = F_y_B[num_species];
        F_y_HLL[num_species + 1] = F_y_B[num_species + 2];
        for (int si = 0; si < num_species - 1; si++)
        {
            F_y_HLL[num_species + 2 + si] = F_y_B[num_species + 4 + si];
        }
    }
    
    if (s_y_T < double(0))
    {
        for (int si = 0; si < num_species; si++)
        {
            F_y_HLL[si] = F_y_T[si];
        }
        F_y_HLL[num_species] = F_y_T[num_species];
        F_y_HLL[num_species + 1] = F_y_T[num_species + 2];
        for (int si = 0; si < num_species - 1; si++)
        {
            F_y_HLL[num_species + 2 + si] = F_y_T[num_species + 4 + si];
        }
    }
    
    if (s_y_star > double(0))
    {
        Chi_y_star_BT = (s_y_B - v_y_B)/(s_y_B - s_y_star);
        
        for (int si = 0; si < num_species; si++)
        {
            Q_y_star_BT[si] = Chi_y_star_BT*Q_y_B[si][idx];
        }
        Q_y_star_BT[num_species] = Chi_y_star_BT*Q_y_B[num_species][idx];
        Q_y_star_BT[num_species + 1] = Chi_y_star_BT*rho_y_B[idx]*s_y_star;
        Q_y_star_BT[num_species + 2] = Chi_y_star_BT*Q_y_B[num_species + 2][idx];
        Q_y_star_BT[num_species + 3] = Chi_y_star_BT*(Q_y_B[num_species + 3][idx] +
            (s_y_star - v_y_B)*(rho_y_B[idx]*s_y_star + p_y_B[idx]/(s_y_B - v_y_B)));
        for (int si = 0; si < num_species - 1; si++)
        {
            Q_y_star_BT[num_species + 4 + si] = Chi_y_star_BT*Q_y_B[num_species + 4 + si][idx];
        }
        
        for (int ei = 0; ei < num_eqn; ei++)
        {
            F_y_HLLC[ei] = F_y_B[ei] + s_y_minus*(Q_y_star_BT[ei] - Q_y_B[ei][idx]);
        }
    }
    else
    {
        Chi_y_star_BT = (s_y_T - v_y_T)/(s_y_T - s_y_star);
        
        for (int si = 0; si < num_species; si++)
        {
            Q_y_star_BT[si] = Chi_y_star_BT*Q_y_T[si][idx];
        }
        Q_y_star_BT[num_species] = Chi_y_star_BT*Q_y_T[num_species][idx];
        Q_y_star_BT[num_species + 1] = Chi_y_star_BT*rho_y_T[idx]*s_y_star;
        Q_y_star_BT[num_species + 2] = Chi_y_star_BT*Q_y_T[num_species + 2][idx];
        Q_y_star_BT[num_species + 3] = Chi_y_star_BT*(Q_y_T[num_species + 3][idx] +
            (s_y_star - v_y_T)*(rho_y_T[idx]*s_y_star + p_y_T[idx]/(s_y_T - v_y_T)));
        for (int si = 0; si < num_species - 1; si++)
        {
            Q_y_star_BT[num_species + 4 + si] = Chi_y_star_BT*Q_y_T[num_species + 4 + si][idx];
        }
        
        for (int ei = 0; ei < num_eqn; ei++)
        {
            F_y_HLLC[ei] = F_y_T[ei] + s_y_plus*(Q_y_star_BT[ei] - Q_y_T[ei][idx]);
        }
    }
    
    /*
     * Calulate the weights beta for hybridization.
     */
    
    const double u_y_diff = u_y_T - u_y_B;
    const double v_y_diff = v_y_T - v_y_B;
    const double w_y_diff = w_y_T - w_y_B;
    const double vel_mag = sqrt(u_y_diff*u_y_diff + v_y_diff*v_y_diff + w_y_diff*w_y_diff);
    
    double alpha_1, alpha_2;
    if (vel_mag < EPSILON)
    {
        alpha_1 = double(1);
        alpha_2 = double(0);
    }
    else
    {
        alpha_1 = fabs(v_y_diff)/vel_mag;
        alpha_2 = sqrt(double(1) - alpha_1*alpha_1);
    }
    
    const double beta_1 = double(1)/double(2)*(double(1) + alpha_1/(alpha_1 + alpha_2));
    const double beta_2 = double(1) - beta_1;
    
    for (int si = 0; si < num_species; si++)
    {
        F_y[si][idx_flux] = beta_1*F_y_HLLC[si] + beta_2*F_y_HLL[si];
    }
    F_y[num_species][idx_flux] = beta_1*F_y_HLLC[num_species] + beta_2*F_y_HLL[num_species];
    F_y[num_species + 1][idx_flux] = F_y_HLLC[num_species + 1];
    F_y[num_species + 2][idx_flux] = beta_1*F_y_HLLC[num_species + 2] + beta_2*F_y_HLL[num_species + 1];
    F_y[num_species + 3][idx_flux] = F_y_HLLC[num_species + 3];
    for (int si = 0; si < num_species - 1; si++)
    {
        F_y[num_species + 4 + si][idx_flux] = beta_1*F_y_HLLC[num_species + 4 + si] +
            beta_2*F_y_HLL[num_species + 2 + si];
    }
}


/*
 * Compute the local convective flux in the z-direction from conservative variables with
 * 3D HLLC-HLL Riemann solver.
 */
static inline __attribute__((always_inline)) void
computeLocalConvectiveFluxInZDirectionFromConservativeVariablesHLLC_HLL3D(
    double** F_z,
    double** Q_z_B,
    double** Q_z_F,
    double* rho_z_B,
    double* rho_z_F,
    double* p_z_B,
    double* p_z_F,
    double* c_z_B,
    double* c_z_F,
    double& w_z_B,
    double& w_z_F,
    double& s_z_minus,
    double& s_z_plus,
    double& s_z_star,
    double& Chi_z_star_BF,
    const int& idx_flux,
    const int& idx,
    const int& num_species,
    const int& num_eqn)
{
    w_z_B = Q_z_B[num_species + 2][idx]/rho_z_B[idx];
    w_z_F = Q_z_F[num_species + 2][idx]/rho_z_F[idx];
    
    const double u_z_B = Q_z_B[num_species][idx]/rho_z_B[idx];
    const double u_z_F = Q_z_F[num_species][idx]/rho_z_F[idx];
    
    const double v_z_B = Q_z_B[num_species + 1][idx]/rho_z_B[idx];
    const double v_z_F = Q_z_F[num_species + 1][idx]/rho_z_F[idx];
   
    const double w_z_average = double(1)/double(2)*(w_z_B + w_z_F);
    const double c_z_average = double(1)/double(2)*(c_z_B[idx] + c_z_F[idx]);
    
    const double s_z_B = fmin(w_z_average - c_z_average, w_z_B - c_z_B[idx]);
    const double s_z_F = fmax(w_z_average + c_z_average, w_z_F + c_z_F[idx]);
    
    s_z_minus = fmin(double(0), s_z_B);
    s_z_plus  = fmax(double(0), s_z_F);
    
    s_z_star = (p_z_F[idx] - p_z_B[idx] +
        Q_z_B[num_species + 2][idx]*(s_z_B - w_z_B) - Q_z_F[num_species + 2][idx]*(s_z_F - w_z_F))/
        (rho_z_B[idx]*(s_z_B - w_z_B) - rho_z_F[idx]*(s_z_F - w_z_F));
    
    double F_z_B[num_eqn];
    double F_z_F[num_eqn];
    double F_z_HLL[2*num_species + 1];
    double F_z_HLLC[num_eqn];
    double Q_z_star_BF[num_eqn];
    
    for (int si = 0; si < num_species; si++)
    {
        F_z_B[si] = w_z_B*Q_z_B[si][idx];
    }
    F_z_B[num_species] = w_z_B*Q_z_B[num_species][idx];
    F_z_B[num_species + 1] = w_z_B*Q_z_B[num_species + 1][idx];
    F_z_B[num_species + 2] = w_z_B*Q_z_B[num_species + 2][idx] + p_z_B[idx];
    F_z_B[num_species + 3] = w_z_B*(Q_z_B[num_species + 3][idx] + p_z_B[idx]);
    for (int si = 0; si < num_species - 1; si++)
    {
        F_z_B[num_species + 4 + si] = w_z_B*(Q_z_B[num_species + 4 + si][idx]);
    }
    
    for (int si = 0; si < num_species; si++)
    {
        F_z_F[si] = w_z_F*Q_z_F[si][idx];
    }
    F_z_F[num_species] = w_z_F*Q_z_F[num_species][idx];
    F_z_F[num_species + 1] = w_z_F*Q_z_F[num_species + 1][idx];
    F_z_F[num_species + 2] = w_z_F*Q_z_F[num_species + 2][idx] + p_z_F[idx];
    F_z_F[num_species + 3] = w_z_F*(Q_z_F[num_species + 3][idx] + p_z_F[idx]);
    for (int si = 0; si < num_species - 1; si++)
    {
        F_z_F[num_species + 4 + si] = w_z_F*(Q_z_F[num_species + 4 + si][idx]);
    }
    
    for (int si = 0; si < num_species; si++)
    {
        F_z_HLL[si] = (s_z_F*F_z_B[si] - s_z_B*F_z_F[si] + s_z_F*s_z_B*(Q_z_F[si][idx] - Q_z_B[si][idx]))/
            (s_z_F - s_z_B);
    }
    F_z_HLL[num_species] = (s_z_F*F_z_B[num_species] - s_z_B*F_z_F[num_species] + s_z_F*s_z_B*
        (Q_z_F[num_species][idx] - Q_z_B[num_species][idx]))/(s_z_F - s_z_B);
    F_z_HLL[num_species + 1] = (s_z_F*F_z_B[num_species + 1] - s_z_B*F_z_F[num_species + 1] + s_z_F*s_z_B*
        (Q_z_F[num_species + 1][idx] - Q_z_B[num_species + 1][idx]))/(s_z_F - s_z_B);
    for (int si = 0; si < num_species - 1; si++)
    {
        F_z_HLL[num_species + 2 + si] = (s_z_F*F_z_B[num_species + 4 + si] - s_z_B*F_z_F[num_species + 4 + si] +
            s_z_F*s_z_B*(Q_z_F[num_species + 4 + si][idx] - Q_z_B[num_species + 4 + si][idx]))/(s_z_F - s_z_B);
    }
    
    if (s_z_B > double(0))
    {
        for (int si = 0; si < num_species; si++)
        {
            F_z_HLL[si] = F_z_B[si];
        }
        F_z_HLL[num_species] = F_z_B[num_species];
        F_z_HLL[num_species + 1] = F_z_B[num_species + 1];
        for (int si = 0; si < num_species - 1; si++)
        {
            F_z_HLL[num_species + 2 + si] = F_z_B[num_species + 4 + si];
        }
    }
    
    if (s_z_F < double(0))
    {
        for (int si = 0; si < num_species; si++)
        {
            F_z_HLL[si] = F_z_F[si];
        }
        F_z_HLL[num_species] = F_z_F[num_species];
        F_z_HLL[num_species + 1] = F_z_F[num_species + 1];
        for (int si = 0; si < num_species - 1; si++)
        {
            F_z_HLL[num_species + 2 + si] = F_z_F[num_species + 4 + si];
        }
    }
    
    if (s_z_star > double(0))
    {
        Chi_z_star_BF = (s_z_B - w_z_B)/(s_z_B - s_z_star);
        
        for (int si = 0; si < num_species; si++)
        {
            Q_z_star_BF[si] = Chi_z_star_BF*Q_z_B[si][idx];
        }
        Q_z_star_BF[num_species] = Chi_z_star_BF*Q_z_B[num_species][idx];
        Q_z_star_BF[num_species + 1] = Chi_z_star_BF*Q_z_B[num_species + 1][idx];
        Q_z_star_BF[num_species + 2] = Chi_z_star_BF*rho_z_B[idx]*s_z_star;
        Q_z_star_BF[num_species + 3] = Chi_z_star_BF*(Q_z_B[num_species + 3][idx] +
            (s_z_star - w_z_B)*(rho_z_B[idx]*s_z_star + p_z_B[idx]/(s_z_B - w_z_B)));
        for (int si = 0; si < num_species - 1; si++)
        {
            Q_z_star_BF[num_species + 4 + si] = Chi_z_star_BF*Q_z_B[num_species + 4 + si][idx];
        }
        
        for (int ei = 0; ei < num_eqn; ei++)
        {
            F_z_HLLC[ei] = F_z_B[ei] + s_z_minus*(Q_z_star_BF[ei] - Q_z_B[ei][idx]);
        }
    }
    else
    {
        Chi_z_star_BF = (s_z_F - w_z_F)/(s_z_F - s_z_star);
        
        for (int si = 0; si < num_species; si++)
        {
            Q_z_star_BF[si] = Chi_z_star_BF*Q_z_F[si][idx];
        }
        Q_z_star_BF[num_species] = Chi_z_star_BF*Q_z_F[num_species][idx];
        Q_z_star_BF[num_species + 1] = Chi_z_star_BF*Q_z_F[num_species + 1][idx];
        Q_z_star_BF[num_species + 2] = Chi_z_star_BF*rho_z_F[idx]*s_z_star;
        Q_z_star_BF[num_species + 3] = Chi_z_star_BF*(Q_z_F[num_species + 3][idx] +
            (s_z_star - w_z_F)*(rho_z_F[idx]*s_z_star + p_z_F[idx]/(s_z_F - w_z_F)));
        for (int si = 0; si < num_species - 1; si++)
        {
            Q_z_star_BF[num_species + 4 + si] = Chi_z_star_BF*Q_z_F[num_species + 4 + si][idx];
        }
        
        for (int ei = 0; ei < num_eqn; ei++)
        {
            F_z_HLLC[ei] = F_z_F[ei] + s_z_plus*(Q_z_star_BF[ei] - Q_z_F[ei][idx]);
        }
    }
    
    /*
     * Calulate the weights beta for hybridization.
     */
    
    const double u_z_diff = u_z_F - u_z_B;
    const double v_z_diff = v_z_F - v_z_B;
    const double w_z_diff = w_z_F - w_z_B;
    const double vel_mag = sqrt(u_z_diff*u_z_diff + v_z_diff*v_z_diff + w_z_diff*w_z_diff);
    
    double alpha_1, alpha_2;
    if (vel_mag < EPSILON)
    {
        alpha_1 = double(1);
        alpha_2 = double(0);
    }
    else
    {
        alpha_1 = fabs(w_z_diff)/vel_mag;
        alpha_2 = sqrt(double(1) - alpha_1*alpha_1);
    }
    
    const double beta_1 = double(1)/double(2)*(double(1) + alpha_1/(alpha_1 + alpha_2));
    const double beta_2 = double(1) - beta_1;
    
    for (int si = 0; si < num_species; si++)
    {
        F_z[si][idx_flux] = beta_1*F_z_HLLC[si] + beta_2*F_z_HLL[si];
    }
    F_z[num_species][idx_flux] = beta_1*F_z_HLLC[num_species] + beta_2*F_z_HLL[num_species];
    F_z[num_species + 1][idx_flux] = beta_1*F_z_HLLC[num_species + 1] + beta_2*F_z_HLL[num_species + 1];
    F_z[num_species + 2][idx_flux] = F_z_HLLC[num_species + 2];
    F_z[num_species + 3][idx_flux] = F_z_HLLC[num_species + 3];
    for (int si = 0; si < num_species - 1; si++)
    {
        F_z[num_species + 4 + si][idx_flux] = beta_1*F_z_HLLC[num_species + 4 + si] +
            beta_2*F_z_HLL[num_species + 2 + si];
    }
}


/*
 * Compute the local convective flux in the x-direction from primitive variables with
 * 1D HLLC-HLL Riemann solver.
 */
static inline __attribute__((always_inline)) void
computeLocalConvectiveFluxInXDirectionFromPrimitiveVariablesHLLC_HLL1D(
    double** F_x,
    double** V_x_L,
    double** V_x_R,
    double* rho_x_L,
    double* rho_x_R,
    double* c_x_L,
    double* c_x_R,
    double* epsilon_x_L,
    double* epsilon_x_R,
    double& s_x_minus,
    double& s_x_plus,
    double& s_x_star,
    double& Chi_x_star_LR,
    const int& idx_flux,
    const int& idx,
    const int& num_species,
    const int& num_eqn)
{
    const double u_x_average = double(1)/double(2)*(V_x_L[num_species][idx] + V_x_R[num_species][idx]);
    const double c_x_average = double(1)/double(2)*(c_x_L[idx] + c_x_R[idx]);
    
    const double s_x_L = fmin(u_x_average - c_x_average, V_x_L[num_species][idx] - c_x_L[idx]);
    const double s_x_R = fmax(u_x_average + c_x_average, V_x_R[num_species][idx] + c_x_R[idx]);
    
    s_x_minus = fmin(double(0), s_x_L);
    s_x_plus  = fmax(double(0), s_x_R);
    
    s_x_star = (V_x_R[num_species + 1][idx] - V_x_L[num_species + 1][idx] +
        rho_x_L[idx]*V_x_L[num_species][idx]*(s_x_L - V_x_L[num_species][idx]) -
        rho_x_R[idx]*V_x_R[num_species][idx]*(s_x_R - V_x_R[num_species][idx]))/
        (rho_x_L[idx]*(s_x_L - V_x_L[num_species][idx]) -
        rho_x_R[idx]*(s_x_R - V_x_R[num_species][idx]));
    
    double Q_x_LR[num_eqn];
    double Q_x_star_LR[num_eqn];
    double F_x_LR[num_eqn];
    
    if (s_x_star > double(0))
    {
        Chi_x_star_LR = (s_x_L - V_x_L[num_species][idx])/(s_x_L - s_x_star);
        
        for (int si = 0; si < num_species; si++)
        {
            Q_x_LR[si] = V_x_L[si][idx];
        }
        Q_x_LR[num_species] = rho_x_L[idx]*V_x_L[num_species][idx];
        Q_x_LR[num_species + 1] = rho_x_L[idx]*(epsilon_x_L[idx] +
            double(1)/double(2)*V_x_L[num_species][idx]*V_x_L[num_species][idx]);
        for (int si = 0; si < num_species - 1; si++)
        {
            Q_x_LR[num_species + 2 + si] = V_x_L[num_species + 2 + si][idx];
        }
        
        for (int si = 0; si < num_species; si++)
        {
            Q_x_star_LR[si] = Chi_x_star_LR*V_x_L[si][idx];
        }
        Q_x_star_LR[num_species] = Chi_x_star_LR*rho_x_L[idx]*s_x_star;
        Q_x_star_LR[num_species + 1] = Chi_x_star_LR*(Q_x_LR[num_species + 1] +
            (s_x_star - V_x_L[num_species][idx])*(rho_x_L[idx]*s_x_star + V_x_L[num_species + 1][idx]/
            (s_x_L - V_x_L[num_species][idx])));
        for (int si = 0; si < num_species - 1; si++)
        {
            Q_x_star_LR[num_species + 2 + si] = Chi_x_star_LR*V_x_L[num_species + 2 + si][idx];
        }
        
        for (int si = 0; si < num_species; si++)
        {
            F_x_LR[si] = V_x_L[num_species][idx]*V_x_L[si][idx];
        }
        F_x_LR[num_species] = V_x_L[num_species][idx]*Q_x_LR[num_species] + V_x_L[num_species + 1][idx];
        F_x_LR[num_species + 1] = V_x_L[num_species][idx]*(Q_x_LR[num_species + 1] + V_x_L[num_species + 1][idx]);
        for (int si = 0; si < num_species - 1; si++)
        {
            F_x_LR[num_species + 2 + si] = V_x_L[num_species][idx]*V_x_L[num_species + 2 + si][idx];
        }
        
        for (int ei = 0; ei < num_eqn; ei++)
        {
            F_x[ei][idx_flux] = F_x_LR[ei] + s_x_minus*(Q_x_star_LR[ei] - Q_x_LR[ei]);
        }
    }
    else
    {
        Chi_x_star_LR = (s_x_R - V_x_R[num_species][idx])/(s_x_R - s_x_star);
        
        for (int si = 0; si < num_species; si++)
        {
            Q_x_LR[si] = V_x_R[si][idx];
        }
        Q_x_LR[num_species] = rho_x_R[idx]*V_x_R[num_species][idx];
        Q_x_LR[num_species + 1] = rho_x_R[idx]*(epsilon_x_R[idx] +
            double(1)/double(2)*V_x_R[num_species][idx]*V_x_R[num_species][idx]);
        for (int si = 0; si < num_species - 1; si++)
        {
            Q_x_LR[num_species + 2 + si] = V_x_R[num_species + 2 + si][idx];
        }
        
        for (int si = 0; si < num_species; si++)
        {
            Q_x_star_LR[si] = Chi_x_star_LR*V_x_R[si][idx];
        }
        Q_x_star_LR[num_species] = Chi_x_star_LR*rho_x_R[idx]*s_x_star;
        Q_x_star_LR[num_species + 1] = Chi_x_star_LR*(Q_x_LR[num_species + 1] +
            (s_x_star - V_x_R[num_species][idx])*(rho_x_R[idx]*s_x_star + V_x_R[num_species + 1][idx]/
            (s_x_R - V_x_R[num_species][idx])));
        for (int si = 0; si < num_species - 1; si++)
        {
            Q_x_star_LR[num_species + 2 + si] = Chi_x_star_LR*V_x_R[num_species + 2 + si][idx];
        }
        
        for (int si = 0; si < num_species; si++)
        {
            F_x_LR[si] = V_x_R[num_species][idx]*V_x_R[si][idx];
        }
        F_x_LR[num_species] = V_x_R[num_species][idx]*Q_x_LR[num_species] + V_x_R[num_species + 1][idx];
        F_x_LR[num_species + 1] = V_x_R[num_species][idx]*(Q_x_LR[num_species + 1] + V_x_R[num_species + 1][idx]);
        for (int si = 0; si < num_species - 1; si++)
        {
            F_x_LR[num_species + 2 + si] = V_x_R[num_species][idx]*V_x_R[num_species + 2 + si][idx];
        }
        
        for (int ei = 0; ei < num_eqn; ei++)
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
    double** F_x,
    double** V_x_L,
    double** V_x_R,
    double* rho_x_L,
    double* rho_x_R,
    double* c_x_L,
    double* c_x_R,
    double* epsilon_x_L,
    double* epsilon_x_R,
    double& s_x_minus,
    double& s_x_plus,
    double& s_x_star,
    double& Chi_x_star_LR,
    const int& idx_flux,
    const int& idx,
    const int& num_species,
    const int& num_eqn)
{
    const double u_x_average = double(1)/double(2)*(V_x_L[num_species][idx] + V_x_R[num_species][idx]);
    const double c_x_average = double(1)/double(2)*(c_x_L[idx] + c_x_R[idx]);
    
    const double s_x_L = fmin(u_x_average - c_x_average, V_x_L[num_species][idx] - c_x_L[idx]);
    const double s_x_R = fmax(u_x_average + c_x_average, V_x_R[num_species][idx] + c_x_R[idx]);
    
    s_x_minus = fmin(double(0), s_x_L);
    s_x_plus  = fmax(double(0), s_x_R);
    
    s_x_star = (V_x_R[num_species + 2][idx] - V_x_L[num_species + 2][idx] +
        rho_x_L[idx]*V_x_L[num_species][idx]*(s_x_L - V_x_L[num_species][idx]) -
        rho_x_R[idx]*V_x_R[num_species][idx]*(s_x_R - V_x_R[num_species][idx]))/
        (rho_x_L[idx]*(s_x_L - V_x_L[num_species][idx]) -
        rho_x_R[idx]*(s_x_R - V_x_R[num_species][idx]));
    
    double Q_x_L[num_eqn];
    double Q_x_R[num_eqn];
    double F_x_L[num_eqn];
    double F_x_R[num_eqn];
    double F_x_HLL[2*num_species];
    double F_x_HLLC[num_eqn];
    double Q_x_star_LR[num_eqn];
    
    for (int si = 0; si < num_species; si++)
    {
        Q_x_L[si] = V_x_L[si][idx];
    }
    Q_x_L[num_species] = rho_x_L[idx]*V_x_L[num_species][idx];
    Q_x_L[num_species + 1] = rho_x_L[idx]*V_x_L[num_species + 1][idx];
    Q_x_L[num_species + 2] = rho_x_L[idx]*(epsilon_x_L[idx] +
        double(1)/double(2)*(V_x_L[num_species][idx]*V_x_L[num_species][idx] +
        V_x_L[num_species + 1][idx]*V_x_L[num_species + 1][idx]));
    for (int si = 0; si < num_species - 1; si++)
    {
        Q_x_L[num_species + 3 + si] = V_x_L[num_species + 3 + si][idx];
    }
    
    for (int si = 0; si < num_species; si++)
    {
        Q_x_R[si] = V_x_R[si][idx];
    }
    Q_x_R[num_species] = rho_x_R[idx]*V_x_R[num_species][idx];
    Q_x_R[num_species + 1] = rho_x_R[idx]*V_x_R[num_species + 1][idx];
    Q_x_R[num_species + 2] = rho_x_R[idx]*(epsilon_x_R[idx] +
        double(1)/double(2)*(V_x_R[num_species][idx]*V_x_R[num_species][idx] +
        V_x_R[num_species + 1][idx]*V_x_R[num_species + 1][idx]));
    for (int si = 0; si < num_species - 1; si++)
    {
        Q_x_R[num_species + 3 + si] = V_x_R[num_species + 3 + si][idx];
    }
    
    for (int si = 0; si < num_species; si++)
    {
        F_x_L[si] = V_x_L[num_species][idx]*V_x_L[si][idx];
    }
    F_x_L[num_species] = V_x_L[num_species][idx]*Q_x_L[num_species] + V_x_L[num_species + 2][idx];
    F_x_L[num_species + 1] = V_x_L[num_species][idx]*Q_x_L[num_species + 1];
    F_x_L[num_species + 2] = V_x_L[num_species][idx]*(Q_x_L[num_species + 2] + V_x_L[num_species + 2][idx]);
    for (int si = 0; si < num_species - 1; si++)
    {
        F_x_L[num_species + 3 + si] = V_x_L[num_species][idx]*V_x_L[num_species + 3 + si][idx];
    }
    
    for (int si = 0; si < num_species; si++)
    {
        F_x_R[si] = V_x_R[num_species][idx]*V_x_R[si][idx];
    }
    F_x_R[num_species] = V_x_R[num_species][idx]*Q_x_R[num_species] + V_x_R[num_species + 2][idx];
    F_x_R[num_species + 1] = V_x_R[num_species][idx]*Q_x_R[num_species + 1];
    F_x_R[num_species + 2] = V_x_R[num_species][idx]*(Q_x_R[num_species + 2] + V_x_R[num_species + 2][idx]);
    for (int si = 0; si < num_species - 1; si++)
    {
        F_x_R[num_species + 3 + si] = V_x_R[num_species][idx]*V_x_R[num_species + 3 + si][idx];
    }
    
    for (int si = 0; si < num_species; si++)
    {
        F_x_HLL[si] = (s_x_R*F_x_L[si] - s_x_L*F_x_R[si] + s_x_R*s_x_L*(Q_x_R[si] - Q_x_L[si]))/
            (s_x_R - s_x_L);
    }
    F_x_HLL[num_species] = (s_x_R*F_x_L[num_species + 1] - s_x_L*F_x_R[num_species + 1] + s_x_R*s_x_L*
        (Q_x_R[num_species + 1] - Q_x_L[num_species + 1]))/(s_x_R - s_x_L);
    for (int si = 0; si < num_species - 1; si++)
    {
        F_x_HLL[num_species + 1 + si] = (s_x_R*F_x_L[num_species + 3 + si] - s_x_L*F_x_R[num_species + 3 + si] +
            s_x_R*s_x_L*(Q_x_R[num_species + 3 + si] - Q_x_L[num_species + 3 + si]))/(s_x_R - s_x_L);
    }
    
    if (s_x_L > double(0))
    {
        for (int si = 0; si < num_species; si++)
        {
            F_x_HLL[si] = F_x_L[si];
        }
        F_x_HLL[num_species] = F_x_L[num_species + 1];
        for (int si = 0; si < num_species - 1; si++)
        {
            F_x_HLL[num_species + 1 + si] = F_x_L[num_species + 3 + si];
        }
    }
    
    if (s_x_R < double(0))
    {
        for (int si = 0; si < num_species; si++)
        {
            F_x_HLL[si] = F_x_R[si];
        }
        F_x_HLL[num_species] = F_x_R[num_species + 1];
        for (int si = 0; si < num_species - 1; si++)
        {
            F_x_HLL[num_species + 1 + si] = F_x_R[num_species + 3 + si];
        }
    }
    
    if (s_x_star > double(0))
    {
        Chi_x_star_LR = (s_x_L - V_x_L[num_species][idx])/(s_x_L - s_x_star);
        
        for (int si = 0; si < num_species; si++)
        {
            Q_x_star_LR[si] = Chi_x_star_LR*V_x_L[si][idx];
        }
        Q_x_star_LR[num_species] = Chi_x_star_LR*rho_x_L[idx]*s_x_star;
        Q_x_star_LR[num_species + 1] = Chi_x_star_LR*Q_x_L[num_species + 1];
        Q_x_star_LR[num_species + 2] = Chi_x_star_LR*(Q_x_L[num_species + 2] +
            (s_x_star - V_x_L[num_species][idx])*(rho_x_L[idx]*s_x_star + V_x_L[num_species + 2][idx]/
            (s_x_L - V_x_L[num_species][idx])));
        for (int si = 0; si < num_species - 1; si++)
        {
            Q_x_star_LR[num_species + 3 + si] = Chi_x_star_LR*V_x_L[num_species + 3 + si][idx];
        }
        
        for (int ei = 0; ei < num_eqn; ei++)
        {
            F_x_HLLC[ei] = F_x_L[ei] + s_x_minus*(Q_x_star_LR[ei] - Q_x_L[ei]);
        }
    }
    else
    {
        Chi_x_star_LR = (s_x_R - V_x_R[num_species][idx])/(s_x_R - s_x_star);
        
        for (int si = 0; si < num_species; si++)
        {
            Q_x_star_LR[si] = Chi_x_star_LR*V_x_R[si][idx];
        }
        Q_x_star_LR[num_species] = Chi_x_star_LR*rho_x_R[idx]*s_x_star;
        Q_x_star_LR[num_species + 1] = Chi_x_star_LR*Q_x_R[num_species + 1];
        Q_x_star_LR[num_species + 2] = Chi_x_star_LR*(Q_x_R[num_species + 2] +
            (s_x_star - V_x_R[num_species][idx])*(rho_x_R[idx]*s_x_star + V_x_R[num_species + 2][idx]/
            (s_x_R - V_x_R[num_species][idx])));
        for (int si = 0; si < num_species - 1; si++)
        {
            Q_x_star_LR[num_species + 3 + si] = Chi_x_star_LR*V_x_R[num_species + 3 + si][idx];
        }
        
        for (int ei = 0; ei < num_eqn; ei++)
        {
            F_x_HLLC[ei] = F_x_R[ei] + s_x_plus*(Q_x_star_LR[ei] - Q_x_R[ei]);
        }
    }
    
    /*
     * Calulate the weights beta for hybridization.
     */
    
    const double u_x_diff = V_x_R[num_species][idx] - V_x_L[num_species][idx];
    const double v_x_diff = V_x_R[num_species + 1][idx] - V_x_L[num_species + 1][idx];
    const double vel_mag = sqrt(u_x_diff*u_x_diff + v_x_diff*v_x_diff);
    
    double alpha_1, alpha_2;
    if (vel_mag < EPSILON)
    {
        alpha_1 = double(1);
        alpha_2 = double(0);
    }
    else
    {
        alpha_1 = fabs(u_x_diff)/vel_mag;
        alpha_2 = sqrt(double(1) - alpha_1*alpha_1);
    }
    
    const double beta_1 = double(1)/double(2)*(double(1) + alpha_1/(alpha_1 + alpha_2));
    const double beta_2 = double(1) - beta_1;
    
    for (int si = 0; si < num_species; si++)
    {
        F_x[si][idx_flux] = beta_1*F_x_HLLC[si] + beta_2*F_x_HLL[si];
    }
    F_x[num_species][idx_flux] = F_x_HLLC[num_species];
    F_x[num_species + 1][idx_flux] = beta_1*F_x_HLLC[num_species + 1] + beta_2*F_x_HLL[num_species];
    F_x[num_species + 2][idx_flux] = F_x_HLLC[num_species + 2];
    for (int si = 0; si < num_species - 1; si++)
    {
        F_x[num_species + 3 + si][idx_flux] = beta_1*F_x_HLLC[num_species + 3 + si] +
            beta_2*F_x_HLL[num_species + 1 + si];
    }
}


/*
 * Compute the local convective flux in the x-direction from primitive variables with
 * 3D HLLC-HLL Riemann solver.
 */
static inline __attribute__((always_inline)) void
computeLocalConvectiveFluxInXDirectionFromPrimitiveVariablesHLLC_HLL3D(
    double** F_x,
    double** V_x_L,
    double** V_x_R,
    double* rho_x_L,
    double* rho_x_R,
    double* c_x_L,
    double* c_x_R,
    double* epsilon_x_L,
    double* epsilon_x_R,
    double& s_x_minus,
    double& s_x_plus,
    double& s_x_star,
    double& Chi_x_star_LR,
    const int& idx_flux,
    const int& idx,
    const int& num_species,
    const int& num_eqn)
{
    const double u_x_average = double(1)/double(2)*(V_x_L[num_species][idx] + V_x_R[num_species][idx]);
    const double c_x_average = double(1)/double(2)*(c_x_L[idx] + c_x_R[idx]);
    
    const double s_x_L = fmin(u_x_average - c_x_average, V_x_L[num_species][idx] - c_x_L[idx]);
    const double s_x_R = fmax(u_x_average + c_x_average, V_x_R[num_species][idx] + c_x_R[idx]);
    
    s_x_minus = fmin(double(0), s_x_L);
    s_x_plus  = fmax(double(0), s_x_R);
    
    s_x_star = (V_x_R[num_species + 3][idx] - V_x_L[num_species + 3][idx] +
        rho_x_L[idx]*V_x_L[num_species][idx]*(s_x_L - V_x_L[num_species][idx]) -
        rho_x_R[idx]*V_x_R[num_species][idx]*(s_x_R - V_x_R[num_species][idx]))/
        (rho_x_L[idx]*(s_x_L - V_x_L[num_species][idx]) -
        rho_x_R[idx]*(s_x_R - V_x_R[num_species][idx]));
    
    double Q_x_L[num_eqn];
    double Q_x_R[num_eqn];
    double F_x_L[num_eqn];
    double F_x_R[num_eqn];
    double F_x_HLL[2*num_species + 1];
    double F_x_HLLC[num_eqn];
    double Q_x_star_LR[num_eqn];
    
    for (int si = 0; si < num_species; si++)
    {
        Q_x_L[si] = V_x_L[si][idx];
    }
    Q_x_L[num_species] = rho_x_L[idx]*V_x_L[num_species][idx];
    Q_x_L[num_species + 1] = rho_x_L[idx]*V_x_L[num_species + 1][idx];
    Q_x_L[num_species + 2] = rho_x_L[idx]*V_x_L[num_species + 2][idx];
    Q_x_L[num_species + 3] = rho_x_L[idx]*(epsilon_x_L[idx] +
        double(1)/double(2)*(V_x_L[num_species][idx]*V_x_L[num_species][idx] +
        V_x_L[num_species + 1][idx]*V_x_L[num_species + 1][idx] +
        V_x_L[num_species + 2][idx]*V_x_L[num_species + 2][idx]));
    for (int si = 0; si < num_species - 1; si++)
    {
        Q_x_L[num_species + 4 + si] = V_x_L[num_species + 4 + si][idx];
    }
    
    for (int si = 0; si < num_species; si++)
    {
        Q_x_R[si] = V_x_R[si][idx];
    }
    Q_x_R[num_species] = rho_x_R[idx]*V_x_R[num_species][idx];
    Q_x_R[num_species + 1] = rho_x_R[idx]*V_x_R[num_species + 1][idx];
    Q_x_R[num_species + 2] = rho_x_R[idx]*V_x_R[num_species + 2][idx];
    Q_x_R[num_species + 3] = rho_x_R[idx]*(epsilon_x_R[idx] +
        double(1)/double(2)*(V_x_R[num_species][idx]*V_x_R[num_species][idx] +
        V_x_R[num_species + 1][idx]*V_x_R[num_species + 1][idx] +
        V_x_R[num_species + 2][idx]*V_x_R[num_species + 2][idx]));
    for (int si = 0; si < num_species - 1; si++)
    {
        Q_x_R[num_species + 4 + si] = V_x_R[num_species + 4 + si][idx];
    }
    
    for (int si = 0; si < num_species; si++)
    {
        F_x_L[si] = V_x_L[num_species][idx]*V_x_L[si][idx];
    }
    F_x_L[num_species] = V_x_L[num_species][idx]*Q_x_L[num_species] + V_x_L[num_species + 3][idx];
    F_x_L[num_species + 1] = V_x_L[num_species][idx]*Q_x_L[num_species + 1];
    F_x_L[num_species + 2] = V_x_L[num_species][idx]*Q_x_L[num_species + 2];
    F_x_L[num_species + 3] = V_x_L[num_species][idx]*(Q_x_L[num_species + 3] + V_x_L[num_species + 3][idx]);
    for (int si = 0; si < num_species - 1; si++)
    {
        F_x_L[num_species + 4 + si] = V_x_L[num_species][idx]*V_x_L[num_species + 4 + si][idx];
    }
    
    for (int si = 0; si < num_species; si++)
    {
        F_x_R[si] = V_x_R[num_species][idx]*V_x_R[si][idx];
    }
    F_x_R[num_species] = V_x_R[num_species][idx]*Q_x_R[num_species] + V_x_R[num_species + 3][idx];
    F_x_R[num_species + 1] = V_x_R[num_species][idx]*Q_x_R[num_species + 1];
    F_x_R[num_species + 2] = V_x_R[num_species][idx]*Q_x_R[num_species + 2];
    F_x_R[num_species + 3] = V_x_R[num_species][idx]*(Q_x_R[num_species + 3] + V_x_R[num_species + 3][idx]);
    for (int si = 0; si < num_species - 1; si++)
    {
        F_x_R[num_species + 4 + si] = V_x_R[num_species][idx]*V_x_R[num_species + 4 + si][idx];
    }
    
    for (int si = 0; si < num_species; si++)
    {
        F_x_HLL[si] = (s_x_R*F_x_L[si] - s_x_L*F_x_R[si] + s_x_R*s_x_L*(Q_x_R[si] - Q_x_L[si]))/
            (s_x_R - s_x_L);
    }
    F_x_HLL[num_species] = (s_x_R*F_x_L[num_species + 1] - s_x_L*F_x_R[num_species + 1] + s_x_R*s_x_L*
        (Q_x_R[num_species + 1] - Q_x_L[num_species + 1]))/(s_x_R - s_x_L);
    F_x_HLL[num_species + 1] = (s_x_R*F_x_L[num_species + 2] - s_x_L*F_x_R[num_species + 2] + s_x_R*s_x_L*
        (Q_x_R[num_species + 2] - Q_x_L[num_species + 2]))/(s_x_R - s_x_L);
    for (int si = 0; si < num_species - 1; si++)
    {
        F_x_HLL[num_species + 2 + si] = (s_x_R*F_x_L[num_species + 4 + si] - s_x_L*F_x_R[num_species + 4 + si] +
            s_x_R*s_x_L*(Q_x_R[num_species + 4 + si] - Q_x_L[num_species + 4 + si]))/(s_x_R - s_x_L);
    }
    
    if (s_x_L > double(0))
    {
        for (int si = 0; si < num_species; si++)
        {
            F_x_HLL[si] = F_x_L[si];
        }
        F_x_HLL[num_species] = F_x_L[num_species + 1];
        F_x_HLL[num_species + 1] = F_x_L[num_species + 2];
        for (int si = 0; si < num_species - 1; si++)
        {
            F_x_HLL[num_species + 2 + si] = F_x_L[num_species + 4 + si];
        }
    }
    
    if (s_x_R < double(0))
    {
        for (int si = 0; si < num_species; si++)
        {
            F_x_HLL[si] = F_x_R[si];
        }
        F_x_HLL[num_species] = F_x_R[num_species + 1];
        F_x_HLL[num_species + 1] = F_x_R[num_species + 2];
        for (int si = 0; si < num_species - 1; si++)
        {
            F_x_HLL[num_species + 2 + si] = F_x_R[num_species + 4 + si];
        }
    }
    
    if (s_x_star > double(0))
    {
        Chi_x_star_LR = (s_x_L - V_x_L[num_species][idx])/(s_x_L - s_x_star);
        
        for (int si = 0; si < num_species; si++)
        {
            Q_x_star_LR[si] = Chi_x_star_LR*V_x_L[si][idx];
        }
        Q_x_star_LR[num_species] = Chi_x_star_LR*rho_x_L[idx]*s_x_star;
        Q_x_star_LR[num_species + 1] = Chi_x_star_LR*Q_x_L[num_species + 1];
        Q_x_star_LR[num_species + 2] = Chi_x_star_LR*Q_x_L[num_species + 2];
        Q_x_star_LR[num_species + 3] = Chi_x_star_LR*(Q_x_L[num_species + 3] +
            (s_x_star - V_x_L[num_species][idx])*(rho_x_L[idx]*s_x_star + V_x_L[num_species + 3][idx]/
            (s_x_L - V_x_L[num_species][idx])));
        for (int si = 0; si < num_species - 1; si++)
        {
            Q_x_star_LR[num_species + 4 + si] = Chi_x_star_LR*V_x_L[num_species + 4 + si][idx];
        }
        
        for (int ei = 0; ei < num_eqn; ei++)
        {
            F_x_HLLC[ei] = F_x_L[ei] + s_x_minus*(Q_x_star_LR[ei] - Q_x_L[ei]);
        }
    }
    else
    {
        Chi_x_star_LR = (s_x_R - V_x_R[num_species][idx])/(s_x_R - s_x_star);
        
        for (int si = 0; si < num_species; si++)
        {
            Q_x_star_LR[si] = Chi_x_star_LR*V_x_R[si][idx];
        }
        Q_x_star_LR[num_species] = Chi_x_star_LR*rho_x_R[idx]*s_x_star;
        Q_x_star_LR[num_species + 1] = Chi_x_star_LR*Q_x_R[num_species + 1];
        Q_x_star_LR[num_species + 2] = Chi_x_star_LR*Q_x_R[num_species + 2];
        Q_x_star_LR[num_species + 3] = Chi_x_star_LR*(Q_x_R[num_species + 3] +
            (s_x_star - V_x_R[num_species][idx])*(rho_x_R[idx]*s_x_star + V_x_R[num_species + 3][idx]/
            (s_x_R - V_x_R[num_species][idx])));
        for (int si = 0; si < num_species - 1; si++)
        {
            Q_x_star_LR[num_species + 4 + si] = Chi_x_star_LR*V_x_R[num_species + 4 + si][idx];
        }
        
        for (int ei = 0; ei < num_eqn; ei++)
        {
            F_x_HLLC[ei] = F_x_R[ei] + s_x_plus*(Q_x_star_LR[ei] - Q_x_R[ei]);
        }
    }
    
    /*
     * Calulate the weights beta for hybridization.
     */
    
    const double u_x_diff = V_x_R[num_species][idx] - V_x_L[num_species][idx];
    const double v_x_diff = V_x_R[num_species + 1][idx] - V_x_L[num_species + 1][idx];
    const double w_x_diff = V_x_R[num_species + 2][idx] - V_x_L[num_species + 2][idx];
    const double vel_mag = sqrt(u_x_diff*u_x_diff + v_x_diff*v_x_diff + w_x_diff*w_x_diff);
    
    double alpha_1, alpha_2;
    if (vel_mag < EPSILON)
    {
        alpha_1 = double(1);
        alpha_2 = double(0);
    }
    else
    {
        alpha_1 = fabs(u_x_diff)/vel_mag;
        alpha_2 = sqrt(double(1) - alpha_1*alpha_1);
    }
    
    const double beta_1 = double(1)/double(2)*(double(1) + alpha_1/(alpha_1 + alpha_2));
    const double beta_2 = double(1) - beta_1;
    
    for (int si = 0; si < num_species; si++)
    {
        F_x[si][idx_flux] = beta_1*F_x_HLLC[si] + beta_2*F_x_HLL[si];
    }
    F_x[num_species][idx_flux] = F_x_HLLC[num_species];
    F_x[num_species + 1][idx_flux] = beta_1*F_x_HLLC[num_species + 1] + beta_2*F_x_HLL[num_species];
    F_x[num_species + 2][idx_flux] = beta_1*F_x_HLLC[num_species + 2] + beta_2*F_x_HLL[num_species + 1];
    F_x[num_species + 3][idx_flux] = F_x_HLLC[num_species + 3];
    for (int si = 0; si < num_species - 1; si++)
    {
        F_x[num_species + 4 + si][idx_flux] = beta_1*F_x_HLLC[num_species + 4 + si] +
            beta_2*F_x_HLL[num_species + 2 + si];
    }
}


/*
 * Compute the local convective flux in the y-direction from primitive variables with
 * 2D HLLC-HLL Riemann solver.
 */
static inline __attribute__((always_inline)) void
computeLocalConvectiveFluxInYDirectionFromPrimitiveVariablesHLLC_HLL2D(
    double** F_y,
    double** V_y_B,
    double** V_y_T,
    double* rho_y_B,
    double* rho_y_T,
    double* c_y_B,
    double* c_y_T,
    double* epsilon_y_B,
    double* epsilon_y_T,
    double& s_y_minus,
    double& s_y_plus,
    double& s_y_star,
    double& Chi_y_star_BT,
    const int& idx_flux,
    const int& idx,
    const int& num_species,
    const int& num_eqn)
{
    const double v_y_average = double(1)/double(2)*(V_y_B[num_species + 1][idx] + V_y_T[num_species + 1][idx]);
    const double c_y_average = double(1)/double(2)*(c_y_B[idx] + c_y_T[idx]);
    
    const double s_y_B = fmin(v_y_average - c_y_average, V_y_B[num_species + 1][idx] - c_y_B[idx]);
    const double s_y_T = fmax(v_y_average + c_y_average, V_y_T[num_species + 1][idx] + c_y_T[idx]);
    
    s_y_minus = fmin(double(0), s_y_B);
    s_y_plus  = fmax(double(0), s_y_T);
    
    s_y_star = (V_y_T[num_species + 2][idx] - V_y_B[num_species + 2][idx] +
        rho_y_B[idx]*V_y_B[num_species + 1][idx]*(s_y_B - V_y_B[num_species + 1][idx]) -
        rho_y_T[idx]*V_y_T[num_species + 1][idx]*(s_y_T - V_y_T[num_species + 1][idx]))/
        (rho_y_B[idx]*(s_y_B - V_y_B[num_species + 1][idx]) -
        rho_y_T[idx]*(s_y_T - V_y_T[num_species + 1][idx]));
    
    double Q_y_B[num_eqn];
    double Q_y_T[num_eqn];
    double F_y_B[num_eqn];
    double F_y_T[num_eqn];
    double F_y_HLL[2*num_species];
    double F_y_HLLC[num_eqn];
    double Q_y_star_BT[num_eqn];
    
    for (int si = 0; si < num_species; si++)
    {
        Q_y_B[si] = V_y_B[si][idx];
    }
    Q_y_B[num_species] = rho_y_B[idx]*V_y_B[num_species][idx];
    Q_y_B[num_species + 1] = rho_y_B[idx]*V_y_B[num_species + 1][idx];
    Q_y_B[num_species + 2] = rho_y_B[idx]*(epsilon_y_B[idx] +
        double(1)/double(2)*(V_y_B[num_species][idx]*V_y_B[num_species][idx] +
        V_y_B[num_species + 1][idx]*V_y_B[num_species + 1][idx]));
    for (int si = 0; si < num_species - 1; si++)
    {
        Q_y_B[num_species + 3 + si] = V_y_B[num_species + 3 + si][idx];
    }
    
    for (int si = 0; si < num_species; si++)
    {
        Q_y_T[si] = V_y_T[si][idx];
    }
    Q_y_T[num_species] = rho_y_T[idx]*V_y_T[num_species][idx];
    Q_y_T[num_species + 1] = rho_y_T[idx]*V_y_T[num_species + 1][idx];
    Q_y_T[num_species + 2] = rho_y_T[idx]*(epsilon_y_T[idx] +
        double(1)/double(2)*(V_y_T[num_species][idx]*V_y_T[num_species][idx] +
        V_y_T[num_species + 1][idx]*V_y_T[num_species + 1][idx]));
    for (int si = 0; si < num_species - 1; si++)
    {
        Q_y_T[num_species + 3 + si] = V_y_T[num_species + 3 + si][idx];
    }
    
    for (int si = 0; si < num_species; si++)
    {
        F_y_B[si] = V_y_B[num_species + 1][idx]*V_y_B[si][idx];
    }
    F_y_B[num_species] = V_y_B[num_species + 1][idx]*Q_y_B[num_species];
    F_y_B[num_species + 1] = V_y_B[num_species + 1][idx]*Q_y_B[num_species + 1] + V_y_B[num_species + 2][idx];
    F_y_B[num_species + 2] = V_y_B[num_species + 1][idx]*(Q_y_B[num_species + 2] + V_y_B[num_species + 2][idx]);
    for (int si = 0; si < num_species - 1; si++)
    {
        F_y_B[num_species + 3 + si] = V_y_B[num_species + 1][idx]*V_y_B[num_species + 3 + si][idx];
    }
    
    for (int si = 0; si < num_species; si++)
    {
        F_y_T[si] = V_y_T[num_species + 1][idx]*V_y_T[si][idx];
    }
    F_y_T[num_species] = V_y_T[num_species + 1][idx]*Q_y_T[num_species];
    F_y_T[num_species + 1] = V_y_T[num_species + 1][idx]*Q_y_T[num_species + 1] + V_y_T[num_species + 2][idx];
    F_y_T[num_species + 2] = V_y_T[num_species + 1][idx]*(Q_y_T[num_species + 2] + V_y_T[num_species + 2][idx]);
    for (int si = 0; si < num_species - 1; si++)
    {
        F_y_T[num_species + 3 + si] = V_y_T[num_species + 1][idx]*V_y_T[num_species + 3 + si][idx];
    }
    
    for (int si = 0; si < num_species; si++)
    {
        F_y_HLL[si] = (s_y_T*F_y_B[si] - s_y_B*F_y_T[si] + s_y_T*s_y_B*(Q_y_T[si] - Q_y_B[si]))/
            (s_y_T - s_y_B);
    }
    F_y_HLL[num_species] = (s_y_T*F_y_B[num_species] - s_y_B*F_y_T[num_species] + s_y_T*s_y_B*
        (Q_y_T[num_species] - Q_y_B[num_species]))/(s_y_T - s_y_B);
    for (int si = 0; si < num_species - 1; si++)
    {
        F_y_HLL[num_species + 1 + si] = (s_y_T*F_y_B[num_species + 3 + si] - s_y_B*F_y_T[num_species + 3 + si] +
            s_y_T*s_y_B*(Q_y_T[num_species + 3 + si] - Q_y_B[num_species + 3 + si]))/(s_y_T - s_y_B);
    }
    
    if (s_y_B > double(0))
    {
        for (int si = 0; si < num_species; si++)
        {
            F_y_HLL[si] = F_y_B[si];
        }
        F_y_HLL[num_species] = F_y_B[num_species];
        for (int si = 0; si < num_species - 1; si++)
        {
            F_y_HLL[num_species + 1 + si] = F_y_B[num_species + 3 + si];
        }
    }
    
    if (s_y_T < double(0))
    {
        for (int si = 0; si < num_species; si++)
        {
            F_y_HLL[si] = F_y_T[si];
        }
        F_y_HLL[num_species] = F_y_T[num_species];
        for (int si = 0; si < num_species - 1; si++)
        {
            F_y_HLL[num_species + 1 + si] = F_y_T[num_species + 3 + si];
        }
    }
    
    if (s_y_star > double(0))
    {
        Chi_y_star_BT = (s_y_B - V_y_B[num_species + 1][idx])/(s_y_B - s_y_star);
        
        for (int si = 0; si < num_species; si++)
        {
            Q_y_star_BT[si] = Chi_y_star_BT*V_y_B[si][idx];
        }
        Q_y_star_BT[num_species] = Chi_y_star_BT*Q_y_B[num_species];
        Q_y_star_BT[num_species + 1] = Chi_y_star_BT*rho_y_B[idx]*s_y_star;
        Q_y_star_BT[num_species + 2] = Chi_y_star_BT*(Q_y_B[num_species + 2] +
            (s_y_star - V_y_B[num_species + 1][idx])*(rho_y_B[idx]*s_y_star + V_y_B[num_species + 2][idx]/
            (s_y_B - V_y_B[num_species + 1][idx])));
        for (int si = 0; si < num_species - 1; si++)
        {
            Q_y_star_BT[num_species + 3 + si] = Chi_y_star_BT*V_y_B[num_species + 3 + si][idx];
        }
        
        for (int ei = 0; ei < num_eqn; ei++)
        {
            F_y_HLLC[ei] = F_y_B[ei] + s_y_minus*(Q_y_star_BT[ei] - Q_y_B[ei]);
        }
    }
    else
    {
        Chi_y_star_BT = (s_y_T - V_y_T[num_species + 1][idx])/(s_y_T - s_y_star);
        
        for (int si = 0; si < num_species; si++)
        {
            Q_y_star_BT[si] = Chi_y_star_BT*V_y_T[si][idx];
        }
        Q_y_star_BT[num_species] = Chi_y_star_BT*Q_y_T[num_species];
        Q_y_star_BT[num_species + 1] = Chi_y_star_BT*rho_y_T[idx]*s_y_star;
        Q_y_star_BT[num_species + 2] = Chi_y_star_BT*(Q_y_T[num_species + 2] +
            (s_y_star - V_y_T[num_species + 1][idx])*(rho_y_T[idx]*s_y_star + V_y_T[num_species + 2][idx]/
            (s_y_T - V_y_T[num_species + 1][idx])));
        for (int si = 0; si < num_species - 1; si++)
        {
            Q_y_star_BT[num_species + 3 + si] = Chi_y_star_BT*V_y_T[num_species + 3 + si][idx];
        }
        
        for (int ei = 0; ei < num_eqn; ei++)
        {
            F_y_HLLC[ei] = F_y_T[ei] + s_y_plus*(Q_y_star_BT[ei] - Q_y_T[ei]);
        }
    }
    
    /*
     * Calulate the weights beta for hybridization.
     */
    
    const double u_y_diff = V_y_T[num_species][idx] - V_y_B[num_species][idx];
    const double v_y_diff = V_y_T[num_species + 1][idx] - V_y_B[num_species + 1][idx];
    const double vel_mag = sqrt(u_y_diff*u_y_diff + v_y_diff*v_y_diff);
    
    double alpha_1, alpha_2;
    if (vel_mag < EPSILON)
    {
        alpha_1 = double(1);
        alpha_2 = double(0);
    }
    else
    {
        alpha_1 = fabs(v_y_diff)/vel_mag;
        alpha_2 = sqrt(double(1) - alpha_1*alpha_1);
    }
    
    const double beta_1 = double(1)/double(2)*(double(1) + alpha_1/(alpha_1 + alpha_2));
    const double beta_2 = double(1) - beta_1;
    
    for (int si = 0; si < num_species; si++)
    {
        F_y[si][idx_flux] = beta_1*F_y_HLLC[si] + beta_2*F_y_HLL[si];
    }
    F_y[num_species][idx_flux] = beta_1*F_y_HLLC[num_species] + beta_2*F_y_HLL[num_species];
    F_y[num_species + 1][idx_flux] = F_y_HLLC[num_species + 1];
    F_y[num_species + 2][idx_flux] = F_y_HLLC[num_species + 2];
    for (int si = 0; si < num_species - 1; si++)
    {
        F_y[num_species + 3 + si][idx_flux] = beta_1*F_y_HLLC[num_species + 3 + si] +
            beta_2*F_y_HLL[num_species + 1 + si];
    }
}


/*
 * Compute the local convective flux in the y-direction from primitive variables with
 * 3D HLLC-HLL Riemann solver.
 */
static inline __attribute__((always_inline)) void
computeLocalConvectiveFluxInYDirectionFromPrimitiveVariablesHLLC_HLL3D(
    double** F_y,
    double** V_y_B,
    double** V_y_T,
    double* rho_y_B,
    double* rho_y_T,
    double* c_y_B,
    double* c_y_T,
    double* epsilon_y_B,
    double* epsilon_y_T,
    double& s_y_minus,
    double& s_y_plus,
    double& s_y_star,
    double& Chi_y_star_BT,
    const int& idx_flux,
    const int& idx,
    const int& num_species,
    const int& num_eqn)
{
    const double v_y_average = double(1)/double(2)*(V_y_B[num_species + 1][idx] + V_y_T[num_species + 1][idx]);
    const double c_y_average = double(1)/double(2)*(c_y_B[idx] + c_y_T[idx]);
    
    const double s_y_B = fmin(v_y_average - c_y_average, V_y_B[num_species + 1][idx] - c_y_B[idx]);
    const double s_y_T = fmax(v_y_average + c_y_average, V_y_T[num_species + 1][idx] + c_y_T[idx]);
    
    s_y_minus = fmin(double(0), s_y_B);
    s_y_plus  = fmax(double(0), s_y_T);
    
    s_y_star = (V_y_T[num_species + 3][idx] - V_y_B[num_species + 3][idx] +
        rho_y_B[idx]*V_y_B[num_species + 1][idx]*(s_y_B - V_y_B[num_species + 1][idx]) -
        rho_y_T[idx]*V_y_T[num_species + 1][idx]*(s_y_T - V_y_T[num_species + 1][idx]))/
        (rho_y_B[idx]*(s_y_B - V_y_B[num_species + 1][idx]) -
        rho_y_T[idx]*(s_y_T - V_y_T[num_species + 1][idx]));
    
    double Q_y_B[num_eqn];
    double Q_y_T[num_eqn];
    double F_y_B[num_eqn];
    double F_y_T[num_eqn];
    double F_y_HLL[2*num_species + 1];
    double F_y_HLLC[num_eqn];
    double Q_y_star_BT[num_eqn];
    
    for (int si = 0; si < num_species; si++)
    {
        Q_y_B[si] = V_y_B[si][idx];
    }
    Q_y_B[num_species] = rho_y_B[idx]*V_y_B[num_species][idx];
    Q_y_B[num_species + 1] = rho_y_B[idx]*V_y_B[num_species + 1][idx];
    Q_y_B[num_species + 2] = rho_y_B[idx]*V_y_B[num_species + 2][idx];
    Q_y_B[num_species + 3] = rho_y_B[idx]*(epsilon_y_B[idx] +
        double(1)/double(2)*(V_y_B[num_species][idx]*V_y_B[num_species][idx] +
        V_y_B[num_species + 1][idx]*V_y_B[num_species + 1][idx] +
        V_y_B[num_species + 2][idx]*V_y_B[num_species + 2][idx]));
    for (int si = 0; si < num_species - 1; si++)
    {
        Q_y_B[num_species + 4 + si] = V_y_B[num_species + 4 + si][idx];
    }
    
    for (int si = 0; si < num_species; si++)
    {
        Q_y_T[si] = V_y_T[si][idx];
    }
    Q_y_T[num_species] = rho_y_T[idx]*V_y_T[num_species][idx];
    Q_y_T[num_species + 1] = rho_y_T[idx]*V_y_T[num_species + 1][idx];
    Q_y_T[num_species + 2] = rho_y_T[idx]*V_y_T[num_species + 2][idx];
    Q_y_T[num_species + 3] = rho_y_T[idx]*(epsilon_y_T[idx] +
        double(1)/double(2)*(V_y_T[num_species][idx]*V_y_T[num_species][idx] +
        V_y_T[num_species + 1][idx]*V_y_T[num_species + 1][idx] +
        V_y_T[num_species + 2][idx]*V_y_T[num_species + 2][idx]));
    for (int si = 0; si < num_species - 1; si++)
    {
        Q_y_T[num_species + 4 + si] = V_y_T[num_species + 4 + si][idx];
    }
    
    for (int si = 0; si < num_species; si++)
    {
        F_y_B[si] = V_y_B[num_species + 1][idx]*V_y_B[si][idx];
    }
    F_y_B[num_species] = V_y_B[num_species + 1][idx]*Q_y_B[num_species];
    F_y_B[num_species + 1] = V_y_B[num_species + 1][idx]*Q_y_B[num_species + 1] + V_y_B[num_species + 3][idx];
    F_y_B[num_species + 2] = V_y_B[num_species + 1][idx]*Q_y_B[num_species + 2];
    F_y_B[num_species + 3] = V_y_B[num_species + 1][idx]*(Q_y_B[num_species + 3] + V_y_B[num_species + 3][idx]);
    for (int si = 0; si < num_species - 1; si++)
    {
        F_y_B[num_species + 4 + si] = V_y_B[num_species + 1][idx]*V_y_B[num_species + 4 + si][idx];
    }
    
    for (int si = 0; si < num_species; si++)
    {
        F_y_T[si] = V_y_T[num_species + 1][idx]*V_y_T[si][idx];
    }
    F_y_T[num_species] = V_y_T[num_species + 1][idx]*Q_y_T[num_species];
    F_y_T[num_species + 1] = V_y_T[num_species + 1][idx]*Q_y_T[num_species + 1] + V_y_T[num_species + 3][idx];
    F_y_T[num_species + 2] = V_y_T[num_species + 1][idx]*Q_y_T[num_species + 2];
    F_y_T[num_species + 3] = V_y_T[num_species + 1][idx]*(Q_y_T[num_species + 3] + V_y_T[num_species + 3][idx]);
    for (int si = 0; si < num_species - 1; si++)
    {
        F_y_T[num_species + 4 + si] = V_y_T[num_species + 1][idx]*V_y_T[num_species + 4 + si][idx];
    }
    
    for (int si = 0; si < num_species; si++)
    {
        F_y_HLL[si] = (s_y_T*F_y_B[si] - s_y_B*F_y_T[si] + s_y_T*s_y_B*(Q_y_T[si] - Q_y_B[si]))/
            (s_y_T - s_y_B);
    }
    F_y_HLL[num_species] = (s_y_T*F_y_B[num_species] - s_y_B*F_y_T[num_species] + s_y_T*s_y_B*
        (Q_y_T[num_species] - Q_y_B[num_species]))/(s_y_T - s_y_B);
    F_y_HLL[num_species + 1] = (s_y_T*F_y_B[num_species + 2] - s_y_B*F_y_T[num_species + 2] + s_y_T*s_y_B*
        (Q_y_T[num_species + 2] - Q_y_B[num_species + 2]))/(s_y_T - s_y_B);
    for (int si = 0; si < num_species - 1; si++)
    {
        F_y_HLL[num_species + 2 + si] = (s_y_T*F_y_B[num_species + 4 + si] - s_y_B*F_y_T[num_species + 4 + si] +
            s_y_T*s_y_B*(Q_y_T[num_species + 4 + si] - Q_y_B[num_species + 4 + si]))/(s_y_T - s_y_B);
    }
    
    if (s_y_B > double(0))
    {
        for (int si = 0; si < num_species; si++)
        {
            F_y_HLL[si] = F_y_B[si];
        }
        F_y_HLL[num_species] = F_y_B[num_species];
        F_y_HLL[num_species + 1] = F_y_B[num_species + 2];
        for (int si = 0; si < num_species - 1; si++)
        {
            F_y_HLL[num_species + 2 + si] = F_y_B[num_species + 4 + si];
        }
    }
    
    if (s_y_T < double(0))
    {
        for (int si = 0; si < num_species; si++)
        {
            F_y_HLL[si] = F_y_T[si];
        }
        F_y_HLL[num_species] = F_y_T[num_species];
        F_y_HLL[num_species + 1] = F_y_T[num_species + 2];
        for (int si = 0; si < num_species - 1; si++)
        {
            F_y_HLL[num_species + 2 + si] = F_y_T[num_species + 4 + si];
        }
    }
    
    if (s_y_star > double(0))
    {
        Chi_y_star_BT = (s_y_B - V_y_B[num_species + 1][idx])/(s_y_B - s_y_star);
        
        for (int si = 0; si < num_species; si++)
        {
            Q_y_star_BT[si] = Chi_y_star_BT*V_y_B[si][idx];
        }
        Q_y_star_BT[num_species] = Chi_y_star_BT*Q_y_B[num_species];
        Q_y_star_BT[num_species + 1] = Chi_y_star_BT*rho_y_B[idx]*s_y_star;
        Q_y_star_BT[num_species + 2] = Chi_y_star_BT*Q_y_B[num_species + 2];
        Q_y_star_BT[num_species + 3] = Chi_y_star_BT*(Q_y_B[num_species + 3] +
            (s_y_star - V_y_B[num_species + 1][idx])*(rho_y_B[idx]*s_y_star + V_y_B[num_species + 3][idx]/
            (s_y_B - V_y_B[num_species + 1][idx])));
        for (int si = 0; si < num_species - 1; si++)
        {
            Q_y_star_BT[num_species + 4 + si] = Chi_y_star_BT*V_y_B[num_species + 4 + si][idx];
        }
        
        for (int ei = 0; ei < num_eqn; ei++)
        {
            F_y_HLLC[ei] = F_y_B[ei] + s_y_minus*(Q_y_star_BT[ei] - Q_y_B[ei]);
        }
    }
    else
    {
        Chi_y_star_BT = (s_y_T - V_y_T[num_species + 1][idx])/(s_y_T - s_y_star);
        
        for (int si = 0; si < num_species; si++)
        {
            Q_y_star_BT[si] = Chi_y_star_BT*V_y_T[si][idx];
        }
        Q_y_star_BT[num_species] = Chi_y_star_BT*Q_y_T[num_species];
        Q_y_star_BT[num_species + 1] = Chi_y_star_BT*rho_y_T[idx]*s_y_star;
        Q_y_star_BT[num_species + 2] = Chi_y_star_BT*Q_y_T[num_species + 2];
        Q_y_star_BT[num_species + 3] = Chi_y_star_BT*(Q_y_T[num_species + 3] +
            (s_y_star - V_y_T[num_species + 1][idx])*(rho_y_T[idx]*s_y_star + V_y_T[num_species + 3][idx]/
            (s_y_T - V_y_T[num_species + 1][idx])));
        for (int si = 0; si < num_species - 1; si++)
        {
            Q_y_star_BT[num_species + 4 + si] = Chi_y_star_BT*V_y_T[num_species + 4 + si][idx];
        }
        
        for (int ei = 0; ei < num_eqn; ei++)
        {
            F_y_HLLC[ei] = F_y_T[ei] + s_y_plus*(Q_y_star_BT[ei] - Q_y_T[ei]);
        }
    }
    
    /*
     * Calulate the weights beta for hybridization.
     */
    
    const double u_y_diff = V_y_T[num_species][idx] - V_y_B[num_species][idx];
    const double v_y_diff = V_y_T[num_species + 1][idx] - V_y_B[num_species + 1][idx];
    const double w_y_diff = V_y_T[num_species + 2][idx] - V_y_B[num_species + 2][idx];
    const double vel_mag = sqrt(u_y_diff*u_y_diff + v_y_diff*v_y_diff + w_y_diff*w_y_diff);
    
    double alpha_1, alpha_2;
    if (vel_mag < EPSILON)
    {
        alpha_1 = double(1);
        alpha_2 = double(0);
    }
    else
    {
        alpha_1 = fabs(v_y_diff)/vel_mag;
        alpha_2 = sqrt(double(1) - alpha_1*alpha_1);
    }
    
    const double beta_1 = double(1)/double(2)*(double(1) + alpha_1/(alpha_1 + alpha_2));
    const double beta_2 = double(1) - beta_1;
    
    for (int si = 0; si < num_species; si++)
    {
        F_y[si][idx_flux] = beta_1*F_y_HLLC[si] + beta_2*F_y_HLL[si];
    }
    F_y[num_species][idx_flux] = beta_1*F_y_HLLC[num_species] + beta_2*F_y_HLL[num_species];
    F_y[num_species + 1][idx_flux] = F_y_HLLC[num_species + 1];
    F_y[num_species + 2][idx_flux] = beta_1*F_y_HLLC[num_species + 2] + beta_2*F_y_HLL[num_species + 1];
    F_y[num_species + 3][idx_flux] = F_y_HLLC[num_species + 3];
    for (int si = 0; si < num_species - 1; si++)
    {
        F_y[num_species + 4 + si][idx_flux] = beta_1*F_y_HLLC[num_species + 4 + si] +
            beta_2*F_y_HLL[num_species + 2 + si];
    }
}


/*
 * Compute the local convective flux in the z-direction from primitive variables with
 * 3D HLLC-HLL Riemann solver.
 */
static inline __attribute__((always_inline)) void
computeLocalConvectiveFluxInZDirectionFromPrimitiveVariablesHLLC_HLL3D(
    double** F_z,
    double** V_z_B,
    double** V_z_F,
    double* rho_z_B,
    double* rho_z_F,
    double* c_z_B,
    double* c_z_F,
    double* epsilon_z_B,
    double* epsilon_z_F,
    double& s_z_minus,
    double& s_z_plus,
    double& s_z_star,
    double& Chi_z_star_BF,
    const int& idx_flux,
    const int& idx,
    const int& num_species,
    const int& num_eqn)
{
    const double w_z_average = double(1)/double(2)*(V_z_B[num_species + 2][idx] + V_z_F[num_species + 2][idx]);
    const double c_z_average = double(1)/double(2)*(c_z_B[idx] + c_z_F[idx]);
    
    const double s_z_B = fmin(w_z_average - c_z_average, V_z_B[num_species + 2][idx] - c_z_B[idx]);
    const double s_z_F = fmax(w_z_average + c_z_average, V_z_F[num_species + 2][idx] + c_z_F[idx]);
    
    s_z_minus = fmin(double(0), s_z_B);
    s_z_plus  = fmax(double(0), s_z_F);
    
    s_z_star = (V_z_F[num_species + 3][idx] - V_z_B[num_species + 3][idx] +
        rho_z_B[idx]*V_z_B[num_species + 2][idx]*(s_z_B - V_z_B[num_species + 2][idx]) -
        rho_z_F[idx]*V_z_F[num_species + 2][idx]*(s_z_F - V_z_F[num_species + 2][idx]))/
        (rho_z_B[idx]*(s_z_B - V_z_B[num_species + 2][idx]) -
        rho_z_F[idx]*(s_z_F - V_z_F[num_species + 2][idx]));
    
    double Q_z_B[num_eqn];
    double Q_z_F[num_eqn];
    double F_z_B[num_eqn];
    double F_z_F[num_eqn];
    double F_z_HLL[2*num_species + 1];
    double F_z_HLLC[num_eqn];
    double Q_z_star_BF[num_eqn];
    
    for (int si = 0; si < num_species; si++)
    {
        Q_z_B[si] = V_z_B[si][idx];
    }
    Q_z_B[num_species] = rho_z_B[idx]*V_z_B[num_species][idx];
    Q_z_B[num_species + 1] = rho_z_B[idx]*V_z_B[num_species + 1][idx];
    Q_z_B[num_species + 2] = rho_z_B[idx]*V_z_B[num_species + 2][idx];
    Q_z_B[num_species + 3] = rho_z_B[idx]*(epsilon_z_B[idx] +
        double(1)/double(2)*(V_z_B[num_species][idx]*V_z_B[num_species][idx] +
        V_z_B[num_species + 1][idx]*V_z_B[num_species + 1][idx] +
        V_z_B[num_species + 2][idx]*V_z_B[num_species + 2][idx]));
    for (int si = 0; si < num_species - 1; si++)
    {
        Q_z_B[num_species + 4 + si] = V_z_B[num_species + 4 + si][idx];
    }
    
    for (int si = 0; si < num_species; si++)
    {
        Q_z_F[si] = V_z_F[si][idx];
    }
    Q_z_F[num_species] = rho_z_F[idx]*V_z_F[num_species][idx];
    Q_z_F[num_species + 1] = rho_z_F[idx]*V_z_F[num_species + 1][idx];
    Q_z_F[num_species + 2] = rho_z_F[idx]*V_z_F[num_species + 2][idx];
    Q_z_F[num_species + 3] = rho_z_F[idx]*(epsilon_z_F[idx] +
        double(1)/double(2)*(V_z_F[num_species][idx]*V_z_F[num_species][idx] +
        V_z_F[num_species + 1][idx]*V_z_F[num_species + 1][idx] +
        V_z_F[num_species + 2][idx]*V_z_F[num_species + 2][idx]));
    for (int si = 0; si < num_species - 1; si++)
    {
        Q_z_F[num_species + 4 + si] = V_z_F[num_species + 4 + si][idx];
    }
    
    for (int si = 0; si < num_species; si++)
    {
        F_z_B[si] = V_z_B[num_species + 2][idx]*V_z_B[si][idx];
    }
    F_z_B[num_species] = V_z_B[num_species + 2][idx]*Q_z_B[num_species];
    F_z_B[num_species + 1] = V_z_B[num_species + 2][idx]*Q_z_B[num_species + 1];
    F_z_B[num_species + 2] = V_z_B[num_species + 2][idx]*Q_z_B[num_species + 2] + V_z_B[num_species + 3][idx];
    F_z_B[num_species + 3] = V_z_B[num_species + 2][idx]*(Q_z_B[num_species + 3] + V_z_B[num_species + 3][idx]);
    for (int si = 0; si < num_species - 1; si++)
    {
        F_z_B[num_species + 4 + si] = V_z_B[num_species + 2][idx]*V_z_B[num_species + 4 + si][idx];
    }
    
    for (int si = 0; si < num_species; si++)
    {
        F_z_F[si] = V_z_F[num_species + 2][idx]*V_z_F[si][idx];
    }
    F_z_F[num_species] = V_z_F[num_species + 2][idx]*Q_z_F[num_species];
    F_z_F[num_species + 1] = V_z_F[num_species + 2][idx]*Q_z_F[num_species + 1];
    F_z_F[num_species + 2] = V_z_F[num_species + 2][idx]*Q_z_F[num_species + 2] + V_z_F[num_species + 3][idx];
    F_z_F[num_species + 3] = V_z_F[num_species + 2][idx]*(Q_z_F[num_species + 3] + V_z_F[num_species + 3][idx]);
    for (int si = 0; si < num_species - 1; si++)
    {
        F_z_F[num_species + 4 + si] = V_z_F[num_species + 2][idx]*V_z_F[num_species + 4 + si][idx];
    }
    
    for (int si = 0; si < num_species; si++)
    {
        F_z_HLL[si] = (s_z_F*F_z_B[si] - s_z_B*F_z_F[si] + s_z_F*s_z_B*(Q_z_F[si] - Q_z_B[si]))/
            (s_z_F - s_z_B);
    }
    F_z_HLL[num_species] = (s_z_F*F_z_B[num_species] - s_z_B*F_z_F[num_species] + s_z_F*s_z_B*
        (Q_z_F[num_species] - Q_z_B[num_species]))/(s_z_F - s_z_B);
    F_z_HLL[num_species + 1] = (s_z_F*F_z_B[num_species + 1] - s_z_B*F_z_F[num_species + 1] + s_z_F*s_z_B*
        (Q_z_F[num_species + 1] - Q_z_B[num_species + 1]))/(s_z_F - s_z_B);
    for (int si = 0; si < num_species - 1; si++)
    {
        F_z_HLL[num_species + 2 + si] = (s_z_F*F_z_B[num_species + 4 + si] - s_z_B*F_z_F[num_species + 4 + si] +
            s_z_F*s_z_B*(Q_z_F[num_species + 4 + si] - Q_z_B[num_species + 4 + si]))/(s_z_F - s_z_B);
    }
    
    if (s_z_B > double(0))
    {
        for (int si = 0; si < num_species; si++)
        {
            F_z_HLL[si] = F_z_B[si];
        }
        F_z_HLL[num_species] = F_z_B[num_species];
        F_z_HLL[num_species + 1] = F_z_B[num_species + 1];
        for (int si = 0; si < num_species - 1; si++)
        {
            F_z_HLL[num_species + 2 + si] = F_z_B[num_species + 4 + si];
        }
    }
    
    if (s_z_F < double(0))
    {
        for (int si = 0; si < num_species; si++)
        {
            F_z_HLL[si] = F_z_F[si];
        }
        F_z_HLL[num_species] = F_z_F[num_species];
        F_z_HLL[num_species + 1] = F_z_F[num_species + 1];
        for (int si = 0; si < num_species - 1; si++)
        {
            F_z_HLL[num_species + 2 + si] = F_z_F[num_species + 4 + si];
        }
    }
    
    if (s_z_star > double(0))
    {
        Chi_z_star_BF = (s_z_B - V_z_B[num_species + 2][idx])/(s_z_B - s_z_star);
        
        for (int si = 0; si < num_species; si++)
        {
            Q_z_star_BF[si] = Chi_z_star_BF*V_z_B[si][idx];
        }
        Q_z_star_BF[num_species] = Chi_z_star_BF*Q_z_B[num_species];
        Q_z_star_BF[num_species + 1] = Chi_z_star_BF*Q_z_B[num_species + 1];
        Q_z_star_BF[num_species + 2] = Chi_z_star_BF*rho_z_B[idx]*s_z_star;
        Q_z_star_BF[num_species + 3] = Chi_z_star_BF*(Q_z_B[num_species + 3] +
            (s_z_star - V_z_B[num_species + 2][idx])*(rho_z_B[idx]*s_z_star + V_z_B[num_species + 3][idx]/
            (s_z_B - V_z_B[num_species + 2][idx])));
        for (int si = 0; si < num_species - 1; si++)
        {
            Q_z_star_BF[num_species + 4 + si] = Chi_z_star_BF*V_z_B[num_species + 4 + si][idx];
        }
        
        for (int ei = 0; ei < num_eqn; ei++)
        {
            F_z_HLLC[ei] = F_z_B[ei] + s_z_minus*(Q_z_star_BF[ei] - Q_z_B[ei]);
        }
    }
    else
    {
        Chi_z_star_BF = (s_z_F - V_z_F[num_species + 2][idx])/(s_z_F - s_z_star);
        
        for (int si = 0; si < num_species; si++)
        {
            Q_z_star_BF[si] = Chi_z_star_BF*V_z_F[si][idx];
        }
        Q_z_star_BF[num_species] = Chi_z_star_BF*Q_z_F[num_species];
        Q_z_star_BF[num_species + 1] = Chi_z_star_BF*Q_z_F[num_species + 1];
        Q_z_star_BF[num_species + 2] = Chi_z_star_BF*rho_z_F[idx]*s_z_star;
        Q_z_star_BF[num_species + 3] = Chi_z_star_BF*(Q_z_F[num_species + 3] +
            (s_z_star - V_z_F[num_species + 2][idx])*(rho_z_F[idx]*s_z_star + V_z_F[num_species + 3][idx]/
            (s_z_F - V_z_F[num_species + 2][idx])));
        for (int si = 0; si < num_species - 1; si++)
        {
            Q_z_star_BF[num_species + 4 + si] = Chi_z_star_BF*V_z_F[num_species + 4 + si][idx];
        }
        
        for (int ei = 0; ei < num_eqn; ei++)
        {
            F_z_HLLC[ei] = F_z_F[ei] + s_z_plus*(Q_z_star_BF[ei] - Q_z_F[ei]);
        }
    }
    
    /*
     * Calulate the weights beta for hybridization.
     */
    
    const double u_z_diff = V_z_F[num_species][idx] - V_z_B[num_species][idx];
    const double v_z_diff = V_z_F[num_species + 1][idx] - V_z_B[num_species + 1][idx];
    const double w_z_diff = V_z_F[num_species + 2][idx] - V_z_B[num_species + 2][idx];
    const double vel_mag = sqrt(u_z_diff*u_z_diff + v_z_diff*v_z_diff + w_z_diff*w_z_diff);
    
    double alpha_1, alpha_2;
    if (vel_mag < EPSILON)
    {
        alpha_1 = double(1);
        alpha_2 = double(0);
    }
    else
    {
        alpha_1 = fabs(w_z_diff)/vel_mag;
        alpha_2 = sqrt(double(1) - alpha_1*alpha_1);
    }
    
    const double beta_1 = double(1)/double(2)*(double(1) + alpha_1/(alpha_1 + alpha_2));
    const double beta_2 = double(1) - beta_1;
    
    for (int si = 0; si < num_species; si++)
    {
        F_z[si][idx_flux] = beta_1*F_z_HLLC[si] + beta_2*F_z_HLL[si];
    }
    F_z[num_species][idx_flux] = beta_1*F_z_HLLC[num_species] + beta_2*F_z_HLL[num_species];
    F_z[num_species + 1][idx_flux] = beta_1*F_z_HLLC[num_species + 1] + beta_2*F_z_HLL[num_species + 1];
    F_z[num_species + 2][idx_flux] = F_z_HLLC[num_species + 2];
    F_z[num_species + 3][idx_flux] = F_z_HLLC[num_species + 3];
    for (int si = 0; si < num_species - 1; si++)
    {
        F_z[num_species + 4 + si][idx_flux] = beta_1*F_z_HLLC[num_species + 4 + si] +
            beta_2*F_z_HLL[num_species + 2 + si];
    }
}


/*
 * Compute the convective flux and velocity in the x-direction from conservative variables with
 * HLLC Riemann solver.
 */
void
FlowModelRiemannSolverFiveEqnAllaire::computeConvectiveFluxAndVelocityInXDirectionFromConservativeVariablesHLLC_HLL(
    boost::shared_ptr<pdat::SideData<double> > convective_flux,
    boost::shared_ptr<pdat::SideData<double> > velocity,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables_L,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables_R,
    const hier::Box& domain,
    bool compute_velocity) const
{
    boost::shared_ptr<FlowModel> flow_model_tmp = d_flow_model.lock();
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
    
    const boost::shared_ptr<EquationOfStateMixingRules> equation_of_state_mixing_rules =
        flow_model_tmp->getEquationOfStateMixingRules();
    
    /*
     * Get the pointers to the side data of convective flux and conservative variables.
     */
    
    std::vector<double*> F_x;
    F_x.reserve(num_eqn);
    for (int ei = 0; ei < num_eqn; ei++)
    {
        F_x.push_back(convective_flux->getPointer(0, ei));
    }
    
    std::vector<double*> Q_x_L;
    std::vector<double*> Q_x_R;
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
    
    boost::shared_ptr<pdat::SideData<double> > density_x_L(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_variables,
            direction_x));
    
    boost::shared_ptr<pdat::SideData<double> > density_x_R(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_variables,
            direction_x));
    
    boost::shared_ptr<pdat::SideData<double> > internal_energy_x_L(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_variables,
            direction_x));
    
    boost::shared_ptr<pdat::SideData<double> > internal_energy_x_R(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_variables,
            direction_x));
    
    boost::shared_ptr<pdat::SideData<double> > mass_fractions_x_L(
        new pdat::SideData<double>(interior_box, d_num_species, num_ghosts_conservative_variables,
            direction_x));
    
    boost::shared_ptr<pdat::SideData<double> > mass_fractions_x_R(
        new pdat::SideData<double>(interior_box, d_num_species, num_ghosts_conservative_variables,
            direction_x));
    
    boost::shared_ptr<pdat::SideData<double> > volume_fractions_x_L(
        new pdat::SideData<double>(interior_box, d_num_species - 1, num_ghosts_conservative_variables,
            direction_x));
    
    boost::shared_ptr<pdat::SideData<double> > volume_fractions_x_R(
        new pdat::SideData<double>(interior_box, d_num_species - 1, num_ghosts_conservative_variables,
            direction_x));
    
    boost::shared_ptr<pdat::SideData<double> > pressure_x_L(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_variables,
            direction_x));
    
    boost::shared_ptr<pdat::SideData<double> > pressure_x_R(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_variables,
            direction_x));
    
    boost::shared_ptr<pdat::SideData<double> > gruneisen_parameter_x_L(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_variables,
            direction_x));
    
    boost::shared_ptr<pdat::SideData<double> > gruneisen_parameter_x_R(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_variables,
            direction_x));
    
    boost::shared_ptr<pdat::SideData<double> > partial_pressure_partial_partial_densities_x_L(
            new pdat::SideData<double>(interior_box, d_num_species, num_ghosts_conservative_variables,
            direction_x));
    
    boost::shared_ptr<pdat::SideData<double> > partial_pressure_partial_partial_densities_x_R(
            new pdat::SideData<double>(interior_box, d_num_species, num_ghosts_conservative_variables,
            direction_x));
    
    boost::shared_ptr<pdat::SideData<double> > sound_speed_x_L(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_variables,
            direction_x));
    
    boost::shared_ptr<pdat::SideData<double> > sound_speed_x_R(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_variables,
            direction_x));
    
    /*
     * Get the pointers to the temporary data.
     */
    
    double* rho_x_L = density_x_L->getPointer(0, 0);
    double* rho_x_R = density_x_R->getPointer(0, 0);
    
    double* epsilon_x_L = internal_energy_x_L->getPointer(0, 0);
    double* epsilon_x_R = internal_energy_x_R->getPointer(0, 0);
    
    std::vector<double*> Y_x_L;
    std::vector<double*> Y_x_R;
    Y_x_L.reserve(d_num_species);
    Y_x_R.reserve(d_num_species);
    for (int si = 0; si < d_num_species; si++)
    {
        Y_x_L.push_back(mass_fractions_x_L->getPointer(0, si));
        Y_x_R.push_back(mass_fractions_x_R->getPointer(0, si));
    }
    
    std::vector<double*> Z_x_L;
    std::vector<double*> Z_x_R;
    Z_x_L.reserve(d_num_species);
    Z_x_R.reserve(d_num_species);
    for (int si = 0; si < d_num_species - 1; si++)
    {
        Z_x_L.push_back(volume_fractions_x_L->getPointer(0, si));
        Z_x_R.push_back(volume_fractions_x_R->getPointer(0, si));
    }
    
    double* p_x_L = pressure_x_L->getPointer(0, 0);
    double* p_x_R = pressure_x_R->getPointer(0, 0);
    
    double* Gamma_x_L = gruneisen_parameter_x_L->getPointer(0, 0);
    double* Gamma_x_R = gruneisen_parameter_x_R->getPointer(0, 0);
    
    std::vector<double*> Psi_x_L;
    std::vector<double*> Psi_x_R;
    Psi_x_L.reserve(d_num_species);
    Psi_x_R.reserve(d_num_species);
    for (int si = 0; si < d_num_species; si++)
    {
        Psi_x_L.push_back(partial_pressure_partial_partial_densities_x_L->getPointer(0, si));
        Psi_x_R.push_back(partial_pressure_partial_partial_densities_x_R->getPointer(0, si));
    }
    
    double* c_x_L = sound_speed_x_L->getPointer(0, 0);
    double* c_x_R = sound_speed_x_R->getPointer(0, 0);
    
    double u_x_L = double(0);
    double u_x_R = double(0);
    
    double s_x_minus = double(0);
    double s_x_plus  = double(0);
    double s_x_star  = double(0);
    
    double Chi_x_star_LR = double(0);
    
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
         * Compute the density field.
         */
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
        {
            // Compute the linear index.
            const int idx = i + num_ghosts_0_conservative_variables;
            
            rho_x_L[idx] = double(0);
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
            {
                // Compute the linear index.
                const int idx = i + num_ghosts_0_conservative_variables;
                
                rho_x_L[idx] += Q_x_L[si][idx];
            }
        }
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
        {
            // Compute the linear index.
            const int idx = i + num_ghosts_0_conservative_variables;
            
            rho_x_R[idx] = double(0);
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
            {
                // Compute the linear index.
                const int idx = i + num_ghosts_0_conservative_variables;
                
                rho_x_R[idx] += Q_x_R[si][idx];
            }
        }
        
        /*
         * Compute the internal energy field.
         */
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
        {
            // Compute the linear index.
            const int idx = i + num_ghosts_0_conservative_variables;
            
            epsilon_x_L[idx] = (Q_x_L[d_num_species + 1][idx] -
                double(1)/double(2)*Q_x_L[d_num_species][idx]*Q_x_L[d_num_species][idx]/rho_x_L[idx])/
                rho_x_L[idx];
        }
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
        {
            // Compute the linear index.
            const int idx = i + num_ghosts_0_conservative_variables;
            
            epsilon_x_R[idx] = (Q_x_R[d_num_species + 1][idx] -
                double(1)/double(2)*Q_x_R[d_num_species][idx]*Q_x_R[d_num_species][idx]/rho_x_R[idx])/
                rho_x_R[idx];
        }
        
        /*
         * Compute the mass fraction fields.
         */
        
        for (int si = 0; si < d_num_species; si++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
            {
                // Compute the linear index.
                const int idx = i + num_ghosts_0_conservative_variables;
                
                Y_x_L[si][idx] = Q_x_L[si][idx]/rho_x_L[idx];
                Y_x_R[si][idx] = Q_x_R[si][idx]/rho_x_R[idx];
            }
        }
        
        /*
         * Compute the volume fraction fields.
         */
        
        for (int si = 0; si < d_num_species - 1; si++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
            {
                // Compute the linear index.
                const int idx = i + num_ghosts_0_conservative_variables;
                
                Z_x_L[si][idx] = Q_x_L[d_num_species + 2 + si][idx];
                Z_x_R[si][idx] = Q_x_R[d_num_species + 2 + si][idx];
            }
        }
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computePressure(
                pressure_x_L,
                density_x_L,
                internal_energy_x_L,
                mass_fractions_x_L,
                volume_fractions_x_L,
                0,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computePressure(
                pressure_x_R,
                density_x_R,
                internal_energy_x_R,
                mass_fractions_x_R,
                volume_fractions_x_R,
                0,
                domain);
        
        /*
         * Compute the sound speed field.
         */
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computeGruneisenParameter(
                gruneisen_parameter_x_L,
                density_x_L,
                pressure_x_L,
                mass_fractions_x_L,
                volume_fractions_x_L,
                0,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computeGruneisenParameter(
                gruneisen_parameter_x_R,
                density_x_R,
                pressure_x_R,
                mass_fractions_x_R,
                volume_fractions_x_R,
                0,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computePressureDerivativeWithPartialDensities(
                partial_pressure_partial_partial_densities_x_L,
                density_x_L,
                pressure_x_L,
                mass_fractions_x_L,
                volume_fractions_x_L,
                0,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computePressureDerivativeWithPartialDensities(
                partial_pressure_partial_partial_densities_x_R,
                density_x_R,
                pressure_x_R,
                mass_fractions_x_R,
                volume_fractions_x_R,
                0,
                domain);
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
        {
            // Compute the linear index.
            const int idx = i + num_ghosts_0_conservative_variables;
            
            c_x_L[idx] = Gamma_x_L[idx]*p_x_L[idx]/rho_x_L[idx];
            c_x_R[idx] = Gamma_x_R[idx]*p_x_R[idx]/rho_x_R[idx];
        }
        
        for (int si = 0; si < d_num_species; si++)
        {   
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
            {
                // Compute the linear index.
                const int idx = i + num_ghosts_0_conservative_variables;
                
                c_x_L[idx] += Y_x_L[si][idx]*Psi_x_L[si][idx];
                c_x_R[idx] += Y_x_R[si][idx]*Psi_x_R[si][idx];
            }
        }
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
        {
            // Compute the linear index.
            const int idx = i + num_ghosts_0_conservative_variables;
            
            c_x_L[idx] = sqrt(c_x_L[idx]);
            c_x_R[idx] = sqrt(c_x_R[idx]);
        }
        
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
            
            double* u = velocity->getPointer(0, 0);
            
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
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
                    rho_x_L,
                    rho_x_R,
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
                    idx,
                    d_num_species,
                    num_eqn);
                
                if (s_x_star > double(0))
                {
                    u[idx_velocity] = u_x_L + s_x_minus*(Chi_x_star_LR - double(1));
                }
                else
                {
                    u[idx_velocity] = u_x_R + s_x_plus*(Chi_x_star_LR - double(1));
                }
            }
        }
        else
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
            {
                // Compute the linear indices.
                const int idx_flux = i + num_ghosts_0_convective_flux;
                const int idx = i + num_ghosts_0_conservative_variables;
                
                computeLocalConvectiveFluxInXDirectionFromConservativeVariablesHLLC_HLL1D(
                    F_x.data(),
                    Q_x_L.data(),
                    Q_x_R.data(),
                    rho_x_L,
                    rho_x_R,
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
                    idx,
                    d_num_species,
                    num_eqn);
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
         * Compute the density field.
         */
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
            {
                // Compute the linear index.
                const int idx = (i + num_ghosts_0_conservative_variables) +
                    (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables;
                
                rho_x_L[idx] = double(0);
            }
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_conservative_variables) +
                        (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables;
                    
                    rho_x_L[idx] += Q_x_L[si][idx];
                }
            }
        }
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
            {
                // Compute the linear index.
                const int idx = (i + num_ghosts_0_conservative_variables) +
                    (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables;
                
                rho_x_R[idx] = double(0);
            }
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_conservative_variables) +
                        (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables;
                    
                    rho_x_R[idx] += Q_x_R[si][idx];
                }
            }
        }
        
        /*
         * Compute the internal energy field.
         */
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
            {
                // Compute the linear index.
                const int idx = (i + num_ghosts_0_conservative_variables) +
                    (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables;
                
                epsilon_x_L[idx] = (Q_x_L[d_num_species + 2][idx] -
                    double(1)/double(2)*(Q_x_L[d_num_species][idx]*Q_x_L[d_num_species][idx] +
                    Q_x_L[d_num_species + 1][idx]*Q_x_L[d_num_species + 1][idx])/
                    rho_x_L[idx])/rho_x_L[idx];
            }
        }
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
            {
                // Compute the linear index.
                const int idx = (i + num_ghosts_0_conservative_variables) +
                    (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables;
                
                epsilon_x_R[idx] = (Q_x_R[d_num_species + 2][idx] -
                    double(1)/double(2)*(Q_x_R[d_num_species][idx]*Q_x_R[d_num_species][idx] +
                    Q_x_R[d_num_species + 1][idx]*Q_x_R[d_num_species + 1][idx])/
                    rho_x_R[idx])/rho_x_R[idx];
            }
        }
        
        /*
         * Compute the mass fraction fields.
         */
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_conservative_variables) +
                        (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables;
                    
                    Y_x_L[si][idx] = Q_x_L[si][idx]/rho_x_L[idx];
                    Y_x_R[si][idx] = Q_x_R[si][idx]/rho_x_R[idx];
                }
            }
        }
        
        /*
         * Compute the volume fraction fields.
         */
        
        for (int si = 0; si < d_num_species - 1; si++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_conservative_variables) +
                        (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables;
                    
                    Z_x_L[si][idx] = Q_x_L[d_num_species + 3 + si][idx];
                    Z_x_R[si][idx] = Q_x_R[d_num_species + 3 + si][idx];
                }
            }
        }
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computePressure(
                pressure_x_L,
                density_x_L,
                internal_energy_x_L,
                mass_fractions_x_L,
                volume_fractions_x_L,
                0,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computePressure(
                pressure_x_R,
                density_x_R,
                internal_energy_x_R,
                mass_fractions_x_R,
                volume_fractions_x_R,
                0,
                domain);
        
        /*
         * Compute the sound speed field.
         */
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computeGruneisenParameter(
                gruneisen_parameter_x_L,
                density_x_L,
                pressure_x_L,
                mass_fractions_x_L,
                volume_fractions_x_L,
                0,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computeGruneisenParameter(
                gruneisen_parameter_x_R,
                density_x_R,
                pressure_x_R,
                mass_fractions_x_R,
                volume_fractions_x_R,
                0,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computePressureDerivativeWithPartialDensities(
                partial_pressure_partial_partial_densities_x_L,
                density_x_L,
                pressure_x_L,
                mass_fractions_x_L,
                volume_fractions_x_L,
                0,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computePressureDerivativeWithPartialDensities(
                partial_pressure_partial_partial_densities_x_R,
                density_x_R,
                pressure_x_R,
                mass_fractions_x_R,
                volume_fractions_x_R,
                0,
                domain);
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
            {
                // Compute the linear index.
                const int idx = (i + num_ghosts_0_conservative_variables) +
                    (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables;
                
                c_x_L[idx] = Gamma_x_L[idx]*p_x_L[idx]/rho_x_L[idx];
                c_x_R[idx] = Gamma_x_R[idx]*p_x_R[idx]/rho_x_R[idx];
            }
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_conservative_variables) +
                        (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables;
                    
                    c_x_L[idx] += Y_x_L[si][idx]*Psi_x_L[si][idx];
                    c_x_R[idx] += Y_x_R[si][idx]*Psi_x_R[si][idx];
                }
            }
        }
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
            {
                // Compute the linear index.
                const int idx = (i + num_ghosts_0_conservative_variables) +
                    (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables;
                
                c_x_L[idx] = sqrt(c_x_L[idx]);
                c_x_R[idx] = sqrt(c_x_R[idx]);
            }
        }
        
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
            
            double* u = velocity->getPointer(0, 0);
            
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
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
                        rho_x_L,
                        rho_x_R,
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
                        idx,
                        d_num_species,
                        num_eqn);
                    
                    if (s_x_star > double(0))
                    {
                        u[idx_velocity] = u_x_L + s_x_minus*(Chi_x_star_LR - double(1));
                    }
                    else
                    {
                        u[idx_velocity] = u_x_R + s_x_plus*(Chi_x_star_LR - double(1));
                    }
                }
            }
        }
        else
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
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
                        rho_x_L,
                        rho_x_R,
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
                        idx,
                        d_num_species,
                        num_eqn);
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
         * Compute the density field.
         */
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_conservative_variables) +
                        (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables +
                        (k + num_ghosts_2_conservative_variables)*ghostcell_dim_0_conservative_variables*
                            ghostcell_dim_1_conservative_variables;
                    
                    rho_x_L[idx] = double(0);
                }
            }
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                    {
                        // Compute the linear index.
                        const int idx = (i + num_ghosts_0_conservative_variables) +
                            (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables +
                            (k + num_ghosts_2_conservative_variables)*ghostcell_dim_0_conservative_variables*
                                ghostcell_dim_1_conservative_variables;
                        
                        rho_x_L[idx] += Q_x_L[si][idx];
                    }
                }
            }
        }
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_conservative_variables) +
                        (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables +
                        (k + num_ghosts_2_conservative_variables)*ghostcell_dim_0_conservative_variables*
                            ghostcell_dim_1_conservative_variables;
                    
                    rho_x_R[idx] = double(0);
                }
            }
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                    {
                        // Compute the linear index.
                        const int idx = (i + num_ghosts_0_conservative_variables) +
                            (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables +
                            (k + num_ghosts_2_conservative_variables)*ghostcell_dim_0_conservative_variables*
                                ghostcell_dim_1_conservative_variables;
                        
                        rho_x_R[idx] += Q_x_R[si][idx];
                    }
                }
            }
        }
        
        /*
         * Compute the internal energy field.
         */
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_conservative_variables) +
                        (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables +
                        (k + num_ghosts_2_conservative_variables)*ghostcell_dim_0_conservative_variables*
                            ghostcell_dim_1_conservative_variables;
                    
                    epsilon_x_L[idx] = (Q_x_L[d_num_species + 3][idx] -
                        double(1)/double(2)*(Q_x_L[d_num_species][idx]*Q_x_L[d_num_species][idx] +
                        Q_x_L[d_num_species + 1][idx]*Q_x_L[d_num_species + 1][idx] +
                        Q_x_L[d_num_species + 2][idx]*Q_x_L[d_num_species + 2][idx])/
                        rho_x_L[idx])/rho_x_L[idx];
                }
            }
        }
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_conservative_variables) +
                        (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables +
                        (k + num_ghosts_2_conservative_variables)*ghostcell_dim_0_conservative_variables*
                            ghostcell_dim_1_conservative_variables;
                    
                    epsilon_x_R[idx] = (Q_x_R[d_num_species + 3][idx] -
                        double(1)/double(2)*(Q_x_R[d_num_species][idx]*Q_x_R[d_num_species][idx] +
                        Q_x_R[d_num_species + 1][idx]*Q_x_R[d_num_species + 1][idx] +
                        Q_x_R[d_num_species + 2][idx]*Q_x_R[d_num_species + 2][idx])/
                        rho_x_R[idx])/rho_x_R[idx];
                }
            }
        }
        
        /*
         * Compute the mass fraction fields.
         */
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                    {
                        // Compute the linear index.
                        const int idx = (i + num_ghosts_0_conservative_variables) +
                            (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables +
                            (k + num_ghosts_2_conservative_variables)*ghostcell_dim_0_conservative_variables*
                                ghostcell_dim_1_conservative_variables;
                        
                        Y_x_L[si][idx] = Q_x_L[si][idx]/rho_x_L[idx];
                        Y_x_R[si][idx] = Q_x_R[si][idx]/rho_x_R[idx];
                    }
                }
            }
        }
        
        /*
         * Compute the volume fraction fields.
         */
        
        for (int si = 0; si < d_num_species - 1; si++)
        {
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                    {
                        // Compute the linear index.
                        const int idx = (i + num_ghosts_0_conservative_variables) +
                            (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables +
                            (k + num_ghosts_2_conservative_variables)*ghostcell_dim_0_conservative_variables*
                                ghostcell_dim_1_conservative_variables;
                        
                        Z_x_L[si][idx] = Q_x_L[d_num_species + 4 + si][idx];
                        Z_x_R[si][idx] = Q_x_R[d_num_species + 4 + si][idx];
                    }
                }
            }
        }
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computePressure(
                pressure_x_L,
                density_x_L,
                internal_energy_x_L,
                mass_fractions_x_L,
                volume_fractions_x_L,
                0,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computePressure(
                pressure_x_R,
                density_x_R,
                internal_energy_x_R,
                mass_fractions_x_R,
                volume_fractions_x_R,
                0,
                domain);
        
        /*
         * Compute the sound speed field.
         */
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computeGruneisenParameter(
                gruneisen_parameter_x_L,
                density_x_L,
                pressure_x_L,
                mass_fractions_x_L,
                volume_fractions_x_L,
                0,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computeGruneisenParameter(
                gruneisen_parameter_x_R,
                density_x_R,
                pressure_x_R,
                mass_fractions_x_R,
                volume_fractions_x_R,
                0,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computePressureDerivativeWithPartialDensities(
                partial_pressure_partial_partial_densities_x_L,
                density_x_L,
                pressure_x_L,
                mass_fractions_x_L,
                volume_fractions_x_L,
                0,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computePressureDerivativeWithPartialDensities(
                partial_pressure_partial_partial_densities_x_R,
                density_x_R,
                pressure_x_R,
                mass_fractions_x_R,
                volume_fractions_x_R,
                0,
                domain);
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_conservative_variables) +
                        (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables +
                        (k + num_ghosts_2_conservative_variables)*ghostcell_dim_0_conservative_variables*
                            ghostcell_dim_1_conservative_variables;
                    
                    c_x_L[idx] = Gamma_x_L[idx]*p_x_L[idx]/rho_x_L[idx];
                    c_x_R[idx] = Gamma_x_R[idx]*p_x_R[idx]/rho_x_R[idx];
                }
            }
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                    {
                        // Compute the linear index.
                        const int idx = (i + num_ghosts_0_conservative_variables) +
                            (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables +
                            (k + num_ghosts_2_conservative_variables)*ghostcell_dim_0_conservative_variables*
                                ghostcell_dim_1_conservative_variables;
                        
                        c_x_L[idx] += Y_x_L[si][idx]*Psi_x_L[si][idx];
                        c_x_R[idx] += Y_x_R[si][idx]*Psi_x_R[si][idx];
                    }
                }
            }
        }
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_conservative_variables) +
                        (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables +
                        (k + num_ghosts_2_conservative_variables)*ghostcell_dim_0_conservative_variables*
                            ghostcell_dim_1_conservative_variables;
                    
                    c_x_L[idx] = sqrt(c_x_L[idx]);
                    c_x_R[idx] = sqrt(c_x_R[idx]);
                }
            }
        }
        
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
            
            double* u = velocity->getPointer(0, 0);
            
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
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
                            rho_x_L,
                            rho_x_R,
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
                            idx,
                            d_num_species,
                            num_eqn);
                        
                        if (s_x_star > double(0))
                        {
                            u[idx_velocity] = u_x_L + s_x_minus*(Chi_x_star_LR - double(1));
                        }
                        else
                        {
                            u[idx_velocity] = u_x_R + s_x_plus*(Chi_x_star_LR - double(1));
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
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
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
                            rho_x_L,
                            rho_x_R,
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
                            idx,
                            d_num_species,
                            num_eqn);
                    }
                }
            }
        }
    }
}


/*
 * Compute the convective flux and velocity in the y-direction from conservative variables with
 * HLLC Riemann solver.
 */
void
FlowModelRiemannSolverFiveEqnAllaire::computeConvectiveFluxAndVelocityInYDirectionFromConservativeVariablesHLLC_HLL(
    boost::shared_ptr<pdat::SideData<double> > convective_flux,
    boost::shared_ptr<pdat::SideData<double> > velocity,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables_B,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables_T,
    const hier::Box& domain,
    bool compute_velocity) const
{
    boost::shared_ptr<FlowModel> flow_model_tmp = d_flow_model.lock();
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
    
    const boost::shared_ptr<EquationOfStateMixingRules> equation_of_state_mixing_rules =
        flow_model_tmp->getEquationOfStateMixingRules();
    
    /*
     * Get the pointers to the side data of convective flux and conservative variables.
     */
    
    std::vector<double*> F_y;
    F_y.reserve(num_eqn);
    for (int ei = 0; ei < num_eqn; ei++)
    {
        F_y.push_back(convective_flux->getPointer(1, ei));
    }
    
    std::vector<double*> Q_y_B;
    std::vector<double*> Q_y_T;
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
    
    boost::shared_ptr<pdat::SideData<double> > density_y_B(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_variables,
            direction_y));
    
    boost::shared_ptr<pdat::SideData<double> > density_y_T(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_variables,
            direction_y));
    
    boost::shared_ptr<pdat::SideData<double> > internal_energy_y_B(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_variables,
            direction_y));
    
    boost::shared_ptr<pdat::SideData<double> > internal_energy_y_T(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_variables,
            direction_y));
    
    boost::shared_ptr<pdat::SideData<double> > mass_fractions_y_B(
        new pdat::SideData<double>(interior_box, d_num_species, num_ghosts_conservative_variables,
            direction_y));
    
    boost::shared_ptr<pdat::SideData<double> > mass_fractions_y_T(
        new pdat::SideData<double>(interior_box, d_num_species, num_ghosts_conservative_variables,
            direction_y));
    
    boost::shared_ptr<pdat::SideData<double> > volume_fractions_y_B(
        new pdat::SideData<double>(interior_box, d_num_species - 1, num_ghosts_conservative_variables,
            direction_y));
    
    boost::shared_ptr<pdat::SideData<double> > volume_fractions_y_T(
        new pdat::SideData<double>(interior_box, d_num_species - 1, num_ghosts_conservative_variables,
            direction_y));
    
    boost::shared_ptr<pdat::SideData<double> > pressure_y_B(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_variables,
            direction_y));
    
    boost::shared_ptr<pdat::SideData<double> > pressure_y_T(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_variables,
            direction_y));
    
    boost::shared_ptr<pdat::SideData<double> > gruneisen_parameter_y_B(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_variables,
            direction_y));
    
    boost::shared_ptr<pdat::SideData<double> > gruneisen_parameter_y_T(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_variables,
            direction_y));
    
    boost::shared_ptr<pdat::SideData<double> > partial_pressure_partial_partial_densities_y_B(
            new pdat::SideData<double>(interior_box, d_num_species, num_ghosts_conservative_variables,
            direction_y));
    
    boost::shared_ptr<pdat::SideData<double> > partial_pressure_partial_partial_densities_y_T(
            new pdat::SideData<double>(interior_box, d_num_species, num_ghosts_conservative_variables,
            direction_y));
    
    boost::shared_ptr<pdat::SideData<double> > sound_speed_y_B(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_variables,
            direction_y));
    
    boost::shared_ptr<pdat::SideData<double> > sound_speed_y_T(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_variables,
            direction_y));
    
    /*
     * Get the pointers to the temporary data.
     */
    
    double* rho_y_B = density_y_B->getPointer(1, 0);
    double* rho_y_T = density_y_T->getPointer(1, 0);
    
    double* epsilon_y_B = internal_energy_y_B->getPointer(1, 0);
    double* epsilon_y_T = internal_energy_y_T->getPointer(1, 0);
    
    std::vector<double*> Y_y_B;
    std::vector<double*> Y_y_T;
    Y_y_B.reserve(d_num_species);
    Y_y_T.reserve(d_num_species);
    for (int si = 0; si < d_num_species; si++)
    {
        Y_y_B.push_back(mass_fractions_y_B->getPointer(1, si));
        Y_y_T.push_back(mass_fractions_y_T->getPointer(1, si));
    }
    
    std::vector<double*> Z_y_B;
    std::vector<double*> Z_y_T;
    Z_y_B.reserve(d_num_species);
    Z_y_T.reserve(d_num_species);
    for (int si = 0; si < d_num_species - 1; si++)
    {
        Z_y_B.push_back(volume_fractions_y_B->getPointer(1, si));
        Z_y_T.push_back(volume_fractions_y_T->getPointer(1, si));
    }
    
    double* p_y_B = pressure_y_B->getPointer(1, 0);
    double* p_y_T = pressure_y_T->getPointer(1, 0);
    
    double* Gamma_y_B = gruneisen_parameter_y_B->getPointer(1, 0);
    double* Gamma_y_T = gruneisen_parameter_y_T->getPointer(1, 0);
    
    std::vector<double*> Psi_y_B;
    std::vector<double*> Psi_y_T;
    Psi_y_B.reserve(d_num_species);
    Psi_y_T.reserve(d_num_species);
    for (int si = 0; si < d_num_species; si++)
    {
        Psi_y_B.push_back(partial_pressure_partial_partial_densities_y_B->getPointer(1, si));
        Psi_y_T.push_back(partial_pressure_partial_partial_densities_y_T->getPointer(1, si));
    }
    
    double* c_y_B = sound_speed_y_B->getPointer(1, 0);
    double* c_y_T = sound_speed_y_T->getPointer(1, 0);
    
    double v_y_B = double(0);
    double v_y_T = double(0);
    
    double s_y_minus = double(0);
    double s_y_plus  = double(0);
    double s_y_star  = double(0);
    
    double Chi_y_star_BT = double(0);
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelRiemannSolverFiveEqnAllaire::"
            << "computeConvectiveFluxAndVelocityInYDirectionFromConservativeVariablesHLLC_HLL()\n"
            << "There is no y direction for 1D problem."
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
         * Compute the density field.
         */
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear index.
                const int idx = (i + num_ghosts_0_conservative_variables) +
                    (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables;
                
                rho_y_B[idx] = double(0);
            }
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_conservative_variables) +
                        (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables;
                    
                    rho_y_B[idx] += Q_y_B[si][idx];
                }
            }
        }
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear index.
                const int idx = (i + num_ghosts_0_conservative_variables) +
                    (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables;
                
                rho_y_T[idx] = double(0);
            }
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_conservative_variables) +
                        (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables;
                    
                    rho_y_T[idx] += Q_y_T[si][idx];
                }
            }
        }
        
        /*
         * Compute the internal energy field.
         */
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear index.
                const int idx = (i + num_ghosts_0_conservative_variables) +
                    (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables;
                
                epsilon_y_B[idx] = (Q_y_B[d_num_species + 2][idx] -
                    double(1)/double(2)*(Q_y_B[d_num_species][idx]*Q_y_B[d_num_species][idx] +
                    Q_y_B[d_num_species + 1][idx]*Q_y_B[d_num_species + 1][idx])/
                    rho_y_B[idx])/rho_y_B[idx];
            }
        }
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear index.
                const int idx = (i + num_ghosts_0_conservative_variables) +
                    (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables;
                
                epsilon_y_T[idx] = (Q_y_T[d_num_species + 2][idx] -
                    double(1)/double(2)*(Q_y_T[d_num_species][idx]*Q_y_T[d_num_species][idx] +
                    Q_y_T[d_num_species + 1][idx]*Q_y_T[d_num_species + 1][idx])/
                    rho_y_T[idx])/rho_y_T[idx];
            }
        }
        
        /*
         * Compute the mass fraction fields.
         */
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_conservative_variables) +
                        (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables;
                    
                    Y_y_B[si][idx] = Q_y_B[si][idx]/rho_y_B[idx];
                    Y_y_T[si][idx] = Q_y_T[si][idx]/rho_y_T[idx];
                }
            }
        }
        
        /*
         * Compute the volume fraction fields.
         */
        
        for (int si = 0; si < d_num_species - 1; si++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_conservative_variables) +
                        (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables;
                    
                    Z_y_B[si][idx] = Q_y_B[d_num_species + 3 + si][idx];
                    Z_y_T[si][idx] = Q_y_T[d_num_species + 3 + si][idx];
                }
            }
        }
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computePressure(
                pressure_y_B,
                density_y_B,
                internal_energy_y_B,
                mass_fractions_y_B,
                volume_fractions_y_B,
                1,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computePressure(
                pressure_y_T,
                density_y_T,
                internal_energy_y_T,
                mass_fractions_y_T,
                volume_fractions_y_T,
                1,
                domain);
        
        /*
         * Compute the sound speed field.
         */
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computeGruneisenParameter(
                gruneisen_parameter_y_B,
                density_y_B,
                pressure_y_B,
                mass_fractions_y_B,
                volume_fractions_y_B,
                1,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computeGruneisenParameter(
                gruneisen_parameter_y_T,
                density_y_T,
                pressure_y_T,
                mass_fractions_y_T,
                volume_fractions_y_T,
                1,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computePressureDerivativeWithPartialDensities(
                partial_pressure_partial_partial_densities_y_B,
                density_y_B,
                pressure_y_B,
                mass_fractions_y_B,
                volume_fractions_y_B,
                1,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computePressureDerivativeWithPartialDensities(
                partial_pressure_partial_partial_densities_y_T,
                density_y_T,
                pressure_y_T,
                mass_fractions_y_T,
                volume_fractions_y_T,
                1,
                domain);
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear index.
                const int idx = (i + num_ghosts_0_conservative_variables) +
                    (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables;
                
                c_y_B[idx] = Gamma_y_B[idx]*p_y_B[idx]/rho_y_B[idx];
                c_y_T[idx] = Gamma_y_T[idx]*p_y_T[idx]/rho_y_T[idx];
            }
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_conservative_variables) +
                        (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables;
                    
                    c_y_B[idx] += Y_y_B[si][idx]*Psi_y_B[si][idx];
                    c_y_T[idx] += Y_y_T[si][idx]*Psi_y_T[si][idx];
                }
            }
        }
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear index.
                const int idx = (i + num_ghosts_0_conservative_variables) +
                    (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables;
                
                c_y_B[idx] = sqrt(c_y_B[idx]);
                c_y_T[idx] = sqrt(c_y_T[idx]);
            }
        }
        
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
            
            double* v = velocity->getPointer(1, 1);
            
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
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
                        rho_y_B,
                        rho_y_T,
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
                        idx,
                        d_num_species,
                        num_eqn);
                    
                    if (s_y_star > double(0))
                    {
                        v[idx_velocity] = v_y_B + s_y_minus*(Chi_y_star_BT - double(1));
                    }
                    else
                    {
                        v[idx_velocity] = v_y_T + s_y_plus*(Chi_y_star_BT - double(1));
                    }
                }
            }
        }
        else
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
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
                        rho_y_B,
                        rho_y_T,
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
                        idx,
                        d_num_species,
                        num_eqn);
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
         * Compute the density field.
         */
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_conservative_variables) +
                        (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables +
                        (k + num_ghosts_2_conservative_variables)*ghostcell_dim_0_conservative_variables*
                            ghostcell_dim_1_conservative_variables;
                    
                    rho_y_B[idx] = double(0);
                }
            }
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx = (i + num_ghosts_0_conservative_variables) +
                            (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables +
                            (k + num_ghosts_2_conservative_variables)*ghostcell_dim_0_conservative_variables*
                                ghostcell_dim_1_conservative_variables;
                        
                        rho_y_B[idx] += Q_y_B[si][idx];
                    }
                }
            }
        }
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_conservative_variables) +
                        (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables +
                        (k + num_ghosts_2_conservative_variables)*ghostcell_dim_0_conservative_variables*
                            ghostcell_dim_1_conservative_variables;
                    
                    rho_y_T[idx] = double(0);
                }
            }
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx = (i + num_ghosts_0_conservative_variables) +
                            (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables +
                            (k + num_ghosts_2_conservative_variables)*ghostcell_dim_0_conservative_variables*
                                ghostcell_dim_1_conservative_variables;
                        
                        rho_y_T[idx] += Q_y_T[si][idx];
                    }
                }
            }
        }
        
        /*
         * Compute the internal energy field.
         */
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_conservative_variables) +
                        (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables +
                        (k + num_ghosts_2_conservative_variables)*ghostcell_dim_0_conservative_variables*
                            ghostcell_dim_1_conservative_variables;
                    
                    epsilon_y_B[idx] = (Q_y_B[d_num_species + 3][idx] -
                        double(1)/double(2)*(Q_y_B[d_num_species][idx]*Q_y_B[d_num_species][idx] +
                        Q_y_B[d_num_species + 1][idx]*Q_y_B[d_num_species + 1][idx] +
                        Q_y_B[d_num_species + 2][idx]*Q_y_B[d_num_species + 2][idx])/
                        rho_y_B[idx])/rho_y_B[idx];
                }
            }
        }
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_conservative_variables) +
                        (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables +
                        (k + num_ghosts_2_conservative_variables)*ghostcell_dim_0_conservative_variables*
                            ghostcell_dim_1_conservative_variables;
                    
                    epsilon_y_T[idx] = (Q_y_T[d_num_species + 3][idx] -
                        double(1)/double(2)*(Q_y_T[d_num_species][idx]*Q_y_T[d_num_species][idx] +
                        Q_y_T[d_num_species + 1][idx]*Q_y_T[d_num_species + 1][idx] +
                        Q_y_T[d_num_species + 2][idx]*Q_y_T[d_num_species + 2][idx])/
                        rho_y_T[idx])/rho_y_T[idx];
                }
            }
        }
        
        /*
         * Compute the mass fraction fields.
         */
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx = (i + num_ghosts_0_conservative_variables) +
                            (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables +
                            (k + num_ghosts_2_conservative_variables)*ghostcell_dim_0_conservative_variables*
                                ghostcell_dim_1_conservative_variables;
                        
                        Y_y_B[si][idx] = Q_y_B[si][idx]/rho_y_B[idx];
                        Y_y_T[si][idx] = Q_y_T[si][idx]/rho_y_T[idx];
                    }
                }
            }
        }
        
        /*
         * Compute the volume fraction fields.
         */
        
        for (int si = 0; si < d_num_species - 1; si++)
        {
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx = (i + num_ghosts_0_conservative_variables) +
                            (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables +
                            (k + num_ghosts_2_conservative_variables)*ghostcell_dim_0_conservative_variables*
                                ghostcell_dim_1_conservative_variables;
                        
                        Z_y_B[si][idx] = Q_y_B[d_num_species + 4 + si][idx];
                        Z_y_T[si][idx] = Q_y_T[d_num_species + 4 + si][idx];
                    }
                }
            }
        }
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computePressure(
                pressure_y_B,
                density_y_B,
                internal_energy_y_B,
                mass_fractions_y_B,
                volume_fractions_y_B,
                1,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computePressure(
                pressure_y_T,
                density_y_T,
                internal_energy_y_T,
                mass_fractions_y_T,
                volume_fractions_y_T,
                1,
                domain);
        
        /*
         * Compute the sound speed field.
         */
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computeGruneisenParameter(
                gruneisen_parameter_y_B,
                density_y_B,
                pressure_y_B,
                mass_fractions_y_B,
                volume_fractions_y_B,
                1,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computeGruneisenParameter(
                gruneisen_parameter_y_T,
                density_y_T,
                pressure_y_T,
                mass_fractions_y_T,
                volume_fractions_y_T,
                1,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computePressureDerivativeWithPartialDensities(
                partial_pressure_partial_partial_densities_y_B,
                density_y_B,
                pressure_y_B,
                mass_fractions_y_B,
                volume_fractions_y_B,
                1,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computePressureDerivativeWithPartialDensities(
                partial_pressure_partial_partial_densities_y_T,
                density_y_T,
                pressure_y_T,
                mass_fractions_y_T,
                volume_fractions_y_T,
                1,
                domain);
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_conservative_variables) +
                        (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables +
                        (k + num_ghosts_2_conservative_variables)*ghostcell_dim_0_conservative_variables*
                            ghostcell_dim_1_conservative_variables;
                    
                    c_y_B[idx] = Gamma_y_B[idx]*p_y_B[idx]/rho_y_B[idx];
                    c_y_T[idx] = Gamma_y_T[idx]*p_y_T[idx]/rho_y_T[idx];
                }
            }
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx = (i + num_ghosts_0_conservative_variables) +
                            (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables +
                            (k + num_ghosts_2_conservative_variables)*ghostcell_dim_0_conservative_variables*
                                ghostcell_dim_1_conservative_variables;
                        
                        c_y_B[idx] += Y_y_B[si][idx]*Psi_y_B[si][idx];
                        c_y_T[idx] += Y_y_T[si][idx]*Psi_y_T[si][idx];
                    }
                }
            }
        }
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_conservative_variables) +
                        (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables +
                        (k + num_ghosts_2_conservative_variables)*ghostcell_dim_0_conservative_variables*
                            ghostcell_dim_1_conservative_variables;
                    
                    c_y_B[idx] = sqrt(c_y_B[idx]);
                    c_y_T[idx] = sqrt(c_y_T[idx]);
                }
            }
        }
        
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
            
            double* v = velocity->getPointer(1, 1);
            
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
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
                            rho_y_B,
                            rho_y_T,
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
                            idx,
                            d_num_species,
                            num_eqn);

                        if (s_y_star > double(0))
                        {
                            v[idx_velocity] = v_y_B + s_y_minus*(Chi_y_star_BT - double(1));
                        }
                        else
                        {
                            v[idx_velocity] = v_y_T + s_y_plus*(Chi_y_star_BT - double(1));
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
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
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
                            rho_y_B,
                            rho_y_T,
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
                            idx,
                            d_num_species,
                            num_eqn);
                    }
                }
            }
        }
    }
}


/*
 * Compute the convective flux and velocity in the z-direction from conservative variables with
 * HLLC Riemann solver.
 */
void
FlowModelRiemannSolverFiveEqnAllaire::computeConvectiveFluxAndVelocityInZDirectionFromConservativeVariablesHLLC_HLL(
    boost::shared_ptr<pdat::SideData<double> > convective_flux,
    boost::shared_ptr<pdat::SideData<double> > velocity,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables_B,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables_F,
    const hier::Box& domain,
    bool compute_velocity) const
{
    boost::shared_ptr<FlowModel> flow_model_tmp = d_flow_model.lock();
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
    
    const boost::shared_ptr<EquationOfStateMixingRules> equation_of_state_mixing_rules =
        flow_model_tmp->getEquationOfStateMixingRules();
    
    /*
     * Get the pointers to the side data of convective flux and conservative variables.
     */
    
    std::vector<double*> F_z;
    F_z.reserve(num_eqn);
    for (int ei = 0; ei < num_eqn; ei++)
    {
        F_z.push_back(convective_flux->getPointer(2, ei));
    }
    
    std::vector<double*> Q_z_B;
    std::vector<double*> Q_z_F;
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
    
    boost::shared_ptr<pdat::SideData<double> > density_z_B(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_variables,
            direction_z));
    
    boost::shared_ptr<pdat::SideData<double> > density_z_F(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_variables,
            direction_z));
    
    boost::shared_ptr<pdat::SideData<double> > internal_energy_z_B(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_variables,
            direction_z));
    
    boost::shared_ptr<pdat::SideData<double> > internal_energy_z_F(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_variables,
            direction_z));
    
    boost::shared_ptr<pdat::SideData<double> > mass_fractions_z_B(
        new pdat::SideData<double>(interior_box, d_num_species, num_ghosts_conservative_variables,
            direction_z));
    
    boost::shared_ptr<pdat::SideData<double> > mass_fractions_z_F(
        new pdat::SideData<double>(interior_box, d_num_species, num_ghosts_conservative_variables,
            direction_z));
    
    boost::shared_ptr<pdat::SideData<double> > volume_fractions_z_B(
        new pdat::SideData<double>(interior_box, d_num_species - 1, num_ghosts_conservative_variables,
            direction_z));
    
    boost::shared_ptr<pdat::SideData<double> > volume_fractions_z_F(
        new pdat::SideData<double>(interior_box, d_num_species - 1, num_ghosts_conservative_variables,
            direction_z));
    
    boost::shared_ptr<pdat::SideData<double> > pressure_z_B(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_variables,
            direction_z));
    
    boost::shared_ptr<pdat::SideData<double> > pressure_z_F(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_variables,
            direction_z));
    
    boost::shared_ptr<pdat::SideData<double> > gruneisen_parameter_z_B(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_variables,
            direction_z));
    
    boost::shared_ptr<pdat::SideData<double> > gruneisen_parameter_z_F(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_variables,
            direction_z));
    
    boost::shared_ptr<pdat::SideData<double> > partial_pressure_partial_partial_densities_z_B(
            new pdat::SideData<double>(interior_box, d_num_species, num_ghosts_conservative_variables,
            direction_z));
    
    boost::shared_ptr<pdat::SideData<double> > partial_pressure_partial_partial_densities_z_F(
            new pdat::SideData<double>(interior_box, d_num_species, num_ghosts_conservative_variables,
            direction_z));
    
    boost::shared_ptr<pdat::SideData<double> > sound_speed_z_B(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_variables,
            direction_z));
    
    boost::shared_ptr<pdat::SideData<double> > sound_speed_z_F(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_variables,
            direction_z));
    
    /*
     * Get the pointers to the temporary data.
     */
    
    double* rho_z_B = density_z_B->getPointer(2, 0);
    double* rho_z_F = density_z_F->getPointer(2, 0);
    
    double* epsilon_z_B = internal_energy_z_B->getPointer(2, 0);
    double* epsilon_z_F = internal_energy_z_F->getPointer(2, 0);
    
    std::vector<double*> Y_z_B;
    std::vector<double*> Y_z_F;
    Y_z_B.reserve(d_num_species);
    Y_z_F.reserve(d_num_species);
    for (int si = 0; si < d_num_species; si++)
    {
        Y_z_B.push_back(mass_fractions_z_B->getPointer(2, si));
        Y_z_F.push_back(mass_fractions_z_F->getPointer(2, si));
    }
    
    std::vector<double*> Z_z_B;
    std::vector<double*> Z_z_F;
    Z_z_B.reserve(d_num_species);
    Z_z_F.reserve(d_num_species);
    for (int si = 0; si < d_num_species - 1; si++)
    {
        Z_z_B.push_back(volume_fractions_z_B->getPointer(2, si));
        Z_z_F.push_back(volume_fractions_z_F->getPointer(2, si));
    }
    
    double* p_z_B = pressure_z_B->getPointer(2, 0);
    double* p_z_F = pressure_z_F->getPointer(2, 0);
    
    double* Gamma_z_B = gruneisen_parameter_z_B->getPointer(2, 0);
    double* Gamma_z_F = gruneisen_parameter_z_F->getPointer(2, 0);
    
    std::vector<double*> Psi_z_B;
    std::vector<double*> Psi_z_F;
    Psi_z_B.reserve(d_num_species);
    Psi_z_F.reserve(d_num_species);
    for (int si = 0; si < d_num_species; si++)
    {
        Psi_z_B.push_back(partial_pressure_partial_partial_densities_z_B->getPointer(2, si));
        Psi_z_F.push_back(partial_pressure_partial_partial_densities_z_F->getPointer(2, si));
    }
    
    double* c_z_B = sound_speed_z_B->getPointer(2, 0);
    double* c_z_F = sound_speed_z_F->getPointer(2, 0);
    
    double w_z_B = double(0);
    double w_z_F = double(0);
    
    double s_z_minus = double(0);
    double s_z_plus  = double(0);
    double s_z_star  = double(0);
    
    double Chi_z_star_BF = double(0);
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelRiemannSolverFiveEqnAllaire::"
            << "computeConvectiveFluxAndVelocityInZDirectionFromConservativeVariablesHLLC_HLL()\n"
            << "There is no z direction for 1D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelRiemannSolverFiveEqnAllaire::"
            << "computeConvectiveFluxAndVelocityInZDirectionFromConservativeVariablesHLLC_HLL()\n"
            << "There is no z direction for 2D problem."
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
         * Compute the density field.
         */
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2 + 1; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_conservative_variables) +
                        (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables +
                        (k + num_ghosts_2_conservative_variables)*ghostcell_dim_0_conservative_variables*
                            ghostcell_dim_1_conservative_variables;
                    
                    rho_z_B[idx] = double(0);
                }
            }
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2 + 1; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx = (i + num_ghosts_0_conservative_variables) +
                            (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables +
                            (k + num_ghosts_2_conservative_variables)*ghostcell_dim_0_conservative_variables*
                                ghostcell_dim_1_conservative_variables;
                        
                        rho_z_B[idx] += Q_z_B[si][idx];
                    }
                }
            }
        }
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2 + 1; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_conservative_variables) +
                        (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables +
                        (k + num_ghosts_2_conservative_variables)*ghostcell_dim_0_conservative_variables*
                            ghostcell_dim_1_conservative_variables;
                    
                    rho_z_F[idx] = double(0);
                }
            }
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2 + 1; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx = (i + num_ghosts_0_conservative_variables) +
                            (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables +
                            (k + num_ghosts_2_conservative_variables)*ghostcell_dim_0_conservative_variables*
                                ghostcell_dim_1_conservative_variables;
                        
                        rho_z_F[idx] += Q_z_F[si][idx];
                    }
                }
            }
        }
        
        /*
         * Compute the internal energy field.
         */
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2 + 1; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_conservative_variables) +
                        (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables +
                        (k + num_ghosts_2_conservative_variables)*ghostcell_dim_0_conservative_variables*
                            ghostcell_dim_1_conservative_variables;
                    
                    epsilon_z_B[idx] = (Q_z_B[d_num_species + 3][idx] -
                        double(1)/double(2)*(Q_z_B[d_num_species][idx]*Q_z_B[d_num_species][idx] +
                        Q_z_B[d_num_species + 1][idx]*Q_z_B[d_num_species + 1][idx] +
                        Q_z_B[d_num_species + 2][idx]*Q_z_B[d_num_species + 2][idx])/
                        rho_z_B[idx])/rho_z_B[idx];
                }
            }
        }
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2 + 1; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_conservative_variables) +
                        (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables +
                        (k + num_ghosts_2_conservative_variables)*ghostcell_dim_0_conservative_variables*
                            ghostcell_dim_1_conservative_variables;
                    
                    epsilon_z_F[idx] = (Q_z_F[d_num_species + 3][idx] -
                        double(1)/double(2)*(Q_z_F[d_num_species][idx]*Q_z_F[d_num_species][idx] +
                        Q_z_F[d_num_species + 1][idx]*Q_z_F[d_num_species + 1][idx] +
                        Q_z_F[d_num_species + 2][idx]*Q_z_F[d_num_species + 2][idx])/
                        rho_z_F[idx])/rho_z_F[idx];
                }
            }
        }
        
        /*
         * Compute the mass fraction fields.
         */
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2 + 1; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx = (i + num_ghosts_0_conservative_variables) +
                            (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables +
                            (k + num_ghosts_2_conservative_variables)*ghostcell_dim_0_conservative_variables*
                                ghostcell_dim_1_conservative_variables;
                        
                        Y_z_B[si][idx] = Q_z_B[si][idx]/rho_z_B[idx];
                        Y_z_F[si][idx] = Q_z_F[si][idx]/rho_z_F[idx];
                    }
                }
            }
        }
        
        /*
         * Compute the volume fraction fields.
         */
        
        for (int si = 0; si < d_num_species - 1; si++)
        {
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2 + 1; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx = (i + num_ghosts_0_conservative_variables) +
                            (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables +
                            (k + num_ghosts_2_conservative_variables)*ghostcell_dim_0_conservative_variables*
                                ghostcell_dim_1_conservative_variables;
                        
                        Z_z_B[si][idx] = Q_z_B[d_num_species + 4 + si][idx];
                        Z_z_F[si][idx] = Q_z_F[d_num_species + 4 + si][idx];
                    }
                }
            }
        }
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computePressure(
                pressure_z_B,
                density_z_B,
                internal_energy_z_B,
                mass_fractions_z_B,
                volume_fractions_z_B,
                2,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computePressure(
                pressure_z_F,
                density_z_F,
                internal_energy_z_F,
                mass_fractions_z_F,
                volume_fractions_z_F,
                2,
                domain);
        
        /*
         * Compute the sound speed field.
         */
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computeGruneisenParameter(
                gruneisen_parameter_z_B,
                density_z_B,
                pressure_z_B,
                mass_fractions_z_B,
                volume_fractions_z_B,
                2,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computeGruneisenParameter(
                gruneisen_parameter_z_F,
                density_z_F,
                pressure_z_F,
                mass_fractions_z_F,
                volume_fractions_z_F,
                2,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computePressureDerivativeWithPartialDensities(
                partial_pressure_partial_partial_densities_z_B,
                density_z_B,
                pressure_z_B,
                mass_fractions_z_B,
                volume_fractions_z_B,
                2,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computePressureDerivativeWithPartialDensities(
                partial_pressure_partial_partial_densities_z_F,
                density_z_F,
                pressure_z_F,
                mass_fractions_z_F,
                volume_fractions_z_F,
                2,
                domain);
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2 + 1; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_conservative_variables) +
                        (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables +
                        (k + num_ghosts_2_conservative_variables)*ghostcell_dim_0_conservative_variables*
                            ghostcell_dim_1_conservative_variables;
                    
                    c_z_B[idx] = Gamma_z_B[idx]*p_z_B[idx]/rho_z_B[idx];
                    c_z_F[idx] = Gamma_z_F[idx]*p_z_F[idx]/rho_z_F[idx];
                }
            }
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2 + 1; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx = (i + num_ghosts_0_conservative_variables) +
                            (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables +
                            (k + num_ghosts_2_conservative_variables)*ghostcell_dim_0_conservative_variables*
                                ghostcell_dim_1_conservative_variables;
                        
                        c_z_B[idx] += Y_z_B[si][idx]*Psi_z_B[si][idx];
                        c_z_F[idx] += Y_z_F[si][idx]*Psi_z_F[si][idx];
                    }
                }
            }
        }
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2 + 1; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_conservative_variables) +
                        (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables +
                        (k + num_ghosts_2_conservative_variables)*ghostcell_dim_0_conservative_variables*
                            ghostcell_dim_1_conservative_variables;
                    
                    c_z_B[idx] = sqrt(c_z_B[idx]);
                    c_z_F[idx] = sqrt(c_z_F[idx]);
                }
            }
        }
        
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
            
            double* w = velocity->getPointer(2, 2);
            
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2 + 1; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
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
                            rho_z_B,
                            rho_z_F,
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
                            idx,
                            d_num_species,
                            num_eqn);
                        
                        if (s_z_star > double(0))
                        {
                            w[idx_velocity] = w_z_B + s_z_minus*(Chi_z_star_BF - double(1));
                        }
                        else
                        {
                            w[idx_velocity] = w_z_F + s_z_plus*(Chi_z_star_BF - double(1));
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
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
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
                            rho_z_B,
                            rho_z_F,
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
                            idx,
                            d_num_species,
                            num_eqn);
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
FlowModelRiemannSolverFiveEqnAllaire::computeConvectiveFluxAndVelocityInXDirectionFromPrimitiveVariablesHLLC_HLL(
    boost::shared_ptr<pdat::SideData<double> > convective_flux,
    boost::shared_ptr<pdat::SideData<double> > velocity,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables_L,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables_R,
    const hier::Box& domain,
    bool compute_velocity) const
{
    boost::shared_ptr<FlowModel> flow_model_tmp = d_flow_model.lock();
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
    
    const boost::shared_ptr<EquationOfStateMixingRules> equation_of_state_mixing_rules =
        flow_model_tmp->getEquationOfStateMixingRules();
    
    /*
     * Get the pointers to the side data of convective flux and primitive variables.
     */
    
    std::vector<double*> F_x;
    F_x.reserve(num_eqn);
    for (int ei = 0; ei < num_eqn; ei++)
    {
        F_x.push_back(convective_flux->getPointer(0, ei));
    }
    
    std::vector<double*> V_x_L;
    std::vector<double*> V_x_R;
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
    
    boost::shared_ptr<pdat::SideData<double> > density_x_L(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_primitive_variables,
            direction_x));
    
    boost::shared_ptr<pdat::SideData<double> > density_x_R(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_primitive_variables,
            direction_x));
    
    boost::shared_ptr<pdat::SideData<double> > mass_fractions_x_L(
        new pdat::SideData<double>(interior_box, d_num_species, num_ghosts_primitive_variables,
            direction_x));
    
    boost::shared_ptr<pdat::SideData<double> > mass_fractions_x_R(
        new pdat::SideData<double>(interior_box, d_num_species, num_ghosts_primitive_variables,
            direction_x));
    
    boost::shared_ptr<pdat::SideData<double> > volume_fractions_x_L(
        new pdat::SideData<double>(interior_box, d_num_species - 1, num_ghosts_primitive_variables,
            direction_x));
    
    boost::shared_ptr<pdat::SideData<double> > volume_fractions_x_R(
        new pdat::SideData<double>(interior_box, d_num_species - 1, num_ghosts_primitive_variables,
            direction_x));
    
    boost::shared_ptr<pdat::SideData<double> > gruneisen_parameter_x_L(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_primitive_variables,
            direction_x));
    
    boost::shared_ptr<pdat::SideData<double> > gruneisen_parameter_x_R(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_primitive_variables,
            direction_x));
    
    boost::shared_ptr<pdat::SideData<double> > partial_pressure_partial_partial_densities_x_L(
            new pdat::SideData<double>(interior_box, d_num_species, num_ghosts_primitive_variables,
            direction_x));
    
    boost::shared_ptr<pdat::SideData<double> > partial_pressure_partial_partial_densities_x_R(
            new pdat::SideData<double>(interior_box, d_num_species, num_ghosts_primitive_variables,
            direction_x));
    
    boost::shared_ptr<pdat::SideData<double> > sound_speed_x_L(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_primitive_variables,
            direction_x));
    
    boost::shared_ptr<pdat::SideData<double> > sound_speed_x_R(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_primitive_variables,
            direction_x));
    
    boost::shared_ptr<pdat::SideData<double> > internal_energy_x_L(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_primitive_variables,
            direction_x));
    
    boost::shared_ptr<pdat::SideData<double> > internal_energy_x_R(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_primitive_variables,
            direction_x));
    
    /*
     * Get the pointers to the temporary data.
     */
    
    double* rho_x_L = density_x_L->getPointer(0, 0);
    double* rho_x_R = density_x_R->getPointer(0, 0);
    
    std::vector<double*> Y_x_L;
    std::vector<double*> Y_x_R;
    Y_x_L.reserve(d_num_species);
    Y_x_R.reserve(d_num_species);
    for (int si = 0; si < d_num_species; si++)
    {
        Y_x_L.push_back(mass_fractions_x_L->getPointer(0, si));
        Y_x_R.push_back(mass_fractions_x_R->getPointer(0, si));
    }
    
    std::vector<double*> Z_x_L;
    std::vector<double*> Z_x_R;
    Z_x_L.reserve(d_num_species);
    Z_x_R.reserve(d_num_species);
    for (int si = 0; si < d_num_species - 1; si++)
    {
        Z_x_L.push_back(volume_fractions_x_L->getPointer(0, si));
        Z_x_R.push_back(volume_fractions_x_R->getPointer(0, si));
    }
    
    double* Gamma_x_L = gruneisen_parameter_x_L->getPointer(0, 0);
    double* Gamma_x_R = gruneisen_parameter_x_R->getPointer(0, 0);
    
    std::vector<double*> Psi_x_L;
    std::vector<double*> Psi_x_R;
    Psi_x_L.reserve(d_num_species);
    Psi_x_R.reserve(d_num_species);
    for (int si = 0; si < d_num_species; si++)
    {
        Psi_x_L.push_back(partial_pressure_partial_partial_densities_x_L->getPointer(0, si));
        Psi_x_R.push_back(partial_pressure_partial_partial_densities_x_R->getPointer(0, si));
    }
    
    double* c_x_L = sound_speed_x_L->getPointer(0, 0);
    double* c_x_R = sound_speed_x_R->getPointer(0, 0);
    
    double* epsilon_x_L = internal_energy_x_L->getPointer(0, 0);
    double* epsilon_x_R = internal_energy_x_R->getPointer(0, 0);
    
    double s_x_minus = double(0);
    double s_x_plus  = double(0);
    double s_x_star  = double(0);
    
    double Chi_x_star_LR = double(0);
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and numbers of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int num_ghosts_0_convective_flux = num_ghosts_convective_flux[0];
        const int num_ghosts_0_primitive_variables = num_ghosts_primitive_variables[0];
        
        /*
         * Compute the density field.
         */
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
        {
            // Compute the linear index.
            const int idx = i + num_ghosts_0_primitive_variables;
            
            rho_x_L[idx] = double(0);
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
            {
                // Compute the linear index.
                const int idx = i + num_ghosts_0_primitive_variables;
                
                rho_x_L[idx] += V_x_L[si][idx];
            }
        }
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
        {
            // Compute the linear index.
            const int idx = i + num_ghosts_0_primitive_variables;
            
            rho_x_R[idx] = double(0);
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
            {
                // Compute the linear index.
                const int idx = i + num_ghosts_0_primitive_variables;
                
                rho_x_R[idx] += V_x_R[si][idx];
            }
        }
        
        /*
         * Compute the mass fraction fields.
         */
        
        for (int si = 0; si < d_num_species; si++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
            {
                // Compute the linear index.
                const int idx = i + num_ghosts_0_primitive_variables;
                
                Y_x_L[si][idx] = V_x_L[si][idx]/rho_x_L[idx];
                Y_x_R[si][idx] = V_x_R[si][idx]/rho_x_R[idx];
            }
        }
        
        /*
         * Compute the volume fraction fields.
         */
        
        for (int si = 0; si < d_num_species - 1; si++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
            {
                // Compute the linear index.
                const int idx = i + num_ghosts_0_primitive_variables;
                
                Z_x_L[si][idx] = V_x_L[d_num_species + 2 + si][idx];
                Z_x_R[si][idx] = V_x_R[d_num_species + 2 + si][idx];
            }
        }
        
        /*
         * Compute the sound speed field.
         */
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computeGruneisenParameter(
                gruneisen_parameter_x_L,
                density_x_L,
                primitive_variables_L[d_num_species + 1],
                mass_fractions_x_L,
                volume_fractions_x_L,
                0,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computeGruneisenParameter(
                gruneisen_parameter_x_R,
                density_x_R,
                primitive_variables_R[d_num_species + 1],
                mass_fractions_x_R,
                volume_fractions_x_R,
                0,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computePressureDerivativeWithPartialDensities(
                partial_pressure_partial_partial_densities_x_L,
                density_x_L,
                primitive_variables_L[d_num_species + 1],
                mass_fractions_x_L,
                volume_fractions_x_L,
                0,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computePressureDerivativeWithPartialDensities(
                partial_pressure_partial_partial_densities_x_R,
                density_x_R,
                primitive_variables_R[d_num_species + 1],
                mass_fractions_x_R,
                volume_fractions_x_R,
                0,
                domain);
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
        {
            // Compute the linear index.
            const int idx = i + num_ghosts_0_primitive_variables;
            
            c_x_L[idx] = Gamma_x_L[idx]*V_x_L[d_num_species + 1][idx]/rho_x_L[idx];
            c_x_R[idx] = Gamma_x_R[idx]*V_x_R[d_num_species + 1][idx]/rho_x_R[idx];
        }
        
        for (int si = 0; si < d_num_species; si++)
        {   
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
            {
                // Compute the linear index.
                const int idx = i + num_ghosts_0_primitive_variables;
                
                c_x_L[idx] += Y_x_L[si][idx]*Psi_x_L[si][idx];
                c_x_R[idx] += Y_x_R[si][idx]*Psi_x_R[si][idx];
            }
        }
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
        {
            // Compute the linear index.
            const int idx = i + num_ghosts_0_primitive_variables;
            
            c_x_L[idx] = sqrt(c_x_L[idx]);
            c_x_R[idx] = sqrt(c_x_R[idx]);
        }
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computeInternalEnergy(
                internal_energy_x_L,
                density_x_L,
                primitive_variables_L[d_num_species + 1],
                mass_fractions_x_L,
                volume_fractions_x_L,
                0,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computeInternalEnergy(
                internal_energy_x_R,
                density_x_R,
                primitive_variables_R[d_num_species + 1],
                mass_fractions_x_R,
                volume_fractions_x_R,
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
            
            double* u = velocity->getPointer(0, 0);
            
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
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
                    rho_x_L,
                    rho_x_R,
                    c_x_L,
                    c_x_R,
                    epsilon_x_L,
                    epsilon_x_R,
                    s_x_minus,
                    s_x_plus,
                    s_x_star,
                    Chi_x_star_LR,
                    idx_flux,
                    idx,
                    d_num_species,
                    num_eqn);
                
                if (s_x_star > double(0))
                {
                    u[idx_velocity] = V_x_L[d_num_species][idx] + s_x_minus*(Chi_x_star_LR - double(1));
                }
                else
                {
                    u[idx_velocity] = V_x_R[d_num_species][idx] + s_x_plus*(Chi_x_star_LR - double(1));
                }
            }
        }
        else
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
            {
                // Compute the linear indices.
                const int idx_flux = i + num_ghosts_0_convective_flux;
                const int idx = i + num_ghosts_0_primitive_variables;
                
                computeLocalConvectiveFluxInXDirectionFromPrimitiveVariablesHLLC_HLL1D(
                    F_x.data(),
                    V_x_L.data(),
                    V_x_R.data(),
                    rho_x_L,
                    rho_x_R,
                    c_x_L,
                    c_x_R,
                    epsilon_x_L,
                    epsilon_x_R,
                    s_x_minus,
                    s_x_plus,
                    s_x_star,
                    Chi_x_star_LR,
                    idx_flux,
                    idx,
                    d_num_species,
                    num_eqn);
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
        
        /*
         * Compute the density field.
         */
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
            {
                // Compute the linear index.
                const int idx = (i + num_ghosts_0_primitive_variables) +
                    (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables;
                
                rho_x_L[idx] = double(0);
            }
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_primitive_variables) +
                        (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables;
                    
                    rho_x_L[idx] += V_x_L[si][idx];
                }
            }
        }
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
            {
                // Compute the linear index.
                const int idx = (i + num_ghosts_0_primitive_variables) +
                    (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables;
                
                rho_x_R[idx] = double(0);
            }
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_primitive_variables) +
                        (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables;
                    
                    rho_x_R[idx] += V_x_R[si][idx];
                }
            }
        }
        
        /*
         * Compute the mass fraction fields.
         */
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_primitive_variables) +
                        (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables;
                    
                    Y_x_L[si][idx] = V_x_L[si][idx]/rho_x_L[idx];
                    Y_x_R[si][idx] = V_x_R[si][idx]/rho_x_R[idx];
                }
            }
        }
        
        /*
         * Compute the volume fraction fields.
         */
        
        for (int si = 0; si < d_num_species - 1; si++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_primitive_variables) +
                        (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables;
                    
                    Z_x_L[si][idx] = V_x_L[d_num_species + 3 + si][idx];
                    Z_x_R[si][idx] = V_x_R[d_num_species + 3 + si][idx];
                }
            }
        }
        
        /*
         * Compute the sound speed field.
         */
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computeGruneisenParameter(
                gruneisen_parameter_x_L,
                density_x_L,
                primitive_variables_L[d_num_species + 2],
                mass_fractions_x_L,
                volume_fractions_x_L,
                0,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computeGruneisenParameter(
                gruneisen_parameter_x_R,
                density_x_R,
                primitive_variables_R[d_num_species + 2],
                mass_fractions_x_R,
                volume_fractions_x_R,
                0,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computePressureDerivativeWithPartialDensities(
                partial_pressure_partial_partial_densities_x_L,
                density_x_L,
                primitive_variables_L[d_num_species + 2],
                mass_fractions_x_L,
                volume_fractions_x_L,
                0,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computePressureDerivativeWithPartialDensities(
                partial_pressure_partial_partial_densities_x_R,
                density_x_R,
                primitive_variables_R[d_num_species + 2],
                mass_fractions_x_R,
                volume_fractions_x_R,
                0,
                domain);
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
            {
                // Compute the linear index.
                const int idx = (i + num_ghosts_0_primitive_variables) +
                    (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables;
                
                c_x_L[idx] = Gamma_x_L[idx]*V_x_L[d_num_species + 2][idx]/rho_x_L[idx];
                c_x_R[idx] = Gamma_x_R[idx]*V_x_R[d_num_species + 2][idx]/rho_x_R[idx];
            }
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_primitive_variables) +
                        (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables;
                    
                    c_x_L[idx] += Y_x_L[si][idx]*Psi_x_L[si][idx];
                    c_x_R[idx] += Y_x_R[si][idx]*Psi_x_R[si][idx];
                }
            }
        }
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
            {
                // Compute the linear index.
                const int idx = (i + num_ghosts_0_primitive_variables) +
                    (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables;
                
                c_x_L[idx] = sqrt(c_x_L[idx]);
                c_x_R[idx] = sqrt(c_x_R[idx]);
            }
        }
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computeInternalEnergy(
                internal_energy_x_L,
                density_x_L,
                primitive_variables_L[d_num_species + 2],
                mass_fractions_x_L,
                volume_fractions_x_L,
                0,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computeInternalEnergy(
                internal_energy_x_R,
                density_x_R,
                primitive_variables_R[d_num_species + 2],
                mass_fractions_x_R,
                volume_fractions_x_R,
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
            
            double* u = velocity->getPointer(0, 0);
            
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
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
                        rho_x_L,
                        rho_x_R,
                        c_x_L,
                        c_x_R,
                        epsilon_x_L,
                        epsilon_x_R,
                        s_x_minus,
                        s_x_plus,
                        s_x_star,
                        Chi_x_star_LR,
                        idx_flux,
                        idx,
                        d_num_species,
                        num_eqn);
                    
                    if (s_x_star > double(0))
                    {
                        u[idx_velocity] = V_x_L[d_num_species][idx] + s_x_minus*(Chi_x_star_LR - double(1));
                    }
                    else
                    {
                        u[idx_velocity] = V_x_R[d_num_species][idx] + s_x_plus*(Chi_x_star_LR - double(1));
                    }
                }
            }
        }
        else
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
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
                        rho_x_L,
                        rho_x_R,
                        c_x_L,
                        c_x_R,
                        epsilon_x_L,
                        epsilon_x_R,
                        s_x_minus,
                        s_x_plus,
                        s_x_star,
                        Chi_x_star_LR,
                        idx_flux,
                        idx,
                        d_num_species,
                        num_eqn);
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
        
        /*
         * Compute the density field.
         */
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_primitive_variables) +
                        (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables +
                        (k + num_ghosts_2_primitive_variables)*ghostcell_dim_0_primitive_variables*
                            ghostcell_dim_1_primitive_variables;
                    
                    rho_x_L[idx] = double(0);
                }
            }
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                    {
                        // Compute the linear index.
                        const int idx = (i + num_ghosts_0_primitive_variables) +
                            (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables +
                            (k + num_ghosts_2_primitive_variables)*ghostcell_dim_0_primitive_variables*
                                ghostcell_dim_1_primitive_variables;
                        
                        rho_x_L[idx] += V_x_L[si][idx];
                    }
                }
            }
        }
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_primitive_variables) +
                        (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables +
                        (k + num_ghosts_2_primitive_variables)*ghostcell_dim_0_primitive_variables*
                            ghostcell_dim_1_primitive_variables;
                    
                    rho_x_R[idx] = double(0);
                }
            }
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                    {
                        // Compute the linear index.
                        const int idx = (i + num_ghosts_0_primitive_variables) +
                            (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables +
                            (k + num_ghosts_2_primitive_variables)*ghostcell_dim_0_primitive_variables*
                                ghostcell_dim_1_primitive_variables;
                        
                        rho_x_R[idx] += V_x_R[si][idx];
                    }
                }
            }
        }
        
        /*
         * Compute the mass fraction fields.
         */
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                    {
                        // Compute the linear index.
                        const int idx = (i + num_ghosts_0_primitive_variables) +
                            (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables +
                            (k + num_ghosts_2_primitive_variables)*ghostcell_dim_0_primitive_variables*
                                ghostcell_dim_1_primitive_variables;
                        
                        Y_x_L[si][idx] = V_x_L[si][idx]/rho_x_L[idx];
                        Y_x_R[si][idx] = V_x_R[si][idx]/rho_x_R[idx];
                    }
                }
            }
        }
        
        /*
         * Compute the volume fraction fields.
         */
        
        for (int si = 0; si < d_num_species - 1; si++)
        {
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                    {
                        // Compute the linear index.
                        const int idx = (i + num_ghosts_0_primitive_variables) +
                            (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables +
                            (k + num_ghosts_2_primitive_variables)*ghostcell_dim_0_primitive_variables*
                                ghostcell_dim_1_primitive_variables;
                        
                        Z_x_L[si][idx] = V_x_L[d_num_species + 4 + si][idx];
                        Z_x_R[si][idx] = V_x_R[d_num_species + 4 + si][idx];
                    }
                }
            }
        }
        
        /*
         * Compute the sound speed field.
         */
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computeGruneisenParameter(
                gruneisen_parameter_x_L,
                density_x_L,
                primitive_variables_L[d_num_species + 3],
                mass_fractions_x_L,
                volume_fractions_x_L,
                0,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computeGruneisenParameter(
                gruneisen_parameter_x_R,
                density_x_R,
                primitive_variables_R[d_num_species + 3],
                mass_fractions_x_R,
                volume_fractions_x_R,
                0,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computePressureDerivativeWithPartialDensities(
                partial_pressure_partial_partial_densities_x_L,
                density_x_L,
                primitive_variables_L[d_num_species + 3],
                mass_fractions_x_L,
                volume_fractions_x_L,
                0,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computePressureDerivativeWithPartialDensities(
                partial_pressure_partial_partial_densities_x_R,
                density_x_R,
                primitive_variables_R[d_num_species + 3],
                mass_fractions_x_R,
                volume_fractions_x_R,
                0,
                domain);
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_primitive_variables) +
                        (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables +
                        (k + num_ghosts_2_primitive_variables)*ghostcell_dim_0_primitive_variables*
                            ghostcell_dim_1_primitive_variables;
                    
                    c_x_L[idx] = Gamma_x_L[idx]*V_x_L[d_num_species + 3][idx]/rho_x_L[idx];
                    c_x_R[idx] = Gamma_x_R[idx]*V_x_R[d_num_species + 3][idx]/rho_x_R[idx];
                }
            }
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                    {
                        // Compute the linear index.
                        const int idx = (i + num_ghosts_0_primitive_variables) +
                            (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables +
                            (k + num_ghosts_2_primitive_variables)*ghostcell_dim_0_primitive_variables*
                                ghostcell_dim_1_primitive_variables;
                        
                        c_x_L[idx] += Y_x_L[si][idx]*Psi_x_L[si][idx];
                        c_x_R[idx] += Y_x_R[si][idx]*Psi_x_R[si][idx];
                    }
                }
            }
        }
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_primitive_variables) +
                        (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables +
                        (k + num_ghosts_2_primitive_variables)*ghostcell_dim_0_primitive_variables*
                            ghostcell_dim_1_primitive_variables;
                    
                    c_x_L[idx] = sqrt(c_x_L[idx]);
                    c_x_R[idx] = sqrt(c_x_R[idx]);
                }
            }
        }
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computeInternalEnergy(
                internal_energy_x_L,
                density_x_L,
                primitive_variables_L[d_num_species + 3],
                mass_fractions_x_L,
                volume_fractions_x_L,
                0,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computeInternalEnergy(
                internal_energy_x_R,
                density_x_R,
                primitive_variables_R[d_num_species + 3],
                mass_fractions_x_R,
                volume_fractions_x_R,
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
            
            double* u = velocity->getPointer(0, 0);
            
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
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
                            rho_x_L,
                            rho_x_R,
                            c_x_L,
                            c_x_R,
                            epsilon_x_L,
                            epsilon_x_R,
                            s_x_minus,
                            s_x_plus,
                            s_x_star,
                            Chi_x_star_LR,
                            idx_flux,
                            idx,
                            d_num_species,
                            num_eqn);
                        
                        if (s_x_star > double(0))
                        {
                            u[idx_velocity] = V_x_L[d_num_species][idx] + s_x_minus*(Chi_x_star_LR - double(1));
                        }
                        else
                        {
                            u[idx_velocity] = V_x_R[d_num_species][idx] + s_x_plus*(Chi_x_star_LR - double(1));
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
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
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
                            rho_x_L,
                            rho_x_R,
                            c_x_L,
                            c_x_R,
                            epsilon_x_L,
                            epsilon_x_R,
                            s_x_minus,
                            s_x_plus,
                            s_x_star,
                            Chi_x_star_LR,
                            idx_flux,
                            idx,
                            d_num_species,
                            num_eqn);
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
FlowModelRiemannSolverFiveEqnAllaire::computeConvectiveFluxAndVelocityInYDirectionFromPrimitiveVariablesHLLC_HLL(
    boost::shared_ptr<pdat::SideData<double> > convective_flux,
    boost::shared_ptr<pdat::SideData<double> > velocity,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables_B,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables_T,
    const hier::Box& domain,
    bool compute_velocity) const
{
    boost::shared_ptr<FlowModel> flow_model_tmp = d_flow_model.lock();
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
    
    const boost::shared_ptr<EquationOfStateMixingRules> equation_of_state_mixing_rules =
        flow_model_tmp->getEquationOfStateMixingRules();
    
    /*
     * Get the pointers to the side data of convective flux and primitive variables.
     */
    
    std::vector<double*> F_y;
    F_y.reserve(num_eqn);
    for (int ei = 0; ei < num_eqn; ei++)
    {
        F_y.push_back(convective_flux->getPointer(1, ei));
    }
    
    std::vector<double*> V_y_B;
    std::vector<double*> V_y_T;
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
    
    boost::shared_ptr<pdat::SideData<double> > density_y_B(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_primitive_variables,
            direction_y));
    
    boost::shared_ptr<pdat::SideData<double> > density_y_T(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_primitive_variables,
            direction_y));
    
    boost::shared_ptr<pdat::SideData<double> > mass_fractions_y_B(
        new pdat::SideData<double>(interior_box, d_num_species, num_ghosts_primitive_variables,
            direction_y));
    
    boost::shared_ptr<pdat::SideData<double> > mass_fractions_y_T(
        new pdat::SideData<double>(interior_box, d_num_species, num_ghosts_primitive_variables,
            direction_y));
    
    boost::shared_ptr<pdat::SideData<double> > volume_fractions_y_B(
        new pdat::SideData<double>(interior_box, d_num_species - 1, num_ghosts_primitive_variables,
            direction_y));
    
    boost::shared_ptr<pdat::SideData<double> > volume_fractions_y_T(
        new pdat::SideData<double>(interior_box, d_num_species - 1, num_ghosts_primitive_variables,
            direction_y));
    
    boost::shared_ptr<pdat::SideData<double> > gruneisen_parameter_y_B(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_primitive_variables,
            direction_y));
    
    boost::shared_ptr<pdat::SideData<double> > gruneisen_parameter_y_T(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_primitive_variables,
            direction_y));
    
    boost::shared_ptr<pdat::SideData<double> > partial_pressure_partial_partial_densities_y_B(
            new pdat::SideData<double>(interior_box, d_num_species, num_ghosts_primitive_variables,
            direction_y));
    
    boost::shared_ptr<pdat::SideData<double> > partial_pressure_partial_partial_densities_y_T(
            new pdat::SideData<double>(interior_box, d_num_species, num_ghosts_primitive_variables,
            direction_y));
    
    boost::shared_ptr<pdat::SideData<double> > sound_speed_y_B(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_primitive_variables,
            direction_y));
    
    boost::shared_ptr<pdat::SideData<double> > sound_speed_y_T(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_primitive_variables,
            direction_y));
    
    boost::shared_ptr<pdat::SideData<double> > internal_energy_y_B(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_primitive_variables,
            direction_y));
    
    boost::shared_ptr<pdat::SideData<double> > internal_energy_y_T(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_primitive_variables,
            direction_y));
    
    /*
     * Get the pointers to the temporary data.
     */
    
    double* rho_y_B = density_y_B->getPointer(1, 0);
    double* rho_y_T = density_y_T->getPointer(1, 0);
    
    std::vector<double*> Y_y_B;
    std::vector<double*> Y_y_T;
    Y_y_B.reserve(d_num_species);
    Y_y_T.reserve(d_num_species);
    for (int si = 0; si < d_num_species; si++)
    {
        Y_y_B.push_back(mass_fractions_y_B->getPointer(1, si));
        Y_y_T.push_back(mass_fractions_y_T->getPointer(1, si));
    }
    
    std::vector<double*> Z_y_B;
    std::vector<double*> Z_y_T;
    Z_y_B.reserve(d_num_species);
    Z_y_T.reserve(d_num_species);
    for (int si = 0; si < d_num_species - 1; si++)
    {
        Z_y_B.push_back(volume_fractions_y_B->getPointer(1, si));
        Z_y_T.push_back(volume_fractions_y_T->getPointer(1, si));
    }
    
    double* Gamma_y_B = gruneisen_parameter_y_B->getPointer(1, 0);
    double* Gamma_y_T = gruneisen_parameter_y_T->getPointer(1, 0);
    
    std::vector<double*> Psi_y_B;
    std::vector<double*> Psi_y_T;
    Psi_y_B.reserve(d_num_species);
    Psi_y_T.reserve(d_num_species);
    for (int si = 0; si < d_num_species; si++)
    {
        Psi_y_B.push_back(partial_pressure_partial_partial_densities_y_B->getPointer(1, si));
        Psi_y_T.push_back(partial_pressure_partial_partial_densities_y_T->getPointer(1, si));
    }
    
    double* c_y_B = sound_speed_y_B->getPointer(1, 0);
    double* c_y_T = sound_speed_y_T->getPointer(1, 0);
    
    double* epsilon_y_B = internal_energy_y_B->getPointer(1, 0);
    double* epsilon_y_T = internal_energy_y_T->getPointer(1, 0);
    
    double s_y_minus = double(0);
    double s_y_plus  = double(0);
    double s_y_star  = double(0);
    
    double Chi_y_star_BT = double(0);
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelRiemannSolverFiveEqnAllaire::"
            << "computeConvectiveFluxAndVelocityInYDirectionFromPrimitiveVariablesHLLC_HLL()\n"
            << "There is no y direction for 1D problem."
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
        
        /*
         * Compute the density field.
         */
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear index.
                const int idx = (i + num_ghosts_0_primitive_variables) +
                    (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables;
                
                rho_y_B[idx] = double(0);
            }
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_primitive_variables) +
                        (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables;
                    
                    rho_y_B[idx] += V_y_B[si][idx];
                }
            }
        }
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear index.
                const int idx = (i + num_ghosts_0_primitive_variables) +
                    (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables;
                
                rho_y_T[idx] = double(0);
            }
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_primitive_variables) +
                        (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables;
                    
                    rho_y_T[idx] += V_y_T[si][idx];
                }
            }
        }
        
        /*
         * Compute the mass fraction fields.
         */
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_primitive_variables) +
                        (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables;
                    
                    Y_y_B[si][idx] = V_y_B[si][idx]/rho_y_B[idx];
                    Y_y_T[si][idx] = V_y_T[si][idx]/rho_y_T[idx];
                }
            }
        }
        
        /*
         * Compute the volume fraction fields.
         */
        
        for (int si = 0; si < d_num_species - 1; si++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_primitive_variables) +
                        (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables;
                    
                    Z_y_B[si][idx] = V_y_B[d_num_species + 3 + si][idx];
                    Z_y_T[si][idx] = V_y_T[d_num_species + 3 + si][idx];
                }
            }
        }
        
        /*
         * Compute the sound speed field.
         */
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computeGruneisenParameter(
                gruneisen_parameter_y_B,
                density_y_B,
                primitive_variables_B[d_num_species + 2],
                mass_fractions_y_B,
                volume_fractions_y_B,
                1,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computeGruneisenParameter(
                gruneisen_parameter_y_T,
                density_y_T,
                primitive_variables_T[d_num_species + 2],
                mass_fractions_y_T,
                volume_fractions_y_T,
                1,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computePressureDerivativeWithPartialDensities(
                partial_pressure_partial_partial_densities_y_B,
                density_y_B,
                primitive_variables_B[d_num_species + 2],
                mass_fractions_y_B,
                volume_fractions_y_B,
                1,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computePressureDerivativeWithPartialDensities(
                partial_pressure_partial_partial_densities_y_T,
                density_y_T,
                primitive_variables_T[d_num_species + 2],
                mass_fractions_y_T,
                volume_fractions_y_T,
                1,
                domain);
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear index.
                const int idx = (i + num_ghosts_0_primitive_variables) +
                    (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables;
                
                c_y_B[idx] = Gamma_y_B[idx]*V_y_B[d_num_species + 2][idx]/rho_y_B[idx];
                c_y_T[idx] = Gamma_y_T[idx]*V_y_T[d_num_species + 2][idx]/rho_y_T[idx];
            }
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_primitive_variables) +
                        (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables;
                    
                    c_y_B[idx] += Y_y_B[si][idx]*Psi_y_B[si][idx];
                    c_y_T[idx] += Y_y_T[si][idx]*Psi_y_T[si][idx];
                }
            }
        }
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear index.
                const int idx = (i + num_ghosts_0_primitive_variables) +
                    (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables;
                
                c_y_B[idx] = sqrt(c_y_B[idx]);
                c_y_T[idx] = sqrt(c_y_T[idx]);
            }
        }
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computeInternalEnergy(
                internal_energy_y_B,
                density_y_B,
                primitive_variables_B[d_num_species + 2],
                mass_fractions_y_B,
                volume_fractions_y_B,
                1,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computeInternalEnergy(
                internal_energy_y_T,
                density_y_T,
                primitive_variables_T[d_num_species + 2],
                mass_fractions_y_T,
                volume_fractions_y_T,
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
            
            double* v = velocity->getPointer(1, 1);
            
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
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
                        rho_y_B,
                        rho_y_T,
                        c_y_B,
                        c_y_T,
                        epsilon_y_B,
                        epsilon_y_T,
                        s_y_minus,
                        s_y_plus,
                        s_y_star,
                        Chi_y_star_BT,
                        idx_flux,
                        idx,
                        d_num_species,
                        num_eqn);
                    
                    if (s_y_star > double(0))
                    {
                        v[idx_velocity] = V_y_B[d_num_species + 1][idx] + s_y_minus*(Chi_y_star_BT - double(1));
                    }
                    else
                    {
                        v[idx_velocity] = V_y_T[d_num_species + 1][idx] + s_y_plus*(Chi_y_star_BT - double(1));
                    }
                }
            }
        }
        else
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
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
                        rho_y_B,
                        rho_y_T,
                        c_y_B,
                        c_y_T,
                        epsilon_y_B,
                        epsilon_y_T,
                        s_y_minus,
                        s_y_plus,
                        s_y_star,
                        Chi_y_star_BT,
                        idx_flux,
                        idx,
                        d_num_species,
                        num_eqn);
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
        
        /*
         * Compute the density field.
         */
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_primitive_variables) +
                        (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables +
                        (k + num_ghosts_2_primitive_variables)*ghostcell_dim_0_primitive_variables*
                            ghostcell_dim_1_primitive_variables;
                    
                    rho_y_B[idx] = double(0);
                }
            }
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx = (i + num_ghosts_0_primitive_variables) +
                            (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables +
                            (k + num_ghosts_2_primitive_variables)*ghostcell_dim_0_primitive_variables*
                                ghostcell_dim_1_primitive_variables;
                        
                        rho_y_B[idx] += V_y_B[si][idx];
                    }
                }
            }
        }
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_primitive_variables) +
                        (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables +
                        (k + num_ghosts_2_primitive_variables)*ghostcell_dim_0_primitive_variables*
                            ghostcell_dim_1_primitive_variables;
                    
                    rho_y_T[idx] = double(0);
                }
            }
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx = (i + num_ghosts_0_primitive_variables) +
                            (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables +
                            (k + num_ghosts_2_primitive_variables)*ghostcell_dim_0_primitive_variables*
                                ghostcell_dim_1_primitive_variables;
                        
                        rho_y_T[idx] += V_y_T[si][idx];
                    }
                }
            }
        }
        
        /*
         * Compute the mass fraction fields.
         */
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx = (i + num_ghosts_0_primitive_variables) +
                            (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables +
                            (k + num_ghosts_2_primitive_variables)*ghostcell_dim_0_primitive_variables*
                                ghostcell_dim_1_primitive_variables;
                        
                        Y_y_B[si][idx] = V_y_B[si][idx]/rho_y_B[idx];
                        Y_y_T[si][idx] = V_y_T[si][idx]/rho_y_T[idx];
                    }
                }
            }
        }
        
        /*
         * Compute the volume fraction fields.
         */
        
        for (int si = 0; si < d_num_species - 1; si++)
        {
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx = (i + num_ghosts_0_primitive_variables) +
                            (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables +
                            (k + num_ghosts_2_primitive_variables)*ghostcell_dim_0_primitive_variables*
                                ghostcell_dim_1_primitive_variables;
                        
                        Z_y_B[si][idx] = V_y_B[d_num_species + 4 + si][idx];
                        Z_y_T[si][idx] = V_y_T[d_num_species + 4 + si][idx];
                    }
                }
            }
        }
        
        /*
         * Compute the sound speed field.
         */
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computeGruneisenParameter(
                gruneisen_parameter_y_B,
                density_y_B,
                primitive_variables_B[d_num_species + 3],
                mass_fractions_y_B,
                volume_fractions_y_B,
                1,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computeGruneisenParameter(
                gruneisen_parameter_y_T,
                density_y_T,
                primitive_variables_T[d_num_species + 3],
                mass_fractions_y_T,
                volume_fractions_y_T,
                1,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computePressureDerivativeWithPartialDensities(
                partial_pressure_partial_partial_densities_y_B,
                density_y_B,
                primitive_variables_B[d_num_species + 3],
                mass_fractions_y_B,
                volume_fractions_y_B,
                1,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computePressureDerivativeWithPartialDensities(
                partial_pressure_partial_partial_densities_y_T,
                density_y_T,
                primitive_variables_T[d_num_species + 3],
                mass_fractions_y_T,
                volume_fractions_y_T,
                1,
                domain);
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_primitive_variables) +
                        (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables +
                        (k + num_ghosts_2_primitive_variables)*ghostcell_dim_0_primitive_variables*
                            ghostcell_dim_1_primitive_variables;
                    
                    c_y_B[idx] = Gamma_y_B[idx]*V_y_B[d_num_species + 3][idx]/rho_y_B[idx];
                    c_y_T[idx] = Gamma_y_T[idx]*V_y_T[d_num_species + 3][idx]/rho_y_T[idx];
                }
            }
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx = (i + num_ghosts_0_primitive_variables) +
                            (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables +
                            (k + num_ghosts_2_primitive_variables)*ghostcell_dim_0_primitive_variables*
                                ghostcell_dim_1_primitive_variables;
                        
                        c_y_B[idx] += Y_y_B[si][idx]*Psi_y_B[si][idx];
                        c_y_T[idx] += Y_y_T[si][idx]*Psi_y_T[si][idx];
                    }
                }
            }
        }
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_primitive_variables) +
                        (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables +
                        (k + num_ghosts_2_primitive_variables)*ghostcell_dim_0_primitive_variables*
                            ghostcell_dim_1_primitive_variables;
                    
                    c_y_B[idx] = sqrt(c_y_B[idx]);
                    c_y_T[idx] = sqrt(c_y_T[idx]);
                }
            }
        }
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computeInternalEnergy(
                internal_energy_y_B,
                density_y_B,
                primitive_variables_B[d_num_species + 3],
                mass_fractions_y_B,
                volume_fractions_y_B,
                1,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computeInternalEnergy(
                internal_energy_y_T,
                density_y_T,
                primitive_variables_T[d_num_species + 3],
                mass_fractions_y_T,
                volume_fractions_y_T,
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
            
            double* v = velocity->getPointer(1, 1);
            
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
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
                            rho_y_B,
                            rho_y_T,
                            c_y_B,
                            c_y_T,
                            epsilon_y_B,
                            epsilon_y_T,
                            s_y_minus,
                            s_y_plus,
                            s_y_star,
                            Chi_y_star_BT,
                            idx_flux,
                            idx,
                            d_num_species,
                            num_eqn);
                        
                        if (s_y_star > double(0))
                        {
                            v[idx_velocity] = V_y_B[d_num_species + 1][idx] + s_y_minus*(Chi_y_star_BT - double(1));
                        }
                        else
                        {
                            v[idx_velocity] = V_y_T[d_num_species + 1][idx] + s_y_plus*(Chi_y_star_BT - double(1));
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
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
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
                            rho_y_B,
                            rho_y_T,
                            c_y_B,
                            c_y_T,
                            epsilon_y_B,
                            epsilon_y_T,
                            s_y_minus,
                            s_y_plus,
                            s_y_star,
                            Chi_y_star_BT,
                            idx_flux,
                            idx,
                            d_num_species,
                            num_eqn);
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
FlowModelRiemannSolverFiveEqnAllaire::computeConvectiveFluxAndVelocityInZDirectionFromPrimitiveVariablesHLLC_HLL(
    boost::shared_ptr<pdat::SideData<double> > convective_flux,
    boost::shared_ptr<pdat::SideData<double> > velocity,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables_B,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables_F,
    const hier::Box& domain,
    bool compute_velocity) const
{
    boost::shared_ptr<FlowModel> flow_model_tmp = d_flow_model.lock();
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
    
    const boost::shared_ptr<EquationOfStateMixingRules> equation_of_state_mixing_rules =
        flow_model_tmp->getEquationOfStateMixingRules();
    
    /*
     * Get the pointers to the side data of convective flux and primitive variables.
     */
    
    std::vector<double*> F_z;
    F_z.reserve(num_eqn);
    for (int ei = 0; ei < num_eqn; ei++)
    {
        F_z.push_back(convective_flux->getPointer(2, ei));
    }
    
    std::vector<double*> V_z_B;
    std::vector<double*> V_z_F;
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
    
    boost::shared_ptr<pdat::SideData<double> > density_z_B(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_primitive_variables,
            direction_z));
    
    boost::shared_ptr<pdat::SideData<double> > density_z_F(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_primitive_variables,
            direction_z));
    
    boost::shared_ptr<pdat::SideData<double> > mass_fractions_z_B(
        new pdat::SideData<double>(interior_box, d_num_species, num_ghosts_primitive_variables,
            direction_z));
    
    boost::shared_ptr<pdat::SideData<double> > mass_fractions_z_F(
        new pdat::SideData<double>(interior_box, d_num_species, num_ghosts_primitive_variables,
            direction_z));
    
    boost::shared_ptr<pdat::SideData<double> > volume_fractions_z_B(
        new pdat::SideData<double>(interior_box, d_num_species - 1, num_ghosts_primitive_variables,
            direction_z));
    
    boost::shared_ptr<pdat::SideData<double> > volume_fractions_z_F(
        new pdat::SideData<double>(interior_box, d_num_species - 1, num_ghosts_primitive_variables,
            direction_z));
    
    boost::shared_ptr<pdat::SideData<double> > gruneisen_parameter_z_B(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_primitive_variables,
            direction_z));
    
    boost::shared_ptr<pdat::SideData<double> > gruneisen_parameter_z_F(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_primitive_variables,
            direction_z));
    
    boost::shared_ptr<pdat::SideData<double> > partial_pressure_partial_partial_densities_z_B(
            new pdat::SideData<double>(interior_box, d_num_species, num_ghosts_primitive_variables,
            direction_z));
    
    boost::shared_ptr<pdat::SideData<double> > partial_pressure_partial_partial_densities_z_F(
            new pdat::SideData<double>(interior_box, d_num_species, num_ghosts_primitive_variables,
            direction_z));
    
    boost::shared_ptr<pdat::SideData<double> > sound_speed_z_B(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_primitive_variables,
            direction_z));
    
    boost::shared_ptr<pdat::SideData<double> > sound_speed_z_F(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_primitive_variables,
            direction_z));
    
    boost::shared_ptr<pdat::SideData<double> > internal_energy_z_B(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_primitive_variables,
            direction_z));
    
    boost::shared_ptr<pdat::SideData<double> > internal_energy_z_F(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_primitive_variables,
            direction_z));
    
    /*
     * Get the pointers to the temporary data.
     */
    
    double* rho_z_B = density_z_B->getPointer(2, 0);
    double* rho_z_F = density_z_F->getPointer(2, 0);
    
    std::vector<double*> Y_z_B;
    std::vector<double*> Y_z_F;
    Y_z_B.reserve(d_num_species);
    Y_z_F.reserve(d_num_species);
    for (int si = 0; si < d_num_species; si++)
    {
        Y_z_B.push_back(mass_fractions_z_B->getPointer(2, si));
        Y_z_F.push_back(mass_fractions_z_F->getPointer(2, si));
    }
    
    std::vector<double*> Z_z_B;
    std::vector<double*> Z_z_F;
    Z_z_B.reserve(d_num_species);
    Z_z_F.reserve(d_num_species);
    for (int si = 0; si < d_num_species - 1; si++)
    {
        Z_z_B.push_back(volume_fractions_z_B->getPointer(2, si));
        Z_z_F.push_back(volume_fractions_z_F->getPointer(2, si));
    }
    
    double* Gamma_z_B = gruneisen_parameter_z_B->getPointer(2, 0);
    double* Gamma_z_F = gruneisen_parameter_z_F->getPointer(2, 0);
    
    std::vector<double*> Psi_z_B;
    std::vector<double*> Psi_z_F;
    Psi_z_B.reserve(d_num_species);
    Psi_z_F.reserve(d_num_species);
    for (int si = 0; si < d_num_species; si++)
    {
        Psi_z_B.push_back(partial_pressure_partial_partial_densities_z_B->getPointer(2, si));
        Psi_z_F.push_back(partial_pressure_partial_partial_densities_z_F->getPointer(2, si));
    }
    
    double* c_z_B = sound_speed_z_B->getPointer(2, 0);
    double* c_z_F = sound_speed_z_F->getPointer(2, 0);
    
    double* epsilon_z_B = internal_energy_z_B->getPointer(2, 0);
    double* epsilon_z_F = internal_energy_z_F->getPointer(2, 0);
    
    double s_z_minus = double(0);
    double s_z_plus  = double(0);
    double s_z_star  = double(0);
    
    double Chi_z_star_BF = double(0);
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelRiemannSolverFiveEqnAllaire::"
            << "computeConvectiveFluxAndVelocityInZDirectionFromPrimitiveVariablesHLLC_HLL()\n"
            << "There is no z direction for 1D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelRiemannSolverFiveEqnAllaire::"
            << "computeConvectiveFluxAndVelocityInZDirectionFromPrimitiveVariablesHLLC_HLL()\n"
            << "There is no z direction for 2D problem."
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
        
        /*
         * Compute the density field.
         */
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2 + 1; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_primitive_variables) +
                        (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables +
                        (k + num_ghosts_2_primitive_variables)*ghostcell_dim_0_primitive_variables*
                            ghostcell_dim_1_primitive_variables;
                    
                    rho_z_B[idx] = double(0);
                }
            }
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2 + 1; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx = (i + num_ghosts_0_primitive_variables) +
                            (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables +
                            (k + num_ghosts_2_primitive_variables)*ghostcell_dim_0_primitive_variables*
                                ghostcell_dim_1_primitive_variables;
                        
                        rho_z_B[idx] += V_z_B[si][idx];
                    }
                }
            }
        }
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2 + 1; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_primitive_variables) +
                        (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables +
                        (k + num_ghosts_2_primitive_variables)*ghostcell_dim_0_primitive_variables*
                            ghostcell_dim_1_primitive_variables;
                    
                    rho_z_F[idx] = double(0);
                }
            }
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2 + 1; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx = (i + num_ghosts_0_primitive_variables) +
                            (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables +
                            (k + num_ghosts_2_primitive_variables)*ghostcell_dim_0_primitive_variables*
                                ghostcell_dim_1_primitive_variables;
                        
                        rho_z_F[idx] += V_z_F[si][idx];
                    }
                }
            }
        }
        
        /*
         * Compute the mass fraction fields.
         */
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2 + 1; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx = (i + num_ghosts_0_primitive_variables) +
                            (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables +
                            (k + num_ghosts_2_primitive_variables)*ghostcell_dim_0_primitive_variables*
                                ghostcell_dim_1_primitive_variables;
                        
                        Y_z_B[si][idx] = V_z_B[si][idx]/rho_z_B[idx];
                        Y_z_F[si][idx] = V_z_F[si][idx]/rho_z_F[idx];
                    }
                }
            }
        }
        
        /*
         * Compute the volume fraction fields.
         */
        
        for (int si = 0; si < d_num_species - 1; si++)
        {
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2 + 1; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx = (i + num_ghosts_0_primitive_variables) +
                            (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables +
                            (k + num_ghosts_2_primitive_variables)*ghostcell_dim_0_primitive_variables*
                                ghostcell_dim_1_primitive_variables;
                        
                        Z_z_B[si][idx] = V_z_B[d_num_species + 4 + si][idx];
                        Z_z_F[si][idx] = V_z_F[d_num_species + 4 + si][idx];
                    }
                }
            }
        }
        
        /*
         * Compute the sound speed field.
         */
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computeGruneisenParameter(
                gruneisen_parameter_z_B,
                density_z_B,
                primitive_variables_B[d_num_species + 3],
                mass_fractions_z_B,
                volume_fractions_z_B,
                2,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computeGruneisenParameter(
                gruneisen_parameter_z_F,
                density_z_F,
                primitive_variables_F[d_num_species + 3],
                mass_fractions_z_F,
                volume_fractions_z_F,
                2,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computePressureDerivativeWithPartialDensities(
                partial_pressure_partial_partial_densities_z_B,
                density_z_B,
                primitive_variables_B[d_num_species + 3],
                mass_fractions_z_B,
                volume_fractions_z_B,
                2,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computePressureDerivativeWithPartialDensities(
                partial_pressure_partial_partial_densities_z_F,
                density_z_F,
                primitive_variables_F[d_num_species + 3],
                mass_fractions_z_F,
                volume_fractions_z_F,
                2,
                domain);
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2 + 1; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_primitive_variables) +
                        (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables +
                        (k + num_ghosts_2_primitive_variables)*ghostcell_dim_0_primitive_variables*
                            ghostcell_dim_1_primitive_variables;
                    
                    c_z_B[idx] = Gamma_z_B[idx]*V_z_B[d_num_species + 3][idx]/rho_z_B[idx];
                    c_z_F[idx] = Gamma_z_F[idx]*V_z_F[d_num_species + 3][idx]/rho_z_F[idx];
                }
            }
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2 + 1; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx = (i + num_ghosts_0_primitive_variables) +
                            (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables +
                            (k + num_ghosts_2_primitive_variables)*ghostcell_dim_0_primitive_variables*
                                ghostcell_dim_1_primitive_variables;
                        
                        c_z_B[idx] += Y_z_B[si][idx]*Psi_z_B[si][idx];
                        c_z_F[idx] += Y_z_F[si][idx]*Psi_z_F[si][idx];
                    }
                }
            }
        }
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2 + 1; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + num_ghosts_0_primitive_variables) +
                        (j + num_ghosts_1_primitive_variables)*ghostcell_dim_0_primitive_variables +
                        (k + num_ghosts_2_primitive_variables)*ghostcell_dim_0_primitive_variables*
                            ghostcell_dim_1_primitive_variables;
                    
                    c_z_B[idx] = sqrt(c_z_B[idx]);
                    c_z_F[idx] = sqrt(c_z_F[idx]);
                }
            }
        }
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computeInternalEnergy(
                internal_energy_z_B,
                density_z_B,
                primitive_variables_B[d_num_species + 3],
                mass_fractions_z_B,
                volume_fractions_z_B,
                2,
                domain);
        
        flow_model_tmp->getEquationOfStateMixingRules()->
            computeInternalEnergy(
                internal_energy_z_F,
                density_z_F,
                primitive_variables_F[d_num_species + 3],
                mass_fractions_z_F,
                volume_fractions_z_F,
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
            
            double* w = velocity->getPointer(2, 2);
            
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2 + 1; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
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
                            rho_z_B,
                            rho_z_F,
                            c_z_B,
                            c_z_F,
                            epsilon_z_B,
                            epsilon_z_F,
                            s_z_minus,
                            s_z_plus,
                            s_z_star,
                            Chi_z_star_BF,
                            idx_flux,
                            idx,
                            d_num_species,
                            num_eqn);
                        
                        if (s_z_star > double(0))
                        {
                            w[idx_velocity] = V_z_B[d_num_species + 2][idx] + s_z_minus*(Chi_z_star_BF - double(1));
                        }
                        else
                        {
                            w[idx_velocity] = V_z_F[d_num_species + 2][idx] + s_z_plus*(Chi_z_star_BF - double(1));
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
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
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
                            rho_z_B,
                            rho_z_F,
                            c_z_B,
                            c_z_F,
                            epsilon_z_B,
                            epsilon_z_F,
                            s_z_minus,
                            s_z_plus,
                            s_z_star,
                            Chi_z_star_BF,
                            idx_flux,
                            idx,
                            d_num_species,
                            num_eqn);
                    }
                }
            }
        }
    }
}
