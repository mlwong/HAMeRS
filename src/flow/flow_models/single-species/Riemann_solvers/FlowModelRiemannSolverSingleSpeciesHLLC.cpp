#include "flow/flow_models/single-species/FlowModelRiemannSolverSingleSpecies.hpp"

/*
 * Compute the local convective flux in the x-direction from conservative variables with
 * 1D HLLC Riemann solver.
 */
static inline __attribute__((always_inline)) void
computeLocalConvectiveFluxInXDirectionFromConservativeVariablesHLLC1D(
    double** F_x,
    double** Q_x_L,
    double** Q_x_R,
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
    const int& idx)
{
    u_x_L = Q_x_L[1][idx]/Q_x_L[0][idx];
    u_x_R = Q_x_R[1][idx]/Q_x_R[0][idx];
    
    const double u_x_average = double(1)/double(2)*(u_x_L + u_x_R);
    const double c_x_average = double(1)/double(2)*(c_x_L[idx] + c_x_R[idx]);
    
    const double s_x_L = fmin(u_x_average - c_x_average, u_x_L - c_x_L[idx]);
    const double s_x_R = fmax(u_x_average + c_x_average, u_x_R + c_x_R[idx]);
    
    s_x_minus = fmin(double(0), s_x_L);
    s_x_plus  = fmax(double(0), s_x_R);
    
    s_x_star = (p_x_R[idx] - p_x_L[idx] +
        Q_x_L[1][idx]*(s_x_L - u_x_L) - Q_x_R[1][idx]*(s_x_R - u_x_R))/
        (Q_x_L[0][idx]*(s_x_L - u_x_L) - Q_x_R[0][idx]*(s_x_R - u_x_R));
    
    double Q_x_star_LR[3];
    double F_x_LR[3];
    
    if (s_x_star > double(0))
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
 * 2D HLLC Riemann solver.
 */
static inline __attribute__((always_inline)) void
computeLocalConvectiveFluxInXDirectionFromConservativeVariablesHLLC2D(
    double** F_x,
    double** Q_x_L,
    double** Q_x_R,
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
    const int& idx)
{
    u_x_L = Q_x_L[1][idx]/Q_x_L[0][idx];
    u_x_R = Q_x_R[1][idx]/Q_x_R[0][idx];
    
    const double u_x_average = double(1)/double(2)*(u_x_L + u_x_R);
    const double c_x_average = double(1)/double(2)*(c_x_L[idx] + c_x_R[idx]);
    
    const double s_x_L = fmin(u_x_average - c_x_average, u_x_L - c_x_L[idx]);
    const double s_x_R = fmax(u_x_average + c_x_average, u_x_R + c_x_R[idx]);
    
    s_x_minus = fmin(double(0), s_x_L);
    s_x_plus  = fmax(double(0), s_x_R);
    
    s_x_star = (p_x_R[idx] - p_x_L[idx] +
        Q_x_L[1][idx]*(s_x_L - u_x_L) - Q_x_R[1][idx]*(s_x_R - u_x_R))/
        (Q_x_L[0][idx]*(s_x_L - u_x_L) - Q_x_R[0][idx]*(s_x_R - u_x_R));
    
    double Q_x_star_LR[4];
    double F_x_LR[4];
    
    if (s_x_star > double(0))
    {
        Chi_x_star_LR = (s_x_L - u_x_L)/(s_x_L - s_x_star);
        
        Q_x_star_LR[0] = Chi_x_star_LR*Q_x_L[0][idx];
        Q_x_star_LR[1] = Chi_x_star_LR*Q_x_L[0][idx]*s_x_star;
        Q_x_star_LR[2] = Chi_x_star_LR*Q_x_L[2][idx];
        Q_x_star_LR[3] = Chi_x_star_LR*(Q_x_L[3][idx] + (s_x_star - u_x_L)*(Q_x_L[0][idx]*s_x_star +
            p_x_L[idx]/(s_x_L - u_x_L)));
        
        F_x_LR[0] = Q_x_L[1][idx];
        F_x_LR[1] = u_x_L*Q_x_L[1][idx] + p_x_L[idx];
        F_x_LR[2] = u_x_L*Q_x_L[2][idx];
        F_x_LR[3] = u_x_L*(Q_x_L[3][idx] + p_x_L[idx]);
        
        for (int ei = 0; ei < 4; ei++)
        {
            F_x[ei][idx_flux] = F_x_LR[ei] + s_x_minus*(Q_x_star_LR[ei] - Q_x_L[ei][idx]);
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
        
        F_x_LR[0] = Q_x_R[1][idx];
        F_x_LR[1] = u_x_R*Q_x_R[1][idx] + p_x_R[idx];
        F_x_LR[2] = u_x_R*Q_x_R[2][idx];
        F_x_LR[3] = u_x_R*(Q_x_R[3][idx] + p_x_R[idx]);
        
        for (int ei = 0; ei < 4; ei++)
        {
            F_x[ei][idx_flux] = F_x_LR[ei] + s_x_plus*(Q_x_star_LR[ei] - Q_x_R[ei][idx]);
        }
    }
}


/*
 * Compute the local convective flux in the x-direction from conservative variables with
 * 3D HLLC Riemann solver.
 */
static inline __attribute__((always_inline)) void
computeLocalConvectiveFluxInXDirectionFromConservativeVariablesHLLC3D(
    double** F_x,
    double** Q_x_L,
    double** Q_x_R,
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
    const int& idx)
{
    u_x_L = Q_x_L[1][idx]/Q_x_L[0][idx];
    u_x_R = Q_x_R[1][idx]/Q_x_R[0][idx];
    
    const double u_x_average = double(1)/double(2)*(u_x_L + u_x_R);
    const double c_x_average = double(1)/double(2)*(c_x_L[idx] + c_x_R[idx]);
    
    const double s_x_L = fmin(u_x_average - c_x_average, u_x_L - c_x_L[idx]);
    const double s_x_R = fmax(u_x_average + c_x_average, u_x_R + c_x_R[idx]);
    
    s_x_minus = fmin(double(0), s_x_L);
    s_x_plus  = fmax(double(0), s_x_R);
    
    s_x_star = (p_x_R[idx] - p_x_L[idx] +
        Q_x_L[1][idx]*(s_x_L - u_x_L) - Q_x_R[1][idx]*(s_x_R - u_x_R))/
        (Q_x_L[0][idx]*(s_x_L - u_x_L) - Q_x_R[0][idx]*(s_x_R - u_x_R));
    
    double Q_x_star_LR[5];
    double F_x_LR[5];
    
    if (s_x_star > double(0))
    {
        Chi_x_star_LR = (s_x_L - u_x_L)/(s_x_L - s_x_star);
        
        Q_x_star_LR[0] = Chi_x_star_LR*Q_x_L[0][idx];
        Q_x_star_LR[1] = Chi_x_star_LR*Q_x_L[0][idx]*s_x_star;
        Q_x_star_LR[2] = Chi_x_star_LR*Q_x_L[2][idx];
        Q_x_star_LR[3] = Chi_x_star_LR*Q_x_L[3][idx];
        Q_x_star_LR[4] = Chi_x_star_LR*(Q_x_L[4][idx] + (s_x_star - u_x_L)*(Q_x_L[0][idx]*s_x_star +
            p_x_L[idx]/(s_x_L - u_x_L)));
        
        F_x_LR[0] = Q_x_L[1][idx];
        F_x_LR[1] = u_x_L*Q_x_L[1][idx] + p_x_L[idx];
        F_x_LR[2] = u_x_L*Q_x_L[2][idx];
        F_x_LR[3] = u_x_L*Q_x_L[3][idx];
        F_x_LR[4] = u_x_L*(Q_x_L[4][idx] + p_x_L[idx]);
        
        for (int ei = 0; ei < 5; ei++)
        {
            F_x[ei][idx_flux] = F_x_LR[ei] + s_x_minus*(Q_x_star_LR[ei] - Q_x_L[ei][idx]);
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
        
        F_x_LR[0] = Q_x_R[1][idx];
        F_x_LR[1] = u_x_R*Q_x_R[1][idx] + p_x_R[idx];
        F_x_LR[2] = u_x_R*Q_x_R[2][idx];
        F_x_LR[3] = u_x_R*Q_x_R[3][idx];
        F_x_LR[4] = u_x_R*(Q_x_R[4][idx] + p_x_R[idx]);
        
        for (int ei = 0; ei < 5; ei++)
        {
            F_x[ei][idx_flux] = F_x_LR[ei] + s_x_plus*(Q_x_star_LR[ei] - Q_x_R[ei][idx]);
        }
    }
}


/*
 * Compute the local convective flux in the y-direction from conservative variables with
 * 2D HLLC Riemann solver.
 */
static inline __attribute__((always_inline)) void
computeLocalConvectiveFluxInYDirectionFromConservativeVariablesHLLC2D(
    double** F_y,
    double** Q_y_B,
    double** Q_y_T,
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
    const int& idx)
{
    v_y_B = Q_y_B[2][idx]/Q_y_B[0][idx];
    v_y_T = Q_y_T[2][idx]/Q_y_T[0][idx];
    
    const double v_y_average = double(1)/double(2)*(v_y_B + v_y_T);
    const double c_y_average = double(1)/double(2)*(c_y_B[idx] + c_y_T[idx]);
    
    const double s_y_B = fmin(v_y_average - c_y_average, v_y_B - c_y_B[idx]);
    const double s_y_T = fmax(v_y_average + c_y_average, v_y_T + c_y_T[idx]);
    
    s_y_minus = fmin(double(0), s_y_B);
    s_y_plus  = fmax(double(0), s_y_T);
    
    s_y_star = (p_y_T[idx] - p_y_B[idx] +
        Q_y_B[2][idx]*(s_y_B - v_y_B) - Q_y_T[2][idx]*(s_y_T - v_y_T))/
        (Q_y_B[0][idx]*(s_y_B - v_y_B) - Q_y_T[0][idx]*(s_y_T - v_y_T));
    
    double Q_y_star_BT[4];
    double F_y_BT[4];
    
    if (s_y_star > double(0))
    {
        Chi_y_star_BT = (s_y_B - v_y_B)/(s_y_B - s_y_star);
        
        Q_y_star_BT[0] = Chi_y_star_BT*Q_y_B[0][idx];
        Q_y_star_BT[1] = Chi_y_star_BT*Q_y_B[1][idx];
        Q_y_star_BT[2] = Chi_y_star_BT*Q_y_B[0][idx]*s_y_star;
        Q_y_star_BT[3] = Chi_y_star_BT*(Q_y_B[3][idx] + (s_y_star - v_y_B)*(Q_y_B[0][idx]*s_y_star +
            p_y_B[idx]/(s_y_B - v_y_B)));
        
        F_y_BT[0] = Q_y_B[2][idx];
        F_y_BT[1] = v_y_B*Q_y_B[1][idx];
        F_y_BT[2] = v_y_B*Q_y_B[2][idx] + p_y_B[idx];
        F_y_BT[3] = v_y_B*(Q_y_B[3][idx] + p_y_B[idx]);
        
        for (int ei = 0; ei < 4; ei++)
        {
            F_y[ei][idx_flux] = F_y_BT[ei] + s_y_minus*(Q_y_star_BT[ei] - Q_y_B[ei][idx]);
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
        
        F_y_BT[0] = Q_y_T[2][idx];
        F_y_BT[1] = v_y_T*Q_y_T[1][idx];
        F_y_BT[2] = v_y_T*Q_y_T[2][idx] + p_y_T[idx];
        F_y_BT[3] = v_y_T*(Q_y_T[3][idx] + p_y_T[idx]);
        
        for (int ei = 0; ei < 4; ei++)
        {
            F_y[ei][idx_flux] = F_y_BT[ei] + s_y_plus*(Q_y_star_BT[ei] - Q_y_T[ei][idx]);
        }
    }
}


/*
 * Compute the local convective flux in the y-direction from conservative variables with
 * 3D HLLC Riemann solver.
 */
static inline __attribute__((always_inline)) void
computeLocalConvectiveFluxInYDirectionFromConservativeVariablesHLLC3D(
    double** F_y,
    double** Q_y_B,
    double** Q_y_T,
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
    const int& idx)
{
    v_y_B = Q_y_B[2][idx]/Q_y_B[0][idx];
    v_y_T = Q_y_T[2][idx]/Q_y_T[0][idx];
    
    const double v_y_average = double(1)/double(2)*(v_y_B + v_y_T);
    const double c_y_average = double(1)/double(2)*(c_y_B[idx] + c_y_T[idx]);
    
    const double s_y_B = fmin(v_y_average - c_y_average, v_y_B - c_y_B[idx]);
    const double s_y_T = fmax(v_y_average + c_y_average, v_y_T + c_y_T[idx]);
    
    s_y_minus = fmin(double(0), s_y_B);
    s_y_plus  = fmax(double(0), s_y_T);
    
    s_y_star = (p_y_T[idx] - p_y_B[idx] +
        Q_y_B[2][idx]*(s_y_B - v_y_B) - Q_y_T[2][idx]*(s_y_T - v_y_T))/
        (Q_y_B[0][idx]*(s_y_B - v_y_B) - Q_y_T[0][idx]*(s_y_T - v_y_T));
    
    double Q_y_star_BT[5];
    double F_y_BT[5];
    
    if (s_y_star > double(0))
    {
        Chi_y_star_BT = (s_y_B - v_y_B)/(s_y_B - s_y_star);
        
        Q_y_star_BT[0] = Chi_y_star_BT*Q_y_B[0][idx];
        Q_y_star_BT[1] = Chi_y_star_BT*Q_y_B[1][idx];
        Q_y_star_BT[2] = Chi_y_star_BT*Q_y_B[0][idx]*s_y_star;
        Q_y_star_BT[3] = Chi_y_star_BT*Q_y_B[3][idx];
        Q_y_star_BT[4] = Chi_y_star_BT*(Q_y_B[4][idx] + (s_y_star - v_y_B)*(Q_y_B[0][idx]*s_y_star +
            p_y_B[idx]/(s_y_B - v_y_B)));
        
        F_y_BT[0] = Q_y_B[2][idx];
        F_y_BT[1] = v_y_B*Q_y_B[1][idx];
        F_y_BT[2] = v_y_B*Q_y_B[2][idx] + p_y_B[idx];
        F_y_BT[3] = v_y_B*Q_y_B[3][idx];
        F_y_BT[4] = v_y_B*(Q_y_B[4][idx] + p_y_B[idx]);
        
        for (int ei = 0; ei < 5; ei++)
        {
            F_y[ei][idx_flux] = F_y_BT[ei] + s_y_minus*(Q_y_star_BT[ei] - Q_y_B[ei][idx]);
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
        
        F_y_BT[0] = Q_y_T[2][idx];
        F_y_BT[1] = v_y_T*Q_y_T[1][idx];
        F_y_BT[2] = v_y_T*Q_y_T[2][idx] + p_y_T[idx];
        F_y_BT[3] = v_y_T*Q_y_T[3][idx];
        F_y_BT[4] = v_y_T*(Q_y_T[4][idx] + p_y_T[idx]);
        
        for (int ei = 0; ei < 5; ei++)
        {
            F_y[ei][idx_flux] = F_y_BT[ei] + s_y_plus*(Q_y_star_BT[ei] - Q_y_T[ei][idx]);
        }
    }
}


/*
 * Compute the local convective flux in the z-direction from conservative variables with
 * 3D HLLC Riemann solver.
 */
static inline __attribute__((always_inline)) void
computeLocalConvectiveFluxInZDirectionFromConservativeVariablesHLLC3D(
    double** F_z,
    double** Q_z_B,
    double** Q_z_F,
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
    const int& idx)
{
    w_z_B = Q_z_B[3][idx]/Q_z_B[0][idx];
    w_z_F = Q_z_F[3][idx]/Q_z_F[0][idx];
    
    const double w_z_average = double(1)/double(2)*(w_z_B + w_z_F);
    const double c_z_average = double(1)/double(2)*(c_z_B[idx] + c_z_F[idx]);
    
    const double s_z_B = fmin(w_z_average - c_z_average, w_z_B - c_z_B[idx]);
    const double s_z_F = fmax(w_z_average + c_z_average, w_z_F + c_z_F[idx]);
    
    s_z_minus = fmin(double(0), s_z_B);
    s_z_plus  = fmax(double(0), s_z_F);
    
    s_z_star = (p_z_F[idx] - p_z_B[idx] +
        Q_z_B[3][idx]*(s_z_B - w_z_B) - Q_z_F[3][idx]*(s_z_F - w_z_F))/
        (Q_z_B[0][idx]*(s_z_B - w_z_B) - Q_z_F[0][idx]*(s_z_F - w_z_F));
    
    double Q_z_star_BF[5];
    double F_z_BF[5];
    
    if (s_z_star > double(0))
    {
        Chi_z_star_BF = (s_z_B - w_z_B)/(s_z_B - s_z_star);
        
        Q_z_star_BF[0] = Chi_z_star_BF*Q_z_B[0][idx];
        Q_z_star_BF[1] = Chi_z_star_BF*Q_z_B[1][idx];
        Q_z_star_BF[2] = Chi_z_star_BF*Q_z_B[2][idx];
        Q_z_star_BF[3] = Chi_z_star_BF*Q_z_B[0][idx]*s_z_star;
        Q_z_star_BF[4] = Chi_z_star_BF*(Q_z_B[4][idx] + (s_z_star - w_z_B)*(Q_z_B[0][idx]*s_z_star +
            p_z_B[idx]/(s_z_B - w_z_B)));
        
        F_z_BF[0] = Q_z_B[3][idx];
        F_z_BF[1] = w_z_B*Q_z_B[1][idx];
        F_z_BF[2] = w_z_B*Q_z_B[2][idx];
        F_z_BF[3] = w_z_B*Q_z_B[3][idx] + p_z_B[idx];
        F_z_BF[4] = w_z_B*(Q_z_B[4][idx] + p_z_B[idx]);
        
        for (int ei = 0; ei < 5; ei++)
        {
            F_z[ei][idx_flux] = F_z_BF[ei] + s_z_minus*(Q_z_star_BF[ei] - Q_z_B[ei][idx]);
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
        
        F_z_BF[0] = Q_z_F[3][idx];
        F_z_BF[1] = w_z_F*Q_z_F[1][idx];
        F_z_BF[2] = w_z_F*Q_z_F[2][idx];
        F_z_BF[3] = w_z_F*Q_z_F[3][idx] + p_z_F[idx];
        F_z_BF[4] = w_z_F*(Q_z_F[4][idx] + p_z_F[idx]);
        
        for (int ei = 0; ei < 5; ei++)
        {
            F_z[ei][idx_flux] = F_z_BF[ei] + s_z_plus*(Q_z_star_BF[ei] - Q_z_F[ei][idx]);
        }
    }
}


/*
 * Compute the local convective flux in the x-direction from primitive variables with
 * 1D HLLC Riemann solver.
 */
static inline __attribute__((always_inline)) void
computeLocalConvectiveFluxInXDirectionFromPrimitiveVariablesHLLC1D(
    double** F_x,
    double** V_x_L,
    double** V_x_R,
    double* c_x_L,
    double* c_x_R,
    double* epsilon_x_L,
    double* epsilon_x_R,
    double& s_x_minus,
    double& s_x_plus,
    double& s_x_star,
    double& Chi_x_star_LR,
    const int& idx_flux,
    const int& idx)
{
    const double u_x_average = double(1)/double(2)*(V_x_L[1][idx] + V_x_R[1][idx]);
    const double c_x_average = double(1)/double(2)*(c_x_L[idx] + c_x_R[idx]);
    
    const double s_x_L = fmin(u_x_average - c_x_average, V_x_L[1][idx] - c_x_L[idx]);
    const double s_x_R = fmax(u_x_average + c_x_average, V_x_R[1][idx] + c_x_R[idx]);
    
    s_x_minus = fmin(double(0), s_x_L);
    s_x_plus  = fmax(double(0), s_x_R);
    
    s_x_star = (V_x_R[2][idx] - V_x_L[2][idx] + V_x_L[0][idx]*V_x_L[1][idx]*(s_x_L - V_x_L[1][idx]) -
        V_x_R[0][idx]*V_x_R[1][idx]*(s_x_R - V_x_R[1][idx]))/
        (V_x_L[0][idx]*(s_x_L - V_x_L[1][idx]) - V_x_R[0][idx]*(s_x_R - V_x_R[1][idx]));
    
    double Q_x_LR[3];
    double Q_x_star_LR[3];
    double F_x_LR[3];
    
    if (s_x_star > double(0))
    {
        Chi_x_star_LR = (s_x_L - V_x_L[1][idx])/(s_x_L - s_x_star);
        
        Q_x_LR[0] = V_x_L[0][idx];
        Q_x_LR[1] = V_x_L[0][idx]*V_x_L[1][idx];
        Q_x_LR[2] = V_x_L[0][idx]*(epsilon_x_L[idx] + double(1)/double(2)*V_x_L[1][idx]*V_x_L[1][idx]);
        
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
        Q_x_LR[2] = V_x_R[0][idx]*(epsilon_x_R[idx] + double(1)/double(2)*V_x_R[1][idx]*V_x_R[1][idx]);
        
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
 * 2D HLLC Riemann solver.
 */
static inline __attribute__((always_inline)) void
computeLocalConvectiveFluxInXDirectionFromPrimitiveVariablesHLLC2D(
    double** F_x,
    double** V_x_L,
    double** V_x_R,
    double* c_x_L,
    double* c_x_R,
    double* epsilon_x_L,
    double* epsilon_x_R,
    double& s_x_minus,
    double& s_x_plus,
    double& s_x_star,
    double& Chi_x_star_LR,
    const int& idx_flux,
    const int& idx)
{
    const double u_x_average = double(1)/double(2)*(V_x_L[1][idx] + V_x_R[1][idx]);
    const double c_x_average = double(1)/double(2)*(c_x_L[idx] + c_x_R[idx]);
    
    const double s_x_L = fmin(u_x_average - c_x_average, V_x_L[1][idx] - c_x_L[idx]);
    const double s_x_R = fmax(u_x_average + c_x_average, V_x_R[1][idx] + c_x_R[idx]);
    
    s_x_minus = fmin(double(0), s_x_L);
    s_x_plus  = fmax(double(0), s_x_R);
    
    s_x_star = (V_x_R[3][idx] - V_x_L[3][idx] + V_x_L[0][idx]*V_x_L[1][idx]*(s_x_L - V_x_L[1][idx]) -
        V_x_R[0][idx]*V_x_R[1][idx]*(s_x_R - V_x_R[1][idx]))/
        (V_x_L[0][idx]*(s_x_L - V_x_L[1][idx]) - V_x_R[0][idx]*(s_x_R - V_x_R[1][idx]));
    
    double Q_x_LR[4];
    double Q_x_star_LR[4];
    double F_x_LR[4];
    
    if (s_x_star > double(0))
    {
        Chi_x_star_LR = (s_x_L - V_x_L[1][idx])/(s_x_L - s_x_star);
        
        Q_x_LR[0] = V_x_L[0][idx];
        Q_x_LR[1] = V_x_L[0][idx]*V_x_L[1][idx];
        Q_x_LR[2] = V_x_L[0][idx]*V_x_L[2][idx];
        Q_x_LR[3] = V_x_L[0][idx]*(epsilon_x_L[idx] + double(1)/double(2)*(V_x_L[1][idx]*V_x_L[1][idx] +
            V_x_L[2][idx]*V_x_L[2][idx]));
        
        Q_x_star_LR[0] = Chi_x_star_LR*V_x_L[0][idx];
        Q_x_star_LR[1] = Chi_x_star_LR*V_x_L[0][idx]*s_x_star;
        Q_x_star_LR[2] = Chi_x_star_LR*Q_x_LR[2];
        Q_x_star_LR[3] = Chi_x_star_LR*(Q_x_LR[3] + (s_x_star - V_x_L[1][idx])*(V_x_L[0][idx]*s_x_star +
            V_x_L[3][idx]/(s_x_L - V_x_L[1][idx])));
        
        F_x_LR[0] = Q_x_LR[1];
        F_x_LR[1] = Q_x_LR[1]*V_x_L[1][idx] + V_x_L[3][idx];
        F_x_LR[2] = Q_x_LR[1]*V_x_L[2][idx];
        F_x_LR[3] = V_x_L[1][idx]*(Q_x_LR[3] + V_x_L[3][idx]);
        
        for (int ei = 0; ei < 4; ei++)
        {
            F_x[ei][idx_flux] = F_x_LR[ei] + s_x_minus*(Q_x_star_LR[ei] - Q_x_LR[ei]);
        }
    }
    else
    {
        Chi_x_star_LR = (s_x_R - V_x_R[1][idx])/(s_x_R - s_x_star);
        
        Q_x_LR[0] = V_x_R[0][idx];
        Q_x_LR[1] = V_x_R[0][idx]*V_x_R[1][idx];
        Q_x_LR[2] = V_x_R[0][idx]*V_x_R[2][idx];
        Q_x_LR[3] = V_x_R[0][idx]*(epsilon_x_R[idx] + double(1)/double(2)*(V_x_R[1][idx]*V_x_R[1][idx] +
            V_x_R[2][idx]*V_x_R[2][idx]));
        
        Q_x_star_LR[0] = Chi_x_star_LR*V_x_R[0][idx];
        Q_x_star_LR[1] = Chi_x_star_LR*V_x_R[0][idx]*s_x_star;
        Q_x_star_LR[2] = Chi_x_star_LR*Q_x_LR[2];
        Q_x_star_LR[3] = Chi_x_star_LR*(Q_x_LR[3] + (s_x_star - V_x_R[1][idx])*(V_x_R[0][idx]*s_x_star +
            V_x_R[3][idx]/(s_x_R - V_x_R[1][idx])));
        
        F_x_LR[0] = Q_x_LR[1];
        F_x_LR[1] = Q_x_LR[1]*V_x_R[1][idx] + V_x_R[3][idx];
        F_x_LR[2] = Q_x_LR[1]*V_x_R[2][idx];
        F_x_LR[3] = V_x_R[1][idx]*(Q_x_LR[3] + V_x_R[3][idx]);
        
        for (int ei = 0; ei < 4; ei++)
        {
            F_x[ei][idx_flux] = F_x_LR[ei] + s_x_plus*(Q_x_star_LR[ei] - Q_x_LR[ei]);
        }
    }
}


/*
 * Compute the local convective flux in the x-direction from primitive variables with
 * 3D HLLC Riemann solver.
 */
static inline __attribute__((always_inline)) void
computeLocalConvectiveFluxInXDirectionFromPrimitiveVariablesHLLC3D(
    double** F_x,
    double** V_x_L,
    double** V_x_R,
    double* c_x_L,
    double* c_x_R,
    double* epsilon_x_L,
    double* epsilon_x_R,
    double& s_x_minus,
    double& s_x_plus,
    double& s_x_star,
    double& Chi_x_star_LR,
    const int& idx_flux,
    const int& idx)
{
    const double u_x_average = double(1)/double(2)*(V_x_L[1][idx] + V_x_R[1][idx]);
    const double c_x_average = double(1)/double(2)*(c_x_L[idx] + c_x_R[idx]);
    
    const double s_x_L = fmin(u_x_average - c_x_average, V_x_L[1][idx] - c_x_L[idx]);
    const double s_x_R = fmax(u_x_average + c_x_average, V_x_R[1][idx] + c_x_R[idx]);
    
    s_x_minus = fmin(double(0), s_x_L);
    s_x_plus  = fmax(double(0), s_x_R);
    
    s_x_star = (V_x_R[4][idx] - V_x_L[4][idx] + V_x_L[0][idx]*V_x_L[1][idx]*(s_x_L - V_x_L[1][idx]) -
        V_x_R[0][idx]*V_x_R[1][idx]*(s_x_R - V_x_R[1][idx]))/
        (V_x_L[0][idx]*(s_x_L - V_x_L[1][idx]) - V_x_R[0][idx]*(s_x_R - V_x_R[1][idx]));
    
    double Q_x_LR[5];
    double Q_x_star_LR[5];
    double F_x_LR[5];
    
    if (s_x_star > double(0))
    {
        Chi_x_star_LR = (s_x_L - V_x_L[1][idx])/(s_x_L - s_x_star);
        
        Q_x_LR[0] = V_x_L[0][idx];
        Q_x_LR[1] = V_x_L[0][idx]*V_x_L[1][idx];
        Q_x_LR[2] = V_x_L[0][idx]*V_x_L[2][idx];
        Q_x_LR[3] = V_x_L[0][idx]*V_x_L[3][idx];
        Q_x_LR[4] = V_x_L[0][idx]*(epsilon_x_L[idx] + double(1)/double(2)*(V_x_L[1][idx]*V_x_L[1][idx] +
            V_x_L[2][idx]*V_x_L[2][idx] + V_x_L[3][idx]*V_x_L[3][idx]));
        
        Q_x_star_LR[0] = Chi_x_star_LR*V_x_L[0][idx];
        Q_x_star_LR[1] = Chi_x_star_LR*V_x_L[0][idx]*s_x_star;
        Q_x_star_LR[2] = Chi_x_star_LR*Q_x_LR[2];
        Q_x_star_LR[3] = Chi_x_star_LR*Q_x_LR[3];
        Q_x_star_LR[4] = Chi_x_star_LR*(Q_x_LR[4] + (s_x_star - V_x_L[1][idx])*(V_x_L[0][idx]*s_x_star +
            V_x_L[4][idx]/(s_x_L - V_x_L[1][idx])));
        
        F_x_LR[0] = Q_x_LR[1];
        F_x_LR[1] = Q_x_LR[1]*V_x_L[1][idx] + V_x_L[4][idx];
        F_x_LR[2] = Q_x_LR[1]*V_x_L[2][idx];
        F_x_LR[3] = Q_x_LR[1]*V_x_L[3][idx];
        F_x_LR[4] = V_x_L[1][idx]*(Q_x_LR[4] + V_x_L[4][idx]);
        
        for (int ei = 0; ei < 5; ei++)
        {
            F_x[ei][idx_flux] = F_x_LR[ei] + s_x_minus*(Q_x_star_LR[ei] - Q_x_LR[ei]);
        }
    }
    else
    {
        Chi_x_star_LR = (s_x_R - V_x_R[1][idx])/(s_x_R - s_x_star);
        
        Q_x_LR[0] = V_x_R[0][idx];
        Q_x_LR[1] = V_x_R[0][idx]*V_x_R[1][idx];
        Q_x_LR[2] = V_x_R[0][idx]*V_x_R[2][idx];
        Q_x_LR[3] = V_x_R[0][idx]*V_x_R[3][idx];
        Q_x_LR[4] = V_x_R[0][idx]*(epsilon_x_R[idx] + double(1)/double(2)*(V_x_R[1][idx]*V_x_R[1][idx] +
            V_x_R[2][idx]*V_x_R[2][idx] + V_x_R[3][idx]*V_x_R[3][idx]));
        
        Q_x_star_LR[0] = Chi_x_star_LR*V_x_R[0][idx];
        Q_x_star_LR[1] = Chi_x_star_LR*V_x_R[0][idx]*s_x_star;
        Q_x_star_LR[2] = Chi_x_star_LR*Q_x_LR[2];
        Q_x_star_LR[3] = Chi_x_star_LR*Q_x_LR[3];
        Q_x_star_LR[4] = Chi_x_star_LR*(Q_x_LR[4] + (s_x_star - V_x_R[1][idx])*(V_x_R[0][idx]*s_x_star +
            V_x_R[4][idx]/(s_x_R - V_x_R[1][idx])));
        
        F_x_LR[0] = Q_x_LR[1];
        F_x_LR[1] = Q_x_LR[1]*V_x_R[1][idx] + V_x_R[4][idx];
        F_x_LR[2] = Q_x_LR[1]*V_x_R[2][idx];
        F_x_LR[3] = Q_x_LR[1]*V_x_R[3][idx];
        F_x_LR[4] = V_x_R[1][idx]*(Q_x_LR[4] + V_x_R[4][idx]);
        
        for (int ei = 0; ei < 5; ei++)
        {
            F_x[ei][idx_flux] = F_x_LR[ei] + s_x_plus*(Q_x_star_LR[ei] - Q_x_LR[ei]);
        }
    }
}


/*
 * Compute the local convective flux in the y-direction from primitive variables with
 * 2D HLLC Riemann solver.
 */
static inline __attribute__((always_inline)) void
computeLocalConvectiveFluxInYDirectionFromPrimitiveVariablesHLLC2D(
    double** F_y,
    double** V_y_B,
    double** V_y_T,
    double* c_y_B,
    double* c_y_T,
    double* epsilon_y_B,
    double* epsilon_y_T,
    double& s_y_minus,
    double& s_y_plus,
    double& s_y_star,
    double& Chi_y_star_BT,
    const int& idx_flux,
    const int& idx)
{
    const double v_y_average = double(1)/double(2)*(V_y_B[2][idx] + V_y_T[2][idx]);
    const double c_y_average = double(1)/double(2)*(c_y_B[idx] + c_y_T[idx]);
    
    const double s_y_B = fmin(v_y_average - c_y_average, V_y_B[2][idx] - c_y_B[idx]);
    const double s_y_T = fmax(v_y_average + c_y_average, V_y_T[2][idx] + c_y_T[idx]);
    
    s_y_minus = fmin(double(0), s_y_B);
    s_y_plus  = fmax(double(0), s_y_T);
    
    s_y_star = (V_y_T[3][idx] - V_y_B[3][idx] + V_y_B[0][idx]*V_y_B[2][idx]*(s_y_B - V_y_B[2][idx]) -
        V_y_T[0][idx]*V_y_T[2][idx]*(s_y_T - V_y_T[2][idx]))/
        (V_y_B[0][idx]*(s_y_B - V_y_B[2][idx]) - V_y_T[0][idx]*(s_y_T - V_y_T[2][idx]));
    
    double Q_y_BT[4];
    double Q_y_star_BT[4];
    double F_y_BT[4];
    
    if (s_y_star > double(0))
    {
        Chi_y_star_BT = (s_y_B - V_y_B[2][idx])/(s_y_B - s_y_star);
        
        Q_y_BT[0] = V_y_B[0][idx];
        Q_y_BT[1] = V_y_B[0][idx]*V_y_B[1][idx];
        Q_y_BT[2] = V_y_B[0][idx]*V_y_B[2][idx];
        Q_y_BT[3] = V_y_B[0][idx]*(epsilon_y_B[idx] + double(1)/double(2)*(V_y_B[1][idx]*V_y_B[1][idx] +
            V_y_B[2][idx]*V_y_B[2][idx]));
        
        Q_y_star_BT[0] = Chi_y_star_BT*V_y_B[0][idx];
        Q_y_star_BT[1] = Chi_y_star_BT*Q_y_BT[1];
        Q_y_star_BT[2] = Chi_y_star_BT*V_y_B[0][idx]*s_y_star;
        Q_y_star_BT[3] = Chi_y_star_BT*(Q_y_BT[3] + (s_y_star - V_y_B[2][idx])*(V_y_B[0][idx]*s_y_star +
            V_y_B[3][idx]/(s_y_B - V_y_B[2][idx])));
        
        F_y_BT[0] = Q_y_BT[2];
        F_y_BT[1] = Q_y_BT[2]*V_y_B[1][idx];
        F_y_BT[2] = Q_y_BT[2]*V_y_B[2][idx] + V_y_B[3][idx];
        F_y_BT[3] = V_y_B[2][idx]*(Q_y_BT[3] + V_y_B[3][idx]);
        
        for (int ei = 0; ei < 4; ei++)
        {
            F_y[ei][idx_flux] = F_y_BT[ei] + s_y_minus*(Q_y_star_BT[ei] - Q_y_BT[ei]);
        }
    }
    else
    {
        Chi_y_star_BT = (s_y_T - V_y_T[2][idx])/(s_y_T - s_y_star);
        
        Q_y_BT[0] = V_y_T[0][idx];
        Q_y_BT[1] = V_y_T[0][idx]*V_y_T[1][idx];
        Q_y_BT[2] = V_y_T[0][idx]*V_y_T[2][idx];
        Q_y_BT[3] = V_y_T[0][idx]*(epsilon_y_T[idx] + double(1)/double(2)*(V_y_T[1][idx]*V_y_T[1][idx] +
            V_y_T[2][idx]*V_y_T[2][idx]));
        
        Q_y_star_BT[0] = Chi_y_star_BT*V_y_T[0][idx];
        Q_y_star_BT[1] = Chi_y_star_BT*Q_y_BT[1];
        Q_y_star_BT[2] = Chi_y_star_BT*V_y_T[0][idx]*s_y_star;
        Q_y_star_BT[3] = Chi_y_star_BT*(Q_y_BT[3] + (s_y_star - V_y_T[2][idx])*(V_y_T[0][idx]*s_y_star +
            V_y_T[3][idx]/(s_y_T - V_y_T[2][idx])));
        
        F_y_BT[0] = Q_y_BT[2];
        F_y_BT[1] = Q_y_BT[2]*V_y_T[1][idx];
        F_y_BT[2] = Q_y_BT[2]*V_y_T[2][idx] + V_y_T[3][idx];
        F_y_BT[3] = V_y_T[2][idx]*(Q_y_BT[3] + V_y_T[3][idx]);
        
        for (int ei = 0; ei < 4; ei++)
        {
            F_y[ei][idx_flux] = F_y_BT[ei] + s_y_plus*(Q_y_star_BT[ei] - Q_y_BT[ei]);
        }
    }
}


/*
 * Compute the local convective flux in the y-direction from primitive variables with
 * 3D HLLC Riemann solver.
 */
static inline __attribute__((always_inline)) void
computeLocalConvectiveFluxInYDirectionFromPrimitiveVariablesHLLC3D(
    double** F_y,
    double** V_y_B,
    double** V_y_T,
    double* c_y_B,
    double* c_y_T,
    double* epsilon_y_B,
    double* epsilon_y_T,
    double& s_y_minus,
    double& s_y_plus,
    double& s_y_star,
    double& Chi_y_star_BT,
    const int& idx_flux,
    const int& idx)
{
    const double v_y_average = double(1)/double(2)*(V_y_B[2][idx] + V_y_T[2][idx]);
    const double c_y_average = double(1)/double(2)*(c_y_B[idx] + c_y_T[idx]);
    
    const double s_y_B = fmin(v_y_average - c_y_average, V_y_B[2][idx] - c_y_B[idx]);
    const double s_y_T = fmax(v_y_average + c_y_average, V_y_T[2][idx] + c_y_T[idx]);
    
    s_y_minus = fmin(double(0), s_y_B);
    s_y_plus  = fmax(double(0), s_y_T);
    
    s_y_star = (V_y_T[4][idx] - V_y_B[4][idx] + V_y_B[0][idx]*V_y_B[2][idx]*(s_y_B - V_y_B[2][idx]) -
        V_y_T[0][idx]*V_y_T[2][idx]*(s_y_T - V_y_T[2][idx]))/
        (V_y_B[0][idx]*(s_y_B - V_y_B[2][idx]) - V_y_T[0][idx]*(s_y_T - V_y_T[2][idx]));
    
    double Q_y_BT[5];
    double Q_y_star_BT[5];
    double F_y_BT[5];
    
    if (s_y_star > double(0))
    {
        Chi_y_star_BT = (s_y_B - V_y_B[2][idx])/(s_y_B - s_y_star);
        
        Q_y_BT[0] = V_y_B[0][idx];
        Q_y_BT[1] = V_y_B[0][idx]*V_y_B[1][idx];
        Q_y_BT[2] = V_y_B[0][idx]*V_y_B[2][idx];
        Q_y_BT[3] = V_y_B[0][idx]*V_y_B[3][idx];
        Q_y_BT[4] = V_y_B[0][idx]*(epsilon_y_B[idx] + double(1)/double(2)*(V_y_B[1][idx]*V_y_B[1][idx] +
            V_y_B[2][idx]*V_y_B[2][idx] + V_y_B[3][idx]*V_y_B[3][idx]));
        
        Q_y_star_BT[0] = Chi_y_star_BT*V_y_B[0][idx];
        Q_y_star_BT[1] = Chi_y_star_BT*Q_y_BT[1];
        Q_y_star_BT[2] = Chi_y_star_BT*V_y_B[0][idx]*s_y_star;
        Q_y_star_BT[3] = Chi_y_star_BT*Q_y_BT[3];
        Q_y_star_BT[4] = Chi_y_star_BT*(Q_y_BT[4] + (s_y_star - V_y_B[2][idx])*(V_y_B[0][idx]*s_y_star +
            V_y_B[4][idx]/(s_y_B - V_y_B[2][idx])));
        
        F_y_BT[0] = Q_y_BT[2];
        F_y_BT[1] = Q_y_BT[2]*V_y_B[1][idx];
        F_y_BT[2] = Q_y_BT[2]*V_y_B[2][idx] + V_y_B[4][idx];
        F_y_BT[3] = Q_y_BT[2]*V_y_B[3][idx];
        F_y_BT[4] = V_y_B[2][idx]*(Q_y_BT[4] + V_y_B[4][idx]);
        
        for (int ei = 0; ei < 5; ei++)
        {
            F_y[ei][idx_flux] = F_y_BT[ei] + s_y_minus*(Q_y_star_BT[ei] - Q_y_BT[ei]);
        }
    }
    else
    {
        Chi_y_star_BT = (s_y_T - V_y_T[2][idx])/(s_y_T - s_y_star);
        
        Q_y_BT[0] = V_y_T[0][idx];
        Q_y_BT[1] = V_y_T[0][idx]*V_y_T[1][idx];
        Q_y_BT[2] = V_y_T[0][idx]*V_y_T[2][idx];
        Q_y_BT[3] = V_y_T[0][idx]*V_y_T[3][idx];
        Q_y_BT[4] = V_y_T[0][idx]*(epsilon_y_T[idx] + double(1)/double(2)*(V_y_T[1][idx]*V_y_T[1][idx] +
            V_y_T[2][idx]*V_y_T[2][idx] + V_y_T[3][idx]*V_y_T[3][idx]));
        
        Q_y_star_BT[0] = Chi_y_star_BT*V_y_T[0][idx];
        Q_y_star_BT[1] = Chi_y_star_BT*Q_y_BT[1];
        Q_y_star_BT[2] = Chi_y_star_BT*V_y_T[0][idx]*s_y_star;
        Q_y_star_BT[3] = Chi_y_star_BT*Q_y_BT[3];
        Q_y_star_BT[4] = Chi_y_star_BT*(Q_y_BT[4] + (s_y_star - V_y_T[2][idx])*(V_y_T[0][idx]*s_y_star +
            V_y_T[4][idx]/(s_y_T - V_y_T[2][idx])));
        
        F_y_BT[0] = Q_y_BT[2];
        F_y_BT[1] = Q_y_BT[2]*V_y_T[1][idx];
        F_y_BT[2] = Q_y_BT[2]*V_y_T[2][idx] + V_y_T[4][idx];
        F_y_BT[3] = Q_y_BT[2]*V_y_T[3][idx];
        F_y_BT[4] = V_y_T[2][idx]*(Q_y_BT[4] + V_y_T[4][idx]);
        
        for (int ei = 0; ei < 5; ei++)
        {
            F_y[ei][idx_flux] = F_y_BT[ei] + s_y_plus*(Q_y_star_BT[ei] - Q_y_BT[ei]);
        }
    }
}


/*
 * Compute the local convective flux in the z-direction from primitive variables with
 * 3D HLLC Riemann solver.
 */
static inline __attribute__((always_inline)) void
computeLocalConvectiveFluxInZDirectionFromPrimitiveVariablesHLLC3D(
    double** F_z,
    double** V_z_B,
    double** V_z_F,
    double* c_z_B,
    double* c_z_F,
    double* epsilon_z_B,
    double* epsilon_z_F,
    double& s_z_minus,
    double& s_z_plus,
    double& s_z_star,
    double& Chi_z_star_BF,
    const int& idx_flux,
    const int& idx)
{
    const double w_z_average = double(1)/double(2)*(V_z_B[3][idx] + V_z_F[3][idx]);
    const double c_z_average = double(1)/double(2)*(c_z_B[idx] + c_z_F[idx]);
    
    const double s_z_B = fmin(w_z_average - c_z_average, V_z_B[3][idx] - c_z_B[idx]);
    const double s_z_F = fmax(w_z_average + c_z_average, V_z_F[3][idx] + c_z_F[idx]);
    
    s_z_minus = fmin(double(0), s_z_B);
    s_z_plus  = fmax(double(0), s_z_F);
    
    s_z_star = (V_z_F[4][idx] - V_z_B[4][idx] + V_z_B[0][idx]*V_z_B[3][idx]*(s_z_B - V_z_B[3][idx]) -
        V_z_F[0][idx]*V_z_F[3][idx]*(s_z_F - V_z_F[3][idx]))/
        (V_z_B[0][idx]*(s_z_B - V_z_B[3][idx]) - V_z_F[0][idx]*(s_z_F - V_z_F[3][idx]));
    
    double Q_z_BF[5];
    double Q_z_star_BF[5];
    double F_z_BF[5];
    
    if (s_z_star > double(0))
    {
        Chi_z_star_BF = (s_z_B - V_z_B[3][idx])/(s_z_B - s_z_star);
        
        Q_z_BF[0] = V_z_B[0][idx];
        Q_z_BF[1] = V_z_B[0][idx]*V_z_B[1][idx];
        Q_z_BF[2] = V_z_B[0][idx]*V_z_B[2][idx];
        Q_z_BF[3] = V_z_B[0][idx]*V_z_B[3][idx];
        Q_z_BF[4] = V_z_B[0][idx]*(epsilon_z_B[idx] + double(1)/double(2)*(V_z_B[1][idx]*V_z_B[1][idx] +
            V_z_B[2][idx]*V_z_B[2][idx] + V_z_B[3][idx]*V_z_B[3][idx]));
        
        Q_z_star_BF[0] = Chi_z_star_BF*V_z_B[0][idx];
        Q_z_star_BF[1] = Chi_z_star_BF*Q_z_BF[1];
        Q_z_star_BF[2] = Chi_z_star_BF*Q_z_BF[2];
        Q_z_star_BF[3] = Chi_z_star_BF*V_z_B[0][idx]*s_z_star;
        Q_z_star_BF[4] = Chi_z_star_BF*(Q_z_BF[4] + (s_z_star - V_z_B[3][idx])*(V_z_B[0][idx]*s_z_star +
            V_z_B[4][idx]/(s_z_B - V_z_B[3][idx])));
        
        F_z_BF[0] = Q_z_BF[3];
        F_z_BF[1] = Q_z_BF[3]*V_z_B[1][idx];
        F_z_BF[2] = Q_z_BF[3]*V_z_B[2][idx];
        F_z_BF[3] = Q_z_BF[3]*V_z_B[3][idx] + V_z_B[4][idx];
        F_z_BF[4] = V_z_B[3][idx]*(Q_z_BF[4] + V_z_B[4][idx]);
        
        for (int ei = 0; ei < 5; ei++)
        {
            F_z[ei][idx_flux] = F_z_BF[ei] + s_z_minus*(Q_z_star_BF[ei] - Q_z_BF[ei]);
        }
    }
    else
    {
        Chi_z_star_BF = (s_z_F - V_z_F[3][idx])/(s_z_F - s_z_star);
        
        Q_z_BF[0] = V_z_F[0][idx];
        Q_z_BF[1] = V_z_F[0][idx]*V_z_F[1][idx];
        Q_z_BF[2] = V_z_F[0][idx]*V_z_F[2][idx];
        Q_z_BF[3] = V_z_F[0][idx]*V_z_F[3][idx];
        Q_z_BF[4] = V_z_F[0][idx]*(epsilon_z_F[idx] + double(1)/double(2)*(V_z_F[1][idx]*V_z_F[1][idx] +
            V_z_F[2][idx]*V_z_F[2][idx] + V_z_F[3][idx]*V_z_F[3][idx]));
        
        Q_z_star_BF[0] = Chi_z_star_BF*V_z_F[0][idx];
        Q_z_star_BF[1] = Chi_z_star_BF*Q_z_BF[1];
        Q_z_star_BF[2] = Chi_z_star_BF*Q_z_BF[2];
        Q_z_star_BF[3] = Chi_z_star_BF*V_z_F[0][idx]*s_z_star;
        Q_z_star_BF[4] = Chi_z_star_BF*(Q_z_BF[4] + (s_z_star - V_z_F[3][idx])*(V_z_F[0][idx]*s_z_star +
            V_z_F[4][idx]/(s_z_F - V_z_F[3][idx])));
        
        F_z_BF[0] = Q_z_BF[3];
        F_z_BF[1] = Q_z_BF[3]*V_z_F[1][idx];
        F_z_BF[2] = Q_z_BF[3]*V_z_F[2][idx];
        F_z_BF[3] = Q_z_BF[3]*V_z_F[3][idx] + V_z_F[4][idx];
        F_z_BF[4] = V_z_F[3][idx]*(Q_z_BF[4] + V_z_F[4][idx]);
        
        for (int ei = 0; ei < 5; ei++)
        {
            F_z[ei][idx_flux] = F_z_BF[ei] + s_z_plus*(Q_z_star_BF[ei] - Q_z_BF[ei]);
        }
    }
}


/*
 * Compute the convective flux and velocity in the x-direction from conservative variables with
 * HLLC Riemann solver.
 */
void
FlowModelRiemannSolverSingleSpecies::computeConvectiveFluxAndVelocityInXDirectionFromConservativeVariablesHLLC(
    HAMERS_SHARED_PTR<pdat::SideData<double> > convective_flux,
    HAMERS_SHARED_PTR<pdat::SideData<double> > velocity,
    const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& conservative_variables_L,
    const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& conservative_variables_R,
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
    
    std::vector<double> thermo_properties;
    std::vector<double*> thermo_properties_ptr;
    std::vector<const double*> thermo_properties_const_ptr;
    
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
    
    HAMERS_SHARED_PTR<pdat::SideData<double> > internal_energy_x_L(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_variables,
            direction_x));
    
    HAMERS_SHARED_PTR<pdat::SideData<double> > internal_energy_x_R(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_variables,
            direction_x));
    
    HAMERS_SHARED_PTR<pdat::SideData<double> > pressure_x_L(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_variables,
            direction_x));
    
    HAMERS_SHARED_PTR<pdat::SideData<double> > pressure_x_R(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_variables,
            direction_x));
    
    HAMERS_SHARED_PTR<pdat::SideData<double> > sound_speed_x_L(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_variables,
            direction_x));
    
    HAMERS_SHARED_PTR<pdat::SideData<double> > sound_speed_x_R(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_variables,
            direction_x));
    
    /*
     * Get the pointers to the temporary data.
     */
    
    double* epsilon_x_L = internal_energy_x_L->getPointer(0, 0);
    double* epsilon_x_R = internal_energy_x_R->getPointer(0, 0);
    
    double* p_x_L = pressure_x_L->getPointer(0, 0);
    double* p_x_R = pressure_x_R->getPointer(0, 0);
    
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
         * Compute the internal energy field.
         */
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
        {
            // Compute the linear index.
            const int idx = i + num_ghosts_0_conservative_variables;
            
            epsilon_x_L[idx] = (Q_x_L[2][idx] -
                double(1)/double(2)*Q_x_L[1][idx]*Q_x_L[1][idx]/Q_x_L[0][idx])/Q_x_L[0][idx];
        }
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
        {
            // Compute the linear index.
            const int idx = i + num_ghosts_0_conservative_variables;
            
            epsilon_x_R[idx] = (Q_x_R[2][idx] -
                double(1)/double(2)*Q_x_R[1][idx]*Q_x_R[1][idx]/Q_x_R[0][idx])/Q_x_R[0][idx];
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
                
                computeLocalConvectiveFluxInXDirectionFromConservativeVariablesHLLC1D(
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
                
                computeLocalConvectiveFluxInXDirectionFromConservativeVariablesHLLC1D(
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
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
            {
                // Compute the linear index.
                const int idx = (i + num_ghosts_0_conservative_variables) +
                    (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables;
                
                epsilon_x_L[idx] = (Q_x_L[3][idx] -
                    double(1)/double(2)*(Q_x_L[1][idx]*Q_x_L[1][idx] + Q_x_L[2][idx]*Q_x_L[2][idx])/
                    Q_x_L[0][idx])/Q_x_L[0][idx];
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
                
                epsilon_x_R[idx] = (Q_x_R[3][idx] -
                    double(1)/double(2)*(Q_x_R[1][idx]*Q_x_R[1][idx] + Q_x_R[2][idx]*Q_x_R[2][idx])/
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
                    
                    computeLocalConvectiveFluxInXDirectionFromConservativeVariablesHLLC2D(
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
                    
                    computeLocalConvectiveFluxInXDirectionFromConservativeVariablesHLLC2D(
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
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                {
                    // Compute the linear indices.
                    const int idx = (i + num_ghosts_0_conservative_variables) +
                        (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables +
                        (k + num_ghosts_2_conservative_variables)*ghostcell_dim_0_conservative_variables*
                            ghostcell_dim_1_conservative_variables;
                    
                    epsilon_x_L[idx] = (Q_x_L[4][idx] -
                        double(1)/double(2)*(Q_x_L[1][idx]*Q_x_L[1][idx] + Q_x_L[2][idx]*Q_x_L[2][idx] +
                        Q_x_L[3][idx]*Q_x_L[3][idx])/Q_x_L[0][idx])/Q_x_L[0][idx];
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
                    // Compute the linear indices.
                    const int idx = (i + num_ghosts_0_conservative_variables) +
                        (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables +
                        (k + num_ghosts_2_conservative_variables)*ghostcell_dim_0_conservative_variables*
                            ghostcell_dim_1_conservative_variables;
                    
                    epsilon_x_R[idx] = (Q_x_R[4][idx] -
                        double(1)/double(2)*(Q_x_R[1][idx]*Q_x_R[1][idx] + Q_x_R[2][idx]*Q_x_R[2][idx] +
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
                        
                        computeLocalConvectiveFluxInXDirectionFromConservativeVariablesHLLC3D(
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
                        
                        computeLocalConvectiveFluxInXDirectionFromConservativeVariablesHLLC3D(
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
 * HLLC Riemann solver.
 */
void
FlowModelRiemannSolverSingleSpecies::computeConvectiveFluxAndVelocityInYDirectionFromConservativeVariablesHLLC(
    HAMERS_SHARED_PTR<pdat::SideData<double> > convective_flux,
    HAMERS_SHARED_PTR<pdat::SideData<double> > velocity,
    const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& conservative_variables_B,
    const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& conservative_variables_T,
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
    
    std::vector<double> thermo_properties;
    std::vector<double*> thermo_properties_ptr;
    std::vector<const double*> thermo_properties_const_ptr;
    
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
    
    HAMERS_SHARED_PTR<pdat::SideData<double> > internal_energy_y_B(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_variables,
            direction_y));
    
    HAMERS_SHARED_PTR<pdat::SideData<double> > internal_energy_y_T(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_variables,
            direction_y));
    
    HAMERS_SHARED_PTR<pdat::SideData<double> > pressure_y_B(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_variables,
            direction_y));
    
    HAMERS_SHARED_PTR<pdat::SideData<double> > pressure_y_T(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_variables,
            direction_y));
    
    HAMERS_SHARED_PTR<pdat::SideData<double> > sound_speed_y_B(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_variables,
            direction_y));
    
    HAMERS_SHARED_PTR<pdat::SideData<double> > sound_speed_y_T(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_variables,
            direction_y));
    
    /*
     * Get the pointers to the temporary data.
     */
    
    double* epsilon_y_B = internal_energy_y_B->getPointer(1, 0);
    double* epsilon_y_T = internal_energy_y_T->getPointer(1, 0);
    
    double* p_y_B = pressure_y_B->getPointer(1, 0);
    double* p_y_T = pressure_y_T->getPointer(1, 0);
    
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
            << ": FlowModelRiemannSolverSingleSpecies::"
            << "computeConvectiveFluxAndVelocityInYDirectionFromConservativeVariablesHLLC()\n"
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
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear index.
                const int idx = (i + num_ghosts_0_conservative_variables) +
                    (j + num_ghosts_1_conservative_variables)*ghostcell_dim_0_conservative_variables;
                
                epsilon_y_B[idx] = (Q_y_B[3][idx] -
                    double(1)/double(2)*(Q_y_B[1][idx]*Q_y_B[1][idx] + Q_y_B[2][idx]*Q_y_B[2][idx])/
                    Q_y_B[0][idx])/Q_y_B[0][idx];
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
                
                epsilon_y_T[idx] = (Q_y_T[3][idx] -
                    double(1)/double(2)*(Q_y_T[1][idx]*Q_y_T[1][idx] + Q_y_T[2][idx]*Q_y_T[2][idx])/
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
                    
                    computeLocalConvectiveFluxInYDirectionFromConservativeVariablesHLLC2D(
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
                    
                    computeLocalConvectiveFluxInYDirectionFromConservativeVariablesHLLC2D(
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
                    
                    epsilon_y_B[idx] = (Q_y_B[4][idx] -
                        double(1)/double(2)*(Q_y_B[1][idx]*Q_y_B[1][idx] + Q_y_B[2][idx]*Q_y_B[2][idx] +
                        Q_y_B[3][idx]*Q_y_B[3][idx])/Q_y_B[0][idx])/Q_y_B[0][idx];
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
                    
                    epsilon_y_T[idx] = (Q_y_T[4][idx] -
                        double(1)/double(2)*(Q_y_T[1][idx]*Q_y_T[1][idx] + Q_y_T[2][idx]*Q_y_T[2][idx] +
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
                        
                        computeLocalConvectiveFluxInYDirectionFromConservativeVariablesHLLC3D(
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
                        
                        computeLocalConvectiveFluxInYDirectionFromConservativeVariablesHLLC3D(
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
 * HLLC Riemann solver.
 */
void
FlowModelRiemannSolverSingleSpecies::computeConvectiveFluxAndVelocityInZDirectionFromConservativeVariablesHLLC(
    HAMERS_SHARED_PTR<pdat::SideData<double> > convective_flux,
    HAMERS_SHARED_PTR<pdat::SideData<double> > velocity,
    const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& conservative_variables_B,
    const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& conservative_variables_F,
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
    
    std::vector<double> thermo_properties;
    std::vector<double*> thermo_properties_ptr;
    std::vector<const double*> thermo_properties_const_ptr;
    
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
    
    HAMERS_SHARED_PTR<pdat::SideData<double> > internal_energy_z_B(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_variables,
            direction_z));
    
    HAMERS_SHARED_PTR<pdat::SideData<double> > internal_energy_z_F(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_variables,
            direction_z));
    
    HAMERS_SHARED_PTR<pdat::SideData<double> > pressure_z_B(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_variables,
            direction_z));
    
    HAMERS_SHARED_PTR<pdat::SideData<double> > pressure_z_F(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_variables,
            direction_z));
    
    HAMERS_SHARED_PTR<pdat::SideData<double> > sound_speed_z_B(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_variables,
            direction_z));
    
    HAMERS_SHARED_PTR<pdat::SideData<double> > sound_speed_z_F(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_variables,
            direction_z));
    
    /*
     * Get the pointers to the temporary data.
     */
    
    double* epsilon_z_B = internal_energy_z_B->getPointer(2, 0);
    double* epsilon_z_F = internal_energy_z_F->getPointer(2, 0);
    
    double* p_z_B = pressure_z_B->getPointer(2, 0);
    double* p_z_F = pressure_z_F->getPointer(2, 0);
    
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
            << ": FlowModelRiemannSolverSingleSpecies::"
            << "computeConvectiveFluxAndVelocityInZDirectionFromConservativeVariablesHLLC()\n"
            << "There is no z direction for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelRiemannSolverSingleSpecies::"
            << "computeConvectiveFluxAndVelocityInZDirectionFromConservativeVariablesHLLC()\n"
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
                    
                    epsilon_z_B[idx] = (Q_z_B[4][idx] -
                        double(1)/double(2)*(Q_z_B[1][idx]*Q_z_B[1][idx] + Q_z_B[2][idx]*Q_z_B[2][idx] +
                        Q_z_B[3][idx]*Q_z_B[3][idx])/Q_z_B[0][idx])/Q_z_B[0][idx];
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
                    
                    epsilon_z_F[idx] = (Q_z_F[4][idx] -
                        double(1)/double(2)*(Q_z_F[1][idx]*Q_z_F[1][idx] + Q_z_F[2][idx]*Q_z_F[2][idx] +
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
                        
                        computeLocalConvectiveFluxInZDirectionFromConservativeVariablesHLLC3D(
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
                        
                        computeLocalConvectiveFluxInZDirectionFromConservativeVariablesHLLC3D(
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
 * HLLC Riemann solver.
 */
void
FlowModelRiemannSolverSingleSpecies::computeConvectiveFluxAndVelocityInXDirectionFromPrimitiveVariablesHLLC(
    HAMERS_SHARED_PTR<pdat::SideData<double> > convective_flux,
    HAMERS_SHARED_PTR<pdat::SideData<double> > velocity,
    const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& primitive_variables_L,
    const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& primitive_variables_R,
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
    
    std::vector<double> thermo_properties;
    std::vector<double*> thermo_properties_ptr;
    std::vector<const double*> thermo_properties_const_ptr;
    
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
    
    HAMERS_SHARED_PTR<pdat::SideData<double> > sound_speed_x_L(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_primitive_variables,
            direction_x));
    
    HAMERS_SHARED_PTR<pdat::SideData<double> > sound_speed_x_R(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_primitive_variables,
            direction_x));
    
    HAMERS_SHARED_PTR<pdat::SideData<double> > internal_energy_x_L(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_primitive_variables,
            direction_x));
    
    HAMERS_SHARED_PTR<pdat::SideData<double> > internal_energy_x_R(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_primitive_variables,
            direction_x));
    
    /*
     * Get the pointers to the temporary data.
     */
    
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
                
                computeLocalConvectiveFluxInXDirectionFromPrimitiveVariablesHLLC1D(
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
                
                if (s_x_star > double(0))
                {
                    u[idx_velocity] = V_x_L[1][idx] + s_x_minus*(Chi_x_star_LR - double(1));
                }
                else
                {
                    u[idx_velocity] = V_x_R[1][idx] + s_x_plus*(Chi_x_star_LR - double(1));
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
                
                computeLocalConvectiveFluxInXDirectionFromPrimitiveVariablesHLLC1D(
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
                    
                    computeLocalConvectiveFluxInXDirectionFromPrimitiveVariablesHLLC2D(
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
                    
                    if (s_x_star > double(0))
                    {
                        u[idx_velocity] = V_x_L[1][idx] + s_x_minus*(Chi_x_star_LR - double(1));
                    }
                    else
                    {
                        u[idx_velocity] = V_x_R[1][idx] + s_x_plus*(Chi_x_star_LR - double(1));
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
                    
                    computeLocalConvectiveFluxInXDirectionFromPrimitiveVariablesHLLC2D(
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
                        
                        computeLocalConvectiveFluxInXDirectionFromPrimitiveVariablesHLLC3D(
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
                        
                        if (s_x_star > double(0))
                        {
                            u[idx_velocity] = V_x_L[1][idx] + s_x_minus*(Chi_x_star_LR - double(1));
                        }
                        else
                        {
                            u[idx_velocity] = V_x_R[1][idx] + s_x_plus*(Chi_x_star_LR - double(1));
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
                        
                        computeLocalConvectiveFluxInXDirectionFromPrimitiveVariablesHLLC3D(
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
 * HLLC Riemann solver.
 */
void
FlowModelRiemannSolverSingleSpecies::computeConvectiveFluxAndVelocityInYDirectionFromPrimitiveVariablesHLLC(
    HAMERS_SHARED_PTR<pdat::SideData<double> > convective_flux,
    HAMERS_SHARED_PTR<pdat::SideData<double> > velocity,
    const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& primitive_variables_B,
    const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& primitive_variables_T,
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
    
    std::vector<double> thermo_properties;
    std::vector<double*> thermo_properties_ptr;
    std::vector<const double*> thermo_properties_const_ptr;
    
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
    
    HAMERS_SHARED_PTR<pdat::SideData<double> > sound_speed_y_B(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_primitive_variables,
            direction_y));
    
    HAMERS_SHARED_PTR<pdat::SideData<double> > sound_speed_y_T(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_primitive_variables,
            direction_y));
    
    HAMERS_SHARED_PTR<pdat::SideData<double> > internal_energy_y_B(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_primitive_variables,
            direction_y));

    HAMERS_SHARED_PTR<pdat::SideData<double> > internal_energy_y_T(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_primitive_variables,
            direction_y));
    
    /*
     * Get the pointers to the temporary data.
     */
    
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
            << ": FlowModelRiemannSolverSingleSpecies::"
            << "computeConvectiveFluxAndVelocityInYDirectionFromPrimitiveVariablesHLLC()\n"
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
                    
                    computeLocalConvectiveFluxInYDirectionFromPrimitiveVariablesHLLC2D(
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
                    
                    if (s_y_star > double(0))
                    {
                        v[idx_velocity] = V_y_B[2][idx] + s_y_minus*(Chi_y_star_BT - double(1));
                    }
                    else
                    {
                        v[idx_velocity] = V_y_T[2][idx] + s_y_plus*(Chi_y_star_BT - double(1));
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
                    
                    computeLocalConvectiveFluxInYDirectionFromPrimitiveVariablesHLLC2D(
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
                        
                        computeLocalConvectiveFluxInYDirectionFromPrimitiveVariablesHLLC3D(
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
                        
                        if (s_y_star > double(0))
                        {
                            v[idx_velocity] = V_y_B[2][idx] + s_y_minus*(Chi_y_star_BT - double(1));
                        }
                        else
                        {
                            v[idx_velocity] = V_y_T[2][idx] + s_y_plus*(Chi_y_star_BT - double(1));
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
                        
                        computeLocalConvectiveFluxInYDirectionFromPrimitiveVariablesHLLC3D(
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
 * HLLC Riemann solver.
 */
void
FlowModelRiemannSolverSingleSpecies::computeConvectiveFluxAndVelocityInZDirectionFromPrimitiveVariablesHLLC(
    HAMERS_SHARED_PTR<pdat::SideData<double> > convective_flux,
    HAMERS_SHARED_PTR<pdat::SideData<double> > velocity,
    const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& primitive_variables_B,
    const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& primitive_variables_F,
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
    
    std::vector<double> thermo_properties;
    std::vector<double*> thermo_properties_ptr;
    std::vector<const double*> thermo_properties_const_ptr;
    
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
    
    HAMERS_SHARED_PTR<pdat::SideData<double> > sound_speed_z_B(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_primitive_variables,
            direction_z));
    
    HAMERS_SHARED_PTR<pdat::SideData<double> > sound_speed_z_F(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_primitive_variables,
            direction_z));
    
    HAMERS_SHARED_PTR<pdat::SideData<double> > internal_energy_z_B(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_primitive_variables,
            direction_z));

    HAMERS_SHARED_PTR<pdat::SideData<double> > internal_energy_z_F(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_primitive_variables,
            direction_z));
    
    /*
     * Get the pointers to the temporary data.
     */
    
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
            << ": FlowModelRiemannSolverSingleSpecies::"
            << "computeConvectiveFluxAndVelocityInZDirectionFromPrimitiveVariablesHLLC()\n"
            << "There is no z direction for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelRiemannSolverSingleSpecies::"
            << "computeConvectiveFluxAndVelocityInZDirectionFromPrimitiveVariablesHLLC()\n"
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
                        
                        computeLocalConvectiveFluxInZDirectionFromPrimitiveVariablesHLLC3D(
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
                        
                        if (s_z_star > double(0))
                        {
                            w[idx_velocity] = V_z_B[3][idx] + s_z_minus*(Chi_z_star_BF - double(1));
                        }
                        else
                        {
                            w[idx_velocity] = V_z_F[3][idx] + s_z_plus*(Chi_z_star_BF - double(1));
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
                        
                        computeLocalConvectiveFluxInZDirectionFromPrimitiveVariablesHLLC3D(
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
