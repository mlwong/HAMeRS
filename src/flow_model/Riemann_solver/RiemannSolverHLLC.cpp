#include "flow_model/Riemann_solver/RiemannSolverHLLC.hpp"

#include <cmath>

/*
 * Compute the flux at the intercell face for single-species flow model from conservative variables.
 */
void
RiemannSolverHLLC::computeIntercellFluxForSingleSpeciesFromConservativeVariables(
    std::vector<double*>& flux_intercell,
    std::vector<double>& conservative_variables_minus,
    std::vector<double>& conservative_variables_plus,
    DIRECTION direction)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(flux_intercell.size()) == d_num_eqn);
    TBOX_ASSERT(static_cast<int>(conservative_variables_minus.size()) == d_num_eqn);
    TBOX_ASSERT(static_cast<int>(conservative_variables_plus.size()) == d_num_eqn);
#endif
    
    switch (direction)
    {
        case X_DIRECTION:
        {
            computeIntercellFluxForSingleSpeciesInXDirectionFromConservativeVariables(
                flux_intercell,
                conservative_variables_minus,
                conservative_variables_plus);
            
            break;
        }
        case Y_DIRECTION:
        {
            computeIntercellFluxForSingleSpeciesInYDirectionFromConservativeVariables(
                flux_intercell,
                conservative_variables_minus,
                conservative_variables_plus);
            
            break;
        }
        case Z_DIRECTION:
        {
            computeIntercellFluxForSingleSpeciesInZDirectionFromConservativeVariables(
                flux_intercell,
                conservative_variables_minus,
                conservative_variables_plus);
            
            break;
        }
        default:
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Unknown direction."
                << std::endl);
        }
    }
}


/*
 * Compute the flux at the intercell face for single-species flow model from primitive variables.
 */
void
RiemannSolverHLLC::computeIntercellFluxForSingleSpeciesFromPrimitiveVariables(
    std::vector<double*>& flux_intercell,
    std::vector<double>& primitive_variables_minus,
    std::vector<double>& primitive_variables_plus,
    DIRECTION direction)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(flux_intercell.size()) == d_num_eqn);
    TBOX_ASSERT(static_cast<int>(primitive_variables_minus.size()) == d_num_eqn);
    TBOX_ASSERT(static_cast<int>(primitive_variables_plus.size()) == d_num_eqn);
#endif
    
    switch (direction)
    {
        case X_DIRECTION:
        {
            computeIntercellFluxForSingleSpeciesInXDirectionFromPrimitiveVariables(
                flux_intercell,
                primitive_variables_minus,
                primitive_variables_plus);
            
            break;
        }
        case Y_DIRECTION:
        {
            computeIntercellFluxForSingleSpeciesInYDirectionFromPrimitiveVariables(
                flux_intercell,
                primitive_variables_minus,
                primitive_variables_plus);
            
            break;
        }
        case Z_DIRECTION:
        {
            computeIntercellFluxForSingleSpeciesInZDirectionFromPrimitiveVariables(
                flux_intercell,
                primitive_variables_minus,
                primitive_variables_plus);
            
            break;
        }
        default:
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Unknown direction."
                << std::endl);
        }
    }
}


/*
 * Compute the fluxes at the intercell face for single-species flow model.
 */
void
RiemannSolverHLLC::computeIntercellFluxForSingleSpecies(
    std::vector<double*>& flux_intercell,
    const double* const density_L,
    const double* const density_R,
    const std::vector<const double*>& momentum_L,
    const std::vector<const double*>& momentum_R,
    const double* const total_energy_L,
    const double* const total_energy_R,
    DIRECTION direction)
{
    if (d_dim == tbox::Dimension(1))
    {
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(static_cast<int>(flux_intercell.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(momentum_L.size()) == 1);
        TBOX_ASSERT(static_cast<int>(momentum_R.size()) == 1);
        
        for (int fi = 0; fi < static_cast<int>(flux_intercell.size()); fi++)
        {
            TBOX_ASSERT(flux_intercell[fi]);
        }
        TBOX_ASSERT(density_L);
        TBOX_ASSERT(density_R);
        for (int mi = 0; mi < static_cast<int>(momentum_L.size()); mi++)
        {
            TBOX_ASSERT(momentum_L[mi]);
        }
        for (int mi = 0; mi < static_cast<int>(momentum_R.size()); mi++)
        {
            TBOX_ASSERT(momentum_R[mi]);
        }
        TBOX_ASSERT(total_energy_L);
        TBOX_ASSERT(total_energy_R);
#endif
        
        switch (direction)
        {
            case X_DIRECTION:
            {
                std::vector<double*> F_x_intercell = flux_intercell;
                
                const double& rho_L = *density_L;
                const double& rho_R = *density_R;
                
                const double& rho_u_L = *(momentum_L[0]);
                const double& rho_u_R = *(momentum_R[0]);
                
                const std::vector<const double*> m_L = momentum_L;
                const std::vector<const double*> m_R = momentum_R;
                
                const double& u_L = rho_u_L/rho_L;
                const double& u_R = rho_u_R/rho_R;
                
                const double& E_L = *total_energy_L;
                const double& E_R = *total_energy_R;
                
                const double p_L = d_equation_of_state->getPressure(
                    &rho_L,
                    m_L,
                    &E_L);
                
                const double p_R = d_equation_of_state->getPressure(
                    &rho_R,
                    m_R,
                    &E_R);
                
                const double c_L = d_equation_of_state->getSoundSpeedWithPressure(
                    &rho_L,
                    &p_L);
                
                const double c_R = d_equation_of_state->getSoundSpeedWithPressure(
                    &rho_R,
                    &p_R);
                
                const double u_average = 0.5*(u_L + u_R);
                const double c_average = 0.5*(c_L + c_R);
                
                const double s_L = fmin(u_average - c_average, u_L - c_L);
                const double s_R = fmax(u_average + c_average, u_R + c_R);
                
                const double s_minus = fmin(0.0, s_L);
                const double s_plus  = fmax(0.0, s_R);
                
                const double s_star =
                    (p_R - p_L + rho_L*u_L*(s_L - u_L) - rho_R*u_R*(s_R - u_R))/(rho_L*(s_L - u_L) -
                        rho_R*(s_R - u_R));
                
                if (s_star > 0)
                {
                    const double Chi_star_L = (s_L - u_L)/(s_L - s_star);
                    
                    std::vector<double> Q_star_L;
                    Q_star_L.push_back(Chi_star_L*rho_L);
                    Q_star_L.push_back(Chi_star_L*rho_L*s_star);
                    Q_star_L.push_back(Chi_star_L*(E_L + (s_star - u_L)*(rho_L*s_star + p_L/(s_L - u_L))));
                    
                    std::vector<const double*> Q_L;
                    Q_L.push_back(&rho_L);
                    Q_L.push_back(m_L[0]);
                    Q_L.push_back(&E_L);
                    
                    std::vector<double> F_x_L;
                    F_x_L.push_back(rho_u_L);
                    F_x_L.push_back(rho_u_L*u_L + p_L);
                    F_x_L.push_back(u_L*(E_L + p_L));
                    
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        *(F_x_intercell[ei]) = F_x_L[ei] + s_minus*(Q_star_L[ei] - *(Q_L[ei]));
                    }
                }
                else
                {
                    const double Chi_star_R = (s_R - u_R)/(s_R - s_star);
                    
                    std::vector<double> Q_star_R;
                    Q_star_R.push_back(Chi_star_R*rho_R);
                    Q_star_R.push_back(Chi_star_R*rho_R*s_star);
                    Q_star_R.push_back(Chi_star_R*(E_R + (s_star - u_R)*(rho_R*s_star + p_R/(s_R - u_R))));
                    
                    std::vector<const double*> Q_R;
                    Q_R.push_back(&rho_R);
                    Q_R.push_back(m_R[0]);
                    Q_R.push_back(&E_R);
                    
                    std::vector<double> F_x_R;
                    F_x_R.push_back(rho_u_R);
                    F_x_R.push_back(rho_u_R*u_R + p_R);
                    F_x_R.push_back(u_R*(E_R + p_R));
                    
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        *(F_x_intercell[ei]) = F_x_R[ei] + s_plus*(Q_star_R[ei] - *(Q_R[ei]));
                    }
                }
                
                break;
            }
            case Y_DIRECTION:
            {
                TBOX_ERROR(d_object_name
                           << ": "
                           << "There is no y direction for 1D problem."
                           << std::endl);
                break;
            }
            case Z_DIRECTION:
            {
                TBOX_ERROR(d_object_name
                           << ": "
                           << "There is no z direction for 1D problem."
                           << std::endl);
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                           << ": "
                           << "Unknown direction."
                           << std::endl);
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(static_cast<int>(flux_intercell.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(momentum_L.size()) == 2);
        TBOX_ASSERT(static_cast<int>(momentum_R.size()) == 2);
        
        for (int fi = 0; fi < static_cast<int>(flux_intercell.size()); fi++)
        {
            TBOX_ASSERT(flux_intercell[fi]);
        }
        TBOX_ASSERT(density_L);
        TBOX_ASSERT(density_R);
        for (int mi = 0; mi < static_cast<int>(momentum_L.size()); mi++)
        {
            TBOX_ASSERT(momentum_L[mi]);
        }
        for (int mi = 0; mi < static_cast<int>(momentum_R.size()); mi++)
        {
            TBOX_ASSERT(momentum_R[mi]);
        }
        TBOX_ASSERT(total_energy_L);
        TBOX_ASSERT(total_energy_R);
#endif
        
        switch (direction)
        {
            case X_DIRECTION:
            {
                std::vector<double*> F_x_intercell = flux_intercell;
                
                const double& rho_L = *density_L;
                const double& rho_R = *density_R;
                
                const double& rho_u_L = *(momentum_L[0]);
                const double& rho_u_R = *(momentum_R[0]);
                
                const double& rho_v_L = *(momentum_L[1]);
                const double& rho_v_R = *(momentum_R[1]);
                
                const double& u_L = rho_u_L/rho_L;
                const double& u_R = rho_u_R/rho_R;
                
                const double& v_L = rho_v_L/rho_L;
                const double& v_R = rho_v_R/rho_R;
                
                const std::vector<const double*> m_L = momentum_L;
                const std::vector<const double*> m_R = momentum_R;
                
                const double& E_L = *total_energy_L;
                const double& E_R = *total_energy_R;
                
                const double p_L = d_equation_of_state->getPressure(
                    &rho_L,
                    m_L,
                    &E_L);
                
                const double p_R = d_equation_of_state->getPressure(
                    &rho_R,
                    m_R,
                    &E_R);
                
                const double c_L = d_equation_of_state->getSoundSpeedWithPressure(
                    &rho_L,
                    &p_L);
                
                const double c_R = d_equation_of_state->getSoundSpeedWithPressure(
                    &rho_R,
                    &p_R);
                
                const double u_average = 0.5*(u_L + u_R);
                const double c_average = 0.5*(c_L + c_R);
                
                const double s_L = fmin(u_average - c_average, u_L - c_L);
                const double s_R = fmax(u_average + c_average, u_R + c_R);
                
                const double s_minus = fmin(0.0, s_L);
                const double s_plus  = fmax(0.0, s_R);
                
                const double s_star =
                    (p_R - p_L + rho_L*u_L*(s_L - u_L) - rho_R*u_R*(s_R - u_R))/(rho_L*(s_L - u_L) -
                        rho_R*(s_R - u_R));
                
                if (s_star > 0)
                {
                    const double Chi_star_L = (s_L - u_L)/(s_L - s_star);
                    
                    std::vector<double> Q_star_L;
                    Q_star_L.push_back(Chi_star_L*rho_L);
                    Q_star_L.push_back(Chi_star_L*rho_L*s_star);
                    Q_star_L.push_back(Chi_star_L*rho_v_L);
                    Q_star_L.push_back(Chi_star_L*(E_L + (s_star - u_L)*(rho_L*s_star + p_L/(s_L - u_L))));
                    
                    std::vector<const double*> Q_L;
                    Q_L.push_back(&rho_L);
                    Q_L.push_back(m_L[0]);
                    Q_L.push_back(m_L[1]);
                    Q_L.push_back(&E_L);
                    
                    std::vector<double> F_x_L;
                    F_x_L.push_back(rho_u_L);
                    F_x_L.push_back(rho_u_L*u_L + p_L);
                    F_x_L.push_back(rho_u_L*v_L);
                    F_x_L.push_back(u_L*(E_L + p_L));
                    
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        *(F_x_intercell[ei]) = F_x_L[ei] + s_minus*(Q_star_L[ei] - *(Q_L[ei]));
                    }
                }
                else
                {
                    const double Chi_star_R = (s_R - u_R)/(s_R - s_star);
                    
                    std::vector<double> Q_star_R;
                    Q_star_R.push_back(Chi_star_R*rho_R);
                    Q_star_R.push_back(Chi_star_R*rho_R*s_star);
                    Q_star_R.push_back(Chi_star_R*rho_v_R);
                    Q_star_R.push_back(Chi_star_R*(E_R + (s_star - u_R)*(rho_R*s_star + p_R/(s_R - u_R))));
                    
                    std::vector<const double*> Q_R;
                    Q_R.push_back(&rho_R);
                    Q_R.push_back(m_R[0]);
                    Q_R.push_back(m_R[1]);
                    Q_R.push_back(&E_R);
                    
                    std::vector<double> F_x_R;
                    F_x_R.push_back(rho_u_R);
                    F_x_R.push_back(rho_u_R*u_R + p_R);
                    F_x_R.push_back(rho_u_R*v_R);
                    F_x_R.push_back(u_R*(E_R + p_R));
                    
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        *(F_x_intercell[ei]) = F_x_R[ei] + s_plus*(Q_star_R[ei] - *(Q_R[ei]));
                    }
                }
                
                break;
            }
            case Y_DIRECTION:
            {
                std::vector<double*> F_y_intercell = flux_intercell;
                
                const double& rho_B = *density_L;
                const double& rho_T = *density_R;
                
                const double& rho_u_B = *(momentum_L[0]);
                const double& rho_u_T = *(momentum_R[0]);
                
                const double& rho_v_B = *(momentum_L[1]);
                const double& rho_v_T = *(momentum_R[1]);
                
                const std::vector<const double*> m_B = momentum_L;
                const std::vector<const double*> m_T = momentum_R;
                
                const double& u_B = rho_u_B/rho_B;
                const double& u_T = rho_u_T/rho_T;
                
                const double& v_B = rho_v_B/rho_B;
                const double& v_T = rho_v_T/rho_T;
                
                const double& E_B = *total_energy_L;
                const double& E_T = *total_energy_R;
                
                const double p_B = d_equation_of_state->getPressure(
                    &rho_B,
                    m_B,
                    &E_B);
                
                const double p_T = d_equation_of_state->getPressure(
                    &rho_T,
                    m_T,
                    &E_T);
                
                const double c_B = d_equation_of_state->getSoundSpeedWithPressure(
                    &rho_B,
                    &p_B);
                
                const double c_T = d_equation_of_state->getSoundSpeedWithPressure(
                    &rho_T,
                    &p_T);
                
                const double v_average = 0.5*(v_B + v_T);
                const double c_average = 0.5*(c_B + c_T);
                
                const double s_B = fmin(v_average - c_average, v_B - c_B);
                const double s_T = fmax(v_average + c_average, v_T + c_T);
                
                const double s_minus = fmin(0.0, s_B);
                const double s_plus  = fmax(0.0, s_T);
                
                const double s_star =
                    (p_T - p_B + rho_B*v_B*(s_B - v_B) - rho_T*v_T*(s_T - v_T))/(rho_B*(s_B - v_B) -
                        rho_T*(s_T - v_T));
                
                if (s_star > 0)
                {
                    const double Chi_star_B = (s_B - v_B)/(s_B - s_star);
                    
                    std::vector<double> Q_star_B;
                    Q_star_B.push_back(Chi_star_B*rho_B);
                    Q_star_B.push_back(Chi_star_B*rho_u_B);
                    Q_star_B.push_back(Chi_star_B*rho_B*s_star);
                    Q_star_B.push_back(Chi_star_B*(E_B + (s_star - v_B)*(rho_B*s_star + p_B/(s_B - v_B))));
                    
                    std::vector<const double*> Q_B;
                    Q_B.push_back(&rho_B);
                    Q_B.push_back(m_B[0]);
                    Q_B.push_back(m_B[1]);
                    Q_B.push_back(&E_B);
                    
                    std::vector<double> F_y_B;
                    F_y_B.push_back(rho_v_B);
                    F_y_B.push_back(rho_v_B*u_B);
                    F_y_B.push_back(rho_v_B*v_B + p_B);
                    F_y_B.push_back(v_B*(E_B + p_B));
                    
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        *(F_y_intercell[ei]) = F_y_B[ei] + s_minus*(Q_star_B[ei] - *(Q_B[ei]));
                    }
                }
                else
                {
                    const double Chi_star_T = (s_T - v_T)/(s_T - s_star);
                    
                    std::vector<double> Q_star_T;
                    Q_star_T.push_back(Chi_star_T*rho_T);
                    Q_star_T.push_back(Chi_star_T*rho_u_T);
                    Q_star_T.push_back(Chi_star_T*rho_T*s_star);
                    Q_star_T.push_back(Chi_star_T*(E_T + (s_star - v_T)*(rho_T*s_star + p_T/(s_T - v_T))));
                    
                    std::vector<const double*> Q_T;
                    Q_T.push_back(&rho_T);
                    Q_T.push_back(m_T[0]);
                    Q_T.push_back(m_T[1]);
                    Q_T.push_back(&E_T);
                    
                    std::vector<double> F_y_T;
                    F_y_T.push_back(rho_v_T);
                    F_y_T.push_back(rho_v_T*u_T);
                    F_y_T.push_back(rho_v_T*v_T + p_T);
                    F_y_T.push_back(v_T*(E_T + p_T));
                    
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        *(F_y_intercell[ei]) = F_y_T[ei] + s_plus*(Q_star_T[ei] - *(Q_T[ei]));
                    }
                }
                
                break;
            }
            case Z_DIRECTION:
            {
                TBOX_ERROR(d_object_name
                           << ": "
                           << "There is no z direction for 1D problem."
                           << std::endl);
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                           << ": "
                           << "Unknown direction."
                           << std::endl);
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(static_cast<int>(flux_intercell.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(momentum_L.size()) == 3);
        TBOX_ASSERT(static_cast<int>(momentum_R.size()) == 3);
        
        for (int fi = 0; fi < static_cast<int>(flux_intercell.size()); fi++)
        {
            TBOX_ASSERT(flux_intercell[fi]);
        }
        TBOX_ASSERT(density_L);
        TBOX_ASSERT(density_R);
        for (int mi = 0; mi < static_cast<int>(momentum_L.size()); mi++)
        {
            TBOX_ASSERT(momentum_L[mi]);
        }
        for (int mi = 0; mi < static_cast<int>(momentum_R.size()); mi++)
        {
            TBOX_ASSERT(momentum_R[mi]);
        }
        TBOX_ASSERT(total_energy_L);
        TBOX_ASSERT(total_energy_R);
#endif
        
        switch (direction)
        {
            case X_DIRECTION:
            {
                std::vector<double*> F_x_intercell = flux_intercell;
                
                const double& rho_L = *density_L;
                const double& rho_R = *density_R;
                
                const double& rho_u_L = *(momentum_L[0]);
                const double& rho_u_R = *(momentum_R[0]);
                
                const double& rho_v_L = *(momentum_L[1]);
                const double& rho_v_R = *(momentum_R[1]);
                
                const double& rho_w_L = *(momentum_L[2]);
                const double& rho_w_R = *(momentum_R[2]);
                
                const std::vector<const double*> m_L = momentum_L;
                const std::vector<const double*> m_R = momentum_R;
                
                const double& u_L = rho_u_L/rho_L;
                const double& u_R = rho_u_R/rho_R;
                
                const double& v_L = rho_v_L/rho_L;
                const double& v_R = rho_v_R/rho_R;
                
                const double& w_L = rho_w_L/rho_L;
                const double& w_R = rho_w_R/rho_R;
                
                const double& E_L = *total_energy_L;
                const double& E_R = *total_energy_R;
                
                const double p_L = d_equation_of_state->getPressure(
                    &rho_L,
                    m_L,
                    &E_L);
                
                const double p_R = d_equation_of_state->getPressure(
                    &rho_R,
                    m_R,
                    &E_R);
                
                const double c_L = d_equation_of_state->getSoundSpeedWithPressure(
                    &rho_L,
                    &p_L);
                
                const double c_R = d_equation_of_state->getSoundSpeedWithPressure(
                    &rho_R,
                    &p_R);
                
                const double u_average = 0.5*(u_L + u_R);
                const double c_average = 0.5*(c_L + c_R);
                
                const double s_L = fmin(u_average - c_average, u_L - c_L);
                const double s_R = fmax(u_average + c_average, u_R + c_R);
                
                const double s_minus = fmin(0.0, s_L);
                const double s_plus  = fmax(0.0, s_R);
                
                const double s_star =
                    (p_R - p_L + rho_L*u_L*(s_L - u_L) - rho_R*u_R*(s_R - u_R))/(rho_L*(s_L - u_L) -
                        rho_R*(s_R - u_R));
                
                if (s_star > 0)
                {
                    const double Chi_star_L = (s_L - u_L)/(s_L - s_star);
                    
                    std::vector<double> Q_star_L;
                    Q_star_L.push_back(Chi_star_L*rho_L);
                    Q_star_L.push_back(Chi_star_L*rho_L*s_star);
                    Q_star_L.push_back(Chi_star_L*rho_v_L);
                    Q_star_L.push_back(Chi_star_L*rho_w_L);
                    Q_star_L.push_back(Chi_star_L*(E_L + (s_star - u_L)*(rho_L*s_star + p_L/(s_L - u_L))));
                    
                    std::vector<const double*> Q_L;
                    Q_L.push_back(&rho_L);
                    Q_L.push_back(m_L[0]);
                    Q_L.push_back(m_L[1]);
                    Q_L.push_back(m_L[2]);
                    Q_L.push_back(&E_L);
                    
                    std::vector<double> F_x_L;
                    F_x_L.push_back(rho_u_L);
                    F_x_L.push_back(rho_u_L*u_L + p_L);
                    F_x_L.push_back(rho_u_L*v_L);
                    F_x_L.push_back(rho_u_L*w_L);
                    F_x_L.push_back(u_L*(E_L + p_L));
                    
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        *(F_x_intercell[ei]) = F_x_L[ei] + s_minus*(Q_star_L[ei] - *(Q_L[ei]));
                    }
                }
                else
                {
                    const double Chi_star_R = (s_R - u_R)/(s_R - s_star);
                    
                    std::vector<double> Q_star_R;
                    Q_star_R.push_back(Chi_star_R*rho_R);
                    Q_star_R.push_back(Chi_star_R*rho_R*s_star);
                    Q_star_R.push_back(Chi_star_R*rho_v_R);
                    Q_star_R.push_back(Chi_star_R*rho_w_R);
                    Q_star_R.push_back(Chi_star_R*(E_R + (s_star - u_R)*(rho_R*s_star + p_R/(s_R - u_R))));
                    
                    std::vector<const double*> Q_R;
                    Q_R.push_back(&rho_R);
                    Q_R.push_back(m_R[0]);
                    Q_R.push_back(m_R[1]);
                    Q_R.push_back(m_R[2]);
                    Q_R.push_back(&E_R);
                    
                    std::vector<double> F_x_R;
                    F_x_R.push_back(rho_u_R);
                    F_x_R.push_back(rho_u_R*u_R + p_R);
                    F_x_R.push_back(rho_u_R*v_R);
                    F_x_R.push_back(rho_u_R*w_R);
                    F_x_R.push_back(u_R*(E_R + p_R));
                    
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        *(F_x_intercell[ei]) = F_x_R[ei] + s_plus*(Q_star_R[ei] - *(Q_R[ei]));
                    }
                }
                
                break;
                
            }
            case Y_DIRECTION:
            {
                std::vector<double*> F_y_intercell = flux_intercell;
                
                const double& rho_B = *density_L;
                const double& rho_T = *density_R;
                
                const double& rho_u_B = *(momentum_L[0]);
                const double& rho_u_T = *(momentum_R[0]);
                
                const double& rho_v_B = *(momentum_L[1]);
                const double& rho_v_T = *(momentum_R[1]);
                
                const double& rho_w_B = *(momentum_L[2]);
                const double& rho_w_T = *(momentum_R[2]);
                
                const std::vector<const double*> m_B = momentum_L;
                const std::vector<const double*> m_T = momentum_R;
                
                const double& u_B = rho_u_B/rho_B;
                const double& u_T = rho_u_T/rho_T;
                
                const double& v_B = rho_v_B/rho_B;
                const double& v_T = rho_v_T/rho_T;
                
                const double& w_B = rho_w_B/rho_B;
                const double& w_T = rho_w_T/rho_T;
                
                const double& E_B = *total_energy_L;
                const double& E_T = *total_energy_R;
                
                const double p_B = d_equation_of_state->getPressure(
                    &rho_B,
                    m_B,
                    &E_B);
                
                const double p_T = d_equation_of_state->getPressure(
                    &rho_T,
                    m_T,
                    &E_T);
                
                const double c_B = d_equation_of_state->getSoundSpeedWithPressure(
                    &rho_B,
                    &p_B);
                
                const double c_T = d_equation_of_state->getSoundSpeedWithPressure(
                    &rho_T,
                    &p_T);
                
                const double v_average = 0.5*(v_B + v_T);
                const double c_average = 0.5*(c_B + c_T);
                
                const double s_B = fmin(v_average - c_average, v_B - c_B);
                const double s_T = fmax(v_average + c_average, v_T + c_T);
                
                const double s_minus = fmin(0.0, s_B);
                const double s_plus  = fmax(0.0, s_T);
                
                const double s_star =
                    (p_T - p_B + rho_B*v_B*(s_B - v_B) - rho_T*v_T*(s_T - v_T))/(rho_B*(s_B - v_B) -
                        rho_T*(s_T - v_T));
                
                if (s_star > 0)
                {
                    const double Chi_star_B = (s_B - v_B)/(s_B - s_star);
                    
                    std::vector<double> Q_star_B;
                    Q_star_B.push_back(Chi_star_B*rho_B);
                    Q_star_B.push_back(Chi_star_B*rho_u_B);
                    Q_star_B.push_back(Chi_star_B*rho_B*s_star);
                    Q_star_B.push_back(Chi_star_B*rho_w_B);
                    Q_star_B.push_back(Chi_star_B*(E_B + (s_star - v_B)*(rho_B*s_star + p_B/(s_B - v_B))));
                    
                    std::vector<const double*> Q_B;
                    Q_B.push_back(&rho_B);
                    Q_B.push_back(m_B[0]);
                    Q_B.push_back(m_B[1]);
                    Q_B.push_back(m_B[2]);
                    Q_B.push_back(&E_B);
                    
                    std::vector<double> F_y_B;
                    F_y_B.push_back(rho_v_B);
                    F_y_B.push_back(rho_v_B*u_B);
                    F_y_B.push_back(rho_v_B*v_B + p_B);
                    F_y_B.push_back(rho_v_B*w_B);
                    F_y_B.push_back(v_B*(E_B + p_B));
                    
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        *(F_y_intercell[ei]) = F_y_B[ei] + s_minus*(Q_star_B[ei] - *(Q_B[ei]));
                    }
                }
                else
                {
                    const double Chi_star_T = (s_T - v_T)/(s_T - s_star);
                    
                    std::vector<double> Q_star_T;
                    Q_star_T.push_back(Chi_star_T*rho_T);
                    Q_star_T.push_back(Chi_star_T*rho_u_T);
                    Q_star_T.push_back(Chi_star_T*rho_T*s_star);
                    Q_star_T.push_back(Chi_star_T*rho_w_T);
                    Q_star_T.push_back(Chi_star_T*(E_T + (s_star - v_T)*(rho_T*s_star + p_T/(s_T - v_T))));
                    
                    std::vector<const double*> Q_T;
                    Q_T.push_back(&rho_T);
                    Q_T.push_back(m_T[0]);
                    Q_T.push_back(m_T[1]);
                    Q_T.push_back(m_T[2]);
                    Q_T.push_back(&E_T);
                    
                    std::vector<double> F_y_T;
                    F_y_T.push_back(rho_v_T);
                    F_y_T.push_back(rho_v_T*u_T);
                    F_y_T.push_back(rho_v_T*v_T + p_T);
                    F_y_T.push_back(rho_v_T*w_T);
                    F_y_T.push_back(v_T*(E_T + p_T));
                    
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        *(F_y_intercell[ei]) = F_y_T[ei] + s_plus*(Q_star_T[ei] - *(Q_T[ei]));
                    }
                }
                
                break;
            }
            case Z_DIRECTION:
            {
                std::vector<double*> F_z_intercell = flux_intercell;
                
                const double& rho_B = *density_L;
                const double& rho_F = *density_R;
                
                const double& rho_u_B = *(momentum_L[0]);
                const double& rho_u_F = *(momentum_R[0]);
                
                const double& rho_v_B = *(momentum_L[1]);
                const double& rho_v_F = *(momentum_R[1]);
                
                const double& rho_w_B = *(momentum_L[2]);
                const double& rho_w_F = *(momentum_R[2]);
                
                const std::vector<const double*> m_B = momentum_L;
                const std::vector<const double*> m_F = momentum_R;
                
                const double& u_B = rho_u_B/rho_B;
                const double& u_F = rho_u_F/rho_F;
                
                const double& v_B = rho_v_B/rho_B;
                const double& v_F = rho_v_F/rho_F;
                
                const double& w_B = rho_w_B/rho_B;
                const double& w_F = rho_w_F/rho_F;
                
                const double& E_B = *total_energy_L;
                const double& E_F = *total_energy_R;
                
                const double p_B = d_equation_of_state->getPressure(
                    &rho_B,
                    m_B,
                    &E_B);
                
                const double p_F = d_equation_of_state->getPressure(
                    &rho_F,
                    m_F,
                    &E_F);
                
                const double c_B = d_equation_of_state->getSoundSpeedWithPressure(
                    &rho_B,
                    &p_B);
                
                const double c_F = d_equation_of_state->getSoundSpeedWithPressure(
                    &rho_F,
                    &p_F);
                
                const double w_average = 0.5*(w_B + w_F);
                const double c_average = 0.5*(c_B + c_F);
                
                const double s_B = fmin(w_average - c_average, w_B - c_B);
                const double s_F = fmax(w_average + c_average, w_F + c_F);
                
                const double s_minus = fmin(0.0, s_B);
                const double s_plus  = fmax(0.0, s_F);
                
                const double s_star =
                    (p_F - p_B + rho_B*w_B*(s_B - w_B) - rho_F*w_F*(s_F - w_F))/(rho_B*(s_B - w_B) -
                        rho_F*(s_F - w_F));
                
                if (s_star > 0)
                {
                    const double Chi_star_B = (s_B - w_B)/(s_B - s_star);
                    
                    std::vector<double> Q_star_B;
                    Q_star_B.push_back(Chi_star_B*rho_B);
                    Q_star_B.push_back(Chi_star_B*rho_u_B);
                    Q_star_B.push_back(Chi_star_B*rho_v_B);
                    Q_star_B.push_back(Chi_star_B*rho_B*s_star);
                    Q_star_B.push_back(Chi_star_B*(E_B + (s_star - w_B)*(rho_B*s_star + p_B/(s_B - w_B))));
                    
                    std::vector<const double*> Q_B;
                    Q_B.push_back(&rho_B);
                    Q_B.push_back(m_B[0]);
                    Q_B.push_back(m_B[1]);
                    Q_B.push_back(m_B[2]);
                    Q_B.push_back(&E_B);
                    
                    std::vector<double> F_z_B;
                    F_z_B.push_back(rho_w_B);
                    F_z_B.push_back(rho_w_B*u_B);
                    F_z_B.push_back(rho_w_B*v_B);
                    F_z_B.push_back(rho_w_B*w_B + p_B);
                    F_z_B.push_back(w_B*(E_B + p_B));
                    
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        *(F_z_intercell[ei]) = F_z_B[ei] + s_minus*(Q_star_B[ei] - *(Q_B[ei]));
                    }
                }
                else
                {
                    const double Chi_star_F = (s_F - w_F)/(s_F - s_star);
                    
                    std::vector<double> Q_star_F;
                    Q_star_F.push_back(Chi_star_F*rho_F);
                    Q_star_F.push_back(Chi_star_F*rho_u_F);
                    Q_star_F.push_back(Chi_star_F*rho_v_F);
                    Q_star_F.push_back(Chi_star_F*rho_F*s_star);
                    Q_star_F.push_back(Chi_star_F*(E_F + (s_star - w_F)*(rho_F*s_star + p_F/(s_F - w_F))));
                    
                    std::vector<const double*> Q_F;
                    Q_F.push_back(&rho_F);
                    Q_F.push_back(m_F[0]);
                    Q_F.push_back(m_F[1]);
                    Q_F.push_back(m_F[2]);
                    Q_F.push_back(&E_F);
                    
                    std::vector<double> F_z_F;
                    F_z_F.push_back(rho_w_F);
                    F_z_F.push_back(rho_w_F*u_F);
                    F_z_F.push_back(rho_w_F*v_F);
                    F_z_F.push_back(rho_w_F*w_F + p_F);
                    F_z_F.push_back(w_F*(E_F + p_F));
                    
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        *(F_z_intercell[ei]) = F_z_F[ei] + s_plus*(Q_star_F[ei] - *(Q_F[ei]));
                    }
                }
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                           << ": "
                           << "Unknown direction."
                           << std::endl);
            }
        }
    }
}


/*
 * Compute the fluxes and velocities at the intercell face
 * for four-equation multi-species flow model by Shyue.
 */
void
RiemannSolverHLLC::computeIntercellFluxAndVelocityForFourEqnShyue(
    std::vector<double*>& flux_intercell,
    std::vector<double*>& velocity_intercell,
    const double* const density_L,
    const double* const density_R,
    const std::vector<const double*>& momentum_L,
    const std::vector<const double*>& momentum_R,
    const double* const total_energy_L,
    const double* const total_energy_R,
    const std::vector<const double*>& mass_fraction_L,
    const std::vector<const double*>& mass_fraction_R,
    DIRECTION direction)
{
    if (d_dim == tbox::Dimension(1))
    {
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(static_cast<int>(flux_intercell.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(velocity_intercell.size()) == 1);
        TBOX_ASSERT(static_cast<int>(momentum_L.size()) == 1);
        TBOX_ASSERT(static_cast<int>(momentum_R.size()) == 1);
        TBOX_ASSERT(static_cast<int>(mass_fraction_L.size()) == d_num_species - 1);
        TBOX_ASSERT(static_cast<int>(mass_fraction_R.size()) == d_num_species - 1);
        
        for (int fi = 0; fi < static_cast<int>(flux_intercell.size()); fi++)
        {
            TBOX_ASSERT(flux_intercell[fi]);
        }
        for (int vi = 0; vi < static_cast<int>(velocity_intercell.size()); vi++)
        {
            TBOX_ASSERT(velocity_intercell[vi]);
        }
        TBOX_ASSERT(density_L);
        TBOX_ASSERT(density_R);
        for (int mi = 0; mi < static_cast<int>(momentum_L.size()); mi++)
        {
            TBOX_ASSERT(momentum_L[mi]);
        }
        for (int mi = 0; mi < static_cast<int>(momentum_R.size()); mi++)
        {
            TBOX_ASSERT(momentum_R[mi]);
        }
        TBOX_ASSERT(total_energy_L);
        TBOX_ASSERT(total_energy_R);
        for (int mi = 0; mi < static_cast<int>(mass_fraction_L.size()); mi++)
        {
            TBOX_ASSERT(mass_fraction_L[mi]);
        }
        for (int mi = 0; mi < static_cast<int>(mass_fraction_R.size()); mi++)
        {
            TBOX_ASSERT(mass_fraction_R[mi]);
        }
#endif
        
        switch (direction)
        {
            case X_DIRECTION:
            {
                std::vector<double*> F_x_intercell = flux_intercell;
                
                double& u_intercell = *(velocity_intercell[0]);
                
                const double& rho_L = *density_L;
                const double& rho_R = *density_R;
                
                const double& rho_u_L = *(momentum_L[0]);
                const double& rho_u_R = *(momentum_R[0]);
                
                const std::vector<const double*> m_L = momentum_L;
                const std::vector<const double*> m_R = momentum_R;
                
                const double& u_L = rho_u_L/rho_L;
                const double& u_R = rho_u_R/rho_R;
                
                const double& E_L = *total_energy_L;
                const double& E_R = *total_energy_R;
                
                const std::vector<const double*> Y_L = mass_fraction_L;
                const std::vector<const double*> Y_R = mass_fraction_R;
                
                const double p_L = d_equation_of_state->
                    getPressureWithMassFraction(
                        &rho_L,
                        m_L,
                        &E_L,
                        Y_L);
                
                const double p_R = d_equation_of_state->
                    getPressureWithMassFraction(
                        &rho_R,
                        m_R,
                        &E_R,
                        Y_R);
                
                const double c_L = d_equation_of_state->
                    getSoundSpeedWithMassFractionAndPressure(
                        &rho_L,
                        Y_L,
                        &p_L);
                
                const double c_R = d_equation_of_state->
                    getSoundSpeedWithMassFractionAndPressure(
                        &rho_R,
                        Y_R,
                        &p_R);
                
                const double u_average = 0.5*(u_L + u_R);
                const double c_average = 0.5*(c_L + c_R);
                
                const double s_L = fmin(u_average - c_average, u_L - c_L);
                const double s_R = fmax(u_average + c_average, u_R + c_R);
                
                const double s_minus = fmin(0.0, s_L);
                const double s_plus  = fmax(0.0, s_R);
                
                const double s_star =
                    (p_R - p_L + rho_L*u_L*(s_L - u_L) - rho_R*u_R*(s_R - u_R))/(rho_L*(s_L - u_L) -
                        rho_R*(s_R - u_R));
                
                if (s_star > 0)
                {
                    const double Chi_star_L = (s_L - u_L)/(s_L - s_star);
                    
                    std::vector<double> Q_star_L;
                    Q_star_L.push_back(Chi_star_L*rho_L);
                    Q_star_L.push_back(Chi_star_L*rho_L*s_star);
                    Q_star_L.push_back(Chi_star_L*(E_L + (s_star - u_L)*(rho_L*s_star + p_L/(s_L - u_L))));
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q_star_L.push_back(Chi_star_L*(*(Y_L[si])));
                    }
                    
                    std::vector<const double*> Q_L;
                    Q_L.push_back(&rho_L);
                    Q_L.push_back(m_L[0]);
                    Q_L.push_back(&E_L);
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q_L.push_back(Y_L[si]);
                    }
                    
                    std::vector<double> F_x_L;
                    F_x_L.push_back(rho_u_L);
                    F_x_L.push_back(rho_u_L*u_L + p_L);
                    F_x_L.push_back(u_L*(E_L + p_L));
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        F_x_L.push_back(u_L*(*(Y_L[si])));
                    }
                    
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        *(F_x_intercell[ei]) = F_x_L[ei] + s_minus*(Q_star_L[ei] - *(Q_L[ei]));
                    }
                    
                    u_intercell = u_L + s_minus*((s_L - u_L)/(s_L - s_star) - 1);
                }
                else
                {
                    const double Chi_star_R = (s_R - u_R)/(s_R - s_star);
                    
                    std::vector<double> Q_star_R;
                    Q_star_R.push_back(Chi_star_R*rho_R);
                    Q_star_R.push_back(Chi_star_R*rho_R*s_star);
                    Q_star_R.push_back(Chi_star_R*(E_R + (s_star - u_R)*(rho_R*s_star + p_R/(s_R - u_R))));
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q_star_R.push_back(Chi_star_R*(*(Y_R[si])));
                    }
                    
                    std::vector<const double*> Q_R;
                    Q_R.push_back(&rho_R);
                    Q_R.push_back(m_R[0]);
                    Q_R.push_back(&E_R);
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q_R.push_back(Y_R[si]);
                    }
                    
                    std::vector<double> F_x_R;
                    F_x_R.push_back(rho_u_R);
                    F_x_R.push_back(rho_u_R*u_R + p_R);
                    F_x_R.push_back(u_R*(E_R + p_R));
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        F_x_R.push_back(u_R*(*(Y_R[si])));
                    }
                    
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        *(F_x_intercell[ei]) = F_x_R[ei] + s_plus*(Q_star_R[ei] - *(Q_R[ei]));
                    }
                    
                    u_intercell = u_R + s_plus*((s_R - u_R)/(s_R - s_star) - 1);
                }
                
                break;  
            }
            case Y_DIRECTION:
            {
                TBOX_ERROR(d_object_name
                           << ": "
                           << "There is no y direction for 1D problem."
                           << std::endl);
                break;
            }
            case Z_DIRECTION:
            {
                TBOX_ERROR(d_object_name
                           << ": "
                           << "There is no z direction for 1D problem."
                           << std::endl);
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                           << ": "
                           << "Unknown direction."
                           << std::endl);
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(static_cast<int>(flux_intercell.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(velocity_intercell.size()) == 2);
        TBOX_ASSERT(static_cast<int>(momentum_L.size()) == 2);
        TBOX_ASSERT(static_cast<int>(momentum_R.size()) == 2);
        TBOX_ASSERT(static_cast<int>(mass_fraction_L.size()) == d_num_species - 1);
        TBOX_ASSERT(static_cast<int>(mass_fraction_R.size()) == d_num_species - 1);
        
        for (int fi = 0; fi < static_cast<int>(flux_intercell.size()); fi++)
        {
            TBOX_ASSERT(flux_intercell[fi]);
        }
        for (int vi = 0; vi < static_cast<int>(velocity_intercell.size()); vi++)
        {
            TBOX_ASSERT(velocity_intercell[vi]);
        }
        TBOX_ASSERT(density_L);
        TBOX_ASSERT(density_R);
        for (int mi = 0; mi < static_cast<int>(momentum_L.size()); mi++)
        {
            TBOX_ASSERT(momentum_L[mi]);
        }
        for (int mi = 0; mi < static_cast<int>(momentum_R.size()); mi++)
        {
            TBOX_ASSERT(momentum_R[mi]);
        }
        TBOX_ASSERT(total_energy_L);
        TBOX_ASSERT(total_energy_R);
        for (int mi = 0; mi < static_cast<int>(mass_fraction_L.size()); mi++)
        {
            TBOX_ASSERT(mass_fraction_L[mi]);
        }
        for (int mi = 0; mi < static_cast<int>(mass_fraction_R.size()); mi++)
        {
            TBOX_ASSERT(mass_fraction_R[mi]);
        }
#endif
        
        switch (direction)
        {
            case X_DIRECTION:
            {
                std::vector<double*> F_x_intercell = flux_intercell;
                
                double& u_intercell = *(velocity_intercell[0]);
                double& v_intercell = *(velocity_intercell[1]);
                
                const double& rho_L = *density_L;
                const double& rho_R = *density_R;
                
                const double& rho_u_L = *(momentum_L[0]);
                const double& rho_u_R = *(momentum_R[0]);
                
                const double& rho_v_L = *(momentum_L[1]);
                const double& rho_v_R = *(momentum_R[1]);
                
                const std::vector<const double*> m_L = momentum_L;
                const std::vector<const double*> m_R = momentum_R;
                
                const double& u_L = rho_u_L/rho_L;
                const double& u_R = rho_u_R/rho_R;
                
                const double& v_L = rho_v_L/rho_L;
                const double& v_R = rho_v_R/rho_R;
                
                const double& E_L = *total_energy_L;
                const double& E_R = *total_energy_R;
                
                const std::vector<const double*> Y_L = mass_fraction_L;
                const std::vector<const double*> Y_R = mass_fraction_R;
                
                const double p_L = d_equation_of_state->
                    getPressureWithMassFraction(
                        &rho_L,
                        m_L,
                        &E_L,
                        Y_L);
                
                const double p_R = d_equation_of_state->
                    getPressureWithMassFraction(
                        &rho_R,
                        m_R,
                        &E_R,
                        Y_R);
                
                const double c_L = d_equation_of_state->
                    getSoundSpeedWithMassFractionAndPressure(
                        &rho_L,
                        Y_L,
                        &p_L);
                
                const double c_R = d_equation_of_state->
                    getSoundSpeedWithMassFractionAndPressure(
                        &rho_R,
                        Y_R,
                        &p_R);
                
                const double u_average = 0.5*(u_L + u_R);
                const double c_average = 0.5*(c_L + c_R);
                
                const double s_L = fmin(u_average - c_average, u_L - c_L);
                const double s_R = fmax(u_average + c_average, u_R + c_R);
                
                const double s_minus = fmin(0.0, s_L);
                const double s_plus  = fmax(0.0, s_R);
                
                const double s_star =
                    (p_R - p_L + rho_L*u_L*(s_L - u_L) - rho_R*u_R*(s_R - u_R))/(rho_L*(s_L - u_L) -
                        rho_R*(s_R - u_R));
                
                if (s_star > 0)
                {
                    const double Chi_star_L = (s_L - u_L)/(s_L - s_star);
                    
                    std::vector<double> Q_star_L;
                    Q_star_L.push_back(Chi_star_L*rho_L);
                    Q_star_L.push_back(Chi_star_L*rho_L*s_star);
                    Q_star_L.push_back(Chi_star_L*rho_v_L);
                    Q_star_L.push_back(Chi_star_L*(E_L + (s_star - u_L)*(rho_L*s_star + p_L/(s_L - u_L))));
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q_star_L.push_back(Chi_star_L*(*(Y_L[si])));
                    }
                    
                    std::vector<const double*> Q_L;
                    Q_L.push_back(&rho_L);
                    Q_L.push_back(m_L[0]);
                    Q_L.push_back(m_L[1]);
                    Q_L.push_back(&E_L);
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q_L.push_back(Y_L[si]);
                    }
                    
                    std::vector<double> F_x_L;
                    F_x_L.push_back(rho_u_L);
                    F_x_L.push_back(rho_u_L*u_L + p_L);
                    F_x_L.push_back(rho_u_L*v_L);
                    F_x_L.push_back(u_L*(E_L + p_L));
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        F_x_L.push_back(u_L*(*(Y_L[si])));
                    }
                    
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        *(F_x_intercell[ei]) = F_x_L[ei] + s_minus*(Q_star_L[ei] - *(Q_L[ei]));
                    }
                    
                    u_intercell = u_L + s_minus*((s_L - u_L)/(s_L - s_star) - 1);
                    v_intercell = v_L;
                }
                else
                {
                    const double Chi_star_R = (s_R - u_R)/(s_R - s_star);
                    
                    std::vector<double> Q_star_R;
                    Q_star_R.push_back(Chi_star_R*rho_R);
                    Q_star_R.push_back(Chi_star_R*rho_R*s_star);
                    Q_star_R.push_back(Chi_star_R*rho_v_R);
                    Q_star_R.push_back(Chi_star_R*(E_R + (s_star - u_R)*(rho_R*s_star + p_R/(s_R - u_R))));
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q_star_R.push_back(Chi_star_R*(*(Y_R[si])));
                    }
                    
                    std::vector<const double*> Q_R;
                    Q_R.push_back(&rho_R);
                    Q_R.push_back(m_R[0]);
                    Q_R.push_back(m_R[1]);
                    Q_R.push_back(&E_R);
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q_R.push_back(Y_R[si]);
                    }
                    
                    std::vector<double> F_x_R;
                    F_x_R.push_back(rho_u_R);
                    F_x_R.push_back(rho_u_R*u_R + p_R);
                    F_x_R.push_back(rho_u_R*v_R);
                    F_x_R.push_back(u_R*(E_R + p_R));
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        F_x_R.push_back(u_R*(*(Y_R[si])));
                    }
                    
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        *(F_x_intercell[ei]) = F_x_R[ei] + s_plus*(Q_star_R[ei] - *(Q_R[ei]));
                    }
                    
                    u_intercell = u_R + s_plus*((s_R - u_R)/(s_R - s_star) - 1);
                    v_intercell = v_R;
                }
                
                break;
            }
            case Y_DIRECTION:
            {
                std::vector<double*> F_y_intercell = flux_intercell;
                
                double& u_intercell = *(velocity_intercell[0]);
                double& v_intercell = *(velocity_intercell[1]);
                
                const double& rho_B = *density_L;
                const double& rho_T = *density_R;
                
                const double& rho_u_B = *(momentum_L[0]);
                const double& rho_u_T = *(momentum_R[0]);
                
                const double& rho_v_B = *(momentum_L[1]);
                const double& rho_v_T = *(momentum_R[1]);
                
                const std::vector<const double*> m_B = momentum_L;
                const std::vector<const double*> m_T = momentum_R;
                
                const double& u_B = rho_u_B/rho_B;
                const double& u_T = rho_u_T/rho_T;
                
                const double& v_B = rho_v_B/rho_B;
                const double& v_T = rho_v_T/rho_T;
                
                const double& E_B = *total_energy_L;
                const double& E_T = *total_energy_R;
                
                const std::vector<const double*> Y_B = mass_fraction_L;
                const std::vector<const double*> Y_T = mass_fraction_R;
                
                const double p_B = d_equation_of_state->
                    getPressureWithMassFraction(
                        &rho_B,
                        m_B,
                        &E_B,
                        Y_B);
                
                const double p_T = d_equation_of_state->
                    getPressureWithMassFraction(
                        &rho_T,
                        m_T,
                        &E_T,
                        Y_T);
                
                const double c_B = d_equation_of_state->
                    getSoundSpeedWithMassFractionAndPressure(
                        &rho_B,
                        Y_B,
                        &p_B);
                
                const double c_T = d_equation_of_state->
                    getSoundSpeedWithMassFractionAndPressure(
                        &rho_T,
                        Y_T,
                        &p_T);
                
                const double v_average = 0.5*(v_B + v_T);
                const double c_average = 0.5*(c_B + c_T);
                
                const double s_B = fmin(v_average - c_average, v_B - c_B);
                const double s_T = fmax(v_average + c_average, v_T + c_T);
                
                const double s_minus = fmin(0.0, s_B);
                const double s_plus  = fmax(0.0, s_T);
                
                const double s_star =
                    (p_T - p_B + rho_B*v_B*(s_B - v_B) - rho_T*v_T*(s_T - v_T))/(rho_B*(s_B - v_B) -
                        rho_T*(s_T - v_T));
                
                if (s_star > 0)
                {
                    const double Chi_star_B = (s_B - v_B)/(s_B - s_star);    
                    
                    std::vector<double> Q_star_B;
                    Q_star_B.push_back(Chi_star_B*rho_B);
                    Q_star_B.push_back(Chi_star_B*rho_u_B);
                    Q_star_B.push_back(Chi_star_B*rho_B*s_star);
                    Q_star_B.push_back(Chi_star_B*(E_B + (s_star - v_B)*(rho_B*s_star + p_B/(s_B - v_B))));
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q_star_B.push_back(Chi_star_B*(*(Y_B[si])));
                    }
                    
                    std::vector<const double*> Q_B;
                    Q_B.push_back(&rho_B);
                    Q_B.push_back(m_B[0]);
                    Q_B.push_back(m_B[1]);
                    Q_B.push_back(&E_B);
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q_B.push_back(Y_B[si]);
                    }
                    
                    std::vector<double> F_y_B;
                    F_y_B.push_back(rho_v_B);
                    F_y_B.push_back(rho_v_B*u_B);
                    F_y_B.push_back(rho_v_B*v_B + p_B);
                    F_y_B.push_back(v_B*(E_B + p_B));
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        F_y_B.push_back(v_B*(*(Y_B[si])));
                    }
                    
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        *(F_y_intercell[ei]) = F_y_B[ei] + s_minus*(Q_star_B[ei] - *(Q_B[ei]));
                    }
                    
                    u_intercell = u_B;
                    v_intercell = v_B + s_minus*((s_B - v_B)/(s_B - s_star) - 1);
                }
                else
                {
                    const double Chi_star_T = (s_T - v_T)/(s_T - s_star); 
                    
                    std::vector<double> Q_star_T;
                    Q_star_T.push_back(Chi_star_T*rho_T);
                    Q_star_T.push_back(Chi_star_T*rho_u_T);
                    Q_star_T.push_back(Chi_star_T*rho_T*s_star);
                    Q_star_T.push_back(Chi_star_T*(E_T + (s_star - v_T)*(rho_T*s_star + p_T/(s_T - v_T))));
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q_star_T.push_back(Chi_star_T*(*(Y_T[si])));
                    }
                    
                    std::vector<const double*> Q_T;
                    Q_T.push_back(&rho_T);
                    Q_T.push_back(m_T[0]);
                    Q_T.push_back(m_T[1]);
                    Q_T.push_back(&E_T);
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q_T.push_back(Y_T[si]);
                    }
                    
                    std::vector<double> F_y_T;
                    F_y_T.push_back(rho_v_T);
                    F_y_T.push_back(rho_v_T*u_T);
                    F_y_T.push_back(rho_v_T*v_T + p_T);
                    F_y_T.push_back(v_T*(E_T + p_T));
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        F_y_T.push_back(v_T*(*(Y_T[si])));
                    }
                    
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        *(F_y_intercell[ei]) = F_y_T[ei] + s_plus*(Q_star_T[ei] - *(Q_T[ei]));
                    }
                    
                    u_intercell = u_T;
                    v_intercell = v_T + s_plus*((s_T - v_T)/(s_T - s_star) - 1);
                }
                
                break;
            }
            case Z_DIRECTION:
            {
                TBOX_ERROR(d_object_name
                           << ": "
                           << "There is no z direction for 1D problem."
                           << std::endl);
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                           << ": "
                           << "Unknown direction."
                           << std::endl);
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(static_cast<int>(flux_intercell.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(velocity_intercell.size()) == 3);
        TBOX_ASSERT(static_cast<int>(momentum_L.size()) == 3);
        TBOX_ASSERT(static_cast<int>(momentum_R.size()) == 3);
        TBOX_ASSERT(static_cast<int>(mass_fraction_L.size()) == d_num_species - 1);
        TBOX_ASSERT(static_cast<int>(mass_fraction_R.size()) == d_num_species - 1);
        
        for (int fi = 0; fi < static_cast<int>(flux_intercell.size()); fi++)
        {
            TBOX_ASSERT(flux_intercell[fi]);
        }
        for (int vi = 0; vi < static_cast<int>(velocity_intercell.size()); vi++)
        {
            TBOX_ASSERT(velocity_intercell[vi]);
        }
        TBOX_ASSERT(density_L);
        TBOX_ASSERT(density_R);
        for (int mi = 0; mi < static_cast<int>(momentum_L.size()); mi++)
        {
            TBOX_ASSERT(momentum_L[mi]);
        }
        for (int mi = 0; mi < static_cast<int>(momentum_R.size()); mi++)
        {
            TBOX_ASSERT(momentum_R[mi]);
        }
        TBOX_ASSERT(total_energy_L);
        TBOX_ASSERT(total_energy_R);
        for (int mi = 0; mi < static_cast<int>(mass_fraction_L.size()); mi++)
        {
            TBOX_ASSERT(mass_fraction_L[mi]);
        }
        for (int mi = 0; mi < static_cast<int>(mass_fraction_R.size()); mi++)
        {
            TBOX_ASSERT(mass_fraction_R[mi]);
        }
#endif
        
        switch (direction)
        {
            case X_DIRECTION:
            {
                std::vector<double*> F_x_intercell = flux_intercell;
                
                double& u_intercell = *(velocity_intercell[0]);
                double& v_intercell = *(velocity_intercell[1]);
                double& w_intercell = *(velocity_intercell[2]);
                
                const double& rho_L = *density_L;
                const double& rho_R = *density_R;
                
                const double& rho_u_L = *(momentum_L[0]);
                const double& rho_u_R = *(momentum_R[0]);
                
                const double& rho_v_L = *(momentum_L[1]);
                const double& rho_v_R = *(momentum_R[1]);
                
                const double& rho_w_L = *(momentum_L[2]);
                const double& rho_w_R = *(momentum_R[2]);
                
                const std::vector<const double*> m_L = momentum_L;
                const std::vector<const double*> m_R = momentum_R;
                
                const double& u_L = rho_u_L/rho_L;
                const double& u_R = rho_u_R/rho_R;
                
                const double& v_L = rho_v_L/rho_L;
                const double& v_R = rho_v_R/rho_R;
                
                const double& w_L = rho_w_L/rho_L;
                const double& w_R = rho_w_R/rho_R;
                
                const double& E_L = *total_energy_L;
                const double& E_R = *total_energy_R;
                
                const std::vector<const double*> Y_L = mass_fraction_L;
                const std::vector<const double*> Y_R = mass_fraction_R;
                
                const double p_L = d_equation_of_state->
                    getPressureWithMassFraction(
                        &rho_L,
                        m_L,
                        &E_L,
                        Y_L);
                
                const double p_R = d_equation_of_state->
                    getPressureWithMassFraction(
                        &rho_R,
                        m_R,
                        &E_R,
                        Y_R);
                
                const double c_L = d_equation_of_state->
                    getSoundSpeedWithMassFractionAndPressure(
                        &rho_L,
                        Y_L,
                        &p_L);
                
                const double c_R = d_equation_of_state->
                    getSoundSpeedWithMassFractionAndPressure(
                        &rho_R,
                        Y_R,
                        &p_R);
                
                const double u_average = 0.5*(u_L + u_R);
                const double c_average = 0.5*(c_L + c_R);
                
                const double s_L = fmin(u_average - c_average, u_L - c_L);
                const double s_R = fmax(u_average + c_average, u_R + c_R);
                
                const double s_minus = fmin(0.0, s_L);
                const double s_plus  = fmax(0.0, s_R);
                
                const double s_star =
                    (p_R - p_L + rho_L*u_L*(s_L - u_L) - rho_R*u_R*(s_R - u_R))/(rho_L*(s_L - u_L) -
                        rho_R*(s_R - u_R));
                
                if (s_star > 0)
                {
                    const double Chi_star_L = (s_L - u_L)/(s_L - s_star);
                    
                    std::vector<double> Q_star_L;
                    Q_star_L.push_back(Chi_star_L*rho_L);
                    Q_star_L.push_back(Chi_star_L*rho_L*s_star);
                    Q_star_L.push_back(Chi_star_L*rho_v_L);
                    Q_star_L.push_back(Chi_star_L*rho_w_L);
                    Q_star_L.push_back(Chi_star_L*(E_L + (s_star - u_L)*(rho_L*s_star + p_L/(s_L - u_L))));
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q_star_L.push_back(Chi_star_L*(*(Y_L[si])));
                    }
                    
                    std::vector<const double*> Q_L;
                    Q_L.push_back(&rho_L);
                    Q_L.push_back(m_L[0]);
                    Q_L.push_back(m_L[1]);
                    Q_L.push_back(m_L[2]);
                    Q_L.push_back(&E_L);
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q_L.push_back(Y_L[si]);
                    }
                    
                    std::vector<double> F_x_L;
                    F_x_L.push_back(rho_u_L);
                    F_x_L.push_back(rho_u_L*u_L + p_L);
                    F_x_L.push_back(rho_u_L*v_L);
                    F_x_L.push_back(rho_u_L*w_L);
                    F_x_L.push_back(u_L*(E_L + p_L));
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        F_x_L.push_back(u_L*(*(Y_L[si])));
                    }
                    
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        *(F_x_intercell[ei]) = F_x_L[ei] + s_minus*(Q_star_L[ei] - *(Q_L[ei]));
                    }
                    
                    u_intercell = u_L + s_minus*((s_L - u_L)/(s_L - s_star) - 1);
                    v_intercell = v_L;
                    w_intercell = w_L;
                }
                else
                {
                    const double Chi_star_R = (s_R - u_R)/(s_R - s_star);
                    
                    std::vector<double> Q_star_R;
                    Q_star_R.push_back(Chi_star_R*rho_R);
                    Q_star_R.push_back(Chi_star_R*rho_R*s_star);
                    Q_star_R.push_back(Chi_star_R*rho_v_R);
                    Q_star_R.push_back(Chi_star_R*rho_w_R);
                    Q_star_R.push_back(Chi_star_R*(E_R + (s_star - u_R)*(rho_R*s_star + p_R/(s_R - u_R))));
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q_star_R.push_back(Chi_star_R*(*(Y_R[si])));
                    }
                    
                    std::vector<const double*> Q_R;
                    Q_R.push_back(&rho_R);
                    Q_R.push_back(m_R[0]);
                    Q_R.push_back(m_R[1]);
                    Q_R.push_back(m_R[2]);
                    Q_R.push_back(&E_R);
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q_R.push_back(Y_R[si]);
                    }
                    
                    std::vector<double> F_x_R;
                    F_x_R.push_back(rho_u_R);
                    F_x_R.push_back(rho_u_R*u_R + p_R);
                    F_x_R.push_back(rho_u_R*v_R);
                    F_x_R.push_back(rho_u_R*w_R);
                    F_x_R.push_back(u_R*(E_R + p_R));
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        F_x_R.push_back(u_R*(*(Y_R[si])));
                    }
                    
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        *(F_x_intercell[ei]) = F_x_R[ei] + s_plus*(Q_star_R[ei] - *(Q_R[ei]));
                    }
                    
                    u_intercell = u_R + s_plus*((s_R - u_R)/(s_R - s_star) - 1);
                    v_intercell = v_R;
                    w_intercell = w_R;
                }
                
                break;
            }
            case Y_DIRECTION:
            {
                std::vector<double*> F_y_intercell = flux_intercell;
                
                double& u_intercell = *(velocity_intercell[0]);
                double& v_intercell = *(velocity_intercell[1]);
                double& w_intercell = *(velocity_intercell[2]);
                
                const double& rho_B = *density_L;
                const double& rho_T = *density_R;
                
                const double& rho_u_B = *(momentum_L[0]);
                const double& rho_u_T = *(momentum_R[0]);
                
                const double& rho_v_B = *(momentum_L[1]);
                const double& rho_v_T = *(momentum_R[1]);
                
                const double& rho_w_B = *(momentum_L[2]);
                const double& rho_w_T = *(momentum_R[2]);
                
                const std::vector<const double*> m_B = momentum_L;
                const std::vector<const double*> m_T = momentum_R;
                
                const double& u_B = rho_u_B/rho_B;
                const double& u_T = rho_u_T/rho_T;
                
                const double& v_B = rho_v_B/rho_B;
                const double& v_T = rho_v_T/rho_T;
                
                const double& w_B = rho_w_B/rho_B;
                const double& w_T = rho_w_T/rho_T;
                
                const double& E_B = *total_energy_L;
                const double& E_T = *total_energy_R;
                
                const std::vector<const double*> Y_B = mass_fraction_L;
                const std::vector<const double*> Y_T = mass_fraction_R;
                
                const double p_B = d_equation_of_state->
                    getPressureWithMassFraction(
                        &rho_B,
                        m_B,
                        &E_B,
                        Y_B);
                
                const double p_T = d_equation_of_state->
                    getPressureWithMassFraction(
                        &rho_T,
                        m_T,
                        &E_T,
                        Y_T);
                
                const double c_B = d_equation_of_state->
                    getSoundSpeedWithMassFractionAndPressure(
                        &rho_B,
                        Y_B,
                        &p_B);
                
                const double c_T = d_equation_of_state->
                    getSoundSpeedWithMassFractionAndPressure(
                        &rho_T,
                        Y_T,
                        &p_T);
                
                const double v_average = 0.5*(v_B + v_T);
                const double c_average = 0.5*(c_B + c_T);
                
                const double s_B = fmin(v_average - c_average, v_B - c_B);
                const double s_T = fmax(v_average + c_average, v_T + c_T);
                
                const double s_minus = fmin(0.0, s_B);
                const double s_plus  = fmax(0.0, s_T);
                
                const double s_star =
                    (p_T - p_B + rho_B*v_B*(s_B - v_B) - rho_T*v_T*(s_T - v_T))/(rho_B*(s_B - v_B) -
                        rho_T*(s_T - v_T));
                
                if (s_star > 0)
                {
                    const double Chi_star_B = (s_B - v_B)/(s_B - s_star);    
                    
                    std::vector<double> Q_star_B;
                    Q_star_B.push_back(Chi_star_B*rho_B);
                    Q_star_B.push_back(Chi_star_B*rho_u_B);
                    Q_star_B.push_back(Chi_star_B*rho_B*s_star);
                    Q_star_B.push_back(Chi_star_B*rho_w_B);
                    Q_star_B.push_back(Chi_star_B*(E_B + (s_star - v_B)*(rho_B*s_star + p_B/(s_B - v_B))));
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q_star_B.push_back(Chi_star_B*(*(Y_B[si])));
                    }
                    
                    std::vector<const double*> Q_B;
                    Q_B.push_back(&rho_B);
                    Q_B.push_back(m_B[0]);
                    Q_B.push_back(m_B[1]);
                    Q_B.push_back(m_B[2]);
                    Q_B.push_back(&E_B);
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q_B.push_back(Y_B[si]);
                    }
                    
                    std::vector<double> F_y_B;
                    F_y_B.push_back(rho_v_B);
                    F_y_B.push_back(rho_v_B*u_B);
                    F_y_B.push_back(rho_v_B*v_B + p_B);
                    F_y_B.push_back(rho_v_B*w_B);
                    F_y_B.push_back(v_B*(E_B + p_B));
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        F_y_B.push_back(v_B*(*(Y_B[si])));
                    }
                    
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        *(F_y_intercell[ei]) = F_y_B[ei] + s_minus*(Q_star_B[ei] - *(Q_B[ei]));
                    }
                    
                    u_intercell = u_B;
                    v_intercell = v_B + s_minus*((s_B - v_B)/(s_B - s_star) - 1);
                    w_intercell = w_B;
                }
                else
                {
                    const double Chi_star_T = (s_T - v_T)/(s_T - s_star); 
                    
                    std::vector<double> Q_star_T;
                    Q_star_T.push_back(Chi_star_T*rho_T);
                    Q_star_T.push_back(Chi_star_T*rho_u_T);
                    Q_star_T.push_back(Chi_star_T*rho_T*s_star);
                    Q_star_T.push_back(Chi_star_T*rho_w_T);
                    Q_star_T.push_back(Chi_star_T*(E_T + (s_star - v_T)*(rho_T*s_star + p_T/(s_T - v_T))));
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q_star_T.push_back(Chi_star_T*(*(Y_T[si])));
                    }
                    
                    std::vector<const double*> Q_T;
                    Q_T.push_back(&rho_T);
                    Q_T.push_back(m_T[0]);
                    Q_T.push_back(m_T[1]);
                    Q_T.push_back(m_T[2]);
                    Q_T.push_back(&E_T);
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q_T.push_back(Y_T[si]);
                    }
                    
                    std::vector<double> F_y_T;
                    F_y_T.push_back(rho_v_T);
                    F_y_T.push_back(rho_v_T*u_T);
                    F_y_T.push_back(rho_v_T*v_T + p_T);
                    F_y_T.push_back(rho_v_T*w_T);
                    F_y_T.push_back(v_T*(E_T + p_T));
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        F_y_T.push_back(v_T*(*(Y_T[si])));
                    }
                    
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        *(F_y_intercell[ei]) = F_y_T[ei] + s_plus*(Q_star_T[ei] - *(Q_T[ei]));
                    }
                    
                    u_intercell = u_T;
                    v_intercell = v_T + s_plus*((s_T - v_T)/(s_T - s_star) - 1);
                    w_intercell = w_T;
                }
                
                break;
            }
            case Z_DIRECTION:
            {
                std::vector<double*> F_z_intercell = flux_intercell;
                
                double& u_intercell = *(velocity_intercell[0]);
                double& v_intercell = *(velocity_intercell[1]);
                double& w_intercell = *(velocity_intercell[2]);
                
                const double& rho_B = *density_L;
                const double& rho_F = *density_R;
                
                const double& rho_u_B = *(momentum_L[0]);
                const double& rho_u_F = *(momentum_R[0]);
                
                const double& rho_v_B = *(momentum_L[1]);
                const double& rho_v_F = *(momentum_R[1]);
                
                const double& rho_w_B = *(momentum_L[2]);
                const double& rho_w_F = *(momentum_R[2]);
                
                const std::vector<const double*> m_B = momentum_L;
                const std::vector<const double*> m_F = momentum_R;
                
                const double& u_B = rho_u_B/rho_B;
                const double& u_F = rho_u_F/rho_F;
                
                const double& v_B = rho_v_B/rho_B;
                const double& v_F = rho_v_F/rho_F;
                
                const double& w_B = rho_w_B/rho_B;
                const double& w_F = rho_w_F/rho_F;
                
                const double& E_B = *total_energy_L;
                const double& E_F = *total_energy_R;
                
                const std::vector<const double*> Y_B = mass_fraction_L;
                const std::vector<const double*> Y_F = mass_fraction_R;
                
                const double p_B = d_equation_of_state->
                    getPressureWithMassFraction(
                        &rho_B,
                        m_B,
                        &E_B,
                        Y_B);
                
                const double p_F = d_equation_of_state->
                    getPressureWithMassFraction(
                        &rho_F,
                        m_F,
                        &E_F,
                        Y_F);
                
                const double c_B = d_equation_of_state->
                    getSoundSpeedWithMassFractionAndPressure(
                        &rho_B,
                        Y_B,
                        &p_B);
                
                const double c_F = d_equation_of_state->
                    getSoundSpeedWithMassFractionAndPressure(
                        &rho_F,
                        Y_F,
                        &p_F);
                
                const double w_average = 0.5*(w_B + w_F);
                const double c_average = 0.5*(c_F + c_F);
                
                const double s_B = fmin(w_average - c_average, w_B - c_B);
                const double s_F = fmax(w_average + c_average, w_F + c_F);
                
                const double s_minus = fmin(0.0, s_B);
                const double s_plus  = fmax(0.0, s_F);
                
                const double s_star =
                    (p_F - p_B + rho_B*w_B*(s_B - w_B) - rho_F*w_F*(s_F - w_F))/(rho_B*(s_B - w_B) -
                        rho_F*(s_F - w_F));
                
                if (s_star > 0)
                {
                    const double Chi_star_B = (s_B - w_B)/(s_B - s_star);    
                    
                    std::vector<double> Q_star_B;
                    Q_star_B.push_back(Chi_star_B*rho_B);
                    Q_star_B.push_back(Chi_star_B*rho_u_B);
                    Q_star_B.push_back(Chi_star_B*rho_v_B);
                    Q_star_B.push_back(Chi_star_B*rho_B*s_star);
                    Q_star_B.push_back(Chi_star_B*(E_B + (s_star - w_B)*(rho_B*s_star + p_B/(s_B - w_B))));
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q_star_B.push_back(Chi_star_B*(*(Y_B[si])));
                    }
                    
                    std::vector<const double*> Q_B;
                    Q_B.push_back(&rho_B);
                    Q_B.push_back(m_B[0]);
                    Q_B.push_back(m_B[1]);
                    Q_B.push_back(m_B[2]);
                    Q_B.push_back(&E_B);
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q_B.push_back(Y_B[si]);
                    }
                    
                    std::vector<double> F_z_B;
                    F_z_B.push_back(rho_w_B);
                    F_z_B.push_back(rho_w_B*u_B);
                    F_z_B.push_back(rho_w_B*v_B);
                    F_z_B.push_back(rho_w_B*w_B + p_B);
                    F_z_B.push_back(w_B*(E_B + p_B));
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        F_z_B.push_back(w_B*(*(Y_B[si])));
                    }
                    
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        *(F_z_intercell[ei]) = F_z_B[ei] + s_minus*(Q_star_B[ei] - *(Q_B[ei]));
                    }
                    
                    u_intercell = u_B;
                    v_intercell = v_B;
                    w_intercell = w_B + s_minus*((s_B - w_B)/(s_B - s_star) - 1);
                }
                else
                {
                    const double Chi_star_F = (s_F - w_F)/(s_F - s_star); 
                    
                    std::vector<double> Q_star_F;
                    Q_star_F.push_back(Chi_star_F*rho_F);
                    Q_star_F.push_back(Chi_star_F*rho_u_F);
                    Q_star_F.push_back(Chi_star_F*rho_v_F);
                    Q_star_F.push_back(Chi_star_F*rho_F*s_star);
                    Q_star_F.push_back(Chi_star_F*(E_F + (s_star - w_F)*(rho_F*s_star + p_F/(s_F - w_F))));
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q_star_F.push_back(Chi_star_F*(*(Y_F[si])));
                    }
                    
                    std::vector<const double*> Q_F;
                    Q_F.push_back(&rho_F);
                    Q_F.push_back(m_F[0]);
                    Q_F.push_back(m_F[1]);
                    Q_F.push_back(m_F[2]);
                    Q_F.push_back(&E_F);
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q_F.push_back(Y_F[si]);
                    }
                    
                    std::vector<double> F_z_F;
                    F_z_F.push_back(rho_w_F);
                    F_z_F.push_back(rho_w_F*u_F);
                    F_z_F.push_back(rho_w_F*v_F);
                    F_z_F.push_back(rho_w_F*w_F + p_F);
                    F_z_F.push_back(w_F*(E_F + p_F));
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        F_z_F.push_back(w_F*(*(Y_F[si])));
                    }
                    
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        *(F_z_intercell[ei]) = F_z_F[ei] + s_plus*(Q_star_F[ei] - *(Q_F[ei]));
                    }
                    
                    u_intercell = u_F;
                    v_intercell = v_F;
                    w_intercell = w_F + s_plus*((s_F - w_F)/(s_F - s_star) - 1);
                }
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                           << ": "
                           << "Unknown direction."
                           << std::endl);
            }
        }
    }
}


/*
 * Compute the fluxes and velocities at the intercell face
 * for five-equation multi-species flow model by Allaire.
 */
void
RiemannSolverHLLC::computeIntercellFluxAndVelocityForFiveEqnAllaire(
    std::vector<double*>& flux_intercell,
    std::vector<double*>& velocity_intercell,
    const std::vector<const double*>& partial_density_L,
    const std::vector<const double*>& partial_density_R,
    const std::vector<const double*>& momentum_L,
    const std::vector<const double*>& momentum_R,
    const double* const total_energy_L,
    const double* const total_energy_R,
    const std::vector<const double*>& volume_fraction_L,
    const std::vector<const double*>& volume_fraction_R,
    DIRECTION direction)
{
    if (d_dim == tbox::Dimension(1))
    {
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(static_cast<int>(flux_intercell.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(velocity_intercell.size()) == 1);
        TBOX_ASSERT(static_cast<int>(partial_density_L.size()) == d_num_species);
        TBOX_ASSERT(static_cast<int>(partial_density_R.size()) == d_num_species);
        TBOX_ASSERT(static_cast<int>(momentum_L.size()) == 1);
        TBOX_ASSERT(static_cast<int>(momentum_R.size()) == 1);
        TBOX_ASSERT(static_cast<int>(volume_fraction_L.size()) == d_num_species - 1);
        TBOX_ASSERT(static_cast<int>(volume_fraction_R.size()) == d_num_species - 1);
        
        for (int fi = 0; fi < static_cast<int>(flux_intercell.size()); fi++)
        {
            TBOX_ASSERT(flux_intercell[fi]);
        }
        for (int vi = 0; vi < static_cast<int>(velocity_intercell.size()); vi++)
        {
            TBOX_ASSERT(velocity_intercell[vi]);
        }
        for (int pi = 0; pi < static_cast<int>(partial_density_L.size()); pi++)
        {
            TBOX_ASSERT(partial_density_L[pi]);
        }
        for (int pi = 0; pi < static_cast<int>(partial_density_R.size()); pi++)
        {
            TBOX_ASSERT(partial_density_R[pi]);
        }
        for (int mi = 0; mi < static_cast<int>(momentum_L.size()); mi++)
        {
            TBOX_ASSERT(momentum_L[mi]);
        }
        for (int mi = 0; mi < static_cast<int>(momentum_R.size()); mi++)
        {
            TBOX_ASSERT(momentum_R[mi]);
        }
        TBOX_ASSERT(total_energy_L);
        TBOX_ASSERT(total_energy_R);
        for (int vi = 0; vi < static_cast<int>(volume_fraction_L.size()); vi++)
        {
            TBOX_ASSERT(volume_fraction_L[vi]);
        }
        for (int vi = 0; vi < static_cast<int>(volume_fraction_R.size()); vi++)
        {
            TBOX_ASSERT(volume_fraction_R[vi]);
        }
#endif
        
        switch (direction)
        {
            case X_DIRECTION:
            {
                std::vector<double*> F_x_intercell = flux_intercell;
                
                double& u_intercell = *(velocity_intercell[0]);
                
                const std::vector<const double*> Z_rho_L = partial_density_L;
                const std::vector<const double*> Z_rho_R = partial_density_R;
                
                const double rho_L = d_equation_of_state->getTotalDensity(
                    partial_density_L);
                
                const double rho_R = d_equation_of_state->getTotalDensity(
                    partial_density_R);
                
                const double& rho_u_L = *(momentum_L[0]);
                const double& rho_u_R = *(momentum_R[0]);
                
                const std::vector<const double*> m_L = momentum_L;
                const std::vector<const double*> m_R = momentum_R;
                
                const double& u_L = rho_u_L/rho_L;
                const double& u_R = rho_u_R/rho_R;
                
                const double& E_L = *total_energy_L;
                const double& E_R = *total_energy_R;
                
                const std::vector<const double*> Z_L = volume_fraction_L;
                const std::vector<const double*> Z_R = volume_fraction_R;
                
                const double p_L = d_equation_of_state->getPressureWithVolumeFraction(
                    &rho_L,
                    m_L,
                    &E_L,
                    Z_L);
                
                const double p_R = d_equation_of_state->getPressureWithVolumeFraction(
                    &rho_R,
                    m_R,
                    &E_R,
                    Z_R);
                
                const double c_L = d_equation_of_state->getSoundSpeedWithVolumeFractionAndPressure(
                    &rho_L,
                    Z_L,
                    &p_L);
                
                const double c_R = d_equation_of_state->getSoundSpeedWithVolumeFractionAndPressure(
                    &rho_R,
                    Z_R,
                    &p_R);
                
                const double u_average = 0.5*(u_L + u_R);
                const double c_average = 0.5*(c_L + c_R);
                
                const double s_L = fmin(u_average - c_average, u_L - c_L);
                const double s_R = fmax(u_average + c_average, u_R + c_R);
                
                const double s_minus = fmin(0.0, s_L);
                const double s_plus  = fmax(0.0, s_R);
                
                const double s_star =
                    (p_R - p_L + rho_L*u_L*(s_L - u_L) - rho_R*u_R*(s_R - u_R))/(rho_L*(s_L - u_L) -
                        rho_R*(s_R - u_R));
                
                if (s_star > 0)
                {
                    const double Chi_star_L = (s_L - u_L)/(s_L - s_star);
                    
                    std::vector<double> Q_star_L;
                    for (int si = 0; si < d_num_species; si++)
                    {
                        Q_star_L.push_back(Chi_star_L*(*(Z_rho_L[si])));
                    }
                    Q_star_L.push_back(Chi_star_L*rho_L*s_star);
                    Q_star_L.push_back(Chi_star_L*(E_L + (s_star - u_L)*(rho_L*s_star + p_L/(s_L - u_L))));
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q_star_L.push_back(Chi_star_L*(*(Z_L[si])));
                    }
                    
                    std::vector<const double*> Q_L;
                    for (int si = 0; si < d_num_species; si++)
                    {
                        Q_L.push_back(Z_rho_L[si]);
                    }
                    Q_L.push_back(m_L[0]);
                    Q_L.push_back(&E_L);
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q_L.push_back(Z_L[si]);
                    }
                    
                    std::vector<double> F_x_L;
                    for (int si = 0; si < d_num_species; si++)
                    {
                        F_x_L.push_back(u_L*(*(Z_rho_L[si])));
                    }
                    F_x_L.push_back(rho_u_L*u_L + p_L);
                    F_x_L.push_back(u_L*(E_L + p_L));
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        F_x_L.push_back(u_L*(*(Z_L[si])));
                    }
                    
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        *(F_x_intercell[ei]) = F_x_L[ei] + s_minus*(Q_star_L[ei] - *(Q_L[ei]));
                    }
                    
                    u_intercell = u_L + s_minus*((s_L - u_L)/(s_L - s_star) - 1);
                }
                else
                {
                    const double Chi_star_R = (s_R - u_R)/(s_R - s_star);
                    
                    std::vector<double> Q_star_R;
                    for (int si = 0; si < d_num_species; si++)
                    {
                        Q_star_R.push_back(Chi_star_R*(*(Z_rho_R[si])));
                    }
                    Q_star_R.push_back(Chi_star_R*rho_R*s_star);
                    Q_star_R.push_back(Chi_star_R*(E_R + (s_star - u_R)*(rho_R*s_star + p_R/(s_R - u_R))));
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q_star_R.push_back(Chi_star_R*(*(Z_R[si])));
                    }
                    
                    std::vector<const double*> Q_R;
                    for (int si = 0; si < d_num_species; si++)
                    {
                        Q_R.push_back(Z_rho_R[si]);
                    }
                    Q_R.push_back(m_R[0]);
                    Q_R.push_back(&E_R);
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q_R.push_back(Z_R[si]);
                    }
                    
                    std::vector<double> F_x_R;
                    for (int si = 0; si < d_num_species; si++)
                    {
                        F_x_R.push_back(u_R*(*(Z_rho_R[si])));
                    }
                    F_x_R.push_back(rho_u_R*u_R + p_R);
                    F_x_R.push_back(u_R*(E_R + p_R));
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        F_x_R.push_back(u_R*(*(Z_R[si])));
                    }
                    
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        *(F_x_intercell[ei]) = F_x_R[ei] + s_plus*(Q_star_R[ei] - *(Q_R[ei]));
                    }
                    
                    u_intercell = u_R + s_plus*((s_R - u_R)/(s_R - s_star) - 1);
                }
                
                break;
            }
            case Y_DIRECTION:
            {
                TBOX_ERROR(d_object_name
                           << ": "
                           << "There is no y direction for 1D problem."
                           << std::endl);
                break;
            }
            case Z_DIRECTION:
            {
                TBOX_ERROR(d_object_name
                           << ": "
                           << "There is no z direction for 1D problem."
                           << std::endl);
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                           << ": "
                           << "Unknown direction."
                           << std::endl);
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(static_cast<int>(flux_intercell.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(velocity_intercell.size()) == 2);
        TBOX_ASSERT(static_cast<int>(partial_density_L.size()) == d_num_species);
        TBOX_ASSERT(static_cast<int>(partial_density_R.size()) == d_num_species);
        TBOX_ASSERT(static_cast<int>(momentum_L.size()) == 2);
        TBOX_ASSERT(static_cast<int>(momentum_R.size()) == 2);
        TBOX_ASSERT(static_cast<int>(volume_fraction_L.size()) == d_num_species - 1);
        TBOX_ASSERT(static_cast<int>(volume_fraction_R.size()) == d_num_species - 1);
        
        for (int fi = 0; fi < static_cast<int>(flux_intercell.size()); fi++)
        {
            TBOX_ASSERT(flux_intercell[fi]);
        }
        for (int vi = 0; vi < static_cast<int>(velocity_intercell.size()); vi++)
        {
            TBOX_ASSERT(velocity_intercell[vi]);
        }
        for (int pi = 0; pi < static_cast<int>(partial_density_L.size()); pi++)
        {
            TBOX_ASSERT(partial_density_L[pi]);
        }
        for (int pi = 0; pi < static_cast<int>(partial_density_R.size()); pi++)
        {
            TBOX_ASSERT(partial_density_R[pi]);
        }
        for (int mi = 0; mi < static_cast<int>(momentum_L.size()); mi++)
        {
            TBOX_ASSERT(momentum_L[mi]);
        }
        for (int mi = 0; mi < static_cast<int>(momentum_R.size()); mi++)
        {
            TBOX_ASSERT(momentum_R[mi]);
        }
        TBOX_ASSERT(total_energy_L);
        TBOX_ASSERT(total_energy_R);
        for (int vi = 0; vi < static_cast<int>(volume_fraction_L.size()); vi++)
        {
            TBOX_ASSERT(volume_fraction_L[vi]);
        }
        for (int vi = 0; vi < static_cast<int>(volume_fraction_R.size()); vi++)
        {
            TBOX_ASSERT(volume_fraction_R[vi]);
        }
#endif
        
        switch (direction)
        {
            case X_DIRECTION:
            {
                std::vector<double*> F_x_intercell = flux_intercell;
                
                double& u_intercell = *(velocity_intercell[0]);
                double& v_intercell = *(velocity_intercell[1]);
                
                const std::vector<const double*> Z_rho_L = partial_density_L;
                const std::vector<const double*> Z_rho_R = partial_density_R;
                
                const double rho_L = d_equation_of_state->getTotalDensity(
                    partial_density_L);
                
                const double rho_R = d_equation_of_state->getTotalDensity(
                    partial_density_R);
                
                const double& rho_u_L = *(momentum_L[0]);
                const double& rho_u_R = *(momentum_R[0]);
                
                const double& rho_v_L = *(momentum_L[1]);
                const double& rho_v_R = *(momentum_R[1]);
                
                const std::vector<const double*> m_L = momentum_L;
                const std::vector<const double*> m_R = momentum_R;
                
                const double& u_L = rho_u_L/rho_L;
                const double& u_R = rho_u_R/rho_R;
                
                const double& v_L = rho_v_L/rho_L;
                const double& v_R = rho_v_R/rho_R;
                
                const double& E_L = *total_energy_L;
                const double& E_R = *total_energy_R;
                
                const std::vector<const double*> Z_L = volume_fraction_L;
                const std::vector<const double*> Z_R = volume_fraction_R;
                
                const double p_L = d_equation_of_state->getPressureWithVolumeFraction(
                    &rho_L,
                    m_L,
                    &E_L,
                    Z_L);
                
                const double p_R = d_equation_of_state->getPressureWithVolumeFraction(
                    &rho_R,
                    m_R,
                    &E_R,
                    Z_R);
                
                const double c_L = d_equation_of_state->getSoundSpeedWithVolumeFractionAndPressure(
                    &rho_L,
                    Z_L,
                    &p_L);
                
                const double c_R = d_equation_of_state->getSoundSpeedWithVolumeFractionAndPressure(
                    &rho_R,
                    Z_R,
                    &p_R);
                
                const double u_average = 0.5*(u_L + u_R);
                const double c_average = 0.5*(c_L + c_R);
                
                const double s_L = fmin(u_average - c_average, u_L - c_L);
                const double s_R = fmax(u_average + c_average, u_R + c_R);
                
                const double s_minus = fmin(0.0, s_L);
                const double s_plus  = fmax(0.0, s_R);
                
                const double s_star =
                    (p_R - p_L + rho_L*u_L*(s_L - u_L) - rho_R*u_R*(s_R - u_R))/(rho_L*(s_L - u_L) -
                        rho_R*(s_R - u_R));
                
                if (s_star > 0)
                {
                    const double Chi_star_L = (s_L - u_L)/(s_L - s_star);
                    
                    std::vector<double> Q_star_L;
                    for (int si = 0; si < d_num_species; si++)
                    {
                        Q_star_L.push_back(Chi_star_L*(*(Z_rho_L[si])));
                    }
                    Q_star_L.push_back(Chi_star_L*rho_L*s_star);
                    Q_star_L.push_back(Chi_star_L*rho_v_L);
                    Q_star_L.push_back(Chi_star_L*(E_L + (s_star - u_L)*(rho_L*s_star + p_L/(s_L - u_L))));
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q_star_L.push_back(Chi_star_L*(*(Z_L[si])));
                    }
                    
                    std::vector<const double*> Q_L;
                    for (int si = 0; si < d_num_species; si++)
                    {
                        Q_L.push_back(Z_rho_L[si]);
                    }
                    Q_L.push_back(m_L[0]);
                    Q_L.push_back(m_L[1]);
                    Q_L.push_back(&E_L);
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q_L.push_back(Z_L[si]);
                    }
                    
                    std::vector<double> F_x_L;
                    for (int si = 0; si < d_num_species; si++)
                    {
                        F_x_L.push_back(u_L*(*(Z_rho_L[si])));
                    }
                    F_x_L.push_back(rho_u_L*u_L + p_L);
                    F_x_L.push_back(rho_u_L*v_L);
                    F_x_L.push_back(u_L*(E_L + p_L));
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        F_x_L.push_back(u_L*(*(Z_L[si])));
                    }
                    
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        *(F_x_intercell[ei]) = F_x_L[ei] + s_minus*(Q_star_L[ei] - *(Q_L[ei]));
                    }
                    
                    u_intercell = u_L + s_minus*((s_L - u_L)/(s_L - s_star) - 1);
                    v_intercell = v_L;
                }
                else
                {
                    const double Chi_star_R = (s_R - u_R)/(s_R - s_star);
                    
                    std::vector<double> Q_star_R;
                    for (int si = 0; si < d_num_species; si++)
                    {
                        Q_star_R.push_back(Chi_star_R*(*(Z_rho_R[si])));
                    }
                    Q_star_R.push_back(Chi_star_R*rho_R*s_star);
                    Q_star_R.push_back(Chi_star_R*rho_v_R);
                    Q_star_R.push_back(Chi_star_R*(E_R + (s_star - u_R)*(rho_R*s_star + p_R/(s_R - u_R))));
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q_star_R.push_back(Chi_star_R*(*(Z_R[si])));
                    }
                    
                    std::vector<const double*> Q_R;
                    for (int si = 0; si < d_num_species; si++)
                    {
                        Q_R.push_back(Z_rho_R[si]);
                    }
                    Q_R.push_back(m_R[0]);
                    Q_R.push_back(m_R[1]);
                    Q_R.push_back(&E_R);
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q_R.push_back(Z_R[si]);
                    }
                    
                    std::vector<double> F_x_R;
                    for (int si = 0; si < d_num_species; si++)
                    {
                        F_x_R.push_back(u_R*(*(Z_rho_R[si])));
                    }
                    F_x_R.push_back(rho_u_R*u_R + p_R);
                    F_x_R.push_back(rho_u_R*v_R);
                    F_x_R.push_back(u_R*(E_R + p_R));
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        F_x_R.push_back(u_R*(*(Z_R[si])));
                    }
                    
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        *(F_x_intercell[ei]) = F_x_R[ei] + s_plus*(Q_star_R[ei] - *(Q_R[ei]));
                    }
                    
                    u_intercell = u_R + s_plus*((s_R - u_R)/(s_R - s_star) - 1);
                    v_intercell = v_R;
                }
                
                break;
            }
            case Y_DIRECTION:
            {
                std::vector<double*> F_y_intercell = flux_intercell;
                
                double& u_intercell = *(velocity_intercell[0]);
                double& v_intercell = *(velocity_intercell[1]);
                
                const std::vector<const double*> Z_rho_B = partial_density_L;
                const std::vector<const double*> Z_rho_T = partial_density_R;
                
                const double rho_B = d_equation_of_state->getTotalDensity(
                    partial_density_L);
                
                const double rho_T = d_equation_of_state->getTotalDensity(
                    partial_density_R);
                
                const double& rho_u_B = *(momentum_L[0]);
                const double& rho_u_T = *(momentum_R[0]);
                
                const double& rho_v_B = *(momentum_L[1]);
                const double& rho_v_T = *(momentum_R[1]);
                
                const std::vector<const double*> m_B = momentum_L;
                const std::vector<const double*> m_T = momentum_R;
                
                const double& u_B = rho_u_B/rho_B;
                const double& u_T = rho_u_T/rho_T;
                
                const double& v_B = rho_v_B/rho_B;
                const double& v_T = rho_v_T/rho_T;
                
                const double& E_B = *total_energy_L;
                const double& E_T = *total_energy_R;
                
                const std::vector<const double*> Z_B = volume_fraction_L;
                const std::vector<const double*> Z_T = volume_fraction_R;
                
                const double p_B = d_equation_of_state->getPressureWithVolumeFraction(
                    &rho_B,
                    m_B,
                    &E_B,
                    Z_B);
                
                const double p_T = d_equation_of_state->getPressureWithVolumeFraction(
                    &rho_T,
                    m_T,
                    &E_T,
                    Z_T);
                
                const double c_B = d_equation_of_state->getSoundSpeedWithVolumeFractionAndPressure(
                    &rho_B,
                    Z_B,
                    &p_B);
                
                const double c_T = d_equation_of_state->getSoundSpeedWithVolumeFractionAndPressure(
                    &rho_T,
                    Z_T,
                    &p_T);
                
                const double v_average = 0.5*(v_B + v_T);
                const double c_average = 0.5*(c_B + c_T);
                
                const double s_B = fmin(v_average - c_average, v_B - c_B);
                const double s_T = fmax(v_average + c_average, v_T + c_T);
                
                const double s_minus = fmin(0.0, s_B);
                const double s_plus  = fmax(0.0, s_T);
                
                const double s_star =
                    (p_T - p_B + rho_B*v_B*(s_B - v_B) - rho_T*v_T*(s_T - v_T))/(rho_B*(s_B - v_B) -
                        rho_T*(s_T - v_T));
                
                if (s_star > 0)
                {
                    const double Chi_star_B = (s_B - v_B)/(s_B - s_star);    
                    
                    std::vector<double> Q_star_B;
                    for (int si = 0; si < d_num_species; si++)
                    {
                        Q_star_B.push_back(Chi_star_B*(*(Z_rho_B[si])));
                    }
                    Q_star_B.push_back(Chi_star_B*rho_u_B);
                    Q_star_B.push_back(Chi_star_B*rho_B*s_star);
                    Q_star_B.push_back(Chi_star_B*(E_B + (s_star - v_B)*(rho_B*s_star + p_B/(s_B - v_B))));
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q_star_B.push_back(Chi_star_B*(*(Z_B[si])));
                    }
                    
                    std::vector<const double*> Q_B;
                    for (int si = 0; si < d_num_species; si++)
                    {
                        Q_B.push_back(Z_rho_B[si]);
                    }
                    Q_B.push_back(m_B[0]);
                    Q_B.push_back(m_B[1]);
                    Q_B.push_back(&E_B);
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q_B.push_back(Z_B[si]);
                    }
                    
                    std::vector<double> F_y_B;
                    for (int si = 0; si < d_num_species; si++)
                    {
                        F_y_B.push_back(v_B*(*(Z_rho_B[si])));
                    }
                    F_y_B.push_back(rho_v_B*u_B);
                    F_y_B.push_back(rho_v_B*v_B + p_B);
                    F_y_B.push_back(v_B*(E_B + p_B));
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        F_y_B.push_back(v_B*(*(Z_B[si])));
                    }
                    
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        *(F_y_intercell[ei]) = F_y_B[ei] + s_minus*(Q_star_B[ei] - *(Q_B[ei]));
                    }
                    
                    u_intercell = u_B;
                    v_intercell = v_B + s_minus*((s_B - v_B)/(s_B - s_star) - 1);
                }
                else
                {
                    const double Chi_star_T = (s_T - v_T)/(s_T - s_star); 
                    
                    std::vector<double> Q_star_T;
                    for (int si = 0; si < d_num_species; si++)
                    {
                        Q_star_T.push_back(Chi_star_T*(*(Z_rho_T[si])));
                    }
                    Q_star_T.push_back(Chi_star_T*rho_u_T);
                    Q_star_T.push_back(Chi_star_T*rho_T*s_star);
                    Q_star_T.push_back(Chi_star_T*(E_T + (s_star - v_T)*(rho_T*s_star + p_T/(s_T - v_T))));
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q_star_T.push_back(Chi_star_T*(*(Z_T[si])));
                    }
                    
                    std::vector<const double*> Q_T;
                    for (int si = 0; si < d_num_species; si++)
                    {
                        Q_T.push_back(Z_rho_T[si]);
                    }
                    Q_T.push_back(m_T[0]);
                    Q_T.push_back(m_T[1]);
                    Q_T.push_back(&E_T);
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q_T.push_back(Z_T[si]);
                    }
                    
                    std::vector<double> F_y_T;
                    for (int si = 0; si < d_num_species; si++)
                    {
                        F_y_T.push_back(v_T*(*(Z_rho_T[si])));
                    }
                    F_y_T.push_back(rho_v_T*u_T);
                    F_y_T.push_back(rho_v_T*v_T + p_T);
                    F_y_T.push_back(v_T*(E_T + p_T));
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        F_y_T.push_back(v_T*(*(Z_T[si])));
                    }
                    
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        *(F_y_intercell[ei]) = F_y_T[ei] + s_plus*(Q_star_T[ei] - *(Q_T[ei]));
                    }
                    
                    u_intercell = u_T;
                    v_intercell = v_T + s_plus*((s_T - v_T)/(s_T - s_star) - 1);
                }
                
                break;
            }
            case Z_DIRECTION:
            {
                TBOX_ERROR(d_object_name
                           << ": "
                           << "There is no z direction for 1D problem."
                           << std::endl);
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                           << ": "
                           << "Unknown direction."
                           << std::endl);
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(static_cast<int>(flux_intercell.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(velocity_intercell.size()) == 3);
        TBOX_ASSERT(static_cast<int>(partial_density_L.size()) == d_num_species);
        TBOX_ASSERT(static_cast<int>(partial_density_R.size()) == d_num_species);
        TBOX_ASSERT(static_cast<int>(momentum_L.size()) == 3);
        TBOX_ASSERT(static_cast<int>(momentum_R.size()) == 3);
        TBOX_ASSERT(static_cast<int>(volume_fraction_L.size()) == d_num_species - 1);
        TBOX_ASSERT(static_cast<int>(volume_fraction_R.size()) == d_num_species - 1);
        
        for (int fi = 0; fi < static_cast<int>(flux_intercell.size()); fi++)
        {
            TBOX_ASSERT(flux_intercell[fi]);
        }
        for (int vi = 0; vi < static_cast<int>(velocity_intercell.size()); vi++)
        {
            TBOX_ASSERT(velocity_intercell[vi]);
        }
        for (int pi = 0; pi < static_cast<int>(partial_density_L.size()); pi++)
        {
            TBOX_ASSERT(partial_density_L[pi]);
        }
        for (int pi = 0; pi < static_cast<int>(partial_density_R.size()); pi++)
        {
            TBOX_ASSERT(partial_density_R[pi]);
        }
        for (int mi = 0; mi < static_cast<int>(momentum_L.size()); mi++)
        {
            TBOX_ASSERT(momentum_L[mi]);
        }
        for (int mi = 0; mi < static_cast<int>(momentum_R.size()); mi++)
        {
            TBOX_ASSERT(momentum_R[mi]);
        }
        TBOX_ASSERT(total_energy_L);
        TBOX_ASSERT(total_energy_R);
        for (int vi = 0; vi < static_cast<int>(volume_fraction_L.size()); vi++)
        {
            TBOX_ASSERT(volume_fraction_L[vi]);
        }
        for (int vi = 0; vi < static_cast<int>(volume_fraction_R.size()); vi++)
        {
            TBOX_ASSERT(volume_fraction_R[vi]);
        }
#endif
        
        switch (direction)
        {
            case X_DIRECTION:
            {
                std::vector<double*> F_x_intercell = flux_intercell;
                
                double& u_intercell = *(velocity_intercell[0]);
                double& v_intercell = *(velocity_intercell[1]);
                double& w_intercell = *(velocity_intercell[2]);
                
                const std::vector<const double*> Z_rho_L = partial_density_L;
                const std::vector<const double*> Z_rho_R = partial_density_R;
                
                const double rho_L = d_equation_of_state->getTotalDensity(
                    partial_density_L);
                
                const double rho_R = d_equation_of_state->getTotalDensity(
                    partial_density_R);
                
                const double& rho_u_L = *(momentum_L[0]);
                const double& rho_u_R = *(momentum_R[0]);
                
                const double& rho_v_L = *(momentum_L[1]);
                const double& rho_v_R = *(momentum_R[1]);
                
                const double& rho_w_L = *(momentum_L[2]);
                const double& rho_w_R = *(momentum_R[2]);
                
                const std::vector<const double*> m_L = momentum_L;
                const std::vector<const double*> m_R = momentum_R;
                
                const double& u_L = rho_u_L/rho_L;
                const double& u_R = rho_u_R/rho_R;
                
                const double& v_L = rho_v_L/rho_L;
                const double& v_R = rho_v_R/rho_R;
                
                const double& w_L = rho_w_L/rho_L;
                const double& w_R = rho_w_R/rho_R;
                
                const double& E_L = *total_energy_L;
                const double& E_R = *total_energy_R;
                
                const std::vector<const double*> Z_L = volume_fraction_L;
                const std::vector<const double*> Z_R = volume_fraction_R;
                
                const double p_L = d_equation_of_state->getPressureWithVolumeFraction(
                    &rho_L,
                    m_L,
                    &E_L,
                    Z_L);
                
                const double p_R = d_equation_of_state->getPressureWithVolumeFraction(
                    &rho_R,
                    m_R,
                    &E_R,
                    Z_R);
                
                const double c_L = d_equation_of_state->getSoundSpeedWithVolumeFractionAndPressure(
                    &rho_L,
                    Z_L,
                    &p_L);
                
                const double c_R = d_equation_of_state->getSoundSpeedWithVolumeFractionAndPressure(
                    &rho_R,
                    Z_R,
                    &p_R);
                
                const double u_average = 0.5*(u_L + u_R);
                const double c_average = 0.5*(c_L + c_R);
                
                const double s_L = fmin(u_average - c_average, u_L - c_L);
                const double s_R = fmax(u_average + c_average, u_R + c_R);
                
                const double s_minus = fmin(0.0, s_L);
                const double s_plus  = fmax(0.0, s_R);
                
                const double s_star =
                    (p_R - p_L + rho_L*u_L*(s_L - u_L) - rho_R*u_R*(s_R - u_R))/(rho_L*(s_L - u_L) -
                        rho_R*(s_R - u_R));
                
                if (s_star > 0)
                {
                    const double Chi_star_L = (s_L - u_L)/(s_L - s_star);
                    
                    std::vector<double> Q_star_L;
                    for (int si = 0; si < d_num_species; si++)
                    {
                        Q_star_L.push_back(Chi_star_L*(*(Z_rho_L[si])));
                    }
                    Q_star_L.push_back(Chi_star_L*rho_L*s_star);
                    Q_star_L.push_back(Chi_star_L*rho_v_L);
                    Q_star_L.push_back(Chi_star_L*rho_w_L);
                    Q_star_L.push_back(Chi_star_L*(E_L + (s_star - u_L)*(rho_L*s_star + p_L/(s_L - u_L))));
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q_star_L.push_back(Chi_star_L*(*(Z_L[si])));
                    }
                    
                    std::vector<const double*> Q_L;
                    for (int si = 0; si < d_num_species; si++)
                    {
                        Q_L.push_back(Z_rho_L[si]);
                    }
                    Q_L.push_back(m_L[0]);
                    Q_L.push_back(m_L[1]);
                    Q_L.push_back(m_L[2]);
                    Q_L.push_back(&E_L);
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q_L.push_back(Z_L[si]);
                    }
                    
                    std::vector<double> F_x_L;
                    for (int si = 0; si < d_num_species; si++)
                    {
                        F_x_L.push_back(u_L*(*(Z_rho_L[si])));
                    }
                    F_x_L.push_back(rho_u_L*u_L + p_L);
                    F_x_L.push_back(rho_u_L*v_L);
                    F_x_L.push_back(rho_u_L*w_L);
                    F_x_L.push_back(u_L*(E_L + p_L));
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        F_x_L.push_back(u_L*(*(Z_L[si])));
                    }
                    
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        *(F_x_intercell[ei]) = F_x_L[ei] + s_minus*(Q_star_L[ei] - *(Q_L[ei]));
                    }
                    
                    u_intercell = u_L + s_minus*((s_L - u_L)/(s_L - s_star) - 1);
                    v_intercell = v_L;
                    w_intercell = w_L;
                }
                else
                {
                    const double Chi_star_R = (s_R - u_R)/(s_R - s_star);
                    
                    std::vector<double> Q_star_R;
                    for (int si = 0; si < d_num_species; si++)
                    {
                        Q_star_R.push_back(Chi_star_R*(*(Z_rho_R[si])));
                    }
                    Q_star_R.push_back(Chi_star_R*rho_R*s_star);
                    Q_star_R.push_back(Chi_star_R*rho_v_R);
                    Q_star_R.push_back(Chi_star_R*rho_w_R);
                    Q_star_R.push_back(Chi_star_R*(E_R + (s_star - u_R)*(rho_R*s_star + p_R/(s_R - u_R))));
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q_star_R.push_back(Chi_star_R*(*(Z_R[si])));
                    }
                    
                    std::vector<const double*> Q_R;
                    for (int si = 0; si < d_num_species; si++)
                    {
                        Q_R.push_back(Z_rho_R[si]);
                    }
                    Q_R.push_back(m_R[0]);
                    Q_R.push_back(m_R[1]);
                    Q_R.push_back(m_R[2]);
                    Q_R.push_back(&E_R);
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q_R.push_back(Z_R[si]);
                    }
                    
                    std::vector<double> F_x_R;
                    for (int si = 0; si < d_num_species; si++)
                    {
                        F_x_R.push_back(u_R*(*(Z_rho_R[si])));
                    }
                    F_x_R.push_back(rho_u_R*u_R + p_R);
                    F_x_R.push_back(rho_u_R*v_R);
                    F_x_R.push_back(rho_u_R*w_R);
                    F_x_R.push_back(u_R*(E_R + p_R));
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        F_x_R.push_back(u_R*(*(Z_R[si])));
                    }
                    
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        *(F_x_intercell[ei]) = F_x_R[ei] + s_plus*(Q_star_R[ei] - *(Q_R[ei]));
                    }
                    
                    u_intercell = u_R + s_plus*((s_R - u_R)/(s_R - s_star) - 1);
                    v_intercell = v_R;
                    w_intercell = w_R;
                }
                
                break;
            }
            case Y_DIRECTION:
            {
                std::vector<double*> F_y_intercell = flux_intercell;
                
                double& u_intercell = *(velocity_intercell[0]);
                double& v_intercell = *(velocity_intercell[1]);
                double& w_intercell = *(velocity_intercell[2]);
                
                const std::vector<const double*> Z_rho_B = partial_density_L;
                const std::vector<const double*> Z_rho_T = partial_density_R;
                
                const double rho_B = d_equation_of_state->getTotalDensity(
                    partial_density_L);
                
                const double rho_T = d_equation_of_state->getTotalDensity(
                    partial_density_R);
                
                const double& rho_u_B = *(momentum_L[0]);
                const double& rho_u_T = *(momentum_R[0]);
                
                const double& rho_v_B = *(momentum_L[1]);
                const double& rho_v_T = *(momentum_R[1]);
                
                const double& rho_w_B = *(momentum_L[2]);
                const double& rho_w_T = *(momentum_R[2]);
                
                const std::vector<const double*> m_B = momentum_L;
                const std::vector<const double*> m_T = momentum_R;
                
                const double& u_B = rho_u_B/rho_B;
                const double& u_T = rho_u_T/rho_T;
                
                const double& v_B = rho_v_B/rho_B;
                const double& v_T = rho_v_T/rho_T;
                
                const double& w_B = rho_w_B/rho_B;
                const double& w_T = rho_w_T/rho_T;
                
                const double& E_B = *total_energy_L;
                const double& E_T = *total_energy_R;
                
                const std::vector<const double*> Z_B = volume_fraction_L;
                const std::vector<const double*> Z_T = volume_fraction_R;
                
                const double p_B = d_equation_of_state->getPressureWithVolumeFraction(
                    &rho_B,
                    m_B,
                    &E_B,
                    Z_B);
                
                const double p_T = d_equation_of_state->getPressureWithVolumeFraction(
                    &rho_T,
                    m_T,
                    &E_T,
                    Z_T);
                
                const double c_B = d_equation_of_state->getSoundSpeedWithVolumeFractionAndPressure(
                    &rho_B,
                    Z_B,
                    &p_B);
                
                const double c_T = d_equation_of_state->getSoundSpeedWithVolumeFractionAndPressure(
                    &rho_T,
                    Z_T,
                    &p_T);
                
                const double v_average = 0.5*(v_B + v_T);
                const double c_average = 0.5*(c_B + c_T);
                
                const double s_B = fmin(v_average - c_average, v_B - c_B);
                const double s_T = fmax(v_average + c_average, v_T + c_T);
                
                const double s_minus = fmin(0.0, s_B);
                const double s_plus  = fmax(0.0, s_T);
                
                const double s_star =
                    (p_T - p_B + rho_B*v_B*(s_B - v_B) - rho_T*v_T*(s_T - v_T))/(rho_B*(s_B - v_B) -
                        rho_T*(s_T - v_T));
                
                if (s_star > 0)
                {
                    const double Chi_star_B = (s_B - v_B)/(s_B - s_star);    
                    
                    std::vector<double> Q_star_B;
                    for (int si = 0; si < d_num_species; si++)
                    {
                        Q_star_B.push_back(Chi_star_B*(*(Z_rho_B[si])));
                    }
                    Q_star_B.push_back(Chi_star_B*rho_u_B);
                    Q_star_B.push_back(Chi_star_B*rho_B*s_star);
                    Q_star_B.push_back(Chi_star_B*rho_w_B);
                    Q_star_B.push_back(Chi_star_B*(E_B + (s_star - v_B)*(rho_B*s_star + p_B/(s_B - v_B))));
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q_star_B.push_back(Chi_star_B*(*(Z_B[si])));
                    }
                    
                    std::vector<const double*> Q_B;
                    for (int si = 0; si < d_num_species; si++)
                    {
                        Q_B.push_back(Z_rho_B[si]);
                    }
                    Q_B.push_back(m_B[0]);
                    Q_B.push_back(m_B[1]);
                    Q_B.push_back(m_B[2]);
                    Q_B.push_back(&E_B);
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q_B.push_back(Z_B[si]);
                    }
                    
                    std::vector<double> F_y_B;
                    for (int si = 0; si < d_num_species; si++)
                    {
                        F_y_B.push_back(v_B*(*(Z_rho_B[si])));
                    }
                    F_y_B.push_back(rho_v_B*u_B);
                    F_y_B.push_back(rho_v_B*v_B + p_B);
                    F_y_B.push_back(rho_v_B*w_B);
                    F_y_B.push_back(v_B*(E_B + p_B));
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        F_y_B.push_back(v_B*(*(Z_B[si])));
                    }
                    
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        *(F_y_intercell[ei]) = F_y_B[ei] + s_minus*(Q_star_B[ei] - *(Q_B[ei]));
                    }
                    
                    u_intercell = u_B;
                    v_intercell = v_B + s_minus*((s_B - v_B)/(s_B - s_star) - 1);
                    w_intercell = w_B;
                }
                else
                {
                    const double Chi_star_T = (s_T - v_T)/(s_T - s_star); 
                    
                    std::vector<double> Q_star_T;
                    for (int si = 0; si < d_num_species; si++)
                    {
                        Q_star_T.push_back(Chi_star_T*(*(Z_rho_T[si])));
                    }
                    Q_star_T.push_back(Chi_star_T*rho_u_T);
                    Q_star_T.push_back(Chi_star_T*rho_T*s_star);
                    Q_star_T.push_back(Chi_star_T*rho_w_T);
                    Q_star_T.push_back(Chi_star_T*(E_T + (s_star - v_T)*(rho_T*s_star + p_T/(s_T - v_T))));
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q_star_T.push_back(Chi_star_T*(*(Z_T[si])));
                    }
                    
                    std::vector<const double*> Q_T;
                    for (int si = 0; si < d_num_species; si++)
                    {
                        Q_T.push_back(Z_rho_T[si]);
                    }
                    Q_T.push_back(m_T[0]);
                    Q_T.push_back(m_T[1]);
                    Q_T.push_back(m_T[2]);
                    Q_T.push_back(&E_T);
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q_T.push_back(Z_T[si]);
                    }
                    
                    std::vector<double> F_y_T;
                    for (int si = 0; si < d_num_species; si++)
                    {
                        F_y_T.push_back(v_T*(*(Z_rho_T[si])));
                    }
                    F_y_T.push_back(rho_v_T*u_T);
                    F_y_T.push_back(rho_v_T*v_T + p_T);
                    F_y_T.push_back(rho_v_T*w_T);
                    F_y_T.push_back(v_T*(E_T + p_T));
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        F_y_T.push_back(v_T*(*(Z_T[si])));
                    }
                    
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        *(F_y_intercell[ei]) = F_y_T[ei] + s_plus*(Q_star_T[ei] - *(Q_T[ei]));
                    }
                    
                    u_intercell = u_T;
                    v_intercell = v_T + s_plus*((s_T - v_T)/(s_T - s_star) - 1);
                    w_intercell = w_T;
                }
                
                break;
            }
            case Z_DIRECTION:
            {
                std::vector<double*> F_z_intercell = flux_intercell;
                
                double& u_intercell = *(velocity_intercell[0]);
                double& v_intercell = *(velocity_intercell[1]);
                double& w_intercell = *(velocity_intercell[2]);
                
                const std::vector<const double*> Z_rho_B = partial_density_L;
                const std::vector<const double*> Z_rho_F = partial_density_R;
                
                const double rho_B = d_equation_of_state->getTotalDensity(
                    partial_density_L);
                
                const double rho_F = d_equation_of_state->getTotalDensity(
                    partial_density_R);
                
                const double& rho_u_B = *(momentum_L[0]);
                const double& rho_u_F = *(momentum_R[0]);
                
                const double& rho_v_B = *(momentum_L[1]);
                const double& rho_v_F = *(momentum_R[1]);
                
                const double& rho_w_B = *(momentum_L[2]);
                const double& rho_w_F = *(momentum_R[2]);
                
                const std::vector<const double*> m_B = momentum_L;
                const std::vector<const double*> m_F = momentum_R;
                
                const double& u_B = rho_u_B/rho_B;
                const double& u_F = rho_u_F/rho_F;
                
                const double& v_B = rho_v_B/rho_B;
                const double& v_F = rho_v_F/rho_F;
                
                const double& w_B = rho_w_B/rho_B;
                const double& w_F = rho_w_F/rho_F;
                
                const double& E_B = *total_energy_L;
                const double& E_F = *total_energy_R;
                
                const std::vector<const double*> Z_B = volume_fraction_L;
                const std::vector<const double*> Z_F = volume_fraction_R;
                
                const double p_B = d_equation_of_state->getPressureWithVolumeFraction(
                    &rho_B,
                    m_B,
                    &E_B,
                    Z_B);
                
                const double p_F = d_equation_of_state->getPressureWithVolumeFraction(
                    &rho_F,
                    m_F,
                    &E_F,
                    Z_F);
                
                const double c_B = d_equation_of_state->getSoundSpeedWithVolumeFractionAndPressure(
                    &rho_B,
                    Z_B,
                    &p_B);
                
                const double c_F = d_equation_of_state->getSoundSpeedWithVolumeFractionAndPressure(
                    &rho_F,
                    Z_F,
                    &p_F);
                
                const double w_average = 0.5*(w_B + w_F);
                const double c_average = 0.5*(c_F + c_F);
                
                const double s_B = fmin(w_average - c_average, w_B - c_B);
                const double s_F = fmax(w_average + c_average, w_F + c_F);
                
                const double s_minus = fmin(0.0, s_B);
                const double s_plus  = fmax(0.0, s_F);
                
                const double s_star =
                    (p_F - p_B + rho_B*w_B*(s_B - w_B) - rho_F*w_F*(s_F - w_F))/(rho_B*(s_B - w_B) -
                        rho_F*(s_F - w_F));
                
                if (s_star > 0)
                {
                    const double Chi_star_B = (s_B - w_B)/(s_B - s_star);    
                    
                    std::vector<double> Q_star_B;
                    for (int si = 0; si < d_num_species; si++)
                    {
                        Q_star_B.push_back(Chi_star_B*(*(Z_rho_B[si])));
                    }
                    Q_star_B.push_back(Chi_star_B*rho_u_B);
                    Q_star_B.push_back(Chi_star_B*rho_v_B);
                    Q_star_B.push_back(Chi_star_B*rho_B*s_star);
                    Q_star_B.push_back(Chi_star_B*(E_B + (s_star - w_B)*(rho_B*s_star + p_B/(s_B - w_B))));
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q_star_B.push_back(Chi_star_B*(*(Z_B[si])));
                    }
                    
                    std::vector<const double*> Q_B;
                    for (int si = 0; si < d_num_species; si++)
                    {
                        Q_B.push_back(Z_rho_B[si]);
                    }
                    Q_B.push_back(m_B[0]);
                    Q_B.push_back(m_B[1]);
                    Q_B.push_back(m_B[2]);
                    Q_B.push_back(&E_B);
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q_B.push_back(Z_B[si]);
                    }
                    
                    std::vector<double> F_z_B;
                    for (int si = 0; si < d_num_species; si++)
                    {
                        F_z_B.push_back(w_B*(*(Z_rho_B[si])));
                    }
                    F_z_B.push_back(rho_w_B*u_B);
                    F_z_B.push_back(rho_w_B*v_B);
                    F_z_B.push_back(rho_w_B*w_B + p_B);
                    F_z_B.push_back(w_B*(E_B + p_B));
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        F_z_B.push_back(w_B*(*(Z_B[si])));
                    }
                    
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        *(F_z_intercell[ei]) = F_z_B[ei] + s_minus*(Q_star_B[ei] - *(Q_B[ei]));
                    }
                    
                    u_intercell = u_B;
                    v_intercell = v_B;
                    w_intercell = w_B + s_minus*((s_B - w_B)/(s_B - s_star) - 1);
                }
                else
                {
                    const double Chi_star_F = (s_F - w_F)/(s_F - s_star); 
                    
                    std::vector<double> Q_star_F;
                    for (int si = 0; si < d_num_species; si++)
                    {
                        Q_star_F.push_back(Chi_star_F*(*(Z_rho_F[si])));
                    }
                    Q_star_F.push_back(Chi_star_F*rho_u_F);
                    Q_star_F.push_back(Chi_star_F*rho_v_F);
                    Q_star_F.push_back(Chi_star_F*rho_F*s_star);
                    Q_star_F.push_back(Chi_star_F*(E_F + (s_star - w_F)*(rho_F*s_star + p_F/(s_F - w_F))));
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q_star_F.push_back(Chi_star_F*(*(Z_F[si])));
                    }
                    
                    std::vector<const double*> Q_F;
                    for (int si = 0; si < d_num_species; si++)
                    {
                        Q_F.push_back(Z_rho_F[si]);
                    }
                    Q_F.push_back(m_F[0]);
                    Q_F.push_back(m_F[1]);
                    Q_F.push_back(m_F[2]);
                    Q_F.push_back(&E_F);
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Q_F.push_back(Z_F[si]);
                    }
                    
                    std::vector<double> F_z_F;
                    for (int si = 0; si < d_num_species; si++)
                    {
                        F_z_F.push_back(w_F*(*(Z_rho_F[si])));
                    }
                    F_z_F.push_back(rho_w_F*u_F);
                    F_z_F.push_back(rho_w_F*v_F);
                    F_z_F.push_back(rho_w_F*w_F + p_F);
                    F_z_F.push_back(w_F*(E_F + p_F));
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        F_z_F.push_back(w_F*(*(Z_F[si])));
                    }
                    
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        *(F_z_intercell[ei]) = F_z_F[ei] + s_plus*(Q_star_F[ei] - *(Q_F[ei]));
                    }
                    
                    u_intercell = u_F;
                    v_intercell = v_F;
                    w_intercell = w_F + s_plus*((s_F - w_F)/(s_F - s_star) - 1);
                }
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                           << ": "
                           << "Unknown direction."
                           << std::endl);
            }
        }
    }
}


/*
 * Compute the flux in the x-direction at the intercell face for single-species flow model
 * from conservative variables.
 */
void
RiemannSolverHLLC::computeIntercellFluxForSingleSpeciesInXDirectionFromConservativeVariables(
    std::vector<double*>& F_x_intercell,
    std::vector<double>& Q_L,
    std::vector<double>& Q_R)
{
    if (d_dim == tbox::Dimension(1))
    {
        const double u_L = Q_L[1]/Q_L[0];
        const double u_R = Q_R[1]/Q_R[0];
        
        std::vector<const double*> m_L(1);
        std::vector<const double*> m_R(1);
        m_L[0] = &Q_L[1];
        m_R[0] = &Q_R[1];
        
        const double p_L = d_equation_of_state->getPressure(
            &Q_L[0],
            m_L,
            &Q_L[2]);
        
        const double p_R = d_equation_of_state->getPressure(
            &Q_R[0],
            m_R,
            &Q_R[2]);
        
        const double c_L = d_equation_of_state->getSoundSpeedWithPressure(
            &Q_L[0],
            &p_L);
        
        const double c_R = d_equation_of_state->getSoundSpeedWithPressure(
            &Q_R[0],
            &p_R);
        
        const double u_average = 0.5*(u_L + u_R);
        const double c_average = 0.5*(c_L + c_R);
        
        const double s_L = fmin(u_average - c_average, u_L - c_L);
        const double s_R = fmax(u_average + c_average, u_R + c_R);
        
        const double s_minus = fmin(0.0, s_L);
        const double s_plus  = fmax(0.0, s_R);
        
        const double s_star =
            (p_R - p_L + Q_L[1]*(s_L - u_L) - Q_R[1]*(s_R - u_R))/(Q_L[0]*(s_L - u_L) -
                Q_R[0]*(s_R - u_R));
        
        if (s_star > 0)
        {
            const double Chi_star_L = (s_L - u_L)/(s_L - s_star);
            
            double Q_star_L[3];
            Q_star_L[0] = Chi_star_L*Q_L[0];
            Q_star_L[1] = Chi_star_L*Q_L[0]*s_star;
            Q_star_L[2] = Chi_star_L*(Q_L[2] + (s_star - u_L)*(Q_L[0]*s_star + p_L/(s_L - u_L)));
            
            double F_x_L[3];
            F_x_L[0] = Q_L[1];
            F_x_L[1] = Q_L[1]*u_L + p_L;
            F_x_L[2] = u_L*(Q_L[2] + p_L);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                (*F_x_intercell[ei]) = F_x_L[ei] + s_minus*(Q_star_L[ei] - Q_L[ei]);
            }
        }
        else
        {
            const double Chi_star_R = (s_R - u_R)/(s_R - s_star);
            
            double Q_star_R[3];
            Q_star_R[0] = Chi_star_R*Q_R[0];
            Q_star_R[1] = Chi_star_R*Q_R[0]*s_star;
            Q_star_R[2] = Chi_star_R*(Q_R[2] + (s_star - u_R)*(Q_R[0]*s_star + p_R/(s_R - u_R)));
            
            double F_x_R[3];
            F_x_R[0] = Q_R[1];
            F_x_R[1] = Q_R[1]*u_R + p_R;
            F_x_R[2] = u_R*(Q_R[2] + p_R);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                (*F_x_intercell[ei]) = F_x_R[ei] + s_plus*(Q_star_R[ei] - Q_R[ei]);
            }
        }
        
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const double u_L = Q_L[1]/Q_L[0];
        const double u_R = Q_R[1]/Q_R[0];
        
        const double v_L = Q_L[2]/Q_L[0];
        const double v_R = Q_R[2]/Q_R[0];
        
        std::vector<const double*> m_L(2);
        std::vector<const double*> m_R(2);
        m_L[0] = &Q_L[1];
        m_R[0] = &Q_R[1];
        m_L[1] = &Q_L[2];
        m_R[1] = &Q_R[2];
        
        const double p_L = d_equation_of_state->getPressure(
            &Q_L[0],
            m_L,
            &Q_L[3]);
        
        const double p_R = d_equation_of_state->getPressure(
            &Q_R[0],
            m_R,
            &Q_R[3]);
        
        const double c_L = d_equation_of_state->getSoundSpeedWithPressure(
            &Q_L[0],
            &p_L);
        
        const double c_R = d_equation_of_state->getSoundSpeedWithPressure(
            &Q_R[0],
            &p_R);
        
        const double u_average = 0.5*(u_L + u_R);
        const double c_average = 0.5*(c_L + c_R);
        
        const double s_L = fmin(u_average - c_average, u_L - c_L);
        const double s_R = fmax(u_average + c_average, u_R + c_R);
        
        const double s_minus = fmin(0.0, s_L);
        const double s_plus  = fmax(0.0, s_R);
        
        const double s_star =
            (p_R - p_L + Q_L[1]*(s_L - u_L) - Q_R[1]*(s_R - u_R))/(Q_L[0]*(s_L - u_L) -
                Q_R[0]*(s_R - u_R));
        
        if (s_star > 0)
        {
            const double Chi_star_L = (s_L - u_L)/(s_L - s_star);
            
            double Q_star_L[4];
            Q_star_L[0] = Chi_star_L*Q_L[0];
            Q_star_L[1] = Chi_star_L*Q_L[0]*s_star;
            Q_star_L[2] = Chi_star_L*Q_L[2];
            Q_star_L[3] = Chi_star_L*(Q_L[3] + (s_star - u_L)*(Q_L[0]*s_star + p_L/(s_L - u_L)));
            
            double F_x_L[4];
            F_x_L[0] = Q_L[1];
            F_x_L[1] = Q_L[1]*u_L + p_L;
            F_x_L[2] = Q_L[1]*v_L;
            F_x_L[3] = u_L*(Q_L[3] + p_L);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                (*F_x_intercell[ei]) = F_x_L[ei] + s_minus*(Q_star_L[ei] - Q_L[ei]);
            }
        }
        else
        {
            const double Chi_star_R = (s_R - u_R)/(s_R - s_star);
            
            double Q_star_R[4];
            Q_star_R[0] = Chi_star_R*Q_R[0];
            Q_star_R[1] = Chi_star_R*Q_R[0]*s_star;
            Q_star_R[2] = Chi_star_R*Q_R[2];
            Q_star_R[3] = Chi_star_R*(Q_R[3] + (s_star - u_R)*(Q_R[0]*s_star + p_R/(s_R - u_R)));
            
            double F_x_R[4];
            F_x_R[0] = Q_R[1];
            F_x_R[1] = Q_R[1]*u_R + p_R;
            F_x_R[2] = Q_R[1]*v_R;
            F_x_R[3] = u_R*(Q_R[3] + p_R);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                (*F_x_intercell[ei]) = F_x_R[ei] + s_plus*(Q_star_R[ei] - Q_R[ei]);
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const double u_L = Q_L[1]/Q_L[0];
        const double u_R = Q_R[1]/Q_R[0];
        
        const double v_L = Q_L[2]/Q_L[0];
        const double v_R = Q_R[2]/Q_R[0];
        
        const double w_L = Q_L[3]/Q_L[0];
        const double w_R = Q_R[3]/Q_R[0];
        
        std::vector<const double*> m_L(3);
        std::vector<const double*> m_R(3);
        m_L[0] = &Q_L[1];
        m_R[0] = &Q_R[1];
        m_L[1] = &Q_L[2];
        m_R[1] = &Q_R[2];
        m_L[2] = &Q_L[3];
        m_R[2] = &Q_R[3];
        
        const double p_L = d_equation_of_state->getPressure(
            &Q_L[0],
            m_L,
            &Q_L[4]);
        
        const double p_R = d_equation_of_state->getPressure(
            &Q_R[0],
            m_R,
            &Q_R[4]);
        
        const double c_L = d_equation_of_state->getSoundSpeedWithPressure(
            &Q_L[0],
            &p_L);
        
        const double c_R = d_equation_of_state->getSoundSpeedWithPressure(
            &Q_R[0],
            &p_R);
        
        const double u_average = 0.5*(u_L + u_R);
        const double c_average = 0.5*(c_L + c_R);
        
        const double s_L = fmin(u_average - c_average, u_L - c_L);
        const double s_R = fmax(u_average + c_average, u_R + c_R);
        
        const double s_minus = fmin(0.0, s_L);
        const double s_plus  = fmax(0.0, s_R);
        
        const double s_star =
            (p_R - p_L + Q_L[1]*(s_L - u_L) - Q_R[1]*(s_R - u_R))/(Q_L[0]*(s_L - u_L) -
                Q_R[0]*(s_R - u_R));
        
        if (s_star > 0)
        {
            const double Chi_star_L = (s_L - u_L)/(s_L - s_star);
            
            double Q_star_L[5];
            Q_star_L[0] = Chi_star_L*Q_L[0];
            Q_star_L[1] = Chi_star_L*Q_L[0]*s_star;
            Q_star_L[2] = Chi_star_L*Q_L[2];
            Q_star_L[3] = Chi_star_L*Q_L[3];
            Q_star_L[4] = Chi_star_L*(Q_L[4] + (s_star - u_L)*(Q_L[0]*s_star + p_L/(s_L - u_L)));
            
            double F_x_L[5];
            F_x_L[0] = Q_L[1];
            F_x_L[1] = Q_L[1]*u_L + p_L;
            F_x_L[2] = Q_L[1]*v_L;
            F_x_L[3] = Q_L[1]*w_L;
            F_x_L[4] = (u_L*(Q_L[4] + p_L));
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                (*F_x_intercell[ei]) = F_x_L[ei] + s_minus*(Q_star_L[ei] - Q_L[ei]);
            }
        }
        else
        {
            const double Chi_star_R = (s_R - u_R)/(s_R - s_star);
            
            double Q_star_R[5];
            Q_star_R[0] = Chi_star_R*Q_R[0];
            Q_star_R[1] = Chi_star_R*Q_R[0]*s_star;
            Q_star_R[2] = Chi_star_R*Q_R[2];
            Q_star_R[3] = Chi_star_R*Q_R[3];
            Q_star_R[4] = Chi_star_R*(Q_R[4] + (s_star - u_R)*(Q_R[0]*s_star + p_R/(s_R - u_R)));
            
            double F_x_R[5];
            F_x_R[0] = Q_R[1];
            F_x_R[1] = Q_R[1]*u_R + p_R;
            F_x_R[2] = Q_R[1]*v_R;
            F_x_R[3] = Q_R[1]*w_R;
            F_x_R[4] = (u_R*(Q_R[4] + p_R));
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                (*F_x_intercell[ei]) = F_x_R[ei] + s_plus*(Q_star_R[ei] - Q_R[ei]);
            }
        }
    }
}


/*
 * Compute the flux in the y-direction at the intercell face for single-species flow model
 * from conservative variables.
 */
void
RiemannSolverHLLC::computeIntercellFluxForSingleSpeciesInYDirectionFromConservativeVariables(
    std::vector<double*>& F_y_intercell,
    std::vector<double>& Q_B,
    std::vector<double>& Q_T)
{
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no y direction for 1D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const double u_B = Q_B[1]/Q_B[0];
        const double u_T = Q_T[1]/Q_T[0];
        
        const double v_B = Q_B[2]/Q_B[0];
        const double v_T = Q_T[2]/Q_T[0];
        
        std::vector<const double*> m_B(2);
        std::vector<const double*> m_T(2);
        m_B[0] = &Q_B[1];
        m_T[0] = &Q_T[1];
        m_B[1] = &Q_B[2];
        m_T[1] = &Q_T[2];
        
        const double p_B = d_equation_of_state->getPressure(
            &Q_B[0],
            m_B,
            &Q_B[3]);
        
        const double p_T = d_equation_of_state->getPressure(
            &Q_T[0],
            m_T,
            &Q_T[3]);
        
        const double c_B = d_equation_of_state->getSoundSpeedWithPressure(
            &Q_B[0],
            &p_B);
        
        const double c_T = d_equation_of_state->getSoundSpeedWithPressure(
            &Q_T[0],
            &p_T);
        
        const double v_average = 0.5*(v_B + v_T);
        const double c_average = 0.5*(c_B + c_T);
        
        const double s_B = fmin(v_average - c_average, v_B - c_B);
        const double s_T = fmax(v_average + c_average, v_T + c_T);
        
        const double s_minus = fmin(0.0, s_B);
        const double s_plus  = fmax(0.0, s_T);
        
        const double s_star =
            (p_T - p_B + Q_B[2]*(s_B - v_B) - Q_T[2]*(s_T - v_T))/(Q_B[0]*(s_B - v_B) -
                Q_T[0]*(s_T - v_T));
        
        if (s_star > 0)
        {
            const double Chi_star_B = (s_B - v_B)/(s_B - s_star);
            
            double Q_star_B[4];
            Q_star_B[0] = Chi_star_B*Q_B[0];
            Q_star_B[1] = Chi_star_B*Q_B[1];
            Q_star_B[2] = Chi_star_B*Q_B[0]*s_star;
            Q_star_B[3] = Chi_star_B*(Q_B[3] + (s_star - v_B)*(Q_B[0]*s_star + p_B/(s_B - v_B)));
            
            double F_y_B[4];
            F_y_B[0] = Q_B[2];
            F_y_B[0] = Q_B[2]*u_B;
            F_y_B[0] = Q_B[2]*v_B + p_B;
            F_y_B[0] = v_B*(Q_B[3] + p_B);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                (*F_y_intercell[ei]) = F_y_B[ei] + s_minus*(Q_star_B[ei] - Q_B[ei]);
            }
        }
        else
        {
            const double Chi_star_T = (s_T - v_T)/(s_T - s_star);
            
            double Q_star_T[4];
            Q_star_T[0] = Chi_star_T*Q_T[0];
            Q_star_T[1] = Chi_star_T*Q_T[1];
            Q_star_T[2] = Chi_star_T*Q_T[0]*s_star;
            Q_star_T[3] = Chi_star_T*(Q_T[3] + (s_star - v_T)*(Q_T[0]*s_star + p_T/(s_T - v_T)));
            
            double F_y_T[4];
            F_y_T[0] = Q_T[2];
            F_y_T[0] = Q_T[2]*u_T;
            F_y_T[0] = Q_T[2]*v_T + p_T;
            F_y_T[0] = v_T*(Q_T[3] + p_T);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                (*F_y_intercell[ei]) = F_y_T[ei] + s_plus*(Q_star_T[ei] - Q_T[ei]);
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const double u_B = Q_B[1]/Q_B[0];
        const double u_T = Q_T[1]/Q_T[0];
        
        const double v_B = Q_B[2]/Q_B[0];
        const double v_T = Q_T[2]/Q_T[0];
        
        const double w_B = Q_B[3]/Q_B[0];
        const double w_T = Q_T[3]/Q_T[0];
        
        std::vector<const double*> m_B(3);
        std::vector<const double*> m_T(3);
        m_B[0] = &Q_B[1];
        m_T[0] = &Q_T[1];
        m_B[1] = &Q_B[2];
        m_T[1] = &Q_T[2];
        m_B[2] = &Q_B[3];
        m_T[2] = &Q_T[3];
        
        const double p_B = d_equation_of_state->getPressure(
            &Q_B[0],
            m_B,
            &Q_B[4]);
        
        const double p_T = d_equation_of_state->getPressure(
            &Q_T[0],
            m_T,
            &Q_T[0]);
        
        const double c_B = d_equation_of_state->getSoundSpeedWithPressure(
            &Q_B[0],
            &p_B);
        
        const double c_T = d_equation_of_state->getSoundSpeedWithPressure(
            &Q_T[0],
            &p_T);
        
        const double v_average = 0.5*(v_B + v_T);
        const double c_average = 0.5*(c_B + c_T);
        
        const double s_B = fmin(v_average - c_average, v_B - c_B);
        const double s_T = fmax(v_average + c_average, v_T + c_T);
        
        const double s_minus = fmin(0.0, s_B);
        const double s_plus  = fmax(0.0, s_T);
        
        const double s_star =
            (p_T - p_B + Q_B[2]*(s_B - v_B) - Q_T[2]*(s_T - v_T))/(Q_B[0]*(s_B - v_B) -
                Q_T[0]*(s_T - v_T));
        
        if (s_star > 0)
        {
            const double Chi_star_B = (s_B - v_B)/(s_B - s_star);
            
            double Q_star_B[5];
            Q_star_B[0] = Chi_star_B*Q_B[0];
            Q_star_B[1] = Chi_star_B*Q_B[1];
            Q_star_B[2] = Chi_star_B*Q_B[0]*s_star;
            Q_star_B[3] = Chi_star_B*Q_B[3];
            Q_star_B[4] = Chi_star_B*(Q_B[4] + (s_star - v_B)*(Q_B[0]*s_star + p_B/(s_B - v_B)));
            
            double F_y_B[5];
            F_y_B[0] = Q_B[2];
            F_y_B[1] = Q_B[2]*u_B;
            F_y_B[2] = Q_B[2]*v_B + p_B;
            F_y_B[3] = Q_B[2]*w_B;
            F_y_B[4] = v_B*(Q_B[4] + p_B);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                (*F_y_intercell[ei]) = F_y_B[ei] + s_minus*(Q_star_B[ei] - Q_B[ei]);
            }
        }
        else
        {
            const double Chi_star_T = (s_T - v_T)/(s_T - s_star);
            
            double Q_star_T[5];
            Q_star_T[0] = Chi_star_T*Q_T[0];
            Q_star_T[1] = Chi_star_T*Q_T[1];
            Q_star_T[2] = Chi_star_T*Q_T[0]*s_star;
            Q_star_T[3] = Chi_star_T*Q_T[3];
            Q_star_T[4] = Chi_star_T*(Q_T[4] + (s_star - v_T)*(Q_T[0]*s_star + p_T/(s_T - v_T)));
            
            double F_y_T[5];
            F_y_T[0] = Q_T[2];
            F_y_T[1] = Q_T[2]*u_T;
            F_y_T[2] = Q_T[2]*v_T + p_T;
            F_y_T[3] = Q_T[2]*w_T;
            F_y_T[4] = v_T*(Q_T[4] + p_T);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                (*F_y_intercell[ei]) = F_y_T[ei] + s_plus*(Q_star_T[ei] - Q_T[ei]);
            }
        }
    }
}


/*
 * Compute the flux in the z-direction at the intercell face for single-species flow model
 * from conservative variables.
 */
void
RiemannSolverHLLC::computeIntercellFluxForSingleSpeciesInZDirectionFromConservativeVariables(
    std::vector<double*>& F_z_intercell,
    std::vector<double>& Q_B,
    std::vector<double>& Q_F)
{
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no z direction for 1D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no z direction for 2D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const double u_B = Q_B[1]/Q_B[0];
        const double u_F = Q_F[1]/Q_F[0];
        
        const double v_B = Q_B[2]/Q_B[0];
        const double v_F = Q_F[2]/Q_F[0];
        
        const double w_B = Q_B[3]/Q_B[0];
        const double w_F = Q_F[3]/Q_F[0];
        
        std::vector<const double*> m_B(3);
        std::vector<const double*> m_F(3);
        m_B[0] = &Q_B[1];
        m_F[0] = &Q_F[1];
        m_B[1] = &Q_B[2];
        m_F[1] = &Q_F[2];
        m_B[2] = &Q_B[3];
        m_F[2] = &Q_F[3];
        
        const double p_B = d_equation_of_state->getPressure(
            &Q_B[0],
            m_B,
            &Q_B[4]);
        
        const double p_F = d_equation_of_state->getPressure(
            &Q_F[0],
            m_F,
            &Q_F[4]);
        
        const double c_B = d_equation_of_state->getSoundSpeedWithPressure(
            &Q_B[0],
            &p_B);
        
        const double c_F = d_equation_of_state->getSoundSpeedWithPressure(
            &Q_F[0],
            &p_F);
        
        const double w_average = 0.5*(w_B + w_F);
        const double c_average = 0.5*(c_B + c_F);
        
        const double s_B = fmin(w_average - c_average, w_B - c_B);
        const double s_F = fmax(w_average + c_average, w_F + c_F);
        
        const double s_minus = fmin(0.0, s_B);
        const double s_plus  = fmax(0.0, s_F);
        
        const double s_star =
            (p_F - p_B + Q_B[3]*(s_B - w_B) - Q_F[3]*(s_F - w_F))/(Q_B[0]*(s_B - w_B) -
                Q_F[0]*(s_F - w_F));
        
        if (s_star > 0)
        {
            const double Chi_star_B = (s_B - w_B)/(s_B - s_star);
            
            double Q_star_B[5];
            Q_star_B[0] = Chi_star_B*Q_B[0];
            Q_star_B[1] = Chi_star_B*Q_B[1];
            Q_star_B[2] = Chi_star_B*Q_B[2];
            Q_star_B[3] = Chi_star_B*Q_B[0]*s_star;
            Q_star_B[4] = Chi_star_B*(Q_B[4] + (s_star - w_B)*(Q_B[0]*s_star + p_B/(s_B - w_B)));
            
            double F_z_B[5];
            F_z_B[0] = Q_B[3];
            F_z_B[1] = Q_B[3]*u_B;
            F_z_B[2] = Q_B[3]*v_B;
            F_z_B[3] = Q_B[3]*w_B + p_B;
            F_z_B[4] = w_B*(Q_B[4] + p_B);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                (*F_z_intercell[ei]) = F_z_B[ei] + s_minus*(Q_star_B[ei] - Q_B[ei]);
            }
        }
        else
        {
            const double Chi_star_F = (s_F - w_F)/(s_F - s_star);
            
            double Q_star_F[5];
            Q_star_F[0] = Chi_star_F*Q_F[0];
            Q_star_F[1] = Chi_star_F*Q_F[1];
            Q_star_F[2] = Chi_star_F*Q_F[2];
            Q_star_F[3] = Chi_star_F*Q_F[0]*s_star;
            Q_star_F[4] = Chi_star_F*(Q_F[4] + (s_star - w_F)*(Q_F[0]*s_star + p_F/(s_F - w_F)));
            
            double F_z_F[5];
            F_z_F[0] = Q_F[3];
            F_z_F[1] = Q_F[3]*u_F;
            F_z_F[2] = Q_F[3]*v_F;
            F_z_F[3] = Q_F[3]*w_F + p_F;
            F_z_F[4] = w_F*(Q_F[4] + p_F);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                (*F_z_intercell[ei]) = F_z_F[ei] + s_plus*(Q_star_F[ei] - Q_F[ei]);
            }
        }
    }
}


/*
 * Compute the flux in the x-direction at the intercell face for single-species flow model
 * from primitive variables.
 */
void
RiemannSolverHLLC::computeIntercellFluxForSingleSpeciesInXDirectionFromPrimitiveVariables(
    std::vector<double*>& F_x_intercell,
    std::vector<double>& V_L,
    std::vector<double>& V_R)
{
    if (d_dim == tbox::Dimension(1))
    {
        const double c_L = d_equation_of_state->getSoundSpeedWithPressure(
            &V_L[0],
            &V_L[2]);
        
        const double c_R = d_equation_of_state->getSoundSpeedWithPressure(
            &V_R[0],
            &V_R[2]);
        
        const double u_average = 0.5*(V_L[1] + V_R[1]);
        const double c_average = 0.5*(c_L + c_R);
        
        const double s_L = fmin(u_average - c_average, V_L[1] - c_L);
        const double s_R = fmax(u_average + c_average, V_R[1] + c_R);
        
        const double s_minus = fmin(0.0, s_L);
        const double s_plus  = fmax(0.0, s_R);
        
        const double s_star =
            (V_R[2] - V_L[2] + V_L[0]*V_L[1]*(s_L - V_L[1]) - V_R[0]*V_R[1]*(s_R - V_R[1]))/(V_L[0]*(s_L - V_L[1]) -
                V_R[0]*(s_R - V_R[1]));
        
        if (s_star > 0)
        {
            const double Chi_star_L = (s_L - V_L[1])/(s_L - s_star);
            
            std::vector<const double*> vel_L(1);
            vel_L[0] = &V_L[1];
            
            std::vector<double> Q_L(3);
            Q_L[0] = V_L[0];
            Q_L[1] = V_L[0]*V_L[1];
            Q_L[2] = d_equation_of_state->getTotalEnergy(
                &V_L[0],
                vel_L,
                &V_L[2]);
            
            std::vector<double> Q_star_L(3);
            Q_star_L[0] = Chi_star_L*V_L[0];
            Q_star_L[1] = Chi_star_L*V_L[0]*s_star;
            Q_star_L[2] = Chi_star_L*(Q_L[2] + (s_star - V_L[1])*(V_L[0]*s_star + V_L[2]/(s_L - V_L[1])));
            
            std::vector<double> F_x_L(3);
            F_x_L[0] = Q_L[1];
            F_x_L[1] = Q_L[1]*V_L[1] + V_L[2];
            F_x_L[2] = V_L[1]*(Q_L[2] + V_L[2]);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                (*F_x_intercell[ei]) = F_x_L[ei] + s_minus*(Q_star_L[ei] - Q_L[ei]);
            }
        }
        else
        {
            const double Chi_star_R = (s_R - V_R[1])/(s_R - s_star);
            
            std::vector<const double*> vel_R(1);
            vel_R[0] = &V_R[1];
            
            std::vector<double> Q_R(3);
            Q_R[0] = V_R[0];
            Q_R[1] = V_R[0]*V_R[1];
            Q_R[2] = d_equation_of_state->getTotalEnergy(
                &V_R[0],
                vel_R,
                &V_R[2]);
            
            std::vector<double> Q_star_R(3);
            Q_star_R[0] = Chi_star_R*V_R[0];
            Q_star_R[1] = Chi_star_R*V_R[0]*s_star;
            Q_star_R[2] = Chi_star_R*(Q_R[2] + (s_star - V_R[1])*(V_R[0]*s_star + V_R[2]/(s_R - V_R[1])));
            
            std::vector<double> F_x_R(3);
            F_x_R[0] = Q_R[1];
            F_x_R[1] = Q_R[1]*V_R[1] + V_R[2];
            F_x_R[2] = V_R[1]*(Q_R[2] + V_R[2]);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                (*F_x_intercell[ei]) = F_x_R[ei] + s_plus*(Q_star_R[ei] - Q_R[ei]);
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const double c_L = d_equation_of_state->getSoundSpeedWithPressure(
            &V_L[0],
            &V_L[3]);
        
        const double c_R = d_equation_of_state->getSoundSpeedWithPressure(
            &V_R[0],
            &V_R[3]);
        
        const double u_average = 0.5*(V_L[1] + V_R[1]);
        const double c_average = 0.5*(c_L + c_R);
        
        const double s_L = fmin(u_average - c_average, V_L[1] - c_L);
        const double s_R = fmax(u_average + c_average, V_R[1] + c_R);
        
        const double s_minus = fmin(0.0, s_L);
        const double s_plus  = fmax(0.0, s_R);
        
        const double s_star =
            (V_R[3] - V_L[3] + V_L[0]*V_L[1]*(s_L - V_L[1]) - V_R[0]*V_R[1]*(s_R - V_R[1]))/(V_L[0]*(s_L - V_L[1]) -
                V_R[0]*(s_R - V_R[1]));
        
        if (s_star > 0)
        {
            const double Chi_star_L = (s_L - V_L[1])/(s_L - s_star);
            
            std::vector<const double*> vel_L(2);
            vel_L[0] = &V_L[1];
            vel_L[1] = &V_L[2];
            
            std::vector<double> Q_L(4);
            Q_L[0] = V_L[0];
            Q_L[1] = V_L[0]*V_L[1];
            Q_L[2] = V_L[0]*V_L[2];
            Q_L[3] = d_equation_of_state->getTotalEnergy(
                &V_L[0],
                vel_L,
                &V_L[3]);
            
            std::vector<double> Q_star_L(4);
            Q_star_L[0] = Chi_star_L*V_L[0];
            Q_star_L[1] = Chi_star_L*V_L[0]*s_star;
            Q_star_L[2] = Chi_star_L*Q_L[2];
            Q_star_L[3] = Chi_star_L*(Q_L[3] + (s_star - V_L[1])*(V_L[0]*s_star + V_L[3]/(s_L - V_L[1])));
            
            std::vector<double> F_x_L(4);
            F_x_L[0] = Q_L[1];
            F_x_L[1] = Q_L[1]*V_L[1] + V_L[3];
            F_x_L[2] = Q_L[1]*V_L[2];
            F_x_L[3] = V_L[1]*(Q_L[3] + V_L[3]);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                (*F_x_intercell[ei]) = F_x_L[ei] + s_minus*(Q_star_L[ei] - Q_L[ei]);
            }
        }
        else
        {
            const double Chi_star_R = (s_R - V_R[1])/(s_R - s_star);
            
            std::vector<const double*> vel_R(2);
            vel_R[0] = &V_R[1];
            vel_R[1] = &V_R[2];
            
            std::vector<double> Q_R(4);
            Q_R[0] = V_R[0];
            Q_R[1] = V_R[0]*V_R[1];
            Q_R[2] = V_R[0]*V_R[2];
            Q_R[3] = d_equation_of_state->getTotalEnergy(
                &V_R[0],
                vel_R,
                &V_R[3]);
            
            std::vector<double> Q_star_R(4);
            Q_star_R[0] = Chi_star_R*V_R[0];
            Q_star_R[1] = Chi_star_R*V_R[0]*s_star;
            Q_star_R[2] = Chi_star_R*Q_R[2];
            Q_star_R[3] = Chi_star_R*(Q_R[3] + (s_star - V_R[1])*(V_R[0]*s_star + V_R[3]/(s_R - V_R[1])));
            
            std::vector<double> F_x_R(4);
            F_x_R[0] = Q_R[1];
            F_x_R[1] = Q_R[1]*V_R[1] + V_R[3];
            F_x_R[2] = Q_R[1]*V_R[2];
            F_x_R[3] = V_R[1]*(Q_R[3] + V_R[3]);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                (*F_x_intercell[ei]) = F_x_R[ei] + s_plus*(Q_star_R[ei] - Q_R[ei]);
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const double c_L = d_equation_of_state->getSoundSpeedWithPressure(
            &V_L[0],
            &V_L[4]);
        
        const double c_R = d_equation_of_state->getSoundSpeedWithPressure(
            &V_R[0],
            &V_R[4]);
        
        const double u_average = 0.5*(V_L[1] + V_R[1]);
        const double c_average = 0.5*(c_L + c_R);
        
        const double s_L = fmin(u_average - c_average, V_L[1] - c_L);
        const double s_R = fmax(u_average + c_average, V_R[1] + c_R);
        
        const double s_minus = fmin(0.0, s_L);
        const double s_plus  = fmax(0.0, s_R);
        
        const double s_star =
            (V_R[4] - V_L[4] + V_L[0]*V_L[1]*(s_L - V_L[1]) - V_R[0]*V_R[1]*(s_R - V_R[1]))/(V_L[0]*(s_L - V_L[1]) -
                V_R[0]*(s_R - V_R[1]));
        
        if (s_star > 0)
        {
            const double Chi_star_L = (s_L - V_L[1])/(s_L - s_star);
            
            std::vector<const double*> vel_L(3);
            vel_L[0] = &V_L[1];
            vel_L[1] = &V_L[2];
            vel_L[2] = &V_L[3];
            
            std::vector<double> Q_L(5);
            Q_L[0] = V_L[0];
            Q_L[1] = V_L[0]*V_L[1];
            Q_L[2] = V_L[0]*V_L[2];
            Q_L[3] = V_L[0]*V_L[3];
            Q_L[4] = d_equation_of_state->getTotalEnergy(
                &V_L[0],
                vel_L,
                &V_L[4]);
            
            std::vector<double> Q_star_L(5);
            Q_star_L[0] = Chi_star_L*V_L[0];
            Q_star_L[1] = Chi_star_L*V_L[0]*s_star;
            Q_star_L[2] = Chi_star_L*Q_L[2];
            Q_star_L[3] = Chi_star_L*Q_L[3];
            Q_star_L[4] = Chi_star_L*(Q_L[4] + (s_star - V_L[1])*(V_L[0]*s_star + V_L[4]/(s_L - V_L[1])));
            
            std::vector<double> F_x_L(5);
            F_x_L[0] = Q_L[1];
            F_x_L[1] = Q_L[1]*V_L[1] + V_L[4];
            F_x_L[2] = Q_L[1]*V_L[2];
            F_x_L[3] = Q_L[1]*V_L[3];
            F_x_L[4] = V_L[1]*(Q_L[4] + V_L[4]);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                (*F_x_intercell[ei]) = F_x_L[ei] + s_minus*(Q_star_L[ei] - Q_L[ei]);
            }
        }
        else
        {
            const double Chi_star_R = (s_R - V_R[1])/(s_R - s_star);
            
            std::vector<const double*> vel_R(3);
            vel_R[0] = &V_R[1];
            vel_R[1] = &V_R[2];
            vel_R[2] = &V_R[3];
            
            std::vector<double> Q_R(5);
            Q_R[0] = V_R[0];
            Q_R[1] = V_R[0]*V_R[1];
            Q_R[2] = V_R[0]*V_R[2];
            Q_R[3] = V_R[0]*V_R[3];
            Q_R[4] = d_equation_of_state->getTotalEnergy(
                &V_R[0],
                vel_R,
                &V_R[4]);
            
            std::vector<double> Q_star_R(5);
            Q_star_R[0] = Chi_star_R*V_R[0];
            Q_star_R[1] = Chi_star_R*V_R[0]*s_star;
            Q_star_R[2] = Chi_star_R*Q_R[2];
            Q_star_R[3] = Chi_star_R*Q_R[3];
            Q_star_R[4] = Chi_star_R*(Q_R[4] + (s_star - V_R[1])*(V_R[0]*s_star + V_R[4]/(s_R - V_R[1])));
            
            std::vector<double> F_x_R(5);
            F_x_R[0] = Q_R[1];
            F_x_R[1] = Q_R[1]*V_R[1] + V_R[4];
            F_x_R[2] = Q_R[1]*V_R[2];
            F_x_R[3] = Q_R[1]*V_R[3];
            F_x_R[4] = V_R[1]*(Q_R[4] + V_R[4]);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                (*F_x_intercell[ei]) = F_x_R[ei] + s_plus*(Q_star_R[ei] - Q_R[ei]);
            }
        }
    }
}


/*
 * Compute the flux in the y-direction at the intercell face for single-species flow model
 * from primitive variables.
 */
void
RiemannSolverHLLC::computeIntercellFluxForSingleSpeciesInYDirectionFromPrimitiveVariables(
    std::vector<double*>& F_y_intercell,
    std::vector<double>& V_B,
    std::vector<double>& V_T)
{
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no y direction for 1D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const double c_B = d_equation_of_state->getSoundSpeedWithPressure(
            &V_B[0],
            &V_B[3]);
        
        const double c_T = d_equation_of_state->getSoundSpeedWithPressure(
            &V_T[0],
            &V_T[3]);
        
        const double v_average = 0.5*(V_B[2] + V_T[2]);
        const double c_average = 0.5*(c_B + c_T);
        
        const double s_B = fmin(v_average - c_average, V_B[2] - c_B);
        const double s_T = fmax(v_average + c_average, V_T[2] + c_T);
        
        const double s_minus = fmin(0.0, s_B);
        const double s_plus  = fmax(0.0, s_T);
        
        const double s_star =
            (V_T[3] - V_B[3] + V_B[0]*V_B[2]*(s_B - V_B[2]) - V_T[0]*V_T[2]*(s_T - V_T[2]))/(V_B[0]*(s_B - V_B[2]) -
                V_T[0]*(s_T - V_T[2]));
        
        if (s_star > 0)
        {
            const double Chi_star_B = (s_B - V_B[2])/(s_B - s_star);
            
            std::vector<const double*> vel_B(2);
            vel_B[0] = &V_B[1];
            vel_B[1] = &V_B[2];
            
            std::vector<double> Q_B(4);
            Q_B[0] = V_B[0];
            Q_B[1] = V_B[0]*V_B[1];
            Q_B[2] = V_B[0]*V_B[2];
            Q_B[3] = d_equation_of_state->getTotalEnergy(
                &V_B[0],
                vel_B,
                &V_B[3]);
            
            std::vector<double> Q_star_B(4);
            Q_star_B[0] = Chi_star_B*V_B[0];
            Q_star_B[1] = Chi_star_B*Q_B[1];
            Q_star_B[2] = Chi_star_B*V_B[0]*s_star;
            Q_star_B[3] = Chi_star_B*(Q_B[3] + (s_star - V_B[2])*(V_B[0]*s_star + V_B[3]/(s_B - V_B[2])));
            
            std::vector<double> F_y_B(4);
            F_y_B[0] = Q_B[2];
            F_y_B[1] = Q_B[2]*V_B[1];
            F_y_B[2] = Q_B[2]*V_B[2] + V_B[3];
            F_y_B[3] = V_B[2]*(Q_B[3] + V_B[3]);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                (*F_y_intercell[ei]) = F_y_B[ei] + s_minus*(Q_star_B[ei] - Q_B[ei]);
            }
        }
        else
        {
            const double Chi_star_T = (s_T - V_T[2])/(s_T - s_star);
            
            std::vector<const double*> vel_T(2);
            vel_T[0] = &V_T[1];
            vel_T[1] = &V_T[2];
            
            std::vector<double> Q_T(4);
            Q_T[0] = V_T[0];
            Q_T[1] = V_T[0]*V_T[1];
            Q_T[2] = V_T[0]*V_T[2];
            Q_T[3] = d_equation_of_state->getTotalEnergy(
                &V_T[0],
                vel_T,
                &V_T[3]);
            
            std::vector<double> Q_star_T(4);
            Q_star_T[0] = Chi_star_T*V_T[0];
            Q_star_T[1] = Chi_star_T*Q_T[1];
            Q_star_T[2] = Chi_star_T*V_T[0]*s_star;
            Q_star_T[3] = Chi_star_T*(Q_T[3] + (s_star - V_T[2])*(V_T[0]*s_star + V_T[3]/(s_T - V_T[2])));
            
            std::vector<double> F_y_T(4);
            F_y_T[0] = Q_T[2];
            F_y_T[1] = Q_T[2]*V_T[1];
            F_y_T[2] = Q_T[2]*V_T[2] + V_T[3];
            F_y_T[3] = V_T[2]*(Q_T[3] + V_T[3]);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                (*F_y_intercell[ei]) = F_y_T[ei] + s_plus*(Q_star_T[ei] - Q_T[ei]);
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const double c_B = d_equation_of_state->getSoundSpeedWithPressure(
            &V_B[0],
            &V_B[4]);
        
        const double c_T = d_equation_of_state->getSoundSpeedWithPressure(
            &V_T[0],
            &V_T[4]);
        
        const double v_average = 0.5*(V_B[2] + V_T[2]);
        const double c_average = 0.5*(c_B + c_T);
        
        const double s_B = fmin(v_average - c_average, V_B[2] - c_B);
        const double s_T = fmax(v_average + c_average, V_T[2] + c_T);
        
        const double s_minus = fmin(0.0, s_B);
        const double s_plus  = fmax(0.0, s_T);
        
        const double s_star =
            (V_T[4] - V_B[4] + V_B[0]*V_B[2]*(s_B - V_B[2]) - V_T[0]*V_T[2]*(s_T - V_T[2]))/(V_B[0]*(s_B - V_B[2]) -
                V_T[0]*(s_T - V_T[2]));
        
        if (s_star > 0)
        {
            const double Chi_star_B = (s_B - V_B[2])/(s_B - s_star);
            
            std::vector<const double*> vel_B(3);
            vel_B[0] = &V_B[1];
            vel_B[1] = &V_B[2];
            vel_B[2] = &V_B[3];
            
            std::vector<double> Q_B(5);
            Q_B[0] = V_B[0];
            Q_B[1] = V_B[0]*V_B[1];
            Q_B[2] = V_B[0]*V_B[2];
            Q_B[3] = V_B[0]*V_B[3];
            Q_B[4] = d_equation_of_state->getTotalEnergy(
                &V_B[0],
                vel_B,
                &V_B[4]);
            
            std::vector<double> Q_star_B(5);
            Q_star_B[0] = Chi_star_B*V_B[0];
            Q_star_B[1] = Chi_star_B*Q_B[1];
            Q_star_B[2] = Chi_star_B*V_B[0]*s_star;
            Q_star_B[3] = Chi_star_B*Q_B[3];
            Q_star_B[4] = Chi_star_B*(Q_B[4] + (s_star - V_B[2])*(V_B[0]*s_star + V_B[4]/(s_B - V_B[2])));
            
            std::vector<double> F_y_B(5);
            F_y_B[0] = Q_B[2];
            F_y_B[1] = Q_B[2]*V_B[1];
            F_y_B[2] = Q_B[2]*V_B[2] + V_B[4];
            F_y_B[3] = Q_B[2]*V_B[3];
            F_y_B[4] = V_B[2]*(Q_B[4] + V_B[4]);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                (*F_y_intercell[ei]) = F_y_B[ei] + s_minus*(Q_star_B[ei] - Q_B[ei]);
            }
        }
        else
        {
            const double Chi_star_T = (s_T - V_T[2])/(s_T - s_star);
            
            std::vector<const double*> vel_T(3);
            vel_T[0] = &V_T[1];
            vel_T[1] = &V_T[2];
            vel_T[2] = &V_T[3];
            
            std::vector<double> Q_T(5);
            Q_T[0] = V_T[0];
            Q_T[1] = V_T[0]*V_T[1];
            Q_T[2] = V_T[0]*V_T[2];
            Q_T[3] = V_T[0]*V_T[3];
            Q_T[4] = d_equation_of_state->getTotalEnergy(
                &V_T[0],
                vel_T,
                &V_T[4]);
            
            std::vector<double> Q_star_T(5);
            Q_star_T[0] = Chi_star_T*V_T[0];
            Q_star_T[1] = Chi_star_T*Q_T[1];
            Q_star_T[2] = Chi_star_T*V_T[0]*s_star;
            Q_star_T[3] = Chi_star_T*Q_T[3];
            Q_star_T[4] = Chi_star_T*(Q_T[4] + (s_star - V_T[2])*(V_T[0]*s_star + V_T[4]/(s_T - V_T[2])));
            
            std::vector<double> F_y_T(5);
            F_y_T[0] = Q_T[2];
            F_y_T[1] = Q_T[2]*V_T[1];
            F_y_T[2] = Q_T[2]*V_T[2] + V_T[4];
            F_y_T[3] = Q_T[2]*V_T[3];
            F_y_T[4] = V_T[2]*(Q_T[4] + V_T[4]);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                (*F_y_intercell[ei]) = F_y_T[ei] + s_plus*(Q_star_T[ei] - Q_T[ei]);
            }
        }
    }
}


/*
 * Compute the flux in the z-direction at the intercell face for single-species flow model
 * from primitive variables.
 */
void
RiemannSolverHLLC::computeIntercellFluxForSingleSpeciesInZDirectionFromPrimitiveVariables(
    std::vector<double*>& F_z_intercell,
    std::vector<double>& V_B,
    std::vector<double>& V_F)
{
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no z direction for 1D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no z direction for 2D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const double c_B = d_equation_of_state->getSoundSpeedWithPressure(
            &V_B[0],
            &V_B[4]);
        
        const double c_F = d_equation_of_state->getSoundSpeedWithPressure(
            &V_F[0],
            &V_F[4]);
        
        const double w_average = 0.5*(V_B[3] + V_F[3]);
        const double c_average = 0.5*(c_B + c_F);
        
        const double s_B = fmin(w_average - c_average, V_B[3] - c_B);
        const double s_F = fmax(w_average + c_average, V_F[3] + c_F);
        
        const double s_minus = fmin(0.0, s_B);
        const double s_plus  = fmax(0.0, s_F);
        
        const double s_star =
            (V_F[4] - V_B[4] + V_B[0]*V_B[3]*(s_B - V_B[3]) - V_F[0]*V_F[3]*(s_F - V_F[3]))/(V_B[0]*(s_B - V_B[3]) -
                V_F[0]*(s_F - V_F[3]));
        
        if (s_star > 0)
        {
            const double Chi_star_B = (s_B - V_B[3])/(s_B - s_star);
            
            std::vector<const double*> vel_B(3);
            vel_B[0] = &V_B[1];
            vel_B[1] = &V_B[2];
            vel_B[2] = &V_B[3];
            
            std::vector<double> Q_B(5);
            Q_B[0] = V_B[0];
            Q_B[1] = V_B[0]*V_B[1];
            Q_B[2] = V_B[0]*V_B[2];
            Q_B[3] = V_B[0]*V_B[3];
            Q_B[4] = d_equation_of_state->getTotalEnergy(
                &V_B[0],
                vel_B,
                &V_B[4]);
            
            std::vector<double> Q_star_B(5);
            Q_star_B[0] = Chi_star_B*V_B[0];
            Q_star_B[1] = Chi_star_B*Q_B[1];
            Q_star_B[2] = Chi_star_B*Q_B[2];
            Q_star_B[3] = Chi_star_B*V_B[0]*s_star;
            Q_star_B[4] = Chi_star_B*(Q_B[4] + (s_star - V_B[3])*(V_B[0]*s_star + V_B[4]/(s_B - V_B[3])));
            
            std::vector<double> F_z_B(5);
            F_z_B[0] = Q_B[3];
            F_z_B[1] = Q_B[3]*V_B[1];
            F_z_B[2] = Q_B[3]*V_B[2];
            F_z_B[3] = Q_B[3]*V_B[3] + V_B[4];
            F_z_B[4] = V_B[3]*(Q_B[4] + V_B[4]);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                (*F_z_intercell[ei]) = F_z_B[ei] + s_minus*(Q_star_B[ei] - Q_B[ei]);
            }
        }
        else
        {
            const double Chi_star_F = (s_F - V_F[3])/(s_F - s_star);
            
            std::vector<const double*> vel_F(3);
            vel_F[0] = &V_F[1];
            vel_F[1] = &V_F[2];
            vel_F[2] = &V_F[3];
            
            std::vector<double> Q_F(5);
            Q_F[0] = V_F[0];
            Q_F[1] = V_F[0]*V_F[1];
            Q_F[2] = V_F[0]*V_F[2];
            Q_F[3] = V_F[0]*V_F[3];
            Q_F[4] = d_equation_of_state->getTotalEnergy(
                &V_F[0],
                vel_F,
                &V_F[4]);
            
            std::vector<double> Q_star_F(5);
            Q_star_F[0] = Chi_star_F*V_F[0];
            Q_star_F[1] = Chi_star_F*Q_F[1];
            Q_star_F[2] = Chi_star_F*Q_F[2];
            Q_star_F[3] = Chi_star_F*V_F[0]*s_star;
            Q_star_F[4] = Chi_star_F*(Q_F[4] + (s_star - V_F[3])*(V_F[0]*s_star + V_F[4]/(s_F - V_F[3])));
            
            std::vector<double> F_z_F(5);
            F_z_F[0] = Q_F[3];
            F_z_F[1] = Q_F[3]*V_F[1];
            F_z_F[2] = Q_F[3]*V_F[2];
            F_z_F[3] = Q_F[3]*V_F[3] + V_F[4];
            F_z_F[4] = V_F[3]*(Q_F[4] + V_F[4]);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                (*F_z_intercell[ei]) = F_z_F[ei] + s_plus*(Q_star_F[ei] - Q_F[ei]);
            }
        }
    }
}
