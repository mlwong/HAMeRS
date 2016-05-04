#include "flow_model/Riemann_solver/RiemannSolverSingleSpeciesHLLC.hpp"

/*
 * Compute the flux at the intercell face from conservative variables.
 */
void
RiemannSolverSingleSpeciesHLLC::computeIntercellFluxFromConservativeVariables(
    std::vector<boost::reference_wrapper<double> >& flux_intercell,
    const std::vector<boost::reference_wrapper<double> >& conservative_variables_minus,
    const std::vector<boost::reference_wrapper<double> >& conservative_variables_plus,
    const DIRECTION& direction)
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(flux_intercell.size()) == d_num_eqn);
    TBOX_ASSERT(static_cast<int>(conservative_variables_minus.size()) == d_num_eqn);
    TBOX_ASSERT(static_cast<int>(conservative_variables_plus.size()) == d_num_eqn);
#endif
    
    switch (direction)
    {
        case X_DIRECTION:
        {
            computeIntercellFluxInXDirectionFromConservativeVariables(
                flux_intercell,
                conservative_variables_minus,
                conservative_variables_plus);
            
            break;
        }
        case Y_DIRECTION:
        {
            computeIntercellFluxInYDirectionFromConservativeVariables(
                flux_intercell,
                conservative_variables_minus,
                conservative_variables_plus);
            
            break;
        }
        case Z_DIRECTION:
        {
            computeIntercellFluxInZDirectionFromConservativeVariables(
                flux_intercell,
                conservative_variables_minus,
                conservative_variables_plus);
            
            break;
        }
        default:
        {
            TBOX_ERROR(d_object_name
                << ": RiemannSolverSingleSpeciesHLLC::"
                << "computeIntercellFluxFromConservativeVariables()\n"
                << "Unknown direction."
                << std::endl);
        }
    }
}


/*
 * Compute the flux at the intercell face from primitive variables.
 */
void
RiemannSolverSingleSpeciesHLLC::computeIntercellFluxFromPrimitiveVariables(
    std::vector<boost::reference_wrapper<double> >& flux_intercell,
    const std::vector<boost::reference_wrapper<double> >& primitive_variables_minus,
    const std::vector<boost::reference_wrapper<double> >& primitive_variables_plus,
    const DIRECTION& direction)
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(flux_intercell.size()) == d_num_eqn);
    TBOX_ASSERT(static_cast<int>(primitive_variables_minus.size()) == d_num_eqn);
    TBOX_ASSERT(static_cast<int>(primitive_variables_plus.size()) == d_num_eqn);
#endif
    
    switch (direction)
    {
        case X_DIRECTION:
        {
            computeIntercellFluxInXDirectionFromPrimitiveVariables(
                flux_intercell,
                primitive_variables_minus,
                primitive_variables_plus);
            
            break;
        }
        case Y_DIRECTION:
        {
            computeIntercellFluxInYDirectionFromPrimitiveVariables(
                flux_intercell,
                primitive_variables_minus,
                primitive_variables_plus);
            
            break;
        }
        case Z_DIRECTION:
        {
            computeIntercellFluxInZDirectionFromPrimitiveVariables(
                flux_intercell,
                primitive_variables_minus,
                primitive_variables_plus);
            
            break;
        }
        default:
        {
            TBOX_ERROR(d_object_name
                << ": RiemannSolverSingleSpeciesHLLC::"
                << "computeIntercellFluxFromPrimitiveVariables()\n"
                << "Unknown direction."
                << std::endl);
        }
    }
}


/*
 * Compute the flux in the x-direction at the intercell face
 * from conservative variables.
 */
void
RiemannSolverSingleSpeciesHLLC::computeIntercellFluxInXDirectionFromConservativeVariables(
    std::vector<boost::reference_wrapper<double> >& F_x_intercell,
    const std::vector<boost::reference_wrapper<double> >& Q_L,
    const std::vector<boost::reference_wrapper<double> >& Q_R)
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(F_x_intercell.size()) == d_num_eqn);
    TBOX_ASSERT(static_cast<int>(Q_L.size()) == d_num_eqn);
    TBOX_ASSERT(static_cast<int>(Q_R.size()) == d_num_eqn);
#endif
    
    if (d_dim == tbox::Dimension(1))
    {
        const double u_L = (Q_L[1].get())/(Q_L[0].get());
        const double u_R = (Q_R[1].get())/(Q_R[0].get());
        
        std::vector<const double*> m_L;
        std::vector<const double*> m_R;
        m_L.reserve(1);
        m_R.reserve(1);
        m_L.push_back(&(Q_L[1].get()));
        m_R.push_back(&(Q_R[1].get()));
        
        const double p_L = d_equation_of_state->getPressure(
            &(Q_L[0].get()),
            m_L,
            &(Q_L[2].get()));
        
        const double p_R = d_equation_of_state->getPressure(
            &(Q_R[0].get()),
            m_R,
            &(Q_R[2].get()));
        
        const double c_L = d_equation_of_state->getSoundSpeedWithPressure(
            &(Q_L[0].get()),
            &p_L);
        
        const double c_R = d_equation_of_state->getSoundSpeedWithPressure(
            &(Q_R[0].get()),
            &p_R);
        
        const double u_average = 0.5*(u_L + u_R);
        const double c_average = 0.5*(c_L + c_R);
        
        const double s_L = fmin(u_average - c_average, u_L - c_L);
        const double s_R = fmax(u_average + c_average, u_R + c_R);
        
        const double s_minus = fmin(0.0, s_L);
        const double s_plus  = fmax(0.0, s_R);
        
        const double s_star =
            (p_R - p_L + (Q_L[1].get())*(s_L - u_L) - (Q_R[1].get())*(s_R - u_R))/((Q_L[0].get())*(s_L - u_L) -
                (Q_R[0].get())*(s_R - u_R));
        
        if (s_star > 0)
        {
            const double Chi_star_L = (s_L - u_L)/(s_L - s_star);
            
            double Q_star_L[d_num_eqn];
            Q_star_L[0] = Chi_star_L*(Q_L[0].get());
            Q_star_L[1] = Chi_star_L*(Q_L[0].get())*s_star;
            Q_star_L[2] = Chi_star_L*((Q_L[2].get()) + (s_star - u_L)*((Q_L[0].get())*s_star + p_L/(s_L - u_L)));
            
            double F_x_L[d_num_eqn];
            F_x_L[0] = (Q_L[1].get());
            F_x_L[1] = u_L*(Q_L[1].get()) + p_L;
            F_x_L[2] = u_L*((Q_L[2].get()) + p_L);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_x_intercell[ei].get() = F_x_L[ei] + s_minus*(Q_star_L[ei] - Q_L[ei]);
            }
        }
        else
        {
            const double Chi_star_R = (s_R - u_R)/(s_R - s_star);
            
            double Q_star_R[d_num_eqn];
            Q_star_R[0] = Chi_star_R*(Q_R[0].get());
            Q_star_R[1] = Chi_star_R*(Q_R[0].get())*s_star;
            Q_star_R[2] = Chi_star_R*((Q_R[2].get()) + (s_star - u_R)*((Q_R[0].get())*s_star + p_R/(s_R - u_R)));
            
            double F_x_R[d_num_eqn];
            F_x_R[0] = (Q_R[1].get());
            F_x_R[1] = u_R*(Q_R[1].get()) + p_R;
            F_x_R[2] = u_R*((Q_R[2].get()) + p_R);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_x_intercell[ei].get() = F_x_R[ei] + s_plus*(Q_star_R[ei] - Q_R[ei]);
            }
        }
        
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const double u_L = (Q_L[1].get())/(Q_L[0].get());
        const double u_R = (Q_R[1].get())/(Q_R[0].get());
        
        std::vector<const double*> m_L;
        std::vector<const double*> m_R;
        m_L.reserve(2);
        m_R.reserve(2);
        m_L.push_back(&(Q_L[1].get()));
        m_R.push_back(&(Q_R[1].get()));
        m_L.push_back(&(Q_L[2].get()));
        m_R.push_back(&(Q_R[2].get()));
        
        const double p_L = d_equation_of_state->getPressure(
            &(Q_L[0].get()),
            m_L,
            &(Q_L[3].get()));
        
        const double p_R = d_equation_of_state->getPressure(
            &(Q_R[0].get()),
            m_R,
            &(Q_R[3].get()));
        
        const double c_L = d_equation_of_state->getSoundSpeedWithPressure(
            &(Q_L[0].get()),
            &p_L);
        
        const double c_R = d_equation_of_state->getSoundSpeedWithPressure(
            &(Q_R[0].get()),
            &p_R);
        
        const double u_average = 0.5*(u_L + u_R);
        const double c_average = 0.5*(c_L + c_R);
        
        const double s_L = fmin(u_average - c_average, u_L - c_L);
        const double s_R = fmax(u_average + c_average, u_R + c_R);
        
        const double s_minus = fmin(0.0, s_L);
        const double s_plus  = fmax(0.0, s_R);
        
        const double s_star =
            (p_R - p_L + (Q_L[1].get())*(s_L - u_L) - (Q_R[1].get())*(s_R - u_R))/((Q_L[0].get())*(s_L - u_L) -
                (Q_R[0].get())*(s_R - u_R));
        
        if (s_star > 0)
        {
            const double Chi_star_L = (s_L - u_L)/(s_L - s_star);
            
            double Q_star_L[d_num_eqn];
            Q_star_L[0] = Chi_star_L*(Q_L[0].get());
            Q_star_L[1] = Chi_star_L*(Q_L[0].get())*s_star;
            Q_star_L[2] = Chi_star_L*(Q_L[2].get());
            Q_star_L[3] = Chi_star_L*((Q_L[3].get()) + (s_star - u_L)*((Q_L[0].get())*s_star + p_L/(s_L - u_L)));
            
            double F_x_L[d_num_eqn];
            F_x_L[0] = (Q_L[1].get());
            F_x_L[1] = u_L*(Q_L[1].get()) + p_L;
            F_x_L[2] = u_L*(Q_L[2].get());
            F_x_L[3] = u_L*((Q_L[3].get()) + p_L);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_x_intercell[ei].get() = F_x_L[ei] + s_minus*(Q_star_L[ei] - Q_L[ei]);
            }
        }
        else
        {
            const double Chi_star_R = (s_R - u_R)/(s_R - s_star);
            
            double Q_star_R[d_num_eqn];
            Q_star_R[0] = Chi_star_R*(Q_R[0].get());
            Q_star_R[1] = Chi_star_R*(Q_R[0].get())*s_star;
            Q_star_R[2] = Chi_star_R*(Q_R[2].get());
            Q_star_R[3] = Chi_star_R*((Q_R[3].get()) + (s_star - u_R)*((Q_R[0].get())*s_star + p_R/(s_R - u_R)));
            
            double F_x_R[d_num_eqn];
            F_x_R[0] = (Q_R[1].get());
            F_x_R[1] = u_R*(Q_R[1].get()) + p_R;
            F_x_R[2] = u_R*(Q_R[2].get());
            F_x_R[3] = u_R*((Q_R[3].get()) + p_R);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_x_intercell[ei].get() = F_x_R[ei] + s_plus*(Q_star_R[ei] - Q_R[ei]);
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const double u_L = (Q_L[1].get())/(Q_L[0].get());
        const double u_R = (Q_R[1].get())/(Q_R[0].get());
        
        std::vector<const double*> m_L;
        std::vector<const double*> m_R;
        m_L.reserve(3);
        m_R.reserve(3);
        m_L.push_back(&(Q_L[1].get()));
        m_R.push_back(&(Q_R[1].get()));
        m_L.push_back(&(Q_L[2].get()));
        m_R.push_back(&(Q_R[2].get()));
        m_L.push_back(&(Q_L[3].get()));
        m_R.push_back(&(Q_R[3].get()));
        
        const double p_L = d_equation_of_state->getPressure(
            &(Q_L[0].get()),
            m_L,
            &(Q_L[4].get()));
        
        const double p_R = d_equation_of_state->getPressure(
            &(Q_R[0].get()),
            m_R,
            &(Q_R[4].get()));
        
        const double c_L = d_equation_of_state->getSoundSpeedWithPressure(
            &(Q_L[0].get()),
            &p_L);
        
        const double c_R = d_equation_of_state->getSoundSpeedWithPressure(
            &(Q_R[0].get()),
            &p_R);
        
        const double u_average = 0.5*(u_L + u_R);
        const double c_average = 0.5*(c_L + c_R);
        
        const double s_L = fmin(u_average - c_average, u_L - c_L);
        const double s_R = fmax(u_average + c_average, u_R + c_R);
        
        const double s_minus = fmin(0.0, s_L);
        const double s_plus  = fmax(0.0, s_R);
        
        const double s_star =
            (p_R - p_L + (Q_L[1].get())*(s_L - u_L) - (Q_R[1].get())*(s_R - u_R))/((Q_L[0].get())*(s_L - u_L) -
                (Q_R[0].get())*(s_R - u_R));
        
        if (s_star > 0)
        {
            const double Chi_star_L = (s_L - u_L)/(s_L - s_star);
            
            double Q_star_L[d_num_eqn];
            Q_star_L[0] = Chi_star_L*(Q_L[0].get());
            Q_star_L[1] = Chi_star_L*(Q_L[0].get())*s_star;
            Q_star_L[2] = Chi_star_L*(Q_L[2].get());
            Q_star_L[3] = Chi_star_L*(Q_L[3].get());
            Q_star_L[4] = Chi_star_L*((Q_L[4].get()) + (s_star - u_L)*((Q_L[0].get())*s_star + p_L/(s_L - u_L)));
            
            double F_x_L[d_num_eqn];
            F_x_L[0] = (Q_L[1].get());
            F_x_L[1] = u_L*(Q_L[1].get()) + p_L;
            F_x_L[2] = u_L*(Q_L[2].get());
            F_x_L[3] = u_L*(Q_L[3].get());
            F_x_L[4] = (u_L*((Q_L[4].get()) + p_L));
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_x_intercell[ei].get() = F_x_L[ei] + s_minus*(Q_star_L[ei] - Q_L[ei]);
            }
        }
        else
        {
            const double Chi_star_R = (s_R - u_R)/(s_R - s_star);
            
            double Q_star_R[d_num_eqn];
            Q_star_R[0] = Chi_star_R*(Q_R[0].get());
            Q_star_R[1] = Chi_star_R*(Q_R[0].get())*s_star;
            Q_star_R[2] = Chi_star_R*(Q_R[2].get());
            Q_star_R[3] = Chi_star_R*(Q_R[3].get());
            Q_star_R[4] = Chi_star_R*((Q_R[4].get()) + (s_star - u_R)*((Q_R[0].get())*s_star + p_R/(s_R - u_R)));
            
            double F_x_R[d_num_eqn];
            F_x_R[0] = (Q_R[1].get());
            F_x_R[1] = u_R*(Q_R[1].get()) + p_R;
            F_x_R[2] = u_R*(Q_R[2].get());
            F_x_R[3] = u_R*(Q_R[3].get());
            F_x_R[4] = (u_R*((Q_R[4].get()) + p_R));
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_x_intercell[ei].get() = F_x_R[ei] + s_plus*(Q_star_R[ei] - Q_R[ei]);
            }
        }
    }
}


/*
 * Compute the flux in the y-direction at the intercell face
 * from conservative variables.
 */
void
RiemannSolverSingleSpeciesHLLC::computeIntercellFluxInYDirectionFromConservativeVariables(
    std::vector<boost::reference_wrapper<double> >& F_y_intercell,
    const std::vector<boost::reference_wrapper<double> >& Q_B,
    const std::vector<boost::reference_wrapper<double> >& Q_T)
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(F_y_intercell.size()) == d_num_eqn);
    TBOX_ASSERT(static_cast<int>(Q_B.size()) == d_num_eqn);
    TBOX_ASSERT(static_cast<int>(Q_T.size()) == d_num_eqn);
#endif
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": RiemannSolverSingleSpeciesHLLC::"
            << "computeIntercellFluxInYDirectionFromConservativeVariables()\n"
            << "There is no y direction for 1D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const double v_B = (Q_B[2].get())/(Q_B[0].get());
        const double v_T = (Q_T[2].get())/(Q_T[0].get());
        
        std::vector<const double*> m_B;
        std::vector<const double*> m_T;
        m_B.reserve(2);
        m_T.reserve(2);
        m_B.push_back(&(Q_B[1].get()));
        m_T.push_back(&(Q_T[1].get()));
        m_B.push_back(&(Q_B[2].get()));
        m_T.push_back(&(Q_T[2].get()));
        
        const double p_B = d_equation_of_state->getPressure(
            &(Q_B[0].get()),
            m_B,
            &(Q_B[3].get()));
        
        const double p_T = d_equation_of_state->getPressure(
            &(Q_T[0].get()),
            m_T,
            &(Q_T[3].get()));
        
        const double c_B = d_equation_of_state->getSoundSpeedWithPressure(
            &(Q_B[0].get()),
            &p_B);
        
        const double c_T = d_equation_of_state->getSoundSpeedWithPressure(
            &(Q_T[0].get()),
            &p_T);
        
        const double v_average = 0.5*(v_B + v_T);
        const double c_average = 0.5*(c_B + c_T);
        
        const double s_B = fmin(v_average - c_average, v_B - c_B);
        const double s_T = fmax(v_average + c_average, v_T + c_T);
        
        const double s_minus = fmin(0.0, s_B);
        const double s_plus  = fmax(0.0, s_T);
        
        const double s_star =
            (p_T - p_B + (Q_B[2].get())*(s_B - v_B) - (Q_T[2].get())*(s_T - v_T))/((Q_B[0].get())*(s_B - v_B) -
                (Q_T[0].get())*(s_T - v_T));
        
        if (s_star > 0)
        {
            const double Chi_star_B = (s_B - v_B)/(s_B - s_star);
            
            double Q_star_B[d_num_eqn];
            Q_star_B[0] = Chi_star_B*(Q_B[0].get());
            Q_star_B[1] = Chi_star_B*(Q_B[1].get());
            Q_star_B[2] = Chi_star_B*(Q_B[0].get())*s_star;
            Q_star_B[3] = Chi_star_B*((Q_B[3].get()) + (s_star - v_B)*((Q_B[0].get())*s_star + p_B/(s_B - v_B)));
            
            double F_y_B[d_num_eqn];
            F_y_B[0] = (Q_B[2].get());
            F_y_B[1] = v_B*(Q_B[1].get());
            F_y_B[2] = v_B*(Q_B[2].get()) + p_B;
            F_y_B[3] = v_B*((Q_B[3].get()) + p_B);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_y_intercell[ei].get() = F_y_B[ei] + s_minus*(Q_star_B[ei] - Q_B[ei]);
            }
        }
        else
        {
            const double Chi_star_T = (s_T - v_T)/(s_T - s_star);
            
            double Q_star_T[d_num_eqn];
            Q_star_T[0] = Chi_star_T*(Q_T[0].get());
            Q_star_T[1] = Chi_star_T*(Q_T[1].get());
            Q_star_T[2] = Chi_star_T*(Q_T[0].get())*s_star;
            Q_star_T[3] = Chi_star_T*((Q_T[3].get()) + (s_star - v_T)*((Q_T[0].get())*s_star + p_T/(s_T - v_T)));
            
            double F_y_T[d_num_eqn];
            F_y_T[0] = (Q_T[2].get());
            F_y_T[1] = v_T*(Q_T[1].get());
            F_y_T[2] = v_T*(Q_T[2].get()) + p_T;
            F_y_T[3] = v_T*((Q_T[3].get()) + p_T);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_y_intercell[ei].get() = F_y_T[ei] + s_plus*(Q_star_T[ei] - Q_T[ei]);
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const double v_B = (Q_B[2].get())/(Q_B[0].get());
        const double v_T = (Q_T[2].get())/(Q_T[0].get());
        
        std::vector<const double*> m_B;
        std::vector<const double*> m_T;
        m_B.reserve(3);
        m_T.reserve(3);
        m_B.push_back(&(Q_B[1].get()));
        m_T.push_back(&(Q_T[1].get()));
        m_B.push_back(&(Q_B[2].get()));
        m_T.push_back(&(Q_T[2].get()));
        m_B.push_back(&(Q_B[3].get()));
        m_T.push_back(&(Q_T[3].get()));
        
        const double p_B = d_equation_of_state->getPressure(
            &(Q_B[0].get()),
            m_B,
            &(Q_B[4].get()));
        
        const double p_T = d_equation_of_state->getPressure(
            &(Q_T[0].get()),
            m_T,
            &(Q_T[0].get()));
        
        const double c_B = d_equation_of_state->getSoundSpeedWithPressure(
            &(Q_B[0].get()),
            &p_B);
        
        const double c_T = d_equation_of_state->getSoundSpeedWithPressure(
            &(Q_T[0].get()),
            &p_T);
        
        const double v_average = 0.5*(v_B + v_T);
        const double c_average = 0.5*(c_B + c_T);
        
        const double s_B = fmin(v_average - c_average, v_B - c_B);
        const double s_T = fmax(v_average + c_average, v_T + c_T);
        
        const double s_minus = fmin(0.0, s_B);
        const double s_plus  = fmax(0.0, s_T);
        
        const double s_star =
            (p_T - p_B + (Q_B[2].get())*(s_B - v_B) - (Q_T[2].get())*(s_T - v_T))/((Q_B[0].get())*(s_B - v_B) -
                (Q_T[0].get())*(s_T - v_T));
        
        if (s_star > 0)
        {
            const double Chi_star_B = (s_B - v_B)/(s_B - s_star);
            
            double Q_star_B[d_num_eqn];
            Q_star_B[0] = Chi_star_B*(Q_B[0].get());
            Q_star_B[1] = Chi_star_B*(Q_B[1].get());
            Q_star_B[2] = Chi_star_B*(Q_B[0].get())*s_star;
            Q_star_B[3] = Chi_star_B*(Q_B[3].get());
            Q_star_B[4] = Chi_star_B*((Q_B[4].get()) + (s_star - v_B)*((Q_B[0].get())*s_star + p_B/(s_B - v_B)));
            
            double F_y_B[d_num_eqn];
            F_y_B[0] = (Q_B[2].get());
            F_y_B[1] = v_B*(Q_B[1].get());
            F_y_B[2] = v_B*(Q_B[2].get()) + p_B;
            F_y_B[3] = v_B*(Q_B[3].get());
            F_y_B[4] = v_B*((Q_B[4].get()) + p_B);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_y_intercell[ei].get() = F_y_B[ei] + s_minus*(Q_star_B[ei] - Q_B[ei]);
            }
        }
        else
        {
            const double Chi_star_T = (s_T - v_T)/(s_T - s_star);
            
            double Q_star_T[d_num_eqn];
            Q_star_T[0] = Chi_star_T*(Q_T[0].get());
            Q_star_T[1] = Chi_star_T*(Q_T[1].get());
            Q_star_T[2] = Chi_star_T*(Q_T[0].get())*s_star;
            Q_star_T[3] = Chi_star_T*(Q_T[3].get());
            Q_star_T[4] = Chi_star_T*((Q_T[4].get()) + (s_star - v_T)*((Q_T[0].get())*s_star + p_T/(s_T - v_T)));
            
            double F_y_T[d_num_eqn];
            F_y_T[0] = (Q_T[2].get());
            F_y_T[1] = v_T*(Q_T[1].get());
            F_y_T[2] = v_T*(Q_T[2].get()) + p_T;
            F_y_T[3] = v_T*(Q_T[3].get());
            F_y_T[4] = v_T*((Q_T[4].get()) + p_T);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_y_intercell[ei].get() = F_y_T[ei] + s_plus*(Q_star_T[ei] - Q_T[ei]);
            }
        }
    }
}


/*
 * Compute the flux in the z-direction at the intercell face
 * from conservative variables.
 */
void
RiemannSolverSingleSpeciesHLLC::computeIntercellFluxInZDirectionFromConservativeVariables(
    std::vector<boost::reference_wrapper<double> >& F_z_intercell,
    const std::vector<boost::reference_wrapper<double> >& Q_B,
    const std::vector<boost::reference_wrapper<double> >& Q_F)
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(F_z_intercell.size()) == d_num_eqn);
    TBOX_ASSERT(static_cast<int>(Q_B.size()) == d_num_eqn);
    TBOX_ASSERT(static_cast<int>(Q_F.size()) == d_num_eqn);
#endif
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": RiemannSolverSingleSpeciesHLLC::"
            << "computeIntercellFluxInZDirectionFromConservativeVariables()\n"
            << "There is no z direction for 1D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": RiemannSolverSingleSpeciesHLLC::"
            << "computeIntercellFluxInZDirectionFromConservativeVariables()\n"
            << "There is no z direction for 2D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const double w_B = (Q_B[3].get())/(Q_B[0].get());
        const double w_F = (Q_F[3].get())/(Q_F[0].get());
        
        std::vector<const double*> m_B;
        std::vector<const double*> m_F;
        m_B.reserve(3);
        m_F.reserve(3);
        m_B.push_back(&(Q_B[1].get()));
        m_F.push_back(&(Q_F[1].get()));
        m_B.push_back(&(Q_B[2].get()));
        m_F.push_back(&(Q_F[2].get()));
        m_B.push_back(&(Q_B[3].get()));
        m_F.push_back(&(Q_F[3].get()));
        
        const double p_B = d_equation_of_state->getPressure(
            &(Q_B[0].get()),
            m_B,
            &(Q_B[4].get()));
        
        const double p_F = d_equation_of_state->getPressure(
            &(Q_F[0].get()),
            m_F,
            &(Q_F[4].get()));
        
        const double c_B = d_equation_of_state->getSoundSpeedWithPressure(
            &(Q_B[0].get()),
            &p_B);
        
        const double c_F = d_equation_of_state->getSoundSpeedWithPressure(
            &(Q_F[0].get()),
            &p_F);
        
        const double w_average = 0.5*(w_B + w_F);
        const double c_average = 0.5*(c_B + c_F);
        
        const double s_B = fmin(w_average - c_average, w_B - c_B);
        const double s_F = fmax(w_average + c_average, w_F + c_F);
        
        const double s_minus = fmin(0.0, s_B);
        const double s_plus  = fmax(0.0, s_F);
        
        const double s_star =
            (p_F - p_B + (Q_B[3].get())*(s_B - w_B) - (Q_F[3].get())*(s_F - w_F))/((Q_B[0].get())*(s_B - w_B) -
                (Q_F[0].get())*(s_F - w_F));
        
        if (s_star > 0)
        {
            const double Chi_star_B = (s_B - w_B)/(s_B - s_star);
            
            double Q_star_B[d_num_eqn];
            Q_star_B[0] = Chi_star_B*(Q_B[0].get());
            Q_star_B[1] = Chi_star_B*(Q_B[1].get());
            Q_star_B[2] = Chi_star_B*(Q_B[2].get());
            Q_star_B[3] = Chi_star_B*(Q_B[0].get())*s_star;
            Q_star_B[4] = Chi_star_B*((Q_B[4].get()) + (s_star - w_B)*((Q_B[0].get())*s_star + p_B/(s_B - w_B)));
            
            double F_z_B[d_num_eqn];
            F_z_B[0] = (Q_B[3].get());
            F_z_B[1] = w_B*(Q_B[1].get());
            F_z_B[2] = w_B*(Q_B[2].get());
            F_z_B[3] = w_B*(Q_B[3].get()) + p_B;
            F_z_B[4] = w_B*((Q_B[4].get()) + p_B);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_z_intercell[ei].get() = F_z_B[ei] + s_minus*(Q_star_B[ei] - Q_B[ei]);
            }
        }
        else
        {
            const double Chi_star_F = (s_F - w_F)/(s_F - s_star);
            
            double Q_star_F[d_num_eqn];
            Q_star_F[0] = Chi_star_F*(Q_F[0].get());
            Q_star_F[1] = Chi_star_F*(Q_F[1].get());
            Q_star_F[2] = Chi_star_F*(Q_F[2].get());
            Q_star_F[3] = Chi_star_F*(Q_F[0].get())*s_star;
            Q_star_F[4] = Chi_star_F*((Q_F[4].get()) + (s_star - w_F)*((Q_F[0].get())*s_star + p_F/(s_F - w_F)));
            
            double F_z_F[d_num_eqn];
            F_z_F[0] = (Q_F[3].get());
            F_z_F[1] = w_F*(Q_F[1].get());
            F_z_F[2] = w_F*(Q_F[2].get());
            F_z_F[3] = w_F*(Q_F[3].get()) + p_F;
            F_z_F[4] = w_F*((Q_F[4].get()) + p_F);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_z_intercell[ei].get() = F_z_F[ei] + s_plus*(Q_star_F[ei] - Q_F[ei]);
            }
        }
    }
}


/*
 * Compute the flux in the x-direction at the intercell face
 * from primitive variables.
 */
void
RiemannSolverSingleSpeciesHLLC::computeIntercellFluxInXDirectionFromPrimitiveVariables(
    std::vector<boost::reference_wrapper<double> >& F_x_intercell,
    const std::vector<boost::reference_wrapper<double> >& V_L,
    const std::vector<boost::reference_wrapper<double> >& V_R)
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(F_x_intercell.size()) == d_num_eqn);
    TBOX_ASSERT(static_cast<int>(V_L.size()) == d_num_eqn);
    TBOX_ASSERT(static_cast<int>(V_R.size()) == d_num_eqn);
#endif
    
    if (d_dim == tbox::Dimension(1))
    {
        const double c_L = d_equation_of_state->getSoundSpeedWithPressure(
            &(V_L[0].get()),
            &(V_L[2].get()));
        
        const double c_R = d_equation_of_state->getSoundSpeedWithPressure(
            &(V_R[0].get()),
            &(V_R[2].get()));
        
        const double u_average = 0.5*((V_L[1].get()) + (V_R[1].get()));
        const double c_average = 0.5*(c_L + c_R);
        
        const double s_L = fmin(u_average - c_average, (V_L[1].get()) - c_L);
        const double s_R = fmax(u_average + c_average, (V_R[1].get()) + c_R);
        
        const double s_minus = fmin(0.0, s_L);
        const double s_plus  = fmax(0.0, s_R);
        
        const double s_star =
            ((V_R[2].get()) - (V_L[2].get()) + (V_L[0].get())*(V_L[1].get())*
                (s_L - (V_L[1].get())) - (V_R[0].get())*(V_R[1].get())*(s_R - (V_R[1].get())))/
                    ((V_L[0].get())*(s_L - (V_L[1].get())) - (V_R[0].get())*(s_R - (V_R[1].get())));
        
        if (s_star > 0)
        {
            const double Chi_star_L = (s_L - (V_L[1].get()))/(s_L - s_star);
            
            std::vector<const double*> vel_L;
            vel_L.reserve(1);
            vel_L.push_back(&(V_L[1].get()));
            
            double Q_L[d_num_eqn];
            Q_L[0] = (V_L[0].get());
            Q_L[1] = (V_L[0].get())*(V_L[1].get());
            Q_L[2] = d_equation_of_state->getTotalEnergy(
                &(V_L[0].get()),
                vel_L,
                &(V_L[2].get()));
            
            double Q_star_L[d_num_eqn];
            Q_star_L[0] = Chi_star_L*(V_L[0].get());
            Q_star_L[1] = Chi_star_L*(V_L[0].get())*s_star;
            Q_star_L[2] = Chi_star_L*(Q_L[2] + (s_star - (V_L[1].get()))*
                ((V_L[0].get())*s_star + (V_L[2].get())/(s_L - (V_L[1].get()))));
            
            double F_x_L[d_num_eqn];
            F_x_L[0] = Q_L[1];
            F_x_L[1] = (V_L[1].get())*Q_L[1] + (V_L[2].get());
            F_x_L[2] = (V_L[1].get())*(Q_L[2] + (V_L[2].get()));
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_x_intercell[ei].get() = F_x_L[ei] + s_minus*(Q_star_L[ei] - Q_L[ei]);
            }
        }
        else
        {
            const double Chi_star_R = (s_R - (V_R[1].get()))/(s_R - s_star);
            
            std::vector<const double*> vel_R;
            vel_R.reserve(1);
            vel_R[0] = &(V_R[1].get());
            
            double Q_R[d_num_eqn];
            Q_R[0] = (V_R[0].get());
            Q_R[1] = (V_R[0].get())*(V_R[1].get());
            Q_R[2] = d_equation_of_state->getTotalEnergy(
                &(V_R[0].get()),
                vel_R,
                &(V_R[2].get()));
            
            double Q_star_R[d_num_eqn];
            Q_star_R[0] = Chi_star_R*(V_R[0].get());
            Q_star_R[1] = Chi_star_R*(V_R[0].get())*s_star;
            Q_star_R[2] = Chi_star_R*(Q_R[2] + (s_star - (V_R[1].get()))*
                ((V_R[0].get())*s_star + (V_R[2].get())/(s_R - (V_R[1].get()))));
            
            double F_x_R[d_num_eqn];
            F_x_R[0] = Q_R[1];
            F_x_R[1] = (V_R[1].get())*Q_R[1] + (V_R[2].get());
            F_x_R[2] = (V_R[1].get())*(Q_R[2] + (V_R[2].get()));
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_x_intercell[ei].get() = F_x_R[ei] + s_plus*(Q_star_R[ei] - Q_R[ei]);
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const double c_L = d_equation_of_state->getSoundSpeedWithPressure(
            &(V_L[0].get()),
            &(V_L[3].get()));
        
        const double c_R = d_equation_of_state->getSoundSpeedWithPressure(
            &(V_R[0].get()),
            &(V_R[3].get()));
        
        const double u_average = 0.5*((V_L[1].get()) + (V_R[1].get()));
        const double c_average = 0.5*(c_L + c_R);
        
        const double s_L = fmin(u_average - c_average, (V_L[1].get()) - c_L);
        const double s_R = fmax(u_average + c_average, (V_R[1].get()) + c_R);
        
        const double s_minus = fmin(0.0, s_L);
        const double s_plus  = fmax(0.0, s_R);
        
        const double s_star =
            ((V_R[3].get()) - (V_L[3].get()) + (V_L[0].get())*(V_L[1].get())*
                (s_L - (V_L[1].get())) - (V_R[0].get())*(V_R[1].get())*(s_R - (V_R[1].get())))/
                    ((V_L[0].get())*(s_L - (V_L[1].get())) - (V_R[0].get())*(s_R - (V_R[1].get())));
        
        if (s_star > 0)
        {
            const double Chi_star_L = (s_L - (V_L[1].get()))/(s_L - s_star);
            
            std::vector<const double*> vel_L;
            vel_L.reserve(2);
            vel_L.push_back(&(V_L[1].get()));
            vel_L.push_back(&(V_L[2].get()));
            
            double Q_L[d_num_eqn];
            Q_L[0] = (V_L[0].get());
            Q_L[1] = (V_L[0].get())*(V_L[1].get());
            Q_L[2] = (V_L[0].get())*(V_L[2].get());
            Q_L[3] = d_equation_of_state->getTotalEnergy(
                &(V_L[0].get()),
                vel_L,
                &(V_L[3].get()));
            
            double Q_star_L[d_num_eqn];
            Q_star_L[0] = Chi_star_L*(V_L[0].get());
            Q_star_L[1] = Chi_star_L*(V_L[0].get())*s_star;
            Q_star_L[2] = Chi_star_L*Q_L[2];
            Q_star_L[3] = Chi_star_L*(Q_L[3] + (s_star - (V_L[1].get()))*
                ((V_L[0].get())*s_star + (V_L[3].get())/(s_L - (V_L[1].get()))));
            
            double F_x_L[d_num_eqn];
            F_x_L[0] = Q_L[1];
            F_x_L[1] = Q_L[1]*(V_L[1].get()) + (V_L[3].get());
            F_x_L[2] = Q_L[1]*(V_L[2].get());
            F_x_L[3] = (V_L[1].get())*(Q_L[3] + (V_L[3].get()));
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_x_intercell[ei].get() = F_x_L[ei] + s_minus*(Q_star_L[ei] - Q_L[ei]);
            }
        }
        else
        {
            const double Chi_star_R = (s_R - (V_R[1].get()))/(s_R - s_star);
            
            std::vector<const double*> vel_R;
            vel_R.reserve(2);
            vel_R.push_back(&(V_R[1].get()));
            vel_R.push_back(&(V_R[2].get()));
            
            double Q_R[d_num_eqn];
            Q_R[0] = (V_R[0].get());
            Q_R[1] = (V_R[0].get())*(V_R[1].get());
            Q_R[2] = (V_R[0].get())*(V_R[2].get());
            Q_R[3] = d_equation_of_state->getTotalEnergy(
                &(V_R[0].get()),
                vel_R,
                &(V_R[3].get()));
            
            double Q_star_R[d_num_eqn];
            Q_star_R[0] = Chi_star_R*(V_R[0].get());
            Q_star_R[1] = Chi_star_R*(V_R[0].get())*s_star;
            Q_star_R[2] = Chi_star_R*Q_R[2];
            Q_star_R[3] = Chi_star_R*(Q_R[3] + (s_star - (V_R[1].get()))*
                ((V_R[0].get())*s_star + (V_R[3].get())/(s_R - (V_R[1].get()))));
            
            double F_x_R[d_num_eqn];
            F_x_R[0] = Q_R[1];
            F_x_R[1] = Q_R[1]*(V_R[1].get()) + (V_R[3].get());
            F_x_R[2] = Q_R[1]*(V_R[2].get());
            F_x_R[3] = (V_R[1].get())*(Q_R[3] + (V_R[3].get()));
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_x_intercell[ei].get() = F_x_R[ei] + s_plus*(Q_star_R[ei] - Q_R[ei]);
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const double c_L = d_equation_of_state->getSoundSpeedWithPressure(
            &(V_L[0].get()),
            &(V_L[4].get()));
        
        const double c_R = d_equation_of_state->getSoundSpeedWithPressure(
            &(V_R[0].get()),
            &(V_R[4].get()));
        
        const double u_average = 0.5*((V_L[1].get()) + (V_R[1].get()));
        const double c_average = 0.5*(c_L + c_R);
        
        const double s_L = fmin(u_average - c_average, (V_L[1].get()) - c_L);
        const double s_R = fmax(u_average + c_average, (V_R[1].get()) + c_R);
        
        const double s_minus = fmin(0.0, s_L);
        const double s_plus  = fmax(0.0, s_R);
        
        const double s_star =
            ((V_R[4].get()) - (V_L[4].get()) + (V_L[0].get())*(V_L[1].get())*
                (s_L - (V_L[1].get())) - (V_R[0].get())*(V_R[1].get())*(s_R - (V_R[1].get())))/
                    ((V_L[0].get())*(s_L - (V_L[1].get())) - (V_R[0].get())*(s_R - (V_R[1].get())));
        
        if (s_star > 0)
        {
            const double Chi_star_L = (s_L - (V_L[1].get()))/(s_L - s_star);
            
            std::vector<const double*> vel_L;
            vel_L.reserve(3);
            vel_L.push_back(&(V_L[1].get()));
            vel_L.push_back(&(V_L[2].get()));
            vel_L.push_back(&(V_L[3].get()));
            
            double Q_L[d_num_eqn];
            Q_L[0] = (V_L[0].get());
            Q_L[1] = (V_L[0].get())*(V_L[1].get());
            Q_L[2] = (V_L[0].get())*(V_L[2].get());
            Q_L[3] = (V_L[0].get())*(V_L[3].get());
            Q_L[4] = d_equation_of_state->getTotalEnergy(
                &(V_L[0].get()),
                vel_L,
                &(V_L[4].get()));
            
            double Q_star_L[d_num_eqn];
            Q_star_L[0] = Chi_star_L*(V_L[0].get());
            Q_star_L[1] = Chi_star_L*(V_L[0].get())*s_star;
            Q_star_L[2] = Chi_star_L*Q_L[2];
            Q_star_L[3] = Chi_star_L*Q_L[3];
            Q_star_L[4] = Chi_star_L*(Q_L[4] + (s_star - (V_L[1].get()))*
                ((V_L[0].get())*s_star + (V_L[4].get())/(s_L - (V_L[1].get()))));
            
            double F_x_L[d_num_eqn];
            F_x_L[0] = Q_L[1];
            F_x_L[1] = Q_L[1]*(V_L[1].get()) + (V_L[4].get());
            F_x_L[2] = Q_L[1]*(V_L[2].get());
            F_x_L[3] = Q_L[1]*(V_L[3].get());
            F_x_L[4] = (V_L[1].get())*(Q_L[4] + (V_L[4].get()));
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_x_intercell[ei].get() = F_x_L[ei] + s_minus*(Q_star_L[ei] - Q_L[ei]);
            }
        }
        else
        {
            const double Chi_star_R = (s_R - (V_R[1].get()))/(s_R - s_star);
            
            std::vector<const double*> vel_R;
            vel_R.reserve(3);
            vel_R.push_back(&(V_R[1].get()));
            vel_R.push_back(&(V_R[2].get()));
            vel_R.push_back(&(V_R[3].get()));
            
            double Q_R[d_num_eqn];
            Q_R[0] = (V_R[0].get());
            Q_R[1] = (V_R[0].get())*(V_R[1].get());
            Q_R[2] = (V_R[0].get())*(V_R[2].get());
            Q_R[3] = (V_R[0].get())*(V_R[3].get());
            Q_R[4] = d_equation_of_state->getTotalEnergy(
                &(V_R[0].get()),
                vel_R,
                &(V_R[4].get()));
            
            double Q_star_R[d_num_eqn];
            Q_star_R[0] = Chi_star_R*(V_R[0].get());
            Q_star_R[1] = Chi_star_R*(V_R[0].get())*s_star;
            Q_star_R[2] = Chi_star_R*Q_R[2];
            Q_star_R[3] = Chi_star_R*Q_R[3];
            Q_star_R[4] = Chi_star_R*(Q_R[4] + (s_star - (V_R[1].get()))*
                ((V_R[0].get())*s_star + (V_R[4].get())/(s_R - (V_R[1].get()))));
            
            double F_x_R[d_num_eqn];
            F_x_R[0] = Q_R[1];
            F_x_R[1] = Q_R[1]*(V_R[1].get()) + (V_R[4].get());
            F_x_R[2] = Q_R[1]*(V_R[2].get());
            F_x_R[3] = Q_R[1]*(V_R[3].get());
            F_x_R[4] = (V_R[1].get())*(Q_R[4] + (V_R[4].get()));
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_x_intercell[ei].get() = F_x_R[ei] + s_plus*(Q_star_R[ei] - Q_R[ei]);
            }
        }
    }
}


/*
 * Compute the flux in the y-direction at the intercell face
 * from primitive variables.
 */
void
RiemannSolverSingleSpeciesHLLC::computeIntercellFluxInYDirectionFromPrimitiveVariables(
    std::vector<boost::reference_wrapper<double> >& F_y_intercell,
    const std::vector<boost::reference_wrapper<double> >& V_B,
    const std::vector<boost::reference_wrapper<double> >& V_T)
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(F_y_intercell.size()) == d_num_eqn);
    TBOX_ASSERT(static_cast<int>(V_B.size()) == d_num_eqn);
    TBOX_ASSERT(static_cast<int>(V_T.size()) == d_num_eqn);
#endif
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": RiemannSolverSingleSpeciesHLLC::"
            << "computeIntercellFluxInYDirectionFromPrimitiveVariables()\n"
            << "There is no y direction for 1D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const double c_B = d_equation_of_state->getSoundSpeedWithPressure(
            &(V_B[0].get()),
            &(V_B[3].get()));
        
        const double c_T = d_equation_of_state->getSoundSpeedWithPressure(
            &(V_T[0].get()),
            &(V_T[3].get()));
        
        const double v_average = 0.5*((V_B[2].get()) + (V_T[2].get()));
        const double c_average = 0.5*(c_B + c_T);
        
        const double s_B = fmin(v_average - c_average, (V_B[2].get()) - c_B);
        const double s_T = fmax(v_average + c_average, (V_T[2].get()) + c_T);
        
        const double s_minus = fmin(0.0, s_B);
        const double s_plus  = fmax(0.0, s_T);
        
        const double s_star =
            ((V_T[3].get()) - (V_B[3].get()) + (V_B[0].get())*(V_B[2].get())*
                (s_B - (V_B[2].get())) - (V_T[0].get())*(V_T[2].get())*(s_T - (V_T[2].get())))/
                    ((V_B[0].get())*(s_B - (V_B[2].get())) - (V_T[0].get())*(s_T - (V_T[2].get())));
        
        if (s_star > 0)
        {
            const double Chi_star_B = (s_B - (V_B[2].get()))/(s_B - s_star);
            
            std::vector<const double*> vel_B;
            vel_B.reserve(2);
            vel_B.push_back(&(V_B[1].get()));
            vel_B.push_back(&(V_B[2].get()));
            
            double Q_B[d_num_eqn];
            Q_B[0] = (V_B[0].get());
            Q_B[1] = (V_B[0].get())*(V_B[1].get());
            Q_B[2] = (V_B[0].get())*(V_B[2].get());
            Q_B[3] = d_equation_of_state->getTotalEnergy(
                &(V_B[0].get()),
                vel_B,
                &(V_B[3].get()));
            
            double Q_star_B[d_num_eqn];
            Q_star_B[0] = Chi_star_B*(V_B[0].get());
            Q_star_B[1] = Chi_star_B*Q_B[1];
            Q_star_B[2] = Chi_star_B*(V_B[0].get())*s_star;
            Q_star_B[3] = Chi_star_B*(Q_B[3] + (s_star - (V_B[2].get()))*
                ((V_B[0].get())*s_star + (V_B[3].get())/(s_B - (V_B[2].get()))));
            
            double F_y_B[d_num_eqn];
            F_y_B[0] = Q_B[2];
            F_y_B[1] = Q_B[2]*(V_B[1].get());
            F_y_B[2] = Q_B[2]*(V_B[2].get()) + (V_B[3].get());
            F_y_B[3] = (V_B[2].get())*(Q_B[3] + (V_B[3].get()));
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_y_intercell[ei].get() = F_y_B[ei] + s_minus*(Q_star_B[ei] - Q_B[ei]);
            }
        }
        else
        {
            const double Chi_star_T = (s_T - (V_T[2].get()))/(s_T - s_star);
            
            std::vector<const double*> vel_T;
            vel_T.reserve(2);
            vel_T.push_back(&(V_T[1].get()));
            vel_T.push_back(&(V_T[2].get()));
            
            double Q_T[d_num_eqn];
            Q_T[0] = (V_T[0].get());
            Q_T[1] = (V_T[0].get())*(V_T[1].get());
            Q_T[2] = (V_T[0].get())*(V_T[2].get());
            Q_T[3] = d_equation_of_state->getTotalEnergy(
                &(V_T[0].get()),
                vel_T,
                &(V_T[3].get()));
            
            double Q_star_T[d_num_eqn];
            Q_star_T[0] = Chi_star_T*(V_T[0].get());
            Q_star_T[1] = Chi_star_T*Q_T[1];
            Q_star_T[2] = Chi_star_T*(V_T[0].get())*s_star;
            Q_star_T[3] = Chi_star_T*(Q_T[3] + (s_star - (V_T[2].get()))*
                ((V_T[0].get())*s_star + (V_T[3].get())/(s_T - (V_T[2].get()))));
            
            double F_y_T[d_num_eqn];
            F_y_T[0] = Q_T[2];
            F_y_T[1] = Q_T[2]*(V_T[1].get());
            F_y_T[2] = Q_T[2]*(V_T[2].get()) + (V_T[3].get());
            F_y_T[3] = (V_T[2].get())*(Q_T[3] + (V_T[3].get()));
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_y_intercell[ei].get() = F_y_T[ei] + s_plus*(Q_star_T[ei] - Q_T[ei]);
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const double c_B = d_equation_of_state->getSoundSpeedWithPressure(
            &(V_B[0].get()),
            &(V_B[4].get()));
        
        const double c_T = d_equation_of_state->getSoundSpeedWithPressure(
            &(V_T[0].get()),
            &(V_T[4].get()));
        
        const double v_average = 0.5*((V_B[2].get()) + (V_T[2].get()));
        const double c_average = 0.5*(c_B + c_T);
        
        const double s_B = fmin(v_average - c_average, (V_B[2].get()) - c_B);
        const double s_T = fmax(v_average + c_average, (V_T[2].get()) + c_T);
        
        const double s_minus = fmin(0.0, s_B);
        const double s_plus  = fmax(0.0, s_T);
        
        const double s_star =
            ((V_T[4].get()) - (V_B[4].get()) + (V_B[0].get())*(V_B[2].get())*
                (s_B - (V_B[2].get())) - (V_T[0].get())*(V_T[2].get())*(s_T - (V_T[2].get())))/
                    ((V_B[0].get())*(s_B - (V_B[2].get())) - (V_T[0].get())*(s_T - (V_T[2].get())));
        
        if (s_star > 0)
        {
            const double Chi_star_B = (s_B - (V_B[2].get()))/(s_B - s_star);
            
            std::vector<const double*> vel_B;
            vel_B.reserve(3);
            vel_B.push_back(&(V_B[1].get()));
            vel_B.push_back(&(V_B[2].get()));
            vel_B.push_back(&(V_B[3].get()));
            
            double Q_B[d_num_eqn];
            Q_B[0] = (V_B[0].get());
            Q_B[1] = (V_B[0].get())*(V_B[1].get());
            Q_B[2] = (V_B[0].get())*(V_B[2].get());
            Q_B[3] = (V_B[0].get())*(V_B[3].get());
            Q_B[4] = d_equation_of_state->getTotalEnergy(
                &(V_B[0].get()),
                vel_B,
                &(V_B[4].get()));
            
            double Q_star_B[d_num_eqn];
            Q_star_B[0] = Chi_star_B*(V_B[0].get());
            Q_star_B[1] = Chi_star_B*Q_B[1];
            Q_star_B[2] = Chi_star_B*(V_B[0].get())*s_star;
            Q_star_B[3] = Chi_star_B*Q_B[3];
            Q_star_B[4] = Chi_star_B*(Q_B[4] + (s_star - (V_B[2].get()))*
                ((V_B[0].get())*s_star + (V_B[4].get())/(s_B - (V_B[2].get()))));
            
            double F_y_B[d_num_eqn];
            F_y_B[0] = Q_B[2];
            F_y_B[1] = Q_B[2]*(V_B[1].get());
            F_y_B[2] = Q_B[2]*(V_B[2].get()) + (V_B[4].get());
            F_y_B[3] = Q_B[2]*(V_B[3].get());
            F_y_B[4] = (V_B[2].get())*(Q_B[4] + (V_B[4].get()));
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_y_intercell[ei].get() = F_y_B[ei] + s_minus*(Q_star_B[ei] - Q_B[ei]);
            }
        }
        else
        {
            const double Chi_star_T = (s_T - (V_T[2].get()))/(s_T - s_star);
            
            std::vector<const double*> vel_T;
            vel_T.reserve(3);
            vel_T.push_back(&(V_T[1].get()));
            vel_T.push_back(&(V_T[2].get()));
            vel_T.push_back(&(V_T[3].get()));
            
            double Q_T[d_num_eqn];
            Q_T[0] = (V_T[0].get());
            Q_T[1] = (V_T[0].get())*(V_T[1].get());
            Q_T[2] = (V_T[0].get())*(V_T[2].get());
            Q_T[3] = (V_T[0].get())*(V_T[3].get());
            Q_T[4] = d_equation_of_state->getTotalEnergy(
                &(V_T[0].get()),
                vel_T,
                &(V_T[4].get()));
            
            double Q_star_T[d_num_eqn];
            Q_star_T[0] = Chi_star_T*(V_T[0].get());
            Q_star_T[1] = Chi_star_T*Q_T[1];
            Q_star_T[2] = Chi_star_T*(V_T[0].get())*s_star;
            Q_star_T[3] = Chi_star_T*Q_T[3];
            Q_star_T[4] = Chi_star_T*(Q_T[4] + (s_star - (V_T[2].get()))*
                ((V_T[0].get())*s_star + (V_T[4].get())/(s_T - (V_T[2].get()))));
            
            double F_y_T[d_num_eqn];
            F_y_T[0] = Q_T[2];
            F_y_T[1] = Q_T[2]*(V_T[1].get());
            F_y_T[2] = Q_T[2]*(V_T[2].get()) + (V_T[4].get());
            F_y_T[3] = Q_T[2]*(V_T[3].get());
            F_y_T[4] = (V_T[2].get())*(Q_T[4] + (V_T[4].get()));
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_y_intercell[ei].get() = F_y_T[ei] + s_plus*(Q_star_T[ei] - Q_T[ei]);
            }
        }
    }
}


/*
 * Compute the flux in the z-direction at the intercell face
 * from primitive variables.
 */
void
RiemannSolverSingleSpeciesHLLC::computeIntercellFluxInZDirectionFromPrimitiveVariables(
    std::vector<boost::reference_wrapper<double> >& F_z_intercell,
    const std::vector<boost::reference_wrapper<double> >& V_B,
    const std::vector<boost::reference_wrapper<double> >& V_F)
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(F_z_intercell.size()) == d_num_eqn);
    TBOX_ASSERT(static_cast<int>(V_B.size()) == d_num_eqn);
    TBOX_ASSERT(static_cast<int>(V_F.size()) == d_num_eqn);
#endif
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": RiemannSolverSingleSpeciesHLLC::"
            << "computeIntercellFluxInZDirectionFromPrimitiveVariables()\n"
            << "There is no z direction for 1D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": RiemannSolverSingleSpeciesHLLC::"
            << "computeIntercellFluxInZDirectionFromPrimitiveVariables()\n"
            << "There is no z direction for 2D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const double c_B = d_equation_of_state->getSoundSpeedWithPressure(
            &(V_B[0].get()),
            &(V_B[4].get()));
        
        const double c_F = d_equation_of_state->getSoundSpeedWithPressure(
            &(V_F[0].get()),
            &(V_F[4].get()));
        
        const double w_average = 0.5*((V_B[3].get()) + (V_F[3].get()));
        const double c_average = 0.5*(c_B + c_F);
        
        const double s_B = fmin(w_average - c_average, (V_B[3].get()) - c_B);
        const double s_F = fmax(w_average + c_average, (V_F[3].get()) + c_F);
        
        const double s_minus = fmin(0.0, s_B);
        const double s_plus  = fmax(0.0, s_F);
        
        const double s_star =
            ((V_F[4].get()) - (V_B[4].get()) + (V_B[0].get())*(V_B[3].get())*
                (s_B - (V_B[3].get())) - (V_F[0].get())*(V_F[3].get())*(s_F - (V_F[3].get())))/
                    ((V_B[0].get())*(s_B - (V_B[3].get())) - (V_F[0].get())*(s_F - (V_F[3].get())));
        
        if (s_star > 0)
        {
            const double Chi_star_B = (s_B - (V_B[3].get()))/(s_B - s_star);
            
            std::vector<const double*> vel_B;
            vel_B.reserve(3);
            vel_B.push_back(&(V_B[1].get()));
            vel_B.push_back(&(V_B[2].get()));
            vel_B.push_back(&(V_B[3].get()));
            
            double Q_B[d_num_eqn];
            Q_B[0] = (V_B[0].get());
            Q_B[1] = (V_B[0].get())*(V_B[1].get());
            Q_B[2] = (V_B[0].get())*(V_B[2].get());
            Q_B[3] = (V_B[0].get())*(V_B[3].get());
            Q_B[4] = d_equation_of_state->getTotalEnergy(
                &(V_B[0].get()),
                vel_B,
                &(V_B[4].get()));
            
            double Q_star_B[d_num_eqn];
            Q_star_B[0] = Chi_star_B*(V_B[0].get());
            Q_star_B[1] = Chi_star_B*Q_B[1];
            Q_star_B[2] = Chi_star_B*Q_B[2];
            Q_star_B[3] = Chi_star_B*(V_B[0].get())*s_star;
            Q_star_B[4] = Chi_star_B*(Q_B[4] + (s_star - (V_B[3].get()))*
                ((V_B[0].get())*s_star + (V_B[4].get())/(s_B - (V_B[3].get()))));
            
            double F_z_B[d_num_eqn];
            F_z_B[0] = Q_B[3];
            F_z_B[1] = Q_B[3]*(V_B[1].get());
            F_z_B[2] = Q_B[3]*(V_B[2].get());
            F_z_B[3] = Q_B[3]*(V_B[3].get()) + (V_B[4].get());
            F_z_B[4] = (V_B[3].get())*(Q_B[4] + (V_B[4].get()));
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_z_intercell[ei].get() = F_z_B[ei] + s_minus*(Q_star_B[ei] - Q_B[ei]);
            }
        }
        else
        {
            const double Chi_star_F = (s_F - (V_F[3].get()))/(s_F - s_star);
            
            std::vector<const double*> vel_F;
            vel_F.reserve(3);
            vel_F.push_back(&(V_F[1].get()));
            vel_F.push_back(&(V_F[2].get()));
            vel_F.push_back(&(V_F[3].get()));
            
            double Q_F[d_num_eqn];
            Q_F[0] = (V_F[0].get());
            Q_F[1] = (V_F[0].get())*(V_F[1].get());
            Q_F[2] = (V_F[0].get())*(V_F[2].get());
            Q_F[3] = (V_F[0].get())*(V_F[3].get());
            Q_F[4] = d_equation_of_state->getTotalEnergy(
                &(V_F[0].get()),
                vel_F,
                &(V_F[4].get()));
            
            double Q_star_F[d_num_eqn];
            Q_star_F[0] = Chi_star_F*(V_F[0].get());
            Q_star_F[1] = Chi_star_F*Q_F[1];
            Q_star_F[2] = Chi_star_F*Q_F[2];
            Q_star_F[3] = Chi_star_F*(V_F[0].get())*s_star;
            Q_star_F[4] = Chi_star_F*(Q_F[4] + (s_star - (V_F[3].get()))*
                ((V_F[0].get())*s_star + (V_F[4].get())/(s_F - (V_F[3].get()))));
            
            double F_z_F[d_num_eqn];
            F_z_F[0] = Q_F[3];
            F_z_F[1] = Q_F[3]*(V_F[1].get());
            F_z_F[2] = Q_F[3]*(V_F[2].get());
            F_z_F[3] = Q_F[3]*(V_F[3].get()) + (V_F[4].get());
            F_z_F[4] = (V_F[3].get())*(Q_F[4] + (V_F[4].get()));
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_z_intercell[ei].get() = F_z_F[ei] + s_plus*(Q_star_F[ei] - Q_F[ei]);
            }
        }
    }
}
