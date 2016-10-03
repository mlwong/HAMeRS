#include "flow/flow_models/Riemann_solvers/RiemannSolverFourEqnConservativeHLLC.hpp"

/*
 * Compute the flux at the intercell face from conservative variables.
 */
void
RiemannSolverFourEqnConservativeHLLC::computeIntercellFluxFromConservativeVariables(
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
                << ": RiemannSolverFourEqnConservativeHLLC::"
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
RiemannSolverFourEqnConservativeHLLC::computeIntercellFluxFromPrimitiveVariables(
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
                << ": RiemannSolverFourEqnConservativeHLLC::"
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
RiemannSolverFourEqnConservativeHLLC::computeIntercellFluxInXDirectionFromConservativeVariables(
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
        std::vector<const double*> rho_Y_L;
        std::vector<const double*> rho_Y_R;
        rho_Y_L.reserve(d_num_species);
        rho_Y_R.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y_L.push_back(&(Q_L[si].get()));
            rho_Y_R.push_back(&(Q_R[si].get()));
        }
        
        const double rho_L = d_equation_of_state_mixing_rules->getMixtureDensity(
            rho_Y_L);
        
        const double rho_R = d_equation_of_state_mixing_rules->getMixtureDensity(
            rho_Y_R);
        
        const double u_L = (Q_L[d_num_species].get())/rho_L;
        const double u_R = (Q_R[d_num_species].get())/rho_R;
        
        std::vector<const double*> m_L;
        std::vector<const double*> m_R;
        m_L.reserve(1);
        m_R.reserve(1);
        m_L.push_back(&(Q_L[d_num_species].get()));
        m_R.push_back(&(Q_R[d_num_species].get()));
        
        std::vector<const double*> vel_L;
        std::vector<const double*> vel_R;
        vel_L.reserve(1);
        vel_R.reserve(1);
        vel_L.push_back(&u_L);
        vel_R.push_back(&u_R);
        
        std::vector<const double*> empty_ptr;
        
        const double p_L = d_equation_of_state_mixing_rules->getPressure(
            rho_Y_L,
            m_L,
            &(Q_L[d_num_species + d_dim.getValue()].get()),
            empty_ptr);
        
        const double p_R = d_equation_of_state_mixing_rules->getPressure(
            rho_Y_R,
            m_R,
            &(Q_R[d_num_species + d_dim.getValue()].get()),
            empty_ptr);
        
        const double c_L = d_equation_of_state_mixing_rules->getSoundSpeed(
            rho_Y_L,
            vel_L,
            &p_L,
            empty_ptr);
        
        const double c_R = d_equation_of_state_mixing_rules->getSoundSpeed(
            rho_Y_R,
            vel_R,
            &p_R,
            empty_ptr);
        
        const double u_average = 0.5*(u_L + u_R);
        const double c_average = 0.5*(c_L + c_R);
        
        const double s_L = fmin(u_average - c_average, u_L - c_L);
        const double s_R = fmax(u_average + c_average, u_R + c_R);
        
        const double s_minus = fmin(0.0, s_L);
        const double s_plus  = fmax(0.0, s_R);
        
        const double s_star =
            (p_R - p_L + (Q_L[d_num_species].get())*(s_L - u_L) - (Q_R[d_num_species].get())*(s_R - u_R))/
                (rho_L*(s_L - u_L) - rho_R*(s_R - u_R));
        
        if (s_star > 0)
        {
            const double Chi_star_L = (s_L - u_L)/(s_L - s_star);
            
            double Q_star_L[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                Q_star_L[si] = Chi_star_L*(Q_L[si].get());
            }
            Q_star_L[d_num_species] = Chi_star_L*rho_L*s_star;
            Q_star_L[d_num_species + 1] = Chi_star_L*((Q_L[d_num_species + 1].get()) +
                (s_star - u_L)*(rho_L*s_star + p_L/(s_L - u_L)));
            
            double F_x_L[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                F_x_L[si] = u_L*(Q_L[si].get());
            }
            F_x_L[d_num_species] = u_L*(Q_L[d_num_species].get()) + p_L;
            F_x_L[d_num_species + 1] = u_L*((Q_L[d_num_species + 1].get()) + p_L);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_x_intercell[ei].get() = F_x_L[ei] + s_minus*(Q_star_L[ei] - Q_L[ei]);
            }
        }
        else
        {
            const double Chi_star_R = (s_R - u_R)/(s_R - s_star);
            
            double Q_star_R[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                Q_star_R[si] = Chi_star_R*(Q_R[si].get());
            }
            Q_star_R[d_num_species] = Chi_star_R*rho_R*s_star;
            Q_star_R[d_num_species + 1] = Chi_star_R*((Q_R[d_num_species + 1].get()) +
                (s_star - u_R)*(rho_R*s_star + p_R/(s_R - u_R)));
            
            double F_x_R[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                F_x_R[si] = u_R*(Q_R[si].get());
            }
            F_x_R[d_num_species] = u_R*(Q_R[d_num_species].get()) + p_R;
            F_x_R[d_num_species + 1] = u_R*((Q_R[d_num_species + 1].get()) + p_R);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_x_intercell[ei].get() = F_x_R[ei] + s_plus*(Q_star_R[ei] - Q_R[ei]);
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        std::vector<const double*> rho_Y_L;
        std::vector<const double*> rho_Y_R;
        rho_Y_L.reserve(d_num_species);
        rho_Y_R.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y_L.push_back(&(Q_L[si].get()));
            rho_Y_R.push_back(&(Q_R[si].get()));
        }
        
        const double rho_L = d_equation_of_state_mixing_rules->getMixtureDensity(
            rho_Y_L);
        
        const double rho_R = d_equation_of_state_mixing_rules->getMixtureDensity(
            rho_Y_R);
        
        const double u_L = (Q_L[d_num_species].get())/rho_L;
        const double u_R = (Q_R[d_num_species].get())/rho_R;
        
        const double v_L = (Q_L[d_num_species + 1].get())/rho_L;
        const double v_R = (Q_R[d_num_species + 1].get())/rho_R;
        
        std::vector<const double*> m_L;
        std::vector<const double*> m_R;
        m_L.reserve(2);
        m_R.reserve(2);
        m_L.push_back(&(Q_L[d_num_species].get()));
        m_R.push_back(&(Q_R[d_num_species].get()));
        m_L.push_back(&(Q_L[d_num_species + 1].get()));
        m_R.push_back(&(Q_R[d_num_species + 1].get()));
        
        std::vector<const double*> vel_L;
        std::vector<const double*> vel_R;
        vel_L.reserve(2);
        vel_R.reserve(2);
        vel_L.push_back(&u_L);
        vel_R.push_back(&u_R);
        vel_L.push_back(&v_L);
        vel_R.push_back(&v_R);
        
        std::vector<const double*> empty_ptr;
        
        const double p_L = d_equation_of_state_mixing_rules->getPressure(
            rho_Y_L,
            m_L,
            &(Q_L[d_num_species + d_dim.getValue()].get()),
            empty_ptr);
        
        const double p_R = d_equation_of_state_mixing_rules->getPressure(
            rho_Y_R,
            m_R,
            &(Q_R[d_num_species + d_dim.getValue()].get()),
            empty_ptr);
        
        const double c_L = d_equation_of_state_mixing_rules->getSoundSpeed(
            rho_Y_L,
            vel_L,
            &p_L,
            empty_ptr);
        
        const double c_R = d_equation_of_state_mixing_rules->getSoundSpeed(
            rho_Y_R,
            vel_R,
            &p_R,
            empty_ptr);
        
        const double u_average = 0.5*(u_L + u_R);
        const double c_average = 0.5*(c_L + c_R);
        
        const double s_L = fmin(u_average - c_average, u_L - c_L);
        const double s_R = fmax(u_average + c_average, u_R + c_R);
        
        const double s_minus = fmin(0.0, s_L);
        const double s_plus  = fmax(0.0, s_R);
        
        const double s_star =
            (p_R - p_L + (Q_L[d_num_species].get())*(s_L - u_L) - (Q_R[d_num_species].get())*(s_R - u_R))/
                (rho_L*(s_L - u_L) - rho_R*(s_R - u_R));
        
        if (s_star > 0)
        {
            const double Chi_star_L = (s_L - u_L)/(s_L - s_star);
            
            double Q_star_L[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                Q_star_L[si] = Chi_star_L*(Q_L[si].get());
            }
            Q_star_L[d_num_species] = Chi_star_L*rho_L*s_star;
            Q_star_L[d_num_species + 1] = Chi_star_L*(Q_L[d_num_species + 1].get());
            Q_star_L[d_num_species + 2] = Chi_star_L*((Q_L[d_num_species + 2].get()) +
                (s_star - u_L)*(rho_L*s_star + p_L/(s_L - u_L)));
            
            double F_x_L[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                F_x_L[si] = u_L*(Q_L[si].get());
            }
            F_x_L[d_num_species] = u_L*(Q_L[d_num_species].get()) + p_L;
            F_x_L[d_num_species + 1] = u_L*(Q_L[d_num_species + 1].get());
            F_x_L[d_num_species + 2] = u_L*((Q_L[d_num_species + 2].get()) + p_L);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_x_intercell[ei].get() = F_x_L[ei] + s_minus*(Q_star_L[ei] - Q_L[ei]);
            }
        }
        else
        {
            const double Chi_star_R = (s_R - u_R)/(s_R - s_star);
            
            double Q_star_R[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                Q_star_R[si] = Chi_star_R*(Q_R[si].get());
            }
            Q_star_R[d_num_species] = Chi_star_R*rho_R*s_star;
            Q_star_R[d_num_species + 1] = Chi_star_R*(Q_R[d_num_species + 1].get());
            Q_star_R[d_num_species + 2] = Chi_star_R*((Q_R[d_num_species + 2].get()) +
                (s_star - u_R)*(rho_R*s_star + p_R/(s_R - u_R)));
            
            double F_x_R[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                F_x_R[si] = u_R*(Q_R[si].get());
            }
            F_x_R[d_num_species] = u_R*(Q_R[d_num_species].get()) + p_R;
            F_x_R[d_num_species + 1] = u_R*(Q_R[d_num_species + 1].get());
            F_x_R[d_num_species + 2] = u_R*((Q_R[d_num_species + 2].get()) + p_R);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_x_intercell[ei].get() = F_x_R[ei] + s_plus*(Q_star_R[ei] - Q_R[ei]);
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        std::vector<const double*> rho_Y_L;
        std::vector<const double*> rho_Y_R;
        rho_Y_L.reserve(d_num_species);
        rho_Y_R.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y_L.push_back(&(Q_L[si].get()));
            rho_Y_R.push_back(&(Q_R[si].get()));
        }
        
        const double rho_L = d_equation_of_state_mixing_rules->getMixtureDensity(
            rho_Y_L);
        
        const double rho_R = d_equation_of_state_mixing_rules->getMixtureDensity(
            rho_Y_R);
        
        const double u_L = (Q_L[d_num_species].get())/rho_L;
        const double u_R = (Q_R[d_num_species].get())/rho_R;
        
        const double v_L = (Q_L[d_num_species + 1].get())/rho_L;
        const double v_R = (Q_R[d_num_species + 1].get())/rho_R;
        
        const double w_L = (Q_L[d_num_species + 2].get())/rho_L;
        const double w_R = (Q_R[d_num_species + 2].get())/rho_R;
        
        std::vector<const double*> m_L;
        std::vector<const double*> m_R;
        m_L.reserve(3);
        m_R.reserve(3);
        m_L.push_back(&(Q_L[d_num_species].get()));
        m_R.push_back(&(Q_R[d_num_species].get()));
        m_L.push_back(&(Q_L[d_num_species + 1].get()));
        m_R.push_back(&(Q_R[d_num_species + 1].get()));
        m_L.push_back(&(Q_L[d_num_species + 2].get()));
        m_R.push_back(&(Q_R[d_num_species + 2].get()));
        
        std::vector<const double*> vel_L;
        std::vector<const double*> vel_R;
        vel_L.reserve(3);
        vel_R.reserve(3);
        vel_L.push_back(&u_L);
        vel_R.push_back(&u_R);
        vel_L.push_back(&v_L);
        vel_R.push_back(&v_R);
        vel_L.push_back(&w_L);
        vel_R.push_back(&w_R);
        
        std::vector<const double*> empty_ptr;
        
        const double p_L = d_equation_of_state_mixing_rules->getPressure(
            rho_Y_L,
            m_L,
            &(Q_L[d_num_species + d_dim.getValue()].get()),
            empty_ptr);
        
        const double p_R = d_equation_of_state_mixing_rules->getPressure(
            rho_Y_R,
            m_R,
            &(Q_R[d_num_species + d_dim.getValue()].get()),
            empty_ptr);
        
        const double c_L = d_equation_of_state_mixing_rules->getSoundSpeed(
            rho_Y_L,
            vel_L,
            &p_L,
            empty_ptr);
        
        const double c_R = d_equation_of_state_mixing_rules->getSoundSpeed(
            rho_Y_R,
            vel_R,
            &p_R,
            empty_ptr);
        
        const double u_average = 0.5*(u_L + u_R);
        const double c_average = 0.5*(c_L + c_R);
        
        const double s_L = fmin(u_average - c_average, u_L - c_L);
        const double s_R = fmax(u_average + c_average, u_R + c_R);
        
        const double s_minus = fmin(0.0, s_L);
        const double s_plus  = fmax(0.0, s_R);
        
        const double s_star =
            (p_R - p_L + (Q_L[d_num_species].get())*(s_L - u_L) - (Q_R[d_num_species].get())*(s_R - u_R))/
                (rho_L*(s_L - u_L) - rho_R*(s_R - u_R));
        
        if (s_star > 0)
        {
            const double Chi_star_L = (s_L - u_L)/(s_L - s_star);
            
            double Q_star_L[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                Q_star_L[si] = Chi_star_L*(Q_L[si].get());
            }
            Q_star_L[d_num_species] = Chi_star_L*rho_L*s_star;
            Q_star_L[d_num_species + 1] = Chi_star_L*(Q_L[d_num_species + 1].get());
            Q_star_L[d_num_species + 2] = Chi_star_L*(Q_L[d_num_species + 2].get());
            Q_star_L[d_num_species + 3] = Chi_star_L*((Q_L[d_num_species + 3].get()) +
                (s_star - u_L)*(rho_L*s_star + p_L/(s_L - u_L)));
            
            double F_x_L[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                F_x_L[si] = u_L*(Q_L[si].get());
            }
            F_x_L[d_num_species] = u_L*(Q_L[d_num_species].get()) + p_L;
            F_x_L[d_num_species + 1] = u_L*(Q_L[d_num_species + 1].get());
            F_x_L[d_num_species + 2] = u_L*(Q_L[d_num_species + 2].get());
            F_x_L[d_num_species + 3] = u_L*((Q_L[d_num_species + 3].get()) + p_L);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_x_intercell[ei].get() = F_x_L[ei] + s_minus*(Q_star_L[ei] - Q_L[ei]);
            }
        }
        else
        {
            const double Chi_star_R = (s_R - u_R)/(s_R - s_star);
            
            double Q_star_R[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                Q_star_R[si] = Chi_star_R*(Q_R[si].get());
            }
            Q_star_R[d_num_species] = Chi_star_R*rho_R*s_star;
            Q_star_R[d_num_species + 1] = Chi_star_R*(Q_R[d_num_species + 1].get());
            Q_star_R[d_num_species + 2] = Chi_star_R*(Q_R[d_num_species + 2].get());
            Q_star_R[d_num_species + 3] = Chi_star_R*((Q_R[d_num_species + 3].get()) +
                (s_star - u_R)*(rho_R*s_star + p_R/(s_R - u_R)));
            
            double F_x_R[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                F_x_R[si] = u_R*(Q_R[si].get());
            }
            F_x_R[d_num_species] = u_R*(Q_R[d_num_species].get()) + p_R;
            F_x_R[d_num_species + 1] = u_R*(Q_R[d_num_species + 1].get());
            F_x_R[d_num_species + 2] = u_R*(Q_R[d_num_species + 2].get());
            F_x_R[d_num_species + 3] = u_R*((Q_R[d_num_species + 3].get()) + p_R);
            
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
RiemannSolverFourEqnConservativeHLLC::computeIntercellFluxInYDirectionFromConservativeVariables(
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
            << ": RiemannSolverFourEqnConservativeHLLC::"
            << "computeIntercellFluxInYDirectionFromConservativeVariables()\n"
            << "There is no y direction for 1D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        std::vector<const double*> rho_Y_B;
        std::vector<const double*> rho_Y_T;
        rho_Y_B.reserve(d_num_species);
        rho_Y_T.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y_B.push_back(&(Q_B[si].get()));
            rho_Y_T.push_back(&(Q_T[si].get()));
        }
        
        const double rho_B = d_equation_of_state_mixing_rules->getMixtureDensity(
            rho_Y_B);
        
        const double rho_T = d_equation_of_state_mixing_rules->getMixtureDensity(
            rho_Y_T);
        
        const double u_B = (Q_B[d_num_species].get())/rho_B;
        const double u_T = (Q_T[d_num_species].get())/rho_T;
        
        const double v_B = (Q_B[d_num_species + 1].get())/rho_B;
        const double v_T = (Q_T[d_num_species + 1].get())/rho_T;
        
        std::vector<const double*> m_B;
        std::vector<const double*> m_T;
        m_B.reserve(2);
        m_T.reserve(2);
        m_B.push_back(&(Q_B[d_num_species].get()));
        m_T.push_back(&(Q_T[d_num_species].get()));
        m_B.push_back(&(Q_B[d_num_species + 1].get()));
        m_T.push_back(&(Q_T[d_num_species + 1].get()));
        
        std::vector<const double*> vel_B;
        std::vector<const double*> vel_T;
        vel_B.reserve(2);
        vel_T.reserve(2);
        vel_B.push_back(&u_B);
        vel_T.push_back(&u_T);
        vel_B.push_back(&v_B);
        vel_T.push_back(&v_T);
        
        std::vector<const double*> empty_ptr;
        
        const double p_B = d_equation_of_state_mixing_rules->getPressure(
            rho_Y_B,
            m_B,
            &(Q_B[d_num_species + d_dim.getValue()].get()),
            empty_ptr);
        
        const double p_T = d_equation_of_state_mixing_rules->getPressure(
            rho_Y_T,
            m_T,
            &(Q_T[d_num_species + d_dim.getValue()].get()),
            empty_ptr);
        
        const double c_B = d_equation_of_state_mixing_rules->getSoundSpeed(
            rho_Y_B,
            vel_B,
            &p_B,
            empty_ptr);
        
        const double c_T = d_equation_of_state_mixing_rules->getSoundSpeed(
            rho_Y_T,
            vel_T,
            &p_T,
            empty_ptr);
        
        const double v_average = 0.5*(v_B + v_T);
        const double c_average = 0.5*(c_B + c_T);
        
        const double s_B = fmin(v_average - c_average, v_B - c_B);
        const double s_T = fmax(v_average + c_average, v_T + c_T);
        
        const double s_minus = fmin(0.0, s_B);
        const double s_plus  = fmax(0.0, s_T);
        
        const double s_star =
            (p_T - p_B + (Q_B[d_num_species + 1].get())*(s_B - v_B) - (Q_T[d_num_species + 1].get())*(s_T - v_T))/
                (rho_B*(s_B - v_B) - rho_T*(s_T - v_T));
        
        if (s_star > 0)
        {
            const double Chi_star_B = (s_B - v_B)/(s_B - s_star);
            
            double Q_star_B[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                Q_star_B[si] = Chi_star_B*(Q_B[si].get());
            }
            Q_star_B[d_num_species] = Chi_star_B*(Q_B[d_num_species].get());
            Q_star_B[d_num_species + 1] = Chi_star_B*rho_B*s_star;
            Q_star_B[d_num_species + 2] = Chi_star_B*((Q_B[d_num_species + 2].get()) +
                (s_star - v_B)*(rho_B*s_star + p_B/(s_B - v_B)));
            
            double F_y_B[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                F_y_B[si] = v_B*(Q_B[si].get());
            }
            F_y_B[d_num_species] = v_B*(Q_B[d_num_species].get());
            F_y_B[d_num_species + 1] = v_B*(Q_B[d_num_species + 1].get()) + p_B;
            F_y_B[d_num_species + 2] = v_B*((Q_B[d_num_species + 2].get()) + p_B);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_y_intercell[ei].get() = F_y_B[ei] + s_minus*(Q_star_B[ei] - Q_B[ei]);
            }
        }
        else
        {
            const double Chi_star_T = (s_T - v_T)/(s_T - s_star);
            
            double Q_star_T[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                Q_star_T[si] = Chi_star_T*(Q_T[si].get());
            }
            Q_star_T[d_num_species] = Chi_star_T*(Q_T[d_num_species].get());
            Q_star_T[d_num_species + 1] = Chi_star_T*rho_T*s_star;
            Q_star_T[d_num_species + 2] = Chi_star_T*((Q_T[d_num_species + 2].get()) +
                (s_star - v_T)*(rho_T*s_star + p_T/(s_T - v_T)));
            
            double F_y_T[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                F_y_T[si] = v_T*(Q_T[si].get());
            }
            F_y_T[d_num_species] = v_T*(Q_T[d_num_species].get());
            F_y_T[d_num_species + 1] = v_T*(Q_T[d_num_species + 1].get()) + p_T;
            F_y_T[d_num_species + 2] = v_T*((Q_T[d_num_species + 2].get()) + p_T);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_y_intercell[ei].get() = F_y_T[ei] + s_plus*(Q_star_T[ei] - Q_T[ei]);
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        std::vector<const double*> rho_Y_B;
        std::vector<const double*> rho_Y_T;
        rho_Y_B.reserve(d_num_species);
        rho_Y_T.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y_B.push_back(&(Q_B[si].get()));
            rho_Y_T.push_back(&(Q_T[si].get()));
        }
        
        const double rho_B = d_equation_of_state_mixing_rules->getMixtureDensity(
            rho_Y_B);
        
        const double rho_T = d_equation_of_state_mixing_rules->getMixtureDensity(
            rho_Y_T);
        
        const double u_B = (Q_B[d_num_species].get())/rho_B;
        const double u_T = (Q_T[d_num_species].get())/rho_T;
        
        const double v_B = (Q_B[d_num_species + 1].get())/rho_B;
        const double v_T = (Q_T[d_num_species + 1].get())/rho_T;
        
        const double w_B = (Q_B[d_num_species + 2].get())/rho_B;
        const double w_T = (Q_T[d_num_species + 2].get())/rho_T;
        
        std::vector<const double*> m_B;
        std::vector<const double*> m_T;
        m_B.reserve(3);
        m_T.reserve(3);
        m_B.push_back(&(Q_B[d_num_species].get()));
        m_T.push_back(&(Q_T[d_num_species].get()));
        m_B.push_back(&(Q_B[d_num_species + 1].get()));
        m_T.push_back(&(Q_T[d_num_species + 1].get()));
        m_B.push_back(&(Q_B[d_num_species + 2].get()));
        m_T.push_back(&(Q_T[d_num_species + 2].get()));
        
        std::vector<const double*> vel_B;
        std::vector<const double*> vel_T;
        vel_B.reserve(3);
        vel_T.reserve(3);
        vel_B.push_back(&u_B);
        vel_T.push_back(&u_T);
        vel_B.push_back(&v_B);
        vel_T.push_back(&v_T);
        vel_B.push_back(&w_B);
        vel_T.push_back(&w_T);
        
        std::vector<const double*> empty_ptr;
        
        const double p_B = d_equation_of_state_mixing_rules->getPressure(
            rho_Y_B,
            m_B,
            &(Q_B[d_num_species + d_dim.getValue()].get()),
            empty_ptr);
        
        const double p_T = d_equation_of_state_mixing_rules->getPressure(
            rho_Y_T,
            m_T,
            &(Q_T[d_num_species + d_dim.getValue()].get()),
            empty_ptr);
        
        const double c_B = d_equation_of_state_mixing_rules->getSoundSpeed(
            rho_Y_B,
            vel_B,
            &p_B,
            empty_ptr);
        
        const double c_T = d_equation_of_state_mixing_rules->getSoundSpeed(
            rho_Y_T,
            vel_T,
            &p_T,
            empty_ptr);
        
        const double v_average = 0.5*(v_B + v_T);
        const double c_average = 0.5*(c_B + c_T);
        
        const double s_B = fmin(v_average - c_average, v_B - c_B);
        const double s_T = fmax(v_average + c_average, v_T + c_T);
        
        const double s_minus = fmin(0.0, s_B);
        const double s_plus  = fmax(0.0, s_T);
        
        const double s_star =
            (p_T - p_B + (Q_B[d_num_species + 1].get())*(s_B - v_B) - (Q_T[d_num_species + 1].get())*(s_T - v_T))/
                (rho_B*(s_B - v_B) - rho_T*(s_T - v_T));
        
        if (s_star > 0)
        {
            const double Chi_star_B = (s_B - v_B)/(s_B - s_star);
            
            double Q_star_B[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                Q_star_B[si] = Chi_star_B*(Q_B[si].get());
            }
            Q_star_B[d_num_species] = Chi_star_B*(Q_B[d_num_species].get());
            Q_star_B[d_num_species + 1] = Chi_star_B*rho_B*s_star;
            Q_star_B[d_num_species + 2] = Chi_star_B*(Q_B[d_num_species + 2].get());
            Q_star_B[d_num_species + 3] = Chi_star_B*((Q_B[d_num_species + 3].get()) +
                (s_star - v_B)*(rho_B*s_star + p_B/(s_B - v_B)));
            
            double F_y_B[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                F_y_B[si] = v_B*(Q_B[si].get());
            }
            F_y_B[d_num_species] = v_B*(Q_B[d_num_species].get());
            F_y_B[d_num_species + 1] = v_B*(Q_B[d_num_species + 1].get()) + p_B;
            F_y_B[d_num_species + 2] = v_B*(Q_B[d_num_species + 2].get());
            F_y_B[d_num_species + 3] = v_B*((Q_B[d_num_species + 3].get()) + p_B);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_y_intercell[ei].get() = F_y_B[ei] + s_minus*(Q_star_B[ei] - Q_B[ei]);
            }
        }
        else
        {
            const double Chi_star_T = (s_T - v_T)/(s_T - s_star);
            
            double Q_star_T[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                Q_star_T[si] = Chi_star_T*(Q_T[si].get());
            }
            Q_star_T[d_num_species] = Chi_star_T*(Q_T[d_num_species].get());
            Q_star_T[d_num_species + 1] = Chi_star_T*rho_T*s_star;
            Q_star_T[d_num_species + 2] = Chi_star_T*(Q_T[d_num_species + 2].get());
            Q_star_T[d_num_species + 3] = Chi_star_T*((Q_T[d_num_species + 3].get()) +
                (s_star - v_T)*(rho_T*s_star + p_T/(s_T - v_T)));
            
            double F_y_T[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                F_y_T[si] = v_T*(Q_T[si].get());
            }
            F_y_T[d_num_species] = v_T*(Q_T[d_num_species].get());
            F_y_T[d_num_species + 1] = v_T*(Q_T[d_num_species + 1].get()) + p_T;
            F_y_T[d_num_species + 2] = v_T*(Q_T[d_num_species + 2].get());
            F_y_T[d_num_species + 3] = v_T*((Q_T[d_num_species + 3].get()) + p_T);
            
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
RiemannSolverFourEqnConservativeHLLC::computeIntercellFluxInZDirectionFromConservativeVariables(
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
            << ": RiemannSolverFourEqnConservativeHLLC::"
            << "computeIntercellFluxInZDirectionFromConservativeVariables()\n"
            << "There is no z direction for 1D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": RiemannSolverFourEqnConservativeHLLC::"
            << "computeIntercellFluxInZDirectionFromConservativeVariables()\n"
            << "There is no z direction for 2D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        std::vector<const double*> rho_Y_B;
        std::vector<const double*> rho_Y_F;
        rho_Y_B.reserve(d_num_species);
        rho_Y_F.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y_B.push_back(&(Q_B[si].get()));
            rho_Y_F.push_back(&(Q_F[si].get()));
        }
        
        const double rho_B = d_equation_of_state_mixing_rules->getMixtureDensity(
            rho_Y_B);
        
        const double rho_F = d_equation_of_state_mixing_rules->getMixtureDensity(
            rho_Y_F);
        
        const double u_B = (Q_B[d_num_species].get())/rho_B;
        const double u_F = (Q_F[d_num_species].get())/rho_F;
        
        const double v_B = (Q_B[d_num_species + 1].get())/rho_B;
        const double v_F = (Q_F[d_num_species + 1].get())/rho_F;
        
        const double w_B = (Q_B[d_num_species + 2].get())/rho_B;
        const double w_F = (Q_F[d_num_species + 2].get())/rho_F;
        
        std::vector<const double*> m_B;
        std::vector<const double*> m_F;
        m_B.reserve(3);
        m_F.reserve(3);
        m_B.push_back(&(Q_B[d_num_species].get()));
        m_F.push_back(&(Q_F[d_num_species].get()));
        m_B.push_back(&(Q_B[d_num_species + 1].get()));
        m_F.push_back(&(Q_F[d_num_species + 1].get()));
        m_B.push_back(&(Q_B[d_num_species + 2].get()));
        m_F.push_back(&(Q_F[d_num_species + 2].get()));
        
        std::vector<const double*> vel_B;
        std::vector<const double*> vel_F;
        vel_B.reserve(3);
        vel_F.reserve(3);
        vel_B.push_back(&u_B);
        vel_F.push_back(&u_F);
        vel_B.push_back(&v_B);
        vel_F.push_back(&v_F);
        vel_B.push_back(&w_B);
        vel_F.push_back(&w_F);
        
        std::vector<const double*> empty_ptr;
        
        const double p_B = d_equation_of_state_mixing_rules->getPressure(
            rho_Y_B,
            m_B,
            &(Q_B[d_num_species + d_dim.getValue()].get()),
            empty_ptr);
        
        const double p_F = d_equation_of_state_mixing_rules->getPressure(
            rho_Y_F,
            m_F,
            &(Q_F[d_num_species + d_dim.getValue()].get()),
            empty_ptr);
        
        const double c_B = d_equation_of_state_mixing_rules->getSoundSpeed(
            rho_Y_B,
            vel_B,
            &p_B,
            empty_ptr);
        
        const double c_F = d_equation_of_state_mixing_rules->getSoundSpeed(
            rho_Y_F,
            vel_F,
            &p_F,
            empty_ptr);
        
        const double w_average = 0.5*(w_B + w_F);
        const double c_average = 0.5*(c_B + c_F);
        
        const double s_B = fmin(w_average - c_average, w_B - c_B);
        const double s_F = fmax(w_average + c_average, w_F + c_F);
        
        const double s_minus = fmin(0.0, s_B);
        const double s_plus  = fmax(0.0, s_F);
        
        const double s_star =
            (p_F - p_B + (Q_B[d_num_species + 2].get())*(s_B - w_B) - (Q_F[d_num_species + 2].get())*(s_F - w_F))/
                (rho_B*(s_B - w_B) - rho_F*(s_F - w_F));
        
        if (s_star > 0)
        {
            const double Chi_star_B = (s_B - w_B)/(s_B - s_star);
            
            double Q_star_B[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                Q_star_B[si] = Chi_star_B*(Q_B[si].get());
            }
            Q_star_B[d_num_species] = Chi_star_B*(Q_B[d_num_species].get());
            Q_star_B[d_num_species + 1] = Chi_star_B*(Q_B[d_num_species + 1].get());
            Q_star_B[d_num_species + 2] = Chi_star_B*rho_B*s_star;
            Q_star_B[d_num_species + 3] = Chi_star_B*((Q_B[d_num_species + 3].get()) +
                (s_star - w_B)*(rho_B*s_star + p_B/(s_B - w_B)));
            
            double F_z_B[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                F_z_B[si] = w_B*(Q_B[si].get());
            }
            F_z_B[d_num_species] = w_B*(Q_B[d_num_species].get());
            F_z_B[d_num_species + 1] = w_B*(Q_B[d_num_species + 1].get());
            F_z_B[d_num_species + 2] = w_B*(Q_B[d_num_species + 2].get()) + p_B;
            F_z_B[d_num_species + 3] = w_B*((Q_B[d_num_species + 3].get()) + p_B);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_z_intercell[ei].get() = F_z_B[ei] + s_minus*(Q_star_B[ei] - Q_B[ei]);
            }
        }
        else
        {
            const double Chi_star_F = (s_F - w_F)/(s_F - s_star);
            
            double Q_star_F[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                Q_star_F[si] = Chi_star_F*(Q_F[si].get());
            }
            Q_star_F[d_num_species] = Chi_star_F*(Q_F[d_num_species].get());
            Q_star_F[d_num_species + 1] = Chi_star_F*(Q_F[d_num_species + 1].get());
            Q_star_F[d_num_species + 2] = Chi_star_F*rho_F*s_star;
            Q_star_F[d_num_species + 3] = Chi_star_F*((Q_F[d_num_species + 3].get()) +
                (s_star - w_F)*(rho_F*s_star + p_F/(s_F - w_F)));
            
            double F_z_F[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                F_z_F[si] = w_F*(Q_F[si].get());
            }
            F_z_F[d_num_species] = w_F*(Q_F[d_num_species].get());
            F_z_F[d_num_species + 1] = w_F*(Q_F[d_num_species + 1].get());
            F_z_F[d_num_species + 2] = w_F*(Q_F[d_num_species + 2].get()) + p_F;
            F_z_F[d_num_species + 3] = w_F*((Q_F[d_num_species + 3].get()) + p_F);
            
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
RiemannSolverFourEqnConservativeHLLC::computeIntercellFluxInXDirectionFromPrimitiveVariables(
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
        std::vector<const double*> rho_Y_L;
        std::vector<const double*> rho_Y_R;
        rho_Y_L.reserve(d_num_species);
        rho_Y_R.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y_L.push_back(&(V_L[si].get()));
            rho_Y_R.push_back(&(V_R[si].get()));
        }
        
        std::vector<const double*> vel_L;
        std::vector<const double*> vel_R;
        vel_L.reserve(1);
        vel_R.reserve(1);
        vel_L.push_back(&(V_L[d_num_species].get()));
        vel_R.push_back(&(V_R[d_num_species].get()));
        
        const double rho_L = d_equation_of_state_mixing_rules->getMixtureDensity(
            rho_Y_L);
        
        const double rho_R = d_equation_of_state_mixing_rules->getMixtureDensity(
            rho_Y_R);
        
        std::vector<const double*> empty_ptr;
        
        const double c_L = d_equation_of_state_mixing_rules->getSoundSpeed(
            rho_Y_L,
            vel_L,
            &(V_L[d_num_species + d_dim.getValue()].get()),
            empty_ptr);
        
        const double c_R = d_equation_of_state_mixing_rules->getSoundSpeed(
            rho_Y_R,
            vel_R,
            &(V_R[d_num_species + d_dim.getValue()].get()),
            empty_ptr);
        
        const double u_average = 0.5*((V_L[d_num_species].get()) + (V_R[d_num_species].get()));
        const double c_average = 0.5*(c_L + c_R);
        
        const double s_L = fmin(u_average - c_average, (V_L[d_num_species].get()) - c_L);
        const double s_R = fmax(u_average + c_average, (V_R[d_num_species].get()) + c_R);
        
        const double s_minus = fmin(0.0, s_L);
        const double s_plus  = fmax(0.0, s_R);
        
        const double s_star =
            ((V_R[d_num_species + d_dim.getValue()].get()) - (V_L[d_num_species + d_dim.getValue()].get()) +
                rho_L*(V_L[d_num_species].get())*(s_L - (V_L[d_num_species].get())) -
                    rho_R*(V_R[d_num_species].get())*(s_R - (V_R[d_num_species].get())))/
                        (rho_L*(s_L - (V_L[d_num_species].get())) - rho_R*(s_R - (V_R[d_num_species].get())));
        
        if (s_star > 0)
        {
            const double Chi_star_L = (s_L - (V_L[d_num_species].get()))/(s_L - s_star);
            
            std::vector<const double*> vel_L;
            vel_L.reserve(1);
            vel_L.push_back(&(V_L[d_num_species].get()));
            
            double Q_L[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                Q_L[si] = (V_L[si].get());
            }
            Q_L[d_num_species] = rho_L*(V_L[d_num_species].get());
            Q_L[d_num_species + d_dim.getValue()] = d_equation_of_state_mixing_rules->getTotalEnergy(
                rho_Y_L,
                vel_L,
                &(V_L[d_num_species + d_dim.getValue()].get()),
                empty_ptr);
            
            double Q_star_L[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                Q_star_L[si] = Chi_star_L*(V_L[si].get());
            }
            Q_star_L[d_num_species] = Chi_star_L*rho_L*s_star;
            Q_star_L[d_num_species + d_dim.getValue()] = Chi_star_L*(Q_L[d_num_species + d_dim.getValue()] +
                (s_star - (V_L[d_num_species].get()))*(rho_L*s_star +
                        (V_L[d_num_species + d_dim.getValue()].get())/(s_L - (V_L[d_num_species].get()))));
            
            double F_x_L[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                F_x_L[si] = (V_L[d_num_species].get())*(V_L[si].get());
            }
            F_x_L[d_num_species] = (V_L[d_num_species].get())*Q_L[d_num_species] +
                (V_L[d_num_species + d_dim.getValue()].get());
            F_x_L[d_num_species + d_dim.getValue()] = (V_L[d_num_species].get())*(Q_L[d_num_species + d_dim.getValue()] +
                (V_L[d_num_species + d_dim.getValue()].get()));
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_x_intercell[ei].get() = F_x_L[ei] + s_minus*(Q_star_L[ei] - Q_L[ei]);
            }
        }
        else
        {
            const double Chi_star_R = (s_R - (V_R[d_num_species].get()))/(s_R - s_star);
            
            std::vector<const double*> vel_R;
            vel_R.reserve(1);
            vel_R.push_back(&(V_R[d_num_species].get()));
            
            double Q_R[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                Q_R[si] = (V_R[si].get());
            }
            Q_R[d_num_species] = rho_R*(V_R[d_num_species].get());
            Q_R[d_num_species + d_dim.getValue()] = d_equation_of_state_mixing_rules->getTotalEnergy(
                rho_Y_R,
                vel_R,
                &(V_R[d_num_species + d_dim.getValue()].get()),
                empty_ptr);
            
            double Q_star_R[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                Q_star_R[si] = Chi_star_R*(V_R[si].get());
            }
            Q_star_R[d_num_species] = Chi_star_R*rho_R*s_star;
            Q_star_R[d_num_species + d_dim.getValue()] = Chi_star_R*(Q_R[d_num_species + d_dim.getValue()] +
                (s_star - (V_R[d_num_species].get()))*(rho_R*s_star +
                        (V_R[d_num_species + d_dim.getValue()].get())/(s_R - (V_R[d_num_species].get()))));
            
            double F_x_R[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                F_x_R[si] = (V_R[d_num_species].get())*(V_R[si].get());
            }
            F_x_R[d_num_species] = (V_R[d_num_species].get())*Q_R[d_num_species] +
                (V_R[d_num_species + d_dim.getValue()].get());
            F_x_R[d_num_species + d_dim.getValue()] = (V_R[d_num_species].get())*(Q_R[d_num_species + d_dim.getValue()] +
                (V_R[d_num_species + d_dim.getValue()].get()));
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_x_intercell[ei].get() = F_x_R[ei] + s_plus*(Q_star_R[ei] - Q_R[ei]);
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        std::vector<const double*> rho_Y_L;
        std::vector<const double*> rho_Y_R;
        rho_Y_L.reserve(d_num_species);
        rho_Y_R.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y_L.push_back(&(V_L[si].get()));
            rho_Y_R.push_back(&(V_R[si].get()));
        }
        
        std::vector<const double*> vel_L;
        std::vector<const double*> vel_R;
        vel_L.reserve(2);
        vel_R.reserve(2);
        vel_L.push_back(&(V_L[d_num_species].get()));
        vel_R.push_back(&(V_R[d_num_species].get()));
        vel_L.push_back(&(V_L[d_num_species + 1].get()));
        vel_R.push_back(&(V_R[d_num_species + 1].get()));
        
        const double rho_L = d_equation_of_state_mixing_rules->getMixtureDensity(
            rho_Y_L);
        
        const double rho_R = d_equation_of_state_mixing_rules->getMixtureDensity(
            rho_Y_R);
        
        std::vector<const double*> empty_ptr;
        
        const double c_L = d_equation_of_state_mixing_rules->getSoundSpeed(
            rho_Y_L,
            vel_L,
            &(V_L[d_num_species + d_dim.getValue()].get()),
            empty_ptr);
        
        const double c_R = d_equation_of_state_mixing_rules->getSoundSpeed(
            rho_Y_R,
            vel_R,
            &(V_R[d_num_species + d_dim.getValue()].get()),
            empty_ptr);
        
        const double u_average = 0.5*((V_L[d_num_species].get()) + (V_R[d_num_species].get()));
        const double c_average = 0.5*(c_L + c_R);
        
        const double s_L = fmin(u_average - c_average, (V_L[d_num_species].get()) - c_L);
        const double s_R = fmax(u_average + c_average, (V_R[d_num_species].get()) + c_R);
        
        const double s_minus = fmin(0.0, s_L);
        const double s_plus  = fmax(0.0, s_R);
        
        const double s_star =
            ((V_R[d_num_species + d_dim.getValue()].get()) - (V_L[d_num_species + d_dim.getValue()].get()) +
                rho_L*(V_L[d_num_species].get())*(s_L - (V_L[d_num_species].get())) -
                    rho_R*(V_R[d_num_species].get())*(s_R - (V_R[d_num_species].get())))/
                        (rho_L*(s_L - (V_L[d_num_species].get())) - rho_R*(s_R - (V_R[d_num_species].get())));
        
        if (s_star > 0)
        {
            const double Chi_star_L = (s_L - (V_L[d_num_species].get()))/(s_L - s_star);
            
            std::vector<const double*> vel_L;
            vel_L.reserve(2);
            vel_L.push_back(&(V_L[d_num_species].get()));
            vel_L.push_back(&(V_L[d_num_species + 1].get()));
            
            double Q_L[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                Q_L[si] = (V_L[si].get());
            }
            Q_L[d_num_species] = rho_L*(V_L[d_num_species].get());
            Q_L[d_num_species + 1] = rho_L*(V_L[d_num_species + 1].get());
            Q_L[d_num_species + d_dim.getValue()] = d_equation_of_state_mixing_rules->getTotalEnergy(
                rho_Y_L,
                vel_L,
                &(V_L[d_num_species + d_dim.getValue()].get()),
                empty_ptr);
            
            double Q_star_L[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                Q_star_L[si] = Chi_star_L*(V_L[si].get());
            }
            Q_star_L[d_num_species] = Chi_star_L*rho_L*s_star;
            Q_star_L[d_num_species + 1] = Chi_star_L*Q_L[d_num_species + 1];
            Q_star_L[d_num_species + d_dim.getValue()] = Chi_star_L*(Q_L[d_num_species + d_dim.getValue()] +
                (s_star - (V_L[d_num_species].get()))*(rho_L*s_star +
                        (V_L[d_num_species + d_dim.getValue()].get())/(s_L - (V_L[d_num_species].get()))));
            
            double F_x_L[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                F_x_L[si] = (V_L[d_num_species].get())*(V_L[si].get());
            }
            F_x_L[d_num_species] = (V_L[d_num_species].get())*Q_L[d_num_species] +
                (V_L[d_num_species + d_dim.getValue()].get());
            F_x_L[d_num_species + 1] = (V_L[d_num_species].get())*Q_L[d_num_species + 1];
            F_x_L[d_num_species + d_dim.getValue()] = (V_L[d_num_species].get())*(Q_L[d_num_species + d_dim.getValue()] +
                (V_L[d_num_species + d_dim.getValue()].get()));
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_x_intercell[ei].get() = F_x_L[ei] + s_minus*(Q_star_L[ei] - Q_L[ei]);
            }
        }
        else
        {
            const double Chi_star_R = (s_R - (V_R[d_num_species].get()))/(s_R - s_star);
            
            std::vector<const double*> vel_R;
            vel_R.reserve(2);
            vel_R.push_back(&(V_R[d_num_species].get()));
            vel_R.push_back(&(V_R[d_num_species + 1].get()));
            
            double Q_R[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                Q_R[si] = (V_R[si].get());
            }
            Q_R[d_num_species] = rho_R*(V_R[d_num_species].get());
            Q_R[d_num_species + 1] = rho_R*(V_R[d_num_species + 1].get());
            Q_R[d_num_species + d_dim.getValue()] = d_equation_of_state_mixing_rules->getTotalEnergy(
                rho_Y_R,
                vel_R,
                &(V_R[d_num_species + d_dim.getValue()].get()),
                empty_ptr);
            
            double Q_star_R[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                Q_star_R[si] = Chi_star_R*(V_R[si].get());
            }
            Q_star_R[d_num_species] = Chi_star_R*rho_R*s_star;
            Q_star_R[d_num_species + 1] = Chi_star_R*Q_R[d_num_species + 1];
            Q_star_R[d_num_species + d_dim.getValue()] = Chi_star_R*(Q_R[d_num_species + d_dim.getValue()] +
                (s_star - (V_R[d_num_species].get()))*(rho_R*s_star +
                        (V_R[d_num_species + d_dim.getValue()].get())/(s_R - (V_R[d_num_species].get()))));
            
            double F_x_R[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                F_x_R[si] = (V_R[d_num_species].get())*(V_R[si].get());
            }
            F_x_R[d_num_species] = (V_R[d_num_species].get())*Q_R[d_num_species] +
                (V_R[d_num_species + d_dim.getValue()].get());
            F_x_R[d_num_species + 1] = (V_R[d_num_species].get())*Q_R[d_num_species + 1];
            F_x_R[d_num_species + d_dim.getValue()] = (V_R[d_num_species].get())*(Q_R[d_num_species + d_dim.getValue()] +
                (V_R[d_num_species + d_dim.getValue()].get()));
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_x_intercell[ei].get() = F_x_R[ei] + s_plus*(Q_star_R[ei] - Q_R[ei]);
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        std::vector<const double*> rho_Y_L;
        std::vector<const double*> rho_Y_R;
        rho_Y_L.reserve(d_num_species);
        rho_Y_R.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y_L.push_back(&(V_L[si].get()));
            rho_Y_R.push_back(&(V_R[si].get()));
        }
        
        std::vector<const double*> vel_L;
        std::vector<const double*> vel_R;
        vel_L.reserve(3);
        vel_R.reserve(3);
        vel_L.push_back(&(V_L[d_num_species].get()));
        vel_R.push_back(&(V_R[d_num_species].get()));
        vel_L.push_back(&(V_L[d_num_species + 1].get()));
        vel_R.push_back(&(V_R[d_num_species + 1].get()));
        vel_L.push_back(&(V_L[d_num_species + 2].get()));
        vel_R.push_back(&(V_R[d_num_species + 2].get()));
        
        const double rho_L = d_equation_of_state_mixing_rules->getMixtureDensity(
            rho_Y_L);
        
        const double rho_R = d_equation_of_state_mixing_rules->getMixtureDensity(
            rho_Y_R);
        
        std::vector<const double*> empty_ptr;
        
        const double c_L = d_equation_of_state_mixing_rules->getSoundSpeed(
            rho_Y_L,
            vel_L,
            &(V_L[d_num_species + d_dim.getValue()].get()),
            empty_ptr);
        
        const double c_R = d_equation_of_state_mixing_rules->getSoundSpeed(
            rho_Y_R,
            vel_R,
            &(V_R[d_num_species + d_dim.getValue()].get()),
            empty_ptr);
        
        const double u_average = 0.5*((V_L[d_num_species].get()) + (V_R[d_num_species].get()));
        const double c_average = 0.5*(c_L + c_R);
        
        const double s_L = fmin(u_average - c_average, (V_L[d_num_species].get()) - c_L);
        const double s_R = fmax(u_average + c_average, (V_R[d_num_species].get()) + c_R);
        
        const double s_minus = fmin(0.0, s_L);
        const double s_plus  = fmax(0.0, s_R);
        
        const double s_star =
            ((V_R[d_num_species + d_dim.getValue()].get()) - (V_L[d_num_species + d_dim.getValue()].get()) +
                rho_L*(V_L[d_num_species].get())*(s_L - (V_L[d_num_species].get())) -
                    rho_R*(V_R[d_num_species].get())*(s_R - (V_R[d_num_species].get())))/
                        (rho_L*(s_L - (V_L[d_num_species].get())) - rho_R*(s_R - (V_R[d_num_species].get())));
        
        if (s_star > 0)
        {
            const double Chi_star_L = (s_L - (V_L[d_num_species].get()))/(s_L - s_star);
            
            std::vector<const double*> vel_L;
            vel_L.reserve(3);
            vel_L.push_back(&(V_L[d_num_species].get()));
            vel_L.push_back(&(V_L[d_num_species + 1].get()));
            vel_L.push_back(&(V_L[d_num_species + 2].get()));
            
            double Q_L[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                Q_L[si] = (V_L[si].get());
            }
            Q_L[d_num_species] = rho_L*(V_L[d_num_species].get());
            Q_L[d_num_species + 1] = rho_L*(V_L[d_num_species + 1].get());
            Q_L[d_num_species + 2] = rho_L*(V_L[d_num_species + 2].get());
            Q_L[d_num_species + d_dim.getValue()] = d_equation_of_state_mixing_rules->getTotalEnergy(
                rho_Y_L,
                vel_L,
                &(V_L[d_num_species + d_dim.getValue()].get()),
                empty_ptr);
            
            double Q_star_L[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                Q_star_L[si] = Chi_star_L*(V_L[si].get());
            }
            Q_star_L[d_num_species] = Chi_star_L*rho_L*s_star;
            Q_star_L[d_num_species + 1] = Chi_star_L*Q_L[d_num_species + 1];
            Q_star_L[d_num_species + 2] = Chi_star_L*Q_L[d_num_species + 2];
            Q_star_L[d_num_species + d_dim.getValue()] = Chi_star_L*(Q_L[d_num_species + d_dim.getValue()] +
                (s_star - (V_L[d_num_species].get()))*(rho_L*s_star +
                        (V_L[d_num_species + d_dim.getValue()].get())/(s_L - (V_L[d_num_species].get()))));
            
            double F_x_L[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                F_x_L[si] = (V_L[d_num_species].get())*(V_L[si].get());
            }
            F_x_L[d_num_species] = (V_L[d_num_species].get())*Q_L[d_num_species] +
                (V_L[d_num_species + d_dim.getValue()].get());
            F_x_L[d_num_species + 1] = (V_L[d_num_species].get())*Q_L[d_num_species + 1];
            F_x_L[d_num_species + 2] = (V_L[d_num_species].get())*Q_L[d_num_species + 2];
            F_x_L[d_num_species + d_dim.getValue()] = (V_L[d_num_species].get())*(Q_L[d_num_species + d_dim.getValue()] +
                (V_L[d_num_species + d_dim.getValue()].get()));
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_x_intercell[ei].get() = F_x_L[ei] + s_minus*(Q_star_L[ei] - Q_L[ei]);
            }
        }
        else
        {
            const double Chi_star_R = (s_R - (V_R[d_num_species].get()))/(s_R - s_star);
            
            std::vector<const double*> vel_R;
            vel_R.reserve(3);
            vel_R.push_back(&(V_R[d_num_species].get()));
            vel_R.push_back(&(V_R[d_num_species + 1].get()));
            vel_R.push_back(&(V_R[d_num_species + 2].get()));
            
            double Q_R[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                Q_R[si] = (V_R[si].get());
            }
            Q_R[d_num_species] = rho_R*(V_R[d_num_species].get());
            Q_R[d_num_species + 1] = rho_R*(V_R[d_num_species + 1].get());
            Q_R[d_num_species + 2] = rho_R*(V_R[d_num_species + 2].get());
            Q_R[d_num_species + d_dim.getValue()] = d_equation_of_state_mixing_rules->getTotalEnergy(
                rho_Y_R,
                vel_R,
                &(V_R[d_num_species + d_dim.getValue()].get()),
                empty_ptr);
            
            double Q_star_R[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                Q_star_R[si] = Chi_star_R*(V_R[si].get());
            }
            Q_star_R[d_num_species] = Chi_star_R*rho_R*s_star;
            Q_star_R[d_num_species + 1] = Chi_star_R*Q_R[d_num_species + 1];
            Q_star_R[d_num_species + 2] = Chi_star_R*Q_R[d_num_species + 2];
            Q_star_R[d_num_species + d_dim.getValue()] = Chi_star_R*(Q_R[d_num_species + d_dim.getValue()] +
                (s_star - (V_R[d_num_species].get()))*(rho_R*s_star +
                        (V_R[d_num_species + d_dim.getValue()].get())/(s_R - (V_R[d_num_species].get()))));
            
            double F_x_R[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                F_x_R[si] = (V_R[d_num_species].get())*(V_R[si].get());
            }
            F_x_R[d_num_species] = (V_R[d_num_species].get())*Q_R[d_num_species] +
                (V_R[d_num_species + d_dim.getValue()].get());
            F_x_R[d_num_species + 1] = (V_R[d_num_species].get())*Q_R[d_num_species + 1];
            F_x_R[d_num_species + 2] = (V_R[d_num_species].get())*Q_R[d_num_species + 2];
            F_x_R[d_num_species + d_dim.getValue()] = (V_R[d_num_species].get())*(Q_R[d_num_species + d_dim.getValue()] +
                (V_R[d_num_species + d_dim.getValue()].get()));
            
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
RiemannSolverFourEqnConservativeHLLC::computeIntercellFluxInYDirectionFromPrimitiveVariables(
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
            << ": RiemannSolverFourEqnConservativeHLLC::"
            << "computeIntercellFluxInYDirectionFromPrimitiveVariables()\n"
            << "There is no y direction for 1D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        std::vector<const double*> rho_Y_B;
        std::vector<const double*> rho_Y_T;
        rho_Y_B.reserve(d_num_species);
        rho_Y_T.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y_B.push_back(&(V_B[si].get()));
            rho_Y_T.push_back(&(V_T[si].get()));
        }
        
        std::vector<const double*> vel_B;
        std::vector<const double*> vel_T;
        vel_B.reserve(2);
        vel_T.reserve(2);
        vel_B.push_back(&(V_B[d_num_species].get()));
        vel_T.push_back(&(V_T[d_num_species].get()));
        vel_B.push_back(&(V_B[d_num_species + 1].get()));
        vel_T.push_back(&(V_T[d_num_species + 1].get()));
        
        const double rho_B = d_equation_of_state_mixing_rules->getMixtureDensity(
            rho_Y_B);
        
        const double rho_T = d_equation_of_state_mixing_rules->getMixtureDensity(
            rho_Y_T);
        
        std::vector<const double*> empty_ptr;
        
        const double c_B = d_equation_of_state_mixing_rules->getSoundSpeed(
            rho_Y_B,
            vel_B,
            &(V_B[d_num_species + d_dim.getValue()].get()),
            empty_ptr);
        
        const double c_T = d_equation_of_state_mixing_rules->getSoundSpeed(
            rho_Y_T,
            vel_T,
            &(V_T[d_num_species + d_dim.getValue()].get()),
            empty_ptr);
        
        const double v_average = 0.5*((V_B[d_num_species + 1].get()) + (V_T[d_num_species + 1].get()));
        const double c_average = 0.5*(c_B + c_T);
        
        const double s_B = fmin(v_average - c_average, (V_B[d_num_species + 1].get()) - c_B);
        const double s_T = fmax(v_average + c_average, (V_T[d_num_species + 1].get()) + c_T);
        
        const double s_minus = fmin(0.0, s_B);
        const double s_plus  = fmax(0.0, s_T);
        
        const double s_star =
            ((V_T[d_num_species + d_dim.getValue()].get()) - (V_B[d_num_species + d_dim.getValue()].get()) +
                rho_B*(V_B[d_num_species + 1].get())*(s_B - (V_B[d_num_species + 1].get())) -
                    rho_T*(V_T[d_num_species + 1].get())*(s_T - (V_T[d_num_species + 1].get())))/
                        (rho_B*(s_B - (V_B[d_num_species + 1].get())) - rho_T*(s_T - (V_T[d_num_species + 1].get())));
        
        if (s_star > 0)
        {
            const double Chi_star_B = (s_B - (V_B[d_num_species + 1].get()))/(s_B - s_star);
            
            std::vector<const double*> vel_B;
            vel_B.reserve(2);
            vel_B.push_back(&(V_B[d_num_species].get()));
            vel_B.push_back(&(V_B[d_num_species + 1].get()));
            
            double Q_B[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                Q_B[si] = (V_B[si].get());
            }
            Q_B[d_num_species] = rho_B*(V_B[d_num_species].get());
            Q_B[d_num_species + 1] = rho_B*(V_B[d_num_species + 1].get());
            Q_B[d_num_species + d_dim.getValue()] = d_equation_of_state_mixing_rules->getTotalEnergy(
                rho_Y_B,
                vel_B,
                &(V_B[d_num_species + d_dim.getValue()].get()),
                empty_ptr);
            
            double Q_star_B[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                Q_star_B[si] = Chi_star_B*(V_B[si].get());
            }
            Q_star_B[d_num_species] = Chi_star_B*Q_B[d_num_species];
            Q_star_B[d_num_species + 1] = Chi_star_B*rho_B*s_star;
            Q_star_B[d_num_species + d_dim.getValue()] = Chi_star_B*(Q_B[d_num_species + d_dim.getValue()] +
                (s_star - (V_B[d_num_species + 1].get()))*(rho_B*s_star +
                        (V_B[d_num_species + d_dim.getValue()].get())/(s_B - (V_B[d_num_species + 1].get()))));
            
            double F_y_B[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                F_y_B[si] = (V_B[d_num_species + 1].get())*(V_B[si].get());
            }
            F_y_B[d_num_species] = (V_B[d_num_species + 1].get())*Q_B[d_num_species];
            F_y_B[d_num_species + 1] = (V_B[d_num_species + 1].get())*Q_B[d_num_species + 1] +
                (V_B[d_num_species + d_dim.getValue()].get());
            F_y_B[d_num_species + d_dim.getValue()] = (V_B[d_num_species + 1].get())*(Q_B[d_num_species + d_dim.getValue()] +
                (V_B[d_num_species + d_dim.getValue()].get()));
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_y_intercell[ei].get() = F_y_B[ei] + s_minus*(Q_star_B[ei] - Q_B[ei]);
            }
        }
        else
        {
            const double Chi_star_T = (s_T - (V_T[d_num_species + 1].get()))/(s_T - s_star);
            
            std::vector<const double*> vel_T;
            vel_T.reserve(2);
            vel_T.push_back(&(V_T[d_num_species].get()));
            vel_T.push_back(&(V_T[d_num_species + 1].get()));
            
            double Q_T[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                Q_T[si] = (V_T[si].get());
            }
            Q_T[d_num_species] = rho_T*(V_T[d_num_species].get());
            Q_T[d_num_species + 1] = rho_T*(V_T[d_num_species + 1].get());
            Q_T[d_num_species + d_dim.getValue()] = d_equation_of_state_mixing_rules->getTotalEnergy(
                rho_Y_T,
                vel_T,
                &(V_T[d_num_species + d_dim.getValue()].get()),
                empty_ptr);
            
            double Q_star_T[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                Q_star_T[si] = Chi_star_T*(V_T[si].get());
            }
            Q_star_T[d_num_species] = Chi_star_T*Q_T[d_num_species];
            Q_star_T[d_num_species + 1] = Chi_star_T*rho_T*s_star;
            Q_star_T[d_num_species + d_dim.getValue()] = Chi_star_T*(Q_T[d_num_species + d_dim.getValue()] +
                (s_star - (V_T[d_num_species + 1].get()))*(rho_T*s_star +
                        (V_T[d_num_species + d_dim.getValue()].get())/(s_T - (V_T[d_num_species + 1].get()))));
            
            double F_y_T[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                F_y_T[si] = (V_T[d_num_species + 1].get())*(V_T[si].get());
            }
            F_y_T[d_num_species] = (V_T[d_num_species + 1].get())*Q_T[d_num_species];
            F_y_T[d_num_species + 1] = (V_T[d_num_species + 1].get())*Q_T[d_num_species + 1] +
                (V_T[d_num_species + d_dim.getValue()].get());
            F_y_T[d_num_species + d_dim.getValue()] = (V_T[d_num_species + 1].get())*(Q_T[d_num_species + d_dim.getValue()] +
                (V_T[d_num_species + d_dim.getValue()].get()));
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_y_intercell[ei].get() = F_y_T[ei] + s_plus*(Q_star_T[ei] - Q_T[ei]);
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        std::vector<const double*> rho_Y_B;
        std::vector<const double*> rho_Y_T;
        rho_Y_B.reserve(d_num_species);
        rho_Y_T.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y_B.push_back(&(V_B[si].get()));
            rho_Y_T.push_back(&(V_T[si].get()));
        }
        
        std::vector<const double*> vel_B;
        std::vector<const double*> vel_T;
        vel_B.reserve(3);
        vel_T.reserve(3);
        vel_B.push_back(&(V_B[d_num_species].get()));
        vel_T.push_back(&(V_T[d_num_species].get()));
        vel_B.push_back(&(V_B[d_num_species + 1].get()));
        vel_T.push_back(&(V_T[d_num_species + 1].get()));
        vel_B.push_back(&(V_B[d_num_species + 2].get()));
        vel_T.push_back(&(V_T[d_num_species + 2].get()));
        
        const double rho_B = d_equation_of_state_mixing_rules->getMixtureDensity(
            rho_Y_B);
        
        const double rho_T = d_equation_of_state_mixing_rules->getMixtureDensity(
            rho_Y_T);
        
        std::vector<const double*> empty_ptr;
        
        const double c_B = d_equation_of_state_mixing_rules->getSoundSpeed(
            rho_Y_B,
            vel_B,
            &(V_B[d_num_species + d_dim.getValue()].get()),
            empty_ptr);
        
        const double c_T = d_equation_of_state_mixing_rules->getSoundSpeed(
            rho_Y_T,
            vel_T,
            &(V_T[d_num_species + d_dim.getValue()].get()),
            empty_ptr);
        
        const double v_average = 0.5*((V_B[d_num_species + 1].get()) + (V_T[d_num_species + 1].get()));
        const double c_average = 0.5*(c_B + c_T);
        
        const double s_B = fmin(v_average - c_average, (V_B[d_num_species + 1].get()) - c_B);
        const double s_T = fmax(v_average + c_average, (V_T[d_num_species + 1].get()) + c_T);
        
        const double s_minus = fmin(0.0, s_B);
        const double s_plus  = fmax(0.0, s_T);
        
        const double s_star =
            ((V_T[d_num_species + d_dim.getValue()].get()) - (V_B[d_num_species + d_dim.getValue()].get()) +
                rho_B*(V_B[d_num_species + 1].get())*(s_B - (V_B[d_num_species + 1].get())) -
                    rho_T*(V_T[d_num_species + 1].get())*(s_T - (V_T[d_num_species + 1].get())))/
                        (rho_B*(s_B - (V_B[d_num_species + 1].get())) - rho_T*(s_T - (V_T[d_num_species + 1].get())));
        
        if (s_star > 0)
        {
            const double Chi_star_B = (s_B - (V_B[d_num_species + 1].get()))/(s_B - s_star);
            
            std::vector<const double*> vel_B;
            vel_B.reserve(3);
            vel_B.push_back(&(V_B[d_num_species].get()));
            vel_B.push_back(&(V_B[d_num_species + 1].get()));
            vel_B.push_back(&(V_B[d_num_species + 2].get()));
            
            double Q_B[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                Q_B[si] = (V_B[si].get());
            }
            Q_B[d_num_species] = rho_B*(V_B[d_num_species].get());
            Q_B[d_num_species + 1] = rho_B*(V_B[d_num_species + 1].get());
            Q_B[d_num_species + 2] = rho_B*(V_B[d_num_species + 2].get());
            Q_B[d_num_species + d_dim.getValue()] = d_equation_of_state_mixing_rules->getTotalEnergy(
                rho_Y_B,
                vel_B,
                &(V_B[d_num_species + d_dim.getValue()].get()),
                empty_ptr);
            
            double Q_star_B[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                Q_star_B[si] = Chi_star_B*(V_B[si].get());
            }
            Q_star_B[d_num_species] = Chi_star_B*Q_B[d_num_species];
            Q_star_B[d_num_species + 1] = Chi_star_B*rho_B*s_star;
            Q_star_B[d_num_species + 2] = Chi_star_B*Q_B[d_num_species + 2];
            Q_star_B[d_num_species + d_dim.getValue()] = Chi_star_B*(Q_B[d_num_species + d_dim.getValue()] +
                (s_star - (V_B[d_num_species + 1].get()))*(rho_B*s_star +
                        (V_B[d_num_species + d_dim.getValue()].get())/(s_B - (V_B[d_num_species + 1].get()))));
            
            double F_y_B[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                F_y_B[si] = (V_B[d_num_species + 1].get())*(V_B[si].get());
            }
            F_y_B[d_num_species] = (V_B[d_num_species + 1].get())*Q_B[d_num_species];
            F_y_B[d_num_species + 1] = (V_B[d_num_species + 1].get())*Q_B[d_num_species + 1] +
                (V_B[d_num_species + d_dim.getValue()].get());
            F_y_B[d_num_species + 2] = (V_B[d_num_species + 1].get())*Q_B[d_num_species + 2];
            F_y_B[d_num_species + d_dim.getValue()] = (V_B[d_num_species + 1].get())*(Q_B[d_num_species + d_dim.getValue()] +
                (V_B[d_num_species + d_dim.getValue()].get()));
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_y_intercell[ei].get() = F_y_B[ei] + s_minus*(Q_star_B[ei] - Q_B[ei]);
            }
        }
        else
        {
            const double Chi_star_T = (s_T - (V_T[d_num_species + 1].get()))/(s_T - s_star);
            
            std::vector<const double*> vel_T;
            vel_T.reserve(3);
            vel_T.push_back(&(V_T[d_num_species].get()));
            vel_T.push_back(&(V_T[d_num_species + 1].get()));
            vel_T.push_back(&(V_T[d_num_species + 2].get()));
            
            double Q_T[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                Q_T[si] = (V_T[si].get());
            }
            Q_T[d_num_species] = rho_T*(V_T[d_num_species].get());
            Q_T[d_num_species + 1] = rho_T*(V_T[d_num_species + 1].get());
            Q_T[d_num_species + 2] = rho_T*(V_T[d_num_species + 2].get());
            Q_T[d_num_species + d_dim.getValue()] = d_equation_of_state_mixing_rules->getTotalEnergy(
                rho_Y_T,
                vel_T,
                &(V_T[d_num_species + d_dim.getValue()].get()),
                empty_ptr);
            
            double Q_star_T[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                Q_star_T[si] = Chi_star_T*(V_T[si].get());
            }
            Q_star_T[d_num_species] = Chi_star_T*Q_T[d_num_species];
            Q_star_T[d_num_species + 1] = Chi_star_T*rho_T*s_star;
            Q_star_T[d_num_species + 2] = Chi_star_T*Q_T[d_num_species + 2];
            Q_star_T[d_num_species + d_dim.getValue()] = Chi_star_T*(Q_T[d_num_species + d_dim.getValue()] +
                (s_star - (V_T[d_num_species + 1].get()))*(rho_T*s_star +
                        (V_T[d_num_species + d_dim.getValue()].get())/(s_T - (V_T[d_num_species + 1].get()))));
            
            double F_y_T[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                F_y_T[si] = (V_T[d_num_species + 1].get())*(V_T[si].get());
            }
            F_y_T[d_num_species] = (V_T[d_num_species + 1].get())*Q_T[d_num_species];
            F_y_T[d_num_species + 1] = (V_T[d_num_species + 1].get())*Q_T[d_num_species + 1] +
                (V_T[d_num_species + d_dim.getValue()].get());
            F_y_T[d_num_species + 2] = (V_T[d_num_species + 1].get())*Q_T[d_num_species + 2];
            F_y_T[d_num_species + d_dim.getValue()] = (V_T[d_num_species + 1].get())*(Q_T[d_num_species + d_dim.getValue()] +
                (V_T[d_num_species + d_dim.getValue()].get()));
            
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
RiemannSolverFourEqnConservativeHLLC::computeIntercellFluxInZDirectionFromPrimitiveVariables(
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
            << ": RiemannSolverFourEqnConservativeHLLC::"
            << "computeIntercellFluxInZDirectionFromPrimitiveVariables()\n"
            << "There is no z direction for 1D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": RiemannSolverFourEqnConservativeHLLC::"
            << "computeIntercellFluxInZDirectionFromPrimitiveVariables()\n"
            << "There is no z direction for 2D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        std::vector<const double*> rho_Y_B;
        std::vector<const double*> rho_Y_F;
        rho_Y_B.reserve(d_num_species);
        rho_Y_F.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y_B.push_back(&(V_B[si].get()));
            rho_Y_F.push_back(&(V_F[si].get()));
        }
        
        std::vector<const double*> vel_B;
        std::vector<const double*> vel_F;
        vel_B.reserve(3);
        vel_F.reserve(3);
        vel_B.push_back(&(V_B[d_num_species].get()));
        vel_F.push_back(&(V_F[d_num_species].get()));
        vel_B.push_back(&(V_B[d_num_species + 1].get()));
        vel_F.push_back(&(V_F[d_num_species + 1].get()));
        vel_B.push_back(&(V_B[d_num_species + 2].get()));
        vel_F.push_back(&(V_F[d_num_species + 2].get()));
        
        const double rho_B = d_equation_of_state_mixing_rules->getMixtureDensity(
            rho_Y_B);
        
        const double rho_F = d_equation_of_state_mixing_rules->getMixtureDensity(
            rho_Y_F);
        
        std::vector<const double*> empty_ptr;
        
        const double c_B = d_equation_of_state_mixing_rules->getSoundSpeed(
            rho_Y_B,
            vel_B,
            &(V_B[d_num_species + d_dim.getValue()].get()),
            empty_ptr);
        
        const double c_F = d_equation_of_state_mixing_rules->getSoundSpeed(
            rho_Y_F,
            vel_F,
            &(V_F[d_num_species + d_dim.getValue()].get()),
            empty_ptr);
        
        const double w_average = 0.5*((V_B[d_num_species + 2].get()) + (V_F[d_num_species + 2].get()));
        const double c_average = 0.5*(c_B + c_F);
        
        const double s_B = fmin(w_average - c_average, (V_B[d_num_species + 2].get()) - c_B);
        const double s_F = fmax(w_average + c_average, (V_F[d_num_species + 2].get()) + c_F);
        
        const double s_minus = fmin(0.0, s_B);
        const double s_plus  = fmax(0.0, s_F);
        
        const double s_star =
            ((V_F[d_num_species + d_dim.getValue()].get()) - (V_B[d_num_species + d_dim.getValue()].get()) +
                rho_B*(V_B[d_num_species + 2].get())*(s_B - (V_B[d_num_species + 2].get())) -
                    rho_F*(V_F[d_num_species + 2].get())*(s_F - (V_F[d_num_species + 2].get())))/
                        (rho_B*(s_B - (V_B[d_num_species + 2].get())) - rho_F*(s_F - (V_F[d_num_species + 2].get())));
        
        if (s_star > 0)
        {
            const double Chi_star_B = (s_B - (V_B[d_num_species + 2].get()))/(s_B - s_star);
            
            std::vector<const double*> vel_B;
            vel_B.reserve(3);
            vel_B.push_back(&(V_B[d_num_species].get()));
            vel_B.push_back(&(V_B[d_num_species + 1].get()));
            vel_B.push_back(&(V_B[d_num_species + 2].get()));
            
            double Q_B[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                Q_B[si] = (V_B[si].get());
            }
            Q_B[d_num_species] = rho_B*(V_B[d_num_species].get());
            Q_B[d_num_species + 1] = rho_B*(V_B[d_num_species + 1].get());
            Q_B[d_num_species + 2] = rho_B*(V_B[d_num_species + 2].get());
            Q_B[d_num_species + d_dim.getValue()] = d_equation_of_state_mixing_rules->getTotalEnergy(
                rho_Y_B,
                vel_B,
                &(V_B[d_num_species + d_dim.getValue()].get()),
                empty_ptr);
            
            double Q_star_B[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                Q_star_B[si] = Chi_star_B*(V_B[si].get());
            }
            Q_star_B[d_num_species] = Chi_star_B*Q_B[d_num_species];
            Q_star_B[d_num_species + 1] = Chi_star_B*Q_B[d_num_species + 1];
            Q_star_B[d_num_species + 2] = Chi_star_B*rho_B*s_star;
            Q_star_B[d_num_species + d_dim.getValue()] = Chi_star_B*(Q_B[d_num_species + d_dim.getValue()] +
                (s_star - (V_B[d_num_species + 2].get()))*(rho_B*s_star +
                        (V_B[d_num_species + d_dim.getValue()].get())/(s_B - (V_B[d_num_species + 2].get()))));
            
            double F_z_B[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                F_z_B[si] = (V_B[d_num_species + 2].get())*(V_B[si].get());
            }
            F_z_B[d_num_species] = (V_B[d_num_species + 2].get())*Q_B[d_num_species];
            F_z_B[d_num_species + 1] = (V_B[d_num_species + 2].get())*Q_B[d_num_species + 1];
            F_z_B[d_num_species + 2] = (V_B[d_num_species + 2].get())*Q_B[d_num_species + 2] +
                (V_B[d_num_species + d_dim.getValue()].get());
            F_z_B[d_num_species + d_dim.getValue()] = (V_B[d_num_species + 2].get())*(Q_B[d_num_species + d_dim.getValue()] +
                (V_B[d_num_species + d_dim.getValue()].get()));
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_z_intercell[ei].get() = F_z_B[ei] + s_minus*(Q_star_B[ei] - Q_B[ei]);
            }
        }
        else
        {
            const double Chi_star_F = (s_F - (V_F[d_num_species + 2].get()))/(s_F - s_star);
            
            std::vector<const double*> vel_F;
            vel_F.reserve(3);
            vel_F.push_back(&(V_F[d_num_species].get()));
            vel_F.push_back(&(V_F[d_num_species + 1].get()));
            vel_F.push_back(&(V_F[d_num_species + 2].get()));
            
            double Q_F[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                Q_F[si] = (V_F[si].get());
            }
            Q_F[d_num_species] = rho_F*(V_F[d_num_species].get());
            Q_F[d_num_species + 1] = rho_F*(V_F[d_num_species + 1].get());
            Q_F[d_num_species + 2] = rho_F*(V_F[d_num_species + 2].get());
            Q_F[d_num_species + d_dim.getValue()] = d_equation_of_state_mixing_rules->getTotalEnergy(
                rho_Y_F,
                vel_F,
                &(V_F[d_num_species + d_dim.getValue()].get()),
                empty_ptr);
            
            double Q_star_F[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                Q_star_F[si] = Chi_star_F*(V_F[si].get());
            }
            Q_star_F[d_num_species] = Chi_star_F*Q_F[d_num_species];
            Q_star_F[d_num_species + 1] = Chi_star_F*Q_F[d_num_species + 1];
            Q_star_F[d_num_species + 2] = Chi_star_F*rho_F*s_star;
            Q_star_F[d_num_species + d_dim.getValue()] = Chi_star_F*(Q_F[d_num_species + d_dim.getValue()] +
                (s_star - (V_F[d_num_species + 2].get()))*(rho_F*s_star +
                        (V_F[d_num_species + d_dim.getValue()].get())/(s_F - (V_F[d_num_species + 2].get()))));
            
            double F_z_F[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                F_z_F[si] = (V_F[d_num_species + 2].get())*(V_F[si].get());
            }
            F_z_F[d_num_species] = (V_F[d_num_species + 2].get())*Q_F[d_num_species];
            F_z_F[d_num_species + 1] = (V_F[d_num_species + 2].get())*Q_F[d_num_species + 1];
            F_z_F[d_num_species + 2] = (V_F[d_num_species + 2].get())*Q_F[d_num_species + 2] +
                (V_F[d_num_species + d_dim.getValue()].get());
            F_z_F[d_num_species + d_dim.getValue()] = (V_F[d_num_species + 2].get())*(Q_F[d_num_species + d_dim.getValue()] +
                (V_F[d_num_species + d_dim.getValue()].get()));
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_z_intercell[ei].get() = F_z_F[ei] + s_plus*(Q_star_F[ei] - Q_F[ei]);
            }
        }
    }
}
