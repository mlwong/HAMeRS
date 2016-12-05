#include "flow/flow_models/four-eqn_conservative/Riemann_solvers/RiemannSolverFourEqnConservativeHLLC_HLL.hpp"

#define EPSILON 1e-40

/*
 * Compute the flux at the intercell face from conservative variables.
 */
void
RiemannSolverFourEqnConservativeHLLC_HLL::computeIntercellFluxFromConservativeVariables(
    std::vector<boost::reference_wrapper<double> >& flux_intercell,
    const std::vector<boost::reference_wrapper<double> >& conservative_variables_minus,
    const std::vector<boost::reference_wrapper<double> >& conservative_variables_plus,
    const DIRECTION::TYPE& direction)
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(flux_intercell.size()) == d_num_eqn);
    TBOX_ASSERT(static_cast<int>(conservative_variables_minus.size()) == d_num_eqn);
    TBOX_ASSERT(static_cast<int>(conservative_variables_plus.size()) == d_num_eqn);
#endif
    
    switch (direction)
    {
        case DIRECTION::X_DIRECTION:
        {
            computeIntercellFluxInXDirectionFromConservativeVariables(
                flux_intercell,
                conservative_variables_minus,
                conservative_variables_plus);
            
            break;
        }
        case DIRECTION::Y_DIRECTION:
        {
            computeIntercellFluxInYDirectionFromConservativeVariables(
                flux_intercell,
                conservative_variables_minus,
                conservative_variables_plus);
            
            break;
        }
        case DIRECTION::Z_DIRECTION:
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
                << ": RiemannSolverFourEqnConservativeHLLC_HLL::"
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
RiemannSolverFourEqnConservativeHLLC_HLL::computeIntercellFluxFromPrimitiveVariables(
    std::vector<boost::reference_wrapper<double> >& flux_intercell,
    const std::vector<boost::reference_wrapper<double> >& primitive_variables_minus,
    const std::vector<boost::reference_wrapper<double> >& primitive_variables_plus,
    const DIRECTION::TYPE& direction)
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(flux_intercell.size()) == d_num_eqn);
    TBOX_ASSERT(static_cast<int>(primitive_variables_minus.size()) == d_num_eqn);
    TBOX_ASSERT(static_cast<int>(primitive_variables_plus.size()) == d_num_eqn);
#endif
    
    switch (direction)
    {
        case DIRECTION::X_DIRECTION:
        {
            computeIntercellFluxInXDirectionFromPrimitiveVariables(
                flux_intercell,
                primitive_variables_minus,
                primitive_variables_plus);
            
            break;
        }
        case DIRECTION::Y_DIRECTION:
        {
            computeIntercellFluxInYDirectionFromPrimitiveVariables(
                flux_intercell,
                primitive_variables_minus,
                primitive_variables_plus);
            
            break;
        }
        case DIRECTION::Z_DIRECTION:
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
                << ": RiemannSolverFourEqnConservativeHLLC_HLL::"
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
RiemannSolverFourEqnConservativeHLLC_HLL::computeIntercellFluxInXDirectionFromConservativeVariables(
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
        
        const double epsilon_L =
            ((Q_L[d_num_species + d_dim.getValue()].get()) -
                0.5*(Q_L[d_num_species].get())*(Q_L[d_num_species].get())/rho_L)/rho_L;
        
        const double epsilon_R =
            ((Q_R[d_num_species + d_dim.getValue()].get()) -
                0.5*(Q_R[d_num_species].get())*(Q_R[d_num_species].get())/rho_R)/rho_R;
        
        /*
         * Compute the mass fractions.
         */
        double Y_L[d_num_species];
        double Y_R[d_num_species];
        for (int si = 0; si < d_num_species; si++)
        {
            Y_L[si] = Q_L[si]/rho_L;
            Y_R[si] = Q_R[si]/rho_R;
        }
        
        /*
         * Get the pointers to the mass fractions.
         */
        std::vector<const double*> Y_L_ptr;
        std::vector<const double*> Y_R_ptr;
        Y_L_ptr.reserve(d_num_species);
        Y_R_ptr.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            Y_L_ptr.push_back(&Y_L[si]);
            Y_R_ptr.push_back(&Y_R[si]);
        }
        
        const double p_L = d_equation_of_state_mixing_rules->getPressure(
            &rho_L,
            &epsilon_L,
            Y_L_ptr);
        
        const double p_R = d_equation_of_state_mixing_rules->getPressure(
            &rho_R,
            &epsilon_R,
            Y_R_ptr);
        
        const double c_L = d_equation_of_state_mixing_rules->getSoundSpeed(
            &rho_L,
            &p_L,
            Y_L_ptr);
        
        const double c_R = d_equation_of_state_mixing_rules->getSoundSpeed(
            &rho_R,
            &p_R,
            Y_R_ptr);
        
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
        double F_x_intercell_HLLC[d_num_eqn];
        double F_x_intercell_HLL[d_num_eqn];
        
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
        
        const double epsilon_L =
            ((Q_L[d_num_species + d_dim.getValue()].get()) -
                0.5*((Q_L[d_num_species].get())*(Q_L[d_num_species].get()) +
                (Q_L[d_num_species + 1].get())*(Q_L[d_num_species + 1].get()))/
                rho_L)/rho_L;
        
        const double epsilon_R =
            ((Q_R[d_num_species + d_dim.getValue()].get()) -
                0.5*((Q_R[d_num_species].get())*(Q_R[d_num_species].get()) +
                (Q_R[d_num_species + 1].get())*(Q_R[d_num_species + 1].get()))/
                rho_R)/rho_R;
        
        /*
         * Compute the mass fractions.
         */
        double Y_L[d_num_species];
        double Y_R[d_num_species];
        for (int si = 0; si < d_num_species; si++)
        {
            Y_L[si] = Q_L[si]/rho_L;
            Y_R[si] = Q_R[si]/rho_R;
        }
        
        /*
         * Get the pointers to the mass fractions.
         */
        std::vector<const double*> Y_L_ptr;
        std::vector<const double*> Y_R_ptr;
        Y_L_ptr.reserve(d_num_species);
        Y_R_ptr.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            Y_L_ptr.push_back(&Y_L[si]);
            Y_R_ptr.push_back(&Y_R[si]);
        }
        
        const double p_L = d_equation_of_state_mixing_rules->getPressure(
            &rho_L,
            &epsilon_L,
            Y_L_ptr);
        
        const double p_R = d_equation_of_state_mixing_rules->getPressure(
            &rho_R,
            &epsilon_R,
            Y_R_ptr);
        
        const double c_L = d_equation_of_state_mixing_rules->getSoundSpeed(
            &rho_L,
            &p_L,
            Y_L_ptr);
        
        const double c_R = d_equation_of_state_mixing_rules->getSoundSpeed(
            &rho_R,
            &p_R,
            Y_R_ptr);
        
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
            /*
             * Compute the HLLC flux.
             */
            
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
                F_x_intercell_HLLC[ei] = F_x_L[ei] + s_minus*(Q_star_L[ei] - Q_L[ei]);
            }
            
            /*
             * Compute the HLL flux.
             */
            
            if (s_L >= 0)
            {
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    F_x_intercell_HLL[ei] = F_x_L[ei];
                }
            }
            else
            {
                double F_x_R[d_num_eqn];
                for (int si = 0; si < d_num_species; si++)
                {
                    F_x_R[si] = u_R*(Q_R[si].get());
                }
                F_x_R[d_num_species] = u_R*(Q_R[d_num_species].get()) + p_R;
                F_x_R[d_num_species + 1] = u_R*(Q_R[d_num_species + 1].get());
                F_x_R[d_num_species + 2] = u_R*((Q_R[d_num_species + 2].get()) + p_R);
                
                if (s_R <= 0)
                {
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        F_x_intercell_HLL[ei] = F_x_R[ei];
                    }
                }
                else
                {
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        F_x_intercell_HLL[ei] = (s_R*F_x_L[ei] - s_L*F_x_R[ei] +
                            s_R*s_L*(Q_R[ei] - Q_L[ei]))/(s_R - s_L);
                    }
                }
            }
        }
        else
        {
            /*
             * Compute the HLLC flux.
             */
            
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
                F_x_intercell_HLLC[ei] = F_x_R[ei] + s_plus*(Q_star_R[ei] - Q_R[ei]);
            }
            
            /*
             * Compute the HLL flux.
             */
            
            if (s_R <= 0)
            {
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    F_x_intercell_HLL[ei] = F_x_R[ei];
                }
            }
            else
            {
                double F_x_L[d_num_eqn];
                for (int si = 0; si < d_num_species; si++)
                {
                    F_x_L[si] = u_L*(Q_L[si].get());
                }
                F_x_L[d_num_species] = u_L*(Q_L[d_num_species].get()) + p_L;
                F_x_L[d_num_species + 1] = u_L*(Q_L[d_num_species + 1].get());
                F_x_L[d_num_species + 2] = u_L*((Q_L[d_num_species + 2].get()) + p_L);
                
                if (s_L >= 0)
                {
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        F_x_intercell_HLL[ei] = F_x_L[ei];
                    }
                }
                else
                {
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        F_x_intercell_HLL[ei] = (s_R*F_x_L[ei] - s_L*F_x_R[ei] +
                            s_R*s_L*(Q_R[ei] - Q_L[ei]))/(s_R - s_L);
                    }
                }
            }
        }
        
        /*
         * Calulate the weights beta for hybridization.
         */
       
        double alpha1, alpha2;
        
        double vel_mag = sqrt(pow(u_R - u_L, 2) + pow(v_R - v_L, 2));
        
        if (vel_mag < EPSILON)
        {
            alpha1 = 1.0;
            alpha2 = 0.0;
        }
        else
        {
            alpha1 = fabs(u_R - u_L)/vel_mag;
            alpha2 = sqrt(1.0 - alpha1*alpha1);
        }
        
        double beta1 = 0.5 + 0.5*alpha1/(alpha1 + alpha2);
        double beta2 = 1.0 - beta1;
        
        /*
         * Compute the HLLC-HLL flux.
         */
        
        for (int si = 0; si < d_num_species; si++)
        {
            F_x_intercell[si].get() = beta1*F_x_intercell_HLLC[si] + beta2*F_x_intercell_HLL[si];
        }
        F_x_intercell[d_num_species].get() = F_x_intercell_HLLC[d_num_species];
        F_x_intercell[d_num_species + 1].get() = beta1*F_x_intercell_HLLC[d_num_species + 1] +
            beta2*F_x_intercell_HLL[d_num_species + 1];
        F_x_intercell[d_num_species + 2].get() = F_x_intercell_HLLC[d_num_species + 2];
    }
    else if (d_dim == tbox::Dimension(3))
    {
        double F_x_intercell_HLLC[d_num_eqn];
        double F_x_intercell_HLL[d_num_eqn];
        
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
        
        const double epsilon_L =
            ((Q_L[d_num_species + d_dim.getValue()].get()) -
                0.5*((Q_L[d_num_species].get())*(Q_L[d_num_species].get()) +
                (Q_L[d_num_species + 1].get())*(Q_L[d_num_species + 1].get()) +
                (Q_L[d_num_species + 2].get())*(Q_L[d_num_species + 2].get()))/
                rho_L)/rho_L;
        
        const double epsilon_R =
            ((Q_R[d_num_species + d_dim.getValue()].get()) -
                0.5*((Q_R[d_num_species].get())*(Q_R[d_num_species].get()) +
                (Q_R[d_num_species + 1].get())*(Q_R[d_num_species + 1].get()) +
                (Q_R[d_num_species + 2].get())*(Q_R[d_num_species + 2].get()))/
                rho_R)/rho_R;
        
        /*
         * Compute the mass fractions.
         */
        double Y_L[d_num_species];
        double Y_R[d_num_species];
        for (int si = 0; si < d_num_species; si++)
        {
            Y_L[si] = Q_L[si]/rho_L;
            Y_R[si] = Q_R[si]/rho_R;
        }
        
        /*
         * Get the pointers to the mass fractions.
         */
        std::vector<const double*> Y_L_ptr;
        std::vector<const double*> Y_R_ptr;
        Y_L_ptr.reserve(d_num_species);
        Y_R_ptr.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            Y_L_ptr.push_back(&Y_L[si]);
            Y_R_ptr.push_back(&Y_R[si]);
        }
        
        const double p_L = d_equation_of_state_mixing_rules->getPressure(
            &rho_L,
            &epsilon_L,
            Y_L_ptr);
        
        const double p_R = d_equation_of_state_mixing_rules->getPressure(
            &rho_R,
            &epsilon_R,
            Y_R_ptr);
        
        const double c_L = d_equation_of_state_mixing_rules->getSoundSpeed(
            &rho_L,
            &p_L,
            Y_L_ptr);
        
        const double c_R = d_equation_of_state_mixing_rules->getSoundSpeed(
            &rho_R,
            &p_R,
            Y_R_ptr);
        
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
            /*
             * Compute the HLLC flux.
             */
            
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
                F_x_intercell_HLLC[ei] = F_x_L[ei] + s_minus*(Q_star_L[ei] - Q_L[ei]);
            }
            
            /*
             * Compute the HLL flux.
             */
            
            if (s_L >= 0)
            {
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    F_x_intercell_HLL[ei] = F_x_L[ei];
                }
            }
            else
            {
                double F_x_R[d_num_eqn];
                for (int si = 0; si < d_num_species; si++)
                {
                    F_x_R[si] = u_R*(Q_R[si].get());
                }
                F_x_R[d_num_species] = u_R*(Q_R[d_num_species].get()) + p_R;
                F_x_R[d_num_species + 1] = u_R*(Q_R[d_num_species + 1].get());
                F_x_R[d_num_species + 2] = u_R*(Q_R[d_num_species + 2].get());
                F_x_R[d_num_species + 3] = u_R*((Q_R[d_num_species + 3].get()) + p_R);
                
                if (s_R <= 0)
                {
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        F_x_intercell_HLL[ei] = F_x_R[ei];
                    }
                }
                else
                {
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        F_x_intercell_HLL[ei] = (s_R*F_x_L[ei] - s_L*F_x_R[ei] +
                            s_R*s_L*(Q_R[ei] - Q_L[ei]))/(s_R - s_L);
                    }
                }
            }
        }
        else
        {
            /*
             * Compute the HLLC flux.
             */
            
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
                F_x_intercell_HLLC[ei] = F_x_R[ei] + s_plus*(Q_star_R[ei] - Q_R[ei]);
            }
            
            /*
             * Compute the HLL flux.
             */
            
            if (s_R <= 0)
            {
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    F_x_intercell_HLL[ei] = F_x_R[ei];
                }
            }
            else
            {
                double F_x_L[d_num_eqn];
                for (int si = 0; si < d_num_species; si++)
                {
                    F_x_L[si] = u_L*(Q_L[si].get());
                }
                F_x_L[d_num_species] = u_L*(Q_L[d_num_species].get()) + p_L;
                F_x_L[d_num_species + 1] = u_L*(Q_L[d_num_species + 1].get());
                F_x_L[d_num_species + 2] = u_L*(Q_L[d_num_species + 2].get());
                F_x_L[d_num_species + 3] = u_L*((Q_L[d_num_species + 3].get()) + p_L);
                
                if (s_L >= 0)
                {
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        F_x_intercell_HLL[ei] = F_x_L[ei];
                    }
                }
                else
                {
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        F_x_intercell_HLL[ei] = (s_R*F_x_L[ei] - s_L*F_x_R[ei] +
                            s_R*s_L*(Q_R[ei] - Q_L[ei]))/(s_R - s_L);
                    }
                }
            }
        }
        
        /*
         * Calulate the weights beta for hybridization.
         */
        
        double alpha1, alpha2;
        
        double vel_mag = sqrt(pow(u_R - u_L, 2) + pow(v_R - v_L, 2) + pow(w_R - w_L, 2));
        
        if (vel_mag < EPSILON)
        {
            alpha1 = 1.0;
            alpha2 = 0.0;
        }
        else
        {
            alpha1 = fabs(u_R - u_L)/vel_mag;
            alpha2 = sqrt(1.0 - alpha1*alpha1);
        }
        
        double beta1 = 0.5 + 0.5*alpha1/(alpha1 + alpha2);
        double beta2 = 1.0 - beta1;
        
        /*
         * Compute the HLLC-HLL flux.
         */
        
        for (int si = 0; si < d_num_species; si++)
        {
            F_x_intercell[si].get() = beta1*F_x_intercell_HLLC[si] + beta2*F_x_intercell_HLL[si];
        }
        F_x_intercell[d_num_species].get() = F_x_intercell_HLLC[d_num_species];
        F_x_intercell[d_num_species + 1].get() = beta1*F_x_intercell_HLLC[d_num_species + 1] +
            beta2*F_x_intercell_HLL[d_num_species + 1];
        F_x_intercell[d_num_species + 2].get() = beta1*F_x_intercell_HLLC[d_num_species + 2] +
            beta2*F_x_intercell_HLL[d_num_species + 2];
        F_x_intercell[d_num_species + 3].get() = F_x_intercell_HLLC[d_num_species + 3];
    }
}


/*
 * Compute the flux in the y-direction at the intercell face
 * from conservative variables.
 */
void
RiemannSolverFourEqnConservativeHLLC_HLL::computeIntercellFluxInYDirectionFromConservativeVariables(
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
            << ": RiemannSolverFourEqnConservativeHLLC_HLL::"
            << "computeIntercellFluxInYDirectionFromConservativeVariables()\n"
            << "There is no y direction for 1D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        double F_y_intercell_HLLC[d_num_eqn];
        double F_y_intercell_HLL[d_num_eqn];
        
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
        
        const double epsilon_B =
            ((Q_B[d_num_species + d_dim.getValue()].get()) -
                0.5*((Q_B[d_num_species].get())*(Q_B[d_num_species].get()) +
                (Q_B[d_num_species + 1].get())*(Q_B[d_num_species + 1].get()))/
                rho_B)/rho_B;
        
        const double epsilon_T =
            ((Q_T[d_num_species + d_dim.getValue()].get()) -
                0.5*((Q_T[d_num_species].get())*(Q_T[d_num_species].get()) +
                (Q_T[d_num_species + 1].get())*(Q_T[d_num_species + 1].get()))/
                rho_T)/rho_T;
        
        /*
         * Compute the mass fractions.
         */
        double Y_B[d_num_species];
        double Y_T[d_num_species];
        for (int si = 0; si < d_num_species; si++)
        {
            Y_B[si] = Q_B[si]/rho_B;
            Y_T[si] = Q_T[si]/rho_T;
        }
        
        /*
         * Get the pointers to the mass fractions.
         */
        std::vector<const double*> Y_B_ptr;
        std::vector<const double*> Y_T_ptr;
        Y_B_ptr.reserve(d_num_species);
        Y_T_ptr.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            Y_B_ptr.push_back(&Y_B[si]);
            Y_T_ptr.push_back(&Y_T[si]);
        }
        
        const double p_B = d_equation_of_state_mixing_rules->getPressure(
            &rho_B,
            &epsilon_B,
            Y_B_ptr);
        
        const double p_T = d_equation_of_state_mixing_rules->getPressure(
            &rho_T,
            &epsilon_T,
            Y_T_ptr);
        
        const double c_B = d_equation_of_state_mixing_rules->getSoundSpeed(
            &rho_B,
            &p_B,
            Y_B_ptr);
        
        const double c_T = d_equation_of_state_mixing_rules->getSoundSpeed(
            &rho_T,
            &p_T,
            Y_T_ptr);
        
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
            /*
             * Compute the HLLC flux.
             */
            
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
                F_y_intercell_HLLC[ei] = F_y_B[ei] + s_minus*(Q_star_B[ei] - Q_B[ei]);
            }
            
            /*
             * Compute the HLL flux.
             */
            
            if (s_B >= 0)
            {
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    F_y_intercell_HLL[ei] = F_y_B[ei];
                }
            }
            else
            {
                double F_y_T[d_num_eqn];
                for (int si = 0; si < d_num_species; si++)
                {
                    F_y_T[si] = v_T*(Q_T[si].get());
                }
                F_y_T[d_num_species] = v_T*(Q_T[d_num_species].get());
                F_y_T[d_num_species + 1] = v_T*(Q_T[d_num_species + 1].get()) + p_T;
                F_y_T[d_num_species + 2] = v_T*((Q_T[d_num_species + 2].get()) + p_T);
                
                if (s_T <= 0)
                {
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        F_y_intercell_HLL[ei] = F_y_T[ei];
                    }
                }
                else
                {
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        F_y_intercell_HLL[ei] = (s_T*F_y_B[ei] - s_B*F_y_T[ei] +
                            s_T*s_B*(Q_T[ei] - Q_B[ei]))/(s_T - s_B);
                    }
                }
            }
        }
        else
        {
            /*
             * Compute the HLLC flux.
             */
            
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
                F_y_intercell_HLLC[ei] = F_y_T[ei] + s_plus*(Q_star_T[ei] - Q_T[ei]);
            }
            
            /*
             * Compute the HLL flux.
             */
            
            if (s_T <= 0)
            {
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    F_y_intercell_HLL[ei] = F_y_T[ei];
                }
            }
            else
            {
                double F_y_B[d_num_eqn];
                for (int si = 0; si < d_num_species; si++)
                {
                    F_y_B[si] = v_B*(Q_B[si].get());
                }
                F_y_B[d_num_species] = v_B*(Q_B[d_num_species].get());
                F_y_B[d_num_species + 1] = v_B*(Q_B[d_num_species + 1].get()) + p_B;
                F_y_B[d_num_species + 2] = v_B*((Q_B[d_num_species + 2].get()) + p_B);
                
                if (s_B >= 0)
                {
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        F_y_intercell_HLL[ei] = F_y_B[ei];
                    }
                }
                else
                {
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        F_y_intercell_HLL[ei] = (s_T*F_y_B[ei] - s_B*F_y_T[ei] +
                            s_T*s_B*(Q_T[ei] - Q_B[ei]))/(s_T - s_B);
                    }
                }
            }
        }
        
        /*
         * Calulate the weights beta for hybridization.
         */
        
        double alpha1, alpha2;
        
        double vel_mag = sqrt(pow(u_T - u_B, 2) + pow(v_T - v_B, 2));
        
        if (vel_mag < EPSILON)
        {
            alpha1 = 1.0;
            alpha2 = 0.0;
        }
        else
        {
            alpha1 = fabs(v_T - v_B)/vel_mag;
            alpha2 = sqrt(1.0 - alpha1*alpha1);
        }
        
        double beta1 = 0.5 + 0.5*alpha1/(alpha1 + alpha2);
        double beta2 = 1.0 - beta1;
        
        /*
         * Compute the HLLC-HLL flux.
         */
        
        for (int si = 0; si < d_num_species; si++)
        {
            F_y_intercell[si].get() = beta1*F_y_intercell_HLLC[si] + beta2*F_y_intercell_HLL[si];
        }
        F_y_intercell[d_num_species].get() = beta1*F_y_intercell_HLLC[d_num_species] +
            beta2*F_y_intercell_HLL[d_num_species];
        F_y_intercell[d_num_species + 1].get() = F_y_intercell_HLLC[d_num_species + 1];
        F_y_intercell[d_num_species + 2].get() = F_y_intercell_HLLC[d_num_species + 2];
    }
    else if (d_dim == tbox::Dimension(3))
    {
        double F_y_intercell_HLLC[d_num_eqn];
        double F_y_intercell_HLL[d_num_eqn];
        
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
        
        const double epsilon_B =
            ((Q_B[d_num_species + d_dim.getValue()].get()) -
                0.5*((Q_B[d_num_species].get())*(Q_B[d_num_species].get()) +
                (Q_B[d_num_species + 1].get())*(Q_B[d_num_species + 1].get()) +
                (Q_B[d_num_species + 2].get())*(Q_B[d_num_species + 2].get()))/
                rho_B)/rho_B;
        
        const double epsilon_T =
            ((Q_T[d_num_species + d_dim.getValue()].get()) -
                0.5*((Q_T[d_num_species].get())*(Q_T[d_num_species].get()) +
                (Q_T[d_num_species + 1].get())*(Q_T[d_num_species + 1].get()) +
                (Q_T[d_num_species + 2].get())*(Q_T[d_num_species + 2].get()))/
                rho_T)/rho_T;
        
        /*
         * Compute the mass fractions.
         */
        double Y_B[d_num_species];
        double Y_T[d_num_species];
        for (int si = 0; si < d_num_species; si++)
        {
            Y_B[si] = Q_B[si]/rho_B;
            Y_T[si] = Q_T[si]/rho_T;
        }
        
        /*
         * Get the pointers to the mass fractions.
         */
        std::vector<const double*> Y_B_ptr;
        std::vector<const double*> Y_T_ptr;
        Y_B_ptr.reserve(d_num_species);
        Y_T_ptr.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            Y_B_ptr.push_back(&Y_B[si]);
            Y_T_ptr.push_back(&Y_T[si]);
        }
        
        const double p_B = d_equation_of_state_mixing_rules->getPressure(
            &rho_B,
            &epsilon_B,
            Y_B_ptr);
        
        const double p_T = d_equation_of_state_mixing_rules->getPressure(
            &rho_T,
            &epsilon_T,
            Y_T_ptr);
        
        const double c_B = d_equation_of_state_mixing_rules->getSoundSpeed(
            &rho_B,
            &p_B,
            Y_B_ptr);
        
        const double c_T = d_equation_of_state_mixing_rules->getSoundSpeed(
            &rho_T,
            &p_T,
            Y_T_ptr);
        
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
            /*
             * Compute the HLLC flux.
             */
            
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
                F_y_intercell_HLLC[ei] = F_y_B[ei] + s_minus*(Q_star_B[ei] - Q_B[ei]);
            }
            
            /*
             * Compute the HLL flux.
             */
            
            if (s_B >= 0)
            {
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    F_y_intercell_HLL[ei] = F_y_B[ei];
                }
            }
            else
            {
                double F_y_T[d_num_eqn];
                for (int si = 0; si < d_num_species; si++)
                {
                    F_y_T[si] = v_T*(Q_T[si].get());
                }
                F_y_T[d_num_species] = v_T*(Q_T[d_num_species].get());
                F_y_T[d_num_species + 1] = v_T*(Q_T[d_num_species + 1].get()) + p_T;
                F_y_T[d_num_species + 2] = v_T*(Q_T[d_num_species + 2].get());
                F_y_T[d_num_species + 3] = v_T*((Q_T[d_num_species + 3].get()) + p_T);
                
                if (s_T <= 0)
                {
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        F_y_intercell_HLL[ei] = F_y_T[ei];
                    }
                }
                else
                {
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        F_y_intercell_HLL[ei] = (s_T*F_y_B[ei] - s_B*F_y_T[ei] +
                            s_T*s_B*(Q_T[ei] - Q_B[ei]))/(s_T - s_B);
                    }
                }
            }
        }
        else
        {
            /*
             * Compute the HLLC flux.
             */
            
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
               F_y_intercell_HLLC[ei] = F_y_T[ei] + s_plus*(Q_star_T[ei] - Q_T[ei]);
            }
            
            /*
             * Compute the HLL flux.
             */
            
            if (s_T <= 0)
            {
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    F_y_intercell_HLL[ei] = F_y_T[ei];
                }
            }
            else
            {
                double F_y_B[d_num_eqn];
                for (int si = 0; si < d_num_species; si++)
                {
                    F_y_B[si] = v_B*(Q_B[si].get());
                }
                F_y_B[d_num_species] = v_B*(Q_B[d_num_species].get());
                F_y_B[d_num_species + 1] = v_B*(Q_B[d_num_species + 1].get()) + p_B;
                F_y_B[d_num_species + 2] = v_B*(Q_B[d_num_species + 2].get());
                F_y_B[d_num_species + 3] = v_B*((Q_B[d_num_species + 3].get()) + p_B);
                
                if (s_B >= 0)
                {
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        F_y_intercell_HLL[ei] = F_y_B[ei];
                    }
                }
                else
                {
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        F_y_intercell_HLL[ei] = (s_T*F_y_B[ei] - s_B*F_y_T[ei] +
                            s_T*s_B*(Q_T[ei] - Q_B[ei]))/(s_T - s_B);
                    }
                }
            }
        }
        
        /*
         * Calulate the weights beta for hybridization.
         */
        
        double alpha1, alpha2;
        
        double vel_mag = sqrt(pow(u_T - u_B, 2) + pow(v_T - v_B, 2) + pow(w_T - w_B, 2));
        
        if (vel_mag < EPSILON)
        {
           alpha1 = 1.0;
           alpha2 = 0.0;
        }
        else
        {
           alpha1 = fabs(v_T - v_B)/vel_mag;
           alpha2 = sqrt(1.0 - alpha1*alpha1);
        }
        
        double beta1 = 0.5 + 0.5*alpha1/(alpha1 + alpha2);
        double beta2 = 1.0 - beta1;
        
        /*
         * Compute the HLLC-HLL flux.
         */
        
        for (int si = 0; si < d_num_species; si++)
        {
            F_y_intercell[si].get() = beta1*F_y_intercell_HLLC[si] + beta2*F_y_intercell_HLL[si];
        }
        F_y_intercell[d_num_species].get() = beta1*F_y_intercell_HLLC[d_num_species] +
            beta2*F_y_intercell_HLL[d_num_species];
        F_y_intercell[d_num_species + 1].get() = F_y_intercell_HLLC[d_num_species + 1];
        F_y_intercell[d_num_species + 2].get() = beta1*F_y_intercell_HLLC[d_num_species + 2] +
            beta2*F_y_intercell_HLL[d_num_species + 2];
        F_y_intercell[d_num_species + 3].get() = F_y_intercell_HLLC[d_num_species + 3];
    }
}


/*
 * Compute the flux in the z-direction at the intercell face
 * from conservative variables.
 */
void
RiemannSolverFourEqnConservativeHLLC_HLL::computeIntercellFluxInZDirectionFromConservativeVariables(
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
            << ": RiemannSolverFourEqnConservativeHLLC_HLL::"
            << "computeIntercellFluxInZDirectionFromConservativeVariables()\n"
            << "There is no z direction for 1D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": RiemannSolverFourEqnConservativeHLLC_HLL::"
            << "computeIntercellFluxInZDirectionFromConservativeVariables()\n"
            << "There is no z direction for 2D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        double F_z_intercell_HLLC[d_num_eqn];
        double F_z_intercell_HLL[d_num_eqn];
        
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
        
        const double epsilon_B =
            ((Q_B[d_num_species + d_dim.getValue()].get()) -
                0.5*((Q_B[d_num_species].get())*(Q_B[d_num_species].get()) +
                (Q_B[d_num_species + 1].get())*(Q_B[d_num_species + 1].get()) +
                (Q_B[d_num_species + 2].get())*(Q_B[d_num_species + 2].get()))/
                rho_B)/rho_B;
        
        const double epsilon_F =
            ((Q_F[d_num_species + d_dim.getValue()].get()) -
                0.5*((Q_F[d_num_species].get())*(Q_F[d_num_species].get()) +
                (Q_F[d_num_species + 1].get())*(Q_F[d_num_species + 1].get()) +
                (Q_F[d_num_species + 2].get())*(Q_F[d_num_species + 2].get()))/
                rho_F)/rho_F;
        
        /*
         * Compute the mass fractions.
         */
        double Y_B[d_num_species];
        double Y_F[d_num_species];
        for (int si = 0; si < d_num_species; si++)
        {
            Y_B[si] = Q_B[si]/rho_B;
            Y_F[si] = Q_F[si]/rho_F;
        }
        
        /*
         * Get the pointers to the mass fractions.
         */
        std::vector<const double*> Y_B_ptr;
        std::vector<const double*> Y_F_ptr;
        Y_B_ptr.reserve(d_num_species);
        Y_F_ptr.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            Y_B_ptr.push_back(&Y_B[si]);
            Y_F_ptr.push_back(&Y_F[si]);
        }
        
        const double p_B = d_equation_of_state_mixing_rules->getPressure(
            &rho_B,
            &epsilon_B,
            Y_B_ptr);
        
        const double p_F = d_equation_of_state_mixing_rules->getPressure(
            &rho_F,
            &epsilon_F,
            Y_F_ptr);
        
        const double c_B = d_equation_of_state_mixing_rules->getSoundSpeed(
            &rho_B,
            &p_B,
            Y_B_ptr);
        
        const double c_F = d_equation_of_state_mixing_rules->getSoundSpeed(
            &rho_F,
            &p_F,
            Y_F_ptr);
        
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
            /*
             * Compute the HLLC flux.
             */
            
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
                F_z_intercell_HLLC[ei] = F_z_B[ei] + s_minus*(Q_star_B[ei] - Q_B[ei]);
            }
            
            /*
             * Compute the HLL flux.
             */
            
            if (s_B >= 0)
            {
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    F_z_intercell_HLL[ei] = F_z_B[ei];
                }
            }
            else
            {
                double F_z_F[d_num_eqn];
                for (int si = 0; si < d_num_species; si++)
                {
                    F_z_F[si] = w_F*(Q_F[si].get());
                }
                F_z_F[d_num_species] = w_F*(Q_F[d_num_species].get());
                F_z_F[d_num_species + 1] = w_F*(Q_F[d_num_species + 1].get());
                F_z_F[d_num_species + 2] = w_F*(Q_F[d_num_species + 2].get()) + p_F;
                F_z_F[d_num_species + 3] = w_F*((Q_F[d_num_species + 3].get()) + p_F);
                
                if (s_F <= 0)
                {
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        F_z_intercell_HLL[ei] = F_z_F[ei];
                    }
                }
                else
                {
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        F_z_intercell_HLL[ei] = (s_F*F_z_B[ei] - s_B*F_z_F[ei] +
                            s_F*s_B*(Q_F[ei] - Q_B[ei]))/(s_F - s_B);
                    }
                }
            }
        }
        else
        {
            /*
             * Compute the HLLC flux.
             */
            
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
                F_z_intercell_HLLC[ei] = F_z_F[ei] + s_plus*(Q_star_F[ei] - Q_F[ei]);
            }
            
            /*
             * Compute the HLL flux.
             */
            
            if (s_F <= 0)
            {
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    F_z_intercell_HLL[ei] = F_z_F[ei];
                }
            }
            else
            {
                double F_z_B[d_num_eqn];
                for (int si = 0; si < d_num_species; si++)
                {
                    F_z_B[si] = w_B*(Q_B[si].get());
                }
                F_z_B[d_num_species] = w_B*(Q_B[d_num_species].get());
                F_z_B[d_num_species + 1] = w_B*(Q_B[d_num_species + 1].get());
                F_z_B[d_num_species + 2] = w_B*(Q_B[d_num_species + 2].get()) + p_B;
                F_z_B[d_num_species + 3] = w_B*((Q_B[d_num_species + 3].get()) + p_B);
                
                if (s_B >= 0)
                {
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        F_z_intercell_HLL[ei] = F_z_B[ei];
                    }
                }
                else
                {
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        F_z_intercell_HLL[ei] = (s_F*F_z_B[ei] - s_B*F_z_F[ei] +
                            s_F*s_B*(Q_F[ei] - Q_B[ei]))/(s_F - s_B);
                    }
                }
            }
        }
        
        /*
         * Calulate the weights beta for hybridization.
         */
        
        double alpha1, alpha2;
        
        double vel_mag = sqrt(pow(u_F - u_B, 2) + pow(v_F - v_B, 2) + pow(w_F - w_B, 2));
        
        if (vel_mag < EPSILON)
        {
            alpha1 = 1.0;
            alpha2 = 0.0;
        }
        else
        {
            alpha1 = fabs(w_F - w_B)/vel_mag;
            alpha2 = sqrt(1.0 - alpha1*alpha1);
        }
        
        double beta1 = 0.5 + 0.5*alpha1/(alpha1 + alpha2);
        double beta2 = 1.0 - beta1;
        
        /*
         * Compute the HLLC-HLL flux.
         */
        
        for (int si = 0; si < d_num_species; si++)
        {
            F_z_intercell[si].get() = beta1*F_z_intercell_HLLC[si] + beta2*F_z_intercell_HLL[si];
        }
        F_z_intercell[d_num_species].get() = beta1*F_z_intercell_HLLC[d_num_species] +
            beta2*F_z_intercell_HLL[d_num_species];
        F_z_intercell[d_num_species + 1].get() = beta1*F_z_intercell_HLLC[d_num_species + 1] +
            beta2*F_z_intercell_HLL[d_num_species + 1];
        F_z_intercell[d_num_species + 2].get() = F_z_intercell_HLLC[d_num_species + 2];
        F_z_intercell[d_num_species + 3].get() = F_z_intercell_HLLC[d_num_species + 3];
    }
}


/*
 * Compute the flux in the x-direction at the intercell face
 * from primitive variables.
 */
void
RiemannSolverFourEqnConservativeHLLC_HLL::computeIntercellFluxInXDirectionFromPrimitiveVariables(
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
        
        const double rho_L = d_equation_of_state_mixing_rules->getMixtureDensity(
            rho_Y_L);
        
        const double rho_R = d_equation_of_state_mixing_rules->getMixtureDensity(
            rho_Y_R);
        
        /*
         * Compute the mass fractions.
         */
        double Y_L[d_num_species];
        double Y_R[d_num_species];
        for (int si = 0; si < d_num_species; si++)
        {
            Y_L[si] = V_L[si]/rho_L;
            Y_R[si] = V_R[si]/rho_R;
        }
        
        /*
         * Get the pointers to the mass fractions.
         */
        std::vector<const double*> Y_L_ptr;
        std::vector<const double*> Y_R_ptr;
        Y_L_ptr.reserve(d_num_species);
        Y_R_ptr.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            Y_L_ptr.push_back(&Y_L[si]);
            Y_R_ptr.push_back(&Y_R[si]);
        }
        
        const double c_L = d_equation_of_state_mixing_rules->getSoundSpeed(
            &rho_L,
            &(V_L[d_num_species + d_dim.getValue()].get()),
            Y_L_ptr);
        
        const double c_R = d_equation_of_state_mixing_rules->getSoundSpeed(
            &rho_R,
            &(V_R[d_num_species + d_dim.getValue()].get()),
            Y_R_ptr);
        
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
            
            const double epsilon_L = d_equation_of_state_mixing_rules->getInternalEnergy(
                &rho_L,
                &(V_L[d_num_species + d_dim.getValue()].get()),
                Y_L_ptr);
            
            double Q_L[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                Q_L[si] = (V_L[si].get());
            }
            Q_L[d_num_species] = rho_L*(V_L[d_num_species].get());
            Q_L[d_num_species + d_dim.getValue()] = rho_L*(epsilon_L +
                0.5*(V_L[d_num_species].get())*(V_L[d_num_species].get()));
            
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
            
            const double epsilon_R = d_equation_of_state_mixing_rules->getInternalEnergy(
                &rho_R,
                &(V_R[d_num_species + d_dim.getValue()].get()),
                Y_R_ptr);
            
            double Q_R[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                Q_R[si] = (V_R[si].get());
            }
            Q_R[d_num_species] = rho_R*(V_R[d_num_species].get());
            Q_R[d_num_species + d_dim.getValue()] = rho_R*(epsilon_R +
                0.5*(V_R[d_num_species].get())*(V_R[d_num_species].get()));
            
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
        double F_x_intercell_HLLC[d_num_eqn];
        double F_x_intercell_HLL[d_num_eqn];
        
        std::vector<const double*> rho_Y_L;
        std::vector<const double*> rho_Y_R;
        rho_Y_L.reserve(d_num_species);
        rho_Y_R.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y_L.push_back(&(V_L[si].get()));
            rho_Y_R.push_back(&(V_R[si].get()));
        }
        
        const double rho_L = d_equation_of_state_mixing_rules->getMixtureDensity(
            rho_Y_L);
        
        const double rho_R = d_equation_of_state_mixing_rules->getMixtureDensity(
            rho_Y_R);
        
        /*
         * Compute the mass fractions.
         */
        double Y_L[d_num_species];
        double Y_R[d_num_species];
        for (int si = 0; si < d_num_species; si++)
        {
            Y_L[si] = V_L[si]/rho_L;
            Y_R[si] = V_R[si]/rho_R;
        }
        
        /*
         * Get the pointers to the mass fractions.
         */
        std::vector<const double*> Y_L_ptr;
        std::vector<const double*> Y_R_ptr;
        Y_L_ptr.reserve(d_num_species);
        Y_R_ptr.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            Y_L_ptr.push_back(&Y_L[si]);
            Y_R_ptr.push_back(&Y_R[si]);
        }
        
        const double c_L = d_equation_of_state_mixing_rules->getSoundSpeed(
            &rho_L,
            &(V_L[d_num_species + d_dim.getValue()].get()),
            Y_L_ptr);
        
        const double c_R = d_equation_of_state_mixing_rules->getSoundSpeed(
            &rho_R,
            &(V_R[d_num_species + d_dim.getValue()].get()),
            Y_R_ptr);
        
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
            /*
             * Compute the HLLC flux.
             */
            
            const double Chi_star_L = (s_L - (V_L[d_num_species].get()))/(s_L - s_star);
            
            const double epsilon_L = d_equation_of_state_mixing_rules->getInternalEnergy(
                &rho_L,
                &(V_L[d_num_species + d_dim.getValue()].get()),
                Y_L_ptr);
            
            double Q_L[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                Q_L[si] = (V_L[si].get());
            }
            Q_L[d_num_species] = rho_L*(V_L[d_num_species].get());
            Q_L[d_num_species + 1] = rho_L*(V_L[d_num_species + 1].get());
            Q_L[d_num_species + d_dim.getValue()] = rho_L*(epsilon_L +
                0.5*((V_L[d_num_species].get())*(V_L[d_num_species].get()) +
                (V_L[d_num_species + 1].get())*(V_L[d_num_species + 1].get())));
            
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
                F_x_intercell_HLLC[ei] = F_x_L[ei] + s_minus*(Q_star_L[ei] - Q_L[ei]);
            }
            
            /*
             * Compute the HLL flux.
             */
            
            if (s_L >= 0)
            {
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    F_x_intercell_HLL[ei] = F_x_L[ei];
                }
            }
            else
            {
                const double epsilon_R = d_equation_of_state_mixing_rules->getInternalEnergy(
                    &rho_R,
                    &(V_R[d_num_species + d_dim.getValue()].get()),
                    Y_R_ptr);
                
                const double E_R = rho_R*(epsilon_R +
                    0.5*((V_R[d_num_species].get())*(V_R[d_num_species].get()) +
                    (V_R[d_num_species + 1].get())*(V_R[d_num_species + 1].get())));
                
                if (s_R <= 0)
                {
                    for (int si = 0; si < d_num_species; si++)
                    {
                        F_x_intercell_HLL[si] = (V_R[d_num_species].get())*(V_R[si].get());
                    }
                    F_x_intercell_HLL[d_num_species] = (V_R[d_num_species].get())*rho_R*(V_R[d_num_species].get()) +
                        (V_R[d_num_species + d_dim.getValue()].get());
                    F_x_intercell_HLL[d_num_species + 1] = (V_R[d_num_species].get())*rho_R*(V_R[d_num_species + 1].get());
                    F_x_intercell_HLL[d_num_species + d_dim.getValue()] = (V_R[d_num_species].get())*(E_R +
                        (V_R[d_num_species + d_dim.getValue()].get()));
                }
                else
                {
                    double Q_R[d_num_eqn];
                    for (int si = 0; si < d_num_species; si++)
                    {
                        Q_R[si] = (V_R[si].get());
                    }
                    Q_R[d_num_species] = rho_R*(V_R[d_num_species].get());
                    Q_R[d_num_species + 1] = rho_R*(V_R[d_num_species + 1].get());
                    Q_R[d_num_species + d_dim.getValue()] = E_R;
                    
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
                        F_x_intercell_HLL[ei] = (s_R*F_x_L[ei] - s_L*F_x_R[ei] +
                            s_R*s_L*(Q_R[ei] - Q_L[ei]))/(s_R - s_L);
                    }
                }
            }
        }
        else
        {
            /*
             * Compute the HLLC flux.
             */
            
            const double Chi_star_R = (s_R - (V_R[d_num_species].get()))/(s_R - s_star);
            
            const double epsilon_R = d_equation_of_state_mixing_rules->getInternalEnergy(
                &rho_R,
                &(V_R[d_num_species + d_dim.getValue()].get()),
                Y_R_ptr);
            
            double Q_R[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                Q_R[si] = (V_R[si].get());
            }
            Q_R[d_num_species] = rho_R*(V_R[d_num_species].get());
            Q_R[d_num_species + 1] = rho_R*(V_R[d_num_species + 1].get());
            Q_R[d_num_species + d_dim.getValue()] = rho_R*(epsilon_R +
                0.5*((V_R[d_num_species].get())*(V_R[d_num_species].get()) +
                (V_R[d_num_species + 1].get())*(V_R[d_num_species + 1].get())));
            
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
                F_x_intercell_HLLC[ei] = F_x_R[ei] + s_plus*(Q_star_R[ei] - Q_R[ei]);
            }
            
            /*
             * Compute the HLL flux.
             */
            
            if (s_R <= 0)
            {
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    F_x_intercell_HLL[ei] = F_x_R[ei];
                }
            }
            else
            {
                const double epsilon_L = d_equation_of_state_mixing_rules->getInternalEnergy(
                    &rho_L,
                    &(V_L[d_num_species + d_dim.getValue()].get()),
                    Y_L_ptr);
                
                const double E_L = rho_L*(epsilon_L +
                    0.5*((V_L[d_num_species].get())*(V_L[d_num_species].get()) +
                    (V_L[d_num_species + 1].get())*(V_L[d_num_species + 1].get())));
                
                if (s_L >= 0)
                {
                    for (int si = 0; si < d_num_species; si++)
                    {
                        F_x_intercell_HLL[si] = (V_L[d_num_species].get())*(V_L[si].get());
                    }
                    F_x_intercell_HLL[d_num_species] = (V_L[d_num_species].get())*rho_L*(V_L[d_num_species].get()) +
                        (V_L[d_num_species + d_dim.getValue()].get());
                    F_x_intercell_HLL[d_num_species + 1] = (V_L[d_num_species].get())*rho_L*(V_L[d_num_species + 1].get());
                    F_x_intercell_HLL[d_num_species + d_dim.getValue()] = (V_L[d_num_species].get())*(E_L +
                        (V_L[d_num_species + d_dim.getValue()].get()));
                }
                else
                {
                    double Q_L[d_num_eqn];
                    for (int si = 0; si < d_num_species; si++)
                    {
                        Q_L[si] = (V_L[si].get());
                    }
                    Q_L[d_num_species] = rho_L*(V_L[d_num_species].get());
                    Q_L[d_num_species + 1] = rho_L*(V_L[d_num_species + 1].get());
                    Q_L[d_num_species + d_dim.getValue()] = E_L;
                    
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
                        F_x_intercell_HLL[ei] = (s_R*F_x_L[ei] - s_L*F_x_R[ei] +
                            s_R*s_L*(Q_R[ei] - Q_L[ei]))/(s_R - s_L);
                    }
                }
            }
        }
        
        /*
         * Calulate the weights beta for hybridization.
         */
       
        double alpha1, alpha2;
        
        double vel_mag = sqrt(pow((V_R[d_num_species].get()) - (V_L[d_num_species].get()), 2) +
            pow((V_R[d_num_species + 1].get()) - (V_L[d_num_species + 1].get()), 2));
        
        if (vel_mag < EPSILON)
        {
            alpha1 = 1.0;
            alpha2 = 0.0;
        }
        else
        {
            alpha1 = fabs((V_R[d_num_species].get()) - (V_L[d_num_species].get()))/vel_mag;
            alpha2 = sqrt(1.0 - alpha1*alpha1);
        }
        
        double beta1 = 0.5 + 0.5*alpha1/(alpha1 + alpha2);
        double beta2 = 1.0 - beta1;
        
        /*
         * Compute the HLLC-HLL flux.
         */
        
        for (int si = 0; si < d_num_species; si++)
        {
            F_x_intercell[si].get() = beta1*F_x_intercell_HLLC[si] + beta2*F_x_intercell_HLL[si];
        }
        F_x_intercell[d_num_species].get() = F_x_intercell_HLLC[d_num_species];
        F_x_intercell[d_num_species + 1].get() = beta1*F_x_intercell_HLLC[d_num_species + 1] +
            beta2*F_x_intercell_HLL[d_num_species + 1];
        F_x_intercell[d_num_species + 2].get() = F_x_intercell_HLLC[d_num_species + 2];
    }
    else if (d_dim == tbox::Dimension(3))
    {
        double F_x_intercell_HLLC[d_num_eqn];
        double F_x_intercell_HLL[d_num_eqn];
        
        std::vector<const double*> rho_Y_L;
        std::vector<const double*> rho_Y_R;
        rho_Y_L.reserve(d_num_species);
        rho_Y_R.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y_L.push_back(&(V_L[si].get()));
            rho_Y_R.push_back(&(V_R[si].get()));
        }
        
        const double rho_L = d_equation_of_state_mixing_rules->getMixtureDensity(
            rho_Y_L);
        
        const double rho_R = d_equation_of_state_mixing_rules->getMixtureDensity(
            rho_Y_R);
        
        /*
         * Compute the mass fractions.
         */
        double Y_L[d_num_species];
        double Y_R[d_num_species];
        for (int si = 0; si < d_num_species; si++)
        {
            Y_L[si] = V_L[si]/rho_L;
            Y_R[si] = V_R[si]/rho_R;
        }
        
        /*
         * Get the pointers to the mass fractions.
         */
        std::vector<const double*> Y_L_ptr;
        std::vector<const double*> Y_R_ptr;
        Y_L_ptr.reserve(d_num_species);
        Y_R_ptr.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            Y_L_ptr.push_back(&Y_L[si]);
            Y_R_ptr.push_back(&Y_R[si]);
        }
        
        const double c_L = d_equation_of_state_mixing_rules->getSoundSpeed(
            &rho_L,
            &(V_L[d_num_species + d_dim.getValue()].get()),
            Y_L_ptr);
        
        const double c_R = d_equation_of_state_mixing_rules->getSoundSpeed(
            &rho_R,
            &(V_R[d_num_species + d_dim.getValue()].get()),
            Y_R_ptr);
        
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
            /*
             * Compute the HLLC flux.
             */
            
            const double Chi_star_L = (s_L - (V_L[d_num_species].get()))/(s_L - s_star);
            
            const double epsilon_L = d_equation_of_state_mixing_rules->getInternalEnergy(
                &rho_L,
                &(V_L[d_num_species + d_dim.getValue()].get()),
                Y_L_ptr);
            
            double Q_L[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                Q_L[si] = (V_L[si].get());
            }
            Q_L[d_num_species] = rho_L*(V_L[d_num_species].get());
            Q_L[d_num_species + 1] = rho_L*(V_L[d_num_species + 1].get());
            Q_L[d_num_species + 2] = rho_L*(V_L[d_num_species + 2].get());
            Q_L[d_num_species + d_dim.getValue()] = rho_L*(epsilon_L +
                0.5*((V_L[d_num_species].get())*(V_L[d_num_species].get()) +
                (V_L[d_num_species + 1].get())*(V_L[d_num_species + 1].get()) +
                (V_L[d_num_species + 2].get())*(V_L[d_num_species + 2].get())));
            
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
                F_x_intercell_HLLC[ei] = F_x_L[ei] + s_minus*(Q_star_L[ei] - Q_L[ei]);
            }
            
            /*
             * Compute the HLL flux.
             */
            
            if (s_L >= 0)
            {
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    F_x_intercell_HLL[ei] = F_x_L[ei];
                }
            }
            else
            {
                const double epsilon_R = d_equation_of_state_mixing_rules->getInternalEnergy(
                    &rho_R,
                    &(V_R[d_num_species + d_dim.getValue()].get()),
                    Y_R_ptr);
                
                const double E_R = rho_R*(epsilon_R +
                    0.5*((V_R[d_num_species].get())*(V_R[d_num_species].get()) +
                    (V_R[d_num_species + 1].get())*(V_R[d_num_species + 1].get()) +
                    (V_R[d_num_species + 2].get())*(V_R[d_num_species + 2].get())));
                
                if (s_R <= 0)
                {
                    for (int si = 0; si < d_num_species; si++)
                    {
                        F_x_intercell_HLL[si] = (V_R[d_num_species].get())*(V_R[si].get());
                    }
                    F_x_intercell_HLL[d_num_species] = (V_R[d_num_species].get())*rho_R*(V_R[d_num_species].get()) +
                        (V_R[d_num_species + d_dim.getValue()].get());
                    F_x_intercell_HLL[d_num_species + 1] = (V_R[d_num_species].get())*rho_R*(V_R[d_num_species + 1].get());
                    F_x_intercell_HLL[d_num_species + 2] = (V_R[d_num_species].get())*rho_R*(V_R[d_num_species + 2].get());
                    F_x_intercell_HLL[d_num_species + d_dim.getValue()] = (V_R[d_num_species].get())*(E_R +
                        (V_R[d_num_species + d_dim.getValue()].get()));
                }
                else
                {
                    double Q_R[d_num_eqn];
                    for (int si = 0; si < d_num_species; si++)
                    {
                        Q_R[si] = (V_R[si].get());
                    }
                    Q_R[d_num_species] = rho_R*(V_R[d_num_species].get());
                    Q_R[d_num_species + 1] = rho_R*(V_R[d_num_species + 1].get());
                    Q_R[d_num_species + 2] = rho_R*(V_R[d_num_species + 2].get());
                    Q_R[d_num_species + d_dim.getValue()] = E_R;
                    
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
                        F_x_intercell_HLL[ei] = (s_R*F_x_L[ei] - s_L*F_x_R[ei] +
                            s_R*s_L*(Q_R[ei] - Q_L[ei]))/(s_R - s_L);
                    }
                }
            }
        }
        else
        {
            /*
             * Compute the HLLC flux.
             */
            
            const double Chi_star_R = (s_R - (V_R[d_num_species].get()))/(s_R - s_star);
            
            const double epsilon_R = d_equation_of_state_mixing_rules->getInternalEnergy(
                &rho_R,
                &(V_R[d_num_species + d_dim.getValue()].get()),
                Y_R_ptr);
            
            double Q_R[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                Q_R[si] = (V_R[si].get());
            }
            Q_R[d_num_species] = rho_R*(V_R[d_num_species].get());
            Q_R[d_num_species + 1] = rho_R*(V_R[d_num_species + 1].get());
            Q_R[d_num_species + 2] = rho_R*(V_R[d_num_species + 2].get());
            Q_R[d_num_species + d_dim.getValue()] = rho_R*(epsilon_R +
                0.5*((V_R[d_num_species].get())*(V_R[d_num_species].get()) +
                (V_R[d_num_species + 1].get())*(V_R[d_num_species + 1].get()) +
                (V_R[d_num_species + 2].get())*(V_R[d_num_species + 2].get())));
            
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
                F_x_intercell_HLLC[ei] = F_x_R[ei] + s_plus*(Q_star_R[ei] - Q_R[ei]);
            }
            
            /*
             * Compute the HLL flux.
             */
            
            if (s_R <= 0)
            {
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    F_x_intercell_HLL[ei] = F_x_R[ei];
                }
            }
            else
            {
                const double epsilon_L = d_equation_of_state_mixing_rules->getInternalEnergy(
                    &rho_L,
                    &(V_L[d_num_species + d_dim.getValue()].get()),
                    Y_L_ptr);
                
                const double E_L = rho_L*(epsilon_L +
                    0.5*((V_L[d_num_species].get())*(V_L[d_num_species].get()) +
                    (V_L[d_num_species + 1].get())*(V_L[d_num_species + 1].get()) +
                    (V_L[d_num_species + 2].get())*(V_L[d_num_species + 2].get())));
                
                if (s_L >= 0)
                {
                    for (int si = 0; si < d_num_species; si++)
                    {
                        F_x_intercell_HLL[si] = (V_L[d_num_species].get())*(V_L[si].get());
                    }
                    F_x_intercell_HLL[d_num_species] = (V_L[d_num_species].get())*rho_L*(V_L[d_num_species].get()) +
                        (V_L[d_num_species + d_dim.getValue()].get());
                    F_x_intercell_HLL[d_num_species + 1] = (V_L[d_num_species].get())*rho_L*(V_L[d_num_species + 1].get());
                    F_x_intercell_HLL[d_num_species + 2] = (V_L[d_num_species].get())*rho_L*(V_L[d_num_species + 2].get());
                    F_x_intercell_HLL[d_num_species + d_dim.getValue()] = (V_L[d_num_species].get())*(E_L +
                        (V_L[d_num_species + d_dim.getValue()].get()));
                }
                else
                {
                    double Q_L[d_num_eqn];
                    for (int si = 0; si < d_num_species; si++)
                    {
                        Q_L[si] = (V_L[si].get());
                    }
                    Q_L[d_num_species] = rho_L*(V_L[d_num_species].get());
                    Q_L[d_num_species + 1] = rho_L*(V_L[d_num_species + 1].get());
                    Q_L[d_num_species + 2] = rho_L*(V_L[d_num_species + 2].get());
                    Q_L[d_num_species + d_dim.getValue()] = E_L;
                    
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
                        F_x_intercell_HLL[ei] = (s_R*F_x_L[ei] - s_L*F_x_R[ei] +
                            s_R*s_L*(Q_R[ei] - Q_L[ei]))/(s_R - s_L);
                    }
                }
            }
        }
        
        /*
         * Calulate the weights beta for hybridization.
         */
        
        double alpha1, alpha2;
        
        double vel_mag = sqrt(pow((V_R[d_num_species].get()) - (V_L[d_num_species].get()), 2) +
            pow((V_R[d_num_species + 1].get()) - (V_L[d_num_species + 1].get()), 2) +
            pow((V_R[d_num_species + 2].get()) - (V_L[d_num_species + 2].get()), 2));
        
        if (vel_mag < EPSILON)
        {
            alpha1 = 1.0;
            alpha2 = 0.0;
        }
        else
        {
            alpha1 = fabs((V_R[d_num_species].get()) - (V_L[d_num_species].get()))/vel_mag;
            alpha2 = sqrt(1.0 - alpha1*alpha1);
        }
        
        double beta1 = 0.5 + 0.5*alpha1/(alpha1 + alpha2);
        double beta2 = 1.0 - beta1;
        
        /*
         * Compute the HLLC-HLL flux.
         */
        
        for (int si = 0; si < d_num_species; si++)
        {
            F_x_intercell[si].get() = beta1*F_x_intercell_HLLC[si] + beta2*F_x_intercell_HLL[si];
        }
        F_x_intercell[d_num_species].get() = F_x_intercell_HLLC[d_num_species];
        F_x_intercell[d_num_species + 1].get() = beta1*F_x_intercell_HLLC[d_num_species + 1] +
            beta2*F_x_intercell_HLL[d_num_species + 1];
        F_x_intercell[d_num_species + 2].get() = beta1*F_x_intercell_HLLC[d_num_species + 2] +
            beta2*F_x_intercell_HLL[d_num_species + 2];
        F_x_intercell[d_num_species + 3].get() = F_x_intercell_HLLC[d_num_species + 3];
    }
}


/*
 * Compute the flux in the y-direction at the intercell face
 * from primitive variables.
 */
void
RiemannSolverFourEqnConservativeHLLC_HLL::computeIntercellFluxInYDirectionFromPrimitiveVariables(
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
            << ": RiemannSolverFourEqnConservativeHLLC_HLL::"
            << "computeIntercellFluxInYDirectionFromPrimitiveVariables()\n"
            << "There is no y direction for 1D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        double F_y_intercell_HLLC[d_num_eqn];
        double F_y_intercell_HLL[d_num_eqn];
        
        std::vector<const double*> rho_Y_B;
        std::vector<const double*> rho_Y_T;
        rho_Y_B.reserve(d_num_species);
        rho_Y_T.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y_B.push_back(&(V_B[si].get()));
            rho_Y_T.push_back(&(V_T[si].get()));
        }
        
        const double rho_B = d_equation_of_state_mixing_rules->getMixtureDensity(
            rho_Y_B);
        
        const double rho_T = d_equation_of_state_mixing_rules->getMixtureDensity(
            rho_Y_T);
        
        /*
         * Compute the mass fractions.
         */
        double Y_B[d_num_species];
        double Y_T[d_num_species];
        for (int si = 0; si < d_num_species; si++)
        {
            Y_B[si] = V_B[si]/rho_B;
            Y_T[si] = V_T[si]/rho_T;
        }
        
        /*
         * Get the pointers to the mass fractions.
         */
        std::vector<const double*> Y_B_ptr;
        std::vector<const double*> Y_T_ptr;
        Y_B_ptr.reserve(d_num_species);
        Y_T_ptr.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            Y_B_ptr.push_back(&Y_B[si]);
            Y_T_ptr.push_back(&Y_T[si]);
        }
        
        const double c_B = d_equation_of_state_mixing_rules->getSoundSpeed(
            &rho_B,
            &(V_B[d_num_species + d_dim.getValue()].get()),
            Y_B_ptr);
        
        const double c_T = d_equation_of_state_mixing_rules->getSoundSpeed(
            &rho_T,
            &(V_T[d_num_species + d_dim.getValue()].get()),
            Y_T_ptr);
        
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
            /*
             * Compute the HLLC flux.
             */
            
            const double Chi_star_B = (s_B - (V_B[d_num_species + 1].get()))/(s_B - s_star);
            
            const double epsilon_B = d_equation_of_state_mixing_rules->getInternalEnergy(
                &rho_B,
                &(V_B[d_num_species + d_dim.getValue()].get()),
                Y_B_ptr);
            
            double Q_B[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                Q_B[si] = (V_B[si].get());
            }
            Q_B[d_num_species] = rho_B*(V_B[d_num_species].get());
            Q_B[d_num_species + 1] = rho_B*(V_B[d_num_species + 1].get());
            Q_B[d_num_species + d_dim.getValue()] = rho_B*(epsilon_B +
                0.5*((V_B[d_num_species].get())*(V_B[d_num_species].get()) +
                (V_B[d_num_species + 1].get())*(V_B[d_num_species + 1].get())));
            
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
                F_y_intercell_HLLC[ei] = F_y_B[ei] + s_minus*(Q_star_B[ei] - Q_B[ei]);
            }
            
            /*
             * Compute the HLL flux.
             */
            
            if (s_B >= 0)
            {
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    F_y_intercell_HLL[ei] = F_y_B[ei];
                }
            }
            else
            {
                const double epsilon_T = d_equation_of_state_mixing_rules->getInternalEnergy(
                    &rho_T,
                    &(V_T[d_num_species + d_dim.getValue()].get()),
                    Y_T_ptr);
                
                const double E_T = rho_T*(epsilon_T +
                    0.5*((V_T[d_num_species].get())*(V_T[d_num_species].get()) +
                    (V_T[d_num_species + 1].get())*(V_T[d_num_species + 1].get())));
                
                if (s_T <= 0)
                {
                    for (int si = 0; si < d_num_species; si++)
                    {
                        F_y_intercell_HLL[si] = (V_T[d_num_species + 1].get())*(V_T[si].get());
                    }
                    F_y_intercell_HLL[d_num_species] = (V_T[d_num_species + 1].get())*rho_T*(V_T[d_num_species].get());
                    F_y_intercell_HLL[d_num_species + 1] = (V_T[d_num_species + 1].get())*rho_T*(V_T[d_num_species + 1].get()) +
                        (V_T[d_num_species + d_dim.getValue()].get());
                    F_y_intercell_HLL[d_num_species + d_dim.getValue()] = (V_T[d_num_species + 1].get())*(E_T +
                        (V_T[d_num_species + d_dim.getValue()].get()));
                }
                else
                {
                    double Q_T[d_num_eqn];
                    for (int si = 0; si < d_num_species; si++)
                    {
                        Q_T[si] = (V_T[si].get());
                    }
                    Q_T[d_num_species] = rho_T*(V_T[d_num_species].get());
                    Q_T[d_num_species + 1] = rho_T*(V_T[d_num_species + 1].get());
                    Q_T[d_num_species + d_dim.getValue()] = E_T;
                    
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
                        F_y_intercell_HLL[ei] = (s_T*F_y_B[ei] - s_B*F_y_T[ei] +
                            s_T*s_B*(Q_T[ei] - Q_B[ei]))/(s_T - s_B);
                    }
                }
            }
        }
        else
        {
            /*
             * Compute the HLLC flux.
             */
            
            const double Chi_star_T = (s_T - (V_T[d_num_species + 1].get()))/(s_T - s_star);
            
            const double epsilon_T = d_equation_of_state_mixing_rules->getInternalEnergy(
                &rho_T,
                &(V_T[d_num_species + d_dim.getValue()].get()),
                Y_T_ptr);
            
            double Q_T[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                Q_T[si] = (V_T[si].get());
            }
            Q_T[d_num_species] = rho_T*(V_T[d_num_species].get());
            Q_T[d_num_species + 1] = rho_T*(V_T[d_num_species + 1].get());
            Q_T[d_num_species + d_dim.getValue()] = rho_T*(epsilon_T +
                0.5*((V_T[d_num_species].get())*(V_T[d_num_species].get()) +
                (V_T[d_num_species + 1].get())*(V_T[d_num_species + 1].get())));
            
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
                F_y_intercell_HLLC[ei] = F_y_T[ei] + s_plus*(Q_star_T[ei] - Q_T[ei]);
            }
            
            /*
             * Compute the HLL flux.
             */
            
            if (s_T <= 0)
            {
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    F_y_intercell_HLL[ei] = F_y_T[ei];
                }
            }
            else
            {
                const double epsilon_B = d_equation_of_state_mixing_rules->getInternalEnergy(
                    &rho_B,
                    &(V_B[d_num_species + d_dim.getValue()].get()),
                    Y_B_ptr);
                
                const double E_B = rho_B*(epsilon_B +
                    0.5*((V_B[d_num_species].get())*(V_B[d_num_species].get()) +
                    (V_B[d_num_species + 1].get())*(V_B[d_num_species + 1].get())));
                
                if (s_B >= 0)
                {
                    for (int si = 0; si < d_num_species; si++)
                    {
                        F_y_intercell_HLL[si] = (V_B[d_num_species + 1].get())*(V_B[si].get());
                    }
                    F_y_intercell_HLL[d_num_species] = (V_B[d_num_species + 1].get())*rho_B*(V_B[d_num_species].get());
                    F_y_intercell_HLL[d_num_species + 1] = (V_B[d_num_species + 1].get())*rho_B*(V_B[d_num_species + 1].get()) +
                        (V_B[d_num_species + d_dim.getValue()].get());
                    F_y_intercell_HLL[d_num_species + d_dim.getValue()] = (V_B[d_num_species + 1].get())*(E_B +
                        (V_B[d_num_species + d_dim.getValue()].get()));
                }
                else
                {
                    double Q_B[d_num_eqn];
                    for (int si = 0; si < d_num_species; si++)
                    {
                        Q_B[si] = (V_B[si].get());
                    }
                    Q_B[d_num_species] = rho_B*(V_B[d_num_species].get());
                    Q_B[d_num_species + 1] = rho_B*(V_B[d_num_species + 1].get());
                    Q_B[d_num_species + d_dim.getValue()] = E_B;
                    
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
                        F_y_intercell_HLL[ei] = (s_T*F_y_B[ei] - s_B*F_y_T[ei] +
                            s_T*s_B*(Q_T[ei] - Q_B[ei]))/(s_T - s_B);
                    }
                }
            }
        }
        
        /*
         * Calulate the weights beta for hybridization.
         */
        
        double alpha1, alpha2;
        
        double vel_mag = sqrt(pow((V_T[d_num_species].get()) - (V_B[d_num_species].get()), 2) +
            pow((V_T[d_num_species + 1].get()) - (V_B[d_num_species + 1].get()), 2));
        
        if (vel_mag < EPSILON)
        {
            alpha1 = 1.0;
            alpha2 = 0.0;
        }
        else
        {
            alpha1 = fabs((V_T[d_num_species + 1].get()) - (V_B[d_num_species + 1].get()))/vel_mag;
            alpha2 = sqrt(1.0 - alpha1*alpha1);
        }
        
        double beta1 = 0.5 + 0.5*alpha1/(alpha1 + alpha2);
        double beta2 = 1.0 - beta1;
        
        /*
         * Compute the HLLC-HLL flux.
         */
        
        for (int si = 0; si < d_num_species; si++)
        {
            F_y_intercell[si].get() = beta1*F_y_intercell_HLLC[si] + beta2*F_y_intercell_HLL[si];
        }
        F_y_intercell[d_num_species].get() = beta1*F_y_intercell_HLLC[d_num_species] +
            beta2*F_y_intercell_HLL[d_num_species];
        F_y_intercell[d_num_species + 1].get() = F_y_intercell_HLLC[d_num_species + 1];
        F_y_intercell[d_num_species + 2].get() = F_y_intercell_HLLC[d_num_species + 2];
    }
    else if (d_dim == tbox::Dimension(3))
    {
        double F_y_intercell_HLLC[d_num_eqn];
        double F_y_intercell_HLL[d_num_eqn];
        
        std::vector<const double*> rho_Y_B;
        std::vector<const double*> rho_Y_T;
        rho_Y_B.reserve(d_num_species);
        rho_Y_T.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y_B.push_back(&(V_B[si].get()));
            rho_Y_T.push_back(&(V_T[si].get()));
        }
        
        const double rho_B = d_equation_of_state_mixing_rules->getMixtureDensity(
            rho_Y_B);
        
        const double rho_T = d_equation_of_state_mixing_rules->getMixtureDensity(
            rho_Y_T);
        
        /*
         * Compute the mass fractions.
         */
        double Y_B[d_num_species];
        double Y_T[d_num_species];
        for (int si = 0; si < d_num_species; si++)
        {
            Y_B[si] = V_B[si]/rho_B;
            Y_T[si] = V_T[si]/rho_T;
        }
        
        /*
         * Get the pointers to the mass fractions.
         */
        std::vector<const double*> Y_B_ptr;
        std::vector<const double*> Y_T_ptr;
        Y_B_ptr.reserve(d_num_species);
        Y_T_ptr.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            Y_B_ptr.push_back(&Y_B[si]);
            Y_T_ptr.push_back(&Y_T[si]);
        }
        
        const double c_B = d_equation_of_state_mixing_rules->getSoundSpeed(
            &rho_B,
            &(V_B[d_num_species + d_dim.getValue()].get()),
            Y_B_ptr);
        
        const double c_T = d_equation_of_state_mixing_rules->getSoundSpeed(
            &rho_T,
            &(V_T[d_num_species + d_dim.getValue()].get()),
            Y_T_ptr);
        
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
            /*
             * Compute the HLLC flux.
             */
            
            const double Chi_star_B = (s_B - (V_B[d_num_species + 1].get()))/(s_B - s_star);
            
            const double epsilon_B = d_equation_of_state_mixing_rules->getInternalEnergy(
                &rho_B,
                &(V_B[d_num_species + d_dim.getValue()].get()),
                Y_B_ptr);
            
            double Q_B[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                Q_B[si] = (V_B[si].get());
            }
            Q_B[d_num_species] = rho_B*(V_B[d_num_species].get());
            Q_B[d_num_species + 1] = rho_B*(V_B[d_num_species + 1].get());
            Q_B[d_num_species + 2] = rho_B*(V_B[d_num_species + 2].get());
            Q_B[d_num_species + d_dim.getValue()] = rho_B*(epsilon_B +
                0.5*((V_B[d_num_species].get())*(V_B[d_num_species].get()) +
                (V_B[d_num_species + 1].get())*(V_B[d_num_species + 1].get()) +
                (V_B[d_num_species + 2].get())*(V_B[d_num_species + 2].get())));
            
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
                F_y_intercell_HLLC[ei] = F_y_B[ei] + s_minus*(Q_star_B[ei] - Q_B[ei]);
            }
            
            /*
             * Compute the HLL flux.
             */
            
            if (s_B >= 0)
            {
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    F_y_intercell_HLL[ei] = F_y_B[ei];
                }
            }
            else
            {
                const double epsilon_T = d_equation_of_state_mixing_rules->getInternalEnergy(
                    &rho_T,
                    &(V_T[d_num_species + d_dim.getValue()].get()),
                    Y_T_ptr);
                
                const double E_T = rho_T*(epsilon_T +
                    0.5*((V_T[d_num_species].get())*(V_T[d_num_species].get()) +
                    (V_T[d_num_species + 1].get())*(V_T[d_num_species + 1].get()) +
                    (V_T[d_num_species + 2].get())*(V_T[d_num_species + 2].get())));
                
                if (s_T <= 0)
                {
                    for (int si = 0; si < d_num_species; si++)
                    {
                        F_y_intercell_HLL[si] = (V_T[d_num_species + 1].get())*(V_T[si].get());
                    }
                    F_y_intercell_HLL[d_num_species] = (V_T[d_num_species + 1].get())*rho_T*(V_T[d_num_species].get());
                    F_y_intercell_HLL[d_num_species + 1] = (V_T[d_num_species + 1].get())*rho_T*(V_T[d_num_species + 1].get()) +
                        (V_T[d_num_species + d_dim.getValue()].get());
                    F_y_intercell_HLL[d_num_species + 2] = (V_T[d_num_species + 1].get())*rho_T*(V_T[d_num_species + 2].get());
                    F_y_intercell_HLL[d_num_species + d_dim.getValue()] = (V_T[d_num_species + 1].get())*(E_T +
                        (V_T[d_num_species + d_dim.getValue()].get()));
                }
                else
                {
                    double Q_T[d_num_eqn];
                    for (int si = 0; si < d_num_species; si++)
                    {
                        Q_T[si] = (V_T[si].get());
                    }
                    Q_T[d_num_species] = rho_T*(V_T[d_num_species].get());
                    Q_T[d_num_species + 1] = rho_T*(V_T[d_num_species + 1].get());
                    Q_T[d_num_species + 2] = rho_T*(V_T[d_num_species + 2].get());
                    Q_T[d_num_species + d_dim.getValue()] = E_T;
                    
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
                        F_y_intercell_HLL[ei] = (s_T*F_y_B[ei] - s_B*F_y_T[ei] +
                            s_T*s_B*(Q_T[ei] - Q_B[ei]))/(s_T - s_B);
                    }
                }
            }
        }
        else
        {
            /*
             * Compute the HLLC flux.
             */
            
            const double Chi_star_T = (s_T - (V_T[d_num_species + 1].get()))/(s_T - s_star);
            
            const double epsilon_T = d_equation_of_state_mixing_rules->getInternalEnergy(
                &rho_T,
                &(V_T[d_num_species + d_dim.getValue()].get()),
                Y_T_ptr);
            
            double Q_T[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                Q_T[si] = (V_T[si].get());
            }
            Q_T[d_num_species] = rho_T*(V_T[d_num_species].get());
            Q_T[d_num_species + 1] = rho_T*(V_T[d_num_species + 1].get());
            Q_T[d_num_species + 2] = rho_T*(V_T[d_num_species + 2].get());
            Q_T[d_num_species + d_dim.getValue()] = rho_T*(V_T[d_num_species + 2].get());
            Q_T[d_num_species + d_dim.getValue()] = rho_T*(epsilon_T +
                0.5*((V_T[d_num_species].get())*(V_T[d_num_species].get()) +
                (V_T[d_num_species + 1].get())*(V_T[d_num_species + 1].get()) +
                (V_T[d_num_species + 2].get())*(V_T[d_num_species + 2].get())));
            
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
                F_y_intercell_HLLC[ei] = F_y_T[ei] + s_plus*(Q_star_T[ei] - Q_T[ei]);
            }
            
            /*
             * Compute the HLL flux.
             */
            
            if (s_T <= 0)
            {
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    F_y_intercell_HLL[ei] = F_y_T[ei];
                }
            }
            else
            {
                const double epsilon_B = d_equation_of_state_mixing_rules->getInternalEnergy(
                    &rho_B,
                    &(V_B[d_num_species + d_dim.getValue()].get()),
                    Y_B_ptr);
                
                const double E_B = rho_B*(epsilon_B +
                    0.5*((V_B[d_num_species].get())*(V_B[d_num_species].get()) +
                    (V_B[d_num_species + 1].get())*(V_B[d_num_species + 1].get()) +
                    (V_B[d_num_species + 2].get())*(V_B[d_num_species + 2].get())));
                
                if (s_B >= 0)
                {
                    for (int si = 0; si < d_num_species; si++)
                    {
                        F_y_intercell_HLL[si] = (V_B[d_num_species + 1].get())*(V_B[si].get());
                    }
                    F_y_intercell_HLL[d_num_species] = (V_B[d_num_species + 1].get())*rho_B*(V_B[d_num_species].get());
                    F_y_intercell_HLL[d_num_species + 1] = (V_B[d_num_species + 1].get())*rho_B*(V_B[d_num_species + 1].get()) +
                        (V_B[d_num_species + d_dim.getValue()].get());
                    F_y_intercell_HLL[d_num_species + 2] = (V_B[d_num_species + 1].get())*rho_B*(V_B[d_num_species + 2].get());
                    F_y_intercell_HLL[d_num_species + d_dim.getValue()] = (V_B[d_num_species + 1].get())*(E_B +
                        (V_B[d_num_species + d_dim.getValue()].get()));
                }
                else
                {
                    double Q_B[d_num_eqn];
                    for (int si = 0; si < d_num_species; si++)
                    {
                        Q_B[si] = (V_B[si].get());
                    }
                    Q_B[d_num_species] = rho_B*(V_B[d_num_species].get());
                    Q_B[d_num_species + 1] = rho_B*(V_B[d_num_species + 1].get());
                    Q_B[d_num_species + 2] = rho_B*(V_B[d_num_species + 2].get());
                    Q_B[d_num_species + d_dim.getValue()] = E_B;
                    
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
                        F_y_intercell_HLL[ei] = (s_T*F_y_B[ei] - s_B*F_y_T[ei] +
                            s_T*s_B*(Q_T[ei] - Q_B[ei]))/(s_T - s_B);
                    }
                }
            }
        }
        
        /*
         * Calulate the weights beta for hybridization.
         */
        
        double alpha1, alpha2;
        
        double vel_mag = sqrt(pow((V_T[d_num_species].get()) - (V_B[d_num_species].get()), 2) +
            pow((V_T[d_num_species + 1].get()) - (V_B[d_num_species + 1].get()), 2) +
            pow((V_T[d_num_species + 2].get()) - (V_B[d_num_species + 2].get()), 2));
        
        if (vel_mag < EPSILON)
        {
            alpha1 = 1.0;
            alpha2 = 0.0;
        }
        else
        {
            alpha1 = fabs((V_T[d_num_species + 1].get()) - (V_B[d_num_species + 1].get()))/vel_mag;
            alpha2 = sqrt(1.0 - alpha1*alpha1);
        }
        
        double beta1 = 0.5 + 0.5*alpha1/(alpha1 + alpha2);
        double beta2 = 1.0 - beta1;
        
        /*
         * Compute the HLLC-HLL flux.
         */
        
        for (int si = 0; si < d_num_species; si++)
        {
            F_y_intercell[si].get() = beta1*F_y_intercell_HLLC[si] + beta2*F_y_intercell_HLL[si];
        }
        F_y_intercell[d_num_species].get() = beta1*F_y_intercell_HLLC[d_num_species] +
            beta2*F_y_intercell_HLL[d_num_species];
        F_y_intercell[d_num_species + 1].get() = F_y_intercell_HLLC[d_num_species + 1];
        F_y_intercell[d_num_species + 2].get() = beta1*F_y_intercell_HLLC[d_num_species + 2] +
            beta2*F_y_intercell_HLL[d_num_species + 2];
        F_y_intercell[d_num_species + 3].get() = F_y_intercell_HLLC[d_num_species + 3];
    }
}


/*
 * Compute the flux in the z-direction at the intercell face
 * from primitive variables.
 */
void
RiemannSolverFourEqnConservativeHLLC_HLL::computeIntercellFluxInZDirectionFromPrimitiveVariables(
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
            << ": RiemannSolverFourEqnConservativeHLLC_HLL::"
            << "computeIntercellFluxInZDirectionFromPrimitiveVariables()\n"
            << "There is no z direction for 1D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": RiemannSolverFourEqnConservativeHLLC_HLL::"
            << "computeIntercellFluxInZDirectionFromPrimitiveVariables()\n"
            << "There is no z direction for 2D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        double F_z_intercell_HLLC[d_num_eqn];
        double F_z_intercell_HLL[d_num_eqn];
        
        std::vector<const double*> rho_Y_B;
        std::vector<const double*> rho_Y_F;
        rho_Y_B.reserve(d_num_species);
        rho_Y_F.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y_B.push_back(&(V_B[si].get()));
            rho_Y_F.push_back(&(V_F[si].get()));
        }
        
        const double rho_B = d_equation_of_state_mixing_rules->getMixtureDensity(
            rho_Y_B);
        
        const double rho_F = d_equation_of_state_mixing_rules->getMixtureDensity(
            rho_Y_F);
        
        /*
         * Compute the mass fractions.
         */
        double Y_B[d_num_species];
        double Y_F[d_num_species];
        for (int si = 0; si < d_num_species; si++)
        {
            Y_B[si] = V_B[si]/rho_B;
            Y_F[si] = V_F[si]/rho_F;
        }
        
        /*
         * Get the pointers to the mass fractions.
         */
        std::vector<const double*> Y_B_ptr;
        std::vector<const double*> Y_F_ptr;
        Y_B_ptr.reserve(d_num_species);
        Y_F_ptr.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            Y_B_ptr.push_back(&Y_B[si]);
            Y_F_ptr.push_back(&Y_F[si]);
        }
        
        const double c_B = d_equation_of_state_mixing_rules->getSoundSpeed(
            &rho_B,
            &(V_B[d_num_species + d_dim.getValue()].get()),
            Y_B_ptr);
        
        const double c_F = d_equation_of_state_mixing_rules->getSoundSpeed(
            &rho_F,
            &(V_F[d_num_species + d_dim.getValue()].get()),
            Y_F_ptr);
        
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
            /*
             * Compute the HLLC flux.
             */
            
            const double Chi_star_B = (s_B - (V_B[d_num_species + 2].get()))/(s_B - s_star);
            
            const double epsilon_B = d_equation_of_state_mixing_rules->getInternalEnergy(
                &rho_B,
                &(V_B[d_num_species + d_dim.getValue()].get()),
                Y_B_ptr);
            
            double Q_B[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                Q_B[si] = (V_B[si].get());
            }
            Q_B[d_num_species] = rho_B*(V_B[d_num_species].get());
            Q_B[d_num_species + 1] = rho_B*(V_B[d_num_species + 1].get());
            Q_B[d_num_species + 2] = rho_B*(V_B[d_num_species + 2].get());
            Q_B[d_num_species + d_dim.getValue()] = rho_B*(epsilon_B +
                0.5*((V_B[d_num_species].get())*(V_B[d_num_species].get()) +
                (V_B[d_num_species + 1].get())*(V_B[d_num_species + 1].get()) +
                (V_B[d_num_species + 2].get())*(V_B[d_num_species + 2].get())));
            
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
                F_z_intercell_HLLC[ei] = F_z_B[ei] + s_minus*(Q_star_B[ei] - Q_B[ei]);
            }
            
            /*
             * Compute the HLL flux.
             */
            
            if (s_B >= 0)
            {
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    F_z_intercell_HLL[ei] = F_z_B[ei];
                }
            }
            else
            {
                const double epsilon_F = d_equation_of_state_mixing_rules->getInternalEnergy(
                    &rho_F,
                    &(V_F[d_num_species + d_dim.getValue()].get()),
                    Y_F_ptr);
                
                const double E_F = rho_F*(epsilon_F +
                0.5*((V_F[d_num_species].get())*(V_F[d_num_species].get()) +
                (V_F[d_num_species + 1].get())*(V_F[d_num_species + 1].get()) +
                (V_F[d_num_species + 2].get())*(V_F[d_num_species + 2].get())));
                
                if (s_F <= 0)
                {
                    for (int si = 0; si < d_num_species; si++)
                    {
                        F_z_intercell_HLL[si] = (V_F[d_num_species + 2].get())*(V_F[si].get());
                    }
                    F_z_intercell_HLL[d_num_species] = (V_F[d_num_species + 2].get())*rho_F*(V_F[d_num_species].get());
                    F_z_intercell_HLL[d_num_species + 1] = (V_F[d_num_species + 2].get())*rho_F*(V_F[d_num_species + 1].get());
                    F_z_intercell_HLL[d_num_species + 2] = (V_F[d_num_species + 2].get())*rho_F*(V_F[d_num_species + 2].get()) +
                        (V_F[d_num_species + d_dim.getValue()].get());
                    F_z_intercell_HLL[d_num_species + d_dim.getValue()] = (V_F[d_num_species + 2].get())*(E_F +
                        (V_F[d_num_species + d_dim.getValue()].get()));
                }
                else
                {
                    double Q_F[d_num_eqn];
                    for (int si = 0; si < d_num_species; si++)
                    {
                        Q_F[si] = (V_F[si].get());
                    }
                    Q_F[d_num_species] = rho_F*(V_F[d_num_species].get());
                    Q_F[d_num_species + 1] = rho_F*(V_F[d_num_species + 1].get());
                    Q_F[d_num_species + 2] = rho_F*(V_F[d_num_species + 2].get());
                    Q_F[d_num_species + d_dim.getValue()] = E_F;
                    
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
                        F_z_intercell_HLL[ei] = (s_F*F_z_B[ei] - s_B*F_z_F[ei] +
                            s_F*s_B*(Q_F[ei] - Q_B[ei]))/(s_F - s_B);
                    }
                }
            }
        }
        else
        {
            /*
             * Compute the HLLC flux.
             */
            
            const double Chi_star_F = (s_F - (V_F[d_num_species + 2].get()))/(s_F - s_star);
            
            const double epsilon_F = d_equation_of_state_mixing_rules->getInternalEnergy(
                &rho_F,
                &(V_F[d_num_species + d_dim.getValue()].get()),
                Y_F_ptr);
            
            double Q_F[d_num_eqn];
            for (int si = 0; si < d_num_species; si++)
            {
                Q_F[si] = (V_F[si].get());
            }
            Q_F[d_num_species] = rho_F*(V_F[d_num_species].get());
            Q_F[d_num_species + 1] = rho_F*(V_F[d_num_species + 1].get());
            Q_F[d_num_species + 2] = rho_F*(V_F[d_num_species + 2].get());
            Q_F[d_num_species + d_dim.getValue()] = rho_F*(epsilon_F +
                0.5*((V_F[d_num_species].get())*(V_F[d_num_species].get()) +
                (V_F[d_num_species + 1].get())*(V_F[d_num_species + 1].get()) +
                (V_F[d_num_species + 2].get())*(V_F[d_num_species + 2].get())));
            
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
                F_z_intercell_HLLC[ei] = F_z_F[ei] + s_plus*(Q_star_F[ei] - Q_F[ei]);
            }
            
            /*
             * Compute the HLL flux.
             */
            
            if (s_F <= 0)
            {
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    F_z_intercell_HLL[ei] = F_z_F[ei];
                }
            }
            else
            {
                const double epsilon_B = d_equation_of_state_mixing_rules->getInternalEnergy(
                    &rho_B,
                    &(V_B[d_num_species + d_dim.getValue()].get()),
                    Y_B_ptr);
                
                const double E_B = rho_B*(epsilon_B +
                    0.5*((V_B[d_num_species].get())*(V_B[d_num_species].get()) +
                    (V_B[d_num_species + 1].get())*(V_B[d_num_species + 1].get()) +
                    (V_B[d_num_species + 2].get())*(V_B[d_num_species + 2].get())));
                
                if (s_B >= 0)
                {
                    for (int si = 0; si < d_num_species; si++)
                    {
                        F_z_intercell_HLL[si] = (V_B[d_num_species + 2].get())*(V_B[si].get());
                    }
                    F_z_intercell_HLL[d_num_species] = (V_B[d_num_species + 2].get())*rho_B*(V_B[d_num_species].get());
                    F_z_intercell_HLL[d_num_species + 1] = (V_B[d_num_species + 2].get())*rho_B*(V_B[d_num_species + 1].get());
                    F_z_intercell_HLL[d_num_species + 2] = (V_B[d_num_species + 2].get())*rho_B*(V_B[d_num_species + 2].get()) +
                        (V_B[d_num_species + d_dim.getValue()].get());
                    F_z_intercell_HLL[d_num_species + d_dim.getValue()] = (V_B[d_num_species + 2].get())*(E_B +
                        (V_B[d_num_species + d_dim.getValue()].get()));
                }
                else
                {
                    double Q_B[d_num_eqn];
                    for (int si = 0; si < d_num_species; si++)
                    {
                        Q_B[si] = (V_B[si].get());
                    }
                    Q_B[d_num_species] = rho_B*(V_B[d_num_species].get());
                    Q_B[d_num_species + 1] = rho_B*(V_B[d_num_species + 1].get());
                    Q_B[d_num_species + 2] = rho_B*(V_B[d_num_species + 2].get());
                    Q_B[d_num_species + d_dim.getValue()] = E_B;
                    
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
                        F_z_intercell_HLL[ei] = (s_F*F_z_B[ei] - s_B*F_z_F[ei] +
                            s_F*s_B*(Q_F[ei] - Q_B[ei]))/(s_F - s_B);
                    }
                }
            }
        }
        
        /*
         * Calulate the weights beta for hybridization.
         */
        
        double alpha1, alpha2;
        
        double vel_mag = sqrt(pow((V_F[d_num_species].get()) - (V_B[d_num_species].get()), 2) +
            pow((V_F[d_num_species + 1].get()) - (V_B[d_num_species + 1].get()), 2) +
            pow((V_F[d_num_species + 2].get()) - (V_B[d_num_species + 2].get()), 2));
        
        if (vel_mag < EPSILON)
        {
            alpha1 = 1.0;
            alpha2 = 0.0;
        }
        else
        {
            alpha1 = fabs((V_F[d_num_species + 2].get()) - (V_B[d_num_species + 2].get()))/vel_mag;
            alpha2 = sqrt(1.0 - alpha1*alpha1);
        }
        
        double beta1 = 0.5 + 0.5*alpha1/(alpha1 + alpha2);
        double beta2 = 1.0 - beta1;
        
        /*
         * Compute the HLLC-HLL flux.
         */
        
        for (int si = 0; si < d_num_species; si++)
        {
            F_z_intercell[si].get() = beta1*F_z_intercell_HLLC[si] + beta2*F_z_intercell_HLL[si];
        }
        F_z_intercell[d_num_species].get() = beta1*F_z_intercell_HLLC[d_num_species] +
            beta2*F_z_intercell_HLL[d_num_species];
        F_z_intercell[d_num_species + 1].get() = beta1*F_z_intercell_HLLC[d_num_species + 1] +
            beta2*F_z_intercell_HLL[d_num_species + 1];
        F_z_intercell[d_num_species + 2].get() = F_z_intercell_HLLC[d_num_species + 2];
        F_z_intercell[d_num_species + 3].get() = F_z_intercell_HLLC[d_num_species + 3];
    }
}
