#include "flow/flow_models/Riemann_solvers/RiemannSolverSingleSpeciesHLLC_HLL.hpp"

#define EPSILON 1e-40

/*
 * Compute the flux at the intercell face from conservative variables.
 */
void
RiemannSolverSingleSpeciesHLLC_HLL::computeIntercellFluxFromConservativeVariables(
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
                << ": RiemannSolverSingleSpeciesHLLC_HLL::"
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
RiemannSolverSingleSpeciesHLLC_HLL::computeIntercellFluxFromPrimitiveVariables(
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
                << ": RiemannSolverSingleSpeciesHLLC_HLL::"
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
RiemannSolverSingleSpeciesHLLC_HLL::computeIntercellFluxInXDirectionFromConservativeVariables(
    std::vector<boost::reference_wrapper<double> >& F_x_intercell,
    const std::vector<boost::reference_wrapper<double> >& Q_L,
    const std::vector<boost::reference_wrapper<double> >& Q_R)
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(F_x_intercell.size()) == d_num_eqn);
    TBOX_ASSERT(static_cast<int>(Q_L.size()) == d_num_eqn);
    TBOX_ASSERT(static_cast<int>(Q_R.size()) == d_num_eqn);
#endif
    
    // Get the thermodynamic properties of the species.
    std::vector<const double*> thermo_properties_ptr;
    thermo_properties_ptr.reserve(static_cast<int> (d_thermo_properties.size()));
    for (int ti = 0; ti < static_cast<int> (d_thermo_properties.size()); ti++)
    {
        thermo_properties_ptr.push_back(&d_thermo_properties[ti]);
    }
    
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
            &(Q_L[2].get()),
            thermo_properties_ptr);
        
        const double p_R = d_equation_of_state->getPressure(
            &(Q_R[0].get()),
            m_R,
            &(Q_R[2].get()),
            thermo_properties_ptr);
        
        const double c_L = d_equation_of_state->getSoundSpeed(
            &(Q_L[0].get()),
            &p_L,
            thermo_properties_ptr);
        
        const double c_R = d_equation_of_state->getSoundSpeed(
            &(Q_R[0].get()),
            &p_R,
            thermo_properties_ptr);
        
        const double u_average = 0.5*(u_L + u_R);
        const double c_average = 0.5*(c_L + c_R);
        
        const double s_L = fmin(u_average - c_average, u_L - c_L);
        const double s_R = fmax(u_average + c_average, u_R + c_R);
        
        const double s_minus = fmin(0.0, s_L);
        const double s_plus  = fmax(0.0, s_R);
        
        const double s_star =
            (p_R - p_L + (Q_L[1].get())*(s_L - u_L) - (Q_R[1].get())*(s_R - u_R))/
                ((Q_L[0].get())*(s_L - u_L) - (Q_R[0].get())*(s_R - u_R));
        
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
        double F_x_intercell_HLLC[d_num_eqn];
        double F_x_intercell_HLL[d_num_eqn];
        
        const double u_L = (Q_L[1].get())/(Q_L[0].get());
        const double u_R = (Q_R[1].get())/(Q_R[0].get());
        
        const double v_L = (Q_L[2].get())/(Q_L[0].get());
        const double v_R = (Q_R[2].get())/(Q_R[0].get());
        
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
            &(Q_L[3].get()),
            thermo_properties_ptr);
        
        const double p_R = d_equation_of_state->getPressure(
            &(Q_R[0].get()),
            m_R,
            &(Q_R[3].get()),
            thermo_properties_ptr);
        
        const double c_L = d_equation_of_state->getSoundSpeed(
            &(Q_L[0].get()),
            &p_L,
            thermo_properties_ptr);
        
        const double c_R = d_equation_of_state->getSoundSpeed(
            &(Q_R[0].get()),
            &p_R,
            thermo_properties_ptr);
        
        const double u_average = 0.5*(u_L + u_R);
        const double c_average = 0.5*(c_L + c_R);
        
        const double s_L = fmin(u_average - c_average, u_L - c_L);
        const double s_R = fmax(u_average + c_average, u_R + c_R);
        
        const double s_minus = fmin(0.0, s_L);
        const double s_plus  = fmax(0.0, s_R);
        
        const double s_star =
            (p_R - p_L + (Q_L[1].get())*(s_L - u_L) - (Q_R[1].get())*(s_R - u_R))/
                ((Q_L[0].get())*(s_L - u_L) - (Q_R[0].get())*(s_R - u_R));
        
        if (s_star > 0)
        {
            /*
             * Compute the HLLC flux.
             */
            
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
                F_x_R[0] = (Q_R[1].get());
                F_x_R[1] = u_R*(Q_R[1].get()) + p_R;
                F_x_R[2] = u_R*(Q_R[2].get());
                F_x_R[3] = u_R*((Q_R[3].get()) + p_R);
                
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
                F_x_L[0] = (Q_L[1].get());
                F_x_L[1] = u_L*(Q_L[1].get()) + p_L;
                F_x_L[2] = u_L*(Q_L[2].get());
                F_x_L[3] = u_L*((Q_L[3].get()) + p_L);
                
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
        
        F_x_intercell[0].get() = beta1*F_x_intercell_HLLC[0] + beta2*F_x_intercell_HLL[0];
        F_x_intercell[1].get() = F_x_intercell_HLLC[1];
        F_x_intercell[2].get() = beta1*F_x_intercell_HLLC[2] + beta2*F_x_intercell_HLL[2];
        F_x_intercell[3].get() = F_x_intercell_HLLC[3];
    }
    else if (d_dim == tbox::Dimension(3))
    {
        double F_x_intercell_HLLC[d_num_eqn];
        double F_x_intercell_HLL[d_num_eqn];
        
        const double u_L = (Q_L[1].get())/(Q_L[0].get());
        const double u_R = (Q_R[1].get())/(Q_R[0].get());
        
        const double v_L = (Q_L[2].get())/(Q_L[0].get());
        const double v_R = (Q_R[2].get())/(Q_R[0].get());
        
        const double w_L = (Q_L[3].get())/(Q_L[0].get());
        const double w_R = (Q_R[3].get())/(Q_R[0].get());
        
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
            &(Q_L[4].get()),
            thermo_properties_ptr);
        
        const double p_R = d_equation_of_state->getPressure(
            &(Q_R[0].get()),
            m_R,
            &(Q_R[4].get()),
            thermo_properties_ptr);
        
        const double c_L = d_equation_of_state->getSoundSpeed(
            &(Q_L[0].get()),
            &p_L,
            thermo_properties_ptr);
        
        const double c_R = d_equation_of_state->getSoundSpeed(
            &(Q_R[0].get()),
            &p_R,
            thermo_properties_ptr);
        
        const double u_average = 0.5*(u_L + u_R);
        const double c_average = 0.5*(c_L + c_R);
        
        const double s_L = fmin(u_average - c_average, u_L - c_L);
        const double s_R = fmax(u_average + c_average, u_R + c_R);
        
        const double s_minus = fmin(0.0, s_L);
        const double s_plus  = fmax(0.0, s_R);
        
        const double s_star =
            (p_R - p_L + (Q_L[1].get())*(s_L - u_L) - (Q_R[1].get())*(s_R - u_R))/
                ((Q_L[0].get())*(s_L - u_L) - (Q_R[0].get())*(s_R - u_R));
        
        if (s_star > 0)
        {
            /*
             * Compute the HLLC flux.
             */
            
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
            F_x_L[4] = u_L*((Q_L[4].get()) + p_L);
            
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
                F_x_R[0] = (Q_R[1].get());
                F_x_R[1] = u_R*(Q_R[1].get()) + p_R;
                F_x_R[2] = u_R*(Q_R[2].get());
                F_x_R[3] = u_R*(Q_R[3].get());
                F_x_R[4] = u_R*((Q_R[4].get()) + p_R);
                
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
            F_x_R[4] = u_R*((Q_R[4].get()) + p_R);
            
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
                F_x_L[0] = (Q_L[1].get());
                F_x_L[1] = u_L*(Q_L[1].get()) + p_L;
                F_x_L[2] = u_L*(Q_L[2].get());
                F_x_L[3] = u_L*(Q_L[3].get());
                F_x_L[4] = u_L*((Q_L[4].get()) + p_L);
                
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
        
        F_x_intercell[0].get() = beta1*F_x_intercell_HLLC[0] + beta2*F_x_intercell_HLL[0];
        F_x_intercell[1].get() = F_x_intercell_HLLC[1];
        F_x_intercell[2].get() = beta1*F_x_intercell_HLLC[2] + beta2*F_x_intercell_HLL[2];
        F_x_intercell[3].get() = beta1*F_x_intercell_HLLC[3] + beta2*F_x_intercell_HLL[3];
        F_x_intercell[4].get() = F_x_intercell_HLLC[4];
    }
}


/*
 * Compute the flux in the y-direction at the intercell face
 * from conservative variables.
 */
void
RiemannSolverSingleSpeciesHLLC_HLL::computeIntercellFluxInYDirectionFromConservativeVariables(
    std::vector<boost::reference_wrapper<double> >& F_y_intercell,
    const std::vector<boost::reference_wrapper<double> >& Q_B,
    const std::vector<boost::reference_wrapper<double> >& Q_T)
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(F_y_intercell.size()) == d_num_eqn);
    TBOX_ASSERT(static_cast<int>(Q_B.size()) == d_num_eqn);
    TBOX_ASSERT(static_cast<int>(Q_T.size()) == d_num_eqn);
#endif
    
    // Get the thermodynamic properties of the species.
    std::vector<const double*> thermo_properties_ptr;
    thermo_properties_ptr.reserve(static_cast<int> (d_thermo_properties.size()));
    for (int ti = 0; ti < static_cast<int> (d_thermo_properties.size()); ti++)
    {
        thermo_properties_ptr.push_back(&d_thermo_properties[ti]);
    }
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": RiemannSolverSingleSpeciesHLLC_HLL::"
            << "computeIntercellFluxInYDirectionFromConservativeVariables()\n"
            << "There is no y direction for 1D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        double F_y_intercell_HLLC[d_num_eqn];
        double F_y_intercell_HLL[d_num_eqn];
        
        const double u_B = (Q_B[1].get())/(Q_B[0].get());
        const double u_T = (Q_T[1].get())/(Q_T[0].get());
        
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
            &(Q_B[3].get()),
            thermo_properties_ptr);
        
        const double p_T = d_equation_of_state->getPressure(
            &(Q_T[0].get()),
            m_T,
            &(Q_T[3].get()),
            thermo_properties_ptr);
        
        const double c_B = d_equation_of_state->getSoundSpeed(
            &(Q_B[0].get()),
            &p_B,
            thermo_properties_ptr);
        
        const double c_T = d_equation_of_state->getSoundSpeed(
            &(Q_T[0].get()),
            &p_T,
            thermo_properties_ptr);
        
        const double v_average = 0.5*(v_B + v_T);
        const double c_average = 0.5*(c_B + c_T);
        
        const double s_B = fmin(v_average - c_average, v_B - c_B);
        const double s_T = fmax(v_average + c_average, v_T + c_T);
        
        const double s_minus = fmin(0.0, s_B);
        const double s_plus  = fmax(0.0, s_T);
        
        const double s_star =
            (p_T - p_B + (Q_B[2].get())*(s_B - v_B) - (Q_T[2].get())*(s_T - v_T))/
                ((Q_B[0].get())*(s_B - v_B) - (Q_T[0].get())*(s_T - v_T));
        
        if (s_star > 0)
        {
            /*
             * Compute the HLLC flux.
             */
            
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
                F_y_T[0] = (Q_T[2].get());
                F_y_T[1] = v_T*(Q_T[1].get());
                F_y_T[2] = v_T*(Q_T[2].get()) + p_T;
                F_y_T[3] = v_T*((Q_T[3].get()) + p_T);
                
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
                F_y_B[0] = (Q_B[2].get());
                F_y_B[1] = v_B*(Q_B[1].get());
                F_y_B[2] = v_B*(Q_B[2].get()) + p_B;
                F_y_B[3] = v_B*((Q_B[3].get()) + p_B);
                
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
        
        F_y_intercell[0].get() = beta1*F_y_intercell_HLLC[0] + beta2*F_y_intercell_HLL[0];
        F_y_intercell[1].get() = beta1*F_y_intercell_HLLC[1] + beta2*F_y_intercell_HLL[1];
        F_y_intercell[2].get() = F_y_intercell_HLLC[2];
        F_y_intercell[3].get() = F_y_intercell_HLLC[3];
    }
    else if (d_dim == tbox::Dimension(3))
    {
        double F_y_intercell_HLLC[d_num_eqn];
        double F_y_intercell_HLL[d_num_eqn];
        
        const double u_B = (Q_B[1].get())/(Q_B[0].get());
        const double u_T = (Q_T[1].get())/(Q_T[0].get());
        
        const double v_B = (Q_B[2].get())/(Q_B[0].get());
        const double v_T = (Q_T[2].get())/(Q_T[0].get());
        
        const double w_B = (Q_B[3].get())/(Q_B[0].get());
        const double w_T = (Q_T[3].get())/(Q_T[0].get());
        
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
            &(Q_B[4].get()),
            thermo_properties_ptr);
        
        const double p_T = d_equation_of_state->getPressure(
            &(Q_T[0].get()),
            m_T,
            &(Q_T[0].get()),
            thermo_properties_ptr);
        
        const double c_B = d_equation_of_state->getSoundSpeed(
            &(Q_B[0].get()),
            &p_B,
            thermo_properties_ptr);
        
        const double c_T = d_equation_of_state->getSoundSpeed(
            &(Q_T[0].get()),
            &p_T,
            thermo_properties_ptr);
        
        const double v_average = 0.5*(v_B + v_T);
        const double c_average = 0.5*(c_B + c_T);
        
        const double s_B = fmin(v_average - c_average, v_B - c_B);
        const double s_T = fmax(v_average + c_average, v_T + c_T);
        
        const double s_minus = fmin(0.0, s_B);
        const double s_plus  = fmax(0.0, s_T);
        
        const double s_star =
            (p_T - p_B + (Q_B[2].get())*(s_B - v_B) - (Q_T[2].get())*(s_T - v_T))/
                ((Q_B[0].get())*(s_B - v_B) - (Q_T[0].get())*(s_T - v_T));
        
        if (s_star > 0)
        {
            /*
             * Compute the HLLC flux.
             */
            
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
                F_y_T[0] = (Q_T[2].get());
                F_y_T[1] = v_T*(Q_T[1].get());
                F_y_T[2] = v_T*(Q_T[2].get()) + p_T;
                F_y_T[3] = v_T*(Q_T[3].get());
                F_y_T[4] = v_T*((Q_T[4].get()) + p_T);
                
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
                F_y_B[0] = (Q_B[2].get());
                F_y_B[1] = v_B*(Q_B[1].get());
                F_y_B[2] = v_B*(Q_B[2].get()) + p_B;
                F_y_B[3] = v_B*(Q_B[3].get());
                F_y_B[4] = v_B*((Q_B[4].get()) + p_B);
                
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
        
        F_y_intercell[0].get() = beta1*F_y_intercell_HLLC[0] + beta2*F_y_intercell_HLL[0];
        F_y_intercell[1].get() = beta1*F_y_intercell_HLLC[1] + beta2*F_y_intercell_HLL[1];
        F_y_intercell[2].get() = F_y_intercell_HLLC[2];
        F_y_intercell[3].get() = beta1*F_y_intercell_HLLC[3] + beta2*F_y_intercell_HLL[3];
        F_y_intercell[4].get() = F_y_intercell_HLLC[4];
    }
}


/*
 * Compute the flux in the z-direction at the intercell face
 * from conservative variables.
 */
void
RiemannSolverSingleSpeciesHLLC_HLL::computeIntercellFluxInZDirectionFromConservativeVariables(
    std::vector<boost::reference_wrapper<double> >& F_z_intercell,
    const std::vector<boost::reference_wrapper<double> >& Q_B,
    const std::vector<boost::reference_wrapper<double> >& Q_F)
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(F_z_intercell.size()) == d_num_eqn);
    TBOX_ASSERT(static_cast<int>(Q_B.size()) == d_num_eqn);
    TBOX_ASSERT(static_cast<int>(Q_F.size()) == d_num_eqn);
#endif
    
    // Get the thermodynamic properties of the species.
    std::vector<const double*> thermo_properties_ptr;
    thermo_properties_ptr.reserve(static_cast<int> (d_thermo_properties.size()));
    for (int ti = 0; ti < static_cast<int> (d_thermo_properties.size()); ti++)
    {
        thermo_properties_ptr.push_back(&d_thermo_properties[ti]);
    }
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": RiemannSolverSingleSpeciesHLLC_HLL::"
            << "computeIntercellFluxInZDirectionFromConservativeVariables()\n"
            << "There is no z direction for 1D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": RiemannSolverSingleSpeciesHLLC_HLL::"
            << "computeIntercellFluxInZDirectionFromConservativeVariables()\n"
            << "There is no z direction for 2D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        double F_z_intercell_HLLC[d_num_eqn];
        double F_z_intercell_HLL[d_num_eqn];
        
        const double u_B = (Q_B[1].get())/(Q_B[0].get());
        const double u_F = (Q_F[1].get())/(Q_F[0].get());
        
        const double v_B = (Q_B[2].get())/(Q_B[0].get());
        const double v_F = (Q_F[2].get())/(Q_F[0].get());
        
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
            &(Q_B[4].get()),
            thermo_properties_ptr);
        
        const double p_F = d_equation_of_state->getPressure(
            &(Q_F[0].get()),
            m_F,
            &(Q_F[4].get()),
            thermo_properties_ptr);
        
        const double c_B = d_equation_of_state->getSoundSpeed(
            &(Q_B[0].get()),
            &p_B,
            thermo_properties_ptr);
        
        const double c_F = d_equation_of_state->getSoundSpeed(
            &(Q_F[0].get()),
            &p_F,
            thermo_properties_ptr);
        
        const double w_average = 0.5*(w_B + w_F);
        const double c_average = 0.5*(c_B + c_F);
        
        const double s_B = fmin(w_average - c_average, w_B - c_B);
        const double s_F = fmax(w_average + c_average, w_F + c_F);
        
        const double s_minus = fmin(0.0, s_B);
        const double s_plus  = fmax(0.0, s_F);
        
        const double s_star =
            (p_F - p_B + (Q_B[3].get())*(s_B - w_B) - (Q_F[3].get())*(s_F - w_F))/
                ((Q_B[0].get())*(s_B - w_B) - (Q_F[0].get())*(s_F - w_F));
        
        if (s_star > 0)
        {
            /*
             * Compute the HLLC flux.
             */
            
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
                F_z_F[0] = (Q_F[3].get());
                F_z_F[1] = w_F*(Q_F[1].get());
                F_z_F[2] = w_F*(Q_F[2].get());
                F_z_F[3] = w_F*(Q_F[3].get()) + p_F;
                F_z_F[4] = w_F*((Q_F[4].get()) + p_F);
                
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
                F_z_B[0] = (Q_B[3].get());
                F_z_B[1] = w_B*(Q_B[1].get());
                F_z_B[2] = w_B*(Q_B[2].get());
                F_z_B[3] = w_B*(Q_B[3].get()) + p_B;
                F_z_B[4] = w_B*((Q_B[4].get()) + p_B);
                
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
        
        F_z_intercell[0].get() = beta1*F_z_intercell_HLLC[0] + beta2*F_z_intercell_HLL[0];
        F_z_intercell[1].get() = beta1*F_z_intercell_HLLC[1] + beta2*F_z_intercell_HLL[1];
        F_z_intercell[2].get() = beta1*F_z_intercell_HLLC[2] + beta2*F_z_intercell_HLL[2];
        F_z_intercell[3].get() = F_z_intercell_HLLC[3];
        F_z_intercell[4].get() = F_z_intercell_HLLC[4];
    }
}


/*
 * Compute the flux in the x-direction at the intercell face
 * from primitive variables.
 */
void
RiemannSolverSingleSpeciesHLLC_HLL::computeIntercellFluxInXDirectionFromPrimitiveVariables(
    std::vector<boost::reference_wrapper<double> >& F_x_intercell,
    const std::vector<boost::reference_wrapper<double> >& V_L,
    const std::vector<boost::reference_wrapper<double> >& V_R)
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(F_x_intercell.size()) == d_num_eqn);
    TBOX_ASSERT(static_cast<int>(V_L.size()) == d_num_eqn);
    TBOX_ASSERT(static_cast<int>(V_R.size()) == d_num_eqn);
#endif
    
    // Get the thermodynamic properties of the species.
    std::vector<const double*> thermo_properties_ptr;
    thermo_properties_ptr.reserve(static_cast<int> (d_thermo_properties.size()));
    for (int ti = 0; ti < static_cast<int> (d_thermo_properties.size()); ti++)
    {
        thermo_properties_ptr.push_back(&d_thermo_properties[ti]);
    }
    
    if (d_dim == tbox::Dimension(1))
    {
        const double c_L = d_equation_of_state->getSoundSpeed(
            &(V_L[0].get()),
            &(V_L[2].get()),
            thermo_properties_ptr);
        
        const double c_R = d_equation_of_state->getSoundSpeed(
            &(V_R[0].get()),
            &(V_R[2].get()),
            thermo_properties_ptr);
        
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
                &(V_L[2].get()),
                thermo_properties_ptr);
            
            double Q_star_L[d_num_eqn];
            Q_star_L[0] = Chi_star_L*(V_L[0].get());
            Q_star_L[1] = Chi_star_L*(V_L[0].get())*s_star;
            Q_star_L[2] = Chi_star_L*(Q_L[2] + (s_star - (V_L[1].get()))*((V_L[0].get())*s_star +
                (V_L[2].get())/(s_L - (V_L[1].get()))));
            
            double F_x_L[d_num_eqn];
            F_x_L[0] = Q_L[1];
            F_x_L[1] = Q_L[1]*(V_L[1].get()) + (V_L[2].get());
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
                &(V_R[2].get()),
                thermo_properties_ptr);
            
            double Q_star_R[d_num_eqn];
            Q_star_R[0] = Chi_star_R*(V_R[0].get());
            Q_star_R[1] = Chi_star_R*(V_R[0].get())*s_star;
            Q_star_R[2] = Chi_star_R*(Q_R[2] + (s_star - (V_R[1].get()))*((V_R[0].get())*s_star +
                (V_R[2].get())/(s_R - (V_R[1].get()))));
            
            double F_x_R[d_num_eqn];
            F_x_R[0] = Q_R[1];
            F_x_R[1] = Q_R[1]*(V_R[1].get()) + (V_R[2].get());
            F_x_R[2] = (V_R[1].get())*(Q_R[2] + (V_R[2].get()));
            
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
        
        const double c_L = d_equation_of_state->getSoundSpeed(
            &(V_L[0].get()),
            &(V_L[3].get()),
            thermo_properties_ptr);
        
        const double c_R = d_equation_of_state->getSoundSpeed(
            &(V_R[0].get()),
            &(V_R[3].get()),
            thermo_properties_ptr);
        
        const double u_average = 0.5*((V_L[1].get()) + (V_R[1].get()));
        const double c_average = 0.5*(c_L + c_R);
        
        const double s_L = fmin(u_average - c_average, (V_L[1].get()) - c_L);
        const double s_R = fmax(u_average + c_average, (V_R[1].get()) + c_R);
        
        const double s_minus = fmin(0.0, s_L);
        const double s_plus  = fmax(0.0, s_R);
        
        const double s_star =
            ((V_R[3].get()) - (V_L[3].get()) + (V_L[0].get())*(V_L[1].get())*
                (s_L - (V_L[1].get())) -(V_R[0].get())*(V_R[1].get())*(s_R - (V_R[1].get())))/
                    ((V_L[0].get())*(s_L - (V_L[1].get())) - (V_R[0].get())*(s_R - (V_R[1].get())));
        
        if (s_star > 0)
        {
            /*
             * Compute the HLLC flux.
             */
            
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
                &(V_L[3].get()),
                thermo_properties_ptr);
            
            double Q_star_L[d_num_eqn];
            Q_star_L[0] = Chi_star_L*(V_L[0].get());
            Q_star_L[1] = Chi_star_L*(V_L[0].get())*s_star;
            Q_star_L[2] = Chi_star_L*Q_L[2];
            Q_star_L[3] = Chi_star_L*(Q_L[3] + (s_star - (V_L[1].get()))*((V_L[0].get())*s_star +
                (V_L[3].get())/(s_L - (V_L[1].get()))));
            
            double F_x_L[d_num_eqn];
            F_x_L[0] = Q_L[1];
            F_x_L[1] = Q_L[1]*(V_L[1].get()) + (V_L[3].get());
            F_x_L[2] = Q_L[1]*(V_L[2].get());
            F_x_L[3] = (V_L[1].get())*(Q_L[3] + (V_L[3].get()));
            
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
                std::vector<const double*> vel_R;
                vel_R.reserve(2);
                vel_R.push_back(&(V_R[1].get()));
                vel_R.push_back(&(V_R[2].get()));
                
                double E_R = d_equation_of_state->getTotalEnergy(
                    &(V_R[0].get()),
                    vel_R,
                    &(V_R[3].get()),
                    thermo_properties_ptr);
                
                if (s_R <= 0)
                {
                    F_x_intercell_HLL[0] = (V_R[0].get())*(V_R[1].get());
                    F_x_intercell_HLL[1] = F_x_intercell_HLL[0]*(V_R[1].get()) + (V_R[3].get());
                    F_x_intercell_HLL[2] = F_x_intercell_HLL[0]*(V_R[2].get());
                    F_x_intercell_HLL[3] = (V_R[1].get())*(E_R + (V_R[3].get()));
                }
                else
                {
                    double Q_R[d_num_eqn];
                    Q_R[0] = (V_R[0].get());
                    Q_R[1] = (V_R[0].get())*(V_R[1].get());
                    Q_R[2] = (V_R[0].get())*(V_R[2].get());
                    Q_R[3] = E_R;
                    
                    double F_x_R[d_num_eqn];
                    F_x_R[0] = Q_R[1];
                    F_x_R[1] = Q_R[1]*(V_R[1].get()) + (V_R[3].get());
                    F_x_R[2] = Q_R[1]*(V_R[2].get());
                    F_x_R[3] = (V_R[1].get())*(Q_R[3] + (V_R[3].get()));
                    
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
                &(V_R[3].get()),
                thermo_properties_ptr);
            
            double Q_star_R[d_num_eqn];
            Q_star_R[0] = Chi_star_R*(V_R[0].get());
            Q_star_R[1] = Chi_star_R*(V_R[0].get())*s_star;
            Q_star_R[2] = Chi_star_R*Q_R[2];
            Q_star_R[3] = Chi_star_R*(Q_R[3] + (s_star - (V_R[1].get()))*((V_R[0].get())*s_star +
                (V_R[3].get())/(s_R - (V_R[1].get()))));
            
            double F_x_R[d_num_eqn];
            F_x_R[0] = Q_R[1];
            F_x_R[1] = Q_R[1]*(V_R[1].get()) + (V_R[3].get());
            F_x_R[2] = Q_R[1]*(V_R[2].get());
            F_x_R[3] = (V_R[1].get())*(Q_R[3] + (V_R[3].get()));
            
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
                std::vector<const double*> vel_L;
                vel_L.reserve(2);
                vel_L.push_back(&(V_L[1].get()));
                vel_L.push_back(&(V_L[2].get()));
                
                double E_L = d_equation_of_state->getTotalEnergy(
                    &(V_L[0].get()),
                    vel_L,
                    &(V_L[3].get()),
                    thermo_properties_ptr);
                
                if (s_L >= 0)
                {
                    F_x_intercell_HLL[0] = (V_L[0].get())*(V_L[1].get());
                    F_x_intercell_HLL[1] = F_x_intercell_HLL[0]*(V_L[1].get()) + (V_L[3].get());
                    F_x_intercell_HLL[2] = F_x_intercell_HLL[0]*(V_L[2].get());
                    F_x_intercell_HLL[3] = (V_L[1].get())*(E_L + (V_L[3].get()));
                }
                else
                {
                    double Q_L[d_num_eqn];
                    Q_L[0] = (V_L[0].get());
                    Q_L[1] = (V_L[0].get())*(V_L[1].get());
                    Q_L[2] = (V_L[0].get())*(V_L[2].get());
                    Q_L[3] = E_L;
                    
                    double F_x_L[d_num_eqn];
                    F_x_L[0] = Q_L[1];
                    F_x_L[1] = Q_L[1]*(V_L[1].get()) + (V_L[3].get());
                    F_x_L[2] = Q_L[1]*(V_L[2].get());
                    F_x_L[3] = (V_L[1].get())*(Q_L[3] + (V_L[3].get()));
                    
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
        
        double vel_mag = sqrt(pow((V_R[1].get()) - (V_L[1].get()), 2) + pow((V_R[2].get()) - (V_L[2].get()), 2));
        
        if (vel_mag < EPSILON)
        {
            alpha1 = 1.0;
            alpha2 = 0.0;
        }
        else
        {
            alpha1 = fabs((V_R[1].get()) - (V_L[1].get()))/vel_mag;
            alpha2 = sqrt(1.0 - alpha1*alpha1);
        }
        
        double beta1 = 0.5 + 0.5*alpha1/(alpha1 + alpha2);
        double beta2 = 1.0 - beta1;
        
        /*
         * Compute the HLLC-HLL flux.
         */
        
        F_x_intercell[0].get()= beta1*F_x_intercell_HLLC[0] + beta2*F_x_intercell_HLL[0];
        F_x_intercell[1].get() = F_x_intercell_HLLC[1];
        F_x_intercell[2].get() = beta1*F_x_intercell_HLLC[2] + beta2*F_x_intercell_HLL[2];
        F_x_intercell[3].get() = F_x_intercell_HLLC[3];
    }
    else if (d_dim == tbox::Dimension(3))
    {
        double F_x_intercell_HLLC[d_num_eqn];
        double F_x_intercell_HLL[d_num_eqn];
        
        const double c_L = d_equation_of_state->getSoundSpeed(
            &(V_L[0].get()),
            &(V_L[4].get()),
            thermo_properties_ptr);
        
        const double c_R = d_equation_of_state->getSoundSpeed(
            &(V_R[0].get()),
            &(V_R[4].get()),
            thermo_properties_ptr);
        
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
            /*
             * Compute the HLLC flux.
             */
            
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
                &(V_L[4].get()),
                thermo_properties_ptr);
            
            double Q_star_L[d_num_eqn];
            Q_star_L[0] = Chi_star_L*(V_L[0].get());
            Q_star_L[1] = Chi_star_L*(V_L[0].get())*s_star;
            Q_star_L[2] = Chi_star_L*Q_L[2];
            Q_star_L[3] = Chi_star_L*Q_L[3];
            Q_star_L[4] = Chi_star_L*(Q_L[4] + (s_star - (V_L[1].get()))*((V_L[0].get())*s_star +
                (V_L[4].get())/(s_L - (V_L[1].get()))));
            
            double F_x_L[d_num_eqn];
            F_x_L[0] = Q_L[1];
            F_x_L[1] = Q_L[1]*(V_L[1].get()) + (V_L[4].get());
            F_x_L[2] = Q_L[1]*(V_L[2].get());
            F_x_L[3] = Q_L[1]*(V_L[3].get());
            F_x_L[4] = (V_L[1].get())*(Q_L[4] + (V_L[4].get()));
            
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
                std::vector<const double*> vel_R;
                vel_R.reserve(3);
                vel_R.push_back(&(V_R[1].get()));
                vel_R.push_back(&(V_R[2].get()));
                vel_R.push_back(&(V_R[3].get()));
                
                double E_R = d_equation_of_state->getTotalEnergy(
                    &(V_R[0].get()),
                    vel_R,
                    &(V_R[4].get()),
                    thermo_properties_ptr);
                
                if (s_R <= 0)
                {
                    F_x_intercell_HLL[0] = (V_R[0].get())*(V_R[1].get());
                    F_x_intercell_HLL[1] = F_x_intercell_HLL[0]*(V_R[1].get()) + (V_R[4].get());
                    F_x_intercell_HLL[2] = F_x_intercell_HLL[0]*(V_R[2].get());
                    F_x_intercell_HLL[3] = F_x_intercell_HLL[0]*(V_R[3].get());
                    F_x_intercell_HLL[4] = (V_R[1].get())*(E_R + (V_R[4].get()));
                }
                else
                {
                    double Q_R[d_num_eqn];
                    Q_R[0] = (V_R[0].get());
                    Q_R[1] = (V_R[0].get())*(V_R[1].get());
                    Q_R[2] = (V_R[0].get())*(V_R[2].get());
                    Q_R[3] = (V_R[0].get())*(V_R[3].get());
                    Q_R[4] = E_R;
                    
                    double F_x_R[d_num_eqn];
                    F_x_R[0] = Q_R[1];
                    F_x_R[1] = Q_R[1]*(V_R[1].get()) + (V_R[4].get());
                    F_x_R[2] = Q_R[1]*(V_R[2].get());
                    F_x_R[3] = Q_R[1]*(V_R[3].get());
                    F_x_R[4] = (V_R[1].get())*(Q_R[4] + (V_R[4].get()));
                    
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
                &(V_R[4].get()),
                thermo_properties_ptr);
            
            double Q_star_R[d_num_eqn];
            Q_star_R[0] = Chi_star_R*(V_R[0].get());
            Q_star_R[1] = Chi_star_R*(V_R[0].get())*s_star;
            Q_star_R[2] = Chi_star_R*Q_R[2];
            Q_star_R[3] = Chi_star_R*Q_R[3];
            Q_star_R[4] = Chi_star_R*(Q_R[4] + (s_star - (V_R[1].get()))*((V_R[0].get())*s_star +
                (V_R[4].get())/(s_R - (V_R[1].get()))));
            
            double F_x_R[d_num_eqn];
            F_x_R[0] = Q_R[1];
            F_x_R[1] = Q_R[1]*(V_R[1].get()) + (V_R[4].get());
            F_x_R[2] = Q_R[1]*(V_R[2].get());
            F_x_R[3] = Q_R[1]*(V_R[3].get());
            F_x_R[4] = (V_R[1].get())*(Q_R[4] + (V_R[4].get()));
            
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
                std::vector<const double*> vel_L;
                vel_L.reserve(3);
                vel_L.push_back(&(V_L[1].get()));
                vel_L.push_back(&(V_L[2].get()));
                vel_L.push_back(&(V_L[3].get()));
                
                double E_L = d_equation_of_state->getTotalEnergy(
                    &(V_L[0].get()),
                    vel_L,
                    &(V_L[4].get()),
                    thermo_properties_ptr);
                
                if (s_L >= 0)
                {
                    F_x_intercell_HLL[0] = (V_L[0].get())*(V_L[1].get());
                    F_x_intercell_HLL[1] = F_x_intercell_HLL[0]*(V_L[1].get()) + (V_L[3].get());
                    F_x_intercell_HLL[2] = F_x_intercell_HLL[0]*(V_L[2].get());
                    F_x_intercell_HLL[3] = F_x_intercell_HLL[0]*(V_L[3].get());
                    F_x_intercell_HLL[4] = (V_L[1].get())*(E_L + (V_L[4].get()));
                }
                else
                {
                    double Q_L[d_num_eqn];
                    Q_L[0] = (V_L[0].get());
                    Q_L[1] = (V_L[0].get())*(V_L[1].get());
                    Q_L[2] = (V_L[0].get())*(V_L[2].get());
                    Q_L[3] = (V_L[0].get())*(V_L[3].get());
                    Q_L[4] = E_L;
                    
                    double F_x_L[d_num_eqn];
                    F_x_L[0] = Q_L[1];
                    F_x_L[1] = Q_L[1]*(V_L[1].get()) + (V_L[4].get());
                    F_x_L[2] = Q_L[1]*(V_L[2].get());
                    F_x_L[3] = Q_L[1]*(V_L[3].get());
                    F_x_L[4] = (V_L[1].get())*(Q_L[4] + (V_L[4].get()));
                    
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
        
        double vel_mag = sqrt(pow((V_R[1].get()) - (V_L[1].get()), 2) + pow((V_R[2].get()) -
            (V_L[2].get()), 2) + pow((V_R[3].get()) - (V_L[3].get()), 2));
        
        if (vel_mag < EPSILON)
        {
            alpha1 = 1.0;
            alpha2 = 0.0;
        }
        else
        {
            alpha1 = fabs((V_R[1].get()) - (V_L[1].get()))/vel_mag;
            alpha2 = sqrt(1.0 - alpha1*alpha1);
        }
        
        double beta1 = 0.5 + 0.5*alpha1/(alpha1 + alpha2);
        double beta2 = 1.0 - beta1;
        
        /*
         * Compute the HLLC-HLL flux.
         */
        
        F_x_intercell[0].get()= beta1*F_x_intercell_HLLC[0] + beta2*F_x_intercell_HLL[0];
        F_x_intercell[1].get() = F_x_intercell_HLLC[1];
        F_x_intercell[2].get() = beta1*F_x_intercell_HLLC[2] + beta2*F_x_intercell_HLL[2];
        F_x_intercell[3].get() = beta1*F_x_intercell_HLLC[3] + beta2*F_x_intercell_HLL[3];
        F_x_intercell[4].get() = F_x_intercell_HLLC[4];
    }
}


/*
 * Compute the flux in the y-direction at the intercell face
 * from primitive variables.
 */
void
RiemannSolverSingleSpeciesHLLC_HLL::computeIntercellFluxInYDirectionFromPrimitiveVariables(
    std::vector<boost::reference_wrapper<double> >& F_y_intercell,
    const std::vector<boost::reference_wrapper<double> >& V_B,
    const std::vector<boost::reference_wrapper<double> >& V_T)
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(F_y_intercell.size()) == d_num_eqn);
    TBOX_ASSERT(static_cast<int>(V_B.size()) == d_num_eqn);
    TBOX_ASSERT(static_cast<int>(V_T.size()) == d_num_eqn);
#endif
    
    // Get the thermodynamic properties of the species.
    std::vector<const double*> thermo_properties_ptr;
    thermo_properties_ptr.reserve(static_cast<int> (d_thermo_properties.size()));
    for (int ti = 0; ti < static_cast<int> (d_thermo_properties.size()); ti++)
    {
        thermo_properties_ptr.push_back(&d_thermo_properties[ti]);
    }
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": RiemannSolverSingleSpeciesHLLC_HLL::"
            << "computeIntercellFluxInYDirectionFromPrimitiveVariables()\n"
            << "There is no y direction for 1D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        double F_y_intercell_HLLC[d_num_eqn];
        double F_y_intercell_HLL[d_num_eqn];
        
        const double c_B = d_equation_of_state->getSoundSpeed(
            &(V_B[0].get()),
            &(V_B[3].get()),
            thermo_properties_ptr);
        
        const double c_T = d_equation_of_state->getSoundSpeed(
            &(V_T[0].get()),
            &(V_T[3].get()),
            thermo_properties_ptr);
        
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
            /*
             * Compute the HLLC flux.
             */
            
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
                &(V_B[3].get()),
                thermo_properties_ptr);
            
            double Q_star_B[d_num_eqn];
            Q_star_B[0] = Chi_star_B*(V_B[0].get());
            Q_star_B[1] = Chi_star_B*Q_B[1];
            Q_star_B[2] = Chi_star_B*(V_B[0].get())*s_star;
            Q_star_B[3] = Chi_star_B*(Q_B[3] + (s_star - (V_B[2].get()))*((V_B[0].get())*s_star + (V_B[3].get())/(s_B - (V_B[2].get()))));
            
            double F_y_B[d_num_eqn];
            F_y_B[0] = Q_B[2];
            F_y_B[1] = Q_B[2]*(V_B[1].get());
            F_y_B[2] = Q_B[2]*(V_B[2].get()) + (V_B[3].get());
            F_y_B[3] = (V_B[2].get())*(Q_B[3] + (V_B[3].get()));
            
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
                std::vector<const double*> vel_T;
                vel_T.reserve(2);
                vel_T.push_back(&(V_T[1].get()));
                vel_T.push_back(&(V_T[2].get()));
                
                double E_T = d_equation_of_state->getTotalEnergy(
                    &(V_T[0].get()),
                    vel_T,
                    &(V_T[3].get()),
                    thermo_properties_ptr);
                
                if (s_T <= 0)
                {
                    F_y_intercell_HLL[0] = (V_T[0].get())*(V_T[2].get());
                    F_y_intercell_HLL[1] = F_y_intercell_HLL[0]*(V_T[1].get());
                    F_y_intercell_HLL[2] = F_y_intercell_HLL[0]*(V_T[2].get()) + (V_T[3].get());
                    F_y_intercell_HLL[3] = (V_T[2].get())*(E_T + (V_T[3].get()));
                }
                else
                {
                    double Q_T[d_num_eqn];
                    Q_T[0] = (V_T[0].get());
                    Q_T[1] = (V_T[0].get())*(V_T[1].get());
                    Q_T[2] = (V_T[0].get())*(V_T[2].get());
                    Q_T[3] = E_T;
                    
                    double F_y_T[d_num_eqn];
                    F_y_T[0] = Q_T[2];
                    F_y_T[1] = Q_T[2]*(V_T[1].get());
                    F_y_T[2] = Q_T[2]*(V_T[2].get()) + (V_T[3].get());
                    F_y_T[3] = (V_T[2].get())*(Q_T[3] + (V_T[3].get()));
                    
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
                &(V_T[3].get()),
                thermo_properties_ptr);
            
            double Q_star_T[d_num_eqn];
            Q_star_T[0] = Chi_star_T*(V_T[0].get());
            Q_star_T[1] = Chi_star_T*Q_T[1];
            Q_star_T[2] = Chi_star_T*(V_T[0].get())*s_star;
            Q_star_T[3] = Chi_star_T*(Q_T[3] + (s_star - (V_T[2].get()))*((V_T[0].get())*s_star +
                (V_T[3].get())/(s_T - (V_T[2].get()))));
            
            double F_y_T[d_num_eqn];
            F_y_T[0] = Q_T[2];
            F_y_T[1] = Q_T[2]*(V_T[1].get());
            F_y_T[2] = Q_T[2]*(V_T[2].get()) + (V_T[3].get());
            F_y_T[3] = (V_T[2].get())*(Q_T[3] + (V_T[3].get()));
            
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
                std::vector<const double*> vel_B;
                vel_B.reserve(2);
                vel_B.push_back(&(V_B[1].get()));
                vel_B.push_back(&(V_B[2].get()));
                
                double E_B = d_equation_of_state->getTotalEnergy(
                    &(V_B[0].get()),
                    vel_B,
                    &(V_B[3].get()),
                    thermo_properties_ptr);
                
                if (s_B >= 0)
                {
                    F_y_intercell_HLL[0] = (V_B[0].get())*(V_B[2].get());
                    F_y_intercell_HLL[1] = F_y_intercell_HLL[0]*(V_B[1].get());
                    F_y_intercell_HLL[2] = F_y_intercell_HLL[0]*(V_B[2].get()) + (V_B[3].get());
                    F_y_intercell_HLL[3] = (V_B[2].get())*(E_B + (V_B[3].get()));
                }
                else
                {
                    double Q_B[d_num_eqn];
                    Q_B[0] = (V_B[0].get());
                    Q_B[1] = (V_B[0].get())*(V_B[1].get());
                    Q_B[2] = (V_B[0].get())*(V_B[2].get());
                    Q_B[3] = E_B;
                    
                    double F_y_B[d_num_eqn];
                    F_y_B[0] = Q_B[2];
                    F_y_B[1] = Q_B[2]*(V_B[1].get());
                    F_y_B[2] = Q_B[2]*(V_B[2].get()) + (V_B[3].get());
                    F_y_B[3] = (V_B[2].get())*(Q_B[3] + (V_B[3].get()));
                    
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
        
        double vel_mag = sqrt(pow((V_T[1].get()) - (V_B[1].get()), 2) + pow((V_T[2].get()) - (V_B[2].get()), 2));
        
        if (vel_mag < EPSILON)
        {
            alpha1 = 1.0;
            alpha2 = 0.0;
        }
        else
        {
            alpha1 = fabs((V_T[2].get()) - (V_B[2].get()))/vel_mag;
            alpha2 = sqrt(1.0 - alpha1*alpha1);
        }
        
        double beta1 = 0.5 + 0.5*alpha1/(alpha1 + alpha2);
        double beta2 = 1.0 - beta1;
        
        /*
         * Compute the HLLC-HLL flux.
         */
        
        F_y_intercell[0].get() = beta1*F_y_intercell_HLLC[0] + beta2*F_y_intercell_HLL[0];
        F_y_intercell[1].get() = beta1*F_y_intercell_HLLC[1] + beta2*F_y_intercell_HLL[1];
        F_y_intercell[2].get() = F_y_intercell_HLLC[2];
        F_y_intercell[3].get() = F_y_intercell_HLLC[3];
    }
    else if (d_dim == tbox::Dimension(3))
    {
        double F_y_intercell_HLLC[d_num_eqn];
        double F_y_intercell_HLL[d_num_eqn];
        
        const double c_B = d_equation_of_state->getSoundSpeed(
            &(V_B[0].get()),
            &(V_B[4].get()),
            thermo_properties_ptr);
        
        const double c_T = d_equation_of_state->getSoundSpeed(
            &(V_T[0].get()),
            &(V_T[4].get()),
            thermo_properties_ptr);
        
        const double v_average = 0.5*((V_B[2].get()) + (V_T[2].get()));
        const double c_average = 0.5*(c_B + c_T);
        
        const double s_B = fmin(v_average - c_average, (V_B[2].get()) - c_B);
        const double s_T = fmax(v_average + c_average, (V_T[2].get()) + c_T);
        
        const double s_minus = fmin(0.0, s_B);
        const double s_plus  = fmax(0.0, s_T);
        
        const double s_star =
            ((V_T[4].get()) - (V_B[4].get()) + (V_B[0].get())*(V_B[2].get())*
                (s_B - (V_B[2].get())) -(V_T[0].get())*(V_T[2].get())*(s_T - (V_T[2].get())))/
                    ((V_B[0].get())*(s_B - (V_B[2].get())) - (V_T[0].get())*(s_T - (V_T[2].get())));
        
        if (s_star > 0)
        {
            /*
             * Compute the HLLC flux.
             */
            
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
                &(V_B[4].get()),
                thermo_properties_ptr);
            
            double Q_star_B[d_num_eqn];
            Q_star_B[0] = Chi_star_B*(V_B[0].get());
            Q_star_B[1] = Chi_star_B*Q_B[1];
            Q_star_B[2] = Chi_star_B*(V_B[0].get())*s_star;
            Q_star_B[3] = Chi_star_B*Q_B[3];
            Q_star_B[4] = Chi_star_B*(Q_B[4] + (s_star - (V_B[2].get()))*((V_B[0].get())*s_star +
                (V_B[4].get())/(s_B - (V_B[2].get()))));
            
            double F_y_B[d_num_eqn];
            F_y_B[0] = Q_B[2];
            F_y_B[1] = Q_B[2]*(V_B[1].get());
            F_y_B[2] = Q_B[2]*(V_B[2].get()) + (V_B[4].get());
            F_y_B[3] = Q_B[2]*(V_B[3].get());
            F_y_B[4] = (V_B[2].get())*(Q_B[4] + (V_B[4].get()));
            
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
                std::vector<const double*> vel_T;
                vel_T.reserve(3);
                vel_T.push_back(&(V_T[1].get()));
                vel_T.push_back(&(V_T[2].get()));
                vel_T.push_back(&(V_T[3].get()));
                
                double E_T = d_equation_of_state->getTotalEnergy(
                    &(V_T[0].get()),
                    vel_T,
                    &(V_T[4].get()),
                    thermo_properties_ptr);
                
                if (s_T <= 0)
                {
                    F_y_intercell_HLL[0] = (V_T[0].get())*(V_T[2].get());
                    F_y_intercell_HLL[1] = F_y_intercell_HLL[0]*(V_T[1].get());
                    F_y_intercell_HLL[2] = F_y_intercell_HLL[0]*(V_T[2].get()) + (V_T[4].get());
                    F_y_intercell_HLL[3] = F_y_intercell_HLL[0]*(V_T[3].get());
                    F_y_intercell_HLL[4] = (V_T[2].get())*(E_T + (V_T[4].get()));
                }
                else
                {
                    double Q_T[d_num_eqn];
                    Q_T[0] = (V_T[0].get());
                    Q_T[1] = (V_T[0].get())*(V_T[1].get());
                    Q_T[2] = (V_T[0].get())*(V_T[2].get());
                    Q_T[3] = (V_T[0].get())*(V_T[3].get());
                    Q_T[4] = E_T;
                    
                    double F_y_T[d_num_eqn];
                    F_y_T[0] = Q_T[2];
                    F_y_T[1] = Q_T[2]*(V_T[1].get());
                    F_y_T[2] = Q_T[2]*(V_T[2].get()) + (V_T[4].get());
                    F_y_T[3] = Q_T[2]*(V_T[3].get());
                    F_y_T[4] = (V_T[2].get())*(Q_T[4] + (V_T[4].get()));
                    
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
                &(V_T[4].get()),
                thermo_properties_ptr);
            
            double Q_star_T[d_num_eqn];
            Q_star_T[0] = Chi_star_T*(V_T[0].get());
            Q_star_T[1] = Chi_star_T*Q_T[1];
            Q_star_T[2] = Chi_star_T*(V_T[0].get())*s_star;
            Q_star_T[3] = Chi_star_T*Q_T[3];
            Q_star_T[4] = Chi_star_T*(Q_T[4] + (s_star - (V_T[2].get()))*((V_T[0].get())*s_star +
                (V_T[4].get())/(s_T - (V_T[2].get()))));
            
            double F_y_T[d_num_eqn];
            F_y_T[0] = Q_T[2];
            F_y_T[1] = Q_T[2]*(V_T[1].get());
            F_y_T[2] = Q_T[2]*(V_T[2].get()) + (V_T[4].get());
            F_y_T[3] = Q_T[2]*(V_T[3].get());
            F_y_T[4] = (V_T[2].get())*(Q_T[4] + (V_T[4].get()));
            
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
                std::vector<const double*> vel_B;
                vel_B.reserve(3);
                vel_B.push_back(&(V_B[1].get()));
                vel_B.push_back(&(V_B[2].get()));
                vel_B.push_back(&(V_B[3].get()));
                
                double E_B = d_equation_of_state->getTotalEnergy(
                    &(V_B[0].get()),
                    vel_B,
                    &(V_B[4].get()),
                    thermo_properties_ptr);
                
                if (s_B >= 0)
                {
                    F_y_intercell_HLL[0] = (V_B[0].get())*(V_B[2].get());
                    F_y_intercell_HLL[1] = F_y_intercell_HLL[0]*(V_B[1].get());
                    F_y_intercell_HLL[2] = F_y_intercell_HLL[0]*(V_B[2].get()) + (V_B[4].get());
                    F_y_intercell_HLL[3] = F_y_intercell_HLL[0]*(V_B[3].get());
                    F_y_intercell_HLL[4] = (V_B[2].get())*(E_B + (V_B[4].get()));
                }
                else
                {
                    double Q_B[d_num_eqn];
                    Q_B[0] = (V_B[0].get());
                    Q_B[1] = (V_B[0].get())*(V_B[1].get());
                    Q_B[2] = (V_B[0].get())*(V_B[2].get());
                    Q_B[3] = (V_B[0].get())*(V_B[3].get());
                    Q_B[4] = E_B;
                    
                    double F_y_B[d_num_eqn];
                    F_y_B[0] = Q_B[2];
                    F_y_B[1] = Q_B[2]*(V_B[1].get());
                    F_y_B[2] = Q_B[2]*(V_B[2].get()) + (V_B[4].get());
                    F_y_B[3] = Q_B[2]*(V_B[3].get());
                    F_y_B[4] = (V_B[2].get())*(Q_B[4] + (V_B[4].get()));
                    
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
        
        double vel_mag = sqrt(pow((V_T[1].get()) - (V_B[1].get()), 2) + pow((V_T[2].get()) -
            (V_B[2].get()), 2) +pow((V_T[3].get()) - (V_B[3].get()), 2));
        
        if (vel_mag < EPSILON)
        {
           alpha1 = 1.0;
           alpha2 = 0.0;
        }
        else
        {
           alpha1 = fabs((V_T[2].get()) - (V_B[2].get()))/vel_mag;
           alpha2 = sqrt(1.0 - alpha1*alpha1);
        }
        
        double beta1 = 0.5 + 0.5*alpha1/(alpha1 + alpha2);
        double beta2 = 1.0 - beta1;
        
        /*
         * Compute the HLLC-HLL flux.
         */
        
        F_y_intercell[0].get() = beta1*F_y_intercell_HLLC[0] + beta2*F_y_intercell_HLL[0];
        F_y_intercell[1].get() = beta1*F_y_intercell_HLLC[1] + beta2*F_y_intercell_HLL[1];
        F_y_intercell[2].get() = F_y_intercell_HLLC[2];
        F_y_intercell[3].get() = beta1*F_y_intercell_HLLC[3] + beta2*F_y_intercell_HLL[3];
        F_y_intercell[4].get() = F_y_intercell_HLLC[4];
    }
}


/*
 * Compute the flux in the z-direction at the intercell face
 * from primitive variables.
 */
void
RiemannSolverSingleSpeciesHLLC_HLL::computeIntercellFluxInZDirectionFromPrimitiveVariables(
    std::vector<boost::reference_wrapper<double> >& F_z_intercell,
    const std::vector<boost::reference_wrapper<double> >& V_B,
    const std::vector<boost::reference_wrapper<double> >& V_F)
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(F_z_intercell.size()) == d_num_eqn);
    TBOX_ASSERT(static_cast<int>(V_B.size()) == d_num_eqn);
    TBOX_ASSERT(static_cast<int>(V_F.size()) == d_num_eqn);
#endif
    
    // Get the thermodynamic properties of the species.
    std::vector<const double*> thermo_properties_ptr;
    thermo_properties_ptr.reserve(static_cast<int> (d_thermo_properties.size()));
    for (int ti = 0; ti < static_cast<int> (d_thermo_properties.size()); ti++)
    {
        thermo_properties_ptr.push_back(&d_thermo_properties[ti]);
    }
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": RiemannSolverSingleSpeciesHLLC_HLL::"
            << "computeIntercellFluxInZDirectionFromPrimitiveVariables()\n"
            << "There is no z direction for 1D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": RiemannSolverSingleSpeciesHLLC_HLL::"
            << "computeIntercellFluxInZDirectionFromPrimitiveVariables()\n"
            << "There is no z direction for 2D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        double F_z_intercell_HLLC[d_num_eqn];
        double F_z_intercell_HLL[d_num_eqn];
        
        const double c_B = d_equation_of_state->getSoundSpeed(
            &(V_B[0].get()),
            &(V_B[4].get()),
            thermo_properties_ptr);
        
        const double c_F = d_equation_of_state->getSoundSpeed(
            &(V_F[0].get()),
            &(V_F[4].get()),
            thermo_properties_ptr);
        
        const double w_average = 0.5*((V_B[3].get()) + (V_F[3].get()));
        const double c_average = 0.5*(c_B + c_F);
        
        const double s_B = fmin(w_average - c_average, (V_B[3].get()) - c_B);
        const double s_F = fmax(w_average + c_average, (V_F[3].get()) + c_F);
        
        const double s_minus = fmin(0.0, s_B);
        const double s_plus  = fmax(0.0, s_F);
        
        const double s_star =
            ((V_F[4].get()) - (V_B[4].get()) + (V_B[0].get())*(V_B[3].get())*
                (s_B - (V_B[3].get())) -(V_F[0].get())*(V_F[3].get())*(s_F - (V_F[3].get())))/
                    ((V_B[0].get())*(s_B - (V_B[3].get())) - (V_F[0].get())*(s_F - (V_F[3].get())));
        
        if (s_star > 0)
        {
            /*
             * Compute the HLLC flux.
             */
            
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
                &(V_B[4].get()),
                thermo_properties_ptr);
            
            double Q_star_B[d_num_eqn];
            Q_star_B[0] = Chi_star_B*(V_B[0].get());
            Q_star_B[1] = Chi_star_B*Q_B[1];
            Q_star_B[2] = Chi_star_B*Q_B[2];
            Q_star_B[3] = Chi_star_B*(V_B[0].get())*s_star;
            Q_star_B[4] = Chi_star_B*(Q_B[4] + (s_star - (V_B[3].get()))*((V_B[0].get())*s_star +
                (V_B[4].get())/(s_B - (V_B[3].get()))));
            
            double F_z_B[d_num_eqn];
            F_z_B[0] = Q_B[3];
            F_z_B[1] = Q_B[3]*(V_B[1].get());
            F_z_B[2] = Q_B[3]*(V_B[2].get());
            F_z_B[3] = Q_B[3]*(V_B[3].get()) + (V_B[4].get());
            F_z_B[4] = (V_B[3].get())*(Q_B[4] + (V_B[4].get()));
            
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
                std::vector<const double*> vel_F;
                vel_F.reserve(3);
                vel_F.push_back(&(V_F[1].get()));
                vel_F.push_back(&(V_F[2].get()));
                vel_F.push_back(&(V_F[3].get()));
                
                double E_F = d_equation_of_state->getTotalEnergy(
                    &(V_F[0].get()),
                    vel_F,
                    &(V_F[4].get()),
                    thermo_properties_ptr);
                
                if (s_F <= 0)
                {
                    F_z_intercell_HLL[0] = (V_F[0].get())*(V_F[3].get());
                    F_z_intercell_HLL[1] = F_z_intercell_HLL[0]*(V_F[1].get());
                    F_z_intercell_HLL[2] = F_z_intercell_HLL[0]*(V_F[2].get());
                    F_z_intercell_HLL[3] = F_z_intercell_HLL[0]*(V_F[3].get()) + (V_F[4].get());
                    F_z_intercell_HLL[4] = (V_F[2].get())*(E_F + (V_F[4].get()));
                }
                else
                {
                    double Q_F[d_num_eqn];
                    Q_F[0] = (V_F[0].get());
                    Q_F[1] = (V_F[0].get())*(V_F[1].get());
                    Q_F[2] = (V_F[0].get())*(V_F[2].get());
                    Q_F[3] = (V_F[0].get())*(V_F[3].get());
                    Q_F[4] = E_F;
                    
                    double F_z_F[d_num_eqn];
                    F_z_F[0] = Q_F[3];
                    F_z_F[1] = Q_F[3]*(V_F[1].get());
                    F_z_F[2] = Q_F[3]*(V_F[2].get());
                    F_z_F[3] = Q_F[3]*(V_F[3].get()) + (V_F[4].get());
                    F_z_F[4] = (V_F[3].get())*(Q_F[4] + (V_F[4].get()));
                    
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
                &(V_F[4].get()),
                thermo_properties_ptr);
            
            double Q_star_F[d_num_eqn];
            Q_star_F[0] = Chi_star_F*(V_F[0].get());
            Q_star_F[1] = Chi_star_F*Q_F[1];
            Q_star_F[2] = Chi_star_F*Q_F[2];
            Q_star_F[3] = Chi_star_F*(V_F[0].get())*s_star;
            Q_star_F[4] = Chi_star_F*(Q_F[4] + (s_star - (V_F[3].get()))*((V_F[0].get())*s_star +
                (V_F[4].get())/(s_F - (V_F[3].get()))));
            
            double F_z_F[d_num_eqn];
            F_z_F[0] = Q_F[3];
            F_z_F[1] = Q_F[3]*(V_F[1].get());
            F_z_F[2] = Q_F[3]*(V_F[2].get());
            F_z_F[3] = Q_F[3]*(V_F[3].get()) + (V_F[4].get());
            F_z_F[4] = (V_F[3].get())*(Q_F[4] + (V_F[4].get()));
            
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
                std::vector<const double*> vel_B;
                vel_B.reserve(3);
                vel_B.push_back(&(V_B[1].get()));
                vel_B.push_back(&(V_B[2].get()));
                vel_B.push_back(&(V_B[3].get()));
                
                double E_B = d_equation_of_state->getTotalEnergy(
                    &(V_B[0].get()),
                    vel_B,
                    &(V_B[4].get()),
                    thermo_properties_ptr);
                
                if (s_B >= 0)
                {
                    F_z_intercell_HLL[0] = (V_B[0].get())*(V_B[3].get());
                    F_z_intercell_HLL[1] = F_z_intercell_HLL[0]*(V_B[1].get());
                    F_z_intercell_HLL[2] = F_z_intercell_HLL[0]*(V_B[2].get());
                    F_z_intercell_HLL[3] = F_z_intercell_HLL[0]*(V_B[3].get()) + (V_B[4].get());
                    F_z_intercell_HLL[4] = (V_B[2].get())*(E_B + (V_B[4].get()));
                }
                else
                {
                    double Q_B[d_num_eqn];
                    Q_B[0] = (V_B[0].get());
                    Q_B[1] = (V_B[0].get())*(V_B[1].get());
                    Q_B[2] = (V_B[0].get())*(V_B[2].get());
                    Q_B[3] = (V_B[0].get())*(V_B[3].get());
                    Q_B[4] = d_equation_of_state->getTotalEnergy(
                        &(V_B[0].get()),
                        vel_B,
                        &(V_B[4].get()),
                        thermo_properties_ptr);
                    
                    double F_z_B[d_num_eqn];
                    F_z_B[0] = Q_B[3];
                    F_z_B[1] = Q_B[3]*(V_B[1].get());
                    F_z_B[2] = Q_B[3]*(V_B[2].get());
                    F_z_B[3] = Q_B[3]*(V_B[3].get()) + (V_B[4].get());
                    F_z_B[4] = (V_B[3].get())*(Q_B[4] + (V_B[4].get()));
                    
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
        
        double vel_mag = sqrt(pow((V_F[1].get()) - (V_B[1].get()), 2) + pow((V_F[2].get()) -
            (V_B[2].get()), 2) + pow((V_F[3].get()) - (V_B[3].get()), 2));
        
        if (vel_mag < EPSILON)
        {
            alpha1 = 1.0;
            alpha2 = 0.0;
        }
        else
        {
            alpha1 = fabs((V_F[3].get()) - (V_B[3].get()))/vel_mag;
            alpha2 = sqrt(1.0 - alpha1*alpha1);
        }
        
        double beta1 = 0.5 + 0.5*alpha1/(alpha1 + alpha2);
        double beta2 = 1.0 - beta1;
        
        /*
         * Compute the HLLC-HLL flux.
         */
        
        F_z_intercell[0].get() = beta1*F_z_intercell_HLLC[0] + beta2*F_z_intercell_HLL[0];
        F_z_intercell[1].get() = beta1*F_z_intercell_HLLC[1] + beta2*F_z_intercell_HLL[1];
        F_z_intercell[2].get() = beta1*F_z_intercell_HLLC[2] + beta2*F_z_intercell_HLL[2];
        F_z_intercell[3].get() = F_z_intercell_HLLC[3];
        F_z_intercell[4].get() = F_z_intercell_HLLC[4];
    }
}
