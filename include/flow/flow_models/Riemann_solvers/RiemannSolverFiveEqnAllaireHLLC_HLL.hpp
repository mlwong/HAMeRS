#ifndef RIEMANN_SOLVER_FIVE_EQN_ALLAIRE_HLLC_HLL_HPP
#define RIEMANN_SOLVER_FIVE_EQN_ALLAIRE_HLLC_HLL_HPP

#include "flow/flow_models/Riemann_solvers/RiemannSolverFiveEqnAllaire.hpp"

using namespace SAMRAI;

class RiemannSolverFiveEqnAllaireHLLC_HLL: public RiemannSolverFiveEqnAllaire
{
    public:
        RiemannSolverFiveEqnAllaireHLLC_HLL(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const int& num_eqn,
            const int& num_species,
            boost::shared_ptr<EquationOfState> equation_of_state):
                RiemannSolverFiveEqnAllaire(
                    object_name,
                    dim,
                    num_eqn,
                    num_species,
                    equation_of_state)
        {}
        
        /*
         * Compute the flux and velocity at the intercell face from conservative variables.
         */
        void
        computeIntercellFluxAndVelocityFromConservativeVariables(
            std::vector<boost::reference_wrapper<double> >& flux_intercell,
            std::vector<boost::reference_wrapper<double> >& velocity_intercell,
            const std::vector<boost::reference_wrapper<double> >& conservative_variables_minus,
            const std::vector<boost::reference_wrapper<double> >& conservative_variables_plus,
            const DIRECTION& direction);
        
        /*
         * Compute the flux and velocity at the intercell face from primitive variables.
         */
        void
        computeIntercellFluxAndVelocityFromPrimitiveVariables(
            std::vector<boost::reference_wrapper<double> >& flux_intercell,
            std::vector<boost::reference_wrapper<double> >& velocity_intercell,
            const std::vector<boost::reference_wrapper<double> >& primitive_variables_minus,
            const std::vector<boost::reference_wrapper<double> >& primitive_variables_plus,
            const DIRECTION& direction);
        
    private:
        /*
         * Compute the flux and velocity in the x-direction at the intercell face
         * from conservative variables.
         */
        void
        computeIntercellFluxAndVelocityInXDirectionFromConservativeVariables(
            std::vector<boost::reference_wrapper<double> >& F_x_intercell,
            std::vector<boost::reference_wrapper<double> >& vel_intercell,
            const std::vector<boost::reference_wrapper<double> >& Q_L,
            const std::vector<boost::reference_wrapper<double> >& Q_R);
        
        /*
         * Compute the flux and velocity in the y-direction at the intercell face
         * from conservative variables.
         */
        void
        computeIntercellFluxAndVelocityInYDirectionFromConservativeVariables(
            std::vector<boost::reference_wrapper<double> >& F_y_intercell,
            std::vector<boost::reference_wrapper<double> >& vel_intercell,
            const std::vector<boost::reference_wrapper<double> >& Q_B,
            const std::vector<boost::reference_wrapper<double> >& Q_T);
        
        /*
         * Compute the flux and velocity in the z-direction at the intercell face
         * from conservative variables.
         */
        void
        computeIntercellFluxAndVelocityInZDirectionFromConservativeVariables(
            std::vector<boost::reference_wrapper<double> >& F_z_intercell,
            std::vector<boost::reference_wrapper<double> >& vel_intercell,
            const std::vector<boost::reference_wrapper<double> >& Q_B,
            const std::vector<boost::reference_wrapper<double> >& Q_F);
        
        /*
         * Compute the flux and velocity in the x-direction at the intercell face
         * from primitive variables.
         */
        void
        computeIntercellFluxAndVelocityInXDirectionFromPrimitiveVariables(
            std::vector<boost::reference_wrapper<double> >& F_x_intercell,
            std::vector<boost::reference_wrapper<double> >& vel_intercell,
            const std::vector<boost::reference_wrapper<double> >& V_L,
            const std::vector<boost::reference_wrapper<double> >& V_R);
        
        /*
         * Compute the flux and velocity in the y-direction at the intercell face
         * from primitive variables.
         */
        void
        computeIntercellFluxAndVelocityInYDirectionFromPrimitiveVariables(
            std::vector<boost::reference_wrapper<double> >& F_y_intercell,
            std::vector<boost::reference_wrapper<double> >& vel_intercell,
            const std::vector<boost::reference_wrapper<double> >& V_B,
            const std::vector<boost::reference_wrapper<double> >& V_T);
        
        /*
         * Compute the flux and velocity in the z-direction at the intercell face
         * from primitive variables.
         */
        void
        computeIntercellFluxAndVelocityInZDirectionFromPrimitiveVariables(
            std::vector<boost::reference_wrapper<double> >& F_z_intercell,
            std::vector<boost::reference_wrapper<double> >& vel_intercell,
            const std::vector<boost::reference_wrapper<double> >& V_B,
            const std::vector<boost::reference_wrapper<double> >& V_F);
        
};

#endif /* RIEMANN_SOLVER_FIVE_EQN_ALLAIRE_HLLC_HLL_HPP */
