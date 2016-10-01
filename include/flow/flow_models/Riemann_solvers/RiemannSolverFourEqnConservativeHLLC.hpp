#ifndef RIEMANN_SOLVER_FOUR_EQN_CONSERVATIVE_HLLC_HPP
#define RIEMANN_SOLVER_FOUR_EQN_CONSERVATIVE_HLLC_HPP

#include "flow/flow_models/Riemann_solvers/RiemannSolverFourEqnConservative.hpp"

using namespace SAMRAI;

class RiemannSolverFourEqnConservativeHLLC: public RiemannSolverFourEqnConservative
{
    public:
        RiemannSolverFourEqnConservativeHLLC(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const int& num_eqn,
            const int& num_species,
            boost::shared_ptr<EquationOfState> equation_of_state,
            boost::shared_ptr<EquationOfStateMixingRules> equation_of_state_mixing_rules):
                RiemannSolverFourEqnConservative(
                    object_name,
                    dim,
                    num_eqn,
                    num_species,
                    equation_of_state,
                    equation_of_state_mixing_rules)
        {}
        
        /*
         * Compute the flux at the intercell face from conservative variables.
         */
        void
        computeIntercellFluxFromConservativeVariables(
            std::vector<boost::reference_wrapper<double> >& flux_intercell,
            const std::vector<boost::reference_wrapper<double> >& conservative_variables_minus,
            const std::vector<boost::reference_wrapper<double> >& conservative_variables_plus,
            const DIRECTION& direction);
        
        /*
         * Compute the flux at the intercell face from primitive variables.
         */
        void
        computeIntercellFluxFromPrimitiveVariables(
            std::vector<boost::reference_wrapper<double> >& flux_intercell,
            const std::vector<boost::reference_wrapper<double> >& primitive_variables_minus,
            const std::vector<boost::reference_wrapper<double> >& primitive_variables_plus,
            const DIRECTION& direction);
        
    private:
        /*
         * Compute the flux in the x-direction at the intercell face
         * from conservative variables.
         */
        void
        computeIntercellFluxInXDirectionFromConservativeVariables(
            std::vector<boost::reference_wrapper<double> >& F_x_intercell,
            const std::vector<boost::reference_wrapper<double> >& Q_L,
            const std::vector<boost::reference_wrapper<double> >& Q_R);
        
        /*
         * Compute the flux in the y-direction at the intercell face
         * from conservative variables.
         */
        void
        computeIntercellFluxInYDirectionFromConservativeVariables(
            std::vector<boost::reference_wrapper<double> >& F_y_intercell,
            const std::vector<boost::reference_wrapper<double> >& Q_B,
            const std::vector<boost::reference_wrapper<double> >& Q_T);
        
        /*
         * Compute the flux in the z-direction at the intercell face
         * from conservative variables.
         */
        void
        computeIntercellFluxInZDirectionFromConservativeVariables(
            std::vector<boost::reference_wrapper<double> >& F_z_intercell,
            const std::vector<boost::reference_wrapper<double> >& Q_B,
            const std::vector<boost::reference_wrapper<double> >& Q_F);
        
        /*
         * Compute the flux in the x-direction at the intercell face
         * from primitive variables.
         */
        void
        computeIntercellFluxInXDirectionFromPrimitiveVariables(
            std::vector<boost::reference_wrapper<double> >& F_x_intercell,
            const std::vector<boost::reference_wrapper<double> >& V_L,
            const std::vector<boost::reference_wrapper<double> >& V_R);
        
        /*
         * Compute the flux in the y-direction at the intercell face
         * from primitive variables.
         */
        void
        computeIntercellFluxInYDirectionFromPrimitiveVariables(
            std::vector<boost::reference_wrapper<double> >& F_y_intercell,
            const std::vector<boost::reference_wrapper<double> >& V_B,
            const std::vector<boost::reference_wrapper<double> >& V_T);
        
        /*
         * Compute the flux in the z-direction at the intercell face
         * from primitive variables.
         */
        void
        computeIntercellFluxInZDirectionFromPrimitiveVariables(
            std::vector<boost::reference_wrapper<double> >& F_z_intercell,
            const std::vector<boost::reference_wrapper<double> >& V_B,
            const std::vector<boost::reference_wrapper<double> >& V_F);
        
};

#endif /* RIEMANN_SOLVER_FOUR_EQN_CONSERVATIVE_HLLC_HPP */
