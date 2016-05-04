#ifndef RIEMANN_SOLVER_HLLC_HPP
#define RIEMANN_SOLVER_HLLC_HPP

#include "flow_model/Riemann_solver/RiemannSolver.hpp"

class RiemannSolverHLLC: public RiemannSolver
{
    public:
        RiemannSolverHLLC(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const int& num_eqn,
            const int& num_species,
            boost::shared_ptr<EquationOfState> equation_of_state):
                RiemannSolver(
                    object_name,
                    dim,
                    num_eqn,
                    num_species,
                    equation_of_state)
        {}
        
        /*
         * Compute the flux at the intercell face for single-species flow model
         * from conservative variables.
         */
        void
        computeIntercellFluxForSingleSpeciesFromConservativeVariables(
            std::vector<double*>& flux_intercell,
            std::vector<double>& conservative_variables_minus,
            std::vector<double>& conservative_variables_plus,
            DIRECTION direction);
        
        /*
         * Compute the flux at the intercell face for single-species flow model
         * from primitive variables.
         */
        void
        computeIntercellFluxForSingleSpeciesFromPrimitiveVariables(
            std::vector<double*>& flux_intercell,
            std::vector<double>& primitive_variables_minus,
            std::vector<double>& primitive_variables_plus,
            DIRECTION direction);
        
        /*
         * Compute the flux and at the intercell face
         * for single-species flow model.
         */
        void
        computeIntercellFluxForSingleSpecies(
            std::vector<double*>& flux_intercell,
            const double* const density_L,
            const double* const density_R,
            const std::vector<const double*>& momentum_L,
            const std::vector<const double*>& momentum_R,
            const double* const total_energy_L,
            const double* const total_energy_R,
            DIRECTION direction);
        
        /*
         * Compute the flux and velocity at the intercell face
         * for four-equation multi-species flow model by Shyue.
         */
        void
        computeIntercellFluxAndVelocityForFourEqnShyue(
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
            DIRECTION direction);
        
        /*
         * Compute the flux and velocity at the intercell face
         * for five-equation multi-species flow model by Allaire.
         */
        void
        computeIntercellFluxAndVelocityForFiveEqnAllaire(
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
            DIRECTION direction);
        
    private:
        /*
         * Compute the flux in the x-direction at the intercell face
         * for single-species flow model from conservative variables.
         */
        void
        computeIntercellFluxForSingleSpeciesInXDirectionFromConservativeVariables(
            std::vector<double*>& F_x_intercell,
            std::vector<double>& Q_L,
            std::vector<double>& Q_R);
        
        /*
         * Compute the flux in the y-direction at the intercell face
         * for single-species flow model from conservative variables.
         */
        void
        computeIntercellFluxForSingleSpeciesInYDirectionFromConservativeVariables(
            std::vector<double*>& F_y_intercell,
            std::vector<double>& Q_B,
            std::vector<double>& Q_T);
        
        /*
         * Compute the flux in the z-direction at the intercell face
         * for single-species flow model from conservative variables.
         */
        void
        computeIntercellFluxForSingleSpeciesInZDirectionFromConservativeVariables(
            std::vector<double*>& F_z_intercell,
            std::vector<double>& Q_B,
            std::vector<double>& Q_F);
        
        /*
         * Compute the flux in the x-direction at the intercell face
         * for single-species flow model from primitive variables.
         */
        void
        computeIntercellFluxForSingleSpeciesInXDirectionFromPrimitiveVariables(
            std::vector<double*>& F_x_intercell,
            std::vector<double>& V_L,
            std::vector<double>& V_R);
        
        /*
         * Compute the flux in the y-direction at the intercell face
         * for single-species flow model from primitive variables.
         */
        void
        computeIntercellFluxForSingleSpeciesInYDirectionFromPrimitiveVariables(
            std::vector<double*>& F_y_intercell,
            std::vector<double>& V_B,
            std::vector<double>& V_T);
        
        /*
         * Compute the flux in the z-direction at the intercell face
         * for single-species flow model from primitive variables.
         */
        void
        computeIntercellFluxForSingleSpeciesInZDirectionFromPrimitiveVariables(
            std::vector<double*>& F_z_intercell,
            std::vector<double>& V_B,
            std::vector<double>& V_F);
        
};

#endif /* RIEMANN_SOLVER_HLLC_HPP */