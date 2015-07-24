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
         * Compute the fluxes and and velocities at the intercell faces
         * for single-species flow model.
         */
        void
        computeIntercellFluxForSingleSpecies(
            std::vector<double*> flux_intercell,
            const double* const density_L,
            const double* const density_R,
            const std::vector<const double*> momentum_L,
            const std::vector<const double*> momentum_R,
            const double* const total_energy_L,
            const double* const total_energy_R,
            DIRECTION direction);
        
        /*
         * Compute the fluxes and and velocities at the intercell faces
         * for four-equation multi-species flow model by Shyue.
         */
        void
        computeIntercellFluxAndVelocityForFourEqnShyue(
            std::vector<double*> flux_intercell,
            std::vector<double*> velocity_intercell,
            const double* const density_L,
            const double* const density_R,
            const std::vector<const double*> momentum_L,
            const std::vector<const double*> momentum_R,
            const double* const total_energy_L,
            const double* const total_energy_R,
            const std::vector<const double*> mass_fraction_L,
            const std::vector<const double*> mass_fraction_R,
            DIRECTION direction);
        
        /*
         * Compute the fluxes and and velocities at the intercell faces
         * for five-equation multi-species flow model by Allaire.
         */
        void
        computeIntercellFluxAndVelocityForFiveEqnAllaire(
            std::vector<double*> flux_intercell,
            std::vector<double*> velocity_intercell,
            const std::vector<const double*> partial_density_L,
            const std::vector<const double*> partial_density_R,
            const std::vector<const double*> momentum_L,
            const std::vector<const double*> momentum_R,
            const double* const total_energy_L,
            const double* const total_energy_R,
            const std::vector<const double*> volume_fraction_L,
            const std::vector<const double*> volume_fraction_R,
            DIRECTION direction);
        
};

#endif /* RIEMANN_SOLVER_HLLC_HPP */