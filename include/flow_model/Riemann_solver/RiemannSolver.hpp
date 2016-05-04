#ifndef RIEMANN_SOLVER_HPP
#define RIEMANN_SOLVER_HPP

#include "SAMRAI/tbox/Dimension.h"

#include "Directions.hpp"
#include "equation_of_state/EquationOfStateIdealGas.hpp"

#include "boost/ptr_container/ptr_vector.hpp"
#include "boost/shared_ptr.hpp"
#include <string>
#include <vector>

using namespace SAMRAI;

class RiemannSolver
{
    public:
        RiemannSolver(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const int& num_eqn,
            const int& num_species,
            boost::shared_ptr<EquationOfState> equation_of_state):
                d_object_name(object_name),
                d_dim(dim),
                d_num_eqn(num_eqn),
                d_num_species(num_species),
                d_equation_of_state(equation_of_state)
        {}
        
        /*
         * Compute the flux at the intercell face for single-species flow model
         * from conservative variables.
         */
        virtual void
        computeIntercellFluxForSingleSpeciesFromConservativeVariables(
            std::vector<double*>& flux_intercell,
            std::vector<double>& conservative_variables_minus,
            std::vector<double>& conservative_variables_plus,
            DIRECTION direction) = 0;
        
        /*
         * Compute the flux at the intercell face for single-species flow model
         * from primitive variables.
         */
        virtual void
        computeIntercellFluxForSingleSpeciesFromPrimitiveVariables(
            std::vector<double*>& flux_intercell,
            std::vector<double>& primitive_variables_minus,
            std::vector<double>& primitive_variables_plus,
            DIRECTION direction) = 0;
        
        /*
         * Compute the flux at the intercell face
         * for single-species flow model.
         */
        virtual void
        computeIntercellFluxForSingleSpecies(
            std::vector<double*>& flux_intercell,
            const double* const density_L,
            const double* const density_R,
            const std::vector<const double*>& momentum_L,
            const std::vector<const double*>& momentum_R,
            const double* const total_energy_L,
            const double* const total_energy_R,
            DIRECTION direction) = 0;
        
        /*
         * Compute the flux and velocity at the intercell face
         * for four-equation multi-species flow model by Shyue.
         */
        virtual void
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
            DIRECTION direction) = 0;
        
        /*
         * Compute the flux and velocity at the intercell face
         * for five-equation multi-species flow model by Allaire.
         */
        virtual void
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
            DIRECTION direction) = 0;
        
    protected:
        /*
         * The object name is used for error/warning reporting.
         */
        const std::string d_object_name;
        
        /*
         * Problem dimension.
         */
        const tbox::Dimension d_dim;
        
        /*
         * Number of equations.
         */
        int d_num_eqn;
        
        /*
         * Number of species.
         */
        int d_num_species;
        
        /*
         * boost::shared_ptr to EquationOfState.
         */
        boost::shared_ptr<EquationOfState> d_equation_of_state;
};

#endif /* RIEMANN_SOLVER_HPP */
