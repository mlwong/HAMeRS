#ifndef RIEMANN_SOLVER_SINGLE_SPECIES_HPP
#define RIEMANN_SOLVER_SINGLE_SPECIES_HPP

#include "SAMRAI/tbox/Dimension.h"

#include "Directions.hpp"
#include "util/equations_of_state/EquationsOfState.hpp"

#include "boost/ref.hpp"
#include "boost/shared_ptr.hpp"
#include <string>
#include <vector>

using namespace SAMRAI;

class RiemannSolverSingleSpecies
{
    public:
        RiemannSolverSingleSpecies(
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
         * Compute the flux at the intercell face from conservative variables.
         */
        virtual void
        computeIntercellFluxFromConservativeVariables(
            std::vector<boost::reference_wrapper<double> >& flux_intercell,
            const std::vector<boost::reference_wrapper<double> >& conservative_variables_minus,
            const std::vector<boost::reference_wrapper<double> >& conservative_variables_plus,
            const DIRECTION& direction) = 0;
        
        /*
         * Compute the flux at the intercell face from primitive variables.
         */
        virtual void
        computeIntercellFluxFromPrimitiveVariables(
            std::vector<boost::reference_wrapper<double> >& flux_intercell,
            const std::vector<boost::reference_wrapper<double> >& primitive_variables_minus,
            const std::vector<boost::reference_wrapper<double> >& primitive_variables_plus,
            const DIRECTION& direction) = 0;
        
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

#endif /* RIEMANN_SOLVER_SINGLE_SPECIES_HPP */
