#ifndef RIEMANN_SOLVER_FOUR_EQN_CONSERVATIVE_HPP
#define RIEMANN_SOLVER_FOUR_EQN_CONSERVATIVE_HPP

#include "SAMRAI/tbox/Dimension.h"

#include "util/Directions.hpp"
#include "util/mixing_rules/equations_of_state/EquationOfStateMixingRules.hpp"

#include "boost/ref.hpp"
#include "boost/shared_ptr.hpp"
#include <string>
#include <vector>

using namespace SAMRAI;

class RiemannSolverFourEqnConservative
{
    public:
        RiemannSolverFourEqnConservative(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const int& num_eqn,
            const int& num_species,
            boost::shared_ptr<EquationOfStateMixingRules> equation_of_state_mixing_rules):
                d_object_name(object_name),
                d_dim(dim),
                d_num_eqn(num_eqn),
                d_num_species(num_species),
                d_equation_of_state_mixing_rules(equation_of_state_mixing_rules)
        {}
        
        /*
         * Compute the flux at the intercell face from conservative variables.
         */
        virtual void
        computeIntercellFluxFromConservativeVariables(
            std::vector<boost::reference_wrapper<double> >& flux_intercell,
            const std::vector<boost::reference_wrapper<double> >& conservative_variables_minus,
            const std::vector<boost::reference_wrapper<double> >& conservative_variables_plus,
            const DIRECTION::TYPE& direction) = 0;
        
        /*
         * Compute the flux at the intercell face from primitive variables.
         */
        virtual void
        computeIntercellFluxFromPrimitiveVariables(
            std::vector<boost::reference_wrapper<double> >& flux_intercell,
            const std::vector<boost::reference_wrapper<double> >& primitive_variables_minus,
            const std::vector<boost::reference_wrapper<double> >& primitive_variables_plus,
            const DIRECTION::TYPE& direction) = 0;
        
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
         * boost::shared_ptr to EquationOfStateMixingRules.
         */
        const boost::shared_ptr<EquationOfStateMixingRules> d_equation_of_state_mixing_rules;
        
};

#endif /* RIEMANN_SOLVER_FOUR_EQN_CONSERVATIVE_HPP */
