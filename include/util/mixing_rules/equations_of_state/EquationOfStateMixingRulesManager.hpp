#ifndef EQUATION_OF_STATE_MIXING_RULES_MANAGER_HPP
#define EQUATION_OF_STATE_MIXING_RULES_MANAGER_HPP

#include "HAMeRS_config.hpp"

#include "util/mixing_rules/equations_of_state/EquationsOfState.hpp"

#include "SAMRAI/tbox/Dimension.h"

#include "boost/shared_ptr.hpp"
#include <string>

using namespace SAMRAI;

class EquationOfStateMixingRulesManager
{
    public:
        EquationOfStateMixingRulesManager(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const int& num_species,
            const MIXING_CLOSURE_MODEL::TYPE& mixing_closure_model,
            const boost::shared_ptr<tbox::Database>& equation_of_state_mixing_rules_db,
            const std::string& equation_of_state_str);
        
        /*
         * Get the type of equation of state.
         */
        const EQN_STATE::TYPE&
        getEquationOfStateType() const
        {
            return d_equation_of_state_type;
        }
        
        /*
         * Get the equation of state mixing rules.
         */
        boost::shared_ptr<EquationOfStateMixingRules>
        getEquationOfStateMixingRules() const
        {
            return d_equation_of_state_mixing_rules;
        }
        
        /*
         * Print all characteristics of equation of state manager.
         */
        void
        printClassData(std::ostream& os) const;
        
    private:
        /*
         * The object name is used for error/warning reporting.
         */
        const std::string d_object_name;
        
        /*
         * Type of equation of state.
         */
        EQN_STATE::TYPE d_equation_of_state_type;
        
        /*
         * boost::shared_ptr to the equation of state.
         */
        boost::shared_ptr<EquationOfStateMixingRules> d_equation_of_state_mixing_rules;
        
};

#endif /* EQUATION_OF_STATE_MIXING_RULES_MANAGER_HPP */
