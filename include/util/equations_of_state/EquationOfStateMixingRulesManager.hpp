#ifndef EQUATION_OF_STATE_MIXING_RULES_MANAGER_HPP
#define EQUATION_OF_STATE_MIXING_RULES_MANAGER_HPP

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/Dimension.h"

#include "util/equations_of_state/EquationsOfState.hpp"

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
            const MIXING_CLOSURE_MODEL& mixing_closure_model,
            const boost::shared_ptr<tbox::Database>& species_db,
            const std::string& equation_of_state_str);
        
        /*
         * Get the label of equation of state.
         */
        const EQUATION_OF_STATE_LABEL&
        getEquationOfStateLabel() const
        {
            return d_equation_of_state_label;
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
         * Label of equation of state.
         */
        EQUATION_OF_STATE_LABEL d_equation_of_state_label;
        
        /*
         * boost::shared_ptr to the equation of state.
         */
        boost::shared_ptr<EquationOfStateMixingRules> d_equation_of_state_mixing_rules;
        
};

#endif /* EQUATION_OF_STATE_MIXING_RULES_MANAGER_HPP */
