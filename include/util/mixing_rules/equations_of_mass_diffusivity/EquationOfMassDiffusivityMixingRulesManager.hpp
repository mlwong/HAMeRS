#ifndef EQUATION_OF_MASS_DIFFUSIVITY_MIXING_RULES_MANAGER_HPP
#define EQUATION_OF_MASS_DIFFUSIVITY_MIXING_RULES_MANAGER_HPP

#include "HAMeRS_config.hpp"

#include "HAMeRS_memory.hpp"

#include "util/mixing_rules/equations_of_mass_diffusivity/EquationsOfMassDiffusivity.hpp"

#include "SAMRAI/tbox/Dimension.h"

#include <string>

using namespace SAMRAI;

class EquationOfMassDiffusivityMixingRulesManager
{
    public:
        EquationOfMassDiffusivityMixingRulesManager(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const int& num_species,
            const MIXING_CLOSURE_MODEL::TYPE& mixing_closure_model,
            const boost::shared_ptr<tbox::Database>& equation_of_mass_diffusivity_mixing_rules_db,
            const std::string& equation_of_mass_diffusivity_str);
        
        /*
         * Get the type of equation of mass diffusivity.
         */
        const EQN_MASS_DIFFUSIVITY::TYPE&
        getEquationOfMassDiffusivityType() const
        {
            return d_equation_of_mass_diffusivity_type;
        }
        
        /*
         * Get the equation of mass diffusivity mixing rules.
         */
        boost::shared_ptr<EquationOfMassDiffusivityMixingRules>
        getEquationOfMassDiffusivityMixingRules() const
        {
            return d_equation_of_mass_diffusivity_mixing_rules;
        }
        
        /*
         * Print all characteristics of equation of mass diffusivity manager.
         */
        void
        printClassData(std::ostream& os) const;
        
    private:
        /*
         * The object name is used for error/warning reporting.
         */
        const std::string d_object_name;
        
        /*
         * Type of equation of mass diffusivity.
         */
        EQN_MASS_DIFFUSIVITY::TYPE d_equation_of_mass_diffusivity_type;
        
        /*
         * boost::shared_ptr to the equation of mass diffusivity.
         */
        boost::shared_ptr<EquationOfMassDiffusivityMixingRules> d_equation_of_mass_diffusivity_mixing_rules;
        
};

#endif /* EQUATION_OF_MASS_DIFFUSIVITY_MIXING_RULES_MANAGER_HPP */
