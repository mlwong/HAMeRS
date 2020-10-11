#ifndef EQUATION_OF_SHEAR_VISCOSITY_MIXING_RULES_MANAGER_HPP
#define EQUATION_OF_SHEAR_VISCOSITY_MIXING_RULES_MANAGER_HPP

#include "HAMeRS_config.hpp"

#include "HAMeRS_memory.hpp"

#include "util/mixing_rules/equations_of_shear_viscosity/EquationsOfShearViscosity.hpp"

#include "SAMRAI/tbox/Dimension.h"

#include <string>

using namespace SAMRAI;

class EquationOfShearViscosityMixingRulesManager
{
    public:
        EquationOfShearViscosityMixingRulesManager(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const int& num_species,
            const MIXING_CLOSURE_MODEL::TYPE& mixing_closure_model,
            const HAMERS_SHARED_PTR<tbox::Database>& equation_of_shear_viscosity_mixing_rules_db,
            const std::string& equation_of_shear_viscosity_str);
        
        /*
         * Get the type of equation of shear viscosity.
         */
        const EQN_SHEAR_VISCOSITY::TYPE&
        getEquationOfShearViscosityType() const
        {
            return d_equation_of_shear_viscosity_type;
        }
        
        /*
         * Get the equation of shear viscosity mixing rules.
         */
        HAMERS_SHARED_PTR<EquationOfShearViscosityMixingRules>
        getEquationOfShearViscosityMixingRules() const
        {
            return d_equation_of_shear_viscosity_mixing_rules;
        }
        
        /*
         * Print all characteristics of equation of shear viscosity manager.
         */
        void
        printClassData(std::ostream& os) const;
        
    private:
        /*
         * The object name is used for error/warning reporting.
         */
        const std::string d_object_name;
        
        /*
         * Type of equation of shear viscosity.
         */
        EQN_SHEAR_VISCOSITY::TYPE d_equation_of_shear_viscosity_type;
        
        /*
         * HAMERS_SHARED_PTR to the equation of shear viscosity.
         */
        HAMERS_SHARED_PTR<EquationOfShearViscosityMixingRules> d_equation_of_shear_viscosity_mixing_rules;
        
};

#endif /* EQUATION_OF_SHEAR_VISCOSITY_MIXING_RULES_MANAGER_HPP */
