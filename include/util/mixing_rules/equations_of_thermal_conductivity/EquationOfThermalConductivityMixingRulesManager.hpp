#ifndef EQUATION_OF_THERMAL_CONDUCTIVITY_MIXING_RULES_MANAGER_HPP
#define EQUATION_OF_THERMAL_CONDUCTIVITY_MIXING_RULES_MANAGER_HPP

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/Dimension.h"

#include "util/mixing_rules/equations_of_thermal_conductivity/EquationsOfThermalConductivity.hpp"

#include "boost/shared_ptr.hpp"
#include <string>

using namespace SAMRAI;

class EquationOfThermalConductivityMixingRulesManager
{
    public:
        EquationOfThermalConductivityMixingRulesManager(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const int& num_species,
            const MIXING_CLOSURE_MODEL& mixing_closure_model,
            const boost::shared_ptr<tbox::Database>& equation_of_thermal_conductivity_mixing_rules_db,
            const std::string& equation_of_thermal_conductivity_str);
        
        /*
         * Get the label of equation of thermal conductivity.
         */
        const EQUATION_OF_THERMAL_CONDUCTIVITY_LABEL&
        getEquationOfThermalConductivityLabel() const
        {
            return d_equation_of_thermal_conductivity_label;
        }
        
        /*
         * Get the equation of thermal conductivity mixing rules.
         */
        boost::shared_ptr<EquationOfThermalConductivityMixingRules>
        getEquationOfThermalConductivityMixingRules() const
        {
            return d_equation_of_thermal_conductivity_mixing_rules;
        }
        
        /*
         * Print all characteristics of equation of thermal conductivity manager.
         */
        void
        printClassData(std::ostream& os) const;
        
    private:
        /*
         * The object name is used for error/warning reporting.
         */
        const std::string d_object_name;
        
        /*
         * Label of equation of thermal conductivity.
         */
        EQUATION_OF_THERMAL_CONDUCTIVITY_LABEL d_equation_of_thermal_conductivity_label;
        
        /*
         * boost::shared_ptr to the equation of thermal conductivity.
         */
        boost::shared_ptr<EquationOfThermalConductivityMixingRules> d_equation_of_thermal_conductivity_mixing_rules;
        
};

#endif /* EQUATION_OF_THERMAL_CONDUCTIVITY_MIXING_RULES_MANAGER_HPP */