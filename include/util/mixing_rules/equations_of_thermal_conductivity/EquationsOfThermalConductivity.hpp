#ifndef EQUATIONS_OF_THERMAL_CONDUCTIVITY_HPP
#define EQUATIONS_OF_THERMAL_CONDUCTIVITY_HPP

#include "util/mixing_rules/equations_of_thermal_conductivity/constant/EquationOfThermalConductivityMixingRulesConstant.hpp"
#include "util/mixing_rules/equations_of_thermal_conductivity/Prandtl/EquationOfThermalConductivityMixingRulesPrandtl.hpp"

#include <map>
#include <string>

namespace EQN_THERMAL_CONDUCTIVITY
{
    enum TYPE { CONSTANT,
                PRANDTL };
}

/*
 * Function to print out enum EQN_THERMAL_CONDUCTIVITY::TYPE value as text.
 */
inline std::ostream& operator<<(std::ostream& os, const EQN_THERMAL_CONDUCTIVITY::TYPE& value)
{
    static std::map<EQN_THERMAL_CONDUCTIVITY::TYPE, std::string> strings;
    
    if (strings.size() == 0)
    {
#define INSERT_ELEMENT(p) strings[p] = #p
        INSERT_ELEMENT(EQN_THERMAL_CONDUCTIVITY::CONSTANT);
        INSERT_ELEMENT(EQN_THERMAL_CONDUCTIVITY::PRANDTL);
#undef INSERT_ELEMENT
    }
    
    return os << strings[value];
}

#endif /* EQUATIONS_OF_THERMAL_CONDUCTIVITY_HPP */
