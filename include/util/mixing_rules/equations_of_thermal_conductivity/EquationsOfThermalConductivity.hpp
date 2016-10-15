#ifndef EQUATIONS_OF_THERMAL_CONDUCTIVITY_HPP
#define EQUATIONS_OF_THERMAL_CONDUCTIVITY_HPP

#include "util/mixing_rules/equations_of_thermal_conductivity/Prandtl/EquationOfThermalConductivityMixingRulesPrandtl.hpp"

#include <map>
#include <string>

enum EQUATION_OF_THERMAL_CONDUCTIVITY_LABEL { PRANDTL };

/*
 * Function to print out enum EQUATION_OF_THERMAL_CONDUCTIVITY_LABEL value as text.
 */
inline std::ostream& operator<<(std::ostream& os, const EQUATION_OF_THERMAL_CONDUCTIVITY_LABEL& value)
{
    static std::map<EQUATION_OF_THERMAL_CONDUCTIVITY_LABEL, std::string> strings;
    
    if (strings.size() == 0)
    {
#define INSERT_ELEMENT(p) strings[p] = #p
        INSERT_ELEMENT(PRANDTL);
#undef INSERT_ELEMENT
    }
    
    return os << strings[value];
}

#endif /* EQUATIONS_OF_THERMAL_CONDUCTIVITY_HPP */