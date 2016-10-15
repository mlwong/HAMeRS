#ifndef EQUATIONS_OF_SHEAR_VISCOSITY_HPP
#define EQUATIONS_OF_SHEAR_VISCOSITY_HPP

#include "util/mixing_rules/equations_of_shear_viscosity/constant/EquationOfShearViscosityMixingRulesConstant.hpp"
#include "util/mixing_rules/equations_of_shear_viscosity/Chapman-Enskog/EquationOfShearViscosityMixingRulesChapmanEnskog.hpp"

#include <map>
#include <string>

enum EQUATION_OF_SHEAR_VISCOSITY_LABEL { CONSTANT,
                                         CHAPMAN_ENSKOG };

/*
 * Function to print out enum EQUATION_OF_SHEAR_VISCOSITY_LABEL value as text.
 */
inline std::ostream& operator<<(std::ostream& os, const EQUATION_OF_SHEAR_VISCOSITY_LABEL& value)
{
    static std::map<EQUATION_OF_SHEAR_VISCOSITY_LABEL, std::string> strings;
    
    if (strings.size() == 0)
    {
#define INSERT_ELEMENT(p) strings[p] = #p
        INSERT_ELEMENT(CONSTANT);
        INSERT_ELEMENT(CHAPMAN_ENSKOG);
#undef INSERT_ELEMENT
    }
    
    return os << strings[value];
}

#endif /* EQUATIONS_OF_SHEAR_VISCOSITY_HPP */