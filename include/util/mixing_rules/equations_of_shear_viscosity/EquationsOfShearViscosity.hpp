#ifndef EQUATIONS_OF_SHEAR_VISCOSITY_HPP
#define EQUATIONS_OF_SHEAR_VISCOSITY_HPP

#include "util/mixing_rules/equations_of_shear_viscosity/constant/EquationOfShearViscosityMixingRulesConstant.hpp"
#include "util/mixing_rules/equations_of_shear_viscosity/Chapman-Enskog/EquationOfShearViscosityMixingRulesChapmanEnskog.hpp"

#include <map>
#include <string>

namespace EQN_SHEAR_VISCOSITY
{
    enum TYPE { CONSTANT,
                CHAPMAN_ENSKOG };
}

/*
 * Function to print out enum EQN_SHEAR_VISCOSITY::TYPE value as text.
 */
inline std::ostream& operator<<(std::ostream& os, const EQN_SHEAR_VISCOSITY::TYPE& value)
{
    static std::map<EQN_SHEAR_VISCOSITY::TYPE, std::string> strings;
    
    if (strings.size() == 0)
    {
#define INSERT_ELEMENT(p) strings[p] = #p
        INSERT_ELEMENT(EQN_SHEAR_VISCOSITY::CONSTANT);
        INSERT_ELEMENT(EQN_SHEAR_VISCOSITY::CHAPMAN_ENSKOG);
#undef INSERT_ELEMENT
    }
    
    return os << strings[value];
}

#endif /* EQUATIONS_OF_SHEAR_VISCOSITY_HPP */
