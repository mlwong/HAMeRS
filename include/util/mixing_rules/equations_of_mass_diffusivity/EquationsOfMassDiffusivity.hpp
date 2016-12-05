#ifndef EQUATIONS_OF_MASS_DIFFUSIVITY_HPP
#define EQUATIONS_OF_MASS_DIFFUSIVITY_HPP

#include "util/mixing_rules/equations_of_mass_diffusivity/constant/EquationOfMassDiffusivityMixingRulesConstant.hpp"
#include "util/mixing_rules/equations_of_mass_diffusivity/Reid/EquationOfMassDiffusivityMixingRulesReid.hpp"

#include <map>
#include <string>

namespace EQN_MASS_DIFFUSIVITY
{
    enum TYPE { CONSTANT,
                REID };
}

/*
 * Function to print out enum EQN_MASS_DIFFUSIVITY::TYPE value as text.
 */
inline std::ostream& operator<<(std::ostream& os, const EQN_MASS_DIFFUSIVITY::TYPE& value)
{
    static std::map<EQN_MASS_DIFFUSIVITY::TYPE, std::string> strings;
    
    if (strings.size() == 0)
    {
#define INSERT_ELEMENT(p) strings[p] = #p
        INSERT_ELEMENT(EQN_MASS_DIFFUSIVITY::CONSTANT);
        INSERT_ELEMENT(EQN_MASS_DIFFUSIVITY::REID);
#undef INSERT_ELEMENT
    }
    
    return os << strings[value];
}

#endif /* EQUATIONS_OF_MASS_DIFFUSIVITY_HPP */
