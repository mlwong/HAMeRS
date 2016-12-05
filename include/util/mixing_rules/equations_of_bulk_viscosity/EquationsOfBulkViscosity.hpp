#ifndef EQUATIONS_OF_BULK_VISCOSITY_HPP
#define EQUATIONS_OF_BULK_VISCOSITY_HPP

#include "util/mixing_rules/equations_of_bulk_viscosity/constant/EquationOfBulkViscosityMixingRulesConstant.hpp"
#include "util/mixing_rules/equations_of_bulk_viscosity/Cramer/EquationOfBulkViscosityMixingRulesCramer.hpp"

#include <map>
#include <string>

namespace EQN_BULK_VISCOSITY
{
    enum TYPE { CONSTANT,
                CRAMER };
}

/*
 * Function to print out enum EQN_BULK_VISCOSITY::TYPE value as text.
 */
inline std::ostream& operator<<(std::ostream& os, const EQN_BULK_VISCOSITY::TYPE& value)
{
    static std::map<EQN_BULK_VISCOSITY::TYPE, std::string> strings;
    
    if (strings.size() == 0)
    {
#define INSERT_ELEMENT(p) strings[p] = #p
        INSERT_ELEMENT(EQN_BULK_VISCOSITY::CONSTANT);
        INSERT_ELEMENT(EQN_BULK_VISCOSITY::CRAMER);
#undef INSERT_ELEMENT
    }
    
    return os << strings[value];
}

#endif /* EQUATIONS_OF_BULK_VISCOSITY_HPP */
