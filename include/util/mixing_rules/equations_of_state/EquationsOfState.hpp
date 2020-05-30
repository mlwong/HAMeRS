#ifndef EQUATIONS_OF_STATE_HPP
#define EQUATIONS_OF_STATE_HPP

#include "util/mixing_rules/equations_of_state/ideal_gas/EquationOfStateIdealGas.hpp"
#include "util/mixing_rules/equations_of_state/ideal_gas/EquationOfStateMixingRulesIdealGas.hpp"
#include "util/mixing_rules/equations_of_state/stiffened_gas/EquationOfStateMixingRulesStiffenedGas.hpp"
#include "util/mixing_rules/equations_of_state/stiffened_gas/EquationOfStateMixingRulesStiffenedGas.hpp"

#include <map>
#include <string>

namespace EQN_STATE
{
    enum TYPE { IDEAL_GAS,
                STIFFENED_GAS };
}

/*
 * Function to print out enum EQN_STATE::TYPE value as text.
 */
inline std::ostream& operator<<(std::ostream& os, const EQN_STATE::TYPE& value)
{
    static std::map<EQN_STATE::TYPE, std::string> strings;
    
    if (strings.size() == 0)
    {
#define INSERT_ELEMENT(p) strings[p] = #p
        INSERT_ELEMENT(EQN_STATE::IDEAL_GAS);
        INSERT_ELEMENT(EQN_STATE::STIFFENED_GAS);
#undef INSERT_ELEMENT
    }
    
    return os << strings[value];
}

#endif /* EQUATIONS_OF_STATE_HPP */
