#ifndef EQUATIONS_OF_STATE_HPP
#define EQUATIONS_OF_STATE_HPP

#include "util/equations_of_state/ideal_gas/EquationOfStateIdealGas.hpp"
#include "util/equations_of_state/ideal_gas/EquationOfStateMixingRulesIdealGas.hpp"

#include <map>
#include <string>

enum EQUATION_OF_STATE_LABEL { IDEAL_GAS };

/*
 * Function to print out enum EQUATION_OF_STATE_LABEL value as text.
 */
inline std::ostream& operator<<(std::ostream& os, const EQUATION_OF_STATE_LABEL& value)
{
    static std::map<EQUATION_OF_STATE_LABEL, std::string> strings;
    
    if (strings.size() == 0)
    {
#define INSERT_ELEMENT(p) strings[p] = #p
        INSERT_ELEMENT(IDEAL_GAS);
#undef INSERT_ELEMENT
    }
    
    return os << strings[value];
}

#endif /* EQUATIONS_OF_STATE_HPP */
