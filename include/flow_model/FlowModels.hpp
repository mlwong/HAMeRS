#ifndef FLOW_MODELS_HPP
#define FLOW_MODELS_HPP

#include <map>
#include <string>

enum FLOW_MODEL { SINGLE_SPECIES,
                  FOUR_EQN_CONSERVATIVE,
                  FOUR_EQN_SHYUE,
                  FIVE_EQN_ALLAIRE };

/*
 * Function to print out enum FLOW_MODEL value as text.
 */
inline std::ostream& operator<<(std::ostream& os, const FLOW_MODEL& value)
{
    static std::map<FLOW_MODEL, std::string> strings;
    
    if (strings.size() == 0)
    {
#define INSERT_ELEMENT(p) strings[p] = #p
        INSERT_ELEMENT(SINGLE_SPECIES);
        INSERT_ELEMENT(FOUR_EQN_CONSERVATIVE);
        INSERT_ELEMENT(FOUR_EQN_SHYUE);
        INSERT_ELEMENT(FIVE_EQN_ALLAIRE);
#undef INSERT_ELEMENT
    }
    
    return os << strings[value];
}

#endif /* FLOW_MODELS_HPP */
