#ifndef MIXING_CLOSURE_MODEL_HPP
#define MIXING_CLOSURE_MODEL_HPP

#include <map>
#include <string>

namespace MIXING_CLOSURE_MODEL
{
    enum TYPE { ISOTHERMAL_AND_ISOBARIC,
                ISOBARIC,
                NO_MODEL };
}

/*
 * Function to print out enum MIXING_CLOSURE_MODEL::TYPE value as text.
 */
inline std::ostream& operator<<(std::ostream& os, const MIXING_CLOSURE_MODEL::TYPE& value)
{
    static std::map<MIXING_CLOSURE_MODEL::TYPE, std::string> strings;
    
    if (strings.size() == 0)
    {
#define INSERT_ELEMENT(p) strings[p] = #p
        INSERT_ELEMENT(MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC);
        INSERT_ELEMENT(MIXING_CLOSURE_MODEL::ISOBARIC);
        INSERT_ELEMENT(MIXING_CLOSURE_MODEL::NO_MODEL);
#undef INSERT_ELEMENT
    }
    
    return os << strings[value];
}

#endif /* MIXING_CLOSURE_MODEL_HPP */
