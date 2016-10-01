#ifndef MIXING_CLOSURE_MODEL_HPP
#define MIXING_CLOSURE_MODEL_HPP

#include <map>
#include <string>

enum MIXING_CLOSURE_MODEL { ISOTHERMAL_AND_ISOBARIC,
                            ISOBARIC,
                            NO_MODEL };

/*
 * Function to print out enum MIXING_CLOSURE_MODEL value as text.
 */
inline std::ostream& operator<<(std::ostream& os, const MIXING_CLOSURE_MODEL& value)
{
    static std::map<MIXING_CLOSURE_MODEL, std::string> strings;
    
    if (strings.size() == 0)
    {
#define INSERT_ELEMENT(p) strings[p] = #p
        INSERT_ELEMENT(ISOTHERMAL_AND_ISOBARIC);
        INSERT_ELEMENT(ISOBARIC);
        INSERT_ELEMENT(NO_MODEL);
#undef INSERT_ELEMENT
    }
    
    return os << strings[value];
}

#endif /* MIXING_CLOSURE_MODEL_HPP */
