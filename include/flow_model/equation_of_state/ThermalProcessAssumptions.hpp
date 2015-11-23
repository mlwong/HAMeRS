#ifndef THERMAL_PROCESS_ASSUMPTIONS_HPP
#define THERMAL_PROCESS_ASSUMPTIONS_HPP

#include <map>
#include <string>

enum THERMAL_PROCESS_ASSUMPTION { ISOTHERMAL,
                                  ISOBARIC,
                                  NO_ASSUMPTION };

/*
 * Function to print out enum THERMAL_PROCESS_ASSUMPTION value as text.
 */
inline std::ostream& operator<<(std::ostream& os, const THERMAL_PROCESS_ASSUMPTION& value)
{
    static std::map<THERMAL_PROCESS_ASSUMPTION, std::string> strings;
    
    if (strings.size() == 0)
    {
#define INSERT_ELEMENT(p) strings[p] = #p
        INSERT_ELEMENT(ISOTHERMAL);
        INSERT_ELEMENT(ISOBARIC);
        INSERT_ELEMENT(NO_ASSUMPTION);
#undef INSERT_ELEMENT
    }
    
    return os << strings[value];
}

#endif /* THERMAL_PROCESS_ASSUMPTIONS_HPP */
