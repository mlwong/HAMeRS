#ifndef BASIC_BOUNDARY_CONDITIONS_HPP
#define BASIC_BOUNDARY_CONDITIONS_HPP

#include <iostream>
#include <map>

//@{
//! @name Definitions for Face, Edge, and Node boundary conditions (see source code for more information):
namespace BDRY_COND
{
    namespace BASIC
    {
        enum TYPE
        {
            FLOW       = 90,
            REFLECT    = 91,
            DIRICHLET  = 92,
            NEUMANN    = 93,
            SYMMETRY   = 94,
            XFLOW      = 900,
            YFLOW      = 901,
            ZFLOW      = 902,
            XREFLECT   = 910,
            YREFLECT   = 911,
            ZREFLECT   = 912,
            XDIRICHLET = 920,
            YDIRICHLET = 921,
            ZDIRICHLET = 922,
            XNEUMANN   = 930,
            YNEUMANN   = 931,
            ZNEUMANN   = 932,
            XSYMMETRY  = 940,
            YSYMMETRY  = 941,
            ZSYMMETRY  = 942
        };
    }
}
//@}


/*
 * Function to print out enum BDRY_COND::BASIC::TYPE value as text.
 */
inline std::ostream& operator<<(std::ostream& output, const BDRY_COND::BASIC::TYPE value)
{
    static std::map<BDRY_COND::BASIC::TYPE, std::string> strings;
    
    if (strings.size() == 0)
    {
#define INSERT_ELEMENT(p) strings[p] = #p
        INSERT_ELEMENT(BDRY_COND::BASIC::FLOW);
        INSERT_ELEMENT(BDRY_COND::BASIC::REFLECT);
        INSERT_ELEMENT(BDRY_COND::BASIC::DIRICHLET);
        INSERT_ELEMENT(BDRY_COND::BASIC::NEUMANN);
        INSERT_ELEMENT(BDRY_COND::BASIC::SYMMETRY);
        INSERT_ELEMENT(BDRY_COND::BASIC::XFLOW);
        INSERT_ELEMENT(BDRY_COND::BASIC::YFLOW);
        INSERT_ELEMENT(BDRY_COND::BASIC::ZFLOW);
        INSERT_ELEMENT(BDRY_COND::BASIC::XREFLECT);
        INSERT_ELEMENT(BDRY_COND::BASIC::YREFLECT);
        INSERT_ELEMENT(BDRY_COND::BASIC::ZREFLECT);
        INSERT_ELEMENT(BDRY_COND::BASIC::XDIRICHLET);
        INSERT_ELEMENT(BDRY_COND::BASIC::YDIRICHLET);
        INSERT_ELEMENT(BDRY_COND::BASIC::ZDIRICHLET);
        INSERT_ELEMENT(BDRY_COND::BASIC::XNEUMANN);
        INSERT_ELEMENT(BDRY_COND::BASIC::YNEUMANN);
        INSERT_ELEMENT(BDRY_COND::BASIC::ZNEUMANN);
        INSERT_ELEMENT(BDRY_COND::BASIC::XSYMMETRY);
        INSERT_ELEMENT(BDRY_COND::BASIC::YSYMMETRY);
        INSERT_ELEMENT(BDRY_COND::BASIC::ZSYMMETRY);
#undef INSERT_ELEMENT
    }
    
    return output << strings[value];
}

#endif /* BASIC_BOUNDARY_CONDITIONS_HPP */
