#ifndef BASIC_BOUNDARY_CONDITIONS_HPP
#define BASIC_BOUNDARY_CONDITIONS_HPP

#include <iostream>
#include <map>

//@{
//! @name Definitions for Face, Edge, and Node boundary conditions (see source code for more information):
namespace BdryCond
{
    namespace Basic
    {
        enum Type
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
 * Function to print out enum BdryCond::Basic::Type value as text.
 */
inline std::ostream& operator<<(std::ostream& output, const BdryCond::Basic::Type value)
{
    static std::map<BdryCond::Basic::Type, std::string> strings;
    
    if (strings.size() == 0)
    {
#define INSERT_ELEMENT(p) strings[p] = #p
        INSERT_ELEMENT(BdryCond::Basic::FLOW);
        INSERT_ELEMENT(BdryCond::Basic::REFLECT);
        INSERT_ELEMENT(BdryCond::Basic::DIRICHLET);
        INSERT_ELEMENT(BdryCond::Basic::NEUMANN);
        INSERT_ELEMENT(BdryCond::Basic::SYMMETRY);
        INSERT_ELEMENT(BdryCond::Basic::XFLOW);
        INSERT_ELEMENT(BdryCond::Basic::YFLOW);
        INSERT_ELEMENT(BdryCond::Basic::ZFLOW);
        INSERT_ELEMENT(BdryCond::Basic::XREFLECT);
        INSERT_ELEMENT(BdryCond::Basic::YREFLECT);
        INSERT_ELEMENT(BdryCond::Basic::ZREFLECT);
        INSERT_ELEMENT(BdryCond::Basic::XDIRICHLET);
        INSERT_ELEMENT(BdryCond::Basic::YDIRICHLET);
        INSERT_ELEMENT(BdryCond::Basic::ZDIRICHLET);
        INSERT_ELEMENT(BdryCond::Basic::XNEUMANN);
        INSERT_ELEMENT(BdryCond::Basic::YNEUMANN);
        INSERT_ELEMENT(BdryCond::Basic::ZNEUMANN);
        INSERT_ELEMENT(BdryCond::Basic::XSYMMETRY);
        INSERT_ELEMENT(BdryCond::Basic::YSYMMETRY);
        INSERT_ELEMENT(BdryCond::Basic::ZSYMMETRY);
#undef INSERT_ELEMENT
    }
    
    return output << strings[value];
}

#endif /* BASIC_BOUNDARY_CONDITIONS_HPP */
