#ifndef CONVECTIVE_FLUX_RECONSTRUCTORS_HPP
#define CONVECTIVE_FLUX_RECONSTRUCTORS_HPP

#include "flow/convective_flux_reconstructors/ConvectiveFluxReconstructorFirstOrderLLF.hpp"
#include "flow/convective_flux_reconstructors/ConvectiveFluxReconstructorFirstOrderHLLC.hpp"
#include "flow/convective_flux_reconstructors/ConvectiveFluxReconstructorWCNS5-JS-HLLC-HLL.hpp"
#include "flow/convective_flux_reconstructors/ConvectiveFluxReconstructorWCNS5-Z-HLLC-HLL.hpp"
#include "flow/convective_flux_reconstructors/ConvectiveFluxReconstructorWCNS6-CU-M2-HLLC-HLL.hpp"
#include "flow/convective_flux_reconstructors/ConvectiveFluxReconstructorWCNS6-HW-HLLC-HLL.hpp"
#include "flow/convective_flux_reconstructors/ConvectiveFluxReconstructorWCNS6-HW-LD-HLLC-HLL.hpp"
#include "flow/convective_flux_reconstructors/ConvectiveFluxReconstructorWCNS6-Test.hpp"

#include <map>
#include <string>

enum CONVECTIVE_FLUX_RECONSTRUCTOR_LABEL { FIRST_ORDER_LLF,
                                           FIRST_ORDER_HLLC,
                                           WCNS5_JS_HLLC_HLL,
                                           WCNS5_Z_HLLC_HLL,
                                           WCNS6_CU_M2_HLLC_HLL,
                                           WCNS6_HW_HLLC_HLL,
                                           WCNS6_HW_LD_HLLC_HLL,
                                           WCNS6_TEST};

/*
 * Function to print out enum CONVECTIVE_FLUX_RECONSTRUCTOR_LABEL value as text.
 */
inline std::ostream& operator<<(std::ostream& os, const CONVECTIVE_FLUX_RECONSTRUCTOR_LABEL& value)
{
    static std::map<CONVECTIVE_FLUX_RECONSTRUCTOR_LABEL, std::string> strings;
    
    if (strings.size() == 0)
    {
#define INSERT_ELEMENT(p) strings[p] = #p
        INSERT_ELEMENT(FIRST_ORDER_LLF);
        INSERT_ELEMENT(FIRST_ORDER_HLLC);
        INSERT_ELEMENT(WCNS5_JS_HLLC_HLL);
        INSERT_ELEMENT(WCNS5_Z_HLLC_HLL);
        INSERT_ELEMENT(WCNS6_CU_M2_HLLC_HLL);
        INSERT_ELEMENT(WCNS6_HW_HLLC_HLL);
        INSERT_ELEMENT(WCNS6_HW_LD_HLLC_HLL);
        INSERT_ELEMENT(WCNS6_TEST);
#undef INSERT_ELEMENT
    }
    
    return os << strings[value];
}

#endif /* CONVECTIVE_FLUX_RECONSTRUCTORS_HPP */
