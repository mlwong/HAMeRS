#ifndef CONVECTIVE_FLUX_RECONSTRUCTORS_HPP
#define CONVECTIVE_FLUX_RECONSTRUCTORS_HPP

#include "flow/convective_flux_reconstructors/first_order/ConvectiveFluxReconstructorFirstOrderLLF.hpp"
#include "flow/convective_flux_reconstructors/first_order/ConvectiveFluxReconstructorFirstOrderHLLC.hpp"
#include "flow/convective_flux_reconstructors/WCNS56/ConvectiveFluxReconstructorWCNS5-JS-HLLC-HLL.hpp"
#include "flow/convective_flux_reconstructors/WCNS56/ConvectiveFluxReconstructorWCNS5-Z-HLLC-HLL.hpp"
#include "flow/convective_flux_reconstructors/WCNS56/ConvectiveFluxReconstructorWCNS6-CU-M2-HLLC-HLL.hpp"
#include "flow/convective_flux_reconstructors/WCNS56/ConvectiveFluxReconstructorWCNS6-LD-HLLC-HLL.hpp"
#include "flow/convective_flux_reconstructors/WCNS56/ConvectiveFluxReconstructorWCNS6-Test.hpp"
#include "flow/convective_flux_reconstructors/central/ConvectiveFluxReconstructorCentral.hpp"
#include "flow/convective_flux_reconstructors/DRP/ConvectiveFluxReconstructorDRP4.hpp"

#include <map>
#include <string>

namespace CONVECTIVE_FLUX_RECONSTRUCTOR
{
    enum TYPE { FIRST_ORDER_LLF,
                FIRST_ORDER_HLLC,
                WCNS5_JS_HLLC_HLL,
                WCNS5_Z_HLLC_HLL,
                WCNS6_CU_M2_HLLC_HLL,
                WCNS6_LD_HLLC_HLL,
                WCNS6_TEST,
                CENTRAL,
                DRP4 };
}

/*
 * Function to print out enum CONVECTIVE_FLUX_RECONSTRUCTOR::TYPE value as text.
 */
inline std::ostream& operator<<(std::ostream& os, const CONVECTIVE_FLUX_RECONSTRUCTOR::TYPE& value)
{
    static std::map<CONVECTIVE_FLUX_RECONSTRUCTOR::TYPE, std::string> strings;
    
    if (strings.size() == 0)
    {
#define INSERT_ELEMENT(p) strings[p] = #p
        INSERT_ELEMENT(CONVECTIVE_FLUX_RECONSTRUCTOR::FIRST_ORDER_LLF);
        INSERT_ELEMENT(CONVECTIVE_FLUX_RECONSTRUCTOR::FIRST_ORDER_HLLC);
        INSERT_ELEMENT(CONVECTIVE_FLUX_RECONSTRUCTOR::WCNS5_JS_HLLC_HLL);
        INSERT_ELEMENT(CONVECTIVE_FLUX_RECONSTRUCTOR::WCNS5_Z_HLLC_HLL);
        INSERT_ELEMENT(CONVECTIVE_FLUX_RECONSTRUCTOR::WCNS6_CU_M2_HLLC_HLL);
        INSERT_ELEMENT(CONVECTIVE_FLUX_RECONSTRUCTOR::WCNS6_LD_HLLC_HLL);
        INSERT_ELEMENT(CONVECTIVE_FLUX_RECONSTRUCTOR::WCNS6_TEST);
        INSERT_ELEMENT(CONVECTIVE_FLUX_RECONSTRUCTOR::CENTRAL);
        INSERT_ELEMENT(CONVECTIVE_FLUX_RECONSTRUCTOR::DRP4);
#undef INSERT_ELEMENT
    }
    
    return os << strings[value];
}

#endif /* CONVECTIVE_FLUX_RECONSTRUCTORS_HPP */
