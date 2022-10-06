#ifndef DIFFUSIVE_FLUX_RECONSTRUCTORS_HPP
#define DIFFUSIVE_FLUX_RECONSTRUCTORS_HPP

#include "flow/diffusive_flux_reconstructors/midpoint/DiffusiveFluxReconstructorMidpointSecondOrder.hpp"
#include "flow/diffusive_flux_reconstructors/midpoint/DiffusiveFluxReconstructorMidpointFourthOrder.hpp"
#include "flow/diffusive_flux_reconstructors/midpoint/DiffusiveFluxReconstructorMidpointSixthOrder.hpp"
#include "flow/diffusive_flux_reconstructors/node/DiffusiveFluxReconstructorNodeSecondOrder.hpp"
#include "flow/diffusive_flux_reconstructors/node/DiffusiveFluxReconstructorNodeFourthOrder.hpp"
#include "flow/diffusive_flux_reconstructors/node/DiffusiveFluxReconstructorNodeSixthOrder.hpp"

#include <map>
#include <string>

namespace DIFFUSIVE_FLUX_RECONSTRUCTOR
{
    enum TYPE { MIDPOINT_SECOND_ORDER,
                MIDPOINT_FOURTH_ORDER,
                MIDPOINT_SIXTH_ORDER,
                SECOND_ORDER,
                FOURTH_ORDER,
                SIXTH_ORDER };
}

/*
 * Function to print out enum DIFFUSIVE_FLUX_RECONSTRUCTOR::TYPE value as text.
 */
inline std::ostream& operator<<(std::ostream& os, const DIFFUSIVE_FLUX_RECONSTRUCTOR::TYPE& value)
{
    static std::map<DIFFUSIVE_FLUX_RECONSTRUCTOR::TYPE, std::string> strings;
    
    if (strings.size() == 0)
    {
#define INSERT_ELEMENT(p) strings[p] = #p
        INSERT_ELEMENT(DIFFUSIVE_FLUX_RECONSTRUCTOR::MIDPOINT_SECOND_ORDER);
        INSERT_ELEMENT(DIFFUSIVE_FLUX_RECONSTRUCTOR::MIDPOINT_FOURTH_ORDER);
        INSERT_ELEMENT(DIFFUSIVE_FLUX_RECONSTRUCTOR::MIDPOINT_SIXTH_ORDER);
        INSERT_ELEMENT(DIFFUSIVE_FLUX_RECONSTRUCTOR::SECOND_ORDER);
        INSERT_ELEMENT(DIFFUSIVE_FLUX_RECONSTRUCTOR::FOURTH_ORDER);
        INSERT_ELEMENT(DIFFUSIVE_FLUX_RECONSTRUCTOR::SIXTH_ORDER);
#undef INSERT_ELEMENT
    }
    
    return os << strings[value];
}

#endif /* DIFFUSIVE_FLUX_RECONSTRUCTORS_HPP */
