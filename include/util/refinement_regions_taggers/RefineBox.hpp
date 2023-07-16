#ifndef REFINE_BOX_HPP
#define REFINE_BOX_HPP

#include "HAMeRS_config.hpp"

#include <vector>

struct RefineBox
{
    const std::vector<Real> x_lo;
    const std::vector<Real> x_hi;
    const int num_refine_levels;
    
    RefineBox(
        const std::vector<Real>& x_lo,
        const std::vector<Real>& x_hi,
        const int num_refine_levels):
            x_lo(x_lo),
            x_hi(x_hi),
            num_refine_levels(num_refine_levels)
            {}
};

#endif /* REFINE_BOX_HPP */
