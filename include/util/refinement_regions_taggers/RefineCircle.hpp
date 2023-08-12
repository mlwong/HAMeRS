#ifndef REFINE_CIRCLE_HPP
#define REFINE_CIRCLE_HPP

#include "HAMeRS_config.hpp"

#include "util/Directions.hpp"

#include <vector>

struct RefineCircle
{
    const std::vector<Real> center_coord;
    const Real radius;
    const DIRECTION::TYPE direction;
    const int num_refine_levels;
    
    RefineCircle(
        const std::vector<Real>& center_coord,
        const Real radius,
        const DIRECTION::TYPE direction,
        const int num_refine_levels):
            center_coord(center_coord),
            radius(radius),
            direction(direction),
            num_refine_levels(num_refine_levels)
            {}
};

#endif /* REFINE_CIRCLE_HPP */