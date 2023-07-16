#ifndef REFINE_TRIANGLE_HPP
#define REFINE_TRIANGLE_HPP

#include "HAMeRS_config.hpp"

#include "util/Directions.hpp"

#include <vector>

struct RefineTriangle
{
    const std::vector<Real> point_coord_0;
    const std::vector<Real> point_coord_1;
    const std::vector<Real> point_coord_2;
    const DIRECTION::TYPE direction;
    const int num_refine_levels;
    
    RefineTriangle(
        const std::vector<Real>& point_coord_0,
        const std::vector<Real>& point_coord_1,
        const std::vector<Real>& point_coord_2,
        const DIRECTION::TYPE direction,
        const int num_refine_levels):
            point_coord_0(point_coord_0),
            point_coord_1(point_coord_1),
            point_coord_2(point_coord_2),
            direction(direction),
            num_refine_levels(num_refine_levels)
            {}
};

#endif /* REFINE_TRIANGLE_HPP */