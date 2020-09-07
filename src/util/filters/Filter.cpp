#include "util/filters/Filter.hpp"

Filter::Filter(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const DIRECTION::TYPE& direction):
        d_object_name(object_name),
        d_dim(dim),
        d_direction(direction),
        d_num_filter_ghosts(hier::IntVector::getZero(d_dim))
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    if (d_direction == DIRECTION::Y_DIRECTION)
    {
        if (d_dim == tbox::Dimension(1))
        {
            TBOX_ERROR(d_object_name
                << ": Filter::Filter()\n"
                << "Filter in y-direction cannot be obtained for problem with diemension less than two."
                << std::endl);
        }
    }
    else if (d_direction == DIRECTION::Z_DIRECTION)
    {
        if (d_dim == tbox::Dimension(1) || d_dim == tbox::Dimension(2))
        {
            TBOX_ERROR(d_object_name
                << ": Filter::Filter()\n"
                << "Filter in z-direction cannot be obtained for problem with diemension less than three"
                << std::endl);
        }
    }
#endif
}
