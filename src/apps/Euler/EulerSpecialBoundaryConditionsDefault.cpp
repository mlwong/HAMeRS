#include "apps/Euler/EulerSpecialBoundaryConditions.hpp"

/*
 * Set the data on the patch physical boundary to some values, depending on the flow problems
 * and flow models.
 */
void
EulerSpecialBoundaryConditions::setSpecialBoundaryConditions(
    hier::Patch& patch,
    const std::vector<boost::shared_ptr<pdat::CellData<double> > >& conservative_variables,
    const double fill_time,
    const hier::IntVector& ghost_width_to_fill)
{
    NULL_USE(patch);
    NULL_USE(conservative_variables);
    NULL_USE(fill_time);
    NULL_USE(ghost_width_to_fill);
}
