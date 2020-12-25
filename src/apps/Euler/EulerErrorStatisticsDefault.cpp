#include "apps/Euler/EulerErrorStatistics.hpp"

/*
 * Set the data on the patch interior to some initial values.
 */
void
EulerErrorStatistics::printErrorStatistics(
    std::ostream& os,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& variable_context) const
{
    NULL_USE(os);
    NULL_USE(patch_hierarchy);
    NULL_USE(variable_context);
}
