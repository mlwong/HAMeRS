#include "apps/Navier-Stokes/NavierStokesErrorStatistics.hpp"

void
NavierStokesErrorStatistics::printErrorStatistics(
    std::ostream& os,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& variable_context,
    const double time) const
{
    NULL_USE(os);
    NULL_USE(patch_hierarchy);
    NULL_USE(variable_context);
    NULL_USE(time);
}
