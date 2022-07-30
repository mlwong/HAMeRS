#include "flow/flow_models/single-species/FlowModelMonitoringStatisticsUtilitiesSingleSpecies.hpp"

/*
 * Compute monitoring statistics.
 */
void
FlowModelMonitoringStatisticsUtilitiesSingleSpecies::computeMonitoringStatistics(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double time)
{
    NULL_USE(patch_hierarchy);
    NULL_USE(data_context);
    NULL_USE(time);
}


/*
 * Output monitoring statistics to screen.
 */
void
FlowModelMonitoringStatisticsUtilitiesSingleSpecies::outputMonitoringStatistics(
    std::ostream& os,
    const double time)
{
    NULL_USE(os);
    NULL_USE(time);
}
