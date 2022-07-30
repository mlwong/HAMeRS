#include "flow/flow_models/four-eqn_conservative/FlowModelMonitoringStatisticsUtilitiesFourEqnConservative.hpp"

/*
 * Output monitoring statistics to screen.
 */
void
FlowModelMonitoringStatisticsUtilitiesFourEqnConservative::computeMonitoringStatistics(
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
FlowModelMonitoringStatisticsUtilitiesFourEqnConservative::outputMonitoringStatistics(
    std::ostream& os,
    const double time)
{
    NULL_USE(os);
    NULL_USE(time);
}
