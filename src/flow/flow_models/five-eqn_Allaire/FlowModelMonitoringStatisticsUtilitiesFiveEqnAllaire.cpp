#include "flow/flow_models/five-eqn_Allaire/FlowModelMonitoringStatisticsUtilitiesFiveEqnAllaire.hpp"

/*
 * Compute monitoring statistics.
 */
void
FlowModelMonitoringStatisticsUtilitiesFiveEqnAllaire::computeMonitoringStatistics(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const int step_num,
    const double time)
{
    NULL_USE(patch_hierarchy);
    NULL_USE(data_context);
    NULL_USE(step_num);
    NULL_USE(time);
}


/*
 * Output monitoring statistics to screen.
 */
void
FlowModelMonitoringStatisticsUtilitiesFiveEqnAllaire::outputMonitoringStatistics(
    std::ostream& os,
    const double time)
{
    NULL_USE(os);
    NULL_USE(time);
}
