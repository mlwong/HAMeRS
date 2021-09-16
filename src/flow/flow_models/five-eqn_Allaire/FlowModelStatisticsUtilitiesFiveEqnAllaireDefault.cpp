#include "flow/flow_models/five-eqn_Allaire/FlowModelStatisticsUtilitiesFiveEqnAllaire.hpp"

/*
 * Output names of statistical quantities to output to a file.
 */
void
FlowModelStatisticsUtilitiesFiveEqnAllaire::outputStatisticalQuantitiesNames(
    const std::string& stat_dump_filename)
{
    NULL_USE(stat_dump_filename);
}

/*
 * Output statisitcal quantities to a file.
 */
void
FlowModelStatisticsUtilitiesFiveEqnAllaire::outputStatisticalQuantities(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double output_time)
{
    NULL_USE(stat_dump_filename);
    NULL_USE(patch_hierarchy);
    NULL_USE(data_context);
    NULL_USE(output_time);
}
