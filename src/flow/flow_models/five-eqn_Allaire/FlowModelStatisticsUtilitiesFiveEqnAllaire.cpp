#include "flow/flow_models/five-eqn_Allaire/FlowModelStatisticsUtilitiesFiveEqnAllaire.hpp"

#include "SAMRAI/hier/FlattenedHierarchy.h"

#include <fstream>

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
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
    NULL_USE(stat_dump_filename);
    NULL_USE(patch_hierarchy);
    NULL_USE(data_context);
}
