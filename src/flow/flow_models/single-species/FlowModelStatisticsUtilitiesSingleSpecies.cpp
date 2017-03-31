#include "flow/flow_models/single-species/FlowModelStatisticsUtilitiesSingleSpecies.hpp"

#include "SAMRAI/hier/FlattenedHierarchy.h"

#include <fstream>

/*
 * Output names of statistical quantities to output to a file.
 */
void
FlowModelStatisticsUtilitiesSingleSpecies::outputStatisticalQuantitiesNames(
    const std::string& filename_statistics)
{
    NULL_USE(filename_statistics);
}


/*
 * Output statisitcal quantities to a file.
 */
void
FlowModelStatisticsUtilitiesSingleSpecies::outputStatisticalQuantities(
    const std::string& filename_statistics,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
    NULL_USE(filename_statistics);
    NULL_USE(patch_hierarchy);
    NULL_USE(data_context);
}
