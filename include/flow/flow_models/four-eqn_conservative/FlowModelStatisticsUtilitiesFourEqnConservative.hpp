#ifndef FLOW_MODEL_STATISTICS_UTILITIES_FOUR_EQN_CONSERVATIVE_HPP
#define FLOW_MODEL_STATISTICS_UTILITIES_FOUR_EQN_CONSERVATIVE_HPP

#include "flow/flow_models/FlowModelStatisticsUtilities.hpp"

class FlowModelStatisticsUtilitiesFourEqnConservative: public FlowModelStatisticsUtilities
{
    public:
        FlowModelStatisticsUtilitiesFourEqnConservative(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_species,
            const boost::shared_ptr<tbox::Database>& flow_model_db):
                FlowModelStatisticsUtilities(
                    object_name,
                    dim,
                    grid_geometry,
                    num_species,
                    flow_model_db)
        {}
        
        /*
         * Output names of statistical quantities to output to a file.
         */
        void
        outputStatisticalQuantitiesNames(
            const std::string& filename_statistics);
        
        /*
         * Output statisitcal quantities to a file.
         */
        void
        outputStatisticalQuantities(
            const std::string& filename_statistics,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
};

#endif /* FLOW_MODEL_STATISTICS_UTILITIES_FOUR_EQN_CONSERVATIVE_HPP */
