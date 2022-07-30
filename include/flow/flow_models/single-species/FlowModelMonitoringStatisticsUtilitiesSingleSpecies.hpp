#ifndef FLOW_MODEL_MONITORING_STATISTICS_UTILITIES_SINGLE_SPECIES_HPP
#define FLOW_MODEL_MONITORING_STATISTICS_UTILITIES_SINGLE_SPECIES_HPP

#include "flow/flow_models/FlowModelMonitoringStatisticsUtilities.hpp"

class FlowModelMonitoringStatisticsUtilitiesSingleSpecies: public FlowModelMonitoringStatisticsUtilities
{
    public:
        FlowModelMonitoringStatisticsUtilitiesSingleSpecies(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_species,
            const HAMERS_SHARED_PTR<tbox::Database>& flow_model_db):
                FlowModelMonitoringStatisticsUtilities(
                    object_name,
                    dim,
                    grid_geometry,
                    num_species,
                    flow_model_db)
        {}
        
        ~FlowModelMonitoringStatisticsUtilitiesSingleSpecies() {}
        
        /*
         * Compute monitoring statistics.
         */
        void
        computeMonitoringStatistics(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double time);
        
        /*
         * Output monitoring statistics to screen.
         */
        void
        outputMonitoringStatistics(
            std::ostream& os,
            const double time);
        
};

#endif /* FLOW_MODEL_MONITORING_STATISTICS_UTILITIES_SINGLE_SPECIES_HPP */
