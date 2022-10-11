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
            const HAMERS_SHARED_PTR<tbox::Database>& flow_model_db);
        
        ~FlowModelMonitoringStatisticsUtilitiesSingleSpecies() {}
        
        /*
         * Compute monitoring statistics.
         */
        void
        computeMonitoringStatistics(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const int step_num,
            const double time);
        
        /*
         * Output names of monitoring statistical quantities to output to a file.
         */
        void
        outputMonitoringStatisticalQuantitiesNames(
            const std::string& monitoring_stat_dump_filename) const;
        
        /*
         * Output monitoring statistics to screen.
         */
        void
        outputMonitoringStatistics(
            std::ostream& os,
            const std::string& monitoring_stat_dump_filename,
            const double time);
        
        /*
         * Get monitoring statistical quantities.
         */
        double getMonitoringStatistics(
            std::string statistics_name) const;
        
        /*
         * Get map of monitoring statistical quantities.
         */
        std::unordered_map<std::string, double> getMonitoringStatisticsMap() const;
        
    private:
        /*
         * Monitoring statistical quantities.
         */
         
        double d_kinetic_energy_avg;
        double d_Mach_num_max;
        
};

#endif /* FLOW_MODEL_MONITORING_STATISTICS_UTILITIES_SINGLE_SPECIES_HPP */
