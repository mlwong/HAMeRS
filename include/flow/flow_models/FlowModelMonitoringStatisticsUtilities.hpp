#ifndef FLOW_MODEL_MONITORING_STATISTICS_UTILITIES_HPP
#define FLOW_MODEL_MONITORING_STATISTICS_UTILITIES_HPP

#include "HAMeRS_config.hpp"

#include "HAMeRS_memory.hpp"

#include "flow/flow_models/FlowModel.hpp"

#include "SAMRAI/geom/CartesianGridGeometry.h"

#include <string>

class FlowModel;

class FlowModelMonitoringStatisticsUtilities
{
    public:
        FlowModelMonitoringStatisticsUtilities(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_species,
            const HAMERS_SHARED_PTR<tbox::Database>& flow_model_db);
        
        virtual ~FlowModelMonitoringStatisticsUtilities() {}
        
        /*
         * Set the weak pointer to the flow model from the parent FlowModel class.
         */
        void setFlowModel(const HAMERS_WEAK_PTR<FlowModel>& flow_model)
        {
            d_flow_model = flow_model;
        }
        
        /*
         * Put the characteristics of the class into the restart database.
         */
        void
        putToRestart(
            const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const;
        
        /*
         * Base function to compute monitoring statistics.
         */
        void
        computeMonitoringStatistics(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const int step_num,
            const double time);
        
        /*
         * Compute monitoring statistics.
         */
        virtual void
        computeMonitoringStatisticsDerived(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const int step_num,
            const double time) = 0;
        
        /*
         * Whether it is the step to output the monitoring statistics.
         */
        bool isStepToOutputMonitoringStatistics(
            const int step_num) const
        {
            if (d_monitoring_time_step_interval > 0 && step_num%d_monitoring_time_step_interval == 0)
            {
                return true;
            }
            return false;
        }
        
        /*
         * Output names of monitoring statistical quantities to output to a file.
         */
        void
        outputMonitoringStatisticalQuantitiesNames(
            const std::string& monitoring_stat_dump_filename) const;
        
        /*
         * Base function to output monitoring statistics.
         */
        void
        outputMonitoringStatistics(
            std::ostream& os,
            const std::string& monitoring_stat_dump_filename,
            const int step_num,
            const double time) const;
        
        /*
         * Output monitoring statistics.
         */
        virtual void
        outputMonitoringStatisticsDerived(
            std::ostream& os,
            std::ofstream& f_out) const = 0;
        
        /*
         * Whether the object has monitoring statistics.
         */
        bool hasMonitoringStatistics() const;
        
        /*
         * Get names of monitoring statistical quantities to output.
         */
        const std::vector<std::string>& getMonitoringStatisticsNames() const
        {
            return d_monitoring_statistics_names;
        }
        
        /*
         * Get monitoring time step interval.
         */
        int getMonitoringTimeStepInterval() const
        {
            return d_monitoring_time_step_interval;
        }
        
        /*
         * Get monitoring statistical quantities.
         */
        virtual Real getMonitoringStatistics(
            std::string statistics_name) const = 0;
        
        /*
         * Get map of monitoring statistical quantities.
         */
        virtual std::unordered_map<std::string, Real> getMonitoringStatisticsMap() const = 0;
        
    protected:
        /*
         * The object name is used for error/warning reporting.
         */
        const std::string d_object_name;
        
        /*
         * Problem dimension.
         */
        const tbox::Dimension d_dim;
        
        /*
         * HAMERS_SHARED_PTR to the grid geometry.
         */
        const HAMERS_SHARED_PTR<geom::CartesianGridGeometry> d_grid_geometry;
        
        /*
         * Number of species.
         */
        const int d_num_species;
        
        /*
         * HAMERS_WEAK_PTR to FlowModel.
         */
        HAMERS_WEAK_PTR<FlowModel> d_flow_model;
        
        /*
         * Names of monitoring statistical quantities to output.
         */
        std::vector<std::string> d_monitoring_statistics_names;
        
        /*
         * The monitoring time step interval.
         */
        int d_monitoring_time_step_interval;
        
        /*
         * Whether to monitor the immersed boundary.
         */
        bool d_monitor_immersed_boundary;
};

#endif /* FLOW_MODEL_MONITORING_STATISTICS_UTILITIES_HPP */
