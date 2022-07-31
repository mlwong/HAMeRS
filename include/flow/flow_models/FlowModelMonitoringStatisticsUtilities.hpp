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
            const HAMERS_SHARED_PTR<tbox::Database>& flow_model_db):
                d_object_name(object_name),
                d_dim(dim),
                d_grid_geometry(grid_geometry),
                d_num_species(num_species),
                d_monitoring_time_step_interval(-1)
        {
            /*
             * Get the monitoring statistics database.
             */
            
            if (flow_model_db->keyExists("monitoring_statistics"))
            {
                d_monitoring_statistics = flow_model_db->getStringVector("monitoring_statistics");
            }
            else if (flow_model_db->keyExists("d_monitoring_statistics"))
            {
                d_monitoring_statistics = flow_model_db->getStringVector("d_monitoring_statistics");
            }
            
            /*
             * Get the monitoring time step interval.
             */
            
            if (flow_model_db->keyExists("monitoring_time_step_interval"))
            {
                d_monitoring_time_step_interval = flow_model_db->getInteger("monitoring_time_step_interval");
            }
            else if (flow_model_db->keyExists("d_monitoring_time_step_interval"))
            {
                d_monitoring_time_step_interval = flow_model_db->getInteger("d_monitoring_time_step_interval");
            }
        }
        
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
            const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const
        {
            if (!d_monitoring_statistics.empty())
            {
                restart_db->putStringVector("d_monitoring_statistics", d_monitoring_statistics);
            }
            
            restart_db->putInteger("d_monitoring_time_step_interval", d_monitoring_time_step_interval);
        }
        
        /*
         * Compute monitoring statistics.
         */
        virtual void
        computeMonitoringStatistics(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const int step_num,
            const double time) = 0;
        
        /*
         * Output monitoring statistics to screen.
         */
        virtual void
        outputMonitoringStatistics(
            std::ostream& os,
            const double time) = 0;
        
        /*
         * Get names of monitoring statistical quantities to output.
         */
        const std::vector<std::string>& getMonitoringStatistics() const
        {
            return d_monitoring_statistics;
        }
        
        /*
         * Get monitoring time step interval.
         */
        int getMonitoringTimeStepInterval() const
        {
            return d_monitoring_time_step_interval;
        }
        
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
        std::vector<std::string> d_monitoring_statistics;
        
        /*
         * The monitoring time step interval.
         */
        int d_monitoring_time_step_interval;
};

#endif /* FLOW_MODEL_MONITORING_STATISTICS_UTILITIES_HPP */
