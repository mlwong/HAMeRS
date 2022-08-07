#include "flow/flow_models/FlowModelMonitoringStatisticsUtilities.hpp"

FlowModelMonitoringStatisticsUtilities::FlowModelMonitoringStatisticsUtilities(
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
        d_monitoring_statistics_names = flow_model_db->getStringVector("monitoring_statistics");
    }
    else if (flow_model_db->keyExists("d_monitoring_statistics_names"))
    {
        d_monitoring_statistics_names = flow_model_db->getStringVector("d_monitoring_statistics_names");
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


/*
 * Put the characteristics of the class into the restart database.
 */
void
FlowModelMonitoringStatisticsUtilities::putToRestart(
    const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const
{
    if (!d_monitoring_statistics_names.empty())
    {
        restart_db->putStringVector("d_monitoring_statistics_names", d_monitoring_statistics_names);
    }
    
    restart_db->putInteger("d_monitoring_time_step_interval", d_monitoring_time_step_interval);
}
