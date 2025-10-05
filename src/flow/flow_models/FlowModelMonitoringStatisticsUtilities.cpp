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
        d_monitoring_time_step_interval(-1),
        d_monitor_immersed_boundary(false)
{
    /*
     * Get the monitoring statistics database.
     */
    
    if (flow_model_db->keyExists("monitoring_statistics_names"))
    {
        d_monitoring_statistics_names = flow_model_db->getStringVector("monitoring_statistics_names");
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
    
    /*
     * Get the monitor immersed boundary flag.
     */
    if (flow_model_db->keyExists("monitor_immersed_boundary"))
    {
        d_monitor_immersed_boundary = flow_model_db->getBool("monitor_immersed_boundary");
    }
    else if (flow_model_db->keyExists("d_monitor_immersed_boundary"))
    {
        d_monitor_immersed_boundary = flow_model_db->getBool("d_monitor_immersed_boundary");
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
    restart_db->putBool("d_monitor_immersed_boundary", d_monitor_immersed_boundary);
}


/*
 * Base function to compute monitoring statistics.
 */
void
FlowModelMonitoringStatisticsUtilities::computeMonitoringStatistics(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const int step_num,
    const double time)
{
    computeMonitoringStatisticsDerived(
        patch_hierarchy,
        data_context,
        step_num,
        time);
}

/*
 * Whether the object has monitoring statistics.
 */
bool
FlowModelMonitoringStatisticsUtilities::hasMonitoringStatistics() const
{
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    bool hasMonitoringStatistics = false;
    
    if (d_monitoring_statistics_names.size() > 0)
    {
        hasMonitoringStatistics = true;
    }
    
    if (flow_model_tmp->useImmersedBoundary() && d_monitor_immersed_boundary)
    {
        hasMonitoringStatistics = true;
    }
    
    return hasMonitoringStatistics;
}


/*
 * Output names of monitoring statistical quantities to output to a file.
 */
void
FlowModelMonitoringStatisticsUtilities::outputMonitoringStatisticalQuantitiesNames(
    const std::string& monitoring_stat_dump_filename) const
{
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    if (mpi.getRank() == 0)
    {
        std::ofstream f_out;
        f_out.open(monitoring_stat_dump_filename.c_str(), std::ios::app);
        
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        for (int si = 0; si < static_cast<int>(d_monitoring_statistics_names.size()); si++)
        {
            // Get the key of the current variable.
            const std::string& statistical_quantity_key = d_monitoring_statistics_names[si];
            f_out << std::setw(25) << statistical_quantity_key;
        }
        
        if (flow_model_tmp->useImmersedBoundary() && d_monitor_immersed_boundary)
        {
            HAMERS_SHARED_PTR<FlowModelImmersedBoundaryMethod> flow_model_immersed_boundary_method =
                flow_model_tmp->getFlowModelImmersedBoundaryMethod();
            
            flow_model_immersed_boundary_method->outputMonitoringStatisticalQuantitiesNames(
                f_out);
        }
        
        f_out.close();
    }
}


/*
 * Base function to output monitoring statistics.
 */
void
FlowModelMonitoringStatisticsUtilities::outputMonitoringStatistics(
    std::ostream& os,
    const std::string& monitoring_stat_dump_filename,
    const int step_num,
    const double time) const
{
    NULL_USE(step_num);
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    std::ofstream f_out;
    
    if (mpi.getRank() == 0)
    {
        f_out.open(monitoring_stat_dump_filename.c_str(), std::ios::app);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output monitoring statistics!"
                << std::endl);
        }
        
        f_out << std::scientific << std::setprecision(16) << std::setw(25) << time;
    }
    
    outputMonitoringStatisticsDerived(os, f_out);
    
    if (mpi.getRank() == 0)
    {
        if (flow_model_tmp->useImmersedBoundary() && d_monitor_immersed_boundary)
        {
            HAMERS_SHARED_PTR<FlowModelImmersedBoundaryMethod> flow_model_immersed_boundary_method =
                flow_model_tmp->getFlowModelImmersedBoundaryMethod();
            
            flow_model_immersed_boundary_method->outputMonitoringStatistics(
                f_out);
        }
        
        f_out << std::endl;
        f_out.close();
    }
}