#include "flow/flow_models/four-eqn_conservative/FlowModelMonitoringStatisticsUtilitiesFourEqnConservative.hpp"

#include "flow/flow_models/MPI_helpers/FlowModelMPIHelperAverage.hpp"
#include "flow/flow_models/MPI_helpers/FlowModelMPIHelperMaxMin.hpp"

FlowModelMonitoringStatisticsUtilitiesFourEqnConservative::FlowModelMonitoringStatisticsUtilitiesFourEqnConservative(
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
            flow_model_db),
        d_kinetic_energy_avg(double(0))
{
    for (int si = 0; si < static_cast<int>(d_monitoring_statistics.size()); si++)
    {
        // Get the key of the current variable.
        std::string statistical_quantity_key = d_monitoring_statistics[si];
        
        if (statistical_quantity_key != "KINETIC_ENERGY_AVG") // &&
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelMonitoringStatisticsUtilitiesFourEqnConservative::"
                << "FlowModelMonitoringStatisticsUtilitiesFourEqnConservative()\n"
                << "Unknown monitoring statistics with variable_key = '" << statistical_quantity_key
                << "' requested."
                << std::endl);
        }
    }
}


/*
 * Output monitoring statistics to screen.
 */
void
FlowModelMonitoringStatisticsUtilitiesFourEqnConservative::computeMonitoringStatistics(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const int step_num,
    const double time)
{
    NULL_USE(time);
    
    if (d_monitoring_time_step_interval > 0 && step_num%d_monitoring_time_step_interval == 0)
    {
        HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
        FlowModelMPIHelperAverage MPI_helper_average = FlowModelMPIHelperAverage(
            "MPI_helper_average",
            d_dim,
            d_grid_geometry,
            patch_hierarchy,
            flow_model_tmp);
    
        for (int si = 0; si < static_cast<int>(d_monitoring_statistics.size()); si++)
        {
            // Get the key of the current variable.
            std::string statistical_quantity_key = d_monitoring_statistics[si];
            
            if (statistical_quantity_key == "KINETIC_ENERGY_AVG")
            {
                std::vector<std::string> quantity_names;
                std::vector<int> component_indices;
                
                double u_sq_avg = double(0);
                double v_sq_avg = double(0);
                double w_sq_avg = double(0);
                
                if (d_dim == tbox::Dimension(1) || d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
                {
                    quantity_names.push_back("VELOCITY");
                    component_indices.push_back(0);
                    quantity_names.push_back("VELOCITY");
                    component_indices.push_back(0);
                    
                    u_sq_avg = MPI_helper_average.getAveragedQuantity(
                        quantity_names,
                        component_indices,
                        data_context);
                    
                    quantity_names.clear();
                    component_indices.clear();
                }
                
                if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
                {
                    quantity_names.push_back("VELOCITY");
                    component_indices.push_back(1);
                    quantity_names.push_back("VELOCITY");
                    component_indices.push_back(1);
                    
                    v_sq_avg = MPI_helper_average.getAveragedQuantity(
                        quantity_names,
                        component_indices,
                        data_context);
                    
                    quantity_names.clear();
                    component_indices.clear();
                }
                
                if (d_dim == tbox::Dimension(3))
                {
                    quantity_names.push_back("VELOCITY");
                    component_indices.push_back(2);
                    quantity_names.push_back("VELOCITY");
                    component_indices.push_back(2);
                    
                    w_sq_avg = MPI_helper_average.getAveragedQuantity(
                        quantity_names,
                        component_indices,
                        data_context);
                    
                    quantity_names.clear();
                    component_indices.clear();
                }
                
                d_kinetic_energy_avg = double(1)/double(2)*(u_sq_avg + v_sq_avg + w_sq_avg);
            }
        }
    }
}


/*
 * Output monitoring statistics to screen.
 */
void
FlowModelMonitoringStatisticsUtilitiesFourEqnConservative::outputMonitoringStatistics(
    std::ostream& os,
    const double time)
{
    NULL_USE(time);
    
    for (int si = 0; si < static_cast<int>(d_monitoring_statistics.size()); si++)
    {
        // Get the key of the current variable.
        std::string statistical_quantity_key = d_monitoring_statistics[si];
        
        if (statistical_quantity_key == "KINETIC_ENERGY_AVG")
        {
            os << "KINETIC_ENERGY_AVG: " << d_kinetic_energy_avg << std::endl;
        }
    }
}
