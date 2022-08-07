#include "flow/flow_models/five-eqn_Allaire/FlowModelMonitoringStatisticsUtilitiesFiveEqnAllaire.hpp"

#include "flow/flow_models/MPI_helpers/FlowModelMPIHelperAverage.hpp"
#include "flow/flow_models/MPI_helpers/FlowModelMPIHelperMaxMin.hpp"

FlowModelMonitoringStatisticsUtilitiesFiveEqnAllaire::FlowModelMonitoringStatisticsUtilitiesFiveEqnAllaire(
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
        d_kinetic_energy_avg(double(0)),
        d_Mach_num_max(double(0))
{
    for (int si = 0; si < static_cast<int>(d_monitoring_statistics.size()); si++)
    {
        // Get the key of the current variable.
        std::string statistical_quantity_key = d_monitoring_statistics[si];
        
        if ((statistical_quantity_key != "KINETIC_ENERGY_AVG") &&
            (statistical_quantity_key != "MACH_NUM_MAX"))
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelMonitoringStatisticsUtilitiesFiveEqnAllaire::"
                << "FlowModelMonitoringStatisticsUtilitiesFiveEqnAllaire()\n"
                << "Unknown monitoring statistics with variable_key = '" << statistical_quantity_key
                << "' requested."
                << std::endl);
        }
    }
}


/*
 * Compute monitoring statistics.
 */
void
FlowModelMonitoringStatisticsUtilitiesFiveEqnAllaire::computeMonitoringStatistics(
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
        
        FlowModelMPIHelperMaxMin MPI_helper_max_min = FlowModelMPIHelperMaxMin(
            "MPI_helper_max_min",
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
            else if (statistical_quantity_key == "MACH_NUM_MAX")
            {
                std::vector<std::string> quantity_names;
                std::vector<int> component_indices;
                std::vector<bool> use_reciprocal;
                
                double u_sq_over_c_sq_max = double(0);
                double v_sq_over_c_sq_max = double(0);
                double w_sq_over_c_sq_max = double(0);
                
                if (d_dim == tbox::Dimension(1) || d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
                {
                    quantity_names.push_back("VELOCITY");
                    component_indices.push_back(0);
                    use_reciprocal.push_back(false);
                    quantity_names.push_back("VELOCITY");
                    component_indices.push_back(0);
                    use_reciprocal.push_back(false);
                    quantity_names.push_back("SOUND_SPEED");
                    component_indices.push_back(0);
                    use_reciprocal.push_back(true);
                    quantity_names.push_back("SOUND_SPEED");
                    component_indices.push_back(0);
                    use_reciprocal.push_back(true);
                    
                    u_sq_over_c_sq_max = MPI_helper_max_min.getMaxQuantity(
                        quantity_names,
                        component_indices,
                        use_reciprocal,
                        data_context);
                    
                    quantity_names.clear();
                    component_indices.clear();
                    use_reciprocal.clear();
                }
                
                if (d_dim == tbox::Dimension(2) || d_dim == tbox::Dimension(3))
                {
                    quantity_names.push_back("VELOCITY");
                    component_indices.push_back(1);
                    use_reciprocal.push_back(false);
                    quantity_names.push_back("VELOCITY");
                    component_indices.push_back(1);
                    use_reciprocal.push_back(false);
                    quantity_names.push_back("SOUND_SPEED");
                    component_indices.push_back(0);
                    use_reciprocal.push_back(true);
                    quantity_names.push_back("SOUND_SPEED");
                    component_indices.push_back(0);
                    use_reciprocal.push_back(true);
                    
                    v_sq_over_c_sq_max = MPI_helper_max_min.getMaxQuantity(
                        quantity_names,
                        component_indices,
                        use_reciprocal,
                        data_context);
                    
                    quantity_names.clear();
                    component_indices.clear();
                    use_reciprocal.clear();
                }
                
                if (d_dim == tbox::Dimension(3))
                {
                    quantity_names.push_back("VELOCITY");
                    component_indices.push_back(2);
                    use_reciprocal.push_back(false);
                    quantity_names.push_back("VELOCITY");
                    component_indices.push_back(2);
                    use_reciprocal.push_back(false);
                    quantity_names.push_back("SOUND_SPEED");
                    component_indices.push_back(0);
                    use_reciprocal.push_back(true);
                    quantity_names.push_back("SOUND_SPEED");
                    component_indices.push_back(0);
                    use_reciprocal.push_back(true);
                    
                    w_sq_over_c_sq_max = MPI_helper_max_min.getMaxQuantity(
                        quantity_names,
                        component_indices,
                        use_reciprocal,
                        data_context);
                    
                    quantity_names.clear();
                    component_indices.clear();
                    use_reciprocal.clear();
                }
                
                d_Mach_num_max = sqrt(fmax(fmax(u_sq_over_c_sq_max, v_sq_over_c_sq_max), w_sq_over_c_sq_max));
            }
        }
    }
}


/*
 * Output monitoring statistics to screen.
 */
void
FlowModelMonitoringStatisticsUtilitiesFiveEqnAllaire::outputMonitoringStatistics(
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
            os << "Avg kinetic energy: " << d_kinetic_energy_avg << std::endl;
        }
        else if (statistical_quantity_key == "MACH_NUM_MAX")
        {
            os << "Max Mach number: " << d_Mach_num_max << std::endl;
        }
    }
}
