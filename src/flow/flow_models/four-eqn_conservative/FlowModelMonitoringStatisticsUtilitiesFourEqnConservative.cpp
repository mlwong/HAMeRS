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
        d_kinetic_energy_avg(Real(0)),
        d_Mach_num_max(Real(0))
{
    for (int si = 0; si < static_cast<int>(d_monitoring_statistics_names.size()); si++)
    {
        // Get the key of the current variable.
        std::string statistical_quantity_key = d_monitoring_statistics_names[si];
        
        if ((statistical_quantity_key != "KINETIC_ENERGY_AVG") &&
            (statistical_quantity_key != "MACH_NUM_MAX"))
        {
            if (!isStatisticsNameValidBase(statistical_quantity_key))
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
}


/*
 * Compute monitoring statistics.
 */
void
FlowModelMonitoringStatisticsUtilitiesFourEqnConservative::computeMonitoringStatisticsDerived(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const int step_num,
    const double time)
{
    NULL_USE(step_num);
    NULL_USE(time);
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
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
    
    for (int si = 0; si < static_cast<int>(d_monitoring_statistics_names.size()); si++)
    {
        // Get the key of the current variable.
        std::string statistical_quantity_key = d_monitoring_statistics_names[si];
        
        if (statistical_quantity_key == "KINETIC_ENERGY_AVG")
        {
            std::vector<std::string> quantity_names;
            std::vector<int> component_indices;
            
            Real u_sq_avg = Real(0);
            Real v_sq_avg = Real(0);
            Real w_sq_avg = Real(0);
            
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
            
            d_kinetic_energy_avg = Real(1)/Real(2)*(u_sq_avg + v_sq_avg + w_sq_avg);
        }
        else if (statistical_quantity_key == "MACH_NUM_MAX")
        {
            std::vector<std::string> quantity_names;
            std::vector<int> component_indices;
            std::vector<bool> use_reciprocal;
            
            Real u_sq_over_c_sq_max = Real(0);
            Real v_sq_over_c_sq_max = Real(0);
            Real w_sq_over_c_sq_max = Real(0);
            
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
            
            d_Mach_num_max = std::sqrt(std::max(std::max(u_sq_over_c_sq_max, v_sq_over_c_sq_max), w_sq_over_c_sq_max));
        }
    }
}


/*
 * Output monitoring statistics.
 */
void
FlowModelMonitoringStatisticsUtilitiesFourEqnConservative::outputMonitoringStatisticsDerived(
    std::ostream& os,
    std::ofstream& f_out) const
{
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    for (int si = 0; si < static_cast<int>(d_monitoring_statistics_names.size()); si++)
    {
        // Get the key of the current variable.
        std::string statistical_quantity_key = d_monitoring_statistics_names[si];
        
        if (statistical_quantity_key == "KINETIC_ENERGY_AVG")
        {
            os << "Avg kinetic energy: " << d_kinetic_energy_avg << std::endl;
            if (mpi.getRank() == 0)
            {
                f_out << std::scientific << std::setprecision(16) << std::setw(25) << d_kinetic_energy_avg;
            }
        }
        else if (statistical_quantity_key == "MACH_NUM_MAX")
        {
            os << "Max Mach number: " << d_Mach_num_max << std::endl;
            if (mpi.getRank() == 0)
            {
                f_out << std::scientific << std::setprecision(16) << std::setw(25) << d_Mach_num_max;
            }
        }
    }
}


/*
 * Get monitoring statistical quantities.
 */
Real
FlowModelMonitoringStatisticsUtilitiesFourEqnConservative::getMonitoringStatisticsDerived(
    std::string statistics_name) const
{
    Real statistical_quantity = Real(0);
    
    if (statistics_name == "KINETIC_ENERGY_AVG")
    {
        statistical_quantity = d_kinetic_energy_avg;
    }
    else if (statistics_name == "MACH_NUM_MAX")
    {
        statistical_quantity = d_Mach_num_max;
    }
    
    return statistical_quantity;
}


/*
 * Get map of monitoring statistical quantities.
 */
void
FlowModelMonitoringStatisticsUtilitiesFourEqnConservative::getMonitoringStatisticsMapDerived(std::unordered_map<std::string, Real>& monitoring_statistics_map) const
{
    monitoring_statistics_map.insert(std::pair<std::string, Real>("KINETIC_ENERGY_AVG", d_kinetic_energy_avg));
    monitoring_statistics_map.insert(std::pair<std::string, Real>("MACH_NUM_MAX", d_Mach_num_max));
}
