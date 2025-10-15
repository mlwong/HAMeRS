#ifndef FLOW_MODEL_MONITORING_STATISTICS_UTILITIES_FIVE_EQN_ALLAIRE_HPP
#define FLOW_MODEL_MONITORING_STATISTICS_UTILITIES_FIVE_EQN_ALLAIRE_HPP

#include "flow/flow_models/FlowModelMonitoringStatisticsUtilities.hpp"

class FlowModelMonitoringStatisticsUtilitiesFiveEqnAllaire: public FlowModelMonitoringStatisticsUtilities
{
    public:
        FlowModelMonitoringStatisticsUtilitiesFiveEqnAllaire(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_species,
            const HAMERS_SHARED_PTR<tbox::Database>& flow_model_db);
        
        ~FlowModelMonitoringStatisticsUtilitiesFiveEqnAllaire() {}
        
    private:
        /*
         * Compute monitoring statistics.
         */
        void computeMonitoringStatisticsDerived(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const int step_num,
            const double time);
        
        /*
         * Output monitoring statistics.
         */
        void outputMonitoringStatisticsDerived(
            std::ostream& os,
            std::ofstream& f_out) const;
        
        /*
         * Get monitoring statistical quantities.
         */
        Real getMonitoringStatisticsDerived(std::string statistics_name) const;
        
        /*
         * Get map of monitoring statistical quantities.
         */
        void getMonitoringStatisticsMapDerived(
            std::unordered_map<std::string, Real>& monitoring_statistics_map) const;
        
        /*
         * Monitoring statistical quantities.
         */
         
        Real d_kinetic_energy_avg;
        Real d_Mach_num_max;
        
};

#endif /* FLOW_MODEL_MONITORING_STATISTICS_UTILITIES_FIVE_EQN_ALLAIRE_HPP */
