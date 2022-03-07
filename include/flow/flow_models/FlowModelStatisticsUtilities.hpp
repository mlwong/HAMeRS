#ifndef FLOW_MODEL_STATISTICS_UTILITIES_HPP
#define FLOW_MODEL_STATISTICS_UTILITIES_HPP

#include "HAMeRS_config.hpp"

#include "HAMeRS_memory.hpp"

#include "flow/flow_models/FlowModel.hpp"
#include "util/derivatives/DerivativeFirstOrder.hpp"
#include "util/ensemble_statistics/EnsembleStatistics.hpp"

#include "SAMRAI/geom/CartesianGridGeometry.h"

#include <string>

class FlowModel;

class FlowModelStatisticsUtilities
{
    public:
        FlowModelStatisticsUtilities(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_species,
            const HAMERS_SHARED_PTR<tbox::Database>& flow_model_db):
                d_object_name(object_name),
                d_dim(dim),
                d_grid_geometry(grid_geometry),
                d_num_species(num_species),
                d_is_ensemble_statistics_initialized(false)
        {
            /*
             * Get the names of statistical quantities to output.
             */
            if (flow_model_db->keyExists("statistical_quantities"))
            {
                d_statistical_quantities = flow_model_db->getStringVector("statistical_quantities");
            }
            else if (flow_model_db->keyExists("d_statistical_quantities"))
            {
                d_statistical_quantities = flow_model_db->getStringVector("d_statistical_quantities");
            }
        }
        
        virtual ~FlowModelStatisticsUtilities() {}
        
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
            if (!d_statistical_quantities.empty())
            {
                restart_db->putStringVector("d_statistical_quantities", d_statistical_quantities);
            }
        }
        
        /*
         * Register the variables required for computing statistics.
         */
        virtual void
        registerVariables(
            RungeKuttaLevelIntegrator* integrator,
            const hier::IntVector& num_ghosts)
        {
            NULL_USE(integrator);
            NULL_USE(num_ghosts);
        }
        
        /*
         * Compute the variables required for computing statistics.
         */
        virtual void
        computeVariables(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
        {
            NULL_USE(patch_hierarchy);
            NULL_USE(data_context);
        }
        
        /*
         * Filter the variables required for computing statistics.
         */
        virtual void
        filterVariables(
            const int level,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
        {
            NULL_USE(level);
            NULL_USE(patch_hierarchy);
            NULL_USE(data_context);
        }
        
        /*
         * Compute statisitcal quantities.
         */
        virtual void
        computeStatisticalQuantities(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double statistics_data_time)
        {
            NULL_USE(patch_hierarchy);
            NULL_USE(data_context);
            NULL_USE(statistics_data_time);
        }
        
        /*
         * Output names of statistical quantities to output to a file.
         */
        virtual void
        outputStatisticalQuantitiesNames(
            const std::string& stat_dump_filename) = 0;
        
        /*
         * Output statisitcal quantities to a file.
         */
        virtual void
        outputStatisticalQuantities(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double output_time) = 0;
        
        /*
         * Get object of storing ensemble statistics.
         */
        HAMERS_SHARED_PTR<EnsembleStatistics>
        getEnsembleStatistics()
        {
            if (d_is_ensemble_statistics_initialized == false)
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "The ensemble statistics object is not initialized yet!"
                    << std::endl);
            }
            
            return d_ensemble_statistics;
        }
        
        /*
         * Set object of storing ensemble statistics.
         */
        void
        setEnsembleStatistics(const HAMERS_SHARED_PTR<EnsembleStatistics> ensemble_statistics)
        {
            d_ensemble_statistics = ensemble_statistics;
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
         * Names of statistical quantities to output.
         */
        std::vector<std::string> d_statistical_quantities;
        
        /*
         * Store ensemble statistics
         */
        HAMERS_SHARED_PTR<EnsembleStatistics> d_ensemble_statistics;
        
        bool d_is_ensemble_statistics_initialized;
};

#endif /* FLOW_MODEL_STATISTICS_UTILITIES_HPP */
