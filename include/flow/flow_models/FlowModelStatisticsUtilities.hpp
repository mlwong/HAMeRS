#ifndef FLOW_MODEL_STATISTICS_UTILITIES_HPP
#define FLOW_MODEL_STATISTICS_UTILITIES_HPP

#include "HAMeRS_config.hpp"

#include "HAMeRS_memory.hpp"

#include "flow/flow_models/FlowModel.hpp"
#include "util/derivatives/DerivativeFirstOrder.hpp"

#include "SAMRAI/geom/CartesianGridGeometry.h"

#include <string>

class FlowModel;

class FlowModelStatisticsUtilities
{
    public:
        FlowModelStatisticsUtilities(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_species,
            const boost::shared_ptr<tbox::Database>& flow_model_db):
                d_object_name(object_name),
                d_dim(dim),
                d_grid_geometry(grid_geometry),
                d_num_species(num_species)
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
        void setFlowModel(const boost::weak_ptr<FlowModel>& flow_model)
        {
            d_flow_model = flow_model;
        }
        
        /*
         * Put the characteristics of the class into the restart database.
         */
        void
        putToRestart(
            const boost::shared_ptr<tbox::Database>& restart_db) const
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
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context)
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
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context)
        {
            NULL_USE(level);
            NULL_USE(patch_hierarchy);
            NULL_USE(data_context);
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
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context) = 0;
        
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
         * boost::shared_ptr to the grid geometry.
         */
        const boost::shared_ptr<geom::CartesianGridGeometry> d_grid_geometry;
        
        /*
         * Number of species.
         */
        const int d_num_species;
        
        /*
         * boost::weak_ptr to FlowModel.
         */
        boost::weak_ptr<FlowModel> d_flow_model;
        
        /*
         * Names of statistical quantities to output.
         */
        std::vector<std::string> d_statistical_quantities;
        
};

#endif /* FLOW_MODEL_STATISTICS_UTILITIES_HPP */
