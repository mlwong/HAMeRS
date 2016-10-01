#ifndef FLOW_MODEL_MANAGER_HPP
#define FLOW_MODEL_MANAGER_HPP

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/appu/VisDerivedDataStrategy.h"
#include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/tbox/Dimension.h"

#include "flow/flow_models/FlowModels.hpp"
#include "flow/initial_conditions/InitialConditions.hpp"

#include "boost/shared_ptr.hpp"
#include <string>
#include <vector>

using namespace SAMRAI;

class FlowModelManager
{
    public:
        FlowModelManager(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_species,
            const boost::shared_ptr<tbox::Database>& flow_model_db,
            const std::string& flow_model_str);
        
        /*
         * Get the label of flow model.
         */
        const FLOW_MODEL_LABEL&
        getFlowModelLabel() const
        {
            return d_flow_model_label;
        }
        
        /*
         * Get the flow model.
         */
        boost::shared_ptr<FlowModel>
        getFlowModel() const
        {
            return d_flow_model;
        }
        
        /*
         * Initialize d_initial_conditions.
         */
        void
        initializeInitialConditions(
            const std::string& project_name,
            boost::shared_ptr<InitialConditions>& initial_conditions);
        
        /*
         * Print all characteristics of flow model manager.
         */
        void
        printClassData(std::ostream& os) const;
        
    private:
        /*
         * Initialize boost::shared_ptr of the variables in d_initial_conditions.
         */
        void
        setVariablesForInitialConditions();
        
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
         * Label of flow model.
         */
        FLOW_MODEL_LABEL d_flow_model_label;
        
        /*
         * Flow model.
         */
        boost::shared_ptr<FlowModel> d_flow_model;
        
        /*
         * Number of species.
         */
        const int d_num_species;
        
        /*
         * boost::shared_ptr to InitialConditions.
         */
        boost::shared_ptr<InitialConditions> d_initial_conditions;
        
};

#endif /* FLOW_MODEL_MANAGER_HPP */
