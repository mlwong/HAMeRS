#ifndef FLOW_MODEL_MANAGER_HPP
#define FLOW_MODEL_MANAGER_HPP

#include "HAMeRS_config.hpp"

#include "HAMeRS_memory.hpp"

#include "flow/flow_models/FlowModels.hpp"

#include "SAMRAI/appu/VisDerivedDataStrategy.h"
#include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/tbox/Dimension.h"

#include <string>
#include <vector>

using namespace SAMRAI;

class FlowModelManager
{
    public:
        FlowModelManager(
            const std::string& object_name,
            const std::string& project_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_species,
            const HAMERS_SHARED_PTR<tbox::Database>& flow_model_db,
            const std::string& flow_model_str);
        
        /*
         * Get the type of flow model.
         */
        const FLOW_MODEL::TYPE&
        getFlowModelType() const
        {
            return d_flow_model_type;
        }
        
        /*
         * Get the flow model.
         */
        HAMERS_SHARED_PTR<FlowModel>
        getFlowModel() const
        {
            return d_flow_model;
        }
        
        /*
         * Print all characteristics of flow model manager.
         */
        void
        printClassData(std::ostream& os) const;
        
    private:
        /*
         * The object name is used for error/warning reporting.
         */
        const std::string d_object_name;
        
        /*
         * Name of the project.
         */
        std::string d_project_name;
        
        /*
         * Problem dimension.
         */
        const tbox::Dimension d_dim;
        
        /*
         * HAMERS_SHARED_PTR to the grid geometry.
         */
        const HAMERS_SHARED_PTR<geom::CartesianGridGeometry> d_grid_geometry;
        
        /*
         * Type of flow model.
         */
        FLOW_MODEL::TYPE d_flow_model_type;
        
        /*
         * Flow model.
         */
        HAMERS_SHARED_PTR<FlowModel> d_flow_model;
        
        /*
         * Number of species.
         */
        const int d_num_species;
        
};

#endif /* FLOW_MODEL_MANAGER_HPP */
