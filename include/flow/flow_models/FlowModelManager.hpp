#ifndef FLOW_MODEL_MANAGER_HPP
#define FLOW_MODEL_MANAGER_HPP

#include "HAMeRS_config.hpp"

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

#include "boost/shared_ptr.hpp"
#include <string>
#include <vector>

using namespace SAMRAI;

class FlowModelManager
{
    public:
        FlowModelManager(
            const std::string& object_name,
            const boost::shared_ptr<tbox::Database>& flow_model_db,
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
         * Create the flow model.
         */
        boost::shared_ptr<FlowModel> createFlowModel(
            const tbox::Dimension& dim,
            const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_species);
        
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
         * Type of flow model.
         */
        FLOW_MODEL::TYPE d_flow_model_type;
        
        /*
         * boost::shared_ptr to the flow model database.
         */
        boost::shared_ptr<tbox::Database> d_flow_model_db;
        
};

#endif /* FLOW_MODEL_MANAGER_HPP */
