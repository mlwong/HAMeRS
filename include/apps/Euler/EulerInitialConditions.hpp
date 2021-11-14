#ifndef EULER_INITIAL_CONDITIONS_HPP
#define EULER_INITIAL_CONDITIONS_HPP

#include "HAMeRS_config.hpp"

#include "HAMeRS_memory.hpp"

#include "flow/flow_models/FlowModels.hpp"

#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/tbox/Dimension.h"

#include <string>
#include <vector>

using namespace SAMRAI;

class EulerInitialConditions
{
    public:
        EulerInitialConditions(
            const std::string& object_name,
            const std::string& project_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const FLOW_MODEL::TYPE& flow_model_type,
            const HAMERS_SHARED_PTR<FlowModel>& flow_model,
            const HAMERS_SHARED_PTR<tbox::Database>& initial_conditions_db):
                d_object_name(object_name),
                d_project_name(project_name),
                d_dim(dim),
                d_grid_geometry(grid_geometry),
                d_flow_model_type(flow_model_type),
                d_flow_model(flow_model),
                d_initial_conditions_db(initial_conditions_db)
        {}
        
        /*
         * Set the data on the patch interior to some initial values, depending on the flow problems
         * and flow models.
         */
        void
        initializeDataOnPatch(
            hier::Patch& patch,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& conservative_variables,
            const double data_time,
            const bool initial_time);
        
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
         * Flow model type.
         */
        const FLOW_MODEL::TYPE d_flow_model_type;
        
        /*
         * Flow model.
         */
        const HAMERS_SHARED_PTR<FlowModel> d_flow_model;
        
        /*
         * Initial conditions database.
         */
        const HAMERS_SHARED_PTR<tbox::Database> d_initial_conditions_db;
        
};

#endif /* EULER_INITIAL_CONDITIONS_HPP */
