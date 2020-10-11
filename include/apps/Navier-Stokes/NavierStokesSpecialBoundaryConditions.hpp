#ifndef NAVIER_STOKES_SPECIAL_BOUNDARY_CONDITIONS_HPP
#define NAVIER_STOKES_SPECIAL_BOUNDARY_CONDITIONS_HPP

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

class NavierStokesSpecialBoundaryConditions
{
    public:
        NavierStokesSpecialBoundaryConditions(
            const std::string& object_name,
            const std::string& project_name,
            const tbox::Dimension& dim,
            const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
            const FLOW_MODEL::TYPE& flow_model_type,
            const boost::shared_ptr<FlowModel>& flow_model):
                d_object_name(object_name),
                d_project_name(project_name),
                d_dim(dim),
                d_grid_geometry(grid_geometry),
                d_flow_model_type(flow_model_type),
                d_flow_model(flow_model)
        {}
        
        /*
         * Set the data on the patch physical boundary to some values, depending on the flow problems
         * and flow models.
         */
        void
        setSpecialBoundaryConditions(
            hier::Patch& patch,
            const std::vector<boost::shared_ptr<pdat::CellData<double> > >& conservative_variables,
            const double fill_time,
            const hier::IntVector& ghost_width_to_fill);
        
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
         * boost::shared_ptr to the grid geometry.
         */
        const boost::shared_ptr<geom::CartesianGridGeometry> d_grid_geometry;
        
        /*
         * Flow model type.
         */
        const FLOW_MODEL::TYPE d_flow_model_type;
        
        /*
         * Flow model.
         */
        const boost::shared_ptr<FlowModel> d_flow_model;
        
};

#endif /* NAVIER_STOKES_SPECIAL_BOUNDARY_CONDITIONS_HPP */
