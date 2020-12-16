#ifndef FLOW_MODEL_HELPER_HPP
#define FLOW_MODEL_HELPER_HPP

#include "flow/flow_models/FlowModel.hpp"

#include <string>

class FlowModelHelper
{
    public:
        FlowModelHelper(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const HAMERS_SHARED_PTR<FlowModel>& flow_model):
                d_object_name(object_name),
                d_dim(dim),
                d_grid_geometry(grid_geometry),
                d_flow_model(flow_model)
        {}
        
        /*
         * Get number of points in the x-direction of the refined domain.
         */
        int
        getRefinedDomainNumberOfPointsX(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const;
        
        /*
         * Get grid spacing in the x-direction of the refined domain.
         */
        double
        getRefinedDomainGridSpacingX(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const;
        
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
         * HAMERS_SHARED_PTR to FlowModel.
         */
        HAMERS_SHARED_PTR<FlowModel> d_flow_model;
        
};

#endif /* FLOW_MODEL_HELPER_HPP */
