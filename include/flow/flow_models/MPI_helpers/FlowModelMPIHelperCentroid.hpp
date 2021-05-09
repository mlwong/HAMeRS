#ifndef FLOW_MODEL_MPI_HELPER_CENTROID_HPP
#define FLOW_MODEL_MPI_HELPER_CENTROID_HPP

#include "flow/flow_models/FlowModel.hpp"

#include "flow/flow_models/MPI_helpers/FlowModelMPIHelper.hpp"

#include <string>

class FlowModelMPIHelperCentroid: public FlowModelMPIHelper
{
    public:
        FlowModelMPIHelperCentroid(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<FlowModel>& flow_model,
            const bool use_diffusive_flux_utilities = false):
                FlowModelMPIHelper(
                    object_name,
                    dim,
                    grid_geometry,
                    patch_hierarchy,
                    flow_model,
                    use_diffusive_flux_utilities)
        {}
        
        /*
         * Compute centroid in x-direction.
         */
        double getCentroidInXDirection(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute centroid in y-direction.
         */
        double getCentroidInYDirection(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute centroid in z-direction.
         */
        double getCentroidInZDirection(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
    private:
        
};

#endif /* FLOW_MODEL_MPI_HELPER_CENTROID_HPP */
