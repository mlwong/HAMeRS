#ifndef FLOW_MODEL_HELPER_HPP
#define FLOW_MODEL_HELPER_HPP

#include "HAMeRS_config.hpp"

#include "HAMeRS_memory.hpp"

#include "flow/flow_models/FlowModel.hpp"
#include "util/MPI_helpers/MPIHelper.hpp"

#include <string>

using namespace SAMRAI;

class FlowModelMPIHelper: public MPIHelper
{
    public:
        FlowModelMPIHelper(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<FlowModel>& flow_model):
                MPIHelper(
                    object_name,
                    dim,
                    grid_geometry,
                    patch_hierarchy),
                d_flow_model(flow_model)
        {}
        
    protected:
        /*
         * HAMERS_SHARED_PTR to FlowModel.
         */
        HAMERS_SHARED_PTR<FlowModel> d_flow_model;
        
};

#endif /* FLOW_MODEL_HELPER_HPP */
