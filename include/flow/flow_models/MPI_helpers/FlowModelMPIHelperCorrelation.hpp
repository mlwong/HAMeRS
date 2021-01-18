#ifndef FLOW_MODEL_HELPER_CORRELATION_HPP
#define FLOW_MODEL_HELPER_CORRELATION_HPP

#include "flow/flow_models/FlowModel.hpp"

#include "flow/flow_models/MPI_helpers/FlowModelMPIHelper.hpp"

#include <string>

class FlowModelMPIHelperCorrelation: public FlowModelMPIHelper
{
    public:
        FlowModelMPIHelperCorrelation(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<FlowModel>& flow_model):
                FlowModelMPIHelper(
                    object_name,
                    dim,
                    grid_geometry,
                    patch_hierarchy,
                    flow_model)
        {}
        
        /*
         * Compute correlation with only x-direction as inhomogeneous direction.
         */
        std::vector<double> getQuantityCorrelationWithInhomogeneousXDirection(
            const std::vector<std::string>& quantity_names,
            const std::vector<int>& component_indices,
            const std::vector<std::vector<double> >& averaged_quantities,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Compute correlation with only x-direction as inhomogeneous direction.
         */
        std::vector<double> getQuantityCorrelationWithInhomogeneousXDirection(
            const std::vector<std::string>& quantity_names,
            const std::vector<int>& component_indices,
            const std::vector<bool>& use_reciprocal,
            const std::vector<std::vector<double> >& averaged_quantities,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
    private:
        
};

#endif /* FLOW_MODEL_HELPER_CORRELATION_HPP */
