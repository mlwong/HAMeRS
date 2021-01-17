#ifndef FLOW_MODEL_HELPER_MAX_MIN_HPP
#define FLOW_MODEL_HELPER_MAX_MIN_HPP

#include "flow/flow_models/FlowModel.hpp"

#include "flow/flow_models/MPI_helpers/FlowModelMPIHelper.hpp"

#include <string>

class FlowModelMPIHelperMaxMin: public FlowModelMPIHelper
{
    public:
        FlowModelMPIHelperMaxMin(
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
         * Compute maximum value with only x-direction as inhomogeneous direction.
         */
        std::vector<double> getMaxQuantityWithInhomogeneousXDirection(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Compute minimum value with only x-direction as inhomogeneous direction.
         */
        std::vector<double> getMinQuantityWithInhomogeneousXDirection(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Compute maximum location within quantity bounds in x-direction.
         */
        double getMaxLocationWithinBoundsInXDirection(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double bound_lo,
            const double bound_hi) const;
        
        /*
         * Compute minimum location within quantity bounds in x-direction.
         */
        double getMinLocationWithinBoundsInXDirection(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double bound_lo,
            const double bound_hi) const;
        
        /*
         * Compute maximum location within quantity bounds in y-direction.
         */
        double getMaxLocationWithinBoundsInYDirection(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double bound_lo,
            const double bound_hi) const;
        
        /*
         * Compute minimum location within quantity bounds in y-direction.
         */
        double getMinLocationWithinBoundsInYDirection(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double bound_lo,
            const double bound_hi) const;
        
        /*
         * Compute maximum location within quantity bounds in z-direction.
         */
        double getMaxLocationWithinBoundsInZDirection(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double bound_lo,
            const double bound_hi) const;
        
        /*
         * Compute minimum location within quantity bounds in z-direction.
         */
        double getMinLocationWithinBoundsInZDirection(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double bound_lo,
            const double bound_hi) const;
        
    private:
        
};

#endif /* FLOW_MODEL_HELPER_MAX_MIN_HPP */
