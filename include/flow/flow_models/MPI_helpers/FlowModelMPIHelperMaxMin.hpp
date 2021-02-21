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
         * Compute maximum value with only y-direction as inhomogeneous direction.
         */
        std::vector<double> getMaxQuantityWithInhomogeneousYDirection(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Compute minimum value with only y-direction as inhomogeneous direction.
         */
        std::vector<double> getMinQuantityWithInhomogeneousYDirection(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Compute maximum value with only z-direction as inhomogeneous direction.
         */
        std::vector<double> getMaxQuantityWithInhomogeneousZDirection(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Compute minimum value with only z-direction as inhomogeneous direction.
         */
        std::vector<double> getMinQuantityWithInhomogeneousZDirection(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Compute maximum location within quantity bounds in x-direction.
         */
        double getMaxLocationWithinQuantityBoundsInXDirection(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double bound_lo,
            const double bound_hi) const;
        
        /*
         * Compute minimum location within quantity bounds in x-direction.
         */
        double getMinLocationWithinQuantityBoundsInXDirection(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double bound_lo,
            const double bound_hi) const;
        
        /*
         * Compute maximum location within quantity bounds in y-direction.
         */
        double getMaxLocationWithinQuantityBoundsInYDirection(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double bound_lo,
            const double bound_hi) const;
        
        /*
         * Compute minimum location within quantity bounds in y-direction.
         */
        double getMinLocationWithinQuantityBoundsInYDirection(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double bound_lo,
            const double bound_hi) const;
        
        /*
         * Compute maximum location within quantity bounds in z-direction.
         */
        double getMaxLocationWithinQuantityBoundsInZDirection(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double bound_lo,
            const double bound_hi) const;
        
        /*
         * Compute minimum location within quantity bounds in z-direction.
         */
        double getMinLocationWithinQuantityBoundsInZDirection(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double bound_lo,
            const double bound_hi) const;
        
        /*
         * Compute maximum value of absolute value of gradient with only x-direction as inhomogeneous direction.
         */
        std::vector<double> getMaxAbsoluteGradientWithInhomogeneousXDirection(
            const std::string quantity_name,
            const int component_idx,
            const int derivative_direction,
            const int num_ghosts_derivative,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
        /*
         * Compute maximum value of magnitude of gradient with only x-direction as inhomogeneous direction.
         */
        std::vector<double> getMaxMagnitudeGradientWithInhomogeneousXDirection(
            const std::string quantity_name,
            const int component_idx,
            const int num_ghosts_derivative,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context) const;
        
    private:
        
};

#endif /* FLOW_MODEL_HELPER_MAX_MIN_HPP */
