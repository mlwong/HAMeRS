#ifndef FLOW_MODEL_MPI_HELPER_MAX_MIN_HPP
#define FLOW_MODEL_MPI_HELPER_MAX_MIN_HPP

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
         * Compute maximum value over the entire domain.
         */
        Real getMaxQuantity(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute maximum reciprocal of value over the entire domain.
         */
        Real getMaxReciprocalOfQuantity(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute maximum value (on product of variables) over the entire domain.
         */
        Real getMaxQuantity(
            const std::vector<std::string>& quantity_names,
            const std::vector<int>& component_indices,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute maximum value (on product of variables) over the entire domain.
         */
        Real getMaxQuantity(
            const std::vector<std::string>& quantity_names,
            const std::vector<int>& component_indices,
            const std::vector<bool>& use_reciprocal,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute minimum value over the entire domain.
         */
        Real getMinQuantity(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute minimum reciprocal of value over the entire domain.
         */
        Real getMinReciprocalOfQuantity(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute minimum value (on product of variables) over the entire domain.
         */
        Real getMinQuantity(
            const std::vector<std::string>& quantity_names,
            const std::vector<int>& component_indices,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute minimum value (on product of variables) over the entire domain.
         */
        Real getMinQuantity(
            const std::vector<std::string>& quantity_names,
            const std::vector<int>& component_indices,
            const std::vector<bool>& use_reciprocal,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute maximum value with only x-direction as inhomogeneous direction.
         */
        std::vector<Real> getMaxQuantityWithInhomogeneousXDirection(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute maximum reciprocal of value with only x-direction as inhomogeneous direction.
         */
        std::vector<Real> getMaxReciprocalOfQuantityWithInhomogeneousXDirection(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute maximum value (on product of variables) with only x-direction as inhomogeneous direction.
         */
        std::vector<Real> getMaxQuantityWithInhomogeneousXDirection(
            const std::vector<std::string>& quantity_names,
            const std::vector<int>& component_indices,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute maximum value (on product of variables) with only x-direction as inhomogeneous direction.
         */
        std::vector<Real> getMaxQuantityWithInhomogeneousXDirection(
            const std::vector<std::string>& quantity_names,
            const std::vector<int>& component_indices,
            const std::vector<bool>& use_reciprocal,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute minimum value with only x-direction as inhomogeneous direction.
         */
        std::vector<Real> getMinQuantityWithInhomogeneousXDirection(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute minimum reciprocal of value with only x-direction as inhomogeneous direction.
         */
        std::vector<Real> getMinReciprocalOfQuantityWithInhomogeneousXDirection(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute minimum value (on product of variables) with only x-direction as inhomogeneous direction.
         */
        std::vector<Real> getMinQuantityWithInhomogeneousXDirection(
            const std::vector<std::string>& quantity_names,
            const std::vector<int>& component_indices,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute minimum value (on product of variables) with only x-direction as inhomogeneous direction.
         */
        std::vector<Real> getMinQuantityWithInhomogeneousXDirection(
            const std::vector<std::string>& quantity_names,
            const std::vector<int>& component_indices,
            const std::vector<bool>& use_reciprocal,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute maximum value with only y-direction as inhomogeneous direction.
         */
        std::vector<Real> getMaxQuantityWithInhomogeneousYDirection(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute maximum reciprocal of value with only y-direction as inhomogeneous direction.
         */
        std::vector<Real> getMaxReciprocalOfQuantityWithInhomogeneousYDirection(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute maximum value (on product of variables) with only y-direction as inhomogeneous direction.
         */
        std::vector<Real> getMaxQuantityWithInhomogeneousYDirection(
            const std::vector<std::string>& quantity_names,
            const std::vector<int>& component_indices,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute maximum value (on product of variables) with only y-direction as inhomogeneous direction.
         */
        std::vector<Real> getMaxQuantityWithInhomogeneousYDirection(
            const std::vector<std::string>& quantity_names,
            const std::vector<int>& component_indices,
            const std::vector<bool>& use_reciprocal,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute minimum value with only y-direction as inhomogeneous direction.
         */
        std::vector<Real> getMinQuantityWithInhomogeneousYDirection(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute minimum reciprocal of value with only y-direction as inhomogeneous direction.
         */
        std::vector<Real> getMinReciprocalOfQuantityWithInhomogeneousYDirection(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute minimum value (on product of variables) with only y-direction as inhomogeneous direction.
         */
        std::vector<Real> getMinQuantityWithInhomogeneousYDirection(
            const std::vector<std::string>& quantity_names,
            const std::vector<int>& component_indices,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute minimum value (on product of variables) with only y-direction as inhomogeneous direction.
         */
        std::vector<Real> getMinQuantityWithInhomogeneousYDirection(
            const std::vector<std::string>& quantity_names,
            const std::vector<int>& component_indices,
            const std::vector<bool>& use_reciprocal,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute maximum value with only z-direction as inhomogeneous direction.
         */
        std::vector<Real> getMaxQuantityWithInhomogeneousZDirection(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute maximum reciprocal of value with only z-direction as inhomogeneous direction.
         */
        std::vector<Real> getMaxReciprocalOfQuantityWithInhomogeneousZDirection(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute maximum value (on product of variables) with only z-direction as inhomogeneous direction.
         */
        std::vector<Real> getMaxQuantityWithInhomogeneousZDirection(
            const std::vector<std::string>& quantity_names,
            const std::vector<int>& component_indices,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute maximum value (on product of variables) with only z-direction as inhomogeneous direction.
         */
        std::vector<Real> getMaxQuantityWithInhomogeneousZDirection(
            const std::vector<std::string>& quantity_names,
            const std::vector<int>& component_indices,
            const std::vector<bool>& use_reciprocal,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute minimum value with only z-direction as inhomogeneous direction.
         */
        std::vector<Real> getMinQuantityWithInhomogeneousZDirection(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute minimum reciprocal of value with only z-direction as inhomogeneous direction.
         */
        std::vector<Real> getMinReciprocalOfQuantityWithInhomogeneousZDirection(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute minimum value (on product of variables) with only z-direction as inhomogeneous direction.
         */
        std::vector<Real> getMinQuantityWithInhomogeneousZDirection(
            const std::vector<std::string>& quantity_names,
            const std::vector<int>& component_indices,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute minimum value (on product of variables) with only z-direction as inhomogeneous direction.
         */
        std::vector<Real> getMinQuantityWithInhomogeneousZDirection(
            const std::vector<std::string>& quantity_names,
            const std::vector<int>& component_indices,
            const std::vector<bool>& use_reciprocal,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute maximum location within quantity bounds in x-direction.
         */
        Real getMaxLocationWithinQuantityBoundsInXDirection(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const Real bound_lo,
            const Real bound_hi);
        
        /*
         * Compute minimum location within quantity bounds in x-direction.
         */
        Real getMinLocationWithinQuantityBoundsInXDirection(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const Real bound_lo,
            const Real bound_hi);
        
        /*
         * Compute maximum location within quantity bounds in y-direction.
         */
        Real getMaxLocationWithinQuantityBoundsInYDirection(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const Real bound_lo,
            const Real bound_hi);
        
        /*
         * Compute minimum location within quantity bounds in y-direction.
         */
        Real getMinLocationWithinQuantityBoundsInYDirection(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const Real bound_lo,
            const Real bound_hi);
        
        /*
         * Compute maximum location within quantity bounds in z-direction.
         */
        Real getMaxLocationWithinQuantityBoundsInZDirection(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const Real bound_lo,
            const Real bound_hi);
        
        /*
         * Compute minimum location within quantity bounds in z-direction.
         */
        Real getMinLocationWithinQuantityBoundsInZDirection(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const Real bound_lo,
            const Real bound_hi);
        
        /*
         * Compute maximum value of absolute value of gradient with only x-direction as inhomogeneous direction.
         */
        std::vector<Real> getMaxAbsoluteGradientWithInhomogeneousXDirection(
            const std::string quantity_name,
            const int component_idx,
            const int derivative_direction,
            const int num_ghosts_derivative,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute maximum value of magnitude of gradient with only x-direction as inhomogeneous direction.
         */
        std::vector<Real> getMaxMagnitudeGradientWithInhomogeneousXDirection(
            const std::string quantity_name,
            const int component_idx,
            const int num_ghosts_derivative,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
    private:
        
};

#endif /* FLOW_MODEL_MPI_HELPER_MAX_MIN_HPP */
