#ifndef FLOW_MODEL_MPI_HELPER_AVERAGE_HPP
#define FLOW_MODEL_MPI_HELPER_AVERAGE_HPP

#include "flow/flow_models/FlowModel.hpp"

#include "flow/flow_models/MPI_helpers/FlowModelMPIHelper.hpp"

#include <string>

class FlowModelMPIHelperAverage: public FlowModelMPIHelper
{
    public:
        FlowModelMPIHelperAverage(
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
         * Compute averaged value over the entire domain.
         */
        double getAveragedQuantity(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute averaged reciprocal of value over the entire domain.
         */
        double getAveragedReciprocalOfQuantity(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute averaged value (on product of variables) over the entire domain.
         */
        double getAveragedQuantity(
            const std::vector<std::string>& quantity_names,
            const std::vector<int>& component_indices,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute averaged value (on product of variables) over the entire domain.
         */
        double getAveragedQuantity(
            const std::vector<std::string>& quantity_names,
            const std::vector<int>& component_indices,
            const std::vector<bool>& use_reciprocal,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute averaged value (on product of variable derivatives) over the entire domain.
         */
        double getAveragedQuantity(
            const std::vector<std::string>& quantity_names,
            const std::vector<int>& component_indices,
            const std::vector<bool>& use_derivative,
            const std::vector<int>& derivative_directions,
            const int num_ghosts_derivative,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute averaged value (on product of variable derivatives) over the entire domain.
         */
        double getAveragedQuantity(
            const std::vector<std::string>& quantity_names,
            const std::vector<int>& component_indices,
            const std::vector<bool>& use_derivative,
            const std::vector<int>& derivative_directions,
            const std::vector<bool>& use_reciprocal,
            const int num_ghosts_derivative,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute averaged derivative of value (on product of variables) over the entire domain.
         */
        double getAveragedDerivativeOfQuantity(
            const std::vector<std::string>& quantity_names,
            const std::vector<int>& component_indices,
            const std::vector<bool>& use_reciprocal,
            const int derivative_direction,
            const int num_ghosts_derivative,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute averaged value with only x-direction as inhomogeneous direction on the coarsest level.
         */
        std::vector<double> getAveragedQuantityWithInhomogeneousXDirectionOnCoarsestLevel(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute averaged reciprocal of value with only x-direction as inhomogeneous direction on the coarsest level.
         */
        std::vector<double> getAveragedReciprocalOfQuantityWithInhomogeneousXDirectionOnCoarsestLevel(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute averaged value (on product of variables) with only x-direction as inhomogeneous direction on
         * the coarsest level.
         */
        std::vector<double> getAveragedQuantityWithInhomogeneousXDirectionOnCoarsestLevel(
            const std::vector<std::string>& quantity_names,
            const std::vector<int>& component_indices,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute averaged value (on product of variables) with only x-direction as inhomogeneous direction on
         * the coarsest level.
         */
        std::vector<double> getAveragedQuantityWithInhomogeneousXDirectionOnCoarsestLevel(
            const std::vector<std::string>& quantity_names,
            const std::vector<int>& component_indices,
            const std::vector<bool>& use_reciprocal,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute averaged value (on product of variable derivatives) with only x direction as inhomogeneous direction
         * on the coarsest level.
         */
        std::vector<double> getAveragedQuantityWithInhomogeneousXDirectionOnCoarsestLevel(
            const std::vector<std::string>& quantity_names,
            const std::vector<int>& component_indices,
            const std::vector<bool>& use_derivative,
            const std::vector<int>& derivative_directions,
            const int num_ghosts_derivative,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute averaged value (on product of variable derivatives) with only x direction as inhomogeneous direction
         * on the coarsest level.
         */
        std::vector<double> getAveragedQuantityWithInhomogeneousXDirectionOnCoarsestLevel(
            const std::vector<std::string>& quantity_names,
            const std::vector<int>& component_indices,
            const std::vector<bool>& use_derivative,
            const std::vector<int>& derivative_directions,
            const std::vector<bool>& use_reciprocal,
            const int num_ghosts_derivative,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute averaged derivative of value (on product of variables) with only x direction as inhomogeneous direction
         * on the coarsest level.
         */
        std::vector<double> getAveragedDerivativeOfQuantityWithInhomogeneousXDirectionOnCoarsestLevel(
            const std::vector<std::string>& quantity_names,
            const std::vector<int>& component_indices,
            const std::vector<bool>& use_reciprocal,
            const int derivative_direction,
            const int num_ghosts_derivative,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute averaged value with only x-direction as inhomogeneous direction.
         */
        std::vector<double> getAveragedQuantityWithInhomogeneousXDirection(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute averaged reciprocal of value with only x-direction as inhomogeneous direction.
         */
        std::vector<double> getAveragedReciprocalOfQuantityWithInhomogeneousXDirection(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute averaged value (on product of variables) with only x-direction as inhomogeneous direction.
         */
        std::vector<double> getAveragedQuantityWithInhomogeneousXDirection(
            const std::vector<std::string>& quantity_names,
            const std::vector<int>& component_indices,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute averaged value (on product of variables) with only x-direction as inhomogeneous direction.
         */
        std::vector<double> getAveragedQuantityWithInhomogeneousXDirection(
            const std::vector<std::string>& quantity_names,
            const std::vector<int>& component_indices,
            const std::vector<bool>& use_reciprocal,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute averaged value (on product of variable derivatives) with only x direction as inhomogeneous direction.
         */
        std::vector<double> getAveragedQuantityWithInhomogeneousXDirection(
            const std::vector<std::string>& quantity_names,
            const std::vector<int>& component_indices,
            const std::vector<bool>& use_derivative,
            const std::vector<int>& derivative_directions,
            const int num_ghosts_derivative,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute averaged value (on product of variable derivatives) with only x direction as inhomogeneous direction.
         */
        std::vector<double> getAveragedQuantityWithInhomogeneousXDirection(
            const std::vector<std::string>& quantity_names,
            const std::vector<int>& component_indices,
            const std::vector<bool>& use_derivative,
            const std::vector<int>& derivative_directions,
            const std::vector<bool>& use_reciprocal,
            const int num_ghosts_derivative,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute averaged derivative of value (on product of variables) with only x direction as inhomogeneous direction.
         */
        std::vector<double> getAveragedDerivativeOfQuantityWithInhomogeneousXDirection(
            const std::vector<std::string>& quantity_names,
            const std::vector<int>& component_indices,
            const std::vector<bool>& use_reciprocal,
            const int derivative_direction,
            const int num_ghosts_derivative,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute averaged value with only y-direction as inhomogeneous direction.
         */
        std::vector<double> getAveragedQuantityWithInhomogeneousYDirection(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute averaged reciprocal of value with only y-direction as inhomogeneous direction.
         */
        std::vector<double> getAveragedReciprocalOfQuantityWithInhomogeneousYDirection(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute averaged value (on product of variables) with only y-direction as inhomogeneous direction.
         */
        std::vector<double> getAveragedQuantityWithInhomogeneousYDirection(
            const std::vector<std::string>& quantity_names,
            const std::vector<int>& component_indices,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute averaged value (on product of variables) with only y-direction as inhomogeneous direction.
         */
        std::vector<double> getAveragedQuantityWithInhomogeneousYDirection(
            const std::vector<std::string>& quantity_names,
            const std::vector<int>& component_indices,
            const std::vector<bool>& use_reciprocal,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute averaged value (on product of variable derivatives) with only y direction as inhomogeneous direction.
         */
        std::vector<double> getAveragedQuantityWithInhomogeneousYDirection(
            const std::vector<std::string>& quantity_names,
            const std::vector<int>& component_indices,
            const std::vector<bool>& use_derivative,
            const std::vector<int>& derivative_directions,
            const int num_ghosts_derivative,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute averaged value (on product of variable derivatives) with only y direction as inhomogeneous direction.
         */
        std::vector<double> getAveragedQuantityWithInhomogeneousYDirection(
            const std::vector<std::string>& quantity_names,
            const std::vector<int>& component_indices,
            const std::vector<bool>& use_derivative,
            const std::vector<int>& derivative_directions,
            const std::vector<bool>& use_reciprocal,
            const int num_ghosts_derivative,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute averaged derivative of value (on product of variables) with only y direction as inhomogeneous direction.
         */
        std::vector<double> getAveragedDerivativeOfQuantityWithInhomogeneousYDirection(
            const std::vector<std::string>& quantity_names,
            const std::vector<int>& component_indices,
            const std::vector<bool>& use_reciprocal,
            const int derivative_direction,
            const int num_ghosts_derivative,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute averaged value with only z-direction as inhomogeneous direction.
         */
        std::vector<double> getAveragedQuantityWithInhomogeneousZDirection(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute averaged reciprocal of value with only z-direction as inhomogeneous direction.
         */
        std::vector<double> getAveragedReciprocalOfQuantityWithInhomogeneousZDirection(
            const std::string quantity_name,
            const int component_idx,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute averaged value (on product of variables) with only z-direction as inhomogeneous direction.
         */
        std::vector<double> getAveragedQuantityWithInhomogeneousZDirection(
            const std::vector<std::string>& quantity_names,
            const std::vector<int>& component_indices,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute averaged value (on product of variables) with only z-direction as inhomogeneous direction.
         */
        std::vector<double> getAveragedQuantityWithInhomogeneousZDirection(
            const std::vector<std::string>& quantity_names,
            const std::vector<int>& component_indices,
            const std::vector<bool>& use_reciprocal,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute averaged value (on product of variable derivatives) with only z direction as inhomogeneous direction.
         */
        std::vector<double> getAveragedQuantityWithInhomogeneousZDirection(
            const std::vector<std::string>& quantity_names,
            const std::vector<int>& component_indices,
            const std::vector<bool>& use_derivative,
            const std::vector<int>& derivative_directions,
            const int num_ghosts_derivative,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute averaged value (on product of variable derivatives) with only z direction as inhomogeneous direction.
         */
        std::vector<double> getAveragedQuantityWithInhomogeneousZDirection(
            const std::vector<std::string>& quantity_names,
            const std::vector<int>& component_indices,
            const std::vector<bool>& use_derivative,
            const std::vector<int>& derivative_directions,
            const std::vector<bool>& use_reciprocal,
            const int num_ghosts_derivative,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Compute averaged derivative of value (on product of variables) with only z direction as inhomogeneous direction.
         */
        std::vector<double> getAveragedDerivativeOfQuantityWithInhomogeneousZDirection(
            const std::vector<std::string>& quantity_names,
            const std::vector<int>& component_indices,
            const std::vector<bool>& use_reciprocal,
            const int derivative_direction,
            const int num_ghosts_derivative,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
    private:
        
};

#endif /* FLOW_MODEL_MPI_HELPER_AVERAGE_HPP */
