#include "flow/flow_models/five-eqn_Allaire/FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire.hpp"

/*
 * Register the required variables for the computation of diffusive fluxes in the registered patch.
 */
void
FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::registerDiffusiveFluxes(const hier::IntVector& num_subghosts)
{
    NULL_USE(num_subghosts);
    
    TBOX_ERROR(d_object_name
        << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::registerDiffusiveFluxes()\n"
        << "Function is not yet implemented!"
        << std::endl);
}


/*
 * Get the variables for the derivatives in the diffusive fluxes.
 */
void
FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::getDiffusiveFluxVariablesForDerivative(
    std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& derivative_var_data,
    std::vector<std::vector<int> >& derivative_var_component_idx,
    const DIRECTION::TYPE& flux_direction,
    const DIRECTION::TYPE& derivative_direction)
{
    NULL_USE(derivative_var_data);
    NULL_USE(derivative_var_component_idx);
    NULL_USE(flux_direction);
    NULL_USE(derivative_direction);
    
    TBOX_ERROR(d_object_name
        << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::getDiffusiveFluxVariablesForDerivativesAtNodes()\n"
        << "Function is not yet implemented!"
        << std::endl);
}


/*
 * Get the diffusivities in the diffusive flux.
 */
void
FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::getDiffusiveFluxDiffusivities(
    std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& diffusivities_data,
    std::vector<std::vector<int> >& diffusivities_component_idx,
    const DIRECTION::TYPE& flux_direction,
    const DIRECTION::TYPE& derivative_direction)
{
    NULL_USE(diffusivities_data);
    NULL_USE(diffusivities_component_idx);
    NULL_USE(flux_direction);
    NULL_USE(derivative_direction);
    
    TBOX_ERROR(d_object_name
        << ": FlowModelDiffusiveFluxUtilitiesFiveEqnAllaire::getDiffusiveFluxDiffusivities()\n"
        << "Function is not yet implemented!"
        << std::endl);
}
