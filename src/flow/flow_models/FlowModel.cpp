#include "flow/flow_models/FlowModel.hpp"

/*
 * Register the required variables for the computation of diffusive fluxes in the registered patch.
 */
void
FlowModel::registerDiffusiveFluxes(const hier::IntVector& num_subghosts)
{
    NULL_USE(num_subghosts);
    
    TBOX_ERROR(d_object_name
        << ": FlowModel::registerDiffusiveFluxes()\n"
        << "Function is not yet implemented!"
        << std::endl);
}


/*
 * Get the variables for the derivatives in the diffusive fluxes.
 */
void
FlowModel::getDiffusiveFluxVariablesForDerivative(
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
        << ": FlowModel::getDiffusiveFluxVariablesForDerivativesAtNodes()\n"
        << "Function is not yet implemented!"
        << std::endl);
}


/*
 * Get the diffusivities in the diffusive flux.
 */
void
FlowModel::getDiffusiveFluxDiffusivities(
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
        << ": FlowModel::getDiffusiveFluxDiffusivities()\n"
        << "Function is not yet implemented!"
        << std::endl);
}


/*
 * Setup the Riemann solver object.
 */
void
FlowModel::setupRiemannSolver()
{
    d_flow_model_riemann_solver->setFlowModel(shared_from_this());
}


/*
 * Setup the statistics utilties object.
 */
void
FlowModel::setupStatisticsUtilities()
{
    d_flow_model_statistics_utilities->setFlowModel(shared_from_this());
}
