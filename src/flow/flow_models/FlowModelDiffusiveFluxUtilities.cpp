#include "flow/flow_models/FlowModelDiffusiveFluxUtilities.hpp"

FlowModelDiffusiveFluxUtilities::FlowModelDiffusiveFluxUtilities(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_species,
    const int& num_eqn):
        d_object_name(object_name),
        d_dim(dim),
        d_grid_geometry(grid_geometry),
        d_num_species(num_species),
        d_num_eqn(num_eqn),
        d_derived_cell_data_computed(false),
        d_num_subghosts_diffusivities(-hier::IntVector::getOne(d_dim)),
        d_subghost_box_diffusivities(hier::Box::getEmptyBox(dim)),
        d_subghostcell_dims_diffusivities(hier::IntVector::getZero(d_dim)),
        d_cell_data_diffusivities_computed(false)
{}


/*
 * Get the variables for the derivatives in the diffusive fluxes.
 */
void
FlowModelDiffusiveFluxUtilities::getCellDataOfDiffusiveFluxVariablesForDerivative(
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
        << ": FlowModelDiffusiveFluxUtilities::getCellDataOfDiffusiveFluxVariablesForDerivativesAtNodes()\n"
        << "Function is not yet implemented!"
        << std::endl);
}


/*
 * Get the diffusivities in the diffusive flux.
 */
void
FlowModelDiffusiveFluxUtilities::getCellDataOfDiffusiveFluxDiffusivities(
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
        << ": FlowModelDiffusiveFluxUtilities::getCellDataOfDiffusiveFluxDiffusivities()\n"
        << "Function is not yet implemented!"
        << std::endl);
}
