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
        d_cell_data_computed_diffusivities(false)
{}


/*
 * Register different derived variables related to this class in the registered patch. The
 * derived variables to be registered are given as entires in a map of the variable name to
 * the number of sub-ghost cells required.
 */
void
FlowModelDiffusiveFluxUtilities::registerDerivedVariables(
    const std::unordered_map<std::string, hier::IntVector>& num_subghosts_of_data)
{
    NULL_USE(num_subghosts_of_data);
    
    TBOX_ERROR(d_object_name
        << ": FlowModelDiffusiveFluxUtilities::registerDerivedVariables()\n"
        << "Function is not yet implemented!"
        << std::endl);
}


/*
 * Register the required variables for the computation of diffusive fluxes in the registered patch.
 */
void
FlowModelDiffusiveFluxUtilities::registerDerivedVariablesForDiffusiveFluxes(
    const hier::IntVector& num_subghosts)
{
    NULL_USE(num_subghosts);
    
    TBOX_ERROR(d_object_name
        << ": FlowModelDiffusiveFluxUtilities::registerDerivedVariablesForDiffusiveFluxes()\n"
        << "Function is not yet implemented!"
        << std::endl);
}


/*
 * Allocate memory for cell data of different registered derived variables related to this
 * class in the registered patch.
 */
void
FlowModelDiffusiveFluxUtilities::allocateMemoryForDerivedCellData()
{
    TBOX_ERROR(d_object_name
        << ": FlowModelDiffusiveFluxUtilities::allocateMemoryForDerivedCellData()\n"
        << "Function is not yet implemented!"
        << std::endl);
}


/*
 * Clear cell data of different derived variables related to this class in the registered patch.
 */
void
FlowModelDiffusiveFluxUtilities::clearCellData()
{
}


/*
 * Compute cell data of different registered derived variables related to this class.
 */
void
FlowModelDiffusiveFluxUtilities::computeDerivedCellData()
{
    TBOX_ERROR(d_object_name
        << ": FlowModelDiffusiveFluxUtilities::computeDerivedCellData()\n"
        << "Function is not yet implemented!"
        << std::endl);
}


/*
 * Get the cell data of one cell variable related to this class in the registered patch.
 */
boost::shared_ptr<pdat::CellData<double> >
FlowModelDiffusiveFluxUtilities::getCellData(const std::string& variable_key)
{
    NULL_USE(variable_key);
    
    TBOX_ERROR(d_object_name
        << ": FlowModelDiffusiveFluxUtilities::getCellData()\n"
        << "Function is not yet implemented!"
        << std::endl);
    
    return nullptr;
}


/*
 * Get the cell data of different cell variables related to this class in the registered patch.
 */
std::vector<boost::shared_ptr<pdat::CellData<double> > >
FlowModelDiffusiveFluxUtilities::getCellData(
    const std::vector<std::string>& variable_keys)
{
    NULL_USE(variable_keys);
    
    TBOX_ERROR(d_object_name
        << ": FlowModelDiffusiveFluxUtilities::getCellData()\n"
        << "Function is not yet implemented!"
        << std::endl);
    
    std::vector<boost::shared_ptr<pdat::CellData<double> > > cell_data(
        static_cast<int>(variable_keys.size()));
    
    return cell_data;
}


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
