#include "flow/flow_models/FlowModelDiffusiveFluxUtilities.hpp"

FlowModelDiffusiveFluxUtilities::FlowModelDiffusiveFluxUtilities(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_species,
    const int& num_eqn,
    const HAMERS_SHARED_PTR<tbox::Database>& flow_model_db):
        d_object_name(object_name),
        d_dim(dim),
        d_grid_geometry(grid_geometry),
        d_num_species(num_species),
        d_num_eqn(num_eqn),
        d_derived_cell_data_computed(false),
        d_num_subghosts_diffusivities(-hier::IntVector::getOne(d_dim)),
        d_subghost_box_diffusivities(hier::Box::getEmptyBox(dim)),
        d_subghostcell_dims_diffusivities(hier::IntVector::getZero(d_dim)),
        d_cell_data_computed_diffusivities(false),
        d_side_data_diffusivities_computed(false),
        d_need_side_diffusivities(false),
        d_use_subgrid_scale_model(false)
{
    if (flow_model_db->keyExists("use_subgrid_scale_model"))
    {
        d_use_subgrid_scale_model = flow_model_db->getBool("use_subgrid_scale_model");
    }
    else if (flow_model_db->keyExists("d_use_subgrid_scale_model"))
    {
        d_use_subgrid_scale_model = flow_model_db->getBool("d_use_subgrid_scale_model");
    }
}


/*
 * Register different derived variables related to this class in the registered patch. The
 * derived variables to be registered are given as entries in a map of the variable name to
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
    const hier::IntVector& num_subghosts,
    const bool need_side_diffusivities)
{
    NULL_USE(num_subghosts);
    NULL_USE(need_side_diffusivities);
    
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
 * Allocate memory for side data of the diffusivities.
 */
void
FlowModelDiffusiveFluxUtilities::allocateMemoryForSideDataOfDiffusiveFluxDiffusivities()
{
    TBOX_ERROR(d_object_name
        << ": FlowModelDiffusiveFluxUtilities::allocateMemoryForSideDataOfDiffusiveFluxDiffusivities()\n"
        << "Function is not yet implemented!"
        << std::endl);
}


/*
 * Clear cell and side data of different derived variables related to this class in the registered patch.
 */
void
FlowModelDiffusiveFluxUtilities::clearCellAndSideData()
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
HAMERS_SHARED_PTR<pdat::CellData<double> >
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
std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >
FlowModelDiffusiveFluxUtilities::getCellData(
    const std::vector<std::string>& variable_keys)
{
    NULL_USE(variable_keys);
    
    TBOX_ERROR(d_object_name
        << ": FlowModelDiffusiveFluxUtilities::getCellData()\n"
        << "Function is not yet implemented!"
        << std::endl);
    
    std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > cell_data(
        static_cast<int>(variable_keys.size()));
    
    return cell_data;
}


/*
 * Get the variables for the derivatives in the diffusive fluxes.
 */
void
FlowModelDiffusiveFluxUtilities::getCellDataOfDiffusiveFluxVariablesForDerivative(
    std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& derivative_var_data,
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
    std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& diffusivities_data,
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


/*
 * Get the cell data that needs interpolation to midpoints for computing side data of diffusivities in the
 * diffusive flux.
 */
void
FlowModelDiffusiveFluxUtilities::getCellDataForInterpolationToSideDataForDiffusiveFluxDiffusivities(
    std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& var_data_for_diffusivities,
    std::vector<int>& var_data_for_diffusivities_component_idx)
{
    NULL_USE(var_data_for_diffusivities);
    NULL_USE(var_data_for_diffusivities_component_idx);
    
    TBOX_ERROR(d_object_name
        << ": FlowModelDiffusiveFluxUtilities::getCellDataForInterpolationToSideDataForDiffusiveFluxDiffusivities()\n"
        << "Function is not yet implemented!"
        << std::endl);
}


/*
 * Compute the side data of the diffusivities in the diffusive flux with the interpolated side data.
 */
void
FlowModelDiffusiveFluxUtilities::computeSideDataOfDiffusiveFluxDiffusivities(
    const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& var_data_for_diffusivities)
{
    NULL_USE(var_data_for_diffusivities);
    
    TBOX_ERROR(d_object_name
        << ": FlowModelDiffusiveFluxUtilities::computeSideDataOfDiffusiveFluxDiffusivities()\n"
        << "Function is not yet implemented!"
        << std::endl);
}


/*
 * Get the side data of the diffusivities in the diffusive fluxa.
 */
void
FlowModelDiffusiveFluxUtilities::getSideDataOfDiffusiveFluxDiffusivities(
    std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > >& diffusivities_data,
    std::vector<std::vector<int> >& diffusivities_component_idx,
    const DIRECTION::TYPE& flux_direction,
    const DIRECTION::TYPE& derivative_direction)
{
    NULL_USE(diffusivities_data);
    NULL_USE(diffusivities_component_idx);
    NULL_USE(flux_direction);
    NULL_USE(derivative_direction);
    
    TBOX_ERROR(d_object_name
        << ": FlowModelDiffusiveFluxUtilities::getSideDataOfDiffusiveFluxDiffusivities()\n"
        << "Function is not yet implemented!"
        << std::endl);
}
