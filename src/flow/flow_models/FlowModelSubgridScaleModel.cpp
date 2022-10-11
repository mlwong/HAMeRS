#include "flow/flow_models/FlowModelSubgridScaleModel.hpp"

FlowModelSubgridScaleModel::FlowModelSubgridScaleModel(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_species,
    const int& num_eqn,
    const HAMERS_SHARED_PTR<tbox::Database>& subgrid_scale_model_db):
        d_object_name(object_name),
        d_dim(dim),
        d_grid_geometry(grid_geometry),
        d_num_species(num_species),
        d_num_eqn(num_eqn)
{
    std::string subgrid_scale_model_str;
    if (subgrid_scale_model_db->keyExists("subgrid_scale_model"))
    {
        subgrid_scale_model_str = subgrid_scale_model_db->getString("subgrid_scale_model");
    }
    else if (subgrid_scale_model_db->keyExists("d_subgrid_scale_model"))
    {
        subgrid_scale_model_str = subgrid_scale_model_db->getString("d_subgrid_scale_model");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSubgridScaleModel::FlowModelSubgridScaleModel()\n"
            << "No key 'subgrid_scale_model'/'d_subgrid_scale_model' found in data for subgrid scale model."
            << std::endl);
    }
    
    if (subgrid_scale_model_str == "VREMAN")
    {
        d_subgrid_scale_model_type = SUBGRID_SCALE_MODEL::VREMAN;
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSubgridScaleModel::FlowModelSubgridScaleModel()\n"
            << "Unknown/unsupported subgrid scale model with string: '" << subgrid_scale_model_str << "'."
            << std::endl);
    }
}


/*
 * Return names of different derived variables required to register.
 */
std::vector<std::string>
FlowModelSubgridScaleModel::getDerivedVariablesToRegister() const
{
    std::vector<std::string> var_to_register;
    return var_to_register;
}


/*
 * Return different derived variables required for interpolation.
 */
void
FlowModelSubgridScaleModel::getDerivedVariablesForInterpolationToSideData(
    std::vector<std::string>& var_to_interpolate,
    std::vector<int>& var_to_interpolate_component_idx) const
{
    NULL_USE(var_to_interpolate);
    NULL_USE(var_to_interpolate_component_idx);
}


/*
 * Get the variables for the derivatives used at computing subgrid scale diffusivity/viscosity at sides.
 */
void
FlowModelSubgridScaleModel::getCellDataOfVariablesForSideDerivativeForSubgridScaleViscosity(
    std::vector<std::string>& derivative_var_data_str,
    std::vector<int>& derivative_var_component_idx,
    const DIRECTION::TYPE& side_direction,
    const DIRECTION::TYPE& derivative_direction)
{
    NULL_USE(derivative_var_data_str);
    NULL_USE(derivative_var_component_idx);
    NULL_USE(side_direction);
    NULL_USE(derivative_direction);
}


/*
 * Modify the side data of the diffusivities/viscosities at sides with subgrid scale diffusivity/viscosity.
 */
void
FlowModelSubgridScaleModel::updateSideDataOfDiffusiveFluxDiffusivities(
    std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& var_data_for_diffusivities,
    const std::map<DIRECTION::TYPE, std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > >& derivatives,
    const DIRECTION::TYPE& side_direction,
    const hier::Patch& patch)
{
    NULL_USE(var_data_for_diffusivities);
    NULL_USE(derivatives);
    NULL_USE(side_direction);
    NULL_USE(patch);
}


/*
 * Put the characteristics of this class into the restart database.
 */
void
FlowModelSubgridScaleModel::putToRestart(
    const HAMERS_SHARED_PTR<tbox::Database>& restart_subgrid_scale_model_db) const
{
    putToRestartBase(restart_subgrid_scale_model_db);
}


/*
 * Put the characteristics of base class into the restart database.
 */
void
FlowModelSubgridScaleModel::putToRestartBase(
    const HAMERS_SHARED_PTR<tbox::Database>& restart_subgrid_scale_model_db) const
{
    if (d_subgrid_scale_model_type == SUBGRID_SCALE_MODEL::VREMAN)
    {
        restart_subgrid_scale_model_db->putString("d_subgrid_scale_model", "VREMAN");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSubgridScaleModel::putToRestartBase()\n"
            << "Unknown/unsupported subgrid scale model."
            << std::endl);
    }
}
