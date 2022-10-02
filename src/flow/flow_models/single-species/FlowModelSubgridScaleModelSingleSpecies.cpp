#include "flow/flow_models/single-species/FlowModelSubgridScaleModelSingleSpecies.hpp"

FlowModelSubgridScaleModelSingleSpecies::FlowModelSubgridScaleModelSingleSpecies(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_species,
    const HAMERS_SHARED_PTR<tbox::Database>& subgrid_scale_model_db):
        FlowModelSubgridScaleModel(
            object_name,
            dim,
            grid_geometry,
            num_species,
            2 + dim.getValue(),
            subgrid_scale_model_db)
{
    d_constant_sgs = double(0.08);
    d_species_Pr_t = double(0.9);
    
    d_constant_sgs = subgrid_scale_model_db->getDoubleWithDefault("constant_sgs",   d_constant_sgs);
    d_constant_sgs = subgrid_scale_model_db->getDoubleWithDefault("d_constant_sgs", d_constant_sgs);
    
    d_species_Pr_t = subgrid_scale_model_db->getDoubleWithDefault("species_Pr_t",   d_species_Pr_t);
    d_species_Pr_t = subgrid_scale_model_db->getDoubleWithDefault("d_species_Pr_t", d_species_Pr_t);
    
    if (subgrid_scale_model_db->keyExists("species_c_p"))
    {
        d_species_c_p = subgrid_scale_model_db->getDouble("species_c_p");
    }
    else if (subgrid_scale_model_db->keyExists("d_species_c_p"))
    {
        d_species_c_p = subgrid_scale_model_db->getDouble("d_species_c_p");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'species_c_p'/'d_species_c_p' found in data for subgrid scale model."
            << std::endl);
    }
}


/*
 * Return names of different derived variables required to register.
 */
std::vector<std::string>
FlowModelSubgridScaleModelSingleSpecies::getDerivedVariablesToRegister() const
{
    std::vector<std::string> var_to_register;
    var_to_register.reserve(2);
    
    var_to_register.push_back("DENSITY");
    var_to_register.push_back("VELOCITY");
    
    return var_to_register;
}


/*
 * Return different derived variables required for interpolation to sides.
 */
void
FlowModelSubgridScaleModelSingleSpecies::getDerivedVariablesForInterpolationToSideData(
    std::vector<std::string>& var_to_interpolate,
    std::vector<int>& var_to_interpolate_component_idx) const
{
    var_to_interpolate.resize(1);
    var_to_interpolate_component_idx.resize(1);
    
    var_to_interpolate.push_back("DENSITY");
    var_to_interpolate_component_idx.push_back(0);
}


/*
 * Get the variables for the derivatives used at computing subgrid scale diffusivity/viscosity at sides.
 */
void
FlowModelSubgridScaleModelSingleSpecies::getCellDataOfVariablesForSideDerivativeForSubgridScaleViscosity(
    std::vector<std::string>& derivative_var_data_str,
    std::vector<int>& derivative_var_component_idx,
    const DIRECTION::TYPE& side_direction,
    const DIRECTION::TYPE& derivative_direction)
{
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSubgridScaleModelSingleSpecies::"
            << "getCellDataOfVariablesForSideDerivativeForSubgridScaleViscosity() "
            << "not implemented for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSubgridScaleModelSingleSpecies::"
            << "getCellDataOfVariablesForSideDerivativeForSubgridScaleViscosity() "
            << "not implemented for two-dimensional problem."
            << std::endl);
    }
    
    NULL_USE(side_direction);
    NULL_USE(derivative_direction);
    
    derivative_var_data_str.resize(3);
    derivative_var_component_idx.resize(3);
    
    derivative_var_data_str[0] = "VELOCITY";
    derivative_var_component_idx[0] = 0;
    
    derivative_var_data_str[1] = "VELOCITY";
    derivative_var_component_idx[1] = 1;
    
    derivative_var_data_str[2] = "VELOCITY";
    derivative_var_component_idx[2] = 2;
}


/*
 * Modify the side data of the diffusivities/viscosities at sides with subgrid scale diffusivity/viscosity.
 */
void
FlowModelSubgridScaleModelSingleSpecies::updateSideDataOfDiffusiveFluxDiffusivities(
    std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& var_data_for_diffusivities,
    const std::map<DIRECTION::TYPE, std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > >& derivatives,
    const DIRECTION::TYPE& side_direction)
{
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSubgridScaleModelSingleSpecies::"
            << "updateSideDataOfDiffusiveFluxDiffusivities() "
            << "not implemented for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSubgridScaleModelSingleSpecies::"
            << "updateSideDataOfDiffusiveFluxDiffusivities()\n"
            << "Not implemented for two-dimensional problem."
            << std::endl);
    }
    
    if (derivatives.find(DIRECTION::X_DIRECTION) == derivatives.end())
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSubgridScaleModelSingleSpecies::"
            << "updateSideDataOfDiffusiveFluxDiffusivities()\n"
            << "Cannot find derivatives in x-direction in the map of derivatives."
            << std::endl);
    }
    
    if (derivatives.find(DIRECTION::Y_DIRECTION) == derivatives.end())
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSubgridScaleModelSingleSpecies::"
            << "updateSideDataOfDiffusiveFluxDiffusivities()\n"
            << "Cannot find derivatives in y-direction in the map of derivatives."
            << std::endl);
    }
    
    if (derivatives.find(DIRECTION::Z_DIRECTION) == derivatives.end())
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSubgridScaleModelSingleSpecies::"
            << "updateSideDataOfDiffusiveFluxDiffusivities()\n"
            << "Cannot find derivatives in z-direction in the map of derivatives."
            << std::endl);
    }
    
    std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > derivatives_x = derivatives.find(DIRECTION::X_DIRECTION)->second;
    std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > derivatives_y = derivatives.find(DIRECTION::Y_DIRECTION)->second;
    std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > derivatives_z = derivatives.find(DIRECTION::Z_DIRECTION)->second;
    
    TBOX_ASSERT(static_cast<int>(derivatives_x.size()) == 3);
    TBOX_ASSERT(static_cast<int>(derivatives_y.size()) == 3);
    TBOX_ASSERT(static_cast<int>(derivatives_z.size()) == 3);
    
    
    
    
}


/*
 * Put the characteristics of this class into the restart database.
 */
void
FlowModelSubgridScaleModelSingleSpecies::putToRestart(
    const HAMERS_SHARED_PTR<tbox::Database>& restart_subgrid_scale_model_db) const
{
    putToRestartBase(restart_subgrid_scale_model_db);
    
    restart_subgrid_scale_model_db->putDouble("d_constant_sgs", d_constant_sgs);
    restart_subgrid_scale_model_db->putDouble("d_species_Pr_t", d_species_Pr_t);
    restart_subgrid_scale_model_db->putDouble("d_species_c_p",  d_species_c_p);
}
