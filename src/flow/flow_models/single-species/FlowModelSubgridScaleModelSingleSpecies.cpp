#include "flow/flow_models/single-species/FlowModelSubgridScaleModelSingleSpecies.hpp"

/*
 * Return names of different derived variables related to this class in the registered patch.
 */
std::vector<std::string>
FlowModelSubgridScaleModelSingleSpecies::getDerivedVariablesToRegister() const
{
    
}


/*
 * Get the variables for the derivatives used at computing subgrid scale diffusivity/viscosity at midpoints.
 */
void
FlowModelSubgridScaleModelSingleSpecies::getCellDataOfVariablesForDerivativeForSubgridScaleViscosity(
    std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& derivative_var_data,
    std::vector<std::vector<int> >& derivative_var_component_idx,
    const DIRECTION::TYPE& midpoint_direction,
    const DIRECTION::TYPE& derivative_direction)
{
    
}


/*
 * Modify the side data of the diffusivities/viscosities with subgrid scale diffusivity/viscosity.
 */
void
FlowModelSubgridScaleModelSingleSpecies::updateSideDataOfDiffusiveFluxDiffusivities(
    std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& var_data_for_diffusivities,
    const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > >& derivatives_x,
    const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > >& derivatives_y,
    const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > > >& derivatives_z)
{
    
}
