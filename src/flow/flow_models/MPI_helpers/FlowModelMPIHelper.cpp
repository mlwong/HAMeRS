#include "flow/flow_models/MPI_helpers/FlowModelMPIHelper.hpp"

/*
 * Register the patch and data in the flow model. Also, set up the diffusive flux utilities object if needed.
 */
void FlowModelMPIHelper::setupFlowModelAndRegisterPatchWithDataContext(
    const hier::Patch& patch,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    if (d_use_diffusive_flux_utilities)
    {
        d_flow_model->setupDiffusiveFluxUtilities();
        d_diffusive_flux_utilities = d_flow_model->getFlowModelDiffusiveFluxUtilities();
    }
    
    d_flow_model->registerPatchWithDataContext(patch, data_context);
}


/*
 * Register different derived variables in the flow model with the registered patch. 
 */
void FlowModelMPIHelper::registerDerivedVariables(
    const std::unordered_map<std::string, hier::IntVector>& num_subghosts_of_data)
{
    d_flow_model->registerDerivedVariables(num_subghosts_of_data);
    
    if (d_use_diffusive_flux_utilities)
    {
        d_diffusive_flux_utilities->registerDerivedVariables(num_subghosts_of_data);
    }
}


/*
 * Allocate memory in the flow model for cell data of different registered derived variables.
 */
void FlowModelMPIHelper::allocateMemoryForDerivedCellData()
{
    d_flow_model->allocateMemoryForDerivedCellData();
    
    if (d_use_diffusive_flux_utilities)
    {
        d_diffusive_flux_utilities->allocateMemoryForDerivedCellData();
    }
}


/*
 * Compute the cell data of different registered derived variables in the flow model with the registered data
 * context.
 */
void FlowModelMPIHelper::computeDerivedCellData()
{
    d_flow_model->computeDerivedCellData();
    
    if (d_use_diffusive_flux_utilities)
    {
        d_diffusive_flux_utilities->computeDerivedCellData();
    }
}


/*
 * Get the cell data of one cell variable in the registered patch.
 */
HAMERS_SHARED_PTR<pdat::CellData<Real> >
FlowModelMPIHelper::getCellData(const std::string& variable_key)
{
    HAMERS_SHARED_PTR<pdat::CellData<Real> > cell_data;
    
    if (d_use_diffusive_flux_utilities)
    {
        cell_data = d_diffusive_flux_utilities->getCellData(variable_key);
    }
    
    if (cell_data == nullptr)
    {
        cell_data = d_flow_model->getCellData(variable_key);
    }
    
    return cell_data;
}


/*
 * Unregister the registered patch in the flow model.
 */
void FlowModelMPIHelper::unregisterPatch()
{
    d_flow_model->unregisterPatch();
    d_diffusive_flux_utilities.reset();
}
