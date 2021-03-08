#include "flow/flow_models/MPI_helpers/FlowModelMPIHelper.hpp"

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


void FlowModelMPIHelper::registerDerivedVariables(
    const std::unordered_map<std::string, hier::IntVector>& num_subghosts_of_data)
{
    d_flow_model->registerDerivedVariables(num_subghosts_of_data);
    
    if (d_use_diffusive_flux_utilities)
    {
        d_diffusive_flux_utilities->registerDerivedVariables(num_subghosts_of_data);
    }
}


void FlowModelMPIHelper::allocateMemoryForDerivedCellData()
{
    d_flow_model->allocateMemoryForDerivedCellData();
    
    if (d_use_diffusive_flux_utilities)
    {
        d_diffusive_flux_utilities->allocateMemoryForDerivedCellData();
    }
}


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
HAMERS_SHARED_PTR<pdat::CellData<double> >
FlowModelMPIHelper::getCellData(const std::string& variable_key)
{
    HAMERS_SHARED_PTR<pdat::CellData<double> > cell_data;
    
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


void FlowModelMPIHelper::unregisterPatch()
{
    d_flow_model->unregisterPatch();
    d_diffusive_flux_utilities.reset();
}
