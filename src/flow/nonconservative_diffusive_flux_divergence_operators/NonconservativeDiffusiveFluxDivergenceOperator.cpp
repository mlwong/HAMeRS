#include "flow/nonconservative_diffusive_flux_divergence_operators/NonconservativeDiffusiveFluxDivergenceOperator.hpp"

#include "SAMRAI/geom/CartesianPatchGeometry.h"

#include <map>

/*
 * Compute the non-conservative diffusive flux divergence on a patch.
 */
void
NonconservativeDiffusiveFluxDivergenceOperator::computeNonconservativeDiffusiveFluxDivergenceOnPatch(
    hier::Patch& patch,
    const HAMERS_SHARED_PTR<pdat::CellVariable<double> >& variable_diffusive_flux_divergence,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double time,
    const double dt,
    const int RK_step_number)
{
    NULL_USE(time);
    NULL_USE(RK_step_number);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(variable_diffusive_flux_divergence);
#endif
    
    d_flow_model->setupDiffusiveFluxUtilities();
    
    HAMERS_SHARED_PTR<FlowModelDiffusiveFluxUtilities> diffusive_flux_utilities =
        d_flow_model->getFlowModelDiffusiveFluxUtilities();
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // diffusive ghost cells.
    hier::Box diff_ghost_box = interior_box;
    diff_ghost_box.grow(d_num_diff_ghosts);
    const hier::IntVector diff_ghostcell_dims = diff_ghost_box.numberCells();
    
    // Get the cell data of diffusive flux divergence.
    HAMERS_SHARED_PTR<pdat::CellData<double> > diffusive_flux_divergence(
        HAMERS_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
            patch.getPatchData(variable_diffusive_flux_divergence, data_context)));
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(diffusive_flux_divergence);
    TBOX_ASSERT(diffusive_flux_divergence->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
#endif
    
    // Initialize the data of diffusive flux to zero.
    diffusive_flux_divergence->fillAll(double(0));
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Register the patch and derived cell variables in the flow model and compute the corresponding cell data.
         */
        
        d_flow_model->registerPatchWithDataContext(patch, data_context);
        
        diffusive_flux_utilities->registerDerivedVariablesForDiffusiveFluxes(d_num_diff_ghosts);
        
        diffusive_flux_utilities->allocateMemoryForDerivedCellData();
        
        diffusive_flux_utilities->computeDerivedCellData();
        
        /*
         * Delcare containers for computing flux derivatives in different directions.
         */
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > var_data_x;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > diffusivities_data_x;
        
        std::vector<std::vector<int> > var_component_idx_x;
        
        std::vector<std::vector<int> > diffusivities_component_idx_x;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > var_derivative_x;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > var_derivative_xx;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > diffusivities_derivative_x;
        
        std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivative_x_computed;
        
        std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivative_xx_computed;
        
        /*
         * Compute the derivatives for diffusive flux in x-direction.
         */
        
        // Get the variables for the terms in x-direction in the diffusive flux in x-direction.
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxVariablesForDerivative(
            var_data_x,
            var_component_idx_x,
            DIRECTION::X_DIRECTION,
            DIRECTION::X_DIRECTION);
        
        // Get the diffusivities for the terms in x-direction in the diffusive flux in x-direction.
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxDiffusivities(
            diffusivities_data_x,
            diffusivities_component_idx_x,
            DIRECTION::X_DIRECTION,
            DIRECTION::X_DIRECTION);
        
        TBOX_ASSERT(static_cast<int>(var_data_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_x.size()) == d_num_eqn);
        
        // Compute the first derivatives of variables in x-direction.
        computeFirstDerivativesInX(
            var_derivative_x,
            derivative_x_computed,
            var_data_x,
            var_component_idx_x,
            patch);
        
        // Compute the first derivatives of diffusivities in x-direction.
        computeFirstDerivativesInX(
            diffusivities_derivative_x,
            derivative_x_computed,
            diffusivities_data_x,
            diffusivities_component_idx_x,
            patch);
        
        // Compute the second derivatives of variables in x-direction.
        computeSecondDerivativesInX(
            var_derivative_xx,
            derivative_xx_computed,
            var_data_x,
            var_component_idx_x,
            patch);
        
        // Add the derivatives to the divergence of diffusive flux.
        addDerivativeToDivergence(
            diffusive_flux_divergence,
            var_derivative_x,
            var_derivative_xx,
            diffusivities_data_x,
            diffusivities_derivative_x,
            var_component_idx_x,
            diffusivities_component_idx_x,
            patch,
            dt);
        
        var_data_x.clear();
        
        diffusivities_data_x.clear();
        
        var_component_idx_x.clear();
        
        diffusivities_component_idx_x.clear();
        
        var_derivative_x.clear();
        diffusivities_derivative_x.clear();
        var_derivative_xx.clear();
        
        /*
         * Unregister the patch and data of all registered derived cell variables in the flow model.
         */
        
        d_flow_model->unregisterPatch();
    }
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Register the patch and derived cell variables in the flow model and compute the corresponding cell data.
         */
        
        d_flow_model->registerPatchWithDataContext(patch, data_context);
        
        diffusive_flux_utilities->registerDerivedVariablesForDiffusiveFluxes(d_num_diff_ghosts);
        
        diffusive_flux_utilities->allocateMemoryForDerivedCellData();
        
        diffusive_flux_utilities->computeDerivedCellData();
        
        /*
         * Delcare containers for computing flux derivatives in different directions.
         */
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > var_data_x;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > var_data_y;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > diffusivities_data_x;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > diffusivities_data_y;
        
        std::vector<std::vector<int> > var_component_idx_x;
        std::vector<std::vector<int> > var_component_idx_y;
        
        std::vector<std::vector<int> > diffusivities_component_idx_x;
        std::vector<std::vector<int> > diffusivities_component_idx_y;
        
        std::vector<std::vector<int> > var_derivative_component_idx_x;
        std::vector<std::vector<int> > var_derivative_component_idx_y;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > var_derivative_x;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > var_derivative_y;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > var_derivative_xx;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > var_derivative_xy;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > var_derivative_yx;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > var_derivative_yy;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > diffusivities_derivative_x;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > diffusivities_derivative_y;
        
        std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivative_x_computed;
        std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivative_y_computed;
        
        std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivative_xx_computed;
        std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivative_xy_computed;
        std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivative_yx_computed;
        std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivative_yy_computed;
        
        /*
         * (1) Compute the derivatives for diffusive flux in x-direction.
         */
        
        // Get the variables for the terms in x-direction in the diffusive flux in x-direction.
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxVariablesForDerivative(
            var_data_x,
            var_component_idx_x,
            DIRECTION::X_DIRECTION,
            DIRECTION::X_DIRECTION);
        
        // Get the diffusivities for the terms in x-direction in the diffusive flux in x-direction.
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxDiffusivities(
            diffusivities_data_x,
            diffusivities_component_idx_x,
            DIRECTION::X_DIRECTION,
            DIRECTION::X_DIRECTION);
        
        TBOX_ASSERT(static_cast<int>(var_data_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_x.size()) == d_num_eqn);
        
        // Compute the first derivatives of variables in x-direction.
        computeFirstDerivativesInX(
            var_derivative_x,
            derivative_x_computed,
            var_data_x,
            var_component_idx_x,
            patch);
        
        // Compute the first derivatives of diffusivities in x-direction.
        computeFirstDerivativesInX(
            diffusivities_derivative_x,
            derivative_x_computed,
            diffusivities_data_x,
            diffusivities_component_idx_x,
            patch);
        
        // Compute the second derivatives of variables in x-direction.
        computeSecondDerivativesInX(
            var_derivative_xx,
            derivative_xx_computed,
            var_data_x,
            var_component_idx_x,
            patch);
        
        // Add the derivatives to the divergence of diffusive flux.
        addDerivativeToDivergence(
            diffusive_flux_divergence,
            var_derivative_x,
            var_derivative_xx,
            diffusivities_data_x,
            diffusivities_derivative_x,
            var_component_idx_x,
            diffusivities_component_idx_x,
            patch,
            dt);
        
        var_data_x.clear();
        
        diffusivities_data_x.clear();
        
        var_component_idx_x.clear();
        
        diffusivities_component_idx_x.clear();
        
        var_derivative_x.clear();
        diffusivities_derivative_x.clear();
        var_derivative_xx.clear();
        
        // Get the variables for the terms in y-direction in the diffusive flux in x-direction.
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxVariablesForDerivative(
            var_data_y,
            var_component_idx_y,
            DIRECTION::X_DIRECTION,
            DIRECTION::Y_DIRECTION);
        
        // Get the diffusivities for the terms in y-direction in the diffusive flux in x-direction.
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxDiffusivities(
            diffusivities_data_y,
            diffusivities_component_idx_y,
            DIRECTION::X_DIRECTION,
            DIRECTION::Y_DIRECTION);
        
        TBOX_ASSERT(static_cast<int>(var_data_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_y.size()) == d_num_eqn);
        
        // Compute the first derivatives of variables in y-direction.
        computeFirstDerivativesInY(
            var_derivative_y,
            derivative_y_computed,
            var_data_y,
            var_component_idx_y,
            patch);
        
        // Compute the first derivatives of diffusivities in x-direction.
        computeFirstDerivativesInX(
            diffusivities_derivative_x,
            derivative_x_computed,
            diffusivities_data_y,
            diffusivities_component_idx_y,
            patch);
        
        // Compute the mixed derivatives of variables.
        var_derivative_component_idx_y.resize(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            var_derivative_component_idx_y[ei].resize(var_derivative_y[ei].size());
            for (int vi = 0; vi < static_cast<int>(var_derivative_component_idx_y[ei].size()); vi++)
            {
                var_derivative_component_idx_y[ei][vi] = 0;
            }
        }
        
        computeFirstDerivativesInX(
            var_derivative_xy,
            derivative_xy_computed,
            var_derivative_y,
            var_derivative_component_idx_y,
            patch);
        
        // Add the derivatives to the divergence of diffusive flux.
        addDerivativeToDivergence(
            diffusive_flux_divergence,
            var_derivative_y,
            var_derivative_xy,
            diffusivities_data_y,
            diffusivities_derivative_x,
            var_component_idx_y,
            diffusivities_component_idx_y,
            patch,
            dt);
        
        var_data_y.clear();
        
        diffusivities_data_y.clear();
        
        var_component_idx_y.clear();
        
        diffusivities_component_idx_y.clear();
        
        var_derivative_component_idx_y.clear();
        
        var_derivative_y.clear();
        diffusivities_derivative_x.clear();
        var_derivative_xy.clear();
        
        /*
         * (2) Compute the derivatives for diffusive flux in y-direction.
         */
        
        // Get the variables for the terms in x-direction in the diffusive flux in y-direction.
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxVariablesForDerivative(
            var_data_x,
            var_component_idx_x,
            DIRECTION::Y_DIRECTION,
            DIRECTION::X_DIRECTION);
        
        // Get the diffusivities for the terms in x-direction in the diffusive flux in y-direction.
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxDiffusivities(
            diffusivities_data_x,
            diffusivities_component_idx_x,
            DIRECTION::Y_DIRECTION,
            DIRECTION::X_DIRECTION);
        
        TBOX_ASSERT(static_cast<int>(var_data_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_x.size()) == d_num_eqn);
        
        // Compute the first derivatives of variables in x-direction.
        computeFirstDerivativesInX(
            var_derivative_x,
            derivative_x_computed,
            var_data_x,
            var_component_idx_x,
            patch);
        
        // Compute the first derivatives of diffusivities in y-direction.
        computeFirstDerivativesInY(
            diffusivities_derivative_y,
            derivative_y_computed,
            diffusivities_data_x,
            diffusivities_component_idx_x,
            patch);
        
        // Compute the mixed derivatives of variables.
        var_derivative_component_idx_x.resize(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            var_derivative_component_idx_x[ei].resize(var_derivative_x[ei].size());
            for (int vi = 0; vi < static_cast<int>(var_derivative_component_idx_x[ei].size()); vi++)
            {
                var_derivative_component_idx_x[ei][vi] = 0;
            }
        }
        
        computeFirstDerivativesInY(
            var_derivative_yx,
            derivative_yx_computed,
            var_derivative_x,
            var_derivative_component_idx_x,
            patch);
        
        // Add the derivatives to the divergence of diffusive flux.
        addDerivativeToDivergence(
            diffusive_flux_divergence,
            var_derivative_x,
            var_derivative_yx,
            diffusivities_data_x,
            diffusivities_derivative_y,
            var_component_idx_x,
            diffusivities_component_idx_x,
            patch,
            dt);
        
        var_data_x.clear();
        
        diffusivities_data_x.clear();
        
        var_component_idx_x.clear();
        
        diffusivities_component_idx_x.clear();
        
        var_derivative_component_idx_x.clear();
        
        var_derivative_x.clear();
        diffusivities_derivative_y.clear();
        var_derivative_yx.clear();
        
        // Get the variables for the terms in y-direction in the diffusive flux in y-direction.
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxVariablesForDerivative(
            var_data_y,
            var_component_idx_y,
            DIRECTION::Y_DIRECTION,
            DIRECTION::Y_DIRECTION);
        
        // Get the diffusivities for the terms in y-direction in the diffusive flux in y-direction.
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxDiffusivities(
            diffusivities_data_y,
            diffusivities_component_idx_y,
            DIRECTION::Y_DIRECTION,
            DIRECTION::Y_DIRECTION);
        
        TBOX_ASSERT(static_cast<int>(var_data_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_y.size()) == d_num_eqn);
        
        // Compute the first derivatives of variables in y-direction.
        computeFirstDerivativesInY(
            var_derivative_y,
            derivative_y_computed,
            var_data_y,
            var_component_idx_y,
            patch);
        
        // Compute the first derivatives of diffusivities in y-direction.
        computeFirstDerivativesInY(
            diffusivities_derivative_y,
            derivative_y_computed,
            diffusivities_data_y,
            diffusivities_component_idx_y,
            patch);
        
        // Compute the second derivatives of variables in y-direction.
        computeSecondDerivativesInY(
            var_derivative_yy,
            derivative_yy_computed,
            var_data_y,
            var_component_idx_y,
            patch);
        
        // Add the derivatives to the divergence of diffusive flux.
        addDerivativeToDivergence(
            diffusive_flux_divergence,
            var_derivative_y,
            var_derivative_yy,
            diffusivities_data_y,
            diffusivities_derivative_y,
            var_component_idx_y,
            diffusivities_component_idx_y,
            patch,
            dt);
        
        var_data_y.clear();
        
        diffusivities_data_y.clear();
        
        var_component_idx_y.clear();
        
        diffusivities_component_idx_y.clear();
        
        var_derivative_y.clear();
        diffusivities_derivative_y.clear();
        var_derivative_yy.clear();
        
        /*
         * Unregister the patch and data of all registered derived cell variables in the flow model.
         */
        
        d_flow_model->unregisterPatch();
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Register the patch and derived cell variables in the flow model and compute the corresponding cell data.
         */
        
        d_flow_model->registerPatchWithDataContext(patch, data_context);
        
        diffusive_flux_utilities->registerDerivedVariablesForDiffusiveFluxes(d_num_diff_ghosts);
        
        diffusive_flux_utilities->allocateMemoryForDerivedCellData();
        
        diffusive_flux_utilities->computeDerivedCellData();
        
        /*
         * Delcare containers for computing flux derivatives in different directions.
         */
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > var_data_x;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > var_data_y;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > var_data_z;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > diffusivities_data_x;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > diffusivities_data_y;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > diffusivities_data_z;
        
        std::vector<std::vector<int> > var_component_idx_x;
        std::vector<std::vector<int> > var_component_idx_y;
        std::vector<std::vector<int> > var_component_idx_z;
        
        std::vector<std::vector<int> > diffusivities_component_idx_x;
        std::vector<std::vector<int> > diffusivities_component_idx_y;
        std::vector<std::vector<int> > diffusivities_component_idx_z;
        
        std::vector<std::vector<int> > var_derivative_component_idx_x;
        std::vector<std::vector<int> > var_derivative_component_idx_y;
        std::vector<std::vector<int> > var_derivative_component_idx_z;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > var_derivative_x;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > var_derivative_y;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > var_derivative_z;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > var_derivative_xx;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > var_derivative_xy;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > var_derivative_xz;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > var_derivative_yx;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > var_derivative_yy;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > var_derivative_yz;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > var_derivative_zx;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > var_derivative_zy;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > var_derivative_zz;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > diffusivities_derivative_x;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > diffusivities_derivative_y;
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > > diffusivities_derivative_z;
        
        std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivative_x_computed;
        std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivative_y_computed;
        std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivative_z_computed;
        
        std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivative_xx_computed;
        std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivative_xy_computed;
        std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivative_xz_computed;
        std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivative_yx_computed;
        std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivative_yy_computed;
        std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivative_yz_computed;
        std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivative_zx_computed;
        std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivative_zy_computed;
        std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivative_zz_computed;
        
        /*
         * (1) Compute the derivatives for diffusive flux in x-direction.
         */
        
        // Get the variables for the terms in x-direction in the diffusive flux in x-direction.
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxVariablesForDerivative(
            var_data_x,
            var_component_idx_x,
            DIRECTION::X_DIRECTION,
            DIRECTION::X_DIRECTION);
        
        // Get the diffusivities for the terms in x-direction in the diffusive flux in x-direction.
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxDiffusivities(
            diffusivities_data_x,
            diffusivities_component_idx_x,
            DIRECTION::X_DIRECTION,
            DIRECTION::X_DIRECTION);
        
        TBOX_ASSERT(static_cast<int>(var_data_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_x.size()) == d_num_eqn);
        
        // Compute the first derivatives of variables in x-direction.
        computeFirstDerivativesInX(
            var_derivative_x,
            derivative_x_computed,
            var_data_x,
            var_component_idx_x,
            patch);
        
        // Compute the first derivatives of diffusivities in x-direction.
        computeFirstDerivativesInX(
            diffusivities_derivative_x,
            derivative_x_computed,
            diffusivities_data_x,
            diffusivities_component_idx_x,
            patch);
        
        // Compute the second derivatives of variables in x-direction.
        computeSecondDerivativesInX(
            var_derivative_xx,
            derivative_xx_computed,
            var_data_x,
            var_component_idx_x,
            patch);
        
        // Add the derivatives to the divergence of diffusive flux.
        addDerivativeToDivergence(
            diffusive_flux_divergence,
            var_derivative_x,
            var_derivative_xx,
            diffusivities_data_x,
            diffusivities_derivative_x,
            var_component_idx_x,
            diffusivities_component_idx_x,
            patch,
            dt);
        
        var_data_x.clear();
        
        diffusivities_data_x.clear();
        
        var_component_idx_x.clear();
        
        diffusivities_component_idx_x.clear();
        
        var_derivative_x.clear();
        diffusivities_derivative_x.clear();
        var_derivative_xx.clear();
        
        // Get the variables for the terms in y-direction in the diffusive flux in x-direction.
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxVariablesForDerivative(
            var_data_y,
            var_component_idx_y,
            DIRECTION::X_DIRECTION,
            DIRECTION::Y_DIRECTION);
        
        // Get the diffusivities for the terms in y-direction in the diffusive flux in x-direction.
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxDiffusivities(
            diffusivities_data_y,
            diffusivities_component_idx_y,
            DIRECTION::X_DIRECTION,
            DIRECTION::Y_DIRECTION);
        
        TBOX_ASSERT(static_cast<int>(var_data_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_y.size()) == d_num_eqn);
        
        // Compute the first derivatives of variables in y-direction.
        computeFirstDerivativesInY(
            var_derivative_y,
            derivative_y_computed,
            var_data_y,
            var_component_idx_y,
            patch);
        
        // Compute the first derivatives of diffusivities in x-direction.
        computeFirstDerivativesInX(
            diffusivities_derivative_x,
            derivative_x_computed,
            diffusivities_data_y,
            diffusivities_component_idx_y,
            patch);
        
        // Compute the mixed derivatives of variables.
        var_derivative_component_idx_y.resize(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            var_derivative_component_idx_y[ei].resize(var_derivative_y[ei].size());
            for (int vi = 0; vi < static_cast<int>(var_derivative_component_idx_y[ei].size()); vi++)
            {
                var_derivative_component_idx_y[ei][vi] = 0;
            }
        }
        
        computeFirstDerivativesInX(
            var_derivative_xy,
            derivative_xy_computed,
            var_derivative_y,
            var_derivative_component_idx_y,
            patch);
        
        // Add the derivatives to the divergence of diffusive flux.
        addDerivativeToDivergence(
            diffusive_flux_divergence,
            var_derivative_y,
            var_derivative_xy,
            diffusivities_data_y,
            diffusivities_derivative_x,
            var_component_idx_y,
            diffusivities_component_idx_y,
            patch,
            dt);
        
        var_data_y.clear();
        
        diffusivities_data_y.clear();
        
        var_component_idx_y.clear();
        
        diffusivities_component_idx_y.clear();
        
        var_derivative_component_idx_y.clear();
        
        var_derivative_y.clear();
        diffusivities_derivative_x.clear();
        var_derivative_xy.clear();
        
        // Get the variables for the terms in z-direction in the diffusive flux in x-direction.
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxVariablesForDerivative(
            var_data_z,
            var_component_idx_z,
            DIRECTION::X_DIRECTION,
            DIRECTION::Z_DIRECTION);
        
        // Get the diffusivities for the terms in z-direction in the diffusive flux in x-direction.
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxDiffusivities(
            diffusivities_data_z,
            diffusivities_component_idx_z,
            DIRECTION::X_DIRECTION,
            DIRECTION::Z_DIRECTION);
        
        TBOX_ASSERT(static_cast<int>(var_data_z.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_z.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_z.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_z.size()) == d_num_eqn);
        
        // Compute the first derivatives of variables in z-direction.
        computeFirstDerivativesInZ(
            var_derivative_z,
            derivative_z_computed,
            var_data_z,
            var_component_idx_z,
            patch);
        
        // Compute the first derivatives of diffusivities in x-direction.
        computeFirstDerivativesInX(
            diffusivities_derivative_x,
            derivative_x_computed,
            diffusivities_data_z,
            diffusivities_component_idx_z,
            patch);
        
        // Compute the mixed derivatives of variables.
        var_derivative_component_idx_z.resize(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            var_derivative_component_idx_z[ei].resize(var_derivative_z[ei].size());
            for (int vi = 0; vi < static_cast<int>(var_derivative_component_idx_z[ei].size()); vi++)
            {
                var_derivative_component_idx_z[ei][vi] = 0;
            }
        }
        
        computeFirstDerivativesInX(
            var_derivative_xz,
            derivative_xz_computed,
            var_derivative_z,
            var_derivative_component_idx_z,
            patch);
        
        // Add the derivatives to the divergence of diffusive flux.
        addDerivativeToDivergence(
            diffusive_flux_divergence,
            var_derivative_z,
            var_derivative_xz,
            diffusivities_data_z,
            diffusivities_derivative_x,
            var_component_idx_z,
            diffusivities_component_idx_z,
            patch,
            dt);
        
        var_data_z.clear();
        
        diffusivities_data_z.clear();
        
        var_component_idx_z.clear();
        
        diffusivities_component_idx_z.clear();
        
        var_derivative_component_idx_z.clear();
        
        var_derivative_z.clear();
        diffusivities_derivative_x.clear();
        var_derivative_xz.clear();
        
        /*
         * (2) Compute the derivatives for diffusive flux in y-direction.
         */
        
        // Get the variables for the terms in x-direction in the diffusive flux in y-direction.
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxVariablesForDerivative(
            var_data_x,
            var_component_idx_x,
            DIRECTION::Y_DIRECTION,
            DIRECTION::X_DIRECTION);
        
        // Get the diffusivities for the terms in x-direction in the diffusive flux in y-direction.
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxDiffusivities(
            diffusivities_data_x,
            diffusivities_component_idx_x,
            DIRECTION::Y_DIRECTION,
            DIRECTION::X_DIRECTION);
        
        TBOX_ASSERT(static_cast<int>(var_data_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_x.size()) == d_num_eqn);
        
        // Compute the first derivatives of variables in x-direction.
        computeFirstDerivativesInX(
            var_derivative_x,
            derivative_x_computed,
            var_data_x,
            var_component_idx_x,
            patch);
        
        // Compute the first derivatives of diffusivities in y-direction.
        computeFirstDerivativesInY(
            diffusivities_derivative_y,
            derivative_y_computed,
            diffusivities_data_x,
            diffusivities_component_idx_x,
            patch);
        
        // Compute the mixed derivatives of variables.
        var_derivative_component_idx_x.resize(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            var_derivative_component_idx_x[ei].resize(var_derivative_x[ei].size());
            for (int vi = 0; vi < static_cast<int>(var_derivative_component_idx_x[ei].size()); vi++)
            {
                var_derivative_component_idx_x[ei][vi] = 0;
            }
        }
        
        computeFirstDerivativesInY(
            var_derivative_yx,
            derivative_yx_computed,
            var_derivative_x,
            var_derivative_component_idx_x,
            patch);
        
        // Add the derivatives to the divergence of diffusive flux.
        addDerivativeToDivergence(
            diffusive_flux_divergence,
            var_derivative_x,
            var_derivative_yx,
            diffusivities_data_x,
            diffusivities_derivative_y,
            var_component_idx_x,
            diffusivities_component_idx_x,
            patch,
            dt);
        
        var_data_x.clear();
        
        diffusivities_data_x.clear();
        
        var_component_idx_x.clear();
        
        diffusivities_component_idx_x.clear();
        
        var_derivative_component_idx_x.clear();
        
        var_derivative_x.clear();
        diffusivities_derivative_y.clear();
        var_derivative_yx.clear();
        
        // Get the variables for the terms in y-direction in the diffusive flux in y-direction.
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxVariablesForDerivative(
            var_data_y,
            var_component_idx_y,
            DIRECTION::Y_DIRECTION,
            DIRECTION::Y_DIRECTION);
        
        // Get the diffusivities for the terms in y-direction in the diffusive flux in y-direction.
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxDiffusivities(
            diffusivities_data_y,
            diffusivities_component_idx_y,
            DIRECTION::Y_DIRECTION,
            DIRECTION::Y_DIRECTION);
        
        TBOX_ASSERT(static_cast<int>(var_data_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_y.size()) == d_num_eqn);
        
        // Compute the first derivatives of variables in y-direction.
        computeFirstDerivativesInY(
            var_derivative_y,
            derivative_y_computed,
            var_data_y,
            var_component_idx_y,
            patch);
        
        // Compute the first derivatives of diffusivities in y-direction.
        computeFirstDerivativesInY(
            diffusivities_derivative_y,
            derivative_y_computed,
            diffusivities_data_y,
            diffusivities_component_idx_y,
            patch);
        
        // Compute the second derivatives of variables in y-direction.
        computeSecondDerivativesInY(
            var_derivative_yy,
            derivative_yy_computed,
            var_data_y,
            var_component_idx_y,
            patch);
        
        // Add the derivatives to the divergence of diffusive flux.
        addDerivativeToDivergence(
            diffusive_flux_divergence,
            var_derivative_y,
            var_derivative_yy,
            diffusivities_data_y,
            diffusivities_derivative_y,
            var_component_idx_y,
            diffusivities_component_idx_y,
            patch,
            dt);
        
        var_data_y.clear();
        
        diffusivities_data_y.clear();
        
        var_component_idx_y.clear();
        
        diffusivities_component_idx_y.clear();
        
        var_derivative_y.clear();
        diffusivities_derivative_y.clear();
        var_derivative_yy.clear();
        
        // Get the variables for the terms in z-direction in the diffusive flux in y-direction.
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxVariablesForDerivative(
            var_data_z,
            var_component_idx_z,
            DIRECTION::Y_DIRECTION,
            DIRECTION::Z_DIRECTION);
        
        // Get the diffusivities for the terms in z-direction in the diffusive flux in y-direction.
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxDiffusivities(
            diffusivities_data_z,
            diffusivities_component_idx_z,
            DIRECTION::Y_DIRECTION,
            DIRECTION::Z_DIRECTION);
        
        TBOX_ASSERT(static_cast<int>(var_data_z.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_z.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_z.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_z.size()) == d_num_eqn);
        
        // Compute the first derivatives of variables in z-direction.
        computeFirstDerivativesInZ(
            var_derivative_z,
            derivative_z_computed,
            var_data_z,
            var_component_idx_z,
            patch);
        
        // Compute the first derivatives of diffusivities in y-direction.
        computeFirstDerivativesInY(
            diffusivities_derivative_y,
            derivative_y_computed,
            diffusivities_data_z,
            diffusivities_component_idx_z,
            patch);
        
        // Compute the mixed derivatives of variables.
        var_derivative_component_idx_z.resize(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            var_derivative_component_idx_z[ei].resize(var_derivative_z[ei].size());
            for (int vi = 0; vi < static_cast<int>(var_derivative_component_idx_z[ei].size()); vi++)
            {
                var_derivative_component_idx_z[ei][vi] = 0;
            }
        }
        
        computeFirstDerivativesInY(
            var_derivative_yz,
            derivative_yz_computed,
            var_derivative_z,
            var_derivative_component_idx_z,
            patch);
        
        // Add the derivatives to the divergence of diffusive flux.
        addDerivativeToDivergence(
            diffusive_flux_divergence,
            var_derivative_z,
            var_derivative_yz,
            diffusivities_data_z,
            diffusivities_derivative_y,
            var_component_idx_z,
            diffusivities_component_idx_z,
            patch,
            dt);
        
        var_data_z.clear();
        
        diffusivities_data_z.clear();
        
        var_component_idx_z.clear();
        
        diffusivities_component_idx_z.clear();
        
        var_derivative_component_idx_z.clear();
        
        var_derivative_z.clear();
        diffusivities_derivative_y.clear();
        var_derivative_yz.clear();
        
        /*
         * (3) Compute the derivatives for diffusive flux in z-direction.
         */
        
        // Get the variables for the terms in x-direction in the diffusive flux in z-direction.
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxVariablesForDerivative(
            var_data_x,
            var_component_idx_x,
            DIRECTION::Z_DIRECTION,
            DIRECTION::X_DIRECTION);
        
        // Get the diffusivities for the terms in x-direction in the diffusive flux in z-direction.
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxDiffusivities(
            diffusivities_data_x,
            diffusivities_component_idx_x,
            DIRECTION::Z_DIRECTION,
            DIRECTION::X_DIRECTION);
        
        TBOX_ASSERT(static_cast<int>(var_data_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_x.size()) == d_num_eqn);
        
        // Compute the first derivatives of variables in x-direction.
        computeFirstDerivativesInX(
            var_derivative_x,
            derivative_x_computed,
            var_data_x,
            var_component_idx_x,
            patch);
        
        // Compute the first derivatives of diffusivities in z-direction.
        computeFirstDerivativesInZ(
            diffusivities_derivative_z,
            derivative_z_computed,
            diffusivities_data_x,
            diffusivities_component_idx_x,
            patch);
        
        // Compute the mixed derivatives of variables.
        var_derivative_component_idx_x.resize(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            var_derivative_component_idx_x[ei].resize(var_derivative_x[ei].size());
            for (int vi = 0; vi < static_cast<int>(var_derivative_component_idx_x[ei].size()); vi++)
            {
                var_derivative_component_idx_x[ei][vi] = 0;
            }
        }
        
        computeFirstDerivativesInZ(
            var_derivative_zx,
            derivative_zx_computed,
            var_derivative_x,
            var_derivative_component_idx_x,
            patch);
        
        // Add the derivatives to the divergence of diffusive flux.
        addDerivativeToDivergence(
            diffusive_flux_divergence,
            var_derivative_x,
            var_derivative_zx,
            diffusivities_data_x,
            diffusivities_derivative_z,
            var_component_idx_x,
            diffusivities_component_idx_x,
            patch,
            dt);
        
        var_data_x.clear();
        
        diffusivities_data_x.clear();
        
        var_component_idx_x.clear();
        
        diffusivities_component_idx_x.clear();
        
        var_derivative_component_idx_x.clear();
        
        var_derivative_x.clear();
        diffusivities_derivative_z.clear();
        var_derivative_zx.clear();
        
        // Get the variables for the terms in y-direction in the diffusive flux in z-direction.
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxVariablesForDerivative(
            var_data_y,
            var_component_idx_y,
            DIRECTION::Z_DIRECTION,
            DIRECTION::Y_DIRECTION);
        
        // Get the diffusivities for the terms in y-direction in the diffusive flux in z-direction.
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxDiffusivities(
            diffusivities_data_y,
            diffusivities_component_idx_y,
            DIRECTION::Z_DIRECTION,
            DIRECTION::Y_DIRECTION);
        
        TBOX_ASSERT(static_cast<int>(var_data_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_y.size()) == d_num_eqn);
        
        // Compute the first derivatives of variables in y-direction.
        computeFirstDerivativesInY(
            var_derivative_y,
            derivative_y_computed,
            var_data_y,
            var_component_idx_y,
            patch);
        
        // Compute the first derivatives of diffusivities in z-direction.
        computeFirstDerivativesInZ(
            diffusivities_derivative_z,
            derivative_z_computed,
            diffusivities_data_y,
            diffusivities_component_idx_y,
            patch);
        
        // Compute the mixed derivatives of variables.
        var_derivative_component_idx_y.resize(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            var_derivative_component_idx_y[ei].resize(var_derivative_y[ei].size());
            for (int vi = 0; vi < static_cast<int>(var_derivative_component_idx_y[ei].size()); vi++)
            {
                var_derivative_component_idx_y[ei][vi] = 0;
            }
        }
        
        computeFirstDerivativesInZ(
            var_derivative_zy,
            derivative_zy_computed,
            var_derivative_y,
            var_derivative_component_idx_y,
            patch);
        
        // Add the derivatives to the divergence of diffusive flux.
        addDerivativeToDivergence(
            diffusive_flux_divergence,
            var_derivative_y,
            var_derivative_zy,
            diffusivities_data_y,
            diffusivities_derivative_z,
            var_component_idx_y,
            diffusivities_component_idx_y,
            patch,
            dt);
        
        var_data_y.clear();
        
        diffusivities_data_y.clear();
        
        var_component_idx_y.clear();
        
        diffusivities_component_idx_y.clear();
        
        var_derivative_component_idx_y.clear();
        
        var_derivative_y.clear();
        diffusivities_derivative_z.clear();
        var_derivative_zy.clear();
        
        // Get the variables for the terms in z-direction in the diffusive flux in z-direction.
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxVariablesForDerivative(
            var_data_z,
            var_component_idx_z,
            DIRECTION::Z_DIRECTION,
            DIRECTION::Z_DIRECTION);
        
        // Get the diffusivities for the terms in z-direction in the diffusive flux in z-direction.
        diffusive_flux_utilities->getCellDataOfDiffusiveFluxDiffusivities(
            diffusivities_data_z,
            diffusivities_component_idx_z,
            DIRECTION::Z_DIRECTION,
            DIRECTION::Z_DIRECTION);
        
        TBOX_ASSERT(static_cast<int>(var_data_z.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_z.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_z.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_z.size()) == d_num_eqn);
        
        // Compute the first derivatives of variables in z-direction.
        computeFirstDerivativesInZ(
            var_derivative_z,
            derivative_z_computed,
            var_data_z,
            var_component_idx_z,
            patch);
        
        // Compute the first derivatives of diffusivities in z-direction.
        computeFirstDerivativesInZ(
            diffusivities_derivative_z,
            derivative_z_computed,
            diffusivities_data_z,
            diffusivities_component_idx_z,
            patch);
        
        // Compute the second derivatives of variables in z-direction.
        computeSecondDerivativesInZ(
            var_derivative_zz,
            derivative_zz_computed,
            var_data_z,
            var_component_idx_z,
            patch);
        
        // Add the derivatives to the divergence of diffusive flux.
        addDerivativeToDivergence(
            diffusive_flux_divergence,
            var_derivative_z,
            var_derivative_zz,
            diffusivities_data_z,
            diffusivities_derivative_z,
            var_component_idx_z,
            diffusivities_component_idx_z,
            patch,
            dt);
        
        var_data_z.clear();
        
        diffusivities_data_z.clear();
        
        var_component_idx_z.clear();
        
        diffusivities_component_idx_z.clear();
        
        var_derivative_z.clear();
        diffusivities_derivative_z.clear();
        var_derivative_zz.clear();
        
        /*
         * Unregister the patch and data of all registered derived cell variables in the flow model.
         */
        
        d_flow_model->unregisterPatch();
    }
}


/*
 * Add derivatives to divergence.
 */
void NonconservativeDiffusiveFluxDivergenceOperator::addDerivativeToDivergence(
    HAMERS_SHARED_PTR<pdat::CellData<double> > & divergence,
    const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& var_first_derivative,
    const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& var_derivative_cross_derivative,
    const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& diffusivities_data,
    const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& diffusivities_first_derivative,
    const std::vector<std::vector<int> >& var_component_idx,
    const std::vector<std::vector<int> >& diffusivities_component_idx,
    const hier::Patch& patch,
    const double dt)
{
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // diffusive ghost cells.
    hier::Box diff_ghost_box = interior_box;
    diff_ghost_box.grow(d_num_diff_ghosts);
    const hier::IntVector diff_ghostcell_dims = diff_ghost_box.numberCells();
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the dimensions and number of ghost cells.
         */
        
        const int interior_dim_0 = interior_dims[0];
        
        const int num_diff_ghosts_0 = d_num_diff_ghosts[0];
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ASSERT(static_cast<int>(diffusivities_component_idx[ei].size()) ==
                        static_cast<int>(var_component_idx[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data[ei].size()) ==
                        static_cast<int>(var_component_idx[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(var_first_derivative[ei].size()) ==
                        static_cast<int>(var_component_idx[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(var_derivative_cross_derivative[ei].size()) ==
                        static_cast<int>(var_component_idx[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_first_derivative[ei].size()) ==
                        static_cast<int>(var_component_idx[ei].size()));
            
            double* nabla_F = divergence->getPointer(ei);
            
            for (int vi = 0; vi < static_cast<int>(var_component_idx[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx[ei][vi];
                
                // Get the pointer to diffusivity.
                double* mu = diffusivities_data[ei][vi]->getPointer(mu_idx);
                
                // Get the pointers to derivatives.
                double* dudx = var_first_derivative[ei][vi]->getPointer(0);
                double* d2udxdy = var_derivative_cross_derivative[ei][vi]->getPointer(0);
                double* dmudy = diffusivities_first_derivative[ei][vi]->getPointer(0);
                
                /*
                 * Get the sub-ghost cell width and ghost box dimensions of the diffusivity.
                 */
                
                hier::IntVector num_subghosts_diffusivity =
                    diffusivities_data[ei][vi]->getGhostCellWidth();
                
                hier::IntVector subghostcell_dims_diffusivity =
                    diffusivities_data[ei][vi]->getGhostBox().numberCells();
                
                const int num_subghosts_0_diffusivity = num_subghosts_diffusivity[0];
                
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_nghost = i;
                    
                    const int idx_diffusivity = i + num_subghosts_0_diffusivity;
                    
                    const int idx_diff = i + num_diff_ghosts_0;
                    
                    nabla_F[idx_nghost] += dt*(dmudy[idx_diff]*dudx[idx_diff] +
                        mu[idx_diffusivity]*d2udxdy[idx_diff]);
                }
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the dimensions and number of ghost cells.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        const int num_diff_ghosts_0 = d_num_diff_ghosts[0];
        const int num_diff_ghosts_1 = d_num_diff_ghosts[1];
        
        const int diff_ghostcell_dim_0 = diff_ghostcell_dims[0];
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ASSERT(static_cast<int>(diffusivities_component_idx[ei].size()) ==
                        static_cast<int>(var_component_idx[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data[ei].size()) ==
                        static_cast<int>(var_component_idx[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(var_first_derivative[ei].size()) ==
                        static_cast<int>(var_component_idx[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(var_derivative_cross_derivative[ei].size()) ==
                        static_cast<int>(var_component_idx[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_first_derivative[ei].size()) ==
                        static_cast<int>(var_component_idx[ei].size()));
            
            double* nabla_F = divergence->getPointer(ei);
            
            for (int vi = 0; vi < static_cast<int>(var_component_idx[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx[ei][vi];
                
                // Get the pointer to diffusivity.
                double* mu = diffusivities_data[ei][vi]->getPointer(mu_idx);
                
                // Get the pointers to derivatives.
                double* dudx = var_first_derivative[ei][vi]->getPointer(0);
                double* d2udxdy = var_derivative_cross_derivative[ei][vi]->getPointer(0);
                double* dmudy = diffusivities_first_derivative[ei][vi]->getPointer(0);
                
                /*
                 * Get the sub-ghost cell width and ghost box dimensions of the diffusivity.
                 */
                
                hier::IntVector num_subghosts_diffusivity =
                    diffusivities_data[ei][vi]->getGhostCellWidth();
                
                hier::IntVector subghostcell_dims_diffusivity =
                    diffusivities_data[ei][vi]->getGhostBox().numberCells();
                
                const int num_subghosts_0_diffusivity = num_subghosts_diffusivity[0];
                const int num_subghosts_1_diffusivity = num_subghosts_diffusivity[1];
                const int subghostcell_dim_0_diffusivity = subghostcell_dims_diffusivity[0];
                
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_nghost = i + j*interior_dim_0;
                        
                        const int idx_diffusivity = (i + num_subghosts_0_diffusivity) +
                            (j + num_subghosts_1_diffusivity)*subghostcell_dim_0_diffusivity;
                        
                        const int idx_diff = (i + num_diff_ghosts_0) +
                            (j + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                        
                        nabla_F[idx_nghost] += dt*(dmudy[idx_diff]*dudx[idx_diff] +
                            mu[idx_diffusivity]*d2udxdy[idx_diff]);
                    }
                }
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the dimensions and number of ghost cells.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        const int interior_dim_2 = interior_dims[2];
        
        const int num_diff_ghosts_0 = d_num_diff_ghosts[0];
        const int num_diff_ghosts_1 = d_num_diff_ghosts[1];
        const int num_diff_ghosts_2 = d_num_diff_ghosts[2];
        
        const int diff_ghostcell_dim_0 = diff_ghostcell_dims[0];
        const int diff_ghostcell_dim_1 = diff_ghostcell_dims[1];
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ASSERT(static_cast<int>(diffusivities_component_idx[ei].size()) ==
                        static_cast<int>(var_component_idx[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data[ei].size()) ==
                        static_cast<int>(var_component_idx[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(var_first_derivative[ei].size()) ==
                        static_cast<int>(var_component_idx[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(var_derivative_cross_derivative[ei].size()) ==
                        static_cast<int>(var_component_idx[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_first_derivative[ei].size()) ==
                        static_cast<int>(var_component_idx[ei].size()));
            
            double* nabla_F = divergence->getPointer(ei);
            
            for (int vi = 0; vi < static_cast<int>(var_component_idx[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx[ei][vi];
                
                // Get the pointer to diffusivity.
                double* mu = diffusivities_data[ei][vi]->getPointer(mu_idx);
                
                // Get the pointers to derivatives.
                double* dudx = var_first_derivative[ei][vi]->getPointer(0);
                double* d2udxdy = var_derivative_cross_derivative[ei][vi]->getPointer(0);
                double* dmudy = diffusivities_first_derivative[ei][vi]->getPointer(0);
                
                /*
                 * Get the sub-ghost cell width and ghost box dimensions of the diffusivity.
                 */
                
                hier::IntVector num_subghosts_diffusivity =
                    diffusivities_data[ei][vi]->getGhostCellWidth();
                
                hier::IntVector subghostcell_dims_diffusivity =
                    diffusivities_data[ei][vi]->getGhostBox().numberCells();
                
                const int num_subghosts_0_diffusivity = num_subghosts_diffusivity[0];
                const int num_subghosts_1_diffusivity = num_subghosts_diffusivity[1];
                const int num_subghosts_2_diffusivity = num_subghosts_diffusivity[2];
                const int subghostcell_dim_0_diffusivity = subghostcell_dims_diffusivity[0];
                const int subghostcell_dim_1_diffusivity = subghostcell_dims_diffusivity[1];
                
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_nghost = i +
                                j*interior_dim_0 +
                                k*interior_dim_0*interior_dim_1;
                            
                            const int idx_diffusivity = (i + num_subghosts_0_diffusivity) +
                                (j + num_subghosts_1_diffusivity)*subghostcell_dim_0_diffusivity +
                                (k + num_subghosts_2_diffusivity)*subghostcell_dim_0_diffusivity*
                                    subghostcell_dim_1_diffusivity;
                            
                            const int idx_diff = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            nabla_F[idx_nghost] += dt*(dmudy[idx_diff]*dudx[idx_diff] +
                                mu[idx_diffusivity]*d2udxdy[idx_diff]);
                        }
                    }
                }
            }
        }
    }
}
