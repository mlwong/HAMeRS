#include "flow/nonconservative_diffusive_flux_divergence_operators/sixth_order/NonconservativeDiffusiveFluxDivergenceOperatorSixthOrder.hpp"

#include "SAMRAI/geom/CartesianPatchGeometry.h"

#include <map>

NonconservativeDiffusiveFluxDivergenceOperatorSixthOrder::NonconservativeDiffusiveFluxDivergenceOperatorSixthOrder(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_eqn,
    const boost::shared_ptr<FlowModel>& flow_model,
    const boost::shared_ptr<tbox::Database>& nonconservative_diffusive_flux_divergence_operator_db):
        NonconservativeDiffusiveFluxDivergenceOperator(
            object_name,
            dim,
            grid_geometry,
            num_eqn,
            flow_model,
            nonconservative_diffusive_flux_divergence_operator_db)
{
    d_num_diff_ghosts = hier::IntVector::getOne(d_dim)*3;
}


/*
 * Print all characteristics of the non-conservative diffusive flux divergence operator class.
 */
void
NonconservativeDiffusiveFluxDivergenceOperatorSixthOrder::printClassData(
    std::ostream& os) const
{
    os << "\nPrint NonconservativeDiffusiveFluxDivergenceOperatorSixthOrder object..."
       << std::endl;
    
    os << std::endl;
    
    os << "NonconservativeDiffusiveFluxDivergenceOperatorSixthOrder: this = "
       << (NonconservativeDiffusiveFluxDivergenceOperatorSixthOrder *)this
       << std::endl;
    os << "d_object_name = "
       << d_object_name
       << std::endl;
}


/*
 * Put the characteristics of the non-conservative diffusive flux divergence operator class into
 * the restart database.
 */
void
NonconservativeDiffusiveFluxDivergenceOperatorSixthOrder::putToRestart(
   const boost::shared_ptr<tbox::Database>& restart_db) const
{
    restart_db->putString("d_nonconservative_diffusive_flux_divergence_operator", "SIXTH_ORDER");
}


/*
 * Compute the non-conservative diffusive flux divergence on a patch.
 */
void
NonconservativeDiffusiveFluxDivergenceOperatorSixthOrder::computeNonconservativeDiffusiveFluxDivergenceOnPatch(
    hier::Patch& patch,
    const boost::shared_ptr<pdat::CellVariable<double> >& variable_diffusive_flux_divergence,
    const boost::shared_ptr<hier::VariableContext>& data_context,
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
    
    boost::shared_ptr<FlowModelDiffusiveFluxUtilities> diffusive_flux_utilities =
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
    boost::shared_ptr<pdat::CellData<double> > diffusive_flux_divergence(
        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
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
        
        d_flow_model->allocateMemoryForDerivedCellData();
        
        d_flow_model->computeDerivedCellData();
        
        /*
         * Delcare containers for computing flux derivatives in different directions.
         */
        
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_data_x;
        
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > diffusivities_data_x;
        
        std::vector<std::vector<int> > var_component_idx_x;
        
        std::vector<std::vector<int> > diffusivities_component_idx_x;
        
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_derivative_x;
        
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_derivative_xx;
        
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > diffusivities_derivative_x;
        
        std::map<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_x_computed;
        
        std::map<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_xx_computed;
        
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
        
        //Compute the first derivatives of variables in x-direction.
        computeFirstDerivativesInX(
            patch,
            var_derivative_x,
            derivative_x_computed,
            var_data_x,
            var_component_idx_x);
        
        // Compute the first derivatives of diffusivities in x-direction.
        computeFirstDerivativesInX(
            patch,
            diffusivities_derivative_x,
            derivative_x_computed,
            diffusivities_data_x,
            diffusivities_component_idx_x);
        
        // Compute the second derivatives of variables in x-direction.
        computeSecondDerivativesInX(
            patch,
            var_derivative_xx,
            derivative_xx_computed,
            var_data_x,
            var_component_idx_x);
        
        // Add the derivatives to the divergence of diffusive flux.
        
        addDerivativeToDivergence(
            patch,
            diffusive_flux_divergence,
            var_derivative_x,
            var_derivative_xx,
            diffusivities_data_x,
            diffusivities_derivative_x,
            var_component_idx_x,
            diffusivities_component_idx_x,
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
        
        d_flow_model->allocateMemoryForDerivedCellData();
        
        d_flow_model->computeDerivedCellData();
        
        /*
         * Delcare containers for computing flux derivatives in different directions.
         */
        
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_data_x;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_data_y;
        
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > diffusivities_data_x;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > diffusivities_data_y;
        
        std::vector<std::vector<int> > var_component_idx_x;
        std::vector<std::vector<int> > var_component_idx_y;
        
        std::vector<std::vector<int> > diffusivities_component_idx_x;
        std::vector<std::vector<int> > diffusivities_component_idx_y;
        
        std::vector<std::vector<int> > var_derivative_component_idx_x;
        std::vector<std::vector<int> > var_derivative_component_idx_y;
        
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_derivative_x;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_derivative_y;
        
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_derivative_xx;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_derivative_xy;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_derivative_yx;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_derivative_yy;
        
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > diffusivities_derivative_x;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > diffusivities_derivative_y;
        
        std::map<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_x_computed;
        std::map<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_y_computed;
        
        std::map<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_xx_computed;
        std::map<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_xy_computed;
        std::map<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_yx_computed;
        std::map<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_yy_computed;
        
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
        
        //Compute the first derivatives of variables in x-direction.
        computeFirstDerivativesInX(
            patch,
            var_derivative_x,
            derivative_x_computed,
            var_data_x,
            var_component_idx_x);
        
        // Compute the first derivatives of diffusivities in x-direction.
        computeFirstDerivativesInX(
            patch,
            diffusivities_derivative_x,
            derivative_x_computed,
            diffusivities_data_x,
            diffusivities_component_idx_x);
        
        // Compute the second derivatives of variables in x-direction.
        computeSecondDerivativesInX(
            patch,
            var_derivative_xx,
            derivative_xx_computed,
            var_data_x,
            var_component_idx_x);
        
        // Add the derivatives to the divergence of diffusive flux.
        
        addDerivativeToDivergence(
            patch,
            diffusive_flux_divergence,
            var_derivative_x,
            var_derivative_xx,
            diffusivities_data_x,
            diffusivities_derivative_x,
            var_component_idx_x,
            diffusivities_component_idx_x,
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
        
        //Compute the first derivatives of variables in y-direction.
        computeFirstDerivativesInY(
            patch,
            var_derivative_y,
            derivative_y_computed,
            var_data_y,
            var_component_idx_y);
        
        // Compute the first derivatives of diffusivities in x-direction.
        computeFirstDerivativesInX(
            patch,
            diffusivities_derivative_x,
            derivative_x_computed,
            diffusivities_data_y,
            diffusivities_component_idx_y);
        
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
            patch,
            var_derivative_xy,
            derivative_xy_computed,
            var_derivative_y,
            var_derivative_component_idx_y);
        
        // Add the derivatives to the divergence of diffusive flux.
        
        addDerivativeToDivergence(
            patch,
            diffusive_flux_divergence,
            var_derivative_y,
            var_derivative_xy,
            diffusivities_data_y,
            diffusivities_derivative_x,
            var_component_idx_y,
            diffusivities_component_idx_y,
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
        
        //Compute the first derivatives of variables in x-direction.
        computeFirstDerivativesInX(
            patch,
            var_derivative_x,
            derivative_x_computed,
            var_data_x,
            var_component_idx_x);
        
        // Compute the first derivatives of diffusivities in y-direction.
        computeFirstDerivativesInY(
            patch,
            diffusivities_derivative_y,
            derivative_y_computed,
            diffusivities_data_x,
            diffusivities_component_idx_x);
        
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
            patch,
            var_derivative_yx,
            derivative_yx_computed,
            var_derivative_x,
            var_derivative_component_idx_x);
        
        // Add the derivatives to the divergence of diffusive flux.
        
        addDerivativeToDivergence(
            patch,
            diffusive_flux_divergence,
            var_derivative_x,
            var_derivative_yx,
            diffusivities_data_x,
            diffusivities_derivative_y,
            var_component_idx_x,
            diffusivities_component_idx_x,
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
        
        //Compute the first derivatives of variables in y-direction.
        computeFirstDerivativesInY(
            patch,
            var_derivative_y,
            derivative_y_computed,
            var_data_y,
            var_component_idx_y);
        
        // Compute the first derivatives of diffusivities in y-direction.
        computeFirstDerivativesInY(
            patch,
            diffusivities_derivative_y,
            derivative_y_computed,
            diffusivities_data_y,
            diffusivities_component_idx_y);
        
        // Compute the second derivatives of variables in y-direction.
        computeSecondDerivativesInY(
            patch,
            var_derivative_yy,
            derivative_yy_computed,
            var_data_y,
            var_component_idx_y);
        
        // Add the derivatives to the divergence of diffusive flux.
        
        addDerivativeToDivergence(
            patch,
            diffusive_flux_divergence,
            var_derivative_y,
            var_derivative_yy,
            diffusivities_data_y,
            diffusivities_derivative_y,
            var_component_idx_y,
            diffusivities_component_idx_y,
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
        
        d_flow_model->allocateMemoryForDerivedCellData();
        
        d_flow_model->computeDerivedCellData();
        
        /*
         * Delcare containers for computing flux derivatives in different directions.
         */
        
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_data_x;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_data_y;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_data_z;
        
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > diffusivities_data_x;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > diffusivities_data_y;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > diffusivities_data_z;
        
        std::vector<std::vector<int> > var_component_idx_x;
        std::vector<std::vector<int> > var_component_idx_y;
        std::vector<std::vector<int> > var_component_idx_z;
        
        std::vector<std::vector<int> > diffusivities_component_idx_x;
        std::vector<std::vector<int> > diffusivities_component_idx_y;
        std::vector<std::vector<int> > diffusivities_component_idx_z;
        
        std::vector<std::vector<int> > var_derivative_component_idx_x;
        std::vector<std::vector<int> > var_derivative_component_idx_y;
        std::vector<std::vector<int> > var_derivative_component_idx_z;
        
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_derivative_x;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_derivative_y;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_derivative_z;
        
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_derivative_xx;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_derivative_xy;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_derivative_xz;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_derivative_yx;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_derivative_yy;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_derivative_yz;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_derivative_zx;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_derivative_zy;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_derivative_zz;
        
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > diffusivities_derivative_x;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > diffusivities_derivative_y;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > diffusivities_derivative_z;
        
        std::map<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_x_computed;
        std::map<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_y_computed;
        std::map<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_z_computed;
        
        std::map<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_xx_computed;
        std::map<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_xy_computed;
        std::map<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_xz_computed;
        std::map<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_yx_computed;
        std::map<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_yy_computed;
        std::map<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_yz_computed;
        std::map<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_zx_computed;
        std::map<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_zy_computed;
        std::map<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_zz_computed;
        
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
        
        //Compute the first derivatives of variables in x-direction.
        computeFirstDerivativesInX(
            patch,
            var_derivative_x,
            derivative_x_computed,
            var_data_x,
            var_component_idx_x);
        
        // Compute the first derivatives of diffusivities in x-direction.
        computeFirstDerivativesInX(
            patch,
            diffusivities_derivative_x,
            derivative_x_computed,
            diffusivities_data_x,
            diffusivities_component_idx_x);
        
        // Compute the second derivatives of variables in x-direction.
        computeSecondDerivativesInX(
            patch,
            var_derivative_xx,
            derivative_xx_computed,
            var_data_x,
            var_component_idx_x);
        
        // Add the derivatives to the divergence of diffusive flux.
        
        addDerivativeToDivergence(
            patch,
            diffusive_flux_divergence,
            var_derivative_x,
            var_derivative_xx,
            diffusivities_data_x,
            diffusivities_derivative_x,
            var_component_idx_x,
            diffusivities_component_idx_x,
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
        
        //Compute the first derivatives of variables in y-direction.
        computeFirstDerivativesInY(
            patch,
            var_derivative_y,
            derivative_y_computed,
            var_data_y,
            var_component_idx_y);
        
        // Compute the first derivatives of diffusivities in x-direction.
        computeFirstDerivativesInX(
            patch,
            diffusivities_derivative_x,
            derivative_x_computed,
            diffusivities_data_y,
            diffusivities_component_idx_y);
        
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
            patch,
            var_derivative_xy,
            derivative_xy_computed,
            var_derivative_y,
            var_derivative_component_idx_y);
        
        // Add the derivatives to the divergence of diffusive flux.
        
        addDerivativeToDivergence(
            patch,
            diffusive_flux_divergence,
            var_derivative_y,
            var_derivative_xy,
            diffusivities_data_y,
            diffusivities_derivative_x,
            var_component_idx_y,
            diffusivities_component_idx_y,
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
        
        //Compute the first derivatives of variables in z-direction.
        computeFirstDerivativesInZ(
            patch,
            var_derivative_z,
            derivative_z_computed,
            var_data_z,
            var_component_idx_z);
        
        // Compute the first derivatives of diffusivities in x-direction.
        computeFirstDerivativesInX(
            patch,
            diffusivities_derivative_x,
            derivative_x_computed,
            diffusivities_data_z,
            diffusivities_component_idx_z);
        
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
            patch,
            var_derivative_xz,
            derivative_xz_computed,
            var_derivative_z,
            var_derivative_component_idx_z);
        
        // Add the derivatives to the divergence of diffusive flux.
        
        addDerivativeToDivergence(
            patch,
            diffusive_flux_divergence,
            var_derivative_z,
            var_derivative_xz,
            diffusivities_data_z,
            diffusivities_derivative_x,
            var_component_idx_z,
            diffusivities_component_idx_z,
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
        
        //Compute the first derivatives of variables in x-direction.
        computeFirstDerivativesInX(
            patch,
            var_derivative_x,
            derivative_x_computed,
            var_data_x,
            var_component_idx_x);
        
        // Compute the first derivatives of diffusivities in y-direction.
        computeFirstDerivativesInY(
            patch,
            diffusivities_derivative_y,
            derivative_y_computed,
            diffusivities_data_x,
            diffusivities_component_idx_x);
        
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
            patch,
            var_derivative_yx,
            derivative_yx_computed,
            var_derivative_x,
            var_derivative_component_idx_x);
        
        // Add the derivatives to the divergence of diffusive flux.
        
        addDerivativeToDivergence(
            patch,
            diffusive_flux_divergence,
            var_derivative_x,
            var_derivative_yx,
            diffusivities_data_x,
            diffusivities_derivative_y,
            var_component_idx_x,
            diffusivities_component_idx_x,
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
        
        //Compute the first derivatives of variables in y-direction.
        computeFirstDerivativesInY(
            patch,
            var_derivative_y,
            derivative_y_computed,
            var_data_y,
            var_component_idx_y);
        
        // Compute the first derivatives of diffusivities in y-direction.
        computeFirstDerivativesInY(
            patch,
            diffusivities_derivative_y,
            derivative_y_computed,
            diffusivities_data_y,
            diffusivities_component_idx_y);
        
        // Compute the second derivatives of variables in y-direction.
        computeSecondDerivativesInY(
            patch,
            var_derivative_yy,
            derivative_yy_computed,
            var_data_y,
            var_component_idx_y);
        
        // Add the derivatives to the divergence of diffusive flux.
        
        addDerivativeToDivergence(
            patch,
            diffusive_flux_divergence,
            var_derivative_y,
            var_derivative_yy,
            diffusivities_data_y,
            diffusivities_derivative_y,
            var_component_idx_y,
            diffusivities_component_idx_y,
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
        
        //Compute the first derivatives of variables in z-direction.
        computeFirstDerivativesInZ(
            patch,
            var_derivative_z,
            derivative_z_computed,
            var_data_z,
            var_component_idx_z);
        
        // Compute the first derivatives of diffusivities in y-direction.
        computeFirstDerivativesInY(
            patch,
            diffusivities_derivative_y,
            derivative_y_computed,
            diffusivities_data_z,
            diffusivities_component_idx_z);
        
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
            patch,
            var_derivative_yz,
            derivative_yz_computed,
            var_derivative_z,
            var_derivative_component_idx_z);
        
        // Add the derivatives to the divergence of diffusive flux.
        
        addDerivativeToDivergence(
            patch,
            diffusive_flux_divergence,
            var_derivative_z,
            var_derivative_yz,
            diffusivities_data_z,
            diffusivities_derivative_y,
            var_component_idx_z,
            diffusivities_component_idx_z,
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
        
        //Compute the first derivatives of variables in x-direction.
        computeFirstDerivativesInX(
            patch,
            var_derivative_x,
            derivative_x_computed,
            var_data_x,
            var_component_idx_x);
        
        // Compute the first derivatives of diffusivities in z-direction.
        computeFirstDerivativesInZ(
            patch,
            diffusivities_derivative_z,
            derivative_z_computed,
            diffusivities_data_x,
            diffusivities_component_idx_x);
        
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
            patch,
            var_derivative_zx,
            derivative_zx_computed,
            var_derivative_x,
            var_derivative_component_idx_x);
        
        // Add the derivatives to the divergence of diffusive flux.
        
        addDerivativeToDivergence(
            patch,
            diffusive_flux_divergence,
            var_derivative_x,
            var_derivative_zx,
            diffusivities_data_x,
            diffusivities_derivative_z,
            var_component_idx_x,
            diffusivities_component_idx_x,
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
        
        //Compute the first derivatives of variables in y-direction.
        computeFirstDerivativesInY(
            patch,
            var_derivative_y,
            derivative_y_computed,
            var_data_y,
            var_component_idx_y);
        
        // Compute the first derivatives of diffusivities in z-direction.
        computeFirstDerivativesInZ(
            patch,
            diffusivities_derivative_z,
            derivative_z_computed,
            diffusivities_data_y,
            diffusivities_component_idx_y);
        
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
            patch,
            var_derivative_zy,
            derivative_zy_computed,
            var_derivative_y,
            var_derivative_component_idx_y);
        
        // Add the derivatives to the divergence of diffusive flux.
        
        addDerivativeToDivergence(
            patch,
            diffusive_flux_divergence,
            var_derivative_y,
            var_derivative_zy,
            diffusivities_data_y,
            diffusivities_derivative_z,
            var_component_idx_y,
            diffusivities_component_idx_y,
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
        
        //Compute the first derivatives of variables in z-direction.
        computeFirstDerivativesInZ(
            patch,
            var_derivative_z,
            derivative_z_computed,
            var_data_z,
            var_component_idx_z);
        
        // Compute the first derivatives of diffusivities in z-direction.
        computeFirstDerivativesInZ(
            patch,
            diffusivities_derivative_z,
            derivative_z_computed,
            diffusivities_data_z,
            diffusivities_component_idx_z);
        
        // Compute the second derivatives of variables in z-direction.
        computeSecondDerivativesInZ(
            patch,
            var_derivative_zz,
            derivative_zz_computed,
            var_data_z,
            var_component_idx_z);
        
        // Add the derivatives to the divergence of diffusive flux.
        
        addDerivativeToDivergence(
            patch,
            diffusive_flux_divergence,
            var_derivative_z,
            var_derivative_zz,
            diffusivities_data_z,
            diffusivities_derivative_z,
            var_component_idx_z,
            diffusivities_component_idx_z,
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
void
NonconservativeDiffusiveFluxDivergenceOperatorSixthOrder::addDerivativeToDivergence(
    hier::Patch& patch,
    boost::shared_ptr<pdat::CellData<double> > & divergence,
    const std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& var_first_derivative,
    const std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& var_derivative_cross_derivative,
    const std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& diffusivities_data,
    const std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& diffusivities_first_derivative,
    const std::vector<std::vector<int> >& var_component_idx,
    const std::vector<std::vector<int> >& diffusivities_component_idx,
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
                
                // Get the pointer to derivatives.
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
                
                // Get the pointer to derivatives.
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
                
                // Get the pointer to derivatives.
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


/*
 * Compute the first derivatives in the x-direction.
 */
void
NonconservativeDiffusiveFluxDivergenceOperatorSixthOrder::computeFirstDerivativesInX(
    hier::Patch& patch,
    std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& derivative_x,
    std::map<double*, boost::shared_ptr<pdat::CellData<double> > >& derivative_x_computed,
    const std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& data_x,
    const std::vector<std::vector<int> >& data_component_idx_x)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    for (int ei = 0; ei < d_num_eqn; ei ++)
    {
        TBOX_ASSERT(static_cast<int>(data_x[ei].size()) ==
                    static_cast<int>(data_component_idx_x[ei].size()));
    }
#endif
    
    derivative_x.resize(d_num_eqn);
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // diffusive ghost cells.
    hier::Box diff_ghost_box = interior_box;
    diff_ghost_box.grow(d_num_diff_ghosts);
    const hier::IntVector diff_ghostcell_dims = diff_ghost_box.numberCells();
    
    // Get the grid spacing.
    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
    const double* const dx = patch_geom->getDx();
    
    const double dx_0 = dx[0];
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the dimensions and number of ghost cells.
         */
        
        const int interior_dim_0 = interior_dims[0];
        
        const int num_diff_ghosts_0 = d_num_diff_ghosts[0];
        
        for (int ei = 0; ei < d_num_eqn; ei ++)
        {
            derivative_x[ei].reserve(static_cast<int>(data_x[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(data_x[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int u_idx = data_component_idx_x[ei][vi];
                
                if (derivative_x_computed.find(data_x[ei][vi]->getPointer(u_idx))
                    == derivative_x_computed.end())
                {
                    // Get the pointer to variable for derivative.
                    double* u = data_x[ei][vi]->getPointer(u_idx);
                    
                    // Declare container to store the derivative.
                    boost::shared_ptr<pdat::CellData<double> > derivative(
                        new pdat::CellData<double>(
                            interior_box, 1, d_num_diff_ghosts));
                    
                    // Get the pointer to the derivative.
                    double* dudx = derivative->getPointer(0);
                    
                    /*
                     * Get the sub-ghost cell width of the variable.
                     */
                    
                    hier::IntVector num_subghosts_data =
                        data_x[ei][vi]->getGhostCellWidth();
                    
                    const int num_subghosts_0_data = num_subghosts_data[0];
                    
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx = i + num_diff_ghosts_0;
                        
                        const int idx_data_LLL = i - 3 + num_subghosts_0_data;
                        const int idx_data_LL  = i - 2 + num_subghosts_0_data;
                        const int idx_data_L   = i - 1 + num_subghosts_0_data;
                        const int idx_data_R   = i + 1 + num_subghosts_0_data;
                        const int idx_data_RR  = i + 2 + num_subghosts_0_data;
                        const int idx_data_RRR = i + 3 + num_subghosts_0_data;
                        
                        dudx[idx] = (double(3)/double(4)*(u[idx_data_R] - u[idx_data_L]) +
                                     double(-3)/double(20)*(u[idx_data_RR] - u[idx_data_LL]) +
                                     double(1)/double(60)*(u[idx_data_RRR] - u[idx_data_LLL]))/
                                        dx_0;
                    }
                    
                    std::pair<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_pair(
                        u,
                        derivative);
                    
                    derivative_x_computed.insert(derivative_pair);
                }
                
                derivative_x[ei].push_back(
                    derivative_x_computed.find(data_x[ei][vi]->getPointer(u_idx))->
                        second);
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
        
        for (int ei = 0; ei < d_num_eqn; ei ++)
        {
            derivative_x[ei].reserve(static_cast<int>(data_x[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(data_x[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int u_idx = data_component_idx_x[ei][vi];
                
                if (derivative_x_computed.find(data_x[ei][vi]->getPointer(u_idx))
                    == derivative_x_computed.end())
                {
                    // Get the pointer to variable for derivative.
                    double* u = data_x[ei][vi]->getPointer(u_idx);
                    
                    // Declare container to store the derivative.
                    boost::shared_ptr<pdat::CellData<double> > derivative(
                        new pdat::CellData<double>(
                            interior_box, 1, d_num_diff_ghosts));
                    
                    // Get the pointer to the derivative.
                    double* dudx = derivative->getPointer(0);
                    
                    /*
                     * Get the sub-ghost cell width and ghost box dimensions of the variable.
                     */
                    
                    hier::IntVector num_subghosts_data =
                        data_x[ei][vi]->getGhostCellWidth();
                    
                    hier::IntVector subghostcell_dims_data =
                        data_x[ei][vi]->getGhostBox().numberCells();
                    
                    const int num_subghosts_0_data = num_subghosts_data[0];
                    const int num_subghosts_1_data = num_subghosts_data[1];
                    const int subghostcell_dim_0_data = subghostcell_dims_data[0];
                    
                    for (int j = -3; j < interior_dim_1 + 3; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                            
                            const int idx_data_LLL = (i - 3 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_LL = (i - 2 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_L = (i - 1 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_R = (i + 1 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_RR = (i + 2 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_RRR = (i + 3 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            dudx[idx] = (double(3)/double(4)*(u[idx_data_R] - u[idx_data_L]) +
                                         double(-3)/double(20)*(u[idx_data_RR] - u[idx_data_LL]) +
                                         double(1)/double(60)*(u[idx_data_RRR] - u[idx_data_LLL]))/
                                            dx_0;
                        }
                    }
                    
                    std::pair<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_pair(
                        u,
                        derivative);
                    
                    derivative_x_computed.insert(derivative_pair);
                }
                
                derivative_x[ei].push_back(
                    derivative_x_computed.find(data_x[ei][vi]->getPointer(u_idx))->
                        second);
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
        
        for (int ei = 0; ei < d_num_eqn; ei ++)
        {
            derivative_x[ei].reserve(static_cast<int>(data_x[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(data_x[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int u_idx = data_component_idx_x[ei][vi];
                
                if (derivative_x_computed.find(data_x[ei][vi]->getPointer(u_idx))
                    == derivative_x_computed.end())
                {
                    // Get the pointer to variable for derivative.
                    double* u = data_x[ei][vi]->getPointer(u_idx);
                    
                    // Declare container to store the derivative.
                    boost::shared_ptr<pdat::CellData<double> > derivative(
                        new pdat::CellData<double>(
                            interior_box, 1, d_num_diff_ghosts));
                    
                    // Get the pointer to the derivative.
                    double* dudx = derivative->getPointer(0);
                    
                    /*
                     * Get the sub-ghost cell width and ghost box dimensions of the variable.
                     */
                    
                    hier::IntVector num_subghosts_data =
                        data_x[ei][vi]->getGhostCellWidth();
                    
                    hier::IntVector subghostcell_dims_data =
                        data_x[ei][vi]->getGhostBox().numberCells();
                    
                    const int num_subghosts_0_data = num_subghosts_data[0];
                    const int num_subghosts_1_data = num_subghosts_data[1];
                    const int num_subghosts_2_data = num_subghosts_data[2];
                    const int subghostcell_dim_0_data = subghostcell_dims_data[0];
                    const int subghostcell_dim_1_data = subghostcell_dims_data[1];
                    
                    for (int k = -3; k < interior_dim_2 + 3; k++)
                    {
                        for (int j = -3; j < interior_dim_1 + 3; j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                // Compute the linear indices.
                                const int idx = (i + num_diff_ghosts_0) +
                                    (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                    (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                        diff_ghostcell_dim_1;
                                
                                const int idx_data_LLL = (i - 3 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_LL = (i - 2 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_L = (i - 1 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_R = (i + 1 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_RR = (i + 2 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_RRR = (i + 3 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                dudx[idx] = (double(3)/double(4)*(u[idx_data_R] - u[idx_data_L]) +
                                             double(-3)/double(20)*(u[idx_data_RR] - u[idx_data_LL]) +
                                             double(1)/double(60)*(u[idx_data_RRR] - u[idx_data_LLL]))/
                                                dx_0;
                            }
                        }
                    }
                    
                    std::pair<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_pair(
                        u,
                        derivative);
                    
                    derivative_x_computed.insert(derivative_pair);
                }
                
                derivative_x[ei].push_back(
                    derivative_x_computed.find(data_x[ei][vi]->getPointer(u_idx))->
                        second);
            }
        }
    }
}


/*
 * Compute the first derivatives in the y-direction.
 */
void
NonconservativeDiffusiveFluxDivergenceOperatorSixthOrder::computeFirstDerivativesInY(
    hier::Patch& patch,
    std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& derivative_y,
    std::map<double*, boost::shared_ptr<pdat::CellData<double> > >& derivative_y_computed,
    const std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& data_y,
    const std::vector<std::vector<int> >& data_component_idx_y)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    for (int ei = 0; ei < d_num_eqn; ei ++)
    {
        TBOX_ASSERT(static_cast<int>(data_y[ei].size()) ==
                    static_cast<int>(data_component_idx_y[ei].size()));
    }
#endif
    
    derivative_y.resize(d_num_eqn);
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // diffusive ghost cells.
    hier::Box diff_ghost_box = interior_box;
    diff_ghost_box.grow(d_num_diff_ghosts);
    const hier::IntVector diff_ghostcell_dims = diff_ghost_box.numberCells();
    
    // Get the grid spacing.
    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
    const double* const dx = patch_geom->getDx();
    
    const double dx_1 = dx[1];
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": NonconservativeDiffusiveFluxDivergenceOperatorSixthOrder::"
            << "computeFirstDerivativesInY()\n"
            << "There isn't y-direction for 1D problem."
            << std::endl);
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
        
        for (int ei = 0; ei < d_num_eqn; ei ++)
        {
            derivative_y[ei].reserve(static_cast<int>(data_y[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(data_y[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int u_idx = data_component_idx_y[ei][vi];
                
                if (derivative_y_computed.find(data_y[ei][vi]->getPointer(u_idx))
                    == derivative_y_computed.end())
                {
                    // Get the pointer to variable for derivative.
                    double* u = data_y[ei][vi]->getPointer(u_idx);
                    
                    // Declare container to store the derivative.
                    boost::shared_ptr<pdat::CellData<double> > derivative(
                        new pdat::CellData<double>(
                            interior_box, 1, d_num_diff_ghosts));
                    
                    // Get the pointer to the derivative.
                    double* dudy = derivative->getPointer(0);
                    
                    /*
                     * Get the sub-ghost cell width and ghost box dimensions of the variable.
                     */
                    
                    hier::IntVector num_subghosts_data =
                        data_y[ei][vi]->getGhostCellWidth();
                    
                    hier::IntVector subghostcell_dims_data =
                        data_y[ei][vi]->getGhostBox().numberCells();
                    
                    const int num_subghosts_0_data = num_subghosts_data[0];
                    const int num_subghosts_1_data = num_subghosts_data[1];
                    const int subghostcell_dim_0_data = subghostcell_dims_data[0];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = -3; i < interior_dim_0 + 3; i++)
                        {
                            // Compute the linear indices.
                            const int idx = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                            
                            const int idx_data_BBB = (i + num_subghosts_0_data) +
                                (j - 3 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_BB = (i + num_subghosts_0_data) +
                                (j - 2 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_B = (i + num_subghosts_0_data) +
                                (j - 1 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_T = (i + num_subghosts_0_data) +
                                (j + 1 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_TT = (i + num_subghosts_0_data) +
                                (j + 2 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_TTT = (i + num_subghosts_0_data) +
                                (j + 3 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            dudy[idx] = (double(3)/double(4)*(u[idx_data_T] - u[idx_data_B]) +
                                         double(-3)/double(20)*(u[idx_data_TT] - u[idx_data_BB]) +
                                         double(1)/double(60)*(u[idx_data_TTT] - u[idx_data_BBB]))/
                                            dx_1;
                        }
                    }
                    
                    std::pair<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_pair(
                        u,
                        derivative);
                    
                    derivative_y_computed.insert(derivative_pair);
                }
                
                derivative_y[ei].push_back(
                    derivative_y_computed.find(data_y[ei][vi]->getPointer(u_idx))->
                        second);
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
        
        for (int ei = 0; ei < d_num_eqn; ei ++)
        {
            derivative_y[ei].reserve(static_cast<int>(data_y[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(data_y[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int u_idx = data_component_idx_y[ei][vi];
                
                if (derivative_y_computed.find(data_y[ei][vi]->getPointer(u_idx))
                    == derivative_y_computed.end())
                {
                    // Get the pointer to variable for derivative.
                    double* u = data_y[ei][vi]->getPointer(u_idx);
                    
                    // Declare container to store the derivative.
                    boost::shared_ptr<pdat::CellData<double> > derivative(
                        new pdat::CellData<double>(
                            interior_box, 1, d_num_diff_ghosts));
                    
                    // Get the pointer to the derivative.
                    double* dudy = derivative->getPointer(0);
                    
                    /*
                     * Get the sub-ghost cell width and ghost box dimensions of the variable.
                     */
                    
                    hier::IntVector num_subghosts_data =
                        data_y[ei][vi]->getGhostCellWidth();
                    
                    hier::IntVector subghostcell_dims_data =
                        data_y[ei][vi]->getGhostBox().numberCells();
                    
                    const int num_subghosts_0_data = num_subghosts_data[0];
                    const int num_subghosts_1_data = num_subghosts_data[1];
                    const int num_subghosts_2_data = num_subghosts_data[2];
                    const int subghostcell_dim_0_data = subghostcell_dims_data[0];
                    const int subghostcell_dim_1_data = subghostcell_dims_data[1];
                    
                    for (int k = -3; k < interior_dim_2 + 3; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = -3; i < interior_dim_0 + 3; i++)
                            {
                                // Compute the linear indices.
                                const int idx = (i + num_diff_ghosts_0) +
                                    (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                    (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                        diff_ghostcell_dim_1;
                                
                                const int idx_data_BBB = (i + num_subghosts_0_data) +
                                    (j - 3 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_BB = (i + num_subghosts_0_data) +
                                    (j - 2 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_B = (i + num_subghosts_0_data) +
                                    (j - 1 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_T = (i + num_subghosts_0_data) +
                                    (j + 1 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_TT = (i + num_subghosts_0_data) +
                                    (j + 2 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_TTT = (i + num_subghosts_0_data) +
                                    (j + 3 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                dudy[idx] = (double(3)/double(4)*(u[idx_data_T] - u[idx_data_B]) +
                                             double(-3)/double(20)*(u[idx_data_TT] - u[idx_data_BB]) +
                                             double(1)/double(60)*(u[idx_data_TTT] - u[idx_data_BBB]))/
                                                dx_1;
                            }
                        }
                    }
                    
                    std::pair<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_pair(
                        u,
                        derivative);
                    
                    derivative_y_computed.insert(derivative_pair);
                }
                
                derivative_y[ei].push_back(
                    derivative_y_computed.find(data_y[ei][vi]->getPointer(u_idx))->
                        second);
            }
        }
    }
}


/*
 * Compute the first derivatives in the z-direction.
 */
void
NonconservativeDiffusiveFluxDivergenceOperatorSixthOrder::computeFirstDerivativesInZ(
    hier::Patch& patch,
    std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& derivative_z,
    std::map<double*, boost::shared_ptr<pdat::CellData<double> > >& derivative_z_computed,
    const std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& data_z,
    const std::vector<std::vector<int> >& data_component_idx_z)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    for (int ei = 0; ei < d_num_eqn; ei ++)
    {
        TBOX_ASSERT(static_cast<int>(data_z[ei].size()) ==
                    static_cast<int>(data_component_idx_z[ei].size()));
    }
#endif
    
    derivative_z.resize(d_num_eqn);
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // diffusive ghost cells.
    hier::Box diff_ghost_box = interior_box;
    diff_ghost_box.grow(d_num_diff_ghosts);
    const hier::IntVector diff_ghostcell_dims = diff_ghost_box.numberCells();
    
    // Get the grid spacing.
    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
    const double* const dx = patch_geom->getDx();
    
    const double dx_2 = dx[2];
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": NonconservativeDiffusiveFluxDivergenceOperatorSixthOrder::"
            << "computeFirstDerivativesInZ()\n"
            << "There isn't z-direction for 1D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": NonconservativeDiffusiveFluxDivergenceOperatorSixthOrder::"
            << "computeFirstDerivativesInZ()\n"
            << "There isn't z-direction for 2D problem."
            << std::endl);
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
        
        for (int ei = 0; ei < d_num_eqn; ei ++)
        {
            derivative_z[ei].reserve(static_cast<int>(data_z[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(data_z[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int u_idx = data_component_idx_z[ei][vi];
                
                if (derivative_z_computed.find(data_z[ei][vi]->getPointer(u_idx))
                    == derivative_z_computed.end())
                {
                    // Get the pointer to variable for derivative.
                    double* u = data_z[ei][vi]->getPointer(u_idx);
                    
                    // Declare container to store the derivative.
                    boost::shared_ptr<pdat::CellData<double> > derivative(
                        new pdat::CellData<double>(
                            interior_box, 1, d_num_diff_ghosts));
                    
                    // Get the pointer to the derivative.
                    double* dudz = derivative->getPointer(0);
                    
                    /*
                     * Get the sub-ghost cell width and ghost box dimensions of the variable.
                     */
                    
                    hier::IntVector num_subghosts_data =
                        data_z[ei][vi]->getGhostCellWidth();
                    
                    hier::IntVector subghostcell_dims_data =
                        data_z[ei][vi]->getGhostBox().numberCells();
                    
                    const int num_subghosts_0_data = num_subghosts_data[0];
                    const int num_subghosts_1_data = num_subghosts_data[1];
                    const int num_subghosts_2_data = num_subghosts_data[2];
                    const int subghostcell_dim_0_data = subghostcell_dims_data[0];
                    const int subghostcell_dim_1_data = subghostcell_dims_data[1];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = -3; j < interior_dim_1 + 3; j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = -3; i < interior_dim_0 + 3; i++)
                            {
                                // Compute the linear indices.
                                const int idx = (i + num_diff_ghosts_0) +
                                    (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                    (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                        diff_ghostcell_dim_1;
                                
                                const int idx_data_BBB = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k - 3 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_BB = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k - 2 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_B = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k - 1 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_F = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + 1 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_FF = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + 2 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_FFF = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + 3 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                dudz[idx] = (double(3)/double(4)*(u[idx_data_F] - u[idx_data_B]) +
                                             double(-3)/double(20)*(u[idx_data_FF] - u[idx_data_BB]) +
                                             double(1)/double(60)*(u[idx_data_FFF] - u[idx_data_BBB]))/
                                                dx_2;
                            }
                        }
                    }
                    
                    std::pair<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_pair(
                        u,
                        derivative);
                    
                    derivative_z_computed.insert(derivative_pair);
                }
                
                derivative_z[ei].push_back(
                    derivative_z_computed.find(data_z[ei][vi]->getPointer(u_idx))->
                        second);
            }
        }
    }
}


/*
 * Compute the second derivatives in the x-direction.
 */
void
NonconservativeDiffusiveFluxDivergenceOperatorSixthOrder::computeSecondDerivativesInX(
    hier::Patch& patch,
    std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& derivative_x,
    std::map<double*, boost::shared_ptr<pdat::CellData<double> > >& derivative_x_computed,
    const std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& data_x,
    const std::vector<std::vector<int> >& data_component_idx_x)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    for (int ei = 0; ei < d_num_eqn; ei ++)
    {
        TBOX_ASSERT(static_cast<int>(data_x[ei].size()) ==
                    static_cast<int>(data_component_idx_x[ei].size()));
    }
#endif
    
    derivative_x.resize(d_num_eqn);
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // diffusive ghost cells.
    hier::Box diff_ghost_box = interior_box;
    diff_ghost_box.grow(d_num_diff_ghosts);
    const hier::IntVector diff_ghostcell_dims = diff_ghost_box.numberCells();
    
    // Get the grid spacing.
    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
    const double* const dx = patch_geom->getDx();
    
    const double dx_sq = dx[0]*dx[0];
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the dimensions and number of ghost cells.
         */
        
        const int interior_dim_0 = interior_dims[0];
        
        const int num_diff_ghosts_0 = d_num_diff_ghosts[0];
        
        for (int ei = 0; ei < d_num_eqn; ei ++)
        {
            derivative_x[ei].reserve(static_cast<int>(data_x[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(data_x[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int u_idx = data_component_idx_x[ei][vi];
                
                if (derivative_x_computed.find(data_x[ei][vi]->getPointer(u_idx))
                    == derivative_x_computed.end())
                {
                    // Get the pointer to variable for derivative.
                    double* u = data_x[ei][vi]->getPointer(u_idx);
                    
                    // Declare container to store the derivative.
                    boost::shared_ptr<pdat::CellData<double> > derivative(
                        new pdat::CellData<double>(
                            interior_box, 1, d_num_diff_ghosts));
                    
                    // Get the pointer to the derivative.
                    double* d2udx2 = derivative->getPointer(0);
                    
                    /*
                     * Get the sub-ghost cell width of the variable.
                     */
                    
                    hier::IntVector num_subghosts_data =
                        data_x[ei][vi]->getGhostCellWidth();
                    
                    const int num_subghosts_0_data = num_subghosts_data[0];
                    
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx = i + num_diff_ghosts_0;
                        
                        const int idx_data_LLL = i - 3 + num_subghosts_0_data;
                        const int idx_data_LL  = i - 2 + num_subghosts_0_data;
                        const int idx_data_L   = i - 1 + num_subghosts_0_data;
                        const int idx_data     = i     + num_subghosts_0_data;
                        const int idx_data_R   = i + 1 + num_subghosts_0_data;
                        const int idx_data_RR  = i + 2 + num_subghosts_0_data;
                        const int idx_data_RRR = i + 3 + num_subghosts_0_data;
                        
                        d2udx2[idx] = (double(-49)/double(18)*u[idx_data] +
                                       double(3)/double(2)*(u[idx_data_L] + u[idx_data_R]) +
                                       double(-3)/double(20)*(u[idx_data_LL] + u[idx_data_RR]) +
                                       double(1)/double(90)*(u[idx_data_LLL] + u[idx_data_RRR]))/
                                        dx_sq;
                    }
                    
                    std::pair<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_pair(
                        u,
                        derivative);
                    
                    derivative_x_computed.insert(derivative_pair);
                }
                
                derivative_x[ei].push_back(
                    derivative_x_computed.find(data_x[ei][vi]->getPointer(u_idx))->
                        second);
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
        
        for (int ei = 0; ei < d_num_eqn; ei ++)
        {
            derivative_x[ei].reserve(static_cast<int>(data_x[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(data_x[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int u_idx = data_component_idx_x[ei][vi];
                
                if (derivative_x_computed.find(data_x[ei][vi]->getPointer(u_idx))
                    == derivative_x_computed.end())
                {
                    // Get the pointer to variable for derivative.
                    double* u = data_x[ei][vi]->getPointer(u_idx);
                    
                    // Declare container to store the derivative.
                    boost::shared_ptr<pdat::CellData<double> > derivative(
                        new pdat::CellData<double>(
                            interior_box, 1, d_num_diff_ghosts));
                    
                    // Get the pointer to the derivative.
                    double* d2udx2 = derivative->getPointer(0);
                    
                    /*
                     * Get the sub-ghost cell width and ghost box dimensions of the variable.
                     */
                    
                    hier::IntVector num_subghosts_data =
                        data_x[ei][vi]->getGhostCellWidth();
                    
                    hier::IntVector subghostcell_dims_data =
                        data_x[ei][vi]->getGhostBox().numberCells();
                    
                    const int num_subghosts_0_data = num_subghosts_data[0];
                    const int num_subghosts_1_data = num_subghosts_data[1];
                    const int subghostcell_dim_0_data = subghostcell_dims_data[0];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                            
                            const int idx_data_LLL = (i - 3 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_LL = (i - 2 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_L = (i - 1 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data = (i + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_R = (i + 1 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_RR = (i + 2 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_RRR = (i + 3 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            d2udx2[idx] = (double(-49)/double(18)*u[idx_data] +
                                           double(3)/double(2)*(u[idx_data_L] + u[idx_data_R]) +
                                           double(-3)/double(20)*(u[idx_data_LL] + u[idx_data_RR]) +
                                           double(1)/double(90)*(u[idx_data_LLL] + u[idx_data_RRR]))/
                                            dx_sq;
                        }
                    }
                    
                    std::pair<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_pair(
                        u,
                        derivative);
                    
                    derivative_x_computed.insert(derivative_pair);
                }
                
                derivative_x[ei].push_back(
                    derivative_x_computed.find(data_x[ei][vi]->getPointer(u_idx))->
                        second);
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
        
        for (int ei = 0; ei < d_num_eqn; ei ++)
        {
            derivative_x[ei].reserve(static_cast<int>(data_x[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(data_x[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int u_idx = data_component_idx_x[ei][vi];
                
                if (derivative_x_computed.find(data_x[ei][vi]->getPointer(u_idx))
                    == derivative_x_computed.end())
                {
                    // Get the pointer to variable for derivative.
                    double* u = data_x[ei][vi]->getPointer(u_idx);
                    
                    // Declare container to store the derivative.
                    boost::shared_ptr<pdat::CellData<double> > derivative(
                        new pdat::CellData<double>(
                            interior_box, 1, d_num_diff_ghosts));
                    
                    // Get the pointer to the derivative.
                    double* d2udx2 = derivative->getPointer(0);
                    
                    /*
                     * Get the sub-ghost cell width and ghost box dimensions of the variable.
                     */
                    
                    hier::IntVector num_subghosts_data =
                        data_x[ei][vi]->getGhostCellWidth();
                    
                    hier::IntVector subghostcell_dims_data =
                        data_x[ei][vi]->getGhostBox().numberCells();
                    
                    const int num_subghosts_0_data = num_subghosts_data[0];
                    const int num_subghosts_1_data = num_subghosts_data[1];
                    const int num_subghosts_2_data = num_subghosts_data[2];
                    const int subghostcell_dim_0_data = subghostcell_dims_data[0];
                    const int subghostcell_dim_1_data = subghostcell_dims_data[1];
                    
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
                                const int idx = (i + num_diff_ghosts_0) +
                                    (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                    (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                        diff_ghostcell_dim_1;
                                
                                const int idx_data_LLL = (i - 3 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_LL = (i - 2 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_L = (i - 1 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_R = (i + 1 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_RR = (i + 2 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_RRR = (i + 3 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                d2udx2[idx] = (double(-49)/double(18)*u[idx_data] +
                                               double(3)/double(2)*(u[idx_data_L] + u[idx_data_R]) +
                                               double(-3)/double(20)*(u[idx_data_LL] + u[idx_data_RR]) +
                                               double(1)/double(90)*(u[idx_data_LLL] + u[idx_data_RRR]))/
                                                dx_sq;
                            }
                        }
                    }
                    
                    std::pair<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_pair(
                        u,
                        derivative);
                    
                    derivative_x_computed.insert(derivative_pair);
                }
                
                derivative_x[ei].push_back(
                    derivative_x_computed.find(data_x[ei][vi]->getPointer(u_idx))->
                        second);
            }
        }
    }
}


/*
 * Compute the second derivatives in the y-direction.
 */
void
NonconservativeDiffusiveFluxDivergenceOperatorSixthOrder::computeSecondDerivativesInY(
    hier::Patch& patch,
    std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& derivative_y,
    std::map<double*, boost::shared_ptr<pdat::CellData<double> > >& derivative_y_computed,
    const std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& data_y,
    const std::vector<std::vector<int> >& data_component_idx_y)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    for (int ei = 0; ei < d_num_eqn; ei ++)
    {
        TBOX_ASSERT(static_cast<int>(data_y[ei].size()) ==
                    static_cast<int>(data_component_idx_y[ei].size()));
    }
#endif
    
    derivative_y.resize(d_num_eqn);
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // diffusive ghost cells.
    hier::Box diff_ghost_box = interior_box;
    diff_ghost_box.grow(d_num_diff_ghosts);
    const hier::IntVector diff_ghostcell_dims = diff_ghost_box.numberCells();
    
    // Get the grid spacing.
    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
    const double* const dx = patch_geom->getDx();
    
    const double dy_sq = dx[1]*dx[1];
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": NonconservativeDiffusiveFluxDivergenceOperatorSixthOrder::"
            << "computeSecondDerivativesInY()\n"
            << "There isn't y-direction for 1D problem."
            << std::endl);
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
        
        for (int ei = 0; ei < d_num_eqn; ei ++)
        {
            derivative_y[ei].reserve(static_cast<int>(data_y[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(data_y[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int u_idx = data_component_idx_y[ei][vi];
                
                if (derivative_y_computed.find(data_y[ei][vi]->getPointer(u_idx))
                    == derivative_y_computed.end())
                {
                    // Get the pointer to variable for derivative.
                    double* u = data_y[ei][vi]->getPointer(u_idx);
                    
                    // Declare container to store the derivative.
                    boost::shared_ptr<pdat::CellData<double> > derivative(
                        new pdat::CellData<double>(
                            interior_box, 1, d_num_diff_ghosts));
                    
                    // Get the pointer to the derivative.
                    double* d2udy2 = derivative->getPointer(0);
                    
                    /*
                     * Get the sub-ghost cell width and ghost box dimensions of the variable.
                     */
                    
                    hier::IntVector num_subghosts_data =
                        data_y[ei][vi]->getGhostCellWidth();
                    
                    hier::IntVector subghostcell_dims_data =
                        data_y[ei][vi]->getGhostBox().numberCells();
                    
                    const int num_subghosts_0_data = num_subghosts_data[0];
                    const int num_subghosts_1_data = num_subghosts_data[1];
                    const int subghostcell_dim_0_data = subghostcell_dims_data[0];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                            
                            const int idx_data_BBB = (i + num_subghosts_0_data) +
                                (j - 3 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_BB = (i + num_subghosts_0_data) +
                                (j - 2 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_B = (i + num_subghosts_0_data) +
                                (j - 1 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data = (i + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_T = (i + num_subghosts_0_data) +
                                (j + 1 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_TT = (i + num_subghosts_0_data) +
                                (j + 2 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_TTT = (i + num_subghosts_0_data) +
                                (j + 3 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            d2udy2[idx] = (double(-49)/double(18)*u[idx_data] +
                                           double(3)/double(2)*(u[idx_data_B] + u[idx_data_T]) +
                                           double(-3)/double(20)*(u[idx_data_BB] + u[idx_data_TT]) +
                                           double(1)/double(90)*(u[idx_data_BBB] + u[idx_data_TTT]))/
                                            dy_sq;
                        }
                    }
                    
                    std::pair<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_pair(
                        u,
                        derivative);
                    
                    derivative_y_computed.insert(derivative_pair);
                }
                
                derivative_y[ei].push_back(
                    derivative_y_computed.find(data_y[ei][vi]->getPointer(u_idx))->
                        second);
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
        
        for (int ei = 0; ei < d_num_eqn; ei ++)
        {
            derivative_y[ei].reserve(static_cast<int>(data_y[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(data_y[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int u_idx = data_component_idx_y[ei][vi];
                
                if (derivative_y_computed.find(data_y[ei][vi]->getPointer(u_idx))
                    == derivative_y_computed.end())
                {
                    // Get the pointer to variable for derivative.
                    double* u = data_y[ei][vi]->getPointer(u_idx);
                    
                    // Declare container to store the derivative.
                    boost::shared_ptr<pdat::CellData<double> > derivative(
                        new pdat::CellData<double>(
                            interior_box, 1, d_num_diff_ghosts));
                    
                    // Get the pointer to the derivative.
                    double* d2udy2 = derivative->getPointer(0);
                    
                    /*
                     * Get the sub-ghost cell width and ghost box dimensions of the variable.
                     */
                    
                    hier::IntVector num_subghosts_data =
                        data_y[ei][vi]->getGhostCellWidth();
                    
                    hier::IntVector subghostcell_dims_data =
                        data_y[ei][vi]->getGhostBox().numberCells();
                    
                    const int num_subghosts_0_data = num_subghosts_data[0];
                    const int num_subghosts_1_data = num_subghosts_data[1];
                    const int num_subghosts_2_data = num_subghosts_data[2];
                    const int subghostcell_dim_0_data = subghostcell_dims_data[0];
                    const int subghostcell_dim_1_data = subghostcell_dims_data[1];
                    
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
                                const int idx = (i + num_diff_ghosts_0) +
                                    (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                    (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                        diff_ghostcell_dim_1;
                                
                                const int idx_data_BBB = (i + num_subghosts_0_data) +
                                    (j - 3 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_BB = (i + num_subghosts_0_data) +
                                    (j - 2 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_B = (i + num_subghosts_0_data) +
                                    (j - 1 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_T = (i + num_subghosts_0_data) +
                                    (j + 1 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_TT = (i + num_subghosts_0_data) +
                                    (j + 2 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_TTT = (i + num_subghosts_0_data) +
                                    (j + 3 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                d2udy2[idx] = (double(-49)/double(18)*u[idx_data] +
                                               double(3)/double(2)*(u[idx_data_B] + u[idx_data_T]) +
                                               double(-3)/double(20)*(u[idx_data_BB] + u[idx_data_TT]) +
                                               double(1)/double(90)*(u[idx_data_BBB] + u[idx_data_TTT]))/
                                                dy_sq;
                            }
                        }
                    }
                    
                    std::pair<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_pair(
                        u,
                        derivative);
                    
                    derivative_y_computed.insert(derivative_pair);
                }
                
                derivative_y[ei].push_back(
                    derivative_y_computed.find(data_y[ei][vi]->getPointer(u_idx))->
                        second);
            }
        }
    }
}


/*
 * Compute the second derivatives in the z-direction.
 */
void
NonconservativeDiffusiveFluxDivergenceOperatorSixthOrder::computeSecondDerivativesInZ(
    hier::Patch& patch,
    std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& derivative_z,
    std::map<double*, boost::shared_ptr<pdat::CellData<double> > >& derivative_z_computed,
    const std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& data_z,
    const std::vector<std::vector<int> >& data_component_idx_z)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    for (int ei = 0; ei < d_num_eqn; ei ++)
    {
        TBOX_ASSERT(static_cast<int>(data_z[ei].size()) ==
                    static_cast<int>(data_component_idx_z[ei].size()));
    }
#endif
    
    derivative_z.resize(d_num_eqn);
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // diffusive ghost cells.
    hier::Box diff_ghost_box = interior_box;
    diff_ghost_box.grow(d_num_diff_ghosts);
    const hier::IntVector diff_ghostcell_dims = diff_ghost_box.numberCells();
    
    // Get the grid spacing.
    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
    const double* const dx = patch_geom->getDx();
    
    const double dz_sq = dx[2]*dx[2];
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": NonconservativeDiffusiveFluxDivergenceOperatorSixthOrder::"
            << "computeSecondDerivativesInZ()\n"
            << "There isn't z-direction for 1D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": NonconservativeDiffusiveFluxDivergenceOperatorSixthOrder::"
            << "computeSecondDerivativesInZ()\n"
            << "There isn't z-direction for 2D problem."
            << std::endl);
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
        
        for (int ei = 0; ei < d_num_eqn; ei ++)
        {
            derivative_z[ei].reserve(static_cast<int>(data_z[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(data_z[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int u_idx = data_component_idx_z[ei][vi];
                
                if (derivative_z_computed.find(data_z[ei][vi]->getPointer(u_idx))
                    == derivative_z_computed.end())
                {
                    // Get the pointer to variable for derivative.
                    double* u = data_z[ei][vi]->getPointer(u_idx);
                    
                    // Declare container to store the derivative.
                    boost::shared_ptr<pdat::CellData<double> > derivative(
                        new pdat::CellData<double>(
                            interior_box, 1, d_num_diff_ghosts));
                    
                    // Get the pointer to the derivative.
                    double* d2udz2 = derivative->getPointer(0);
                    
                    /*
                     * Get the sub-ghost cell width and ghost box dimensions of the variable.
                     */
                    
                    hier::IntVector num_subghosts_data =
                        data_z[ei][vi]->getGhostCellWidth();
                    
                    hier::IntVector subghostcell_dims_data =
                        data_z[ei][vi]->getGhostBox().numberCells();
                    
                    const int num_subghosts_0_data = num_subghosts_data[0];
                    const int num_subghosts_1_data = num_subghosts_data[1];
                    const int num_subghosts_2_data = num_subghosts_data[2];
                    const int subghostcell_dim_0_data = subghostcell_dims_data[0];
                    const int subghostcell_dim_1_data = subghostcell_dims_data[1];
                    
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
                                const int idx = (i + num_diff_ghosts_0) +
                                    (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                    (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                        diff_ghostcell_dim_1;
                                
                                const int idx_data_BBB = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k - 3 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_BB = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k - 2 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_B = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k - 1 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_F = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + 1 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_FF = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + 2 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_FFF = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + 3 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                d2udz2[idx] = (double(-49)/double(18)*u[idx_data] +
                                               double(3)/double(2)*(u[idx_data_B] + u[idx_data_F]) +
                                               double(-3)/double(20)*(u[idx_data_BB] + u[idx_data_FF]) +
                                               double(1)/double(90)*(u[idx_data_BBB] + u[idx_data_FFF]))/
                                                dz_sq;
                            }
                        }
                    }
                    
                    std::pair<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_pair(
                        u,
                        derivative);
                    
                    derivative_z_computed.insert(derivative_pair);
                }
                
                derivative_z[ei].push_back(
                    derivative_z_computed.find(data_z[ei][vi]->getPointer(u_idx))->
                        second);
            }
        }
    }
}
