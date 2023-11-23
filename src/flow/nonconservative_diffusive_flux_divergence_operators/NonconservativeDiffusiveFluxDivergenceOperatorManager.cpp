#include "flow/nonconservative_diffusive_flux_divergence_operators/NonconservativeDiffusiveFluxDivergenceOperatorManager.hpp"

NonconservativeDiffusiveFluxDivergenceOperatorManager::NonconservativeDiffusiveFluxDivergenceOperatorManager(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_eqn,
    const HAMERS_SHARED_PTR<FlowModel>& flow_model,
    const HAMERS_SHARED_PTR<tbox::Database>& nonconservative_diffusive_flux_divergence_operator_db,
    const std::string& nonconservative_diffusive_flux_divergence_operator_str):
        d_object_name(object_name)
{
    
    if (nonconservative_diffusive_flux_divergence_operator_str == "SECOND_ORDER")
    {
        d_nonconservative_diffusive_flux_divergence_operator_type =
            NONCONSERVATIVE_DIFFUSIVE_FLUX_DIVERGENCE_OPERATOR::SECOND_ORDER;
            
        d_noncons_diff_flux_div_op.reset(new NonconservativeDiffusiveFluxDivergenceOperatorSecondOrder(
            "d_nonconservative_diffusive_flux_divergence_operator",
            dim,
            grid_geometry,
            flow_model->getNumberOfEquations(),
            flow_model,
            nonconservative_diffusive_flux_divergence_operator_db));
    }
    else if (nonconservative_diffusive_flux_divergence_operator_str == "SIXTH_ORDER")
    {
        d_nonconservative_diffusive_flux_divergence_operator_type =
            NONCONSERVATIVE_DIFFUSIVE_FLUX_DIVERGENCE_OPERATOR::SIXTH_ORDER;
        
        d_noncons_diff_flux_div_op.reset(new NonconservativeDiffusiveFluxDivergenceOperatorSixthOrder(
            "d_nonconservative_diffusive_flux_divergence_operator",
            dim,
            grid_geometry,
            flow_model->getNumberOfEquations(),
            flow_model,
            nonconservative_diffusive_flux_divergence_operator_db));
    }
    else if (nonconservative_diffusive_flux_divergence_operator_str == "EIGHTH_ORDER")
    {
        d_nonconservative_diffusive_flux_divergence_operator_type =
            NONCONSERVATIVE_DIFFUSIVE_FLUX_DIVERGENCE_OPERATOR::EIGHTH_ORDER;
        
        d_noncons_diff_flux_div_op.reset(new NonconservativeDiffusiveFluxDivergenceOperatorEighthOrder(
            "d_nonconservative_diffusive_flux_divergence_operator",
            dim,
            grid_geometry,
            flow_model->getNumberOfEquations(),
            flow_model,
            nonconservative_diffusive_flux_divergence_operator_db));
    }
    else if (nonconservative_diffusive_flux_divergence_operator_str == "TENTH_ORDER")
    {
        d_nonconservative_diffusive_flux_divergence_operator_type =
            NONCONSERVATIVE_DIFFUSIVE_FLUX_DIVERGENCE_OPERATOR::TENTH_ORDER;
        
        d_noncons_diff_flux_div_op.reset(new NonconservativeDiffusiveFluxDivergenceOperatorTenthOrder(
            "d_nonconservative_diffusive_flux_divergence_operator",
            dim,
            grid_geometry,
            flow_model->getNumberOfEquations(),
            flow_model,
            nonconservative_diffusive_flux_divergence_operator_db));
    }
    else if (nonconservative_diffusive_flux_divergence_operator_str == "TWELFTH_ORDER")
    {
        d_nonconservative_diffusive_flux_divergence_operator_type =
            NONCONSERVATIVE_DIFFUSIVE_FLUX_DIVERGENCE_OPERATOR::TWELFTH_ORDER;
        
        d_noncons_diff_flux_div_op.reset(new NonconservativeDiffusiveFluxDivergenceOperatorTwelfthOrder(
            "d_nonconservative_diffusive_flux_divergence_operator",
            dim,
            grid_geometry,
            flow_model->getNumberOfEquations(),
            flow_model,
            nonconservative_diffusive_flux_divergence_operator_db));
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Unknown nonconservative_diffusive_flux_divergence_operator string = '"
            << nonconservative_diffusive_flux_divergence_operator_str
            << "' found in input."
            << std::endl);
    }
}


/*
 * Print all characteristics of diffusive flux divergence operator manager.
 */
void
NonconservativeDiffusiveFluxDivergenceOperatorManager::printClassData(std::ostream& os) const
{
    os << "\nPrint NonconservativeDiffusiveFluxDivergenceOperatorManager object..."
       << std::endl;
    
    os << std::endl;
    
    os << "NonconservativeDiffusiveFluxDivergenceOperatorManager: this = "
       << (NonconservativeDiffusiveFluxDivergenceOperatorManager *)this
       << std::endl;
    
    os << "d_object_name = "
       << d_object_name
       << std::endl;
    
    os << "d_nonconservative_diffusive_flux_divergence_operator_type = "
       << d_nonconservative_diffusive_flux_divergence_operator_type
       << std::endl;
}
