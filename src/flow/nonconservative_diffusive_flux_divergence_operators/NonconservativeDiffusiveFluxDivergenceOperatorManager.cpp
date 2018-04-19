#include "flow/nonconservative_diffusive_flux_divergence_operators/NonconservativeDiffusiveFluxDivergenceOperatorManager.hpp"

NonconservativeDiffusiveFluxDivergenceOperatorManager::NonconservativeDiffusiveFluxDivergenceOperatorManager(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_eqn,
    const boost::shared_ptr<FlowModel>& flow_model,
    const boost::shared_ptr<tbox::Database>& nonconservative_diffusive_flux_divergence_operator_db,
    const std::string& nonconservative_diffusive_flux_divergence_operator_str):
        d_object_name(object_name)
{
    if (nonconservative_diffusive_flux_divergence_operator_str == "SIXTH_ORDER")
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
