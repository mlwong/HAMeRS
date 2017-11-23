#include "flow/flow_models/five-eqn_Allaire/FlowModelRiemannSolverFiveEqnAllaire.hpp"

/*
 * Compute the convective flux from conservative variables.
 */
void
FlowModelRiemannSolverFiveEqnAllaire::computeConvectiveFluxFromConservativeVariables(
    boost::shared_ptr<pdat::SideData<double> > convective_flux,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables_minus,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables_plus,
    const hier::Box& domain,
    const RIEMANN_SOLVER::TYPE& riemann_solver_type)
{
}


/*
 * Compute the convective flux from primitive variables.
 */
void
FlowModelRiemannSolverFiveEqnAllaire::computeConvectiveFluxFromPrimitiveVariables(
    boost::shared_ptr<pdat::SideData<double> > convective_flux,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables_minus,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables_plus,
    const hier::Box& domain,
    const RIEMANN_SOLVER::TYPE& riemann_solver_type)
{
}
