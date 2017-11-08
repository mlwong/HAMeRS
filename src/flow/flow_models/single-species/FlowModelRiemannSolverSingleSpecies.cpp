#include "flow/flow_models/single-species/FlowModelRiemannSolverSingleSpecies.hpp"

/*
 * Compute the convective flux from conservative variables.
 */
void
FlowModelRiemannSolverSingleSpecies::computeConvectiveFluxFromConservativeVariables(
    boost::shared_ptr<pdat::SideData<double> > convective_flux,
    const std::vector<boost::shared_ptr<pdat::CellData<double> > >& conservative_variables,
    const hier::Box& domain,
    const RIEMANN_SOLVER::TYPE& riemann_solver_type)
{
}


/*
 * Compute the convective flux from primitive variables.
 */
void
FlowModelRiemannSolverSingleSpecies::computeConvectiveFluxFromPrimitiveVariables(
    boost::shared_ptr<pdat::SideData<double> > convective_flux,
    const std::vector<boost::shared_ptr<pdat::CellData<double> > >& primitive_variables,
    const hier::Box& domain,
    const RIEMANN_SOLVER::TYPE& riemann_solver_type)
{
}
