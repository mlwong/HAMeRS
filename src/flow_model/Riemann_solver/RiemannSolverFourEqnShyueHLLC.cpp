#include "flow_model/Riemann_solver/RiemannSolverFourEqnShyueHLLC.hpp"

/*
 * Compute the flux and velocity at the intercell face from conservative variables.
 */
void
RiemannSolverFourEqnShyueHLLC::computeIntercellFluxAndVelocityFromConservativeVariables(
    std::vector<boost::reference_wrapper<double> >& flux_intercell,
    std::vector<boost::reference_wrapper<double> >& velocity_intercell,
    const std::vector<boost::reference_wrapper<double> >& conservative_variables_minus,
    const std::vector<boost::reference_wrapper<double> >& conservative_variables_plus,
    const DIRECTION& direction)
{
    NULL_USE(flux_intercell);
    NULL_USE(velocity_intercell);
    NULL_USE(conservative_variables_minus);
    NULL_USE(conservative_variables_plus);
}


/*
 * Compute the flux and velocity at the intercell face from primitive variables.
 */
void
RiemannSolverFourEqnShyueHLLC::computeIntercellFluxAndVelocityFromPrimitiveVariables(
    std::vector<boost::reference_wrapper<double> >& flux_intercell,
    std::vector<boost::reference_wrapper<double> >& velocity_intercell,
    const std::vector<boost::reference_wrapper<double> >& primitive_variables_minus,
    const std::vector<boost::reference_wrapper<double> >& primitive_variables_plus,
    const DIRECTION& direction)
{
    NULL_USE(flux_intercell);
    NULL_USE(velocity_intercell);
    NULL_USE(primitive_variables_minus);
    NULL_USE(primitive_variables_plus);
}


/*
 * Compute the flux and velocity in the x-direction at the intercell face
 * from conservative variables.
 */
void
RiemannSolverFourEqnShyueHLLC::computeIntercellFluxAndVelocityInXDirectionFromConservativeVariables(
    std::vector<boost::reference_wrapper<double> >& F_x_intercell,
    std::vector<boost::reference_wrapper<double> >& vel_intercell,
    const std::vector<boost::reference_wrapper<double> >& Q_L,
    const std::vector<boost::reference_wrapper<double> >& Q_R)
{
    NULL_USE(F_x_intercell);
    NULL_USE(vel_intercell);
    NULL_USE(Q_L);
    NULL_USE(Q_R);
}


/*
 * Compute the flux and velocity in the y-direction at the intercell face
 * from conservative variables.
 */
void
RiemannSolverFourEqnShyueHLLC::computeIntercellFluxAndVelocityInYDirectionFromConservativeVariables(
    std::vector<boost::reference_wrapper<double> >& F_y_intercell,
    std::vector<boost::reference_wrapper<double> >& vel_intercell,
    const std::vector<boost::reference_wrapper<double> >& Q_B,
    const std::vector<boost::reference_wrapper<double> >& Q_T)
{
    NULL_USE(F_y_intercell);
    NULL_USE(vel_intercell);
    NULL_USE(Q_B);
    NULL_USE(Q_T);
}


/*
 * Compute the flux and velocity in the z-direction at the intercell face
 * from conservative variables.
 */
void
RiemannSolverFourEqnShyueHLLC::computeIntercellFluxAndVelocityInZDirectionFromConservativeVariables(
    std::vector<boost::reference_wrapper<double> >& F_z_intercell,
    std::vector<boost::reference_wrapper<double> >& vel_intercell,
    const std::vector<boost::reference_wrapper<double> >& Q_B,
    const std::vector<boost::reference_wrapper<double> >& Q_F)
{
    NULL_USE(F_z_intercell);
    NULL_USE(vel_intercell);
    NULL_USE(Q_B);
    NULL_USE(Q_F);
}


/*
 * Compute the flux and velocity in the x-direction at the intercell face
 * from primitive variables.
 */
void
RiemannSolverFourEqnShyueHLLC::computeIntercellFluxAndVelocityInXDirectionFromPrimitiveVariables(
    std::vector<boost::reference_wrapper<double> >& F_x_intercell,
    std::vector<boost::reference_wrapper<double> >& vel_intercell,
    const std::vector<boost::reference_wrapper<double> >& V_L,
    const std::vector<boost::reference_wrapper<double> >& V_R)
{
    NULL_USE(F_x_intercell);
    NULL_USE(vel_intercell);
    NULL_USE(V_L);
    NULL_USE(V_R);
}


/*
 * Compute the flux and velocity in the y-direction at the intercell face
 * from primitive variables.
 */
void
RiemannSolverFourEqnShyueHLLC::computeIntercellFluxAndVelocityInYDirectionFromPrimitiveVariables(
    std::vector<boost::reference_wrapper<double> >& F_y_intercell,
    std::vector<boost::reference_wrapper<double> >& vel_intercell,
    const std::vector<boost::reference_wrapper<double> >& V_B,
    const std::vector<boost::reference_wrapper<double> >& V_T)
{
    NULL_USE(F_y_intercell);
    NULL_USE(vel_intercell);
    NULL_USE(V_B);
    NULL_USE(V_T);
}


/*
 * Compute the flux and velocity in the z-direction at the intercell face
 * from primitive variables.
 */
void
RiemannSolverFourEqnShyueHLLC::computeIntercellFluxAndVelocityInZDirectionFromPrimitiveVariables(
    std::vector<boost::reference_wrapper<double> >& F_z_intercell,
    std::vector<boost::reference_wrapper<double> >& vel_intercell,
    const std::vector<boost::reference_wrapper<double> >& V_B,
    const std::vector<boost::reference_wrapper<double> >& V_F)
{
    NULL_USE(F_z_intercell);
    NULL_USE(vel_intercell);
    NULL_USE(V_B);
    NULL_USE(V_F);
}
