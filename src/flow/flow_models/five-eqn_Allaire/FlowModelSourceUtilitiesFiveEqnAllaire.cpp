#include "flow/flow_models/five-eqn_Allaire/FlowModelSourceUtilitiesFiveEqnAllaire.hpp"

void
FlowModelSourceUtilitiesFiveEqnAllaire::computeSourceOnPatch(
    hier::Patch& patch,
    const boost::shared_ptr<pdat::CellVariable<double> >& variable_source,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double time,
    const double dt,
    const int RK_step_number)
{
    NULL_USE(patch);
    NULL_USE(variable_source);
    NULL_USE(data_context);
    NULL_USE(time);
    NULL_USE(dt);
    NULL_USE(RK_step_number);
}
