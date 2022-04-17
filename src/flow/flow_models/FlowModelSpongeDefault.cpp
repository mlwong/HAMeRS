#include "flow/flow_models/FlowModelSponge.hpp"

/*
 * Add the effect of the sponge to the source terms.
 */
void
FlowModelSponge::computeSpongeSourceTermsOnPatch(
    HAMERS_SHARED_PTR<pdat::CellData<double> >& source,
    const hier::Patch& patch,
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& conservative_variables,
    const double time,
    const double dt,
    const int RK_step_number)
{
    NULL_USE(patch);
    NULL_USE(source);
    NULL_USE(conservative_variables);
    NULL_USE(time);
    NULL_USE(dt);
    NULL_USE(RK_step_number);
    
    TBOX_WARNING(d_object_name
        << ": "
        << "No sponge is implemented yet while sponge is used!\n"
        << std::endl);
}
