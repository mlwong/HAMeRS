#include "flow/flow_models/FlowModelSpecialSourceTerms.hpp"

/*
 * Add the effects of the special source terms.
 */
void
FlowModelSpecialSourceTerms::computeSpecialSourceTermsOnPatch(
    HAMERS_SHARED_PTR<pdat::CellData<Real> >& source,
    const hier::Patch& patch,
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& conservative_variables,
    const std::unordered_map<std::string, Real>& monitoring_statistics_map,
    const double time,
    const double dt,
    const int RK_step_number)
{
    NULL_USE(patch);
    NULL_USE(source);
    NULL_USE(conservative_variables);
    NULL_USE(monitoring_statistics_map);
    NULL_USE(time);
    NULL_USE(dt);
    NULL_USE(RK_step_number);
    
    TBOX_WARNING(d_object_name
        << ": "
        << "No special source terms are implemented yet while special source terms are used!\n"
        << std::endl);
}


void
FlowModelSpecialSourceTerms::putToRestart(const HAMERS_SHARED_PTR<tbox::Database>& restart_source_terms_db)
{
    putToRestartBase(restart_source_terms_db);
}
