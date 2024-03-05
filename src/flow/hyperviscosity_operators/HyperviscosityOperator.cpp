#include "flow/hyperviscosity_operators/HyperviscosityOperator.hpp"

/*
 * Print all characteristics of the non-conservative diffusive flux divergence operator class.
 */
void
HyperviscosityOperator::printClassData(
    std::ostream& os) const
{
    os << "\nPrint HyperviscosityOperator object..."
       << std::endl;
    
    os << std::endl;
    
    os << "HyperviscosityOperator: this = "
       << (HyperviscosityOperator *)this
       << std::endl;
    os << "d_object_name = "
       << d_object_name
       << std::endl;
}


/*
 * Put the characteristics of the hyperviscosity operator into the restart database.
 */
void
HyperviscosityOperator::putToRestart(
   const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const
{
    NULL_USE(restart_db);
}


/*
 * Perform the hyperviscosity operator on a patch.
 */
void
HyperviscosityOperator::performHyperviscosityOperatorOnPatch(
    hier::Patch& patch,
    const HAMERS_SHARED_PTR<hier::CoarseFineBoundary> coarse_fine_bdry,
    const HAMERS_SHARED_PTR<pdat::SideVariable<Real> >& variable_convective_flux,
    const HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_source,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double time,
    const double dt,
    const int RK_step_number)
{
    NULL_USE(patch);
    NULL_USE(coarse_fine_bdry);
    NULL_USE(variable_convective_flux);
    NULL_USE(variable_source);
    NULL_USE(data_context);
    NULL_USE(time);
    NULL_USE(dt);
    NULL_USE(RK_step_number);
}
