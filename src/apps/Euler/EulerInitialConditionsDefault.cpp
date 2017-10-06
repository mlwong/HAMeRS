#include "apps/Euler/EulerInitialConditions.hpp"

/*
 * Set the data on the patch interior to some initial values.
 */
void
EulerInitialConditions::initializeDataOnPatch(
    hier::Patch& patch,
    const std::vector<boost::shared_ptr<pdat::CellVariable<double> > >& conservative_variables,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double data_time,
    const bool initial_time)
{
    NULL_USE(data_time);
    
    TBOX_ERROR(d_object_name
        << ": "
        << "No initial conditon is implemented yet!\n"
        << std::endl);
}
