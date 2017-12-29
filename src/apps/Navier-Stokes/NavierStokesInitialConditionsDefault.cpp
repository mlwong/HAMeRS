#include "apps/Navier-Stokes/NavierStokesInitialConditions.hpp"

/*
 * Set the data on the patch interior to some initial values.
 */
void
NavierStokesInitialConditions::initializeDataOnPatch(
    hier::Patch& patch,
    const std::vector<boost::shared_ptr<pdat::CellData<double> > >& conservative_variables,
    const double data_time,
    const bool initial_time)
{
    NULL_USE(patch);
    NULL_USE(conservative_variables);
    NULL_USE(data_time);
    NULL_USE(initial_time);
    
    TBOX_ERROR(d_object_name
        << ": "
        << "No initial conditon is implemented yet!\n"
        << std::endl);
}
