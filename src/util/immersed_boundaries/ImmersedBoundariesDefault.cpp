#include "util/immersed_boundaries/ImmersedBoundaries.hpp"

void
ImmersedBoundaries::setImmersedBoundaryVariablesOnPatch(
    const hier::Patch& patch,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims,
    const HAMERS_SHARED_PTR<pdat::CellData<int> >& data_mask,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_wall_distance,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_surface_normal,
    const double data_time,
    const bool initial_time)
{
    NULL_USE(patch);
    NULL_USE(domain_lo);
    NULL_USE(domain_dims);
    NULL_USE(data_mask);
    NULL_USE(data_wall_distance);
    NULL_USE(data_surface_normal);
    NULL_USE(data_time);
    NULL_USE(initial_time);
}

