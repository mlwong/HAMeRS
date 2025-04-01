#include "util/immersed_boundaries/ImmersedBoundaries.hpp"

void
ImmersedBoundaries::setImmersedBoundaryVariablesOnPatch(
    const hier::Patch& patch,
    const double data_time,
    const bool initial_time,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims,
    const HAMERS_SHARED_PTR<pdat::CellData<int> >& data_mask,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_wall_distance,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_surface_normal)
{
    NULL_USE(patch);
    NULL_USE(data_time);
    NULL_USE(initial_time);
    NULL_USE(domain_lo);
    NULL_USE(domain_dims);
    NULL_USE(data_mask);
    NULL_USE(data_wall_distance);
    NULL_USE(data_surface_normal);
}

void
ImmersedBoundaries::generateSurfaceTriangulation(
    std::vector<std::array<Real, 3> >& nodes,
    std::vector<std::array<int, 3> >& connectivities,
    std::vector<int> component_id)
{
    NULL_USE(nodes);
    NULL_USE(connectivities);
    NULL_USE(component_id);
}