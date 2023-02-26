#include "util/MPI_helpers/MPIHelper.hpp"

MPIHelper::MPIHelper(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy):
        d_object_name(object_name),
        d_dim(dim),
        d_grid_geometry(grid_geometry),
        d_patch_hierarchy(patch_hierarchy),
        d_mpi(tbox::SAMRAI_MPI::getSAMRAIWorld()),
        d_ratio_finest_level_to_coarsest_level(dim),
        d_coarsest_level_dims(dim),
        d_dx_coarsest_level_dims(dim.getValue()),
        d_finest_level_dims(dim),
        d_dx_finest_level_dims(dim.getValue())
{
    /*
     * Compute the refinement ratio from the finest level to the coarsest level.
     */
    
    const int num_levels = d_patch_hierarchy->getNumberOfLevels();
    
    d_ratio_finest_level_to_coarsest_level =
        d_patch_hierarchy->getRatioToCoarserLevel(num_levels - 1);
    for (int li = num_levels - 2; li > 0 ; li--)
    {
        d_ratio_finest_level_to_coarsest_level *= d_patch_hierarchy->getRatioToCoarserLevel(li);
    }
    
    /*
     * Compute the number of cells of physical domain refined to the finest level.
     */
    
    const hier::BoxContainer& physical_domain = d_grid_geometry->getPhysicalDomain();
    const hier::Box& physical_domain_box = physical_domain.front();
    const hier::IntVector& physical_domain_dims = physical_domain_box.numberCells();
    
    d_coarsest_level_dims = physical_domain_dims;
    d_finest_level_dims   = physical_domain_dims*d_ratio_finest_level_to_coarsest_level;
    
    /*
     * Compute grid spacing of the finest refined domain.
     */
    
    const double* dx_tmp = d_grid_geometry->getDx();
    
    for (int di = 0; di < d_dim.getValue(); di++)
    {
        d_dx_coarsest_level_dims[di] = Real(dx_tmp[di]);
        d_dx_finest_level_dims[di]   = Real(dx_tmp[di])/Real(d_ratio_finest_level_to_coarsest_level[di]);
    }
}
