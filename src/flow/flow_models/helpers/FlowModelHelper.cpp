#include "flow/flow_models/helpers/FlowModelHelper.hpp"

/*
 * Get number of points in the x-direction of the refined domain.
 */
int
FlowModelHelper::getRefinedDomainNumberOfPointsX(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const
{
    /*
     * Get the refinement ratio from the finest level to the coarest level.
     */
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    hier::IntVector ratioFinestLevelToCoarestLevel =
        patch_hierarchy->getRatioToCoarserLevel(num_levels - 1);
    for (int li = num_levels - 2; li > 0 ; li--)
    {
        ratioFinestLevelToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(li);
    }
    
    /*
     * Get the number of cells of physical domain refined to the finest level.
     */
    
    const hier::BoxContainer& physical_domain = d_grid_geometry->getPhysicalDomain();
    const hier::Box& physical_domain_box = physical_domain.front();
    const hier::IntVector& physical_domain_dims = physical_domain_box.numberCells();
    const hier::IntVector finest_level_dims = physical_domain_dims*ratioFinestLevelToCoarestLevel;
    
    const int finest_level_dim_0 = finest_level_dims[0];
    return finest_level_dim_0;
}


/*
 * Get grid spacing in the x-direction of the refined domain.
 */
double
FlowModelHelper::getRefinedDomainGridSpacingX(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy) const
{
    /*
     * Get the refinement ratio from the finest level to the coarest level.
     */
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    hier::IntVector ratioFinestLevelToCoarestLevel =
        patch_hierarchy->getRatioToCoarserLevel(num_levels - 1);
    for (int li = num_levels - 2; li > 0 ; li--)
    {
        ratioFinestLevelToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(li);
    }
    
    const double* dx = d_grid_geometry->getDx();
    
    return dx[0]/ratioFinestLevelToCoarestLevel[0];
}
