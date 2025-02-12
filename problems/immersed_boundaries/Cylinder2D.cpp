#include "util/immersed_boundaries/ImmersedBoundaries.hpp"

void
ImmersedBoundaries::setImmersedBoundaryVariablesOnPatch(
    const hier::Patch& patch,
    const double data_time,
    const bool initial_time,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims,
    const HAMERS_SHARED_PTR<pdat::CellData<int> >& data_mask,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_wall_distance,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_surface_normal)
{
    NULL_USE(data_time);
    
    const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
        HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(patch_geom);
#endif
    
    const double* const dx = patch_geom->getDx();
    const double* const patch_xlo = patch_geom->getXLower();
    
    const hier::IntVector num_ghosts = data_mask->getGhostCellWidth();
    const hier::IntVector ghostcell_dims = data_mask->getGhostBox().numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(num_ghosts == data_wall_distance->getGhostCellWidth());
    TBOX_ASSERT(num_ghosts == data_surface_normal->getGhostCellWidth());
#endif
    
    /*
     * Get the pointers to the data.
     */
    int* mask      = data_mask->getPointer(0);
    double* dist   = data_wall_distance->getPointer(0);
    double* norm_0 = data_surface_normal->getPointer(0);
    double* norm_1 = data_surface_normal->getPointer(1);
    
    /*
     * Get the local lower index, numbers of cells in each dimension and numbers of ghost cells.
     */
    
    const int domain_lo_0 = domain_lo[0];
    const int domain_lo_1 = domain_lo[1];
    const int domain_dim_0 = domain_dims[0];
    const int domain_dim_1 = domain_dims[1];
    
    const int num_ghosts_0 = num_ghosts[0];
    const int num_ghosts_1 = num_ghosts[1];
    const int ghostcell_dim_0 = ghostcell_dims[0];
    
    /************************************************
     * Set the immersed boundary variables from here.
     ************************************************/
    
    /*
     * Set the parameters of the cylinder here.
     */
    
    const double half = double(1)/double(2);
    
    /*
     * These will be read from the input file.
     */
    
    Real radius_c = half; //double(20); AFK 
    Real x_c = Real(1); //half;  AFK
    Real y_c = Real(1); //half;  AFK

    if (d_initial_conditions_db != nullptr)
            {
                TBOX_ASSERT(d_initial_conditions_db->keyExists("x_c"));
                TBOX_ASSERT(d_initial_conditions_db->keyExists("y_c"));
                x_c     = d_initial_conditions_db->getReal("x_c");
                y_c     = d_initial_conditions_db->getReal("y_c");
                radius_c = d_initial_conditions_db->getReal("radius");
            }

    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
    {
        HAMERS_PRAGMA_SIMD
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear index.
            const int idx = (i + num_ghosts_0) +
                (j + num_ghosts_1)*ghostcell_dim_0;
            
            // Compute the coordinates.
            double x[2];
            x[0] = patch_xlo[0] + (double(i) + double(1)/double(2))*dx[0]; // x coordinates of the point.
            x[1] = patch_xlo[1] + (double(j) + double(1)/double(2))*dx[1]; // y coordinates of the point.
            
            // Distance from the cylinder center.
            const double radius = sqrt(pow(x[0] - x_c, 2) + pow(x[1] - y_c, 2));
            // Angle between x axis and a line passing through center and current cell.
            const double theta  = atan2(x[1] - y_c, x[0] - x_c);
            
            if (radius < radius_c) // Condition that should be satisfied to be in cylinder.
            {
                double x_p = double(0); // x coordinates on the cylinder where y = x[1].
                double y_p = double(0); // y coordinates on the cylinder where x = x[0].
                
                if (x[0] > x_c)
                {
                    x_p = x_c + sqrt(pow(radius_c, 2) - pow(radius*sin(theta), 2));
                }
                else
                {
                    x_p = x_c - sqrt(pow(radius_c, 2) - pow(radius*sin(theta), 2));
                }
                
                if (x[1] > y_c)
                {
                    y_p = y_c + sqrt(pow(radius_c, 2) - pow(radius*cos(theta), 2));
                }
                else
                {
                    y_p = y_c - sqrt(pow(radius_c, 2) - pow(radius*cos(theta), 2));
                }
                
                if ((fabs(x_p - x[0]) < (double(d_num_immersed_boundary_ghosts[0]))*dx[0]) ||
                    (fabs(y_p - x[1]) < (double(d_num_immersed_boundary_ghosts[1]))*dx[1]))
                {
                    mask[idx]   = int(IB_MASK::IB_GHOST);
                    dist[idx]   = radius_c - radius;
                    norm_0[idx] = (x[0] - x_c)/radius;
                    norm_1[idx] = (x[1] - y_c)/radius;
                }
                else
                {
                    mask[idx]   = int(IB_MASK::BODY);
                    dist[idx]   = double(0);
                    norm_0[idx] = double(0);
                    norm_1[idx] = double(0);
                }
            }
            else
            {
                mask[idx]   = int(IB_MASK::FLUID);
                dist[idx]   = double(0);
                norm_0[idx] = double(0);
                norm_1[idx] = double(0);
            }
        }
    }
}
