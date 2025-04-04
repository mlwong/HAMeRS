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
    double* norm_2 = data_surface_normal->getPointer(2);
    
    /*
     * Get the local lower index, numbers of cells in each dimension and numbers of ghost cells.
     */
    
    const int domain_lo_0 = domain_lo[0];
    const int domain_lo_1 = domain_lo[1];
    const int domain_lo_2 = domain_lo[2];
    const int domain_dim_0 = domain_dims[0];
    const int domain_dim_1 = domain_dims[1];
    const int domain_dim_2 = domain_dims[2];
    
    const int num_ghosts_0 = num_ghosts[0];
    const int num_ghosts_1 = num_ghosts[1];
    const int num_ghosts_2 = num_ghosts[2];
    const int ghostcell_dim_0 = ghostcell_dims[0];
    const int ghostcell_dim_1 = ghostcell_dims[1];
    
    /************************************************
     * Set the immersed boundary variables from here.
     ************************************************/
    
    /*
     * Set the parameters of the 3D cylinder here.
     */
    
    const double half = double(1)/double(2);
    
    /*
     * These will be read from the input file.
     */
    
    Real radius_c = half;
    Real x_c = Real(1);
    Real y_c = Real(1);
    Real z_c = Real(1);
    
    if (d_initial_conditions_db != nullptr)
    {
        TBOX_ASSERT(d_initial_conditions_db->keyExists("x_c"));
        TBOX_ASSERT(d_initial_conditions_db->keyExists("y_c"));
        TBOX_ASSERT(d_initial_conditions_db->keyExists("z_c")); 
        
        x_c = d_initial_conditions_db->getReal("x_c");
        y_c = d_initial_conditions_db->getReal("y_c");
        z_c = d_initial_conditions_db->getReal("z_c");
        
        radius_c = d_initial_conditions_db->getReal("radius");
    }
    
    for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
    {
       for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear index.
                const int idx = (i + num_ghosts_0) +
                    (j + num_ghosts_1)*ghostcell_dim_0 +
                    (k + num_ghosts_2)*ghostcell_dim_0*ghostcell_dim_1;
                
                // Compute the coordinates.
                double x[3];
                x[0] = patch_xlo[0] + (double(i) + double(1)/double(2))*dx[0]; // x coordinates of the point.
                x[1] = patch_xlo[1] + (double(j) + double(1)/double(2))*dx[1]; // y coordinates of the point.
                x[2] = patch_xlo[2] + (double(k) + double(1)/double(2))*dx[2]; // z coordinates of the point.
                
                // Distance from the cylinder center.
                   const double radius = sqrt(pow(x[0] - x_c, 2) + pow(x[1] - y_c, 2));
                // Angle between x axis and a line passing through center and current cell.
                const double theta  = atan2(x[1] - y_c, x[0] - x_c);
                
                // Double check this eqn for cylindrical coord
                const double phi = atan((x[1] - y_c) / x[0] - x_c);
                
                if (radius < radius_c) // Condition that should be satisfied to be in cylinder.
                {
                    double x_p = double(0); // x coordinates on the cylinder where y = x[1], z = x[2].
                    double y_p = double(0); // y coordinates on the cylinder where x = x[0], z = x[2].
                    //double check this
                    double z_p = double(0); // z coordinates on the cylinder where x = x[0], y = x[1].
                    
                    // For checking ghost cell for convective flux.
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
                    //double check this, z_p does not depend on phi in cylindrical coords, x2 + y2 = r2
                    if (x[2] > z_c)
                    {    
                        z_p = x[2];
                    }
                    else
                    {
                        z_p = x[2];
                    }
                    // For checking ghost cell for viscous flux.
                    // 8 diagonal ghost cells for 3D
                    // Check first diagonal ghost cell.
                    // radius_d does not depend on z, only x and y
                    double x_d[3];
                    x_d[0] = patch_xlo[0] + (double(i+1) + double(1)/double(2))*dx[0]; // x coordinates of the point.
                    x_d[1] = patch_xlo[1] + (double(j+1) + double(1)/double(2))*dx[1]; // y coordinates of the point.
                    x_d[2] = patch_xlo[2] + (double(k+1) + double(1)/double(2))*dx[2]; // z coordinates of the point.
                    double radius_d = sqrt(pow(x_d[0] - x_c, 2) + pow(x_d[1] - y_c, 2));
                    bool is_corner_ghost = radius_d > radius_c;    
                    
                    // Check second diagonal ghost cell.
                    x_d[0] = patch_xlo[0] + (double(i+1) + double(1)/double(2))*dx[0]; // x coordinates of the point.
                    x_d[1] = patch_xlo[1] + (double(j+1) + double(1)/double(2))*dx[1]; // y coordinates of the point.
                    x_d[2] = patch_xlo[2] + (double(k-1) + double(1)/double(2))*dx[2]; // z coordinates of the point.
                    radius_d = sqrt(pow(x_d[0] - x_c, 2) + pow(x_d[1] - y_c, 2));
                    is_corner_ghost |= radius_d > radius_c;
                    
                    // Check third diagonal ghost cell.
                    x_d[0] = patch_xlo[0] + (double(i+1) + double(1)/double(2))*dx[0]; // x coordinates of the point.
                    x_d[1] = patch_xlo[1] + (double(j-1) + double(1)/double(2))*dx[1]; // y coordinates of the point.
                    x_d[2] = patch_xlo[2] + (double(k+1) + double(1)/double(2))*dx[2]; // z coordinates of the point.
                    radius_d = sqrt(pow(x_d[0] - x_c, 2) + pow(x_d[1] - y_c, 2));
                    is_corner_ghost |= radius_d > radius_c;
                    
                    // Check fourth diagonal ghost cell.
                    x_d[0] = patch_xlo[0] + (double(i+1) + double(1)/double(2))*dx[0]; // x coordinates of the point.
                    x_d[1] = patch_xlo[1] + (double(j-1) + double(1)/double(2))*dx[1]; // y coordinates of the point.
                    x_d[2] = patch_xlo[2] + (double(k-1) + double(1)/double(2))*dx[2]; // z coordinates of the point.
                    radius_d = sqrt(pow(x_d[0] - x_c, 2) + pow(x_d[1] - y_c, 2));
                    is_corner_ghost |= radius_d > radius_c;
                    
                    // Check fifth diagonal ghost cell.
                    x_d[0] = patch_xlo[0] + (double(i-1) + double(1)/double(2))*dx[0]; // x coordinates of the point.
                    x_d[1] = patch_xlo[1] + (double(j+1) + double(1)/double(2))*dx[1]; // y coordinates of the point.
                    x_d[2] = patch_xlo[2] + (double(k+1) + double(1)/double(2))*dx[2]; // z coordinates of the point.
                    radius_d = sqrt(pow(x_d[0] - x_c, 2) + pow(x_d[1] - y_c, 2));
                    is_corner_ghost |= radius_d > radius_c;
                    // Check sixth diagonal ghost cell.
                    x_d[0] = patch_xlo[0] + (double(i-1) + double(1)/double(2))*dx[0]; // x coordinates of the point.
                    x_d[1] = patch_xlo[1] + (double(j+1) + double(1)/double(2))*dx[1]; // y coordinates of the point.
                    x_d[2] = patch_xlo[2] + (double(k-1) + double(1)/double(2))*dx[2]; // z coordinates of the point.
                    radius_d = sqrt(pow(x_d[0] - x_c, 2) + pow(x_d[1] - y_c, 2));
                    is_corner_ghost |= radius_d > radius_c;
                    
                    // Check seventh diagonal ghost cell.
                    x_d[0] = patch_xlo[0] + (double(i-1) + double(1)/double(2))*dx[0]; // x coordinates of the point.
                    x_d[1] = patch_xlo[1] + (double(j-1) + double(1)/double(2))*dx[1]; // y coordinates of the point.
                    x_d[2] = patch_xlo[2] + (double(k+1) + double(1)/double(2))*dx[2]; // z coordinates of the point.
                    radius_d = sqrt(pow(x_d[0] - x_c, 2) + pow(x_d[1] - y_c, 2));
                    is_corner_ghost |= radius_d > radius_c;
                    
                    // Check eighth diagonal ghost cell.
                    x_d[0] = patch_xlo[0] + (double(i-1) + double(1)/double(2))*dx[0]; // x coordinates of the point.
                    x_d[1] = patch_xlo[1] + (double(j-1) + double(1)/double(2))*dx[1]; // y coordinates of the point.
                    x_d[2] = patch_xlo[2] + (double(k-1) + double(1)/double(2))*dx[2]; // z coordinates of the point.
                    radius_d = sqrt(pow(x_d[0] - x_c, 2) + pow(x_d[1] - y_c, 2));
                    is_corner_ghost |= radius_d > radius_c;
                    
                    if ((fabs(x_p - x[0]) < (double(d_num_immersed_boundary_ghosts[0]))*dx[0]) ||
                        (fabs(y_p - x[1]) < (double(d_num_immersed_boundary_ghosts[1]))*dx[1]) ||
                    //    (fabs(z_p - x[2]) < (double(d_num_immersed_boundary_ghosts[2]))*dx[2]) ||
                        is_corner_ghost)
                    {
                        mask[idx]   = int(IB_MASK::IB_GHOST);
                        dist[idx]   = radius_c - radius;
                        norm_0[idx] = (x[0] - x_c)/radius;
                        norm_1[idx] = (x[1] - y_c)/radius;
                        //check this, normal does not depend on z 
                        norm_2[idx] = double(0);
                    }
                    else
                    {
                        mask[idx]   = int(IB_MASK::BODY);
                        dist[idx]   = double(0);
                        norm_0[idx] = double(0);
                        norm_1[idx] = double(0);
                        norm_2[idx] = double(0);
                    }
                }
                else
                {
                    mask[idx]   = int(IB_MASK::FLUID);
                    dist[idx]   = double(0);
                    norm_0[idx] = double(0);
                    norm_1[idx] = double(0);                        
                    norm_2[idx] = double(0);
                 }
            }
        }
    }
}


void
ImmersedBoundaries::generateSurfaceTriangulation(
    std::vector<std::array<Real, 3> >& nodes,
    std::vector<std::array<int, 3> >& connectivities,
    std::vector<int>& component_ids)
{
    NULL_USE(nodes);
    NULL_USE(connectivities);
    NULL_USE(component_ids);
}
