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
    
    Real radius_c = half;
    Real x_c = Real(1);
    Real y_c = Real(1);

    if (d_initial_conditions_db != nullptr)
    {
        TBOX_ASSERT(d_initial_conditions_db->keyExists("x_c"));
        TBOX_ASSERT(d_initial_conditions_db->keyExists("y_c"));
        
        x_c = d_initial_conditions_db->getReal("x_c");
        y_c = d_initial_conditions_db->getReal("y_c");
        
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
                
                // Ghost cell check boolean variable
                bool is_ghost_cell = false;
                
                // Check cell for convective flux (direct neighbors)
                if ((fabs(x_p - x[0]) < (double(d_num_immersed_boundary_ghosts[0]))*dx[0]) || (fabs(y_p - x[1]) < (double(d_num_immersed_boundary_ghosts[1]))*dx[1]))
                 
                {
                    is_ghost_cell = true;
                }
                
                // Determine maximum ghost layers in x and y directions
                const int max_ghost_x = d_num_immersed_boundary_ghosts[0];
                const int max_ghost_y = d_num_immersed_boundary_ghosts[1];
                double x_d[2];

                // Check diagonal ghost cells for all possible layers
                if (!is_ghost_cell) 
                {
                    for (int gx = -max_ghost_x; gx <= max_ghost_x; gx++) 
                    {
                        for (int gy = -max_ghost_y; gy <= max_ghost_y; gy++) 
                        {
                            // Skip the center cell and Cartesian-aligned cases
                            if (gx == 0 && gy == 0) continue;
                            if (gx == 0 || gy == 0) continue;
                            
                            x_d[0] = patch_xlo[0] + (double(i + gx) + 0.5) * dx[0];
                            x_d[1] = patch_xlo[1] + (double(j + gy) + 0.5) * dx[1];
                            double radius_d = sqrt(pow(x_d[0] - x_c, 2) + pow(x_d[1] - y_c, 2));
                            
                            if (radius_d > radius_c) 
                            {
                                is_ghost_cell = true;
                                break;
                            }
                        }
                        if (is_ghost_cell) break;
                    }

                }

                // Apply ghost cell marking if any condition is met
                if (is_ghost_cell) {
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


void
ImmersedBoundaries::generateSurfaceTriangulation(
    std::vector<std::array<double, 3> >& nodes,
    std::vector<std::array<int, 3> >& connectivities,
    std::vector<std::array<double, 3> >& normal_nodes,
    std::vector<int>& component_ids)
{
    NULL_USE(nodes);
    NULL_USE(connectivities);
    NULL_USE(normal_nodes);
    NULL_USE(component_ids);
}