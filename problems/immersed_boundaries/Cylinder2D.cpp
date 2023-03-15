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
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_surface_normal,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_ip_index,      // AFK 03/14/23
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_ip_corr)      // AFK 03/14/2
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

    double* ip_location_index_0 = data_ip_index->getPointer(0);   // AFK 03/14/23
    double* ip_location_index_1 = data_ip_index->getPointer(1);
    double* ip_ratio_0          = data_ip_corr->getPointer(0);   // AFK 03/14/23
    double* ip_ratio_1          = data_ip_corr->getPointer(1);
    
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
    
    const double radius_c = double(1)/double(20);
    const double x_c = half;
    const double y_c = half;
    
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
            x[0] = patch_xlo[0] + (double(i) + double(1)/double(2))*dx[0]; // local x coordinates 
            x[1] = patch_xlo[1] + (double(j) + double(1)/double(2))*dx[1]; // local y coordinates
            
            const double radius = sqrt(pow(x[0] - x_c, 2) + pow(x[1] - y_c, 2)); // distance from the cylinder center
            const double theta  = atan2(x[1] - y_c, x[0] - x_c);                 // angle between x axis and a line passing through center and current cell

            if (radius < radius_c)      // condition that should be satisfied to be in cylinder.
            {
                double x_p = double(0); // x coordinates on the cylinder where y = x[1].
                double y_p = double(0); // y coordinates on the cylinder where x = x[0].
                
                double x_ip = double(0); // AFK 03/14/23 x coordinates of the image points 
                double y_ip = double(0); // AFK 03/14/23 y coordinates of the image points
 
                double x_ip_BL = double(0); // AFK 03/14/23 x coordinates of the BottomLeft Corner of image points
                double y_ip_BL = double(0); // AFK 03/14/23 y coordinates of the BottomLeft Corner of image points
                
                double d_ip = sqrt(double(2)) + double(0.000001);  // distance from cylinder boundary to the image point sqrt(2 + epsilon)

                if (x[0] > x_c)
                {
                    //x_p = x_c + sqrt(pow(radius_c, 2) - pow(radius*sin(theta), 2));
                    x_p = x_c + sqrt(pow(radius_c, 2) - pow(radius_c*sin(theta), 2));  // AFK 03/14/23 
                }
                else
                {
                    //x_p = x_c - sqrt(pow(radius_c, 2) - pow(radius*sin(theta), 2));
                    x_p = x_c - sqrt(pow(radius_c, 2) - pow(radius_c*sin(theta), 2));   // AFK 03/15/23 
                }
                
                if (x[1] > y_c)
                {
                    //y_p = y_c + sqrt(pow(radius_c, 2) - pow(radius*cos(theta), 2));
                    y_p = y_c + sqrt(pow(radius_c, 2) - pow(radius_c*cos(theta), 2));   // AFK 03/15/23 
                }
                else
                {
                    //y_p = y_c - sqrt(pow(radius_c, 2) - pow(radius*cos(theta), 2));
                    y_p = y_c - sqrt(pow(radius_c, 2) - pow(radius_c*cos(theta), 2)); // AFK 03/15/23  
                }
                
                if ((fabs(x_p - x[0]) < (double(d_num_immersed_boundary_ghosts[0]))*dx[0]) ||
                    (fabs(y_p - x[1]) < (double(d_num_immersed_boundary_ghosts[1]))*dx[1]))
                {
                    mask[idx]   = int(IB_MASK::IB_GHOST);
                    
                    //dist[idx]   = radius;
                    
                    dist[idx]   = radius_c - radius;      // AFK 03/14/23 instead of distance from center storing distance from the boundary
                    norm_0[idx] = (x[0] - x_c)/radius;    
                    norm_1[idx] = (x[1] - y_c)/radius;
                    
                    // AFK 03/14/23 Finding the image point local coordinates

                    if (x[0] > x_c)
                    {
                        x_ip    = x_c + sqrt(pow(radius_c + d_ip, 2) - pow((radius_c + d_ip)*sin(theta), 2));

                    }
                    else
                    { 
                        x_ip    = x_c - sqrt(pow(radius_c + d_ip, 2) - pow((radius_c + d_ip)*sin(theta), 2));

                    }

                    if (x[1] > y_c)
                    {
                        y_ip    = y_c + sqrt(pow(radius_c + d_ip, 2) - pow((radius_c + d_ip)*cos(theta), 2));

                    }
                    else
                    {
                        y_ip    = y_c - sqrt(pow(radius_c + d_ip, 2) - pow((radius_c + d_ip)*cos(theta), 2));
                    
                    }

                    // AFK 03/14/23  Finding the image point indexes and bilinear interpolation coefficients

                    ip_location_index_0[idx]    = floor((x_ip - half * dx[0]) / dx[0]);      // x axis index of the bottom-left cell in the interpolation stencil
                    ip_location_index_1[idx]    = floor((y_ip - half * dx[1]) / dx[1]);      // y axis index of the bottom-left cell in the interpolation stencil
                    
                    x_ip_BL = patch_xlo[0] + (ip_location_index_0[idx] + double(1)/double(2))*dx[0]; // AFK
                    y_ip_BL = patch_xlo[1] + (ip_location_index_1[idx] + double(1)/double(2))*dx[1]; // AFK
  
                    ip_ratio_0[idx]             = (x_ip - x_ip_BL) / dx[0];      // AFK x coefficient in the bilinear interpolation 
                    ip_ratio_1[idx]             = (y_ip - y_ip_BL) / dx[1];      // AFK y coefficient in the bilinear interpolation

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
