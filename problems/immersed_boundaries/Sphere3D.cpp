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
    int* mask    = data_mask->getPointer(0);
    Real* dist   = data_wall_distance->getPointer(0);
    Real* norm_0 = data_surface_normal->getPointer(0);
    Real* norm_1 = data_surface_normal->getPointer(1);
    Real* norm_2 = data_surface_normal->getPointer(2);
    
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
     * Set the parameters of the sphere here.
     */
    
    // These will be read from the input file.
    double x_c      = 1.0;
    double y_c      = 1.0;
    double z_c      = 1.0;
    double radius_c = 0.5;
    
    if (d_initial_conditions_db != nullptr)
    {
        TBOX_ASSERT(d_initial_conditions_db->keyExists("x_c"));
        TBOX_ASSERT(d_initial_conditions_db->keyExists("y_c"));
        TBOX_ASSERT(d_initial_conditions_db->keyExists("z_c"));
        
        x_c      = d_initial_conditions_db->getDouble("x_c");
        y_c      = d_initial_conditions_db->getDouble("y_c");
        z_c      = d_initial_conditions_db->getDouble("z_c");
        radius_c = d_initial_conditions_db->getDouble("radius");
    }
    
    for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++) 
    {
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++) 
        {
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
                
                // Distance from the sphere center.
                const double radius = sqrt(pow(x[0] - x_c, 2) + pow(x[1] - y_c, 2) + pow(x[2] - z_c, 2));
                
                if (radius < radius_c)  // Condition that should be satisfied to be in sphere
                {   
                    double x_p; // x coordinates on the cylinder where y = x[1] and z = x[2].
                    double y_p; // y coordinates on the cylinder where x = x[0] and z = x[2].
                    double z_p; // z coordinates on the cylinder where x = x[0] and y = x[1].
                    
                    // For checking ghost cell for convective flux.
                    if (x[0] > x_c)
                    {
                        x_p = x_c + sqrt(pow((radius_c), 2) - pow(radius, 2) + pow((x[0] - x_c),2));
                    }
                    else
                    {
                        x_p = x_c - sqrt(pow((radius_c), 2) - pow(radius, 2) + pow((x[0] - x_c),2));
                    }
                    
                    if (x[1] > y_c)
                    {
                        y_p = y_c + sqrt(pow((radius_c), 2) - pow(radius, 2) + pow((x[1] - y_c),2));
                    }
                    else
                    {
                        y_p = y_c - sqrt(pow((radius_c), 2) - pow(radius, 2) + pow((x[1] - y_c),2));
                    }
                    if (x[2] > z_c)
                    {
                        z_p = z_c + sqrt(pow((radius_c), 2) - pow(radius, 2) + pow((x[2] - z_c),2));
                    }
                    else
                    {
                        z_p = z_c - sqrt(pow((radius_c), 2) - pow(radius, 2) + pow((x[2] - z_c),2));
                    }
                    
                    // Determine maximum ghost layers in x and y directions
                    const int max_ghost_x = d_num_immersed_boundary_ghosts[0];
                    const int max_ghost_y = d_num_immersed_boundary_ghosts[1];
                    const int max_ghost_z = d_num_immersed_boundary_ghosts[2];
                    
                    if ((max_ghost_x != max_ghost_y) || (max_ghost_x != max_ghost_z) || (max_ghost_y != max_ghost_z))
                    {
                    TBOX_ERROR("num_immersed_boundary_ghosts should have the same value in x, y, and z directions\n");
                    }
                    
                    bool is_ghost_cell   = false;
                    bool is_corner_ghost = false;
                    
                    double x_d[3];
                    
                    for (int gx = 1; gx <= max_ghost_x; gx++) 
                    {
                        if ((fabs(x_p - x[0]) < (double(gx))*dx[0]) || (fabs(y_p - x[1]) < (double(gx))*dx[1]) 
                                                                    || (fabs(z_p - x[2]) < (double(gx))*dx[2])) // Ghost cells excluding corner ghost cells
                        {
                            is_ghost_cell = true;
                            break;
                        }
                        
                        x_d[0] = patch_xlo[0] + (double(i + gx) + double(1)/double(2)) * dx[0];
                        x_d[1] = patch_xlo[1] + (double(j + gx) + double(1)/double(2)) * dx[1];
                        x_d[2] = patch_xlo[2] + (double(k + gx) + double(1)/double(2)) * dx[2];
                        double radius_d_RTF = sqrt(pow(x_d[0] - x_c, 2) + pow(x_d[1] - y_c, 2) + pow(x_d[2] - z_c, 2));
                        
                        x_d[0] = patch_xlo[0] + (double(i - gx) + double(1)/double(2)) * dx[0];
                        x_d[1] = patch_xlo[1] + (double(j + gx) + double(1)/double(2)) * dx[1];
                        x_d[2] = patch_xlo[2] + (double(k + gx) + double(1)/double(2)) * dx[2];
                        double radius_d_LTF = sqrt(pow(x_d[0] - x_c, 2) + pow(x_d[1] - y_c, 2) + pow(x_d[2] - z_c, 2));
                        
                        x_d[0] = patch_xlo[0] + (double(i + gx) + double(1)/double(2)) * dx[0];
                        x_d[1] = patch_xlo[1] + (double(j - gx) + double(1)/double(2)) * dx[1];
                        x_d[2] = patch_xlo[2] + (double(k + gx) + double(1)/double(2)) * dx[2];
                        double radius_d_RBF = sqrt(pow(x_d[0] - x_c, 2) + pow(x_d[1] - y_c, 2) + pow(x_d[2] - z_c, 2));
                        
                        x_d[0] = patch_xlo[0] + (double(i - gx) + double(1)/double(2)) * dx[0];
                        x_d[1] = patch_xlo[1] + (double(j - gx) + double(1)/double(2)) * dx[1];
                        x_d[2] = patch_xlo[2] + (double(k + gx) + double(1)/double(2)) * dx[2];
                        double radius_d_LBF = sqrt(pow(x_d[0] - x_c, 2) + pow(x_d[1] - y_c, 2) + pow(x_d[2] - z_c, 2));
                        
                        x_d[0] = patch_xlo[0] + (double(i + gx) + double(1)/double(2)) * dx[0];
                        x_d[1] = patch_xlo[1] + (double(j + gx) + double(1)/double(2)) * dx[1];
                        x_d[2] = patch_xlo[2] + (double(k - gx) + double(1)/double(2)) * dx[2];
                        double radius_d_RTK = sqrt(pow(x_d[0] - x_c, 2) + pow(x_d[1] - y_c, 2) + pow(x_d[2] - z_c, 2));
                        
                        x_d[0] = patch_xlo[0] + (double(i - gx) + double(1)/double(2)) * dx[0];
                        x_d[1] = patch_xlo[1] + (double(j + gx) + double(1)/double(2)) * dx[1];
                        x_d[2] = patch_xlo[2] + (double(k - gx) + double(1)/double(2)) * dx[2];
                        double radius_d_LTK = sqrt(pow(x_d[0] - x_c, 2) + pow(x_d[1] - y_c, 2) + pow(x_d[2] - z_c, 2));
                        
                        x_d[0] = patch_xlo[0] + (double(i + gx) + double(1)/double(2)) * dx[0];
                        x_d[1] = patch_xlo[1] + (double(j - gx) + double(1)/double(2)) * dx[1];
                        x_d[2] = patch_xlo[2] + (double(k - gx) + double(1)/double(2)) * dx[2];
                        double radius_d_RBK = sqrt(pow(x_d[0] - x_c, 2) + pow(x_d[1] - y_c, 2) + pow(x_d[2] - z_c, 2));
                        
                        x_d[0] = patch_xlo[0] + (double(i - gx) + double(1)/double(2)) * dx[0];
                        x_d[1] = patch_xlo[1] + (double(j - gx) + double(1)/double(2)) * dx[1];
                        x_d[2] = patch_xlo[2] + (double(k - gx) + double(1)/double(2)) * dx[2];
                        double radius_d_LBK = sqrt(pow(x_d[0] - x_c, 2) + pow(x_d[1] - y_c, 2) + pow(x_d[2] - z_c, 2));
                        
                        if ((radius_d_RTF > radius_c) || (radius_d_LTF > radius_c) ||
                            (radius_d_RBF > radius_c) || (radius_d_LBF > radius_c) ||
                            (radius_d_RTK > radius_c) || (radius_d_LTK > radius_c) ||
                            (radius_d_RBK > radius_c) || (radius_d_LBK > radius_c))
                        {
                            is_corner_ghost = true;
                            break;
                        }
                    }
                    
                    if (is_ghost_cell || is_corner_ghost)
                    {
                        dist[idx]   = Real(radius_c - radius);
                        norm_0[idx] = Real((x[0] - x_c)/radius); // cos(theta) * sin(phi);
                        norm_1[idx] = Real((x[1] - y_c)/radius); // sin(theta) * sin(phi);
                        norm_2[idx] = Real((x[2] - z_c)/radius); // cos(phi);
                        
                        // Corner ghost cells required for viscous fluxes
                        if (is_corner_ghost)  
                        {
                            mask[idx]   = int(IB_MASK::IB_GHOST_CORNER);
                        }
                        else // Ghost cells required for convective fluxes
                        {
                            mask[idx]   = int(IB_MASK::IB_GHOST);
                        }
                    }
                    else
                    {
                        mask[idx]   = int(IB_MASK::BODY);
                        dist[idx]   = Real(0);
                        norm_0[idx] = Real(0);
                        norm_1[idx] = Real(0);
                        norm_2[idx] = Real(0);
                    }
                }
                else 
                {
                    mask[idx]   = int(IB_MASK::FLUID);
                    dist[idx]   = Real(0);
                    norm_0[idx] = Real(0);
                    norm_1[idx] = Real(0);
                    norm_2[idx] = Real(0);
                }
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
    double x_cen         = 1.0;
    double y_cen         = 1.0;
    double z_cen         = 1.0;
    double radius_sphere = 0.5;
    double edge_length   = 0.008;
    
    if (d_initial_conditions_db != nullptr)
    {
        x_cen         = d_initial_conditions_db->getDouble("x_c");
        y_cen         = d_initial_conditions_db->getDouble("y_c");
        z_cen         = d_initial_conditions_db->getDouble("z_c");
        radius_sphere = d_initial_conditions_db->getDouble("radius");
        edge_length   = d_initial_conditions_db->getDoubleWithDefault("edge_length", edge_length);
    }
    
    // Start creating the surface mesh.
    
    const int n_theta_init = 7;
    const double edge_length_tol = 1.5;
    const double pi = 2.0*std::asin(1.0);
    
    const double x_start_L = -radius_sphere;
    const double y_start_L = 0.0;
    const double z_start_L = 0.0;
    
    const double x_start_R = radius_sphere;
    const double y_start_R = 0.0;
    const double z_start_R = 0.0;
    
    std::vector<std::array<double, 3>> nodes_L;
    std::vector<std::array<double, 3>> nodes_R;
    std::vector<std::array<int, 3>> connectivity_L;
    std::vector<std::array<int, 3>> connectivity_R;
    // Create point at origin.
    nodes_L.push_back({x_start_L, y_start_L, z_start_L});
    nodes_R.push_back({x_start_R, y_start_R, z_start_R});
    
    std::vector<int> n_thetas;
    n_thetas.push_back(n_theta_init);
    
    // Create n_theta nodes around the origin at a distance of edge_length.
    double shift = edge_length; // std::sqrt(radius_sphere*radius_sphere - radius*radius);
    double radius = std::sqrt(radius_sphere*radius_sphere - (radius_sphere - shift)*(radius_sphere - shift));
    if (radius > edge_length)
    {
        radius = edge_length;
        shift = radius_sphere - std::sqrt(radius_sphere*radius_sphere - radius*radius);
    }
    for (int i = 0; i < n_thetas[0]; ++i)
    {
        const double angle = 2.0 * pi * double(i) / double(n_thetas[0]);
        nodes_L.push_back({x_start_L + shift, y_start_L + radius * std::cos(angle), z_start_L + radius * std::sin(angle)});
        nodes_R.push_back({x_start_R - shift, y_start_R + radius * std::cos(angle), z_start_R + radius * std::sin(angle)});
    }
    // Create the triangles around the origin to form the disc.
    for (int i = 1; i < n_thetas[0]; ++i)
    {
        const int i0 = 1;
        const int i1 = 1 + (i + 0);
        const int i2 = 1 + (i + 1);
        connectivity_L.push_back({i2, i1, i0});
        connectivity_R.push_back({i0, i1, i2});
    }
    connectivity_L.push_back({1 + 1, 1 + n_thetas[0], 1});
    connectivity_R.push_back({1, 1 + n_thetas[0], 1 + 1});
    
    const int n_r = 5;
    int node_offset = 1;
    double angle_shift = 0.0;
    int j = 1;
    shift += edge_length; // at next row.
    double radius_old = radius;
    radius = std::sqrt(radius_sphere*radius_sphere - (radius_sphere - shift)*(radius_sphere - shift));
    if (radius - radius_old > edge_length)
    {
        radius = radius_old + edge_length;
        shift = radius_sphere - std::sqrt(radius_sphere*radius_sphere - radius*radius);
    }
    while (shift < radius_sphere + edge_length)
    {
        shift = std::min(shift, radius_sphere);
        radius = std::sqrt(radius_sphere*radius_sphere - (radius_sphere - shift)*(radius_sphere - shift));
        bool increased_n_theta = false;
        // Increase n_theta if circumference per n_theta is much larger than the edge_length.
        if (2.0*pi*radius/double(n_thetas.back()) > edge_length_tol*edge_length)
        {
            n_thetas.push_back(n_thetas.back()*2);
            increased_n_theta = true;
        }
        else
        {
            n_thetas.push_back(n_thetas.back());
        }
        
        if (!increased_n_theta)
        {
            angle_shift += pi/double(n_thetas.back());
        }
        for (int i = 0; i < n_thetas.back(); ++i)
        {
            double angle = 2.0*pi*double(i)/double(n_thetas.back()) + angle_shift;
            nodes_L.push_back({x_start_L + shift, y_start_L + radius*std::cos(angle), z_start_L + radius*std::sin(angle)});
            nodes_R.push_back({x_start_R - shift, y_start_R + radius*std::cos(angle), z_start_R + radius*std::sin(angle)});
        }
        
        // Create the triangles for different rows.
        if (increased_n_theta)
        {
            const int n_theta_lo = n_thetas[n_thetas.size() - 2];
            const int n_theta_hi = n_thetas.back();
            for (int i = 0; i < n_theta_lo; ++i)
            {
                const int lo_0 = 1 + node_offset + (i + 0);
                const int lo_1 = 1 + node_offset + (i + 1)%n_theta_lo;
                const int hi_0 = 1 + node_offset + n_theta_lo + (2*i + 0);
                const int hi_1 = 1 + node_offset + n_theta_lo + (2*i + 1);
                const int hi_2 = 1 + node_offset + n_theta_lo + (2*i + 2)%n_theta_hi;
                
                connectivity_L.push_back({hi_1, hi_0, lo_0});
                connectivity_L.push_back({lo_1, hi_1, lo_0});
                connectivity_L.push_back({hi_2, hi_1, lo_1});
                
                connectivity_R.push_back({lo_0, hi_0, hi_1});
                connectivity_R.push_back({lo_0, hi_1, lo_1});
                connectivity_R.push_back({lo_1, hi_1, hi_2});
            }
            node_offset += n_theta_lo;
        }
        else
        {
            const int n_theta = n_thetas.back();
            for (int i = 0; i < n_theta; ++i)
            {
                const int lo_0 = 1 + node_offset + (i + 0);
                const int lo_1 = 1 + node_offset + (i + 1)%n_theta;
                const int hi_0 = 1 + node_offset + n_theta + (i + 0);
                const int hi_1 = 1 + node_offset + n_theta + (i + 1)%n_theta;
                
                connectivity_L.push_back({lo_1, hi_0, lo_0});
                connectivity_L.push_back({hi_1, hi_0, lo_1});
                
                connectivity_R.push_back({lo_0, hi_0, lo_1});
                connectivity_R.push_back({lo_1, hi_0, hi_1});
            }
            node_offset += n_theta;
        }
        j++;
        shift += edge_length; // at next row.
        radius_old = radius;
        radius = std::sqrt(radius_sphere*radius_sphere - (radius_sphere - shift)*(radius_sphere - shift));
        if (radius - radius_old > edge_length)
        {
            radius = radius_old + edge_length;
            shift = radius_sphere - std::sqrt(radius_sphere*radius_sphere - radius*radius);
        }
    }
    
    // Merge the two meshes.
    for (int i = 0; i < nodes_L.size(); ++i)
    {
        nodes.push_back(nodes_L[i]);
    }
    for (int i = 0; i < nodes_R.size(); ++i)
    {
        nodes.push_back(nodes_R[i]);
    }
    for (int i = 0; i < nodes.size(); ++i)
    {
        nodes[i][0] += x_cen;
        nodes[i][1] += y_cen;
        nodes[i][2] += z_cen;
        if (std::abs(nodes[i][0]) < HAMERS_EPSILON)
        {
            nodes[i][0] = HAMERS_EPSILON;
        }
        if (std::abs(nodes[i][1]) < HAMERS_EPSILON)
        {
            nodes[i][1] = HAMERS_EPSILON;
        }
        if (std::abs(nodes[i][2]) < HAMERS_EPSILON)
        {
            nodes[i][2] = HAMERS_EPSILON;
        }
        const double radius_node =
            std::sqrt(pow(nodes[i][0] - x_cen, 2) + pow(nodes[i][1] - y_cen, 2) + pow(nodes[i][2] - z_cen, 2));
        const std::array<double, 3> normal_node = {
            (nodes[i][0] - x_cen)/radius_node,
            (nodes[i][1] - y_cen)/radius_node,
            (nodes[i][2] - z_cen)/radius_node};
        normal_nodes.push_back(normal_node);
    }
    
    for (int i = 0; i < connectivity_L.size(); ++i)
    {
        connectivities.push_back(connectivity_L[i]);
    }
    const int n_nodes_L = static_cast<int>(nodes_L.size());
    for (int i = 0; i < connectivity_R.size(); ++i)
    {
        connectivity_R[i][0] += n_nodes_L;
        connectivity_R[i][1] += n_nodes_L;
        connectivity_R[i][2] += n_nodes_L;
    }
    for (int i = 0; i < connectivity_R.size(); ++i)
    {
        connectivities.push_back(connectivity_R[i]);
    }
    
    const int num_centroids = static_cast<int>(connectivities.size());
    
    component_ids.assign(num_centroids, 0);
}