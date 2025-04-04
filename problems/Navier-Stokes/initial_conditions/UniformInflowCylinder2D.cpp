#include "apps/Navier-Stokes/NavierStokesInitialConditions.hpp"
/*
 * Set the data on the patch interior to some initial values.
 */
void
NavierStokesInitialConditions::initializeDataOnPatch(
    hier::Patch& patch,
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& conservative_variables,
    const double data_time,
    const bool initial_time)
{
    NULL_USE(data_time);
    
    if (d_project_name != "2D uniform flow")
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Can only initialize data for 'project_name' = '2D uniform flow'!\n"
            << "'project_name' = '"
            << d_project_name
            << "' is given."
            << std::endl);
    }
    
    if (d_dim != tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Dimension of problem should be 2!"
            << std::endl);
    }
    
    if (d_flow_model_type != FLOW_MODEL::SINGLE_SPECIES && d_flow_model_type != FLOW_MODEL::FIVE_EQN_ALLAIRE)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Flow model should be single-species or five-equation model by Allaire!"
            << std::endl);
    }
    
    if (d_flow_model_type == FLOW_MODEL::FIVE_EQN_ALLAIRE)
    {
        if (d_flow_model->getNumberOfSpecies() != 2)
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Number of species should be 2!"
                << std::endl);
        }
    }
    
    if (initial_time)
    {
        const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
            HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                patch.getPatchGeometry()));
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(patch_geom);
#endif
        
        // Get the dimensions of box that covers the interior of Patch.
        hier::Box patch_box = patch.getBox();
        const hier::IntVector patch_dims = patch_box.numberCells();
        const double* const dx = patch_geom->getDx();            //AFK
        const double* const patch_xlo = patch_geom->getXLower(); //AFK
        
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_cons_var = conservative_variables[0]->getGhostCellWidth();
        
        // Get the dimensions of the ghost cell boxes.
        const hier::Box ghost_box_cons_var = conservative_variables[0]->getGhostBox();
        const hier::IntVector ghostcell_dims_cons_var = ghost_box_cons_var.numberCells();
        
        /*
         * Initialize data for a 2D density wave advection problem.
         */
        
        if (d_flow_model_type == FLOW_MODEL::SINGLE_SPECIES)
        {
            HAMERS_SHARED_PTR<pdat::CellData<Real> > density      = conservative_variables[0];
            HAMERS_SHARED_PTR<pdat::CellData<Real> > momentum     = conservative_variables[1];
            HAMERS_SHARED_PTR<pdat::CellData<Real> > total_energy = conservative_variables[2];
            
            Real* rho   = density->getPointer(0);
            Real* rho_u = momentum->getPointer(0);
            Real* rho_v = momentum->getPointer(1);
            Real* E     = total_energy->getPointer(0);
            
            //Real gamma = Real(7)/Real(5);
            
            // Initial conditions.
            Real rho_inf = Real(1);
            Real u_inf   = Real(1);
            Real v_inf   = Real(1);
            Real p_inf   = Real(1);
            Real gamma   = Real(1);
            Real x_c     = Real(1);
            Real y_c     = Real(1); 
            Real D       = Real(1);
            Real r       = Real(1);
            Real theta   = Real(1);
            Real V_r     = Real(1); 
            Real V_theta = Real(1);
            Real u_ic    = Real(1);
            Real v_ic    = Real(1);
            Real p_ic    = Real(1);
            Real u       = Real(1);
            Real v       = Real(1);
            Real p       = Real(1);
            Real spongeR = Real(1);
            Real spongeL = Real(1);
            Real spongeB = Real(1);
            Real spongeT = Real(1);
            
            Real half = Real(1)/Real(2);
            
            if (d_initial_conditions_db != nullptr)
            {
                TBOX_ASSERT(d_initial_conditions_db->keyExists("rho_inf"));
                TBOX_ASSERT(d_initial_conditions_db->keyExists("u_inf"));
                TBOX_ASSERT(d_initial_conditions_db->keyExists("v_inf"));
                TBOX_ASSERT(d_initial_conditions_db->keyExists("p_inf"));
                
                rho_inf = d_initial_conditions_db->getReal("rho_inf");
                u_inf   = d_initial_conditions_db->getReal("u_inf");
                v_inf   = d_initial_conditions_db->getReal("v_inf");
                p_inf   = d_initial_conditions_db->getReal("p_inf");
                gamma   = d_initial_conditions_db->getReal("gamma");
                x_c     = d_initial_conditions_db->getReal("x_c");
                y_c     = d_initial_conditions_db->getReal("y_c");
                D       = d_initial_conditions_db->getReal("D");
                spongeL = d_initial_conditions_db->getReal("spongeL");
                spongeR = d_initial_conditions_db->getReal("spongeR");
                spongeB = d_initial_conditions_db->getReal("spongeB");
                spongeT = d_initial_conditions_db->getReal("spongeT");
                // std::cout << std::setprecision(17) << u_inf << std::endl;
                // std::cout << v_inf << std::endl;
                // std::cout << rho_inf << std::endl;
                // std::cout << p_inf << std::endl;
            }
            
            for (int j = -num_ghosts_cons_var[1]; j < patch_dims[1] + num_ghosts_cons_var[1]; j++)
            {
                for (int i = -num_ghosts_cons_var[0]; i < patch_dims[0] + num_ghosts_cons_var[0]; i++)
                {
                    // Compute index into linear data array.
                    int idx_cell = (i + num_ghosts_cons_var[0]) +
                        (j + num_ghosts_cons_var[1])*ghostcell_dims_cons_var[0];
                    
                    Real x[2];
                    x[0] = patch_xlo[0] + (Real(i) + half)*Real(dx[0]); // x coordinates of the point.
                    x[1] = patch_xlo[1] + (Real(j) + half)*Real(dx[1]);
                    
                    rho[idx_cell] = rho_inf;
                    
                    r       = std::pow((x[0] - x_c), Real(2)) + std::pow((x[1] - y_c), Real(2));
                    r       = std::pow(r, half);
                    theta   = std::atan((x[1] - y_c)/(x[0] - x_c)) + 30.0*M_PI/180.0; 
                    V_r     =  u_inf*(Real(1) - D*D/(Real(4)*r*r))*std::cos(theta);
                    V_theta = -u_inf*(Real(1) + D*D/(Real(4)*r*r))*std::sin(theta);
                    if (r < D/Real(4))
                    {
                        V_r     = Real(0);
                        V_theta = Real(0);
                    }
                    u_ic = (V_r*std::cos(theta) - V_theta*std::sin(theta));
                    v_ic = (V_r*std::sin(theta) + V_theta*std::cos(theta));
                    p_ic = p_inf + half*rho_inf*(u_inf*u_inf - (V_r*V_r + V_theta*V_theta));
                    
                    if(x[0] < x_c)
                    {
                        u = half * (u_inf + u_ic) + half * (u_ic - u_inf)*erf((x[0]-spongeL)/(D/Real(5)));
                        v = half * (v_inf + v_ic) + half * (v_ic - v_inf)*erf((x[0]-spongeL)/(D/Real(5)));
                        p = half * (p_inf + p_ic) + half * (p_ic - p_inf)*erf((x[0]-spongeL)/(D/Real(5)));
                    }
                    else
                    {
                        u = half * (u_inf + u_ic) - half * (u_ic - u_inf)*erf((x[0]-spongeR)/(D/Real(5)));
                        v = half * (v_inf + v_ic) - half * (v_ic - v_inf)*erf((x[0]-spongeR)/(D/Real(5)));
                        p = half * (p_inf + p_ic) - half * (p_ic - p_inf)*erf((x[0]-spongeR)/(D/Real(5)));
                    }
                    
                    if (x[1] < y_c)
                    {
                        u = half * (u_inf + u_ic) + half * (u_ic - u_inf)*erf((x[1]-spongeB)/(D/Real(5)));
                        v = half * (v_inf + v_ic) + half * (v_ic - v_inf)*erf((x[1]-spongeB)/(D/Real(5)));
                        p = half * (p_inf + p_ic) + half * (p_ic - p_inf)*erf((x[1]-spongeB)/(D/Real(5)));
                    }
                    else
                    {
                        u = half * (u_inf + u_ic) - half * (u_ic - u_inf)*erf((x[1]-spongeT)/(D/Real(5)));
                        v = half * (v_inf + v_ic) - half * (v_ic - v_inf)*erf((x[1]-spongeT)/(D/Real(5)));
                        p = half * (p_inf + p_ic) - half * (p_ic - p_inf)*erf((x[1]-spongeT)/(D/Real(5)));
                    }
                    
                    // rho_u[idx_cell] = rho_inf*u;
                    // rho_v[idx_cell] = rho_inf*v;
                    // E[idx_cell]     = p/(gamma - Real(1)) + half*rho_inf*(u*u + v*v);
                    
                    rho_u[idx_cell] = rho_inf*u_inf;
                    rho_v[idx_cell] = rho_inf*v_inf;
                    E[idx_cell]     = p_inf/(gamma - Real(1)) + half*rho_inf*(u_inf*u_inf + v_inf*v_inf);
                }
            }
        }
        else if (d_flow_model_type == FLOW_MODEL::FIVE_EQN_ALLAIRE)
        {
            HAMERS_SHARED_PTR<pdat::CellData<Real> > partial_density = conservative_variables[0];
            HAMERS_SHARED_PTR<pdat::CellData<Real> > momentum        = conservative_variables[1];
            HAMERS_SHARED_PTR<pdat::CellData<Real> > total_energy    = conservative_variables[2];
            HAMERS_SHARED_PTR<pdat::CellData<Real> > volume_fraction = conservative_variables[3];
            
            Real* Z_rho_1 = partial_density->getPointer(0);
            Real* Z_rho_2 = partial_density->getPointer(1);
            Real* rho_u   = momentum->getPointer(0);
            Real* rho_v   = momentum->getPointer(1);
            Real* E       = total_energy->getPointer(0);
            Real* Z_1     = volume_fraction->getPointer(0);
            Real* Z_2     = volume_fraction->getPointer(1);
            
            // Species 1.
            Real gamma_1 = Real(8)/Real(5); // 1.6
            
            // Species 2.
            Real gamma_2 = Real(7)/Real(5); // 1.4
            
            // Initial conditions.
            Real Z_rho_1_inf = Real(1);
            Real Z_rho_2_inf = Real(1);
            Real u_inf       = Real(1);
            Real v_inf       = Real(1);
            Real p_inf       = Real(1);
            Real Z_1_inf     = Real(1)/Real(2);
            Real Z_2_inf     = Real(1)/Real(2);
            
            if (d_initial_conditions_db != nullptr)
            {
                TBOX_ASSERT(d_initial_conditions_db->keyExists("Z_rho_1_inf"));
                TBOX_ASSERT(d_initial_conditions_db->keyExists("Z_rho_2_inf"));
                TBOX_ASSERT(d_initial_conditions_db->keyExists("u_inf"));
                TBOX_ASSERT(d_initial_conditions_db->keyExists("v_inf"));
                TBOX_ASSERT(d_initial_conditions_db->keyExists("p_inf"));
                TBOX_ASSERT(d_initial_conditions_db->keyExists("Z_1_inf"));
                TBOX_ASSERT(d_initial_conditions_db->keyExists("Z_2_inf"));
                
                Z_rho_1_inf = d_initial_conditions_db->getReal("Z_rho_1_inf");
                Z_rho_2_inf = d_initial_conditions_db->getReal("Z_rho_2_inf");
                u_inf       = d_initial_conditions_db->getReal("u_inf");
                v_inf       = d_initial_conditions_db->getReal("v_inf");
                p_inf       = d_initial_conditions_db->getReal("p_inf");
                Z_1_inf     = d_initial_conditions_db->getReal("Z_1_inf");
                Z_2_inf     = d_initial_conditions_db->getReal("Z_2_inf");
            }
            
            const Real rho_inf = Z_rho_1_inf + Z_rho_2_inf;
            const Real gamma_m = Real(1)/(Z_1_inf/(gamma_1 - Real(1)) + Z_2_inf/(gamma_2 - Real(1))) + Real(1);
            
            for (int j = 0; j < patch_dims[1]; j++)
            {
                for (int i = 0; i < patch_dims[0]; i++)
                {
                    // Compute index into linear data array.
                    int idx_cell = i + j*patch_dims[0];
                    
                    Z_rho_1[idx_cell] = Z_rho_1_inf;
                    Z_rho_2[idx_cell] = Z_rho_2_inf;
                    rho_u[idx_cell]   = rho_inf*u_inf;
                    rho_v[idx_cell]   = rho_inf*v_inf;
                    E[idx_cell]       = p_inf/(gamma_m - Real(1)) + Real(1)/Real(2)*rho_inf*(u_inf*u_inf + v_inf*v_inf);
                    Z_1[idx_cell]     = Z_1_inf;
                    Z_2[idx_cell]     = Z_2_inf;
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
//    const double x_cen = 2.0;
//    const double y_cen = 2.0;
//    const double z_cen = 2.0;
//    
//    // Start creating the surface mesh.
//    const double radius_sphere = 10.0;
//    const double edge_length = 1.0;
//
//    const int n_theta_init = 7;
//    const double edge_length_tol = 1.5;
//    const double pi = 2.0*std::asin(1.0);
//    
//    const double x_start_L = -radius_sphere;
//    const double y_start_L = 0.0;
//    const double z_start_L = 0.0;
//    
//    const double x_start_R = radius_sphere;
//    const double y_start_R = 0.0;
//    const double z_start_R = 0.0;
//    
//    std::vector<std::array<double, 3>> nodes_L;
//    std::vector<std::array<double, 3>> nodes_R;
//    std::vector<std::array<int, 3>> connectivity_L;
//    std::vector<std::array<int, 3>> connectivity_R;
//    // Create point at origin.
//    nodes_L.push_back({x_start_L, y_start_L, z_start_L});
//    nodes_R.push_back({x_start_R, y_start_R, z_start_R});
//    
//    std::vector<int> n_thetas;
//    n_thetas.push_back(n_theta_init);
//    
//    // Create n_theta nodes around the origin at a distance of edge_length.
//    double shift = edge_length; // std::sqrt(radius_sphere*radius_sphere - radius*radius);
//    double radius = std::sqrt(radius_sphere*radius_sphere - (radius_sphere - shift)*(radius_sphere - shift));
//    if (radius > edge_length)
//    {
//        radius = edge_length;
//        shift = radius_sphere - std::sqrt(radius_sphere*radius_sphere - radius*radius);
//    }
//    for (int i = 0; i < n_thetas[0]; ++i)
//    {
//        const double angle = 2.0 * pi * double(i) / double(n_thetas[0]);
//        nodes_L.push_back({x_start_L + shift, y_start_L + radius * std::cos(angle), z_start_L + radius * std::sin(angle)});
//        nodes_R.push_back({x_start_R - shift, y_start_R + radius * std::cos(angle), z_start_R + radius * std::sin(angle)});
//    }
//    // Create the triangles around the origin to form the disc.
//    for (int i = 1; i < n_thetas[0]; ++i)
//    {
//        const int i0 = 1;
//        const int i1 = 1 + (i + 0);
//        const int i2 = 1 + (i + 1);
//        connectivity_L.push_back({i2, i1, i0});
//        connectivity_R.push_back({i0, i1, i2});
//    }
//    connectivity_L.push_back({1 + 1, 1 + n_thetas[0], 1});
//    connectivity_R.push_back({1, 1 + n_thetas[0], 1 + 1});
//    
//    const int n_r = 5;
//    int node_offset = 1;
//    double angle_shift = 0.0;
//    int j = 1;
//    shift += edge_length; // at next row.
//    double radius_old = radius;
//    radius = std::sqrt(radius_sphere*radius_sphere - (radius_sphere - shift)*(radius_sphere - shift));
//    if (radius - radius_old > edge_length)
//    {
//        radius = radius_old + edge_length;
//        shift = radius_sphere - std::sqrt(radius_sphere*radius_sphere - radius*radius);
//    }
//    while (shift < radius_sphere + edge_length)
//    {
//        shift = std::min(shift, radius_sphere);
//        radius = std::sqrt(radius_sphere*radius_sphere - (radius_sphere - shift)*(radius_sphere - shift));
//        bool increased_n_theta = false;
//        // Increase n_theta if circumference per n_theta is much larger than the edge_length.
//        if (2.0*pi*radius/double(n_thetas.back()) > edge_length_tol*edge_length)
//        {
//            n_thetas.push_back(n_thetas.back()*2);
//            increased_n_theta = true;
//        }
//        else
//        {
//            n_thetas.push_back(n_thetas.back());
//        }
//        
//        if (!increased_n_theta)
//        {
//            angle_shift += pi/double(n_thetas.back());
//        }
//        for (int i = 0; i < n_thetas.back(); ++i)
//        {
//            double angle = 2.0*pi*double(i)/double(n_thetas.back()) + angle_shift;
//            nodes_L.push_back({x_start_L + shift, y_start_L + radius*std::cos(angle), z_start_L + radius*std::sin(angle)});
//            nodes_R.push_back({x_start_R - shift, y_start_R + radius*std::cos(angle), z_start_R + radius*std::sin(angle)});
//        }
//        
//        // Create the triangles for different rows.
//        if (increased_n_theta)
//        {
//            const int n_theta_lo = n_thetas[n_thetas.size() - 2];
//            const int n_theta_hi = n_thetas.back();
//            for (int i = 0; i < n_theta_lo; ++i)
//            {
//                const int lo_0 = 1 + node_offset + (i + 0);
//                const int lo_1 = 1 + node_offset + (i + 1)%n_theta_lo;
//                const int hi_0 = 1 + node_offset + n_theta_lo + (2*i + 0);
//                const int hi_1 = 1 + node_offset + n_theta_lo + (2*i + 1);
//                const int hi_2 = 1 + node_offset + n_theta_lo + (2*i + 2)%n_theta_hi;
//                
//                connectivity_L.push_back({hi_1, hi_0, lo_0});
//                connectivity_L.push_back({lo_1, hi_1, lo_0});
//                connectivity_L.push_back({hi_2, hi_1, lo_1});
//                
//                connectivity_R.push_back({lo_0, hi_0, hi_1});
//                connectivity_R.push_back({lo_0, hi_1, lo_1});
//                connectivity_R.push_back({lo_1, hi_1, hi_2});
//            }
//            node_offset += n_theta_lo;
//        }
//        else
//        {
//            const int n_theta = n_thetas.back();
//            for (int i = 0; i < n_theta; ++i)
//            {
//                const int lo_0 = 1 + node_offset + (i + 0);
//                const int lo_1 = 1 + node_offset + (i + 1)%n_theta;
//                const int hi_0 = 1 + node_offset + n_theta + (i + 0);
//                const int hi_1 = 1 + node_offset + n_theta + (i + 1)%n_theta;
//                
//                connectivity_L.push_back({lo_1, hi_0, lo_0});
//                connectivity_L.push_back({hi_1, hi_0, lo_1});
//                
//                connectivity_R.push_back({lo_0, hi_0, lo_1});
//                connectivity_R.push_back({lo_1, hi_0, hi_1});
//            }
//            node_offset += n_theta;
//        }
//        j++;
//        shift += edge_length; // at next row.
//        radius_old = radius;
//        radius = std::sqrt(radius_sphere*radius_sphere - (radius_sphere - shift)*(radius_sphere - shift));
//        if (radius - radius_old > edge_length)
//        {
//            radius = radius_old + edge_length;
//            shift = radius_sphere - std::sqrt(radius_sphere*radius_sphere - radius*radius);
//        }
//    }
//    
//    // Merge the two meshes.
//    for (int i = 0; i < nodes_L.size(); ++i)
//    {
//        nodes.push_back(nodes_L[i]);
//    }
//    for (int i = 0; i < nodes_R.size(); ++i)
//    {
//        nodes.push_back(nodes_R[i]);
//    }
//    for (int i = 0; i < nodes.size(); ++i)
//    {
//        nodes[i][0] += x_cen;
//        nodes[i][1] += y_cen;
//        nodes[i][2] += z_cen;
//    }
//    
//    for (int i = 0; i < connectivity_L.size(); ++i)
//    {
//        connectivities.push_back(connectivity_L[i]);
//    }
//    const int n_nodes_L = static_cast<int>(nodes_L.size());
//    for (int i = 0; i < connectivity_R.size(); ++i)
//    {
//        connectivity_R[i][0] += n_nodes_L;
//        connectivity_R[i][1] += n_nodes_L;
//        connectivity_R[i][2] += n_nodes_L;
//    }
//    for (int i = 0; i < connectivity_R.size(); ++i)
//    {
//        connectivities.push_back(connectivity_R[i]);
//    }
//    
//    const int num_centroids = static_cast<int>(connectivities.size());
//    
//    component_ids.assign(num_centroids, 0);
}