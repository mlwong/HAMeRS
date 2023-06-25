#include "flow/flow_models/single-species/FlowModelImmersedBoundaryMethodSingleSpecies.hpp"

FlowModelImmersedBoundaryMethodSingleSpecies::FlowModelImmersedBoundaryMethodSingleSpecies(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_species,
    const int& num_eqn,
    const HAMERS_SHARED_PTR<ImmersedBoundaries>& immersed_boundaries,
    const HAMERS_SHARED_PTR<tbox::Database>& immersed_boundary_method_db,
    const HAMERS_SHARED_PTR<EquationOfStateMixingRules>& equation_of_state_mixing_rules):
        FlowModelImmersedBoundaryMethod(
            object_name,
            dim,
            grid_geometry,
            num_species,
            num_eqn,
            immersed_boundaries,
            immersed_boundary_method_db,
            equation_of_state_mixing_rules)
{
    /*
     * Read the body density.
     */
    
    if (immersed_boundary_method_db->keyExists("body_density"))
    {
        d_rho_body = immersed_boundary_method_db->getReal("body_density");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelImmersedBoundaryMethodSingleSpecies::FlowModelImmersedBoundaryMethodSingleSpecies()\n"
            << "Required 'body_density' entry from input database missing."
            << std::endl);
    }
    
    /*
     * Read the body velocity.
     */
    
    if (immersed_boundary_method_db->keyExists("body_velocity"))
    {
        d_vel_body = immersed_boundary_method_db->getRealVector("body_velocity");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelImmersedBoundaryMethodSingleSpecies::FlowModelImmersedBoundaryMethodSingleSpecies()\n"
            << "Required 'body_velocity' entry from input database missing."
            << std::endl);
    }
    
    if (static_cast<int>(d_vel_body.size()) != d_dim.getValue())
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelImmersedBoundaryMethodSingleSpecies::FlowModelImmersedBoundaryMethodSingleSpecies()\n"
            << "Size of 'body_velocity' entry from input database does not match problem dimension."
            << std::endl);
    }
    
    /*
     * Read the body pressure.
     */
    
    if (immersed_boundary_method_db->keyExists("body_pressure"))
    {
        d_p_body = immersed_boundary_method_db->getReal("body_pressure");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelImmersedBoundaryMethodSingleSpecies::FlowModelImmersedBoundaryMethodSingleSpecies()\n"
            << "Required 'body_pressure' entry from input database missing."
            << std::endl);
    }
    
    /*
     * Read the problem type. (Euler or NS)
     */
    
    if (immersed_boundary_method_db->keyExists("problem_type"))
    {
        d_type = immersed_boundary_method_db->getString("problem_type");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelImmersedBoundaryMethodSingleSpecies::FlowModelImmersedBoundaryMethodSingleSpecies()\n"
            << "Required 'problem_type' entry from input database missing."
            << std::endl);
    }
    
    std::vector<Real*> thermo_properties_ptr;
    std::vector<const Real*> thermo_properties_const_ptr;
    
    const int num_thermo_properties = d_equation_of_state_mixing_rules->
        getNumberOfSpeciesThermodynamicProperties();
    
    thermo_properties_ptr.reserve(num_thermo_properties);
    thermo_properties_const_ptr.reserve(num_thermo_properties);
    d_thermo_properties.resize(num_thermo_properties);
    
    for (int ti = 0; ti < num_thermo_properties; ti++)
    {
        thermo_properties_ptr.push_back(&d_thermo_properties[ti]);
        thermo_properties_const_ptr.push_back(&d_thermo_properties[ti]);
    }
    
    d_equation_of_state_mixing_rules->getSpeciesThermodynamicProperties(
        thermo_properties_ptr,
        0);
    
    /*
     * Compute the values of conservative variables inside the body.
     */
    
    const Real T_body = d_equation_of_state_mixing_rules->getEquationOfState()->
        getTemperature(
            &d_rho_body,
            &d_p_body,
            thermo_properties_const_ptr);
    
    const Real epsilon_body = d_equation_of_state_mixing_rules->getEquationOfState()->
        getInternalEnergyFromTemperature(
            &d_rho_body,
            &T_body,
            thermo_properties_const_ptr);
    
    d_mom_body.resize(d_dim.getValue());
    
    if (d_dim == tbox::Dimension(1))
    {
        const Real& u_body = d_vel_body[0];
        
        d_mom_body[0] = d_rho_body*u_body;;
        
        d_E_body = d_rho_body*(epsilon_body + Real(1)/Real(2)*
            (u_body*u_body));
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const Real& u_body = d_vel_body[0];
        const Real& v_body = d_vel_body[1];
        
        d_mom_body[0] = d_rho_body*u_body;
        d_mom_body[1] = d_rho_body*v_body;
        
        d_E_body = d_rho_body*(epsilon_body + Real(1)/Real(2)*
            (u_body*u_body + v_body*v_body));
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const Real& u_body = d_vel_body[0];
        const Real& v_body = d_vel_body[1];
        const Real& w_body = d_vel_body[2];
        
        d_mom_body[0] = d_rho_body*u_body;
        d_mom_body[1] = d_rho_body*v_body;
        d_mom_body[2] = d_rho_body*w_body;
        
        d_E_body = d_rho_body*(epsilon_body + Real(1)/Real(2)*
            (u_body*u_body + v_body*v_body + w_body*w_body));
    }
}


/*
 * Set the immersed boundary method ghost cells for the cell data of conservative variables.
 */
void FlowModelImmersedBoundaryMethodSingleSpecies::setConservativeVariablesCellDataImmersedBoundaryGhosts(
    const hier::Patch& patch,
    const double data_time,
    const bool initial_time,
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& conservative_var_data,
    const HAMERS_SHARED_PTR<pdat::CellData<int> >& data_mask,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_wall_distance,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_surface_normal,
    const HAMERS_SHARED_PTR<pdat::CellData<int> >& data_ip_index,          // AFK
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_ip_corr,           // AFK
    const hier::IntVector& offset_cons_var,
    const hier::IntVector& offset_IB,
    const hier::IntVector& ghostcell_dims_cons_var,
    const hier::IntVector& ghostcell_dims_IB,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims)
{
    NULL_USE(data_time);
    
    if (!initial_time)
    {
        /*
         * Get the grid spacings and the lower coordinates.
         */
        
        const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
            HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                patch.getPatchGeometry()));
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(patch_geom);
#endif
        
        const double* const dx = patch_geom->getDx();
        const double* const patch_xlo = patch_geom->getXLower();
        
        /*
         * Get the data of the conservative variables.
         */
        
        const HAMERS_SHARED_PTR<pdat::CellData<Real> > data_density      = conservative_var_data[0];
        const HAMERS_SHARED_PTR<pdat::CellData<Real> > data_momentum     = conservative_var_data[1];
        const HAMERS_SHARED_PTR<pdat::CellData<Real> > data_total_energy = conservative_var_data[2];
        
        /*
         * Get the pointers to the data.
         */
        
        Real* rho = data_density->getPointer(0);
        Real* E   = data_total_energy->getPointer(0);
        
        int* mask = data_mask->getPointer(0);
        Real* dist = data_wall_distance->getPointer(0);  // AFK 03/15/23 distance variable uncommented 
        
        const Real& rho_body = d_rho_body;
        const Real& E_body   = d_E_body;
        
        if (d_dim == tbox::Dimension(1))
        {
            const Real& rho_u_body = d_mom_body[0];
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_dim_0 = domain_dims[0];
            
            const int offset_0_cons_var = offset_cons_var[0];
            const int offset_0_IB = offset_IB[0];
            
            // Get the pointers to the data.
            Real* rho_u = data_momentum->getPointer(0);
            
            // Real* norm_0 = data_surface_normal->getPointer(0);
            
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_cons_var = i + offset_0_cons_var;
                const int idx_IB = i + offset_0_IB;
                
                if (mask[idx_IB] == int(IB_MASK::IB_GHOST))
                {
                    // NEED TO BE CHANGED!!!
                    rho[idx_cons_var]   = rho_body;
                    rho_u[idx_cons_var] = rho_u_body;
                    E[idx_cons_var]     = E_body;
                }
                else if (mask[idx_IB] == int(IB_MASK::BODY))
                {
                    rho[idx_cons_var]   = rho_body;
                    rho_u[idx_cons_var] = rho_u_body;
                    E[idx_cons_var]     = E_body;
                }
            }
        }
        else if (d_dim == tbox::Dimension(2))
        {
            const Real& rho_u_body = d_mom_body[0];
            const Real& rho_v_body = d_mom_body[1];
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_lo_1 = domain_lo[1];
            const int domain_dim_0 = domain_dims[0];
            const int domain_dim_1 = domain_dims[1];
            
            const int offset_0_cons_var = offset_cons_var[0];
            const int offset_1_cons_var = offset_cons_var[1];
            const int ghostcell_dim_0_cons_var = ghostcell_dims_cons_var[0];
            
            const int offset_0_IB = offset_IB[0];
            const int offset_1_IB = offset_IB[1];
            const int ghostcell_dim_0_IB = ghostcell_dims_IB[0];
            
            // Get the pointers to the data.
            Real* rho_u = data_momentum->getPointer(0);
            Real* rho_v = data_momentum->getPointer(1);
            
            Real* norm_0 = data_surface_normal->getPointer(0);  // AFK 03/15/23 norm_0 variable uncommented
            Real* norm_1 = data_surface_normal->getPointer(1);  // AFK 03/15/23 norm_1 variable uncommented
            
            //int* ip_location_index_0 = data_ip_index->getPointer(0);  // AFK 03/15/23 
            //int* ip_location_index_1 = data_ip_index->getPointer(1);  // AFK
            
            //Real* ip_ratio_0          = data_ip_corr->getPointer(0);  // AFK 03/15/23 
            //Real* ip_ratio_1          = data_ip_corr->getPointer(1);  // AFK
            
            Real one  = Real(1);                                    // AFK one
            Real half = Real(1)/Real(2);
            Real d_ip = sqrt(Real(2))*Real(dx[0]) + HAMERS_REAL_EPSILON;  // AFK distance from cylinder boundary to the image point sqrt(2 + epsilon)
                                                                    // ML: assuming isotropic grid cells
            
            Real x_ip = Real(0);  // AFK x and y locations of image points
            Real y_ip = Real(0);
            
            int ip_location_index_0 = 0;  // x index of bottom left in the bilinear interpolation stencil
            int ip_location_index_1 = 0;  // y index of bottom left in the bilinear interpolation stencil
            
            Real x_ip_BL = Real(0);  // AFK x location of the bottom left cell in bilinear interpolation stencil
            Real y_ip_BL = Real(0);  // AFK y location of the bottom left cell in bilinear interpolation stencil
            
            Real ip_ratio_0 = Real(0); // AFK x ratio of the bilinear interpolation
            Real ip_ratio_1 = Real(0); // AFK y ratio of the bilinear interpolation
            
            Real rho_f1 = Real(0);  // AFK variables needed in interpolating density
            Real rho_f2 = Real(0);
            Real rho_ip = Real(0);
            Real rho_gc = Real(0);
            
            Real u_f1 = Real(0);  // AFK variables needed in interpolating u velocity
            Real u_f2 = Real(0);
            Real u_ip = Real(0);
            
            Real v_f1 = Real(0);  // AFK variables needed in interpolating v velocity
            Real v_f2 = Real(0);
            Real v_ip = Real(0);
            Real T_f1 = Real(0);  // AFK variables needed in interpolating v velocity
            Real T_f2 = Real(0);
            Real T_ip = Real(0);
            Real T_gc = Real(0);
            
            Real u_IP_BL = Real(0);  // AFK values of velocity in x axis for cells in the interpolation stencil
            Real u_IP_BR = Real(0);
            Real u_IP_TL = Real(0);
            Real u_IP_TR = Real(0);
            
            Real u_body  = Real(0);  // AFK u velocity of the body
            Real u_gc    = Real(0);  // AFK u velocity of the ghost cell
            
            Real v_IP_BL = Real(0);  // AFK values of velocity in y axis for cells in the interpolation stencil
            Real v_IP_BR = Real(0);
            Real v_IP_TL = Real(0);
            Real v_IP_TR = Real(0);
            
            Real v_body = Real(0);  // AFK v velocity of the body
            Real v_gc   = Real(0);  // AFK u velocity of the ghost cell
            
            Real epsilon_BL = Real(0);  // AFK
            Real epsilon_TL = Real(0);  // AFK
            Real epsilon_BR = Real(0);  // AFK
            Real epsilon_TR = Real(0);  // AFK
            
            Real p_BL = Real(0);  // AFK
            Real p_TL = Real(0);  // AFK
            Real p_BR = Real(0);  // AFK
            Real p_TR = Real(0);  // AFK
            
            Real T_BL = Real(0);  // AFK
            Real T_TL = Real(0);  // AFK
            Real T_BR = Real(0);  // AFK
            Real T_TR = Real(0);  // AFK
            
            Real epsilon_gc = Real(0);  // AFK
            Real E_gc       = Real(0);  // AFK
            
            Real vel_ip_n   =  Real(0);
            Real vel_ip_t   =  Real(0);
            Real vel_gc_n   =  Real(0);
            Real vel_gc_t   =  Real(0);
            Real vel_body_n =  Real(0);
            
            // !!! SHOULD NOT HARD-CODE c_v
            const Real c_v = Real(5) / Real(2);  // AFK isochoric specific heat to get internal energy from temperature
            
            std::vector<Real*> thermo_properties_ptr;
            std::vector<const Real*> thermo_properties_const_ptr;
            
            const int num_thermo_properties = d_equation_of_state_mixing_rules->
                getNumberOfSpeciesThermodynamicProperties();
            
            thermo_properties_ptr.reserve(num_thermo_properties);
            thermo_properties_const_ptr.reserve(num_thermo_properties);
            d_thermo_properties.resize(num_thermo_properties);
            
            for (int ti = 0; ti < num_thermo_properties; ti++)
            {
                thermo_properties_ptr.push_back(&d_thermo_properties[ti]);
                thermo_properties_const_ptr.push_back(&d_thermo_properties[ti]);
            }
            
            d_equation_of_state_mixing_rules->getSpeciesThermodynamicProperties(
                thermo_properties_ptr,
                0);
            
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_cons_var = (i + offset_0_cons_var) +
                        (j + offset_1_cons_var)*ghostcell_dim_0_cons_var;
                    
                    const int idx_IB = (i + offset_0_IB) +
                        (j + offset_1_IB)*ghostcell_dim_0_IB;
                    
                    // Compute the coordinates.
                    Real x[2];
                    x[0] = patch_xlo[0] + (Real(i) + half)*Real(dx[0]); // local x coordinates
                    x[1] = patch_xlo[1] + (Real(j) + half)*Real(dx[1]); // local y coordinates
                    
                    if (mask[idx_IB] == int(IB_MASK::IB_GHOST))
                    {
                        // General computation of image point locations.
                        x_ip = x[0] + (dist[idx_IB] + d_ip)*norm_0[idx_IB];
                        y_ip = x[1] + (dist[idx_IB] + d_ip)*norm_1[idx_IB];
                        
                        ip_location_index_0 = int(floor((x_ip - patch_xlo[0] - half * Real(dx[0])) / Real(dx[0])));
                        ip_location_index_1 = int(floor((y_ip - patch_xlo[1] - half * Real(dx[1])) / Real(dx[1])));
                        
                        x_ip_BL = patch_xlo[0] + (Real(ip_location_index_0) + half) * Real(dx[0]);
                        y_ip_BL = patch_xlo[1] + (Real(ip_location_index_1) + half) * Real(dx[1]);
                        
                        ip_ratio_0 = (x_ip - x_ip_BL) / Real(dx[0]);
                        ip_ratio_1 = (y_ip - y_ip_BL) / Real(dx[1]);
                        
                        // AFK declaration of the indexes of the cells in interpolation
                        const int idx_cons_var_BL  = (ip_location_index_0 + offset_0_cons_var) +
                                                     (ip_location_index_1 + offset_1_cons_var) * ghostcell_dim_0_cons_var;
                        
                        const int idx_cons_var_BR  = idx_cons_var_BL + 1;
                        const int idx_cons_var_TL  = idx_cons_var_BL + ghostcell_dim_0_cons_var;
                        const int idx_cons_var_TR  = idx_cons_var_BL + ghostcell_dim_0_cons_var + 1;
                        
                        // AFK bilinear interpolation to find image point density value
                        rho_f1 = (one - ip_ratio_0) * rho[idx_cons_var_BL] + ip_ratio_0 * rho[idx_cons_var_BR];
                        rho_f2 = (one - ip_ratio_0) * rho[idx_cons_var_TL] + ip_ratio_0 * rho[idx_cons_var_TR];
                        rho_ip = (one - ip_ratio_1) * rho_f1               + ip_ratio_1 * rho_f2;
                        rho_gc = rho_ip;
                        
                        // AFK bilinear interpolation to find image point x velocity value
                        u_IP_BL = rho_u[idx_cons_var_BL] / rho[idx_cons_var_BL];
                        u_IP_BR = rho_u[idx_cons_var_BR] / rho[idx_cons_var_BR];
                        u_IP_TL = rho_u[idx_cons_var_TL] / rho[idx_cons_var_TL];
                        u_IP_TR = rho_u[idx_cons_var_TR] / rho[idx_cons_var_TR];
                        
                        //u_body     = rho_u_body / rho[idx_cons_var];
                        
                        u_f1 = (one - ip_ratio_0) * u_IP_BL + ip_ratio_0 * u_IP_BR;
                        u_f2 = (one - ip_ratio_0) * u_IP_TL + ip_ratio_0 * u_IP_TR;
                        u_ip = (one - ip_ratio_1) * u_f1    + ip_ratio_1 * u_f2;
                        
                        // AFK bilinear interpolation to find image point x velocity value
                        v_IP_BL = rho_v[idx_cons_var_BL] / rho[idx_cons_var_BL];
                        v_IP_BR = rho_v[idx_cons_var_BR] / rho[idx_cons_var_BR];
                        v_IP_TL = rho_v[idx_cons_var_TL] / rho[idx_cons_var_TL];
                        v_IP_TR = rho_v[idx_cons_var_TR] / rho[idx_cons_var_TR];
                        
                        //v_body = rho_v_body / rho[idx_cons_var];
                        
                        v_f1 = (one - ip_ratio_0) * v_IP_BL + ip_ratio_0 * v_IP_BR; 
                        v_f2 = (one - ip_ratio_0) * v_IP_TL + ip_ratio_0 * v_IP_TR;
                        v_ip = (one - ip_ratio_1) * v_f1    + ip_ratio_1 * v_f2;
                        
                        // Compute velocities at the ghost cell for Euler equation (no-slip for normal direction and slip for tangential direction)
                        if(d_type == "Euler")
                        {
                            vel_ip_n =  u_ip * norm_0[idx_IB] + v_ip * norm_1[idx_IB];
                            vel_ip_t = -u_ip * norm_1[idx_IB] + v_ip * norm_0[idx_IB];
                            
                            vel_gc_n = vel_ip_n - ((d_ip + dist[idx_IB]) / d_ip) * (vel_ip_n);
                            vel_gc_t = vel_ip_t;
                            
                            u_gc = vel_gc_n * norm_0[idx_IB] - vel_gc_t * norm_1[idx_IB];
                            v_gc = vel_gc_n * norm_1[idx_IB] + vel_gc_t * norm_0[idx_IB]; 
                        }
                        else
                        {
                            u_gc = u_ip - ((d_ip + dist[idx_IB])/d_ip) * (u_ip);
                            v_gc = v_ip - ((d_ip + dist[idx_IB])/d_ip) * (v_ip);
                        }
                        
                        // Compute velocities at the ghost cell for NS equations (no-slip bcs)
                        // u_gc       = u_ip - ((d_ip + dist[idx_IB])/d_ip) * (u_ip);
                        // v_gc       = v_ip - ((d_ip + dist[idx_IB])/d_ip) * (v_ip);
                        
                        // AFK bilinear interpolation to find image point energy
                        epsilon_BL = E[idx_cons_var_BL] / rho[idx_cons_var_BL] - half*(u_IP_BL*u_IP_BL + v_IP_BL*v_IP_BL);
                        epsilon_TL = E[idx_cons_var_TL] / rho[idx_cons_var_TL] - half*(u_IP_TL*u_IP_TL + v_IP_TL*v_IP_TL);
                        epsilon_BR = E[idx_cons_var_BR] / rho[idx_cons_var_BR] - half*(u_IP_BR*u_IP_BR + v_IP_BR*v_IP_BR);
                        epsilon_TR = E[idx_cons_var_TR] / rho[idx_cons_var_TR] - half*(u_IP_TR*u_IP_TR + v_IP_TR*v_IP_TR);
                        
                        // SHOULD MAKE the following GENERAL using the equation of state object
                        T_BL    = epsilon_BL / c_v;
                        T_TL    = epsilon_TL / c_v;
                        T_BR    = epsilon_BR / c_v;
                        T_TR    = epsilon_TR / c_v;
                        
                        T_f1 = (one - ip_ratio_0) * T_BL + ip_ratio_0 * T_BR;
                        T_f2 = (one - ip_ratio_0) * T_TL + ip_ratio_0 * T_TR;
                        T_ip = (one - ip_ratio_1) * T_f1 + ip_ratio_1 * T_f2;
                        
                        T_gc = T_ip;
                        
                        E_gc = rho_gc * ( c_v * T_gc + half*(u_gc*u_gc + v_gc*v_gc)); // Calculating the energy at the ghost cell value by hard coding the internal energy
                        
                        rho[idx_cons_var]   = rho_gc;
                        rho_u[idx_cons_var] = rho[idx_cons_var] * u_gc;
                        rho_v[idx_cons_var] = rho[idx_cons_var] * v_gc;
                        E[idx_cons_var]     = E_gc;
                    }
                    else if (mask[idx_IB] == int(IB_MASK::BODY))
                    {
                        rho[idx_cons_var]   = rho_body;
                        rho_u[idx_cons_var] = rho_u_body;
                        rho_v[idx_cons_var] = rho_v_body;
                        E[idx_cons_var]     = E_body;
                    }
                }
            }
        }
        else if (d_dim == tbox::Dimension(3))
        {
            const Real& rho_u_body = d_mom_body[0];
            const Real& rho_v_body = d_mom_body[1];
            const Real& rho_w_body = d_mom_body[2];
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_lo_1 = domain_lo[1];
            const int domain_lo_2 = domain_lo[2];
            const int domain_dim_0 = domain_dims[0];
            const int domain_dim_1 = domain_dims[1];
            const int domain_dim_2 = domain_dims[2];
            
            const int offset_0_cons_var = offset_cons_var[0];
            const int offset_1_cons_var = offset_cons_var[1];
            const int offset_2_cons_var = offset_cons_var[2];
            const int ghostcell_dim_0_cons_var = ghostcell_dims_cons_var[0];
            const int ghostcell_dim_1_cons_var = ghostcell_dims_cons_var[1];
            
            const int offset_0_IB = offset_IB[0];
            const int offset_1_IB = offset_IB[1];
            const int offset_2_IB = offset_IB[2];
            const int ghostcell_dim_0_IB = ghostcell_dims_IB[0];
            const int ghostcell_dim_1_IB = ghostcell_dims_IB[1];
            
            // Get the pointers to the data.
            Real* rho_u = data_momentum->getPointer(0);
            Real* rho_v = data_momentum->getPointer(1);
            Real* rho_w = data_momentum->getPointer(2);
            
            // Real* norm_0 = data_surface_normal->getPointer(0);
            // Real* norm_1 = data_surface_normal->getPointer(1);
            // Real* norm_2 = data_surface_normal->getPointer(2);
            
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_cons_var = (i + offset_0_cons_var) +
                            (j + offset_1_cons_var)*ghostcell_dim_0_cons_var +
                            (k + offset_2_cons_var)*ghostcell_dim_0_cons_var*
                                ghostcell_dim_1_cons_var;
                        
                        const int idx_IB = (i + offset_0_IB) +
                            (j + offset_1_IB)*ghostcell_dim_0_IB +
                            (k + offset_2_IB)*ghostcell_dim_0_IB*
                                ghostcell_dim_1_IB;
                        
                        if (mask[idx_IB] == int(IB_MASK::IB_GHOST))
                        {
                            // NEED TO BE CHANGED!!!
                            rho[idx_cons_var]   = rho_body;
                            rho_u[idx_cons_var] = rho_u_body;
                            rho_v[idx_cons_var] = rho_v_body;
                            rho_w[idx_cons_var] = rho_w_body;
                            E[idx_cons_var]     = E_body;
                        }
                        else if (mask[idx_IB] == int(IB_MASK::BODY))
                        {
                            rho[idx_cons_var]   = rho_body;
                            rho_u[idx_cons_var] = rho_u_body;
                            rho_v[idx_cons_var] = rho_v_body;
                            rho_w[idx_cons_var] = rho_w_body;
                            E[idx_cons_var]     = E_body;
                        }
                    }
                }
            }
        }
    }
}
