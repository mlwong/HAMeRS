#include "flow/flow_models/single-species/FlowModelImmersedBoundaryMethodSingleSpecies.hpp"

FlowModelImmersedBoundaryMethodSingleSpecies::FlowModelImmersedBoundaryMethodSingleSpecies(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_species,
    const int& num_eqn,
    const HAMERS_SHARED_PTR<ImmersedBoundaries>& immersed_boundaries,
    const HAMERS_SHARED_PTR<tbox::Database>& immersed_boundary_method_db,
    const HAMERS_SHARED_PTR<EquationOfStateMixingRules>& equation_of_state_mixing_rules,
    const HAMERS_SHARED_PTR<EquationOfShearViscosityMixingRules> equation_of_shear_viscosity_mixing_rules,
    const HAMERS_SHARED_PTR<EquationOfBulkViscosityMixingRules> equation_of_bulk_viscosity_mixing_rules,
    const HAMERS_SHARED_PTR<EquationOfThermalConductivityMixingRules> equation_of_thermal_conductivity_mixing_rules):
        FlowModelImmersedBoundaryMethod(
            object_name,
            dim,
            grid_geometry,
            num_species,
            num_eqn,
            immersed_boundaries,
            immersed_boundary_method_db,
            equation_of_state_mixing_rules),
        d_equation_of_shear_viscosity_mixing_rules(equation_of_shear_viscosity_mixing_rules),
        d_equation_of_bulk_viscosity_mixing_rules(equation_of_bulk_viscosity_mixing_rules),
        d_equation_of_thermal_conductivity_mixing_rules(equation_of_thermal_conductivity_mixing_rules),
        d_surface_triangulation_integrated_F_p_x(0),
        d_surface_triangulation_integrated_F_p_y(0),
        d_surface_triangulation_integrated_F_p_z(0),
        d_surface_triangulation_integrated_F_v_x(0),
        d_surface_triangulation_integrated_F_v_y(0),
        d_surface_triangulation_integrated_F_v_z(0)
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
     * Read the body temperature.
     */
    
    if (immersed_boundary_method_db->keyExists("body_temperature"))
    {
        d_T_body = immersed_boundary_method_db->getReal("body_temperature");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelImmersedBoundaryMethodSingleSpecies::FlowModelImmersedBoundaryMethodSingleSpecies()\n"
            << "Required 'body_temperature' entry from input database missing."
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
    const hier::IntVector& offset_cons_var,
    const hier::IntVector& offset_IB,
    const hier::IntVector& ghostcell_dims_cons_var,
    const hier::IntVector& ghostcell_dims_IB,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims)
{
    NULL_USE(data_time);
    NULL_USE(initial_time);
    
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
    Real* dist = data_wall_distance->getPointer(0);
    
    const Real& rho_body  = d_rho_body;
    const Real& E_body    = d_E_body;
    const Real& T_body    = d_T_body;
    
    // Get the thermodynamic properties of the species.
    std::vector<const Real*> thermo_properties_ptr;
    thermo_properties_ptr.reserve(static_cast<int> (d_thermo_properties.size()));
    for (int ti = 0; ti < static_cast<int> (d_thermo_properties.size()); ti++)
    {
        thermo_properties_ptr.push_back(&d_thermo_properties[ti]);
    }
    
    constexpr Real half = Real(1)/Real(2);
    
    // Distance from cylinder boundary to the image point sqrt(2 + epsilon), isotropic grid cells are assuemd.
    
    if (d_dim == tbox::Dimension(1))
    {
        const Real d_ip = sqrt(Real(1))*Real(dx[0]) + HAMERS_REAL_EPSILON;
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
        if (fabs(dx[0] - dx[1]) > 10.0*std::numeric_limits<double>::epsilon())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelImmersedBoundaryMethodSingleSpecies::setConservativeVariablesCellDataImmersedBoundaryGhosts()\n"
                << "The grid is assumed to be isotropic but the grid spacings are different."
                << std::endl);
        }
        
        const Real dx_inv = Real(1)/Real(dx[0]);
        
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
        
        Real* norm_0 = data_surface_normal->getPointer(0);
        Real* norm_1 = data_surface_normal->getPointer(1);
        
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
                const Real x[2] = {
                    Real(patch_xlo[0]) + (Real(i) + half)*Real(dx[0]),
                    Real(patch_xlo[1]) + (Real(j) + half)*Real(dx[1])
                    };
                
                if ((mask[idx_IB] == int(IB_MASK::IB_GHOST)) || (mask[idx_IB] == int(IB_MASK::IB_GHOST_CORNER)))
                {
                    Real d_ip;
                    // First image point distance is set to sqrt(2)*dx + epsilon for corner ghost cells.
                    if (mask[idx_IB] == int(IB_MASK::IB_GHOST_CORNER))
                    {
                        d_ip   = std::sqrt(Real(2))*Real(dx[0]) + HAMERS_REAL_EPSILON;
                    }
                    // First image point distance is set to dx/maximum(norm) for ghost cells for convective fluxes.
                    else
                    {
                        Real norm_max = std::max(std::abs(norm_0[idx_IB]), std::abs(norm_1[idx_IB]));
                        d_ip  = Real(dx[0]) / norm_max + HAMERS_REAL_EPSILON;
                    }
                    
                    // Second image point distance.
                    Real d_ip2 = d_ip + Real(0.25)*Real(dx[0]);
                    
                    // Coordinates of the image point 1.
                    const Real x_ip  = x[0] + (dist[idx_IB] + d_ip)*norm_0[idx_IB];
                    const Real y_ip  = x[1] + (dist[idx_IB] + d_ip)*norm_1[idx_IB];
                    
                    // Coordinates of the image point 2.
                    const Real x_ip2 = x[0] + (dist[idx_IB] + d_ip2)*norm_0[idx_IB];
                    const Real y_ip2 = x[1] + (dist[idx_IB] + d_ip2)*norm_1[idx_IB];
                    
                    // Get indices of the cells in interpolation for image point 1.
                    int idx_ip_cons_var_BL, idx_ip_cons_var_BR, idx_ip_cons_var_TL, idx_ip_cons_var_TR;
                    int idx_ip_IB_BL, idx_ip_IB_BR, idx_ip_IB_TL, idx_ip_IB_TR;
                    Real x_ip_BL, y_ip_BL;
                    getBilinearInterpolationIndices2D(
                        idx_ip_cons_var_BL,
                        idx_ip_cons_var_BR,
                        idx_ip_cons_var_TL,
                        idx_ip_cons_var_TR,
                        idx_ip_IB_BL,
                        idx_ip_IB_BR, 
                        idx_ip_IB_TL, 
                        idx_ip_IB_TR,
                        x_ip_BL,
                        y_ip_BL,
                        x_ip,
                        y_ip,
                        Real(patch_xlo[0]),
                        Real(patch_xlo[1]),
                        offset_0_cons_var,
                        offset_1_cons_var,
                        ghostcell_dim_0_cons_var,
                        offset_0_IB,
                        offset_1_IB,
                        ghostcell_dim_0_IB,
                        Real(dx[0]),
                        dx_inv);
                    
                    // Checking first image point interpolation stencil to ensure only fluid cell values are used
                    if (mask[idx_ip_IB_BL] != int(IB_MASK::FLUID))
                    {
                        TBOX_ERROR("Error: Bottom-left cell is not FLUID at index:" << idx_ip_IB_BL 
                                    << "\n with mask value: " << mask[idx_ip_IB_BL]
                                    << "\n x_ip: " << x_ip
                                    << "\n y_ip: " << y_ip
                                    << "\n x: " << x[0]
                                    << "\n y: " << x[1]
                                    << "\n d_ip: " << d_ip/Real(dx[0])
                                    << "\n d_gc: " << dist[idx_IB]
                                    << "\n norm_0:" << norm_0[idx_IB]
                                    << "\n norm_1:" << norm_1[idx_IB]
                                    << "\n dx:" << Real(dx[0])
                                    << "\n dx_inv: " << dx_inv
                                    << "\n x_ip_BL:" << x_ip_BL
                                    << "\n y_ip_BL:" << y_ip_BL);
                    }

                    if (mask[idx_ip_IB_TL] != int(IB_MASK::FLUID))
                    {
                        TBOX_ERROR("Error: Top-left cell is not FLUID at index: " << idx_ip_IB_TL 
                                    << "\n with mask value: " << mask[idx_ip_IB_TL]
                                    << "\n x_ip: " << x_ip
                                    << "\n y_ip: " << y_ip
                                    << "\n x: " << x[0]
                                    << "\n y: " << x[1]
                                    << "\n d_ip: " << d_ip/Real(dx[0])
                                    << "\n d_gc: " << dist[idx_IB]
                                    << "\n norm_0: " << norm_0[idx_IB]
                                    << "\n norm_1: " << norm_1[idx_IB]
                                    << "\n dx: " << Real(dx[0])
                                    << "\n dx_inv: " << dx_inv
                                    << "\n x_ip_BL: " << x_ip_BL
                                    << "\n y_ip_BL: " << y_ip_BL);
                    }

                    if (mask[idx_ip_IB_BR] != int(IB_MASK::FLUID))
                    {
                        TBOX_ERROR("Error: Bottom-right cell is not FLUID at index: " << idx_ip_IB_BR 
                                    << "\n with mask value: " << mask[idx_ip_IB_BR]
                                    << "\n x_ip: " << x_ip
                                    << "\n y_ip: " << y_ip
                                    << "\n x: " << x[0]
                                    << "\n y: " << x[1]
                                    << "\n d_ip: " << d_ip/Real(dx[0])
                                    << "\n d_gc: " << dist[idx_IB]
                                    << "\n norm_0: " << norm_0[idx_IB]
                                    << "\n norm_1: " << norm_1[idx_IB]
                                    << "\n dx: " << Real(dx[0])
                                    << "\n dx_inv: " << dx_inv
                                    << "\n x_ip_BL: " << x_ip_BL
                                    << "\n y_ip_BL: " << y_ip_BL);
                    }

                    if (mask[idx_ip_IB_TR] != int(IB_MASK::FLUID))
                    {
                        TBOX_ERROR("Error: Top-right cell is not FLUID at index " << idx_ip_IB_TR 
                                    << "\n with mask value: " << mask[idx_ip_IB_TR]
                                    << "\n x_ip: " << x_ip
                                    << "\n y_ip: " << y_ip
                                    << "\n x: " << x[0]
                                    << "\n y: " << x[1]
                                    << "\n d_ip: " << d_ip/Real(dx[0])
                                    << "\n d_gc: " << dist[idx_IB]
                                    << "\n norm_0:" << norm_0[idx_IB]
                                    << "\n norm_1:" << norm_1[idx_IB]
                                    << "\n dx: " << Real(dx[0])
                                    << "\n dx_inv" << dx_inv
                                    << "\n x_ip_BL: " << x_ip_BL
                                    << "\n y_ip_BL: " << y_ip_BL);
                    }
                    
                    // Get indices of the cells in interpolation for image point 2.
                    int idx_ip2_cons_var_BL, idx_ip2_cons_var_BR, idx_ip2_cons_var_TL, idx_ip2_cons_var_TR;
                    int idx_ip2_IB_BL, idx_ip2_IB_BR, idx_ip2_IB_TL, idx_ip2_IB_TR;
                    Real x_ip2_BL, y_ip2_BL;
                    getBilinearInterpolationIndices2D(
                        idx_ip2_cons_var_BL,
                        idx_ip2_cons_var_BR,
                        idx_ip2_cons_var_TL,
                        idx_ip2_cons_var_TR,
                        idx_ip2_IB_BL,
                        idx_ip2_IB_BR, 
                        idx_ip2_IB_TL, 
                        idx_ip2_IB_TR,
                        x_ip2_BL,
                        y_ip2_BL,
                        x_ip2,
                        y_ip2,
                        Real(patch_xlo[0]),
                        Real(patch_xlo[1]),
                        offset_0_cons_var,
                        offset_1_cons_var,
                        ghostcell_dim_0_cons_var,
                        offset_0_IB,
                        offset_1_IB,
                        ghostcell_dim_0_IB,
                        Real(dx[0]),
                        dx_inv);
                    
                    // Bilinear interpolation to find image point x-component of velocity values.
                    const Real u_ip_BL = rho_u[idx_ip_cons_var_BL]/rho[idx_ip_cons_var_BL];
                    const Real u_ip_BR = rho_u[idx_ip_cons_var_BR]/rho[idx_ip_cons_var_BR];
                    const Real u_ip_TL = rho_u[idx_ip_cons_var_TL]/rho[idx_ip_cons_var_TL];
                    const Real u_ip_TR = rho_u[idx_ip_cons_var_TR]/rho[idx_ip_cons_var_TR];
                    
                    const Real u_ip2_BL = rho_u[idx_ip2_cons_var_BL]/rho[idx_ip2_cons_var_BL];
                    const Real u_ip2_BR = rho_u[idx_ip2_cons_var_BR]/rho[idx_ip2_cons_var_BR];
                    const Real u_ip2_TL = rho_u[idx_ip2_cons_var_TL]/rho[idx_ip2_cons_var_TL];
                    const Real u_ip2_TR = rho_u[idx_ip2_cons_var_TR]/rho[idx_ip2_cons_var_TR];
                    
                    const Real u_ip = bilinearInterpolate2D(
                        u_ip_BL,
                        u_ip_BR,
                        u_ip_TL,
                        u_ip_TR,
                        x_ip,
                        y_ip,
                        x_ip_BL,
                        y_ip_BL,
                        dx_inv);
                    
                    // Bilinear interpolation to find image point y-component of velocity values.
                    const Real v_ip_BL = rho_v[idx_ip_cons_var_BL]/rho[idx_ip_cons_var_BL];
                    const Real v_ip_BR = rho_v[idx_ip_cons_var_BR]/rho[idx_ip_cons_var_BR];
                    const Real v_ip_TL = rho_v[idx_ip_cons_var_TL]/rho[idx_ip_cons_var_TL];
                    const Real v_ip_TR = rho_v[idx_ip_cons_var_TR]/rho[idx_ip_cons_var_TR];
                    
                    const Real v_ip2_BL = rho_v[idx_ip2_cons_var_BL]/rho[idx_ip2_cons_var_BL];
                    const Real v_ip2_BR = rho_v[idx_ip2_cons_var_BR]/rho[idx_ip2_cons_var_BR];
                    const Real v_ip2_TL = rho_v[idx_ip2_cons_var_TL]/rho[idx_ip2_cons_var_TL];
                    const Real v_ip2_TR = rho_v[idx_ip2_cons_var_TR]/rho[idx_ip2_cons_var_TR];
                    
                    const Real v_ip = bilinearInterpolate2D(
                        v_ip_BL,
                        v_ip_BR,
                        v_ip_TL,
                        v_ip_TR,
                        x_ip,
                        y_ip,
                        x_ip_BL,
                        y_ip_BL,
                        dx_inv);
                    
                    // Compute velocity components at the ghost cell depending on the type of immersed boundary condition.
                    Real u_gc = Real(0); // x-component of velocity of the ghost cell
                    Real v_gc = Real(0); // y-component of velocity of the ghost cell
                    
                    if (d_bc_type_velocity == VELOCITY_IBC::SLIP)
                    {
                        const Real u_ip2 = bilinearInterpolate2D(
                            u_ip2_BL,
                            u_ip2_BR,
                            u_ip2_TL,
                            u_ip2_TR,
                            x_ip2,
                            y_ip2,
                            x_ip2_BL,
                            y_ip2_BL,
                            dx_inv);
                        
                        const Real v_ip2 = bilinearInterpolate2D(
                            v_ip2_BL,
                            v_ip2_BR,
                            v_ip2_TL,
                            v_ip2_TR,
                            x_ip2,
                            y_ip2,
                            x_ip2_BL,
                            y_ip2_BL,
                            dx_inv);
                        
                        // Given d_ip and d_ip2, interpolate to the velocity components (u_mirror and v_mirror)
                        // at the mirror image point at dist[idx_IB].
                        const Real diff_ip2_ip = d_ip2 - d_ip;
                        const Real diff_mirror_ip  = dist[idx_IB] - d_ip;
                        const Real diff_ip2_mirror = d_ip2 - dist[idx_IB];
                        
                        // x-component of velocity at the mirror image point.
                        const Real u_mirror = (diff_ip2_mirror*u_ip + diff_mirror_ip*u_ip2)/diff_ip2_ip;
                        // y-component of velocity at the mirror image point.
                        const Real v_mirror = (diff_ip2_mirror*v_ip + diff_mirror_ip*v_ip2)/diff_ip2_ip;
                        
                        // Velocity component normal to the boundary at the mirror image point.
                        const Real vel_mirror_n = dotProduct2D(u_mirror, v_mirror, norm_0[idx_IB], norm_1[idx_IB]);
                        
                        // No-penetration boundary condition.
                        u_gc = u_mirror - Real(2)*vel_mirror_n*norm_0[idx_IB];
                        v_gc = v_mirror - Real(2)*vel_mirror_n*norm_1[idx_IB];
                        
                        // Real vel_ip_n =  u_ip*norm_0[idx_IB] + v_ip*norm_1[idx_IB];
                        // Real vel_ip_t = -u_ip*norm_1[idx_IB] + v_ip*norm_0[idx_IB];
                        
                        // Real vel_gc_n = vel_ip_n - ((d_ip + dist[idx_IB])/d_ip)*(vel_ip_n);
                        // Real vel_gc_t = vel_ip_t;
                        
                        // u_gc = vel_gc_n*norm_0[idx_IB] - vel_gc_t*norm_1[idx_IB];
                        // v_gc = vel_gc_n*norm_1[idx_IB] + vel_gc_t*norm_0[idx_IB];
                    }
                    else if (d_bc_type_velocity == VELOCITY_IBC::NO_SLIP)
                    {
                        u_gc = getGhostValueDirichletBC(
                            Real(0),
                            u_ip,
                            d_ip,
                            dist[idx_IB]);
                        
                        v_gc = getGhostValueDirichletBC(
                            Real(0),
                            v_ip,
                            d_ip,
                            dist[idx_IB]);
                    }
                    
                    // Bilinear interpolation to find specific internal energy for image point 1.
                    const Real epsilon_ip_BL = E[idx_ip_cons_var_BL]/rho[idx_ip_cons_var_BL] - half*(u_ip_BL*u_ip_BL + v_ip_BL*v_ip_BL);
                    const Real epsilon_ip_TL = E[idx_ip_cons_var_TL]/rho[idx_ip_cons_var_TL] - half*(u_ip_TL*u_ip_TL + v_ip_TL*v_ip_TL);
                    const Real epsilon_ip_BR = E[idx_ip_cons_var_BR]/rho[idx_ip_cons_var_BR] - half*(u_ip_BR*u_ip_BR + v_ip_BR*v_ip_BR);
                    const Real epsilon_ip_TR = E[idx_ip_cons_var_TR]/rho[idx_ip_cons_var_TR] - half*(u_ip_TR*u_ip_TR + v_ip_TR*v_ip_TR);
                    
                    // Bilinear interpolation to find pecific internal energy for image point 2.
                    const Real epsilon_ip2_BL = E[idx_ip2_cons_var_BL]/rho[idx_ip2_cons_var_BL] - half*(u_ip2_BL*u_ip2_BL + v_ip2_BL*v_ip2_BL);
                    const Real epsilon_ip2_TL = E[idx_ip2_cons_var_TL]/rho[idx_ip2_cons_var_TL] - half*(u_ip2_TL*u_ip2_TL + v_ip2_TL*v_ip2_TL);
                    const Real epsilon_ip2_BR = E[idx_ip2_cons_var_BR]/rho[idx_ip2_cons_var_BR] - half*(u_ip2_BR*u_ip2_BR + v_ip2_BR*v_ip2_BR);
                    const Real epsilon_ip2_TR = E[idx_ip2_cons_var_TR]/rho[idx_ip2_cons_var_TR] - half*(u_ip2_TR*u_ip2_TR + v_ip2_TR*v_ip2_TR);
                    
                    // Compute the pressure values in the stencils for image point 1.
                    const Real p_ip_BL = d_equation_of_state_mixing_rules->getEquationOfState()->getPressure(&rho[idx_ip_cons_var_BL], &epsilon_ip_BL, thermo_properties_ptr);
                    const Real p_ip_TL = d_equation_of_state_mixing_rules->getEquationOfState()->getPressure(&rho[idx_ip_cons_var_TL], &epsilon_ip_TL, thermo_properties_ptr);
                    const Real p_ip_BR = d_equation_of_state_mixing_rules->getEquationOfState()->getPressure(&rho[idx_ip_cons_var_BR], &epsilon_ip_BR, thermo_properties_ptr);
                    const Real p_ip_TR = d_equation_of_state_mixing_rules->getEquationOfState()->getPressure(&rho[idx_ip_cons_var_TR], &epsilon_ip_TR, thermo_properties_ptr);
                    
                    // Compute the pressure values in the stencils for image point 2.
                    const Real p_ip2_BL = d_equation_of_state_mixing_rules->getEquationOfState()->getPressure(&rho[idx_ip2_cons_var_BL], &epsilon_ip2_BL, thermo_properties_ptr);
                    const Real p_ip2_TL = d_equation_of_state_mixing_rules->getEquationOfState()->getPressure(&rho[idx_ip2_cons_var_TL], &epsilon_ip2_TL, thermo_properties_ptr);
                    const Real p_ip2_BR = d_equation_of_state_mixing_rules->getEquationOfState()->getPressure(&rho[idx_ip2_cons_var_BR], &epsilon_ip2_BR, thermo_properties_ptr);
                    const Real p_ip2_TR = d_equation_of_state_mixing_rules->getEquationOfState()->getPressure(&rho[idx_ip2_cons_var_TR], &epsilon_ip2_TR, thermo_properties_ptr);
                    
                    // Compute the temperature values in the stencils for image point 1.
                    const Real T_ip_BL = d_equation_of_state_mixing_rules->getEquationOfState()->getTemperature(&rho[idx_ip_cons_var_BL], &p_ip_BL, thermo_properties_ptr);
                    const Real T_ip_TL = d_equation_of_state_mixing_rules->getEquationOfState()->getTemperature(&rho[idx_ip_cons_var_TL], &p_ip_TL, thermo_properties_ptr);
                    const Real T_ip_BR = d_equation_of_state_mixing_rules->getEquationOfState()->getTemperature(&rho[idx_ip_cons_var_BR], &p_ip_BR, thermo_properties_ptr);
                    const Real T_ip_TR = d_equation_of_state_mixing_rules->getEquationOfState()->getTemperature(&rho[idx_ip_cons_var_TR], &p_ip_TR, thermo_properties_ptr);
                    
                    // Compute the temperature values in the stencils for image point 2.
                    const Real T_ip2_BL = d_equation_of_state_mixing_rules->getEquationOfState()->getTemperature(&rho[idx_ip2_cons_var_BL], &p_ip2_BL, thermo_properties_ptr);
                    const Real T_ip2_TL = d_equation_of_state_mixing_rules->getEquationOfState()->getTemperature(&rho[idx_ip2_cons_var_TL], &p_ip2_TL, thermo_properties_ptr);
                    const Real T_ip2_BR = d_equation_of_state_mixing_rules->getEquationOfState()->getTemperature(&rho[idx_ip2_cons_var_BR], &p_ip2_BR, thermo_properties_ptr);
                    const Real T_ip2_TR = d_equation_of_state_mixing_rules->getEquationOfState()->getTemperature(&rho[idx_ip2_cons_var_TR], &p_ip2_TR, thermo_properties_ptr);
                    
                    // Bilinear interpolation to find temperature of image point 1.
                    const Real T_ip = bilinearInterpolate2D(
                        T_ip_BL,
                        T_ip_BR,
                        T_ip_TL,
                        T_ip_TR,
                        x_ip,
                        y_ip,
                        x_ip_BL,
                        y_ip_BL,
                        dx_inv);
                    
                    // Bilinear interpolation to find temperature of image point 2.
                    const Real T_ip2 = bilinearInterpolate2D(
                        T_ip2_BL,
                        T_ip2_BR,
                        T_ip2_TL,
                        T_ip2_TR,
                        x_ip2,
                        y_ip2,
                        x_ip2_BL,
                        y_ip2_BL,
                        dx_inv);
                    
                    // Bilinear interpolation to find density of image point 1.
                    const Real rho_ip = bilinearInterpolate2D(
                        rho[idx_ip_cons_var_BL],
                        rho[idx_ip_cons_var_BR],
                        rho[idx_ip_cons_var_TL],
                        rho[idx_ip_cons_var_TR],
                        x_ip,
                        y_ip,
                        x_ip_BL,
                        y_ip_BL,
                        dx_inv);
                    
                    // Bilinear interpolation to find density of image point 2.
                    const Real rho_ip2 = bilinearInterpolate2D(
                        rho[idx_ip2_cons_var_BL],
                        rho[idx_ip2_cons_var_BR],
                        rho[idx_ip2_cons_var_TL],
                        rho[idx_ip2_cons_var_TR],
                        x_ip2,
                        y_ip2,
                        x_ip2_BL,
                        y_ip2_BL,
                        dx_inv);
                    
                    const Real epsilon_ip = d_equation_of_state_mixing_rules->getEquationOfState()->
                        getInternalEnergyFromTemperature(
                            &rho_ip,
                            &T_ip,
                            thermo_properties_ptr);
                        
                    const Real epsilon_ip2 = d_equation_of_state_mixing_rules->getEquationOfState()->
                        getInternalEnergyFromTemperature(
                            &rho_ip2,
                            &T_ip2,
                            thermo_properties_ptr);
                    
                    const Real p_ip = d_equation_of_state_mixing_rules->getEquationOfState()->
                        getPressure(
                            &rho_ip,
                            &epsilon_ip,
                            thermo_properties_ptr);
                    
                    const Real p_ip2 = d_equation_of_state_mixing_rules->getEquationOfState()->
                        getPressure(
                            &rho_ip2,
                            &epsilon_ip2,
                            thermo_properties_ptr);
                    
                    // dP/dn = 0
                    const Real p_gc = getGhostValueNeumannBC(
                            p_ip,
                            p_ip2,
                            d_ip,
                            d_ip2,
                            dist[idx_IB]);
                    
                    Real T_gc = Real(0);
                    if (d_bc_type_temperature == TEMPERATURE_IBC::ADIABATIC)
                    {
                        // dT/dn = 0
                        T_gc = getGhostValueNeumannBC(
                            T_ip,
                            T_ip2,
                            d_ip,
                            d_ip2,
                            dist[idx_IB]);
                    }
                    else if (d_bc_type_temperature == TEMPERATURE_IBC::ISOTHERMAL)
                    {
                        // Iso-thermal boundary condition.
                        T_gc = getGhostValueDirichletBC(
                            T_body,
                            T_ip,
                            d_ip,
                            dist[idx_IB]);
                    }
                    
                    // Compute density using calculated p and T.
                    const Real rho_gc = d_equation_of_state_mixing_rules->getEquationOfState()->
                        getDensity(
                            &p_gc,
                            &T_gc,
                            thermo_properties_ptr);
                    
                    // Compute the total energy at the ghost cell.
                    const Real epsilon_gc = d_equation_of_state_mixing_rules->getEquationOfState()->
                        getInternalEnergyFromTemperature(
                            &rho_gc,
                            &T_gc,
                            thermo_properties_ptr);
                    
                    const Real E_gc = rho_gc*(epsilon_gc + half*(u_gc*u_gc + v_gc*v_gc));
                    
                    rho[idx_cons_var]   = rho_gc;
                    rho_u[idx_cons_var] = rho[idx_cons_var]*u_gc;
                    rho_v[idx_cons_var] = rho[idx_cons_var]*v_gc;
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
        if (fabs(dx[0] - dx[1]) > 10.0*std::numeric_limits<double>::epsilon())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelImmersedBoundaryMethodSingleSpecies::setConservativeVariablesCellDataImmersedBoundaryGhosts()\n"
                << "The grid is assumed to be isotropic but the grid spacings are different."
                << std::endl);
        }
        if (fabs(dx[0] - dx[2]) > 10.0*std::numeric_limits<double>::epsilon())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelImmersedBoundaryMethodSingleSpecies::setConservativeVariablesCellDataImmersedBoundaryGhosts()\n"
                << "The grid is assumed to be isotropic but the grid spacings are different."
                << std::endl);
        }
        
        const Real dx_inv = Real(1)/Real(dx[0]);
        
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
        
        Real* norm_0 = data_surface_normal->getPointer(0);
        Real* norm_1 = data_surface_normal->getPointer(1);
        Real* norm_2 = data_surface_normal->getPointer(2);
        
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
                    
                    // Compute the coordinates.
                    const Real x[3] = {
                        Real(patch_xlo[0]) + (Real(i) + half)*Real(dx[0]),
                        Real(patch_xlo[1]) + (Real(j) + half)*Real(dx[1]),
                        Real(patch_xlo[2]) + (Real(k) + half)*Real(dx[2])
                        };
                    
                    if (mask[idx_IB] == int(IB_MASK::IB_GHOST) || mask[idx_IB] == int(IB_MASK::IB_GHOST_CORNER))
                    {
                        
                        Real d_ip;
                        // First image point distance is set to sqrt(3)*dx + epsilon for corner ghost cells.
                        if (mask[idx_IB] == int(IB_MASK::IB_GHOST_CORNER))
                        {
                            d_ip   = std::sqrt(Real(3))*Real(dx[0]) + HAMERS_REAL_EPSILON;
                        }
                        // First image point distance is set to dx/maximum(norm) for ghost cells for convective fluxes.
                        else
                        {
                            Real norm_max = std::max(std::max(std::abs(norm_0[idx_IB]), std::abs(norm_1[idx_IB])), std::abs(norm_2[idx_IB]));
                            d_ip  = Real(dx[0]) / norm_max + HAMERS_REAL_EPSILON;
                        }
                        
                        // Second image point distance is set to d_ip + 0.25*dx.
                        Real d_ip2 = d_ip + Real(0.25)*Real(dx[0]);
                        // Coordinates of the image point 1.
                        const Real x_ip = x[0] + (dist[idx_IB] + d_ip)*norm_0[idx_IB];
                        const Real y_ip = x[1] + (dist[idx_IB] + d_ip)*norm_1[idx_IB];
                        const Real z_ip = x[2] + (dist[idx_IB] + d_ip)*norm_2[idx_IB];
                        
                        // Coordinates of the image point 2.
                        const Real x_ip2 = x[0] + (dist[idx_IB] + d_ip2)*norm_0[idx_IB];
                        const Real y_ip2 = x[1] + (dist[idx_IB] + d_ip2)*norm_1[idx_IB];
                        const Real z_ip2 = x[2] + (dist[idx_IB] + d_ip2)*norm_2[idx_IB];
                        
                        // Get indices of the cells in interpolation for image point 1.
                        int idx_ip_cons_var_LBK, idx_ip_cons_var_RBK, idx_ip_cons_var_LTK, idx_ip_cons_var_RTK,
                            idx_ip_cons_var_LBF, idx_ip_cons_var_RBF, idx_ip_cons_var_LTF, idx_ip_cons_var_RTF;
                        int idx_ip_IB_LBK, idx_ip_IB_RBK, idx_ip_IB_LTK, idx_ip_IB_RTK,
                            idx_ip_IB_LBF, idx_ip_IB_RBF, idx_ip_IB_LTF, idx_ip_IB_RTF;
                        Real x_ip_LBK, y_ip_LBK, z_ip_LBK;
                        getTrilinearInterpolationIndices3D(
                            idx_ip_cons_var_LBK,
                            idx_ip_cons_var_RBK,
                            idx_ip_cons_var_LTK,
                            idx_ip_cons_var_RTK,
                            idx_ip_cons_var_LBF,
                            idx_ip_cons_var_RBF,
                            idx_ip_cons_var_LTF,
                            idx_ip_cons_var_RTF,
                            idx_ip_IB_LBK,
                            idx_ip_IB_RBK,
                            idx_ip_IB_LTK,
                            idx_ip_IB_RTK,
                            idx_ip_IB_LBF,
                            idx_ip_IB_RBF,
                            idx_ip_IB_LTF,
                            idx_ip_IB_RTF,
                            x_ip_LBK,
                            y_ip_LBK,
                            z_ip_LBK,
                            x_ip,
                            y_ip,
                            z_ip,
                            Real(patch_xlo[0]),
                            Real(patch_xlo[1]),
                            Real(patch_xlo[2]),
                            offset_0_cons_var,
                            offset_1_cons_var,
                            offset_2_cons_var,
                            ghostcell_dim_0_cons_var,
                            ghostcell_dim_1_cons_var,
                            offset_0_IB,
                            offset_1_IB,
                            offset_2_IB,
                            ghostcell_dim_0_IB,
                            ghostcell_dim_1_IB,
                            Real(dx[0]),
                            dx_inv);
                        
                        // Checking first image point interpolation stencil to ensure only fluid cell values are used
                        if (mask[idx_ip_IB_LBK] != int(IB_MASK::FLUID))
                        {
                            TBOX_ERROR("Error: Left-bottom-back cell is not FLUID at index: " << idx_ip_IB_LBK 
                                << "\n with mask value: " << mask[idx_ip_IB_LBK]
                                << "\n x_ip: " << x_ip
                                << "\n y_ip: " << y_ip
                                << "\n x: " << x[0]
                                << "\n y: " << x[1]
                                << "\n d_ip: " << d_ip/Real(dx[0])
                                << "\n d_gc: " << dist[idx_IB]
                                << "\n norm_0: " << norm_0[idx_IB]
                                << "\n norm_1: " << norm_1[idx_IB]
                                << "\n dx: " << Real(dx[0])
                                << "\n dx_inv: " << dx_inv
                                << "\n x_ip_LBK: " << x_ip_LBK
                                << "\n y_ip_LBK: " << y_ip_LBK);
                        }
                        
                        if (mask[idx_ip_IB_RBK] != int(IB_MASK::FLUID))
                        {
                            TBOX_ERROR("Error: Right-bottom-back cell is not FLUID at index: " << idx_ip_IB_RBK 
                                << "\n with mask value: " << mask[idx_ip_IB_RBK]
                                << "\n x_ip: " << x_ip
                                << "\n y_ip: " << y_ip
                                << "\n x: " << x[0]
                                << "\n y: " << x[1]
                                << "\n d_ip: " << d_ip/Real(dx[0])
                                << "\n d_gc: " << dist[idx_IB]
                                << "\n norm_0: " << norm_0[idx_IB]
                                << "\n norm_1: " << norm_1[idx_IB]
                                << "\n dx: " << Real(dx[0])
                                << "\n dx_inv: " << dx_inv
                                << "\n x_ip_LBK: " << x_ip_LBK
                                << "\n y_ip_LBK: " << y_ip_LBK);
                        }
                        
                        if (mask[idx_ip_IB_LTK] != int(IB_MASK::FLUID))
                        {
                            TBOX_ERROR("Error: Left-top-back cell is not FLUID at index: " << idx_ip_IB_LTK 
                                << "\n with mask value: " << mask[idx_ip_IB_LTK]
                                << "\n x_ip: " << x_ip
                                << "\n y_ip: " << y_ip
                                << "\n x: " << x[0]
                                << "\n y: " << x[1]
                                << "\n d_ip: " << d_ip/Real(dx[0])
                                << "\n d_gc: " << dist[idx_IB]
                                << "\n norm_0: " << norm_0[idx_IB]
                                << "\n norm_1: " << norm_1[idx_IB]
                                << "\n dx: " << Real(dx[0])
                                << "\n dx_inv: " << dx_inv
                                << "\n x_ip_LBK: " << x_ip_LBK
                                << "\n y_ip_LBK: " << y_ip_LBK);
                        }
                        
                        if (mask[idx_ip_IB_RTK] != int(IB_MASK::FLUID))
                        {
                            TBOX_ERROR("Error: Right-top-back cell is not FLUID at index: " << idx_ip_IB_RTK 
                                << "\n with mask value: " << mask[idx_ip_IB_RTK]
                                << "\n x_ip: " << x_ip
                                << "\n y_ip: " << y_ip
                                << "\n x: " << x[0]
                                << "\n y: " << x[1]
                                << "\n d_ip: " << d_ip/Real(dx[0])
                                << "\n d_gc: " << dist[idx_IB]
                                << "\n norm_0: " << norm_0[idx_IB]
                                << "\n norm_1: " << norm_1[idx_IB]
                                << "\n dx: " << Real(dx[0])
                                << "\n dx_inv: " << dx_inv
                                << "\n x_ip_LBK: " << x_ip_LBK
                                << "\n y_ip_LBK: " << y_ip_LBK);
                        }
                        
                        // Checking first image point interpolation stencil to ensure only fluid cell values are used
                        if (mask[idx_ip_IB_LBF] != int(IB_MASK::FLUID))
                        {
                            TBOX_ERROR("Error: Left-bottom-front cell is not FLUID at index: " << idx_ip_IB_LBF 
                                << "\n with mask value: " << mask[idx_ip_IB_LBF]
                                << "\n x_ip: " << x_ip
                                << "\n y_ip: " << y_ip
                                << "\n x: " << x[0]
                                << "\n y: " << x[1]
                                << "\n d_ip: " << d_ip/Real(dx[0])
                                << "\n d_gc: " << dist[idx_IB]
                                << "\n norm_0: " << norm_0[idx_IB]
                                << "\n norm_1: " << norm_1[idx_IB]
                                << "\n dx: " << Real(dx[0])
                                << "\n dx_inv: " << dx_inv
                                << "\n x_ip_LBK: " << x_ip_LBK
                                << "\n y_ip_LBK: " << y_ip_LBK);
                        }
                        
                        if (mask[idx_ip_IB_RBF] != int(IB_MASK::FLUID))
                        {
                            TBOX_ERROR("Error: Right-bottom-front cell is not FLUID at index: " << idx_ip_IB_RBF 
                                << "\n with mask value: " << mask[idx_ip_IB_RBF]
                                << "\n x_ip: " << x_ip
                                << "\n y_ip: " << y_ip
                                << "\n x: " << x[0]
                                << "\n y: " << x[1]
                                << "\n d_ip: " << d_ip/Real(dx[0])
                                << "\n d_gc: " << dist[idx_IB]
                                << "\n norm_0: " << norm_0[idx_IB]
                                << "\n norm_1: " << norm_1[idx_IB]
                                << "\n dx: " << Real(dx[0])
                                << "\n dx_inv: " << dx_inv
                                << "\n x_ip_LBK: " << x_ip_LBK
                                << "\n y_ip_LBK: " << y_ip_LBK);
                        }
                        
                        if (mask[idx_ip_IB_LTF] != int(IB_MASK::FLUID))
                        {
                            TBOX_ERROR("Error: Left-top-front cell is not FLUID at index: " << idx_ip_IB_LTF 
                                << "\n with mask value: " << mask[idx_ip_IB_LTF]
                                << "\n x_ip: " << x_ip
                                << "\n y_ip: " << y_ip
                                << "\n x: " << x[0]
                                << "\n y: " << x[1]
                                << "\n d_ip: " << d_ip/Real(dx[0])
                                << "\n d_gc: " << dist[idx_IB]
                                << "\n norm_0: " << norm_0[idx_IB]
                                << "\n norm_1: " << norm_1[idx_IB]
                                << "\n dx: " << Real(dx[0])
                                << "\n dx_inv: " << dx_inv
                                << "\n x_ip_LBK: " << x_ip_LBK
                                << "\n y_ip_LBK: " << y_ip_LBK);
                        }
                        
                        if (mask[idx_ip_IB_RTF] != int(IB_MASK::FLUID))
                        {
                            TBOX_ERROR("Error: Right-top-front cell is not FLUID at index: " << idx_ip_IB_RTF 
                                << "\n with mask value: " << mask[idx_ip_IB_RTF]
                                << "\n x_ip: " << x_ip
                                << "\n y_ip: " << y_ip
                                << "\n x: " << x[0]
                                << "\n y: " << x[1]
                                << "\n d_ip: " << d_ip/Real(dx[0])
                                << "\n d_gc: " << dist[idx_IB]
                                << "\n norm_0: " << norm_0[idx_IB]
                                << "\n norm_1: " << norm_1[idx_IB]
                                << "\n dx: " << Real(dx[0])
                                << "\n dx_inv: " << dx_inv
                                << "\n x_ip_LBK: " << x_ip_LBK
                                << "\n y_ip_LBK: " << y_ip_LBK);
                        }

                        // Get indices of the cells in interpolation for image point 2.
                        int idx_ip2_cons_var_LBK, idx_ip2_cons_var_RBK, idx_ip2_cons_var_LTK, idx_ip2_cons_var_RTK,
                            idx_ip2_cons_var_LBF, idx_ip2_cons_var_RBF, idx_ip2_cons_var_LTF, idx_ip2_cons_var_RTF;
                        int idx_ip2_IB_LBK, idx_ip2_IB_RBK, idx_ip2_IB_LTK, idx_ip2_IB_RTK,
                            idx_ip2_IB_LBF, idx_ip2_IB_RBF, idx_ip2_IB_LTF, idx_ip2_IB_RTF;
                        Real x_ip2_LBK, y_ip2_LBK, z_ip2_LBK;
                        getTrilinearInterpolationIndices3D(
                            idx_ip2_cons_var_LBK,
                            idx_ip2_cons_var_RBK,
                            idx_ip2_cons_var_LTK,
                            idx_ip2_cons_var_RTK,
                            idx_ip2_cons_var_LBF,
                            idx_ip2_cons_var_RBF,
                            idx_ip2_cons_var_LTF,
                            idx_ip2_cons_var_RTF,
                            idx_ip2_IB_LBK,
                            idx_ip2_IB_RBK,
                            idx_ip2_IB_LTK,
                            idx_ip2_IB_RTK,
                            idx_ip2_IB_LBF,
                            idx_ip2_IB_RBF,
                            idx_ip2_IB_LTF,
                            idx_ip2_IB_RTF,
                            x_ip2_LBK,
                            y_ip2_LBK,
                            z_ip2_LBK,
                            x_ip2,
                            y_ip2,
                            z_ip2,
                            Real(patch_xlo[0]),
                            Real(patch_xlo[1]),
                            Real(patch_xlo[2]),
                            offset_0_cons_var,
                            offset_1_cons_var,
                            offset_2_cons_var,
                            ghostcell_dim_0_cons_var,
                            ghostcell_dim_1_cons_var,
                            offset_0_IB,
                            offset_1_IB,
                            offset_2_IB,
                            ghostcell_dim_0_IB,
                            ghostcell_dim_1_IB,
                            Real(dx[0]),
                            dx_inv);
                        
                        // Bilinear interpolation to find image point x-component of velocity values.
                        const Real u_ip_LBK = rho_u[idx_ip_cons_var_LBK]/rho[idx_ip_cons_var_LBK];
                        const Real u_ip_RBK = rho_u[idx_ip_cons_var_RBK]/rho[idx_ip_cons_var_RBK];
                        const Real u_ip_LTK = rho_u[idx_ip_cons_var_LTK]/rho[idx_ip_cons_var_LTK];
                        const Real u_ip_RTK = rho_u[idx_ip_cons_var_RTK]/rho[idx_ip_cons_var_RTK];
                        const Real u_ip_LBF = rho_u[idx_ip_cons_var_LBF]/rho[idx_ip_cons_var_LBF];
                        const Real u_ip_RBF = rho_u[idx_ip_cons_var_RBF]/rho[idx_ip_cons_var_RBF];
                        const Real u_ip_LTF = rho_u[idx_ip_cons_var_LTF]/rho[idx_ip_cons_var_LTF];
                        const Real u_ip_RTF = rho_u[idx_ip_cons_var_RTF]/rho[idx_ip_cons_var_RTF];
                        
                        const Real u_ip2_LBK = rho_u[idx_ip2_cons_var_LBK]/rho[idx_ip2_cons_var_LBK];
                        const Real u_ip2_RBK = rho_u[idx_ip2_cons_var_RBK]/rho[idx_ip2_cons_var_RBK];
                        const Real u_ip2_LTK = rho_u[idx_ip2_cons_var_LTK]/rho[idx_ip2_cons_var_LTK];
                        const Real u_ip2_RTK = rho_u[idx_ip2_cons_var_RTK]/rho[idx_ip2_cons_var_RTK];
                        const Real u_ip2_LBF = rho_u[idx_ip2_cons_var_LBF]/rho[idx_ip2_cons_var_LBF];
                        const Real u_ip2_RBF = rho_u[idx_ip2_cons_var_RBF]/rho[idx_ip2_cons_var_RBF];
                        const Real u_ip2_LTF = rho_u[idx_ip2_cons_var_LTF]/rho[idx_ip2_cons_var_LTF];
                        const Real u_ip2_RTF = rho_u[idx_ip2_cons_var_RTF]/rho[idx_ip2_cons_var_RTF];
                        
                        const Real u_ip = trilinearInterpolate3D(
                            u_ip_LBK,
                            u_ip_RBK,
                            u_ip_LTK,
                            u_ip_RTK,
                            u_ip_LBF,
                            u_ip_RBF,
                            u_ip_LTF,
                            u_ip_RTF,
                            x_ip,
                            y_ip,
                            z_ip,
                            x_ip_LBK,
                            y_ip_LBK,
                            z_ip_LBK,
                            dx_inv);
                        
                        // Bilinear interpolation to find image point y-component of velocity values.
                        const Real v_ip_LBK = rho_v[idx_ip_cons_var_LBK]/rho[idx_ip_cons_var_LBK];
                        const Real v_ip_RBK = rho_v[idx_ip_cons_var_RBK]/rho[idx_ip_cons_var_RBK];
                        const Real v_ip_LTK = rho_v[idx_ip_cons_var_LTK]/rho[idx_ip_cons_var_LTK];
                        const Real v_ip_RTK = rho_v[idx_ip_cons_var_RTK]/rho[idx_ip_cons_var_RTK];
                        const Real v_ip_LBF = rho_v[idx_ip_cons_var_LBF]/rho[idx_ip_cons_var_LBF];
                        const Real v_ip_RBF = rho_v[idx_ip_cons_var_RBF]/rho[idx_ip_cons_var_RBF];
                        const Real v_ip_LTF = rho_v[idx_ip_cons_var_LTF]/rho[idx_ip_cons_var_LTF];
                        const Real v_ip_RTF = rho_v[idx_ip_cons_var_RTF]/rho[idx_ip_cons_var_RTF];
                        
                        const Real v_ip2_LBK = rho_v[idx_ip2_cons_var_LBK]/rho[idx_ip2_cons_var_LBK];
                        const Real v_ip2_RBK = rho_v[idx_ip2_cons_var_RBK]/rho[idx_ip2_cons_var_RBK];
                        const Real v_ip2_LTK = rho_v[idx_ip2_cons_var_LTK]/rho[idx_ip2_cons_var_LTK];
                        const Real v_ip2_RTK = rho_v[idx_ip2_cons_var_RTK]/rho[idx_ip2_cons_var_RTK];
                        const Real v_ip2_LBF = rho_v[idx_ip2_cons_var_LBF]/rho[idx_ip2_cons_var_LBF];
                        const Real v_ip2_RBF = rho_v[idx_ip2_cons_var_RBF]/rho[idx_ip2_cons_var_RBF];
                        const Real v_ip2_LTF = rho_v[idx_ip2_cons_var_LTF]/rho[idx_ip2_cons_var_LTF];
                        const Real v_ip2_RTF = rho_v[idx_ip2_cons_var_RTF]/rho[idx_ip2_cons_var_RTF];
                        
                        const Real v_ip = trilinearInterpolate3D(
                            v_ip_LBK,
                            v_ip_RBK,
                            v_ip_LTK,
                            v_ip_RTK,
                            v_ip_LBF,
                            v_ip_RBF,
                            v_ip_LTF,
                            v_ip_RTF,
                            x_ip,
                            y_ip,
                            z_ip,
                            x_ip_LBK,
                            y_ip_LBK,
                            z_ip_LBK,
                            dx_inv);
                        
                        // Bilinear interpolation to find image point z-component of velocity values.
                        const Real w_ip_LBK = rho_w[idx_ip_cons_var_LBK]/rho[idx_ip_cons_var_LBK];
                        const Real w_ip_RBK = rho_w[idx_ip_cons_var_RBK]/rho[idx_ip_cons_var_RBK];
                        const Real w_ip_LTK = rho_w[idx_ip_cons_var_LTK]/rho[idx_ip_cons_var_LTK];
                        const Real w_ip_RTK = rho_w[idx_ip_cons_var_RTK]/rho[idx_ip_cons_var_RTK];
                        const Real w_ip_LBF = rho_w[idx_ip_cons_var_LBF]/rho[idx_ip_cons_var_LBF];
                        const Real w_ip_RBF = rho_w[idx_ip_cons_var_RBF]/rho[idx_ip_cons_var_RBF];
                        const Real w_ip_LTF = rho_w[idx_ip_cons_var_LTF]/rho[idx_ip_cons_var_LTF];
                        const Real w_ip_RTF = rho_w[idx_ip_cons_var_RTF]/rho[idx_ip_cons_var_RTF];
                        
                        const Real w_ip2_LBK = rho_w[idx_ip2_cons_var_LBK]/rho[idx_ip2_cons_var_LBK];
                        const Real w_ip2_RBK = rho_w[idx_ip2_cons_var_RBK]/rho[idx_ip2_cons_var_RBK];
                        const Real w_ip2_LTK = rho_w[idx_ip2_cons_var_LTK]/rho[idx_ip2_cons_var_LTK];
                        const Real w_ip2_RTK = rho_w[idx_ip2_cons_var_RTK]/rho[idx_ip2_cons_var_RTK];
                        const Real w_ip2_LBF = rho_w[idx_ip2_cons_var_LBF]/rho[idx_ip2_cons_var_LBF];
                        const Real w_ip2_RBF = rho_w[idx_ip2_cons_var_RBF]/rho[idx_ip2_cons_var_RBF];
                        const Real w_ip2_LTF = rho_w[idx_ip2_cons_var_LTF]/rho[idx_ip2_cons_var_LTF];
                        const Real w_ip2_RTF = rho_w[idx_ip2_cons_var_RTF]/rho[idx_ip2_cons_var_RTF];
                        
                        const Real w_ip = trilinearInterpolate3D(
                            w_ip_LBK,
                            w_ip_RBK,
                            w_ip_LTK,
                            w_ip_RTK,
                            w_ip_LBF,
                            w_ip_RBF,
                            w_ip_LTF,
                            w_ip_RTF,
                            x_ip,
                            y_ip,
                            z_ip,
                            x_ip_LBK,
                            y_ip_LBK,
                            z_ip_LBK,
                            dx_inv);
                        
                        // Compute velocity components at the ghost cell depending on the type of immersed boundary condition.
                        Real u_gc = Real(0); // x-component of velocity of the ghost cell
                        Real v_gc = Real(0); // y-component of velocity of the ghost cell
                        Real w_gc = Real(0); // z-component of velocity of the ghost cell
                        
                        if (d_bc_type_velocity == VELOCITY_IBC::SLIP)
                        {
                            const Real u_ip2 = trilinearInterpolate3D(
                                u_ip2_LBK,
                                u_ip2_RBK,
                                u_ip2_LTK,
                                u_ip2_RTK,
                                u_ip2_LBF,
                                u_ip2_RBF,
                                u_ip2_LTF,
                                u_ip2_RTF,
                                x_ip2,
                                y_ip2,
                                z_ip2,
                                x_ip2_LBK,
                                y_ip2_LBK,
                                z_ip2_LBK,
                                dx_inv);
                            
                            const Real v_ip2 = trilinearInterpolate3D(
                                v_ip2_LBK,
                                v_ip2_RBK,
                                v_ip2_LTK,
                                v_ip2_RTK,
                                v_ip2_LBF,
                                v_ip2_RBF,
                                v_ip2_LTF,
                                v_ip2_RTF,
                                x_ip2,
                                y_ip2,
                                z_ip2,
                                x_ip2_LBK,
                                y_ip2_LBK,
                                z_ip2_LBK,
                                dx_inv);
                            
                            const Real w_ip2 = trilinearInterpolate3D(
                                w_ip2_LBK,
                                w_ip2_RBK,
                                w_ip2_LTK,
                                w_ip2_RTK,
                                w_ip2_LBF,
                                w_ip2_RBF,
                                w_ip2_LTF,
                                w_ip2_RTF,
                                x_ip2,
                                y_ip2,
                                z_ip2,
                                x_ip2_LBK,
                                y_ip2_LBK,
                                z_ip2_LBK,
                                dx_inv);
                            
                            // Given d_ip and d_ip2, interpolate to the velocity components (u_mirror and v_mirror)
                            // at the mirror image point at dist[idx_IB].
                            const Real diff_ip2_ip = d_ip2 - d_ip;
                            const Real diff_mirror_ip  = dist[idx_IB] - d_ip;
                            const Real diff_ip2_mirror = d_ip2 - dist[idx_IB];
                            
                            // x-component of velocity at the mirror image point.
                            const Real u_mirror = (diff_ip2_mirror*u_ip + diff_mirror_ip*u_ip2)/diff_ip2_ip;
                            // y-component of velocity at the mirror image point.
                            const Real v_mirror = (diff_ip2_mirror*v_ip + diff_mirror_ip*v_ip2)/diff_ip2_ip;
                            // z-component of velocity at the mirror image point.
                            const Real w_mirror = (diff_ip2_mirror*w_ip + diff_mirror_ip*w_ip2)/diff_ip2_ip;
                            
                            // Velocity component normal to the boundary at the mirror image point.
                            const Real vel_mirror_n = dotProduct3D(u_mirror, v_mirror, w_mirror, norm_0[idx_IB], norm_1[idx_IB], norm_2[idx_IB]);
                            
                            // No-penetration boundary condition.
                            u_gc = u_mirror - Real(2)*vel_mirror_n*norm_0[idx_IB];
                            v_gc = v_mirror - Real(2)*vel_mirror_n*norm_1[idx_IB];
                            w_gc = w_mirror - Real(2)*vel_mirror_n*norm_2[idx_IB];
                        }
                        else if (d_bc_type_velocity == VELOCITY_IBC::NO_SLIP)
                        {
                            u_gc = getGhostValueDirichletBC(
                                Real(0),
                                u_ip,
                                d_ip,
                                dist[idx_IB]);
                            
                            v_gc = getGhostValueDirichletBC(
                                Real(0),
                                v_ip,
                                d_ip,
                                dist[idx_IB]);
                            
                            w_gc = getGhostValueDirichletBC(
                                Real(0),
                                w_ip,
                                d_ip,
                                dist[idx_IB]);
                            
                            if (fabs(u_gc) <= HAMERS_REAL_EPSILON)
                            {
                                u_gc = Real(0);
                            }
                            if (fabs(v_gc) <= HAMERS_REAL_EPSILON)
                            {
                                v_gc = Real(0);
                            }
                            if (fabs(w_gc) <= HAMERS_REAL_EPSILON)
                            {
                                w_gc = Real(0);
                            }
                        }
                        
                        // Bilinear interpolation to find specific internal energy for image point 1.
                        const Real epsilon_ip_LBK = E[idx_ip_cons_var_LBK]/rho[idx_ip_cons_var_LBK] - half*(u_ip_LBK*u_ip_LBK + v_ip_LBK*v_ip_LBK + w_ip_LBK*w_ip_LBK);
                        const Real epsilon_ip_RBK = E[idx_ip_cons_var_RBK]/rho[idx_ip_cons_var_RBK] - half*(u_ip_RBK*u_ip_RBK + v_ip_RBK*v_ip_RBK + w_ip_RBK*w_ip_RBK);
                        const Real epsilon_ip_LTK = E[idx_ip_cons_var_LTK]/rho[idx_ip_cons_var_LTK] - half*(u_ip_LTK*u_ip_LTK + v_ip_LTK*v_ip_LTK + w_ip_LTK*w_ip_LTK);
                        const Real epsilon_ip_RTK = E[idx_ip_cons_var_RTK]/rho[idx_ip_cons_var_RTK] - half*(u_ip_RTK*u_ip_RTK + v_ip_RTK*v_ip_RTK + w_ip_RTK*w_ip_RTK);
                        const Real epsilon_ip_LBF = E[idx_ip_cons_var_LBF]/rho[idx_ip_cons_var_LBF] - half*(u_ip_LBF*u_ip_LBF + v_ip_LBF*v_ip_LBF + w_ip_LBF*w_ip_LBF);
                        const Real epsilon_ip_RBF = E[idx_ip_cons_var_RBF]/rho[idx_ip_cons_var_RBF] - half*(u_ip_RBF*u_ip_RBF + v_ip_RBF*v_ip_RBF + w_ip_RBF*w_ip_RBF);
                        const Real epsilon_ip_LTF = E[idx_ip_cons_var_LTF]/rho[idx_ip_cons_var_LTF] - half*(u_ip_LTF*u_ip_LTF + v_ip_LTF*v_ip_LTF + w_ip_LTF*w_ip_LTF);
                        const Real epsilon_ip_RTF = E[idx_ip_cons_var_RTF]/rho[idx_ip_cons_var_RTF] - half*(u_ip_RTF*u_ip_RTF + v_ip_RTF*v_ip_RTF + w_ip_RTF*w_ip_RTF);
                        
                        // Bilinear interpolation to find specific internal energy for image point 2.
                        const Real epsilon_ip2_LBK = E[idx_ip2_cons_var_LBK]/rho[idx_ip2_cons_var_LBK] - half*(u_ip2_LBK*u_ip2_LBK + v_ip2_LBK*v_ip2_LBK + w_ip2_LBK*w_ip2_LBK);
                        const Real epsilon_ip2_RBK = E[idx_ip2_cons_var_RBK]/rho[idx_ip2_cons_var_RBK] - half*(u_ip2_RBK*u_ip2_RBK + v_ip2_RBK*v_ip2_RBK + w_ip2_RBK*w_ip2_RBK);
                        const Real epsilon_ip2_LTK = E[idx_ip2_cons_var_LTK]/rho[idx_ip2_cons_var_LTK] - half*(u_ip2_LTK*u_ip2_LTK + v_ip2_LTK*v_ip2_LTK + w_ip2_LTK*w_ip2_LTK);
                        const Real epsilon_ip2_RTK = E[idx_ip2_cons_var_RTK]/rho[idx_ip2_cons_var_RTK] - half*(u_ip2_RTK*u_ip2_RTK + v_ip2_RTK*v_ip2_RTK + w_ip2_RTK*w_ip2_RTK);
                        const Real epsilon_ip2_LBF = E[idx_ip2_cons_var_LBF]/rho[idx_ip2_cons_var_LBF] - half*(u_ip2_LBF*u_ip2_LBF + v_ip2_LBF*v_ip2_LBF + w_ip2_LBF*w_ip2_LBF);
                        const Real epsilon_ip2_RBF = E[idx_ip2_cons_var_RBF]/rho[idx_ip2_cons_var_RBF] - half*(u_ip2_RBF*u_ip2_RBF + v_ip2_RBF*v_ip2_RBF + w_ip2_RBF*w_ip2_RBF);
                        const Real epsilon_ip2_LTF = E[idx_ip2_cons_var_LTF]/rho[idx_ip2_cons_var_LTF] - half*(u_ip2_LTF*u_ip2_LTF + v_ip2_LTF*v_ip2_LTF + w_ip2_LTF*w_ip2_LTF);
                        const Real epsilon_ip2_RTF = E[idx_ip2_cons_var_RTF]/rho[idx_ip2_cons_var_RTF] - half*(u_ip2_RTF*u_ip2_RTF + v_ip2_RTF*v_ip2_RTF + w_ip2_RTF*w_ip2_RTF);
                        
                        // Compute the pressure values in the stencils for image point 1.
                        const Real p_ip_LBK = d_equation_of_state_mixing_rules->getEquationOfState()->getPressure(&rho[idx_ip_cons_var_LBK], &epsilon_ip_LBK, thermo_properties_ptr);
                        const Real p_ip_RBK = d_equation_of_state_mixing_rules->getEquationOfState()->getPressure(&rho[idx_ip_cons_var_RBK], &epsilon_ip_RBK, thermo_properties_ptr);
                        const Real p_ip_LTK = d_equation_of_state_mixing_rules->getEquationOfState()->getPressure(&rho[idx_ip_cons_var_LTK], &epsilon_ip_LTK, thermo_properties_ptr);
                        const Real p_ip_RTK = d_equation_of_state_mixing_rules->getEquationOfState()->getPressure(&rho[idx_ip_cons_var_RTK], &epsilon_ip_RTK, thermo_properties_ptr);
                        const Real p_ip_LBF = d_equation_of_state_mixing_rules->getEquationOfState()->getPressure(&rho[idx_ip_cons_var_LBF], &epsilon_ip_LBF, thermo_properties_ptr);
                        const Real p_ip_RBF = d_equation_of_state_mixing_rules->getEquationOfState()->getPressure(&rho[idx_ip_cons_var_RBF], &epsilon_ip_RBF, thermo_properties_ptr);
                        const Real p_ip_LTF = d_equation_of_state_mixing_rules->getEquationOfState()->getPressure(&rho[idx_ip_cons_var_LTF], &epsilon_ip_LTF, thermo_properties_ptr);
                        const Real p_ip_RTF = d_equation_of_state_mixing_rules->getEquationOfState()->getPressure(&rho[idx_ip_cons_var_RTF], &epsilon_ip_RTF, thermo_properties_ptr);
                        
                        // Compute the pressure values in the stencils for image point 2.
                        const Real p_ip2_LBK = d_equation_of_state_mixing_rules->getEquationOfState()->getPressure(&rho[idx_ip2_cons_var_LBK], &epsilon_ip2_LBK, thermo_properties_ptr);
                        const Real p_ip2_RBK = d_equation_of_state_mixing_rules->getEquationOfState()->getPressure(&rho[idx_ip2_cons_var_RBK], &epsilon_ip2_RBK, thermo_properties_ptr);
                        const Real p_ip2_LTK = d_equation_of_state_mixing_rules->getEquationOfState()->getPressure(&rho[idx_ip2_cons_var_LTK], &epsilon_ip2_LTK, thermo_properties_ptr);
                        const Real p_ip2_RTK = d_equation_of_state_mixing_rules->getEquationOfState()->getPressure(&rho[idx_ip2_cons_var_RTK], &epsilon_ip2_RTK, thermo_properties_ptr);
                        const Real p_ip2_LBF = d_equation_of_state_mixing_rules->getEquationOfState()->getPressure(&rho[idx_ip2_cons_var_LBF], &epsilon_ip2_LBF, thermo_properties_ptr);
                        const Real p_ip2_RBF = d_equation_of_state_mixing_rules->getEquationOfState()->getPressure(&rho[idx_ip2_cons_var_RBF], &epsilon_ip2_RBF, thermo_properties_ptr);
                        const Real p_ip2_LTF = d_equation_of_state_mixing_rules->getEquationOfState()->getPressure(&rho[idx_ip2_cons_var_LTF], &epsilon_ip2_LTF, thermo_properties_ptr);
                        const Real p_ip2_RTF = d_equation_of_state_mixing_rules->getEquationOfState()->getPressure(&rho[idx_ip2_cons_var_RTF], &epsilon_ip2_RTF, thermo_properties_ptr);
                        
                        // Compute the temperature values in the stencils for image point 1.
                        const Real T_ip_LBK = d_equation_of_state_mixing_rules->getEquationOfState()->getTemperature(&rho[idx_ip_cons_var_LBK], &p_ip_LBK, thermo_properties_ptr);
                        const Real T_ip_RBK = d_equation_of_state_mixing_rules->getEquationOfState()->getTemperature(&rho[idx_ip_cons_var_RBK], &p_ip_RBK, thermo_properties_ptr);
                        const Real T_ip_LTK = d_equation_of_state_mixing_rules->getEquationOfState()->getTemperature(&rho[idx_ip_cons_var_LTK], &p_ip_LTK, thermo_properties_ptr);
                        const Real T_ip_RTK = d_equation_of_state_mixing_rules->getEquationOfState()->getTemperature(&rho[idx_ip_cons_var_RTK], &p_ip_RTK, thermo_properties_ptr);
                        const Real T_ip_LBF = d_equation_of_state_mixing_rules->getEquationOfState()->getTemperature(&rho[idx_ip_cons_var_LBF], &p_ip_LBF, thermo_properties_ptr);
                        const Real T_ip_RBF = d_equation_of_state_mixing_rules->getEquationOfState()->getTemperature(&rho[idx_ip_cons_var_RBF], &p_ip_RBF, thermo_properties_ptr);
                        const Real T_ip_LTF = d_equation_of_state_mixing_rules->getEquationOfState()->getTemperature(&rho[idx_ip_cons_var_LTF], &p_ip_LTF, thermo_properties_ptr);
                        const Real T_ip_RTF = d_equation_of_state_mixing_rules->getEquationOfState()->getTemperature(&rho[idx_ip_cons_var_RTF], &p_ip_RTF, thermo_properties_ptr);
                        
                        // Compute the temperature values in the stencils for image point 2.
                        const Real T_ip2_LBK = d_equation_of_state_mixing_rules->getEquationOfState()->getTemperature(&rho[idx_ip2_cons_var_LBK], &p_ip2_LBK, thermo_properties_ptr);
                        const Real T_ip2_RBK = d_equation_of_state_mixing_rules->getEquationOfState()->getTemperature(&rho[idx_ip2_cons_var_RBK], &p_ip2_RBK, thermo_properties_ptr);
                        const Real T_ip2_LTK = d_equation_of_state_mixing_rules->getEquationOfState()->getTemperature(&rho[idx_ip2_cons_var_LTK], &p_ip2_LTK, thermo_properties_ptr);
                        const Real T_ip2_RTK = d_equation_of_state_mixing_rules->getEquationOfState()->getTemperature(&rho[idx_ip2_cons_var_RTK], &p_ip2_RTK, thermo_properties_ptr);
                        const Real T_ip2_LBF = d_equation_of_state_mixing_rules->getEquationOfState()->getTemperature(&rho[idx_ip2_cons_var_LBF], &p_ip2_LBF, thermo_properties_ptr);
                        const Real T_ip2_RBF = d_equation_of_state_mixing_rules->getEquationOfState()->getTemperature(&rho[idx_ip2_cons_var_RBF], &p_ip2_RBF, thermo_properties_ptr);
                        const Real T_ip2_LTF = d_equation_of_state_mixing_rules->getEquationOfState()->getTemperature(&rho[idx_ip2_cons_var_LTF], &p_ip2_LTF, thermo_properties_ptr);
                        const Real T_ip2_RTF = d_equation_of_state_mixing_rules->getEquationOfState()->getTemperature(&rho[idx_ip2_cons_var_RTF], &p_ip2_RTF, thermo_properties_ptr);
                        
                        // Trilinear interpolation to find temperature of image point 1.
                        const Real T_ip = trilinearInterpolate3D(
                            T_ip_LBK,
                            T_ip_RBK,
                            T_ip_LTK,
                            T_ip_RTK,
                            T_ip_LBF,
                            T_ip_RBF,
                            T_ip_LTF,
                            T_ip_RTF,
                            x_ip,
                            y_ip,
                            z_ip,
                            x_ip_LBK,
                            y_ip_LBK,
                            z_ip_LBK,
                            dx_inv);
                        
                        // Trilinear interpolation to find temperature of image point 2.
                        const Real T_ip2 = trilinearInterpolate3D(
                            T_ip2_LBK,
                            T_ip2_RBK,
                            T_ip2_LTK,
                            T_ip2_RTK,
                            T_ip2_LBF,
                            T_ip2_RBF,
                            T_ip2_LTF,
                            T_ip2_RTF,
                            x_ip2,
                            y_ip2,
                            z_ip2,
                            x_ip2_LBK,
                            y_ip2_LBK,
                            z_ip2_LBK,
                            dx_inv);
                        
                        // Trilinear interpolation to find density of image point 1.
                        const Real rho_ip = trilinearInterpolate3D(
                            rho[idx_ip_cons_var_LBK],
                            rho[idx_ip_cons_var_RBK],
                            rho[idx_ip_cons_var_LTK],
                            rho[idx_ip_cons_var_RTK],
                            rho[idx_ip_cons_var_LBF],
                            rho[idx_ip_cons_var_RBF],
                            rho[idx_ip_cons_var_LTF],
                            rho[idx_ip_cons_var_RTF],
                            x_ip,
                            y_ip,
                            z_ip,
                            x_ip_LBK,
                            y_ip_LBK,
                            z_ip_LBK,
                            dx_inv);
                        // Trilinear interpolation to find density of image point 2.
                        const Real rho_ip2 = trilinearInterpolate3D(
                            rho[idx_ip2_cons_var_LBK],
                            rho[idx_ip2_cons_var_RBK],
                            rho[idx_ip2_cons_var_LTK],
                            rho[idx_ip2_cons_var_RTK],
                            rho[idx_ip2_cons_var_LBF],
                            rho[idx_ip2_cons_var_RBF],
                            rho[idx_ip2_cons_var_LTF],
                            rho[idx_ip2_cons_var_RTF],
                            x_ip2,
                            y_ip2,
                            z_ip2,
                            x_ip2_LBK,
                            y_ip2_LBK,
                            z_ip2_LBK,
                            dx_inv);
                        
                        const Real epsilon_ip = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getInternalEnergyFromTemperature(
                                &rho_ip,
                                &T_ip,
                                thermo_properties_ptr);
                            
                        const Real epsilon_ip2 = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getInternalEnergyFromTemperature(
                                &rho_ip2,
                                &T_ip2,
                                thermo_properties_ptr);
                        
                        const Real p_ip = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getPressure(
                                &rho_ip,
                                &epsilon_ip,
                                thermo_properties_ptr);
                        
                        const Real p_ip2 = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getPressure(
                                &rho_ip2,
                                &epsilon_ip2,
                                thermo_properties_ptr);
                        
                        // dP/dn = 0
                        const Real p_gc = getGhostValueNeumannBC(
                            p_ip,
                            p_ip2,
                            d_ip,
                            d_ip2,
                            dist[idx_IB]);
                        
                        Real T_gc = Real(0);
                        if (d_bc_type_temperature == TEMPERATURE_IBC::ADIABATIC)
                        {   
                            // dT/dn = 0
                            T_gc = getGhostValueNeumannBC(
                                T_ip,
                                T_ip2,
                                d_ip,
                                d_ip2,
                                dist[idx_IB]);
                        }
                        else if (d_bc_type_temperature == TEMPERATURE_IBC::ISOTHERMAL)
                        {
                            // Iso-thermal boundary condition.
                            T_gc = getGhostValueDirichletBC(
                                T_body,
                                T_ip,
                                d_ip,
                                dist[idx_IB]);
                        }
                        
                        // Compute density using calculated p and T.
                        const Real rho_gc = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getDensity(
                                &p_gc,
                                &T_gc,
                                thermo_properties_ptr);
                        
                        // Compute the total energy at the ghost cell.
                        const Real epsilon_gc = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getInternalEnergyFromTemperature(
                                &rho_gc,
                                &T_gc,
                                thermo_properties_ptr);
                        
                        const Real E_gc = rho_gc*(epsilon_gc + half*(u_gc*u_gc + v_gc*v_gc + w_gc*w_gc));
                        
                        rho[idx_cons_var]   = rho_gc;
                        rho_u[idx_cons_var] = rho_gc*u_gc;
                        rho_v[idx_cons_var] = rho_gc*v_gc;
                        rho_w[idx_cons_var] = rho_gc*w_gc;
                        E[idx_cons_var]     = E_gc;
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


/*
 * Output the surface triangulation with surface data.
 */
void FlowModelImmersedBoundaryMethodSingleSpecies::writeSurfaceTriangulationWithData(
    const std::string& file_name) const
{
#ifdef HAMERS_USE_TECIO
    const SurfaceTriangulation& surface_triangulation = d_immersed_boundaries->getSurfaceTriangulation();
    if (surface_triangulation.nodes.size() == 0)
    {
        TBOX_WARNING(d_object_name
            << ": FlowModelImmersedBoundaryMethodSingleSpecies::writeSurfaceTriangulationWithData()\n"
            << "The surface triangulation is empty."
            << " No surface file will be written."
            << std::endl);
        return;
    }
    
    std::vector<std::string> variable_names;
    std::vector<HAMERS_SHARED_PTR<std::vector<double> > > variable_data;
    
    // Add surface pressure.
    variable_names.push_back("surf_data_pressure");
    variable_data.push_back(d_surface_triangulation_p);
    // Add surface x-component of velocity.
    variable_names.push_back("surf_data_u");
    variable_data.push_back(d_surface_triangulation_u);
    // Add surface y-component of velocity.
    variable_names.push_back("surf_data_v");
    variable_data.push_back(d_surface_triangulation_v);
    // Add surface z-component of velocity.
    variable_names.push_back("surf_data_w");
    variable_data.push_back(d_surface_triangulation_w);
    // Add surface temperature.
    variable_names.push_back("surf_data_temperature");
    variable_data.push_back(d_surface_triangulation_T);
    // Add surface density.
    variable_names.push_back("surf_data_density");
    variable_data.push_back(d_surface_triangulation_rho);
    // Add surface traction x-component.
    variable_names.push_back("surf_data_traction_v_x");
    variable_data.push_back(d_surface_triangulation_tx_v);
    // Add surface traction y-component.
    variable_names.push_back("surf_data_traction_v_y");
    variable_data.push_back(d_surface_triangulation_ty_v);
    // Add surface traction z-component.
    variable_names.push_back("surf_data_traction_v_z");
    variable_data.push_back(d_surface_triangulation_tz_v);
    
    writeSurfaceTriangulationWithDataBase(file_name, variable_names, variable_data);
#endif
}


/*
 * Output names of monitoring statistical quantities to output to a file.
 */
 void FlowModelImmersedBoundaryMethodSingleSpecies::outputMonitoringStatisticalQuantitiesNames(
    std::ofstream& f_out) const
{
#ifdef HAMERS_USE_TECIO
    f_out << "\t" << "F_p_x                ";
    f_out << "\t" << "F_p_y                ";
    f_out << "\t" << "F_p_z                ";
    f_out << "\t" << "F_v_x                ";
    f_out << "\t" << "F_v_y                ";
    f_out << "\t" << "F_v_z                ";
#else
    NULL_USE(f_out);
#endif
}


/*
 * Output monitoring statistics to screen.
 */
void FlowModelImmersedBoundaryMethodSingleSpecies::outputMonitoringStatistics(
    std::ofstream& f_out) const
{
#ifdef HAMERS_USE_TECIO
    f_out << std::scientific << std::setprecision(std::numeric_limits<Real>::digits10)
        << "\t" << d_surface_triangulation_integrated_F_p_x;
    f_out << std::scientific << std::setprecision(std::numeric_limits<Real>::digits10)
        << "\t" << d_surface_triangulation_integrated_F_p_y;
    f_out << std::scientific << std::setprecision(std::numeric_limits<Real>::digits10)
        << "\t" << d_surface_triangulation_integrated_F_p_z;
    f_out << std::scientific << std::setprecision(std::numeric_limits<Real>::digits10)
        << "\t" << d_surface_triangulation_integrated_F_v_x;
    f_out << std::scientific << std::setprecision(std::numeric_limits<Real>::digits10)
        << "\t" << d_surface_triangulation_integrated_F_v_y;
    f_out << std::scientific << std::setprecision(std::numeric_limits<Real>::digits10)
        << "\t" << d_surface_triangulation_integrated_F_v_z;
#else
    NULL_USE(f_out);
#endif
}


/*
 * Compute the data on the surface triangulation.
 */
void FlowModelImmersedBoundaryMethodSingleSpecies::computeSurfaceTriangulationData(
    const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    computeSurfaceTriangulationDataBase(
        grid_geometry,
        patch_hierarchy,
        data_context);
    
    const SurfaceTriangulation& surface_triangulation = d_immersed_boundaries->getSurfaceTriangulation();
    const std::vector<std::array<double, 3> >& nodes = surface_triangulation.nodes;
    const std::vector<std::array<double, 3> >& normal_nodes = surface_triangulation.normal_nodes;
    
    if (nodes.empty())
    {
        return;
    }
    
    const std::vector<std::array<int, 3> >& connectivities = surface_triangulation.connectivities;
    const std::vector<std::array<double, 3> >& normal_centroids = surface_triangulation.normal_centroids;
    const std::vector<double>& areas = surface_triangulation.areas;
    
    const int num_nodes = static_cast<int>(nodes.size());
    const int num_centroids = static_cast<int>(connectivities.size());
    
    d_surface_triangulation_p    = HAMERS_SHARED_PTR<std::vector<double> >(new std::vector<double>(nodes.size(), 0.0));
    d_surface_triangulation_u    = HAMERS_SHARED_PTR<std::vector<double> >(new std::vector<double>(nodes.size(), 0.0));
    d_surface_triangulation_v    = HAMERS_SHARED_PTR<std::vector<double> >(new std::vector<double>(nodes.size(), 0.0));
    d_surface_triangulation_w    = HAMERS_SHARED_PTR<std::vector<double> >(new std::vector<double>(nodes.size(), 0.0));
    d_surface_triangulation_T    = HAMERS_SHARED_PTR<std::vector<double> >(new std::vector<double>(nodes.size(), 0.0));
    d_surface_triangulation_rho  = HAMERS_SHARED_PTR<std::vector<double> >(new std::vector<double>(nodes.size(), 0.0));
    d_surface_triangulation_tx_v = HAMERS_SHARED_PTR<std::vector<double> >(new std::vector<double>(nodes.size(), 0.0));
    d_surface_triangulation_ty_v = HAMERS_SHARED_PTR<std::vector<double> >(new std::vector<double>(nodes.size(), 0.0));
    d_surface_triangulation_tz_v = HAMERS_SHARED_PTR<std::vector<double> >(new std::vector<double>(nodes.size(), 0.0));
    
    double* p_data    = d_surface_triangulation_p->data();
    double* u_data    = d_surface_triangulation_u->data();
    double* v_data    = d_surface_triangulation_v->data();
    double* w_data    = d_surface_triangulation_w->data();
    double* T_data    = d_surface_triangulation_T->data();
    double* rho_data  = d_surface_triangulation_rho->data();
    double* tx_v_data = d_surface_triangulation_tx_v->data();
    double* ty_v_data = d_surface_triangulation_ty_v->data();
    double* tz_v_data = d_surface_triangulation_tz_v->data();
    
    // Get the thermodynamic properties of the species.
    std::vector<const Real*> thermo_properties_ptr;
    thermo_properties_ptr.reserve(static_cast<int> (d_thermo_properties.size()));
    for (int ti = 0; ti < static_cast<int> (d_thermo_properties.size()); ti++)
    {
        thermo_properties_ptr.push_back(&d_thermo_properties[ti]);
    }
    
    /*
     * Get the molecular properties of the species for shear viscosity.
     */
    
    std::vector<Real> molecular_properties_shear_viscosity;
    std::vector<Real*> molecular_properties_shear_viscosity_ptr;
    std::vector<const Real*> molecular_properties_shear_viscosity_const_ptr;
    
    const int num_molecular_properties_shear_viscosity = d_equation_of_shear_viscosity_mixing_rules->
        getNumberOfSpeciesMolecularProperties();
    
    molecular_properties_shear_viscosity.resize(num_molecular_properties_shear_viscosity);
    molecular_properties_shear_viscosity_ptr.reserve(num_molecular_properties_shear_viscosity);
    molecular_properties_shear_viscosity_const_ptr.reserve(num_molecular_properties_shear_viscosity);
    
    for (int ti = 0; ti < num_molecular_properties_shear_viscosity; ti++)
    {
        molecular_properties_shear_viscosity_ptr.push_back(&molecular_properties_shear_viscosity[ti]);
        molecular_properties_shear_viscosity_const_ptr.push_back(&molecular_properties_shear_viscosity[ti]);
    }
    
    d_equation_of_shear_viscosity_mixing_rules->getSpeciesMolecularProperties(
        molecular_properties_shear_viscosity_ptr,
        0);
    
    /*
     * Get the molecular properties of the species for bulk viscosity.
     */
    
    std::vector<Real> molecular_properties_bulk_viscosity;
    std::vector<Real*> molecular_properties_bulk_viscosity_ptr;
    std::vector<const Real*> molecular_properties_bulk_viscosity_const_ptr;
    
    const int num_molecular_properties_bulk_viscosity = d_equation_of_bulk_viscosity_mixing_rules->
        getNumberOfSpeciesMolecularProperties();
        
    molecular_properties_bulk_viscosity.resize(num_molecular_properties_bulk_viscosity);
    molecular_properties_bulk_viscosity_ptr.reserve(num_molecular_properties_bulk_viscosity);
    molecular_properties_bulk_viscosity_const_ptr.reserve(num_molecular_properties_bulk_viscosity);
    
    for (int ti = 0; ti < num_molecular_properties_bulk_viscosity; ti++)
    {
        molecular_properties_bulk_viscosity_ptr.push_back(&molecular_properties_bulk_viscosity[ti]);
        molecular_properties_bulk_viscosity_const_ptr.push_back(&molecular_properties_bulk_viscosity[ti]);
    }
    
    d_equation_of_bulk_viscosity_mixing_rules->getSpeciesMolecularProperties(
        molecular_properties_bulk_viscosity_ptr,
        0);
    
    if (d_dim == tbox::Dimension(1))
    {
        // Do nothing for now.
    }
    else if (d_dim == tbox::Dimension(2))
    {
        // Do nothing for now.
    }
    else if (d_dim == tbox::Dimension(3))
    {
        for (int ni = 0; ni < num_nodes; ni++)
        {
            const double dx_grid   = d_surface_triangulation_dx_grid[ni];
            const Real vec_norm[3] = {Real(normal_nodes[ni][0]), Real(normal_nodes[ni][1]), Real(normal_nodes[ni][2])};
            
            const Real d_ip_1 = Real(d_surface_triangulation_coeff_ip_1*dx_grid);
            const Real d_ip_2 = Real(d_surface_triangulation_coeff_ip_2*dx_grid);
            
            const Real rho_ip_1   = Real(d_surface_triangulation_cons_var_ip_1[0][ni]);
            const Real rho_u_ip_1 = Real(d_surface_triangulation_cons_var_ip_1[1][ni]);
            const Real rho_v_ip_1 = Real(d_surface_triangulation_cons_var_ip_1[2][ni]);
            const Real rho_w_ip_1 = Real(d_surface_triangulation_cons_var_ip_1[3][ni]);
            const Real E_ip_1     = Real(d_surface_triangulation_cons_var_ip_1[4][ni]);
            
            const Real rho_ip_2   = Real(d_surface_triangulation_cons_var_ip_2[0][ni]);
            const Real rho_u_ip_2 = Real(d_surface_triangulation_cons_var_ip_2[1][ni]);
            const Real rho_v_ip_2 = Real(d_surface_triangulation_cons_var_ip_2[2][ni]);
            const Real rho_w_ip_2 = Real(d_surface_triangulation_cons_var_ip_2[3][ni]);
            const Real E_ip_2     = Real(d_surface_triangulation_cons_var_ip_2[4][ni]);
            
            const Real epsilon_ip_1 = (E_ip_1 -
                0.5*(rho_u_ip_1*rho_u_ip_1 + rho_v_ip_1*rho_v_ip_1 + rho_w_ip_1*rho_w_ip_1)/rho_ip_1)/rho_ip_1;
            
            const Real epsilon_ip_2 = (E_ip_2 -
                0.5*(rho_u_ip_2*rho_u_ip_2 + rho_v_ip_2*rho_v_ip_2 + rho_w_ip_2*rho_w_ip_2)/rho_ip_2)/rho_ip_2;
            
            const Real p_ip_1 = d_equation_of_state_mixing_rules->getEquationOfState()->getPressure(
                &rho_ip_1,
                &epsilon_ip_1,
                thermo_properties_ptr);
            
            const Real p_ip_2 = d_equation_of_state_mixing_rules->getEquationOfState()->getPressure(
                &rho_ip_2,
                &epsilon_ip_2,
                thermo_properties_ptr);
            
            const Real T_ip_1 = d_equation_of_state_mixing_rules->getEquationOfState()->getTemperature(
                &rho_ip_1,
                &p_ip_1,
                thermo_properties_ptr);
            
            const Real T_ip_2 = d_equation_of_state_mixing_rules->getEquationOfState()->getTemperature(
                &rho_ip_2,
                &p_ip_2,
                thermo_properties_ptr);
            
            const Real p_surf = getGhostValueNeumannBC(
                p_ip_1,
                p_ip_2,
                d_ip_1,
                d_ip_2,
                Real(0));
            
            Real T_surf = Real(0);
            if (d_bc_type_temperature == TEMPERATURE_IBC::ADIABATIC)
            {
                // dT/dn = 0
                T_surf = getGhostValueNeumannBC(
                    T_ip_1,
                    T_ip_2,
                    d_ip_1,
                    d_ip_2,
                    Real(0));
            }
            else if (d_bc_type_temperature == TEMPERATURE_IBC::ISOTHERMAL)
            {
                // Dirichlet BC.
                T_surf = getGhostValueDirichletBC(
                    d_T_body,
                    T_ip_1,
                    d_ip_1,
                    Real(0));
            }
            
            const Real rho_surf = d_equation_of_state_mixing_rules->getEquationOfState()->getDensity(
                &p_surf,
                &T_surf,
                thermo_properties_ptr);
            
            const Real mu_surf = d_equation_of_shear_viscosity_mixing_rules->getEquationOfShearViscosity()->
                getShearViscosity(
                    &p_surf,
                    &T_surf,
                    molecular_properties_shear_viscosity_const_ptr);
            
            const Real mu_v_surf = d_equation_of_bulk_viscosity_mixing_rules->getEquationOfBulkViscosity()->
                getBulkViscosity(
                    &p_surf,
                    &T_surf,
                    molecular_properties_bulk_viscosity_const_ptr);
            
            const Real u_ip_1 = rho_u_ip_1/rho_ip_1;
            const Real v_ip_1 = rho_v_ip_1/rho_ip_1;
            const Real w_ip_1 = rho_w_ip_1/rho_ip_1;
            const Real u_ip_2 = rho_u_ip_2/rho_ip_2;
            const Real v_ip_2 = rho_v_ip_2/rho_ip_2;
            const Real w_ip_2 = rho_w_ip_2/rho_ip_2;

            Real u_surf = Real(0); // x-component of velocity at the surface
            Real v_surf = Real(0); // y-component of velocity at the surface
            Real w_surf = Real(0); // z-component of velocity at the surface 

            if (d_bc_type_velocity == VELOCITY_IBC::SLIP)
            {
                // Given d_ip and d_ip2, interpolate to the velocity components (u_mirror and v_mirror)
                // at the mirror image point at the surface.
                const Real diff_ip_2_ip_1 = d_ip_2 - d_ip_1;
                const Real diff_mirror_ip_1  = - d_ip_1;
                const Real diff_ip_2_mirror = d_ip_2;

                // x-component of velocity at the mirror image point.
                const Real u_mirror = (diff_ip_2_mirror*u_ip_1 + diff_mirror_ip_1*u_ip_2)/diff_ip_2_ip_1;
                // y-component of velocity at the mirror image point.
                const Real v_mirror = (diff_ip_2_mirror*v_ip_1 + diff_mirror_ip_1*v_ip_2)/diff_ip_2_ip_1;
                // z-component of velocity at the mirror image point.
                const Real w_mirror = (diff_ip_2_mirror*w_ip_1 + diff_mirror_ip_1*w_ip_2)/diff_ip_2_ip_1;
                
                // Velocity component normal to the boundary at the mirror image point.
                const Real vel_mirror_n = dotProduct3D(u_mirror, v_mirror, w_mirror, normal_nodes[ni][0], normal_nodes[ni][1], normal_nodes[ni][2]);
                
                // No-penetration boundary condition.
                u_surf = u_mirror - Real(2)*vel_mirror_n*normal_nodes[ni][0];
                v_surf = v_mirror - Real(2)*vel_mirror_n*normal_nodes[ni][1];
                w_surf = w_mirror - Real(2)*vel_mirror_n*normal_nodes[ni][2];

            }
            else if (d_bc_type_velocity == VELOCITY_IBC::NO_SLIP)
            {
                u_surf = getGhostValueDirichletBC(
                    Real(0),
                    u_ip_1,
                    d_ip_1,
                    Real(0));
                
                v_surf = getGhostValueDirichletBC(
                    Real(0),
                    v_ip_1,
                    d_ip_1,
                    Real(0));
                
                w_surf = getGhostValueDirichletBC(
                    Real(0),
                    w_ip_1,
                    d_ip_1,
                    Real(0));
            }
            
            // Get the vectors in the two tangent directions.
            // A vector orthogonal to (a, b, c) is (-b, a, 0), or (-c, 0, a) or (0, -c, b).
            Real vec_tan_1[3] = {-vec_norm[1], vec_norm[0], Real(0)};
            Real vec_tan_2[3] = {-vec_norm[2], Real(0), vec_norm[0]};
            // Consider special cases when the normal vector is (1, 0, 0) or (0, 1, 0) or (0, 0, 1).
            if (std::abs(vec_norm[0] - Real(1)) < std::numeric_limits<Real>::epsilon())
            {
                vec_tan_1[0] = Real(0);
                vec_tan_1[1] = Real(1);
                vec_tan_1[2] = Real(0);
                
                vec_tan_2[0] = Real(0);
                vec_tan_2[1] = Real(0);
                vec_tan_2[2] = Real(1);
            }
            else if (std::abs(vec_norm[1] - Real(1)) < std::numeric_limits<Real>::epsilon())
            {
                vec_tan_1[0] = Real(1);
                vec_tan_1[1] = Real(0);
                vec_tan_1[2] = Real(0);
                
                vec_tan_2[0] = Real(0);
                vec_tan_2[1] = Real(0);
                vec_tan_2[2] = Real(1);
            }
            else if (std::abs(vec_norm[2] - Real(1)) < std::numeric_limits<Real>::epsilon())
            {
                vec_tan_1[0] = Real(1);
                vec_tan_1[1] = Real(0);
                vec_tan_1[2] = Real(0);
                
                vec_tan_2[0] = Real(0);
                vec_tan_2[1] = Real(1);
                vec_tan_2[2] = Real(0);
            }
            
            // Normalize the tangent vectors.
            const Real norm_tan_1 = sqrt(vec_tan_1[0]*vec_tan_1[0] + vec_tan_1[1]*vec_tan_1[1] + vec_tan_1[2]*vec_tan_1[2]);
            const Real norm_tan_2 = sqrt(vec_tan_2[0]*vec_tan_2[0] + vec_tan_2[1]*vec_tan_2[1] + vec_tan_2[2]*vec_tan_2[2]);
            
            vec_tan_1[0] /= norm_tan_1;
            vec_tan_1[1] /= norm_tan_1;
            vec_tan_1[2] /= norm_tan_1;
            
            vec_tan_2[0] /= norm_tan_2;
            vec_tan_2[1] /= norm_tan_2;
            vec_tan_2[2] /= norm_tan_2;
            
            // Rotate u_ip_1, v_ip_1, w_ip_1 to the normal direction.
            const Real vel_norm_ip_1  = u_ip_1*vec_norm[0]  + v_ip_1*vec_norm[1]  + w_ip_1*vec_norm[2];
            const Real vel_tan_1_ip_1 = u_ip_1*vec_tan_1[0] + v_ip_1*vec_tan_1[1] + w_ip_1*vec_tan_1[2];
            const Real vel_tan_2_ip_1 = u_ip_1*vec_tan_2[0] + v_ip_1*vec_tan_2[1] + w_ip_1*vec_tan_2[2];
            
            // Rotate u_ip_2, v_ip_2, w_ip_2 to the normal direction.
            const Real vel_norm_ip_2  = u_ip_2*vec_norm[0]  + v_ip_2*vec_norm[1]  + w_ip_2*vec_norm[2];
            const Real vel_tan_1_ip_2 = u_ip_2*vec_tan_1[0] + v_ip_2*vec_tan_1[1] + w_ip_2*vec_tan_1[2];
            const Real vel_tan_2_ip_2 = u_ip_2*vec_tan_2[0] + v_ip_2*vec_tan_2[1] + w_ip_2*vec_tan_2[2];
            
            Real ddn_vel_norm_surf  = Real(0);
            Real ddn_vel_tan_1_surf = Real(0);
            Real ddn_vel_tan_2_surf = Real(0);
            
            if (d_bc_type_velocity == VELOCITY_IBC::SLIP)
            {
                ddn_vel_norm_surf = getGradientBC(
                    Real(0),
                    vel_norm_ip_1,
                    vel_norm_ip_2,
                    d_ip_1,
                    d_ip_2);
                
                // ddn_vel_tan_1_surf = Real(0);
                // ddn_vel_tan_2_surf = Real(0);
            }
            else if (d_bc_type_velocity == VELOCITY_IBC::NO_SLIP)
            {
                ddn_vel_norm_surf = getGradientBC(
                    Real(0),
                    vel_norm_ip_1,
                    vel_norm_ip_2,
                    d_ip_1,
                    d_ip_2);
                
                ddn_vel_tan_1_surf = getGradientBC(
                    Real(0),
                    vel_tan_1_ip_1,
                    vel_tan_1_ip_2,
                    d_ip_1,
                    d_ip_2);
                
                ddn_vel_tan_2_surf = getGradientBC(
                    Real(0),
                    vel_tan_2_ip_1,
                    vel_tan_2_ip_2,
                    d_ip_1,
                    d_ip_2);
            }
            
            const Real t_norm_surf  = Real(2)*mu_surf*ddn_vel_norm_surf - (Real(2)/Real(3)*mu_surf - mu_v_surf)*ddn_vel_norm_surf;
            const Real t_tan_1_surf = mu_surf*ddn_vel_tan_1_surf;
            const Real t_tan_2_surf = mu_surf*ddn_vel_tan_2_surf;
            
            const Real tx_v_surf = t_norm_surf*vec_norm[0] + t_tan_1_surf*vec_tan_1[0] + t_tan_2_surf*vec_tan_2[0];
            const Real ty_v_surf = t_norm_surf*vec_norm[1] + t_tan_1_surf*vec_tan_1[1] + t_tan_2_surf*vec_tan_2[1];
            const Real tz_v_surf = t_norm_surf*vec_norm[2] + t_tan_1_surf*vec_tan_1[2] + t_tan_2_surf*vec_tan_2[2];
            
            p_data[ni]    = double(p_surf);
            u_data[ni]    = double(u_surf);
            v_data[ni]    = double(v_surf);
            w_data[ni]    = double(w_surf);
            T_data[ni]    = double(T_surf);
            rho_data[ni]  = double(rho_surf);
            tx_v_data[ni] = double(tx_v_surf);
            ty_v_data[ni] = double(ty_v_surf);
            tz_v_data[ni] = double(tz_v_surf);
        }
        
        d_surface_triangulation_integrated_F_p_x = 0.0;
        d_surface_triangulation_integrated_F_p_y = 0.0;
        d_surface_triangulation_integrated_F_p_z = 0.0;
        d_surface_triangulation_integrated_F_v_x = 0.0;
        d_surface_triangulation_integrated_F_v_y = 0.0;
        d_surface_triangulation_integrated_F_v_z = 0.0;
        
        for (int ci = 0; ci < num_centroids; ci++)
        {
            const int& node_0 = connectivities[ci][0];
            const int& node_1 = connectivities[ci][1];
            const int& node_2 = connectivities[ci][2];
            
            const double& p_0 = p_data[node_0];
            const double& p_1 = p_data[node_1];
            const double& p_2 = p_data[node_2];
            
            const double& tx_v_0 = tx_v_data[node_0];
            const double& tx_v_1 = tx_v_data[node_1];
            const double& tx_v_2 = tx_v_data[node_2];
            const double& ty_v_0 = ty_v_data[node_0];
            const double& ty_v_1 = ty_v_data[node_1];
            const double& ty_v_2 = ty_v_data[node_2];
            const double& tz_v_0 = tz_v_data[node_0];
            const double& tz_v_1 = tz_v_data[node_1];
            const double& tz_v_2 = tz_v_data[node_2];
            
            const double p_centroid = (p_0 + p_1 + p_2)/3.0;
            const double tx_v_centroid = (tx_v_0 + tx_v_1 + tx_v_2)/3.0;
            const double ty_v_centroid = (ty_v_0 + ty_v_1 + ty_v_2)/3.0;
            const double tz_v_centroid = (tz_v_0 + tz_v_1 + tz_v_2)/3.0;
            
            d_surface_triangulation_integrated_F_p_x -= normal_centroids[ci][0]*p_centroid*areas[ci];
            d_surface_triangulation_integrated_F_p_y -= normal_centroids[ci][1]*p_centroid*areas[ci];
            d_surface_triangulation_integrated_F_p_z -= normal_centroids[ci][2]*p_centroid*areas[ci];
            
            d_surface_triangulation_integrated_F_v_x += tx_v_centroid*areas[ci];
            d_surface_triangulation_integrated_F_v_y += ty_v_centroid*areas[ci];
            d_surface_triangulation_integrated_F_v_z += tz_v_centroid*areas[ci];
        }
        
    } // if (d_dim == tbox::Dimension(3))
}
