#include "flow/flow_models/five-eqn_Allaire/FlowModelImmersedBoundaryMethodFiveEqnAllaire.hpp"

FlowModelImmersedBoundaryMethodFiveEqnAllaire::FlowModelImmersedBoundaryMethodFiveEqnAllaire(
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
     * Read the body partial densities.
     */
    
    if (immersed_boundary_method_db->keyExists("body_partial_densities"))
    {
        d_Z_rho_body = immersed_boundary_method_db->getRealVector("body_partial_densities");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelImmersedBoundaryMethodFiveEqnAllaire::FlowModelImmersedBoundaryMethodFiveEqnAllaire()\n"
            << "Required 'body_partial_densities' entry from input database missing."
            << std::endl);
    }
    
    if (static_cast<int>(d_Z_rho_body.size()) != d_num_species)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelImmersedBoundaryMethodFiveEqnAllaire::FlowModelImmersedBoundaryMethodFiveEqnAllaire()\n"
            << "Size of 'body_partial_densities' entry from input database does not the number of species."
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
            << ": FlowModelImmersedBoundaryMethodFiveEqnAllaire::FlowModelImmersedBoundaryMethodFiveEqnAllaire()\n"
            << "Required 'body_velocity' entry from input database missing."
            << std::endl);
    }
    
    if (static_cast<int>(d_vel_body.size()) != d_dim.getValue())
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelImmersedBoundaryMethodFiveEqnAllaire::FlowModelImmersedBoundaryMethodFiveEqnAllaire()\n"
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
            << ": FlowModelImmersedBoundaryMethodFiveEqnAllaire::FlowModelImmersedBoundaryMethodFiveEqnAllaire()\n"
            << "Required 'body_pressure' entry from input database missing."
            << std::endl);
    }
    
    /*
     * Read the body volume fractions.
     */
    
    if (immersed_boundary_method_db->keyExists("body_volume_fractions"))
    {
        d_Z_body = immersed_boundary_method_db->getRealVector("body_volume_fractions");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelImmersedBoundaryMethodFiveEqnAllaire::FlowModelImmersedBoundaryMethodFiveEqnAllaire()\n"
            << "Required 'body_volume_fractions' entry from input database missing."
            << std::endl);
    }
    
    if (static_cast<int>(d_Z_body.size()) != d_num_species)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelImmersedBoundaryMethodFiveEqnAllaire::FlowModelImmersedBoundaryMethodFiveEqnAllaire()\n"
            << "Size of 'body_volume_fractions' entry from input database does not the number of species."
            << std::endl);
    }
    
    Real Z_sum = Real(0);
    for (int si = 0; si < d_num_species; si++)
    {
        Z_sum += d_Z_body[si];
    }
    
    if (std::abs(Z_sum - Real(1)) > HAMERS_EPSILON*Real(1000))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelImmersedBoundaryMethodFiveEqnAllaire::FlowModelImmersedBoundaryMethodFiveEqnAllaire()\n"
            << "'body_volume_fractions' do not sum up to one."
            << std::endl);
    }
    else
    {
        // Make d_Z_body straightly sum up to one.
        d_Z_body[d_num_species - 1] = Real(1);
        for (int si = 0; si < d_num_species - 1; si++)
        {
            d_Z_body[d_num_species - 1] -= d_Z_body[si];
        }
    }
    
    /*
     * Compute the mixture density of the body.
     */
    
    Real rho_body = Real(0);
    for (int si = 0; si < d_num_species; si++)
    {
        rho_body += d_Z_rho_body[si];
    }
    
    /*
     * Get the pointers to the volume fractions of the body.
     */
    
    std::vector<const Real*> Z_body_ptr;
    Z_body_ptr.reserve(d_num_species);
    for (int si = 0; si < d_num_species; si++)
    {
        Z_body_ptr.push_back(&d_Z_body[si]);
    }
    
    /*
     * Compute the total energy and the momentum of the body.
     */
    
    Real epsilon_body = Real(0);
    
    for (int si = 0; si < d_num_species; si++)
    {
        std::vector<Real> species_thermo_properties;
        std::vector<Real*> species_thermo_properties_ptr;
        std::vector<const Real*> species_thermo_properties_const_ptr;
        
        const int num_thermo_properties = d_equation_of_state_mixing_rules->
            getNumberOfSpeciesThermodynamicProperties(si);
        
        species_thermo_properties.resize(num_thermo_properties);
        species_thermo_properties_ptr.reserve(num_thermo_properties);
        species_thermo_properties_const_ptr.reserve(num_thermo_properties);
        
        for (int ti = 0; ti < num_thermo_properties; ti++)
        {
            species_thermo_properties_ptr.push_back(&species_thermo_properties[ti]);
            species_thermo_properties_const_ptr.push_back(&species_thermo_properties[ti]);
        }
        
        d_equation_of_state_mixing_rules->getSpeciesThermodynamicProperties(
            species_thermo_properties_ptr,
            si);
        
        const Real rho_i_body = d_Z_rho_body[si]/d_Z_body[si];
        
        const Real T_i_body = d_equation_of_state_mixing_rules->getEquationOfState(si)->
            getTemperature(
                &rho_i_body,
                &d_p_body,
                species_thermo_properties_const_ptr);
        
        const Real epsilon_i_body = d_equation_of_state_mixing_rules->getEquationOfState(si)->
            getInternalEnergyFromTemperature(
                &rho_i_body,
                &T_i_body,
                species_thermo_properties_const_ptr);
        
        epsilon_body += epsilon_i_body*d_Z_rho_body[si];
    }
    
    epsilon_body /= rho_body;
    
    d_mom_body.resize(d_dim.getValue());
    
    if (d_dim == tbox::Dimension(1))
    {
        const Real& u_body = d_vel_body[0];
        
        d_mom_body[0] = rho_body*u_body;;
        
        d_E_body = rho_body*(epsilon_body + Real(1)/Real(2)*
            (u_body*u_body));
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const Real& u_body = d_vel_body[0];
        const Real& v_body = d_vel_body[1];
        
        d_mom_body[0] = rho_body*u_body;
        d_mom_body[1] = rho_body*v_body;
        
        d_E_body = rho_body*(epsilon_body + Real(1)/Real(2)*
            (u_body*u_body + v_body*v_body));
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const Real& u_body = d_vel_body[0];
        const Real& v_body = d_vel_body[1];
        const Real& w_body = d_vel_body[2];
        
        d_mom_body[0] = rho_body*u_body;
        d_mom_body[1] = rho_body*v_body;
        d_mom_body[2] = rho_body*w_body;
        
        d_E_body = rho_body*(epsilon_body + Real(1)/Real(2)*
            (u_body*u_body + v_body*v_body + w_body*w_body));
    }
}


/*
 * Set the immersed boundary method ghost cells for the cell data of conservative variables.
 */
void FlowModelImmersedBoundaryMethodFiveEqnAllaire::setConservativeVariablesCellDataImmersedBoundaryGhosts(
    const hier::Patch& patch,
    const double data_time,
    const bool initial_time,
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& conservative_var_data,
    const HAMERS_SHARED_PTR<pdat::CellData<int> >& data_mask,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_wall_distance,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_surface_normal,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_ip_index, //AFK
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_ip_corr,  //AFK
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
    
    // const double* const dx = patch_geom->getDx();
    // const double* const patch_xlo = patch_geom->getXLower();
    
    /*
     * Get the data of the conservative variables.
     */
    
    HAMERS_SHARED_PTR<pdat::CellData<Real> > data_partial_densities = conservative_var_data[0];
    HAMERS_SHARED_PTR<pdat::CellData<Real> > data_momentum          = conservative_var_data[1];
    HAMERS_SHARED_PTR<pdat::CellData<Real> > data_total_energy      = conservative_var_data[2];
    HAMERS_SHARED_PTR<pdat::CellData<Real> > data_volume_fractions  = conservative_var_data[3];
    
    /*
     * Get the pointers to the data.
     */
    
    std::vector<Real*> Z_rho;
    Z_rho.reserve(d_num_species);
    for (int si = 0; si < d_num_species; si++)
    {
        Z_rho.push_back(data_partial_densities->getPointer(si));
    }
    Real* E = data_total_energy->getPointer(0);
    std::vector<Real*> Z;
    Z.reserve(d_num_species);
    for (int si = 0; si < d_num_species; si++)
    {
        Z.push_back(data_volume_fractions->getPointer(si));
    }
    
    int* mask = data_mask->getPointer(0);
    // Real* dist = data_wall_distance->getPointer(0);
    
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
        
        for (int si = 0; si < d_num_species; si++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_cons_var = i + offset_0_cons_var;
                const int idx_IB = i + offset_0_IB;
                
                if (mask[idx_IB] == int(IB_MASK::IB_GHOST))
                {
                    // NEED TO BE CHANGED!!!
                    Z_rho[si][idx_cons_var] = d_Z_rho_body[si];
                    Z[si][idx_cons_var]     = d_Z_body[si];
                }
                else if (mask[idx_IB] == int(IB_MASK::BODY))
                {
                    Z_rho[si][idx_cons_var] = d_Z_rho_body[si];
                    Z[si][idx_cons_var]     = d_Z_body[si];
                }
            }
        }
        
        HAMERS_PRAGMA_SIMD
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_cons_var = i + offset_0_cons_var;
            const int idx_IB = i + offset_0_IB;
            
            if (mask[idx_IB] == int(IB_MASK::IB_GHOST))
            {
                // NEED TO BE CHANGED!!!
                rho_u[idx_cons_var] = rho_u_body;
                E[idx_cons_var]     = d_E_body;
            }
            else if (mask[idx_IB] == int(IB_MASK::BODY))
            {
                rho_u[idx_cons_var] = rho_u_body;
                E[idx_cons_var]     = d_E_body;
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
        
        // Real* norm_0 = data_surface_normal->getPointer(0);
        // Real* norm_1 = data_surface_normal->getPointer(1);
        
        for (int si = 0; si < d_num_species; si++)
        {
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
                    
                    if (mask[idx_IB] == int(IB_MASK::IB_GHOST))
                    {
                        // NEED TO BE CHANGED!!!
                        Z_rho[si][idx_cons_var] = d_Z_rho_body[si];
                        Z[si][idx_cons_var]     = d_Z_body[si];
                    }
                    else if (mask[idx_IB] == int(IB_MASK::BODY))
                    {
                        Z_rho[si][idx_cons_var] = d_Z_rho_body[si];
                        Z[si][idx_cons_var]     = d_Z_body[si];
                    }
                }
            }
        }
        
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
                
                if (mask[idx_IB] == int(IB_MASK::IB_GHOST))
                {
                    // NEED TO BE CHANGED!!!
                    rho_u[idx_cons_var] = rho_u_body;
                    rho_v[idx_cons_var] = rho_v_body;
                    E[idx_cons_var]     = d_E_body;
                }
                else if (mask[idx_IB] == int(IB_MASK::BODY))
                {
                    rho_u[idx_cons_var] = rho_u_body;
                    rho_v[idx_cons_var] = rho_v_body;
                    E[idx_cons_var]     = d_E_body;
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
        
        for (int si = 0; si < d_num_species; si++)
        {
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
                            Z_rho[si][idx_cons_var] = d_Z_rho_body[si];
                            Z[si][idx_cons_var]     = d_Z_body[si];
                        }
                        else if (mask[idx_IB] == int(IB_MASK::BODY))
                        {
                            Z_rho[si][idx_cons_var] = d_Z_rho_body[si];
                            Z[si][idx_cons_var]     = d_Z_body[si];
                        }
                    }
                }
            }
        }
        
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
                        rho_u[idx_cons_var] = rho_u_body;
                        rho_v[idx_cons_var] = rho_v_body;
                        rho_w[idx_cons_var] = rho_w_body;
                        E[idx_cons_var]     = d_E_body;
                    }
                    else if (mask[idx_IB] == int(IB_MASK::BODY))
                    {
                        rho_u[idx_cons_var] = rho_u_body;
                        rho_v[idx_cons_var] = rho_v_body;
                        rho_w[idx_cons_var] = rho_w_body;
                        E[idx_cons_var]     = d_E_body;
                    }
                }
            }
        }
    }
}
