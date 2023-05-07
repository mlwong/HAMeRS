#include "flow/flow_models/FlowModelSpecialSourceTerms.hpp"

/*
 * Add the effects of the special source terms.
 */
void
FlowModelSpecialSourceTerms::computeSpecialSourceTermsOnPatch(
    HAMERS_SHARED_PTR<pdat::CellData<double> >& source,
    const hier::Patch& patch,
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& conservative_variables,
    const std::unordered_map<std::string, double>& monitoring_statistics_map,
    const double time,
    const double dt,
    const int RK_step_number)
{
    // Follow Reckinger, Scott J., Daniel Livescu, and Oleg V. Vasilyev.
    // "Comprehensive numerical methodology for direct numerical simulations of compressible Rayleigh–Taylor instability."
    
    if ((d_project_name != "2D discontinuous Rayleigh-Taylor instability") && (d_project_name != "2D smooth Rayleigh-Taylor instability") && (d_project_name != "2D smooth multi-mode Rayleigh-Taylor instability")) 
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Can only initialize data for 'project_name' = '2D discontinuous Rayleigh-Taylor instability' or "
            << "'2D smooth Rayleigh-Taylor instability' or "
            << "'2D smooth multi-mode Rayleigh-Taylor instability'!\n"
            << "'project_name' = '"
            << d_project_name
            << "' is given."
            << std::endl);
    }
    
    if ((d_project_name != "2D discontinuous Rayleigh-Taylor instability") && (d_project_name != "2D smooth Rayleigh-Taylor instability") && (d_project_name != "2D smooth multi-mode Rayleigh-Taylor instability")) 
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Can only initialize data for 'project_name' = '2D discontinuous Rayleigh-Taylor instability' or "
            << "'2D smooth Rayleigh-Taylor instability' or "
            << "'2D smooth multi-mode Rayleigh-Taylor instability'!\n"
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
    
    if (d_special_source_exterior == false)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The 'special_source_exterior' option should be true!"
            << std::endl);
    }
    
    TBOX_ASSERT(d_source_terms_db->keyExists("sponge_rate"));
    
    double sponge_rate = double(0);
    if (d_source_terms_db->keyExists("sponge_rate"))
    {
        sponge_rate = d_source_terms_db->getDouble("sponge_rate");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'sponge_rate' found in data for source terms."
            << std::endl);
    }
    
    TBOX_ASSERT(d_source_terms_db->keyExists("species_mass"));
    
    std::vector<double> species_mass;

    if (d_source_terms_db->keyExists("species_mass"))
    {
        d_source_terms_db->getVector("species_mass", species_mass);
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'species_mass' found in data for source terms."
            << std::endl);
    }

    const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
        HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(patch_geom);
#endif
    
    std::vector<double*> S;
    S.reserve(d_num_eqn);
    for (int si = 0; si < d_num_eqn; si++)
    {
        S.push_back(source->getPointer(si));
    }
    
    /*
     * Get the numbers of ghost cells source and conservative variables.
     */
    const hier::IntVector num_ghosts_source     = source->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_source = source->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_cons_var     = conservative_variables[0]->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_cons_var = conservative_variables[0]->getGhostBox().numberCells();
    
    const double* const dx = patch_geom->getDx();
    const double* const patch_xlo = patch_geom->getXLower();
    
    // Get the dimensions of box that covers the interior of Patch.
    hier::Box patch_box = patch.getBox();
    const hier::IntVector patch_dims = patch_box.numberCells();
    
    /*
     * Initialize data for a 2D Rayleigh-Taylor instability problem (At = 0.04, M = 0.3).
     */
    
    HAMERS_SHARED_PTR<pdat::CellData<double> > partial_density = conservative_variables[0];
    HAMERS_SHARED_PTR<pdat::CellData<double> > momentum        = conservative_variables[1];
    HAMERS_SHARED_PTR<pdat::CellData<double> > total_energy    = conservative_variables[2];
    
    double* rho_Y_0 = partial_density->getPointer(0);
    double* rho_Y_1 = partial_density->getPointer(1);
    double* rho_u   = momentum->getPointer(0);
    double* rho_v   = momentum->getPointer(1);
    double* E       = total_energy->getPointer(0);
    
    const double gamma = double(7)/double(5); // assume both gases have the same ratio of specific heat ratios
    
    double lambda = 701.53278340668; // wavelength of single-mode perturbation
    double eta_0  = 0.01*lambda;      // 1% perturbation
   
    const double p_i = 100000.0; // interface pressure
    const double T_0 = 300.0;    // background temperature
    
    TBOX_ASSERT(d_source_terms_db != nullptr);
    TBOX_ASSERT(d_source_terms_db->keyExists("has_gravity") || d_source_terms_db->keyExists("d_has_gravity"));
    
    std::vector<double> gravity_vector;
    
    if (d_source_terms_db->keyExists("has_gravity"))
    {
        if (d_source_terms_db->keyExists("gravity"))
        {
            d_source_terms_db->getVector("gravity", gravity_vector);
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "No key 'gravity' found in data for source terms."
                << std::endl);
        }
    }
    else if (d_source_terms_db->keyExists("d_has_gravity"))
    {
        if (d_source_terms_db->keyExists("d_gravity"))
        {
            d_source_terms_db->getVector("d_gravity", gravity_vector);
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "No key 'd_gravity' found in data for source terms."
                << std::endl);
        }
    }

    const double g   = gravity_vector[0]; // gravity
    const double W_1 = species_mass[0]; // gravity
    const double W_2 = species_mass[1]; // gravity
    
    
    const double R_u = 8.31446261815324; // universal gas constant
    const double R_1 = R_u/W_1;          // gas constant of heavier gas
    const double R_2 = R_u/W_2;          // gas constant of lighter gas
    
    const double* const domain_xlo = d_grid_geometry->getXLower();
    const double* const domain_xhi = d_grid_geometry->getXUpper();

    if (d_project_name == "2D discontinuous Rayleigh-Taylor instability" ||
        d_project_name == "2D smooth Rayleigh-Taylor instability")
    {
        for (int j = 0; j < patch_dims[1]; j++)
        {
            for (int i = 0; i < patch_dims[0]; i++)
            {
                // Compute the linear indices.
                const int idx_source = (i + num_ghosts_source[0]) +
                    (j + num_ghosts_source[1])*ghostcell_dims_source[0];
                
                const int idx_cons_var = (i + num_ghosts_cons_var[0]) +
                    (j + num_ghosts_cons_var[1])*ghostcell_dims_cons_var[0];
                
        		const int idx_cons_var_L = (i - 1 + num_ghosts_cons_var[0]) +
                    (j + num_ghosts_cons_var[1])*ghostcell_dims_cons_var[0];

        		const int idx_cons_var_R = (i + 1 + num_ghosts_cons_var[0]) +
                    (j + num_ghosts_cons_var[1])*ghostcell_dims_cons_var[0];

        		const int idx_cons_var_B = (i + num_ghosts_cons_var[0]) +
                    (j - 1 + num_ghosts_cons_var[1])*ghostcell_dims_cons_var[0];

        		const int idx_cons_var_T = (i + num_ghosts_cons_var[0]) +
                    (j + 1 + num_ghosts_cons_var[1])*ghostcell_dims_cons_var[0];
		

                // Compute the coordinates.
                double x[2];
                x[0] = patch_xlo[0] + (double(i) + double(1)/double(2))*dx[0];
                x[1] = patch_xlo[1] + (double(j) + double(1)/double(2))*dx[1];
                
        		double x_L;
                x_L = patch_xlo[0] + (double(i) - 1.0 + double(1)/double(2))*dx[0];
		
		        double x_R;
                x_R = patch_xlo[0] + (double(i) + 1.0  + double(1)/double(2))*dx[0];

                // Check whether it is outside the special source box.
                if (!(x[0] >= d_special_source_box_lo[0] && x[0] <= d_special_source_box_hi[0]))
                {
                    double rho_Y_0_ref;
                    double rho_Y_1_ref;
                    double rho_u_ref;
                    double rho_v_ref;
                    double E_ref;
                    double sponge_rate_tot;
                    double xi_b;

        		    double rho_Y_0_ref_L;
                    double rho_Y_1_ref_L;
                    double rho_u_ref_L;
                    double rho_v_ref_L;
                    double E_ref_L;

        		    double rho_Y_0_ref_R;
                    double rho_Y_1_ref_R;
                    double rho_u_ref_R;
                    double rho_v_ref_R;
                    double E_ref_R;

        		    double rho_Y_0_ref_B;
                    double rho_Y_1_ref_B;
                    double rho_u_ref_B;
                    double rho_v_ref_B;
                    double E_ref_B;

        		    double rho_Y_0_ref_T;
                    double rho_Y_1_ref_T;
                    double rho_u_ref_T;
                    double rho_v_ref_T;
                    double E_ref_T;                   
 
                    if (x[0] < 0.0) // left-side
                    {
                        const double rho_ref = p_i/(R_1*T_0)*exp((g*x[0])/(R_1*T_0));
                        rho_Y_0_ref = rho_ref;
	            		const double rho_ref_L = p_i/(R_1*T_0)*exp((g*x_L)/(R_1*T_0));
			            rho_Y_0_ref_L = rho_ref_L;
			            const double rho_ref_R = p_i/(R_1*T_0)*exp((g*x_R)/(R_1*T_0));
                        rho_Y_0_ref_R = rho_ref_R;
                        rho_Y_0_ref_B = rho_ref; // they are at the same height with rho_ref
                        rho_Y_0_ref_T = rho_ref; // they are at the same height with rho_ref
                        rho_Y_1_ref   = 0.0;
			            rho_Y_1_ref_L = 0.0;
			            rho_Y_1_ref_R = 0.0;
			            rho_Y_1_ref_B = 0.0;
			            rho_Y_1_ref_T = 0.0;
                        
                        const double p_ref   = p_i*exp((g*x[0])/(R_1*T_0));
			            const double p_ref_L = p_i*exp((g*x_L)/(R_1*T_0));
			            const double p_ref_R = p_i*exp((g*x_R)/(R_1*T_0));
			            const double p_ref_B = p_ref; // they are at the same height with p_ref
			            const double p_ref_T = p_ref; // they are at the same height with p_ref
                        
                        const double u_ref = 0.0;
                        const double v_ref = 0.0;
                        
                        rho_u_ref         = rho_ref*u_ref;
			            rho_u_ref_L       = rho_ref_L*u_ref;
			            rho_u_ref_R       = rho_ref_R*u_ref;
			            rho_u_ref_B       = rho_ref*u_ref; // they are at the same height with rho_ref
			            rho_u_ref_T       = rho_ref*u_ref; // they are at the same height with rho_ref
                        rho_v_ref         = rho_ref*v_ref;
			            rho_v_ref_L       = rho_ref_L*v_ref;
			            rho_v_ref_R       = rho_ref_R*v_ref;
			            rho_v_ref_B       = rho_ref*v_ref; // they are at the same height with rho_ref
			            rho_v_ref_T       = rho_ref*v_ref; // they are at the same height with rho_ref

                        E_ref             = p_ref/(gamma - double(1)) + double(1)/double(2)*rho_ref*(u_ref*u_ref + v_ref*v_ref);
            			E_ref_L           = p_ref_L/(gamma - double(1)) + double(1)/double(2)*rho_ref_L*(u_ref*u_ref + v_ref*v_ref);
			            E_ref_R           = p_ref_R/(gamma - double(1)) + double(1)/double(2)*rho_ref_R*(u_ref*u_ref + v_ref*v_ref);
			            E_ref_B           = E_ref; // they are at the same height with E_ref
			            E_ref_T           = E_ref; // they are at the same height with E_ref

                        sponge_rate_tot   = (pow((p_ref/rho_ref),0.5))*sponge_rate;
            			const double erf_start_lo  = double(0.75)*(domain_xlo[0]-d_special_source_box_lo[0]) + d_special_source_box_lo[0]; //center of erf is 3/4 of the way into sponge
			            const double erf_offset_lo = double(-0.5) * erf((d_special_source_box_lo[0]-erf_start_lo)/(abs(erf_start_lo)*double(0.35))) + double(0.5); //value of erf at start of sponge
                        xi_b		= double(-0.5) * erf((x[0]-erf_start_lo)/(abs(erf_start_lo)*double(0.35))) + double(0.5) - erf_offset_lo; //subtract value of erf at start of sponge to start xi_b at zero
                        // xi_b            = pow((x[0]-d_special_source_box_lo[0])/(domain_xlo[0]-d_special_source_box_lo[0]),3.0);

                    }
                    else // right-side
                    {
                        const double rho_ref   = p_i/(R_2*T_0)*exp((g*x[0])/(R_2*T_0));
            			rho_Y_1_ref            = rho_ref;
			            const double rho_ref_L = p_i/(R_2*T_0)*exp((g*x_L)/(R_2*T_0));
                        rho_Y_1_ref_L          = rho_ref_L;
                        const double rho_ref_R = p_i/(R_2*T_0)*exp((g*x_R)/(R_2*T_0));
                        rho_Y_1_ref_R          = rho_ref_R;
                        
                        rho_Y_1_ref_B          = rho_ref; // they are at the same height with rho_ref
                        rho_Y_1_ref_T          = rho_ref; // they are at the same height with rho_ref
			            
                        rho_Y_0_ref   = 0.0;
			            rho_Y_0_ref_L = 0.0;
                        rho_Y_0_ref_R = 0.0;
                        rho_Y_0_ref_B = 0.0;
                        rho_Y_0_ref_T = 0.0;
                        
                        const double p_ref = p_i*exp((g*x[0])/(R_2*T_0));
            			const double p_ref_L = p_i*exp((g*x_L)/(R_2*T_0));
                        const double p_ref_R = p_i*exp((g*x_R)/(R_2*T_0));
                        const double p_ref_B = p_ref; // they are at the same height with p_ref
                        const double p_ref_T = p_ref; // they are at the same height with p_ref
                        
                        const double u_ref = 0.0;
                        const double v_ref = 0.0;

            			rho_u_ref         = rho_ref*u_ref;
                        rho_u_ref_L       = rho_ref_L*u_ref;
                        rho_u_ref_R       = rho_ref_R*u_ref;
                        rho_u_ref_B       = rho_ref*u_ref; // they are at the same height with rho_ref
                        rho_u_ref_T       = rho_ref*u_ref; // they are at the same height with rho_ref
                        rho_v_ref         = rho_ref*v_ref;
                        rho_v_ref_L       = rho_ref_L*v_ref;
                        rho_v_ref_R       = rho_ref_R*v_ref;
                        rho_v_ref_B       = rho_ref*v_ref; // they are at the same height with rho_ref
                        rho_v_ref_T       = rho_ref*v_ref; // they are at the same height with rho_ref
			
                        E_ref             = p_ref/(gamma - double(1)) + double(1)/double(2)*rho_ref*(u_ref*u_ref + v_ref*v_ref);
            			E_ref_L           = p_ref_L/(gamma - double(1)) + double(1)/double(2)*rho_ref_L*(u_ref*u_ref + v_ref*v_ref);
                        E_ref_R           = p_ref_R/(gamma - double(1)) + double(1)/double(2)*rho_ref_R*(u_ref*u_ref + v_ref*v_ref);
                        E_ref_B           = E_ref; // they are at the same height with E_ref
                        E_ref_T           = p_ref; // they are at the same height with E_ref

                        sponge_rate_tot = (pow((p_ref/rho_ref),0.5))*sponge_rate;
			            const double erf_start_hi = double(0.75)*(domain_xhi[0]-d_special_source_box_hi[0]) + d_special_source_box_hi[0]; //center of erf is 3/4 of the way into sponge
			            const double erf_offset_hi = double(0.5) * erf((d_special_source_box_hi[0]-erf_start_hi)/(abs(erf_start_hi)*double(0.35))) + double(0.5); //value of erf at start of sponge
                        xi_b            = double(0.5) * erf((x[0]-erf_start_hi)/(abs(erf_start_hi)*double(0.35))) + double(0.5) - erf_offset_hi; //subtract value of erf at start of sponge to start xi_b at zero
                        // xi_b            = pow((x[0]-d_special_source_box_hi[0])/(domain_xhi[0]-d_special_source_box_hi[0]),3.0);
                    }
                    
                    const double rho_Y_0_p = rho_Y_0[idx_cons_var] - rho_Y_0_ref;
                    const double rho_Y_1_p = rho_Y_1[idx_cons_var] - rho_Y_1_ref;
                    const double rho_u_p   = rho_u[idx_cons_var]   - rho_u_ref;
                    const double rho_v_p   = rho_v[idx_cons_var]   - rho_v_ref;
                    const double E_p       = E[idx_cons_var]       - E_ref;
                   
		            const double rho_Y_0_p_L = rho_Y_0[idx_cons_var_L] - rho_Y_0_ref_L;
                    const double rho_Y_1_p_L = rho_Y_1[idx_cons_var_L] - rho_Y_1_ref_L;
                    const double rho_u_p_L   = rho_u[idx_cons_var_L]   - rho_u_ref_L;
                    const double rho_v_p_L   = rho_v[idx_cons_var_L]   - rho_v_ref_L;
                    const double E_p_L       = E[idx_cons_var_L]       - E_ref_L;

        		    const double rho_Y_0_p_R = rho_Y_0[idx_cons_var_R] - rho_Y_0_ref_R;
                    const double rho_Y_1_p_R = rho_Y_1[idx_cons_var_R] - rho_Y_1_ref_R;
                    const double rho_u_p_R   = rho_u[idx_cons_var_R]   - rho_u_ref_R;
                    const double rho_v_p_R   = rho_v[idx_cons_var_R]   - rho_v_ref_R;
                    const double E_p_R       = E[idx_cons_var_R]       - E_ref_R;

        		    const double rho_Y_0_p_B = rho_Y_0[idx_cons_var_B] - rho_Y_0_ref_B;
                    const double rho_Y_1_p_B = rho_Y_1[idx_cons_var_B] - rho_Y_1_ref_B;
                    const double rho_u_p_B   = rho_u[idx_cons_var_B]   - rho_u_ref_B;
                    const double rho_v_p_B   = rho_v[idx_cons_var_B]   - rho_v_ref_B;
                    const double E_p_B       = E[idx_cons_var_B]       - E_ref_B;

        		    const double rho_Y_0_p_T = rho_Y_0[idx_cons_var_T] - rho_Y_0_ref_T;
                    const double rho_Y_1_p_T = rho_Y_1[idx_cons_var_T] - rho_Y_1_ref_T;
                    const double rho_u_p_T   = rho_u[idx_cons_var_T]   - rho_u_ref_T;
                    const double rho_v_p_T   = rho_v[idx_cons_var_T]   - rho_v_ref_T;
                    const double E_p_T       = E[idx_cons_var_T]       - E_ref_T;
		    
                    double sponge_nu       = (dx[0]*dx[0]*dx[1]*dx[1])/(dx[0]*dx[0]+dx[1]*dx[1])/2.0/dt*0.8;
 
                    S[0][idx_source] -= dt*sponge_rate_tot*xi_b*rho_Y_0_p;
                    S[1][idx_source] -= dt*sponge_rate_tot*xi_b*rho_Y_1_p;
                    S[2][idx_source] -= dt*sponge_rate_tot*xi_b*rho_u_p;
                    S[3][idx_source] -= dt*sponge_rate_tot*xi_b*rho_v_p;
                    S[4][idx_source] -= dt*sponge_rate_tot*xi_b*E_p;

		            S[0][idx_source] += dt*sponge_nu*xi_b*( ((rho_Y_0_p_L - 2.0*rho_Y_0_p + rho_Y_0_p_R)/pow(dx[0], 2.0)) + ((rho_Y_0_p_B - 2.0*rho_Y_0_p + rho_Y_0_p_T)/pow(dx[1], 2.0)) );
                    S[1][idx_source] += dt*sponge_nu*xi_b*( ((rho_Y_1_p_L - 2.0*rho_Y_1_p + rho_Y_1_p_R)/pow(dx[0], 2.0)) + ((rho_Y_1_p_B - 2.0*rho_Y_1_p + rho_Y_1_p_T)/pow(dx[1], 2.0)) );
                    S[2][idx_source] += dt*sponge_nu*xi_b*( ((rho_u_p_L - 2.0*rho_u_p + rho_u_p_R)/pow(dx[0], 2.0)) + ((rho_u_p_B - 2.0*rho_u_p + rho_u_p_T)/pow(dx[1], 2.0)) );
                    S[3][idx_source] += dt*sponge_nu*xi_b*( ((rho_v_p_L - 2.0*rho_v_p + rho_v_p_R)/pow(dx[0], 2.0)) + ((rho_v_p_B - 2.0*rho_v_p + rho_v_p_T)/pow(dx[1], 2.0)) );
                    S[4][idx_source] += dt*sponge_nu*xi_b*( ((E_p_L - 2.0*E_p + E_p_R)/pow(dx[0], 2.0)) + ((E_p_B - 2.0*E_p + E_p_T)/pow(dx[1], 2.0)) );
                }
            }
        }
    }
    else if ((d_project_name == "2D smooth multi-mode isopycnic Rayleigh-Taylor instability") ||
             (d_project_name == "2D smooth isopycnic Rayleigh-Taylor instability 3 species"))
    {
        const double delta = 0.02*lambda; // characteristic length of interface.
        const double shift = lambda/4.0;
        double* rho_Y_2    = partial_density->getPointer(2);

        const double W_3 = species_mass[2]; // molecular weight of third fluid if any
        const double R_3 = R_u/W_3;          // gas constant of third gas
        
        const double rho_1 = p_i/(R_1*T_0);
        const double rho_2 = p_i/(R_2*T_0);
        const double rho_3 = p_i/(R_3*T_0);

        for (int j = 0; j < patch_dims[1]; j++)
        {
            for (int i = 0; i < patch_dims[0]; i++)
            {
                // Compute the linear indices.
                const int idx_source = (i + num_ghosts_source[0]) +
                    (j + num_ghosts_source[1])*ghostcell_dims_source[0];
        
                const int idx_cons_var = (i + num_ghosts_cons_var[0]) +
                    (j + num_ghosts_cons_var[1])*ghostcell_dims_cons_var[0];
        
                const int idx_cons_var_L = (i - 1 + num_ghosts_cons_var[0]) +
                    (j + num_ghosts_cons_var[1])*ghostcell_dims_cons_var[0];

                const int idx_cons_var_R = (i + 1 + num_ghosts_cons_var[0]) +
                    (j + num_ghosts_cons_var[1])*ghostcell_dims_cons_var[0];

                const int idx_cons_var_B = (i + num_ghosts_cons_var[0]) +
                    (j - 1 + num_ghosts_cons_var[1])*ghostcell_dims_cons_var[0];

                const int idx_cons_var_T = (i + num_ghosts_cons_var[0]) +
                    (j + 1 + num_ghosts_cons_var[1])*ghostcell_dims_cons_var[0];
                
                // Compute the coordinates.
                double x[2];
                x[0] = patch_xlo[0] + (double(i) + double(1)/double(2))*dx[0];
                x[1] = patch_xlo[1] + (double(j) + double(1)/double(2))*dx[1];
                double x_L;
                x_L = patch_xlo[0] + (double(i) - 1.0 + double(1)/double(2))*dx[0];
		
		        double x_R;
                x_R = patch_xlo[0] + (double(i) + 1.0 + double(1)/double(2))*dx[0];
                    
                // Check whether it is outside the special source box.
                if (!(x[0] >= d_special_source_box_lo[0] && x[0] <= d_special_source_box_hi[0]))
                {
        
                    double rho_Y_0_ref;
                    double rho_Y_1_ref;
                    double rho_Y_2_ref;
                    double rho_u_ref;
                    double rho_v_ref;
                    double E_ref;
                    double sponge_rate_tot;
                    double xi_b;

                    double E_ref_L;

                    double E_ref_R;

                    double E_ref_B;

                    double E_ref_T; 

                    const double Z_2_H = 0.5*(1.0 + erf(((x[0] + shift)/delta))) - 0.5*(1.0 + erf(((x[0] - shift)/delta)));
                    const double Z_3_H = 0.5*(1.0 + erf(((x[0] - shift)/delta)));
                    const double Z_1_H = 1.0 - Z_2_H - Z_3_H;

                    if (x[0] < 0.0) // left-side
                    {
	
            		    //Gideon's Code on 04/19/2023

            		    const double ksi_1 = (x[0] + shift)/delta;
    	                const double ksi_2 = (x[0] - shift)/delta;

    	                const double p = p_i + (g*rho_1*x[0]) + ((0.5*g*delta)*(rho_2 - rho_1)*((ksi_1*erf(ksi_1)) + ((exp(-pow(ksi_1,2.0)))/sqrt(M_PI)))) +
                            ((0.5*g*delta)*(rho_1 - rho_2)*((ksi_2*erf(ksi_2)) + ((exp(-pow(ksi_2,2.0)))/sqrt(M_PI)))) +
  			                (0.5*g*x[0]*(rho_3 - rho_1)) + ((0.5*g*delta*(rho_3 - rho_1))*((exp(-pow(ksi_1,2.0)))/sqrt(M_PI)));

                        const double p_L = p_i + (g*rho_1*x_L) + ((0.5*g*delta)*(rho_2 - rho_1)*((ksi_1*erf(ksi_1)) + ((exp(-pow(ksi_1,2.0)))/sqrt(M_PI)))) +
                            ((0.5*g*delta)*(rho_1 - rho_2)*((ksi_2*erf(ksi_2)) + ((exp(-pow(ksi_2,2.0)))/sqrt(M_PI)))) +
			                (0.5*g*x_L*(rho_3 - rho_1)) + ((0.5*g*delta*(rho_3 - rho_1))*((exp(-pow(ksi_1,2.0)))/sqrt(M_PI)));

                        const double p_R = p_i + (g*rho_1*x_R) + ((0.5*g*delta)*(rho_2 - rho_1)*((ksi_1*erf(ksi_1)) + ((exp(-pow(ksi_1,2.0)))/sqrt(M_PI)))) +
                            ((0.5*g*delta)*(rho_1 - rho_2)*((ksi_2*erf(ksi_2)) + ((exp(-pow(ksi_2,2.0)))/sqrt(M_PI)))) +
			                (0.5*g*x_R*(rho_3 - rho_1)) + ((0.5*g*delta*(rho_3 - rho_1))*((exp(-pow(ksi_1,2.0)))/sqrt(M_PI)));
                    
                        rho_Y_0_ref = rho_1;

                        rho_Y_1_ref = 0.0;
                        rho_Y_2_ref = 0.0;
                    
                        const double u    = 0.0;
                        const double v    = 0.0;
                    
                        rho_u_ref         = rho_1*u;
                        rho_v_ref         = rho_1*v;
                        E_ref             = p/(gamma - double(1)) + double(1)/double(2)*rho_1*(u*u + v*v);
                        E_ref_L           = p_L/(gamma - double(1)) + double(1)/double(2)*rho_1*(u*u + v*v);
                        E_ref_R           = p_R/(gamma - double(1)) + double(1)/double(2)*rho_1*(u*u + v*v);


                        sponge_rate_tot   = (pow((p/rho_1),0.5))*sponge_rate;
            			const double erf_start_lo  = double(0.75)*(domain_xlo[0]-d_special_source_box_lo[0]) + d_special_source_box_lo[0]; //center of erf is 3/4 of the way into sponge
			            const double erf_offset_lo = double(-0.5) * erf((d_special_source_box_lo[0]-erf_start_lo)/(abs(erf_start_lo)*double(0.35))) + double(0.5); //value of erf at start of sponge
                        xi_b		= double(-0.5) * erf((x[0]-erf_start_lo)/(abs(erf_start_lo)*double(0.35))) + double(0.5) - erf_offset_lo; //subtract value of erf at start of sponge to start xi_b at zero
                        // xi_b            = pow((x[0]-d_special_source_box_lo[0])/(domain_xlo[0]-d_special_source_box_lo[0]),3.0);
                    }
                    else // right-side
                    {
                        const double ksi_1 = (x[0] + shift)/delta;
    	                const double ksi_2 = (x[0] - shift)/delta;

    	                const double p = p_i + (g*rho_1*x[0]) + ((0.5*g*delta)*(rho_2 - rho_1)*((ksi_1*erf(ksi_1)) + ((exp(-pow(ksi_1,2.0)))/sqrt(M_PI)))) +
                            ((0.5*g*delta)*(rho_1 - rho_2)*((ksi_2*erf(ksi_2)) + ((exp(-pow(ksi_2,2.0)))/sqrt(M_PI)))) +
  			                (0.5*g*x[0]*(rho_3 - rho_1)) + ((0.5*g*delta*(rho_3 - rho_1))*((exp(-pow(ksi_1,2.0)))/sqrt(M_PI)));

                        const double p_L = p_i + (g*rho_1*x_L) + ((0.5*g*delta)*(rho_2 - rho_1)*((ksi_1*erf(ksi_1)) + ((exp(-pow(ksi_1,2.0)))/sqrt(M_PI)))) +
                            ((0.5*g*delta)*(rho_1 - rho_2)*((ksi_2*erf(ksi_2)) + ((exp(-pow(ksi_2,2.0)))/sqrt(M_PI)))) +
			                (0.5*g*x_L*(rho_3 - rho_1)) + ((0.5*g*delta*(rho_3 - rho_1))*((exp(-pow(ksi_1,2.0)))/sqrt(M_PI)));

                        const double p_R = p_i + (g*rho_1*x_R) + ((0.5*g*delta)*(rho_2 - rho_1)*((ksi_1*erf(ksi_1)) + ((exp(-pow(ksi_1,2.0)))/sqrt(M_PI)))) +
                            ((0.5*g*delta)*(rho_1 - rho_2)*((ksi_2*erf(ksi_2)) + ((exp(-pow(ksi_2,2.0)))/sqrt(M_PI)))) +
			                (0.5*g*x_R*(rho_3 - rho_1)) + ((0.5*g*delta*(rho_3 - rho_1))*((exp(-pow(ksi_1,2.0)))/sqrt(M_PI)));

                        rho_Y_0_ref = rho_1;

                        rho_Y_1_ref = 0.0;
                        rho_Y_2_ref = 0.0;
                    
                        const double u    = 0.0;
                        const double v    = 0.0;
                    
                        rho_u_ref         = rho_1*u;
                        rho_v_ref         = rho_1*v;
                        E_ref             = p/(gamma - double(1)) + double(1)/double(2)*rho_1*(u*u + v*v);
                        E_ref_L           = p_L/(gamma - double(1)) + double(1)/double(2)*rho_1*(u*u + v*v);
                        E_ref_R           = p_R/(gamma - double(1)) + double(1)/double(2)*rho_1*(u*u + v*v);


                        sponge_rate_tot = (pow((p/rho_3),0.5))*sponge_rate;
			            const double erf_start_hi = double(0.75)*(domain_xhi[0]-d_special_source_box_hi[0]) + d_special_source_box_hi[0]; //center of erf is 3/4 of the way into sponge
			            const double erf_offset_hi = double(0.5) * erf((d_special_source_box_hi[0]-erf_start_hi)/(abs(erf_start_hi)*double(0.35))) + double(0.5); //value of erf at start of sponge
                        xi_b            = double(0.5) * erf((x[0]-erf_start_hi)/(abs(erf_start_hi)*double(0.35))) + double(0.5) - erf_offset_hi; //subtract value of erf at start of sponge to start xi_b at zero

                    }

                    const double rho_Y_0_p = rho_Y_0[idx_cons_var] - rho_Y_0_ref;
                    const double rho_Y_1_p = rho_Y_1[idx_cons_var] - rho_Y_1_ref;
                    const double rho_Y_2_p = rho_Y_1[idx_cons_var] - rho_Y_2_ref;
                    const double rho_u_p   = rho_u[idx_cons_var]   - rho_u_ref;
                    const double rho_v_p   = rho_v[idx_cons_var]   - rho_v_ref;
                    const double E_p       = E[idx_cons_var]       - E_ref;
                   
	                const double rho_Y_0_p_L = rho_Y_0[idx_cons_var_L] - rho_Y_0_ref;
                    const double rho_Y_1_p_L = rho_Y_1[idx_cons_var_L] - rho_Y_1_ref;
                    const double rho_Y_2_p_L = rho_Y_2[idx_cons_var_L] - rho_Y_2_ref;
                    const double rho_u_p_L   = rho_u[idx_cons_var_L]   - rho_u_ref;
                    const double rho_v_p_L   = rho_v[idx_cons_var_L]   - rho_v_ref;
                    const double E_p_L       = E[idx_cons_var_L]       - E_ref_L;

           		    const double rho_Y_0_p_R = rho_Y_0[idx_cons_var_R] - rho_Y_0_ref;
                    const double rho_Y_1_p_R = rho_Y_1[idx_cons_var_R] - rho_Y_1_ref;
                    const double rho_Y_2_p_R = rho_Y_2[idx_cons_var_R] - rho_Y_2_ref;
                    const double rho_u_p_R   = rho_u[idx_cons_var_R]   - rho_u_ref;
                    const double rho_v_p_R   = rho_v[idx_cons_var_R]   - rho_v_ref;
                    const double E_p_R       = E[idx_cons_var_R]       - E_ref_R;
      		    
                    const double rho_Y_0_p_B = rho_Y_0[idx_cons_var_B] - rho_Y_0_ref;
                    const double rho_Y_1_p_B = rho_Y_1[idx_cons_var_B] - rho_Y_1_ref;
                    const double rho_Y_2_p_B = rho_Y_2[idx_cons_var_B] - rho_Y_2_ref;
                    const double rho_u_p_B   = rho_u[idx_cons_var_B]   - rho_u_ref;
                    const double rho_v_p_B   = rho_v[idx_cons_var_B]   - rho_v_ref;
                    const double E_p_B       = E[idx_cons_var_B]       - E_ref;

               		const double rho_Y_0_p_T = rho_Y_0[idx_cons_var_T] - rho_Y_0_ref;
                    const double rho_Y_1_p_T = rho_Y_1[idx_cons_var_T] - rho_Y_1_ref;
                    const double rho_Y_2_p_T = rho_Y_2[idx_cons_var_T] - rho_Y_2_ref;
                    const double rho_u_p_T   = rho_u[idx_cons_var_T]   - rho_u_ref;
                    const double rho_v_p_T   = rho_v[idx_cons_var_T]   - rho_v_ref;
                    const double E_p_T       = E[idx_cons_var_T]       - E_ref;


                    double sponge_nu       = (dx[0]*dx[0]*dx[1]*dx[1])/(dx[0]*dx[0]+dx[1]*dx[1])/2.0/dt*0.8;

                    S[0][idx_source] -= dt*sponge_rate_tot*xi_b*rho_Y_0_p;
                    S[1][idx_source] -= dt*sponge_rate_tot*xi_b*rho_Y_1_p;
                    S[2][idx_source] -= dt*sponge_rate_tot*xi_b*rho_Y_2_p;
                    S[3][idx_source] -= dt*sponge_rate_tot*xi_b*rho_u_p;
                    S[4][idx_source] -= dt*sponge_rate_tot*xi_b*rho_v_p;
                    S[5][idx_source] -= dt*sponge_rate_tot*xi_b*E_p;

    		        S[0][idx_source] += dt*sponge_nu*xi_b*( ((rho_Y_0_p_L - 2.0*rho_Y_0_p + rho_Y_0_p_R)/pow(dx[0], 2.0)) + ((rho_Y_0_p_B - 2.0*rho_Y_0_p + rho_Y_0_p_T)/pow(dx[1], 2.0)) );
                    S[1][idx_source] += dt*sponge_nu*xi_b*( ((rho_Y_1_p_L - 2.0*rho_Y_1_p + rho_Y_1_p_R)/pow(dx[0], 2.0)) + ((rho_Y_1_p_B - 2.0*rho_Y_1_p + rho_Y_1_p_T)/pow(dx[1], 2.0)) );
                    S[2][idx_source] += dt*sponge_nu*xi_b*( ((rho_Y_2_p_L - 2.0*rho_Y_2_p + rho_Y_2_p_R)/pow(dx[0], 2.0)) + ((rho_Y_2_p_B - 2.0*rho_Y_2_p + rho_Y_2_p_T)/pow(dx[1], 2.0)) );
                    S[3][idx_source] += dt*sponge_nu*xi_b*( ((rho_u_p_L - 2.0*rho_u_p + rho_u_p_R)/pow(dx[0], 2.0)) + ((rho_u_p_B - 2.0*rho_u_p + rho_u_p_T)/pow(dx[1], 2.0)) );
                    S[4][idx_source] += dt*sponge_nu*xi_b*( ((rho_v_p_L - 2.0*rho_v_p + rho_v_p_R)/pow(dx[0], 2.0)) + ((rho_v_p_B - 2.0*rho_v_p + rho_v_p_T)/pow(dx[1], 2.0)) );
                    S[5][idx_source] += dt*sponge_nu*xi_b*( ((E_p_L - 2.0*E_p + E_p_R)/pow(dx[0], 2.0)) + ((E_p_B - 2.0*E_p + E_p_T)/pow(dx[1], 2.0)) );
                }
            }
        }

    }



}


void
FlowModelSpecialSourceTerms::putToRestart(const HAMERS_SHARED_PTR<tbox::Database>& restart_source_terms_db)
{
    putToRestartBase(restart_source_terms_db);
    
    double sponge_rate = double(0);
    if (d_source_terms_db->keyExists("sponge_rate"))
    {
        sponge_rate = d_source_terms_db->getDouble("sponge_rate");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'sponge_rate' found in data for source terms."
            << std::endl);
    }
    
    restart_source_terms_db->putDouble("sponge_rate", sponge_rate);

    std::vector<double> species_mass;
    if (d_source_terms_db->keyExists("species_mass"))
    {
        d_source_terms_db->getVector("species_mass", species_mass);
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'species_mass' found in data for source terms."
            << std::endl);
    }
    
    restart_source_terms_db->putVector("species_mass", species_mass);
}
