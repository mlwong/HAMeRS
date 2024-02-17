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
    
    if ((d_project_name != "3D discontinuous Rayleigh-Taylor instability") && 
        (d_project_name != "3D smooth Rayleigh-Taylor instability") && 
        (d_project_name != "3D smooth multi-mode Rayleigh-Taylor instability")) 
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Can only initialize data for 'project_name' = '3D discontinuous Rayleigh-Taylor instability' or "
            << "'3D smooth Rayleigh-Taylor instability' or "
            << "'3D smooth multi-mode Rayleigh-Taylor instability'!\n"
            << "'project_name' = '"
            << d_project_name
            << "' is given."
            << std::endl);
    }
    
    if (d_dim != tbox::Dimension(3))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Dimension of problem should be 3!"
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
    
    TBOX_ASSERT(d_source_terms_db->keyExists("sponge_nu_fact"));
    
    double sponge_nu_fact = double(0); // 1.0 could be an input parameter for sponge_nu_fact
    if (d_source_terms_db->keyExists("sponge_nu_fact"))
    {
        sponge_nu_fact = d_source_terms_db->getDouble("sponge_nu_fact");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'sponge_nu_fact' found in data for source terms."
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
     * Initialize data for a 3D Rayleigh-Taylor instability problem (At = 0.04, M = 0.3).
     */
    
    HAMERS_SHARED_PTR<pdat::CellData<double> > partial_density = conservative_variables[0];
    HAMERS_SHARED_PTR<pdat::CellData<double> > momentum        = conservative_variables[1];
    HAMERS_SHARED_PTR<pdat::CellData<double> > total_energy    = conservative_variables[2];
    
    double* rho_Y_0 = partial_density->getPointer(0);
    double* rho_Y_1 = partial_density->getPointer(1);
    double* rho_u   = momentum->getPointer(0);
    double* rho_v   = momentum->getPointer(1);
    double* rho_w   = momentum->getPointer(2);
    double* E       = total_energy->getPointer(0);
    
    const double gamma = double(7)/double(5); // assume both gases have the same ratio of specific heat ratios
    
    double lambda = 701.53278340668; // wavelength of single-mode perturbation
    double eta_0  = 0.02*lambda;
    
    const double p_i = 100000.0; // interface pressure
    const double T_0 = 300.0;    // background temperature
    
    const double sponge_nu = (dx[0]*dx[0]*dx[1]*dx[1]*dx[2]*dx[2])/(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2])/(3.0*dt)*sponge_nu_fact; // It was (2.0/dt) for 2D, ask this.
    
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
    const double W_1 = species_mass[0]; // species mass 1
    const double W_2 = species_mass[1]; // species mass 2
    
    
    const double R_u = 8.31446261815324; // universal gas constant
    const double R_1 = R_u/W_1;          // gas constant of heavier gas
    const double R_2 = R_u/W_2;          // gas constant of lighter gas
    
    const double* const domain_xlo = d_grid_geometry->getXLower();
    const double* const domain_xhi = d_grid_geometry->getXUpper();
    
    if (d_project_name == "3D discontinuous Rayleigh-Taylor instability")
    {
        // for (int j = 0; j < patch_dims[1]; j++) // Where should insert k
        // {
        //     for (int i = 0; i < patch_dims[0]; i++)
        //     {
        //         // Compute the linear indices.
        //         const int idx_source = (i + num_ghosts_source[0]) +
        //             (j + num_ghosts_source[1])*ghostcell_dims_source[0];
                
        //         const int idx_cons_var = (i + num_ghosts_cons_var[0]) +
        //             (j + num_ghosts_cons_var[1])*ghostcell_dims_cons_var[0];
                
        //         const int idx_cons_var_L = (i - 1 + num_ghosts_cons_var[0]) +
        //             (j + num_ghosts_cons_var[1])*ghostcell_dims_cons_var[0];
                
        //         const int idx_cons_var_R = (i + 1 + num_ghosts_cons_var[0]) +
        //             (j + num_ghosts_cons_var[1])*ghostcell_dims_cons_var[0];
                
        //         const int idx_cons_var_B = (i + num_ghosts_cons_var[0]) +
        //             (j - 1 + num_ghosts_cons_var[1])*ghostcell_dims_cons_var[0];
                
        //         const int idx_cons_var_T = (i + num_ghosts_cons_var[0]) +
        //             (j + 1 + num_ghosts_cons_var[1])*ghostcell_dims_cons_var[0];
                
        //         // Compute the coordinates.
        //         double x[2];
        //         x[0] = patch_xlo[0] + (double(i) + double(1)/double(2))*dx[0];
        //         x[1] = patch_xlo[1] + (double(j) + double(1)/double(2))*dx[1];
                
        //         double x_L;
        //         x_L = patch_xlo[0] + (double(i) - 1.0 + double(1)/double(2))*dx[0];
                
        //         double x_R;
        //         x_R = patch_xlo[0] + (double(i) + 1.0  + double(1)/double(2))*dx[0];
                    
        //         // Check whether it is outside the special source box.
        //         if (!(x[0] >= d_special_source_box_lo[0] && x[0] <= d_special_source_box_hi[0]))
        //         {
        //             double rho_Y_0_ref;
        //             double rho_Y_1_ref;
        //             double rho_u_ref;
        //             double rho_v_ref;
        //             double rho_w_ref;
        //             double E_ref;
        //             double sponge_rate_tot;
        //             double xi_b;
                    
        //             double rho_Y_0_ref_L;
        //             double rho_Y_1_ref_L;
        //             double rho_u_ref_L;
        //             double rho_v_ref_L;
        //             double rho_w_ref_L;
        //             double E_ref_L;
                    
        //             double rho_Y_0_ref_R;
        //             double rho_Y_1_ref_R;
        //             double rho_u_ref_R;
        //             double rho_v_ref_R;
        //             double rho_w_ref_R;
        //             double E_ref_R;
                    
        //             double rho_Y_0_ref_B;
        //             double rho_Y_1_ref_B;
        //             double rho_u_ref_B;
        //             double rho_v_ref_B;
        //             double rho_w_ref_B;
        //             double E_ref_B;
                    
        //             double rho_Y_0_ref_T;
        //             double rho_Y_1_ref_T;
        //             double rho_u_ref_T;
        //             double rho_v_ref_T;
        //             double rho_w_ref_T;
        //             double E_ref_T;
                    
        //             if (x[0] < 0.0) // left-side
        //             {
        //                 const double rho_ref = p_i/(R_1*T_0)*exp((g*x[0])/(R_1*T_0));
        //                 rho_Y_0_ref = rho_ref;
        //                 const double rho_ref_L = p_i/(R_1*T_0)*exp((g*x_L)/(R_1*T_0));
        //                 rho_Y_0_ref_L = rho_ref_L;
        //                 const double rho_ref_R = p_i/(R_1*T_0)*exp((g*x_R)/(R_1*T_0));
        //                 rho_Y_0_ref_R = rho_ref_R;
        //                 rho_Y_0_ref_B = rho_ref; // they are at the same height with rho_ref
        //                 rho_Y_0_ref_T = rho_ref; // they are at the same height with rho_ref
        //                 rho_Y_1_ref   = 0.0;
        //                 rho_Y_1_ref_L = 0.0;
        //                 rho_Y_1_ref_R = 0.0;
        //                 rho_Y_1_ref_B = 0.0;
        //                 rho_Y_1_ref_T = 0.0;
                        
        //                 const double p_ref   = p_i*exp((g*x[0])/(R_1*T_0));
        //                 const double p_ref_L = p_i*exp((g*x_L)/(R_1*T_0));
        //                 const double p_ref_R = p_i*exp((g*x_R)/(R_1*T_0));
        //                 const double p_ref_B = p_ref; // they are at the same height with p_ref
        //                 const double p_ref_T = p_ref; // they are at the same height with p_ref
                        
        //                 const double u_ref = 0.0;
        //                 const double v_ref = 0.0;
        //                 const double w_ref = 0.0;
                        
        //                 rho_u_ref         = rho_ref*u_ref;
        //                 rho_u_ref_L       = rho_ref_L*u_ref;
        //                 rho_u_ref_R       = rho_ref_R*u_ref;
        //                 rho_u_ref_B       = rho_ref*u_ref; // they are at the same height with rho_ref
        //                 rho_u_ref_T       = rho_ref*u_ref; // they are at the same height with rho_ref
                        
        //                 rho_v_ref         = rho_ref*v_ref;
        //                 rho_v_ref_L       = rho_ref_L*v_ref;
        //                 rho_v_ref_R       = rho_ref_R*v_ref;
        //                 rho_v_ref_B       = rho_ref*v_ref; // they are at the same height with rho_ref
        //                 rho_v_ref_T       = rho_ref*v_ref; // they are at the same height with rho_ref

        //                 rho_w_ref         = rho_ref*w_ref;
        //                 rho_w_ref_L       = rho_ref_L*w_ref;
        //                 rho_w_ref_R       = rho_ref_R*w_ref;
        //                 rho_w_ref_B       = rho_ref*w_ref; // they are at the same height with rho_ref
        //                 rho_w_ref_T       = rho_ref*w_ref; // they are at the same height with rho_ref
                        
        //                 E_ref             = p_ref/(gamma - double(1)) + double(1)/double(2)*rho_ref*(u_ref*u_ref + v_ref*v_ref + w_ref*w_ref);
        //                 E_ref_L           = p_ref_L/(gamma - double(1)) + double(1)/double(2)*rho_ref_L*(u_ref*u_ref + v_ref*v_ref + w_ref*w_ref);
        //                 E_ref_R           = p_ref_R/(gamma - double(1)) + double(1)/double(2)*rho_ref_R*(u_ref*u_ref + v_ref*v_ref + w_ref*w_ref);
        //                 E_ref_B           = E_ref; // they are at the same height with E_ref
        //                 E_ref_T           = E_ref; // they are at the same height with E_ref
                        
        //                 sponge_rate_tot   = (pow((gamma*p_ref/rho_ref),0.5))*sponge_rate;
        //                 const double erf_start_lo  = double(0.75)*(domain_xlo[0]-d_special_source_box_lo[0]) + d_special_source_box_lo[0]; //center of erf is 3/4 of the way into sponge
        //                 const double erf_offset_lo = double(-0.5) * erf((d_special_source_box_lo[0]-erf_start_lo)/(abs(erf_start_lo)*double(0.35))) + double(0.5); //value of erf at start of sponge
        //                 xi_b        = double(-0.5) * erf((x[0]-erf_start_lo)/(abs(erf_start_lo)*double(0.35))) + double(0.5) - erf_offset_lo; //subtract value of erf at start of sponge to start xi_b at zero
        //                 xi_b        = 1.0;
        //                 // xi_b            = pow((x[0]-d_special_source_box_lo[0])/(domain_xlo[0]-d_special_source_box_lo[0]),3.0);
        //             }
        //             else // right-side
        //             {
        //                 const double rho_ref   = p_i/(R_2*T_0)*exp((g*x[0])/(R_2*T_0));
        //                 rho_Y_1_ref            = rho_ref;
        //                 const double rho_ref_L = p_i/(R_2*T_0)*exp((g*x_L)/(R_2*T_0));
        //                 rho_Y_1_ref_L          = rho_ref_L;
        //                 const double rho_ref_R = p_i/(R_2*T_0)*exp((g*x_R)/(R_2*T_0));
        //                 rho_Y_1_ref_R          = rho_ref_R;
                        
        //                 rho_Y_1_ref_B          = rho_ref; // they are at the same height with rho_ref
        //                 rho_Y_1_ref_T          = rho_ref; // they are at the same height with rho_ref
                        
        //                 rho_Y_0_ref   = 0.0;
        //                 rho_Y_0_ref_L = 0.0;
        //                 rho_Y_0_ref_R = 0.0;
        //                 rho_Y_0_ref_B = 0.0;
        //                 rho_Y_0_ref_T = 0.0;
                        
        //                 const double p_ref = p_i*exp((g*x[0])/(R_2*T_0));
        //                 const double p_ref_L = p_i*exp((g*x_L)/(R_2*T_0));
        //                 const double p_ref_R = p_i*exp((g*x_R)/(R_2*T_0));
        //                 const double p_ref_B = p_ref; // they are at the same height with p_ref
        //                 const double p_ref_T = p_ref; // they are at the same height with p_ref
                        
        //                 const double u_ref = 0.0;
        //                 const double v_ref = 0.0;
        //                 const double w_ref = 0.0;
                        
        //                 rho_u_ref         = rho_ref*u_ref;
        //                 rho_u_ref_L       = rho_ref_L*u_ref;
        //                 rho_u_ref_R       = rho_ref_R*u_ref;
        //                 rho_u_ref_B       = rho_ref*u_ref; // they are at the same height with rho_ref
        //                 rho_u_ref_T       = rho_ref*u_ref; // they are at the same height with rho_ref
        //                 rho_v_ref         = rho_ref*v_ref;
        //                 rho_v_ref_L       = rho_ref_L*v_ref;
        //                 rho_v_ref_R       = rho_ref_R*v_ref;
        //                 rho_v_ref_B       = rho_ref*v_ref; // they are at the same height with rho_ref
        //                 rho_v_ref_T       = rho_ref*v_ref; // they are at the same height with rho_ref
        //                 rho_w_ref         = rho_ref*w_ref;
        //                 rho_w_ref_L       = rho_ref_L*w_ref;
        //                 rho_w_ref_R       = rho_ref_R*w_ref;
        //                 rho_w_ref_B       = rho_ref*w_ref; // they are at the same height with rho_ref
        //                 rho_w_ref_T       = rho_ref*w_ref; // they are at the same height with rho_ref
                        
        //                 E_ref             = p_ref/(gamma - double(1)) + double(1)/double(2)*rho_ref*(u_ref*u_ref + v_ref*v_ref + w_ref*w_ref);
        //                 E_ref_L           = p_ref_L/(gamma - double(1)) + double(1)/double(2)*rho_ref_L*(u_ref*u_ref + v_ref*v_ref + w_ref*w_ref);
        //                 E_ref_R           = p_ref_R/(gamma - double(1)) + double(1)/double(2)*rho_ref_R*(u_ref*u_ref + v_ref*v_ref + w_ref*w_ref);
        //                 E_ref_B           = E_ref; // they are at the same height with E_ref
        //                 E_ref_T           = p_ref; // they are at the same height with E_ref
                    
        //                 sponge_rate_tot = (pow((gamma*p_ref/rho_ref),0.5))*sponge_rate;
        //                 const double erf_start_hi = double(0.75)*(domain_xhi[0]-d_special_source_box_hi[0]) + d_special_source_box_hi[0]; //center of erf is 3/4 of the way into sponge
        //                 const double erf_offset_hi = double(0.5) * erf((d_special_source_box_hi[0]-erf_start_hi)/(abs(erf_start_hi)*double(0.35))) + double(0.5); //value of erf at start of sponge
        //                 xi_b            = double(0.5) * erf((x[0]-erf_start_hi)/(abs(erf_start_hi)*double(0.35))) + double(0.5) - erf_offset_hi; //subtract value of erf at start of sponge to start xi_b at zero
        //                 xi_b = 1.0;
        //                 // xi_b            = pow((x[0]-d_special_source_box_hi[0])/(domain_xhi[0]-d_special_source_box_hi[0]),3.0);
        //             }
                    
        //             const double rho_Y_0_p = rho_Y_0[idx_cons_var] - rho_Y_0_ref;
        //             const double rho_Y_1_p = rho_Y_1[idx_cons_var] - rho_Y_1_ref;
        //             const double rho_u_p   = rho_u[idx_cons_var]   - rho_u_ref;
        //             const double rho_v_p   = rho_v[idx_cons_var]   - rho_v_ref;
        //             const double rho_w_p   = rho_w[idx_cons_var]   - rho_w_ref;
        //             const double E_p       = E[idx_cons_var]       - E_ref;
                   
        //             const double rho_Y_0_p_L = rho_Y_0[idx_cons_var_L] - rho_Y_0_ref_L;
        //             const double rho_Y_1_p_L = rho_Y_1[idx_cons_var_L] - rho_Y_1_ref_L;
        //             const double rho_u_p_L   = rho_u[idx_cons_var_L]   - rho_u_ref_L;
        //             const double rho_v_p_L   = rho_v[idx_cons_var_L]   - rho_v_ref_L;
        //             const double rho_w_p_L   = rho_w[idx_cons_var_L]   - rho_w_ref_L;
        //             const double E_p_L       = E[idx_cons_var_L]       - E_ref_L;
                    
        //             const double rho_Y_0_p_R = rho_Y_0[idx_cons_var_R] - rho_Y_0_ref_R;
        //             const double rho_Y_1_p_R = rho_Y_1[idx_cons_var_R] - rho_Y_1_ref_R;
        //             const double rho_u_p_R   = rho_u[idx_cons_var_R]   - rho_u_ref_R;
        //             const double rho_v_p_R   = rho_v[idx_cons_var_R]   - rho_v_ref_R;
        //             const double rho_w_p_R   = rho_w[idx_cons_var_R]   - rho_w_ref_R;
        //             const double E_p_R       = E[idx_cons_var_R]       - E_ref_R;
                    
        //             const double rho_Y_0_p_B = rho_Y_0[idx_cons_var_B] - rho_Y_0_ref_B;
        //             const double rho_Y_1_p_B = rho_Y_1[idx_cons_var_B] - rho_Y_1_ref_B;
        //             const double rho_u_p_B   = rho_u[idx_cons_var_B]   - rho_u_ref_B;
        //             const double rho_v_p_B   = rho_v[idx_cons_var_B]   - rho_v_ref_B;
        //             const double rho_w_p_B   = rho_w[idx_cons_var_B]   - rho_w_ref_B;
        //             const double E_p_B       = E[idx_cons_var_B]       - E_ref_B;
                    
        //             const double rho_Y_0_p_T = rho_Y_0[idx_cons_var_T] - rho_Y_0_ref_T;
        //             const double rho_Y_1_p_T = rho_Y_1[idx_cons_var_T] - rho_Y_1_ref_T;
        //             const double rho_u_p_T   = rho_u[idx_cons_var_T]   - rho_u_ref_T;
        //             const double rho_v_p_T   = rho_v[idx_cons_var_T]   - rho_v_ref_T;
        //             const double rho_w_p_T   = rho_w[idx_cons_var_T]   - rho_w_ref_T;
        //             const double E_p_T       = E[idx_cons_var_T]       - E_ref_T;
                    
        //             S[0][idx_source] -= dt*sponge_rate_tot*xi_b*rho_Y_0_p;
        //             S[1][idx_source] -= dt*sponge_rate_tot*xi_b*rho_Y_1_p;
        //             S[2][idx_source] -= dt*sponge_rate_tot*xi_b*rho_u_p;
        //             S[3][idx_source] -= dt*sponge_rate_tot*xi_b*rho_v_p;
        //             S[4][idx_source] -= dt*sponge_rate_tot*xi_b*rho_w_p;
        //             S[5][idx_source] -= dt*sponge_rate_tot*xi_b*E_p;
        //         }
        //     }
        // }
    }
    else if (d_project_name == "3D smooth Rayleigh-Taylor instability")
    {
        const double delta = 0.02*lambda; // characteristic length of interface.
        const double shift = 0.0; // location of interface.
        
        // Discretize the domain in x-direction for the approximated integral.
        const int integral_N_x = 10000;
        
        std::vector<double>& integral_vector = d_special_source_vector; // Reference to the class member vector that stores the numerical integral.
        
        if (integral_vector.empty())
        {
            integral_vector.resize(integral_N_x + 3);
            std::ifstream f_in;
            std::string integral_filename = "integral.dat";
            f_in.open(integral_filename, std::ios::in | std::ios::binary);
            f_in.read((char*)&integral_vector[0], sizeof(double)*integral_vector.size());
            f_in.close();
        }
        
        const double x_domain_lo = integral_vector[integral_N_x + 0];
        const double x_domain_hi = integral_vector[integral_N_x + 1];
        const double dx_uniform  = integral_vector[integral_N_x + 2];
        
        const int patch_dims_0 = patch_dims[0];
        const int patch_dims_1 = patch_dims[1];
        const int patch_dims_2 = patch_dims[2];
        const int num_ghosts_source_0 = num_ghosts_source[0];
        const int num_ghosts_source_1 = num_ghosts_source[1];
        const int num_ghosts_source_2 = num_ghosts_source[2];
        const int ghostcell_dims_source_0 = ghostcell_dims_source[0];
        const int num_ghosts_cons_var_0 = num_ghosts_cons_var[0];
        const int num_ghosts_cons_var_1 = num_ghosts_cons_var[1];
        const int num_ghosts_cons_var_2 = num_ghosts_cons_var[2];
        const int ghostcell_dims_cons_var_0 = ghostcell_dims_cons_var[0];
        
        for (int k = 0; k < patch_dims_2; k++)
        {
            for (int j = 0; j < patch_dims_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = 0; i < patch_dims_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_source = (i + num_ghosts_source[0]) +
                        (j + num_ghosts_source[1])*ghostcell_dims_source[0] + 
                        (k + num_ghosts_source[2])*ghostcell_dims_source[0]*ghostcell_dims_source[1];

                    const int idx_cons_var = (i + num_ghosts_cons_var[0]) +
                        (j + num_ghosts_cons_var[1])*ghostcell_dims_cons_var[0] +
                        (k + num_ghosts_cons_var[2])*ghostcell_dims_cons_var[0]*ghostcell_dims_cons_var[1];

                    

                    // Compute the coordinates.
                    double x[3];
                    x[0] = patch_xlo[0] + (double(i) + double(1)/double(2))*dx[0];
                    x[1] = patch_xlo[1] + (double(j) + double(1)/double(2))*dx[1];
                    x[2] = patch_xlo[2] + (double(k) + double(1)/double(2))*dx[2];

                    const double eta = eta_0*cos(2.0*M_PI/lambda*x[1])*cos(2.0*M_PI/lambda*x[2]);
                        
                    double X_2_H = 0.5*(1.0 + erf((x[0] - eta - shift)/delta)); // mass fraction of second species (Y_2)
                    const double R_H   = R_1*(1.0 - X_2_H) + X_2_H*R_2;
        
                    if (!(x[0] >= d_special_source_box_lo[0] && x[0] <= d_special_source_box_hi[0]))
                    {
                        double rho_Y_0_ref;
                        double rho_Y_1_ref;
                        double rho_u_ref;
                        double rho_v_ref;
                        double rho_w_ref;
                        double E_ref;
                        double sponge_rate_tot;
                        double xi_b;
                        
                        const int N_int = 100000; // number of numerical quadrature points
                        const double dx_p = (x[0] - shift)/(N_int - 1.0);
                    
                        double integral = 0.0;
                        double p_H = 0.0;
                        double rho_H = 0.0;
                        
                        
                        // Compute the integral with linear interpolation.
                        
                        const int idx_integral_vector_lo = int(std::floor((x[0] - x_domain_lo)/dx_uniform));
                        const int idx_integral_vector_hi = idx_integral_vector_lo + 1;
                        
                        const double x_integral_vector_lo = double(idx_integral_vector_lo)*dx_uniform + 0.5*dx_uniform + x_domain_lo;
                        const double x_integral_vector_hi = x_integral_vector_lo + dx_uniform;
                        
                        const double weight_lo = ( x_integral_vector_hi - x[0])/dx_uniform;
                        const double weight_hi = (-x_integral_vector_lo + x[0])/dx_uniform;
                        
                        integral = weight_lo*integral_vector[idx_integral_vector_lo] + weight_hi*integral_vector[idx_integral_vector_hi];
                        
                        //for (int ii = 0; ii < N_int; ii++)
                        //{
                        //    const double x_p = shift + ii*dx_p;  //Bug fixed 3.22.2023
                        //    integral += 1.0/(0.5*(R_2 - R_1)*erf((x_p - shift)/(delta)) + 0.5*(R_1 + R_2))*dx_p;
                        //}
                        
                        
                        p_H = p_i*exp(g/T_0*integral);
                        rho_H = p_H/(R_H*T_0);
                        
                        double rho_ref, p_ref;
                        rho_ref = rho_H;
                        p_ref   = p_H;
                        
                        rho_Y_0_ref = rho_ref*(1.0 - X_2_H);
                        rho_Y_1_ref = rho_ref*X_2_H;
                        
                        const double u_ref = 0.0;
                        const double v_ref = 0.0;
                        const double w_ref = 0.0;
                        
                        rho_u_ref = rho_ref*u_ref;
                        rho_v_ref = rho_ref*v_ref;
                        rho_w_ref = rho_ref*w_ref;
                        E_ref     = p_ref/(gamma - double(1)) + double(1)/double(2)*rho_ref*(u_ref*u_ref + v_ref*v_ref + w_ref*w_ref);
                        
                        sponge_rate_tot = sponge_rate*300.0; //300.0 is roughly (gamma*p_i/rho_i)^.05 OLD:(pow((gamma*p_ref/rho_ref),0.5))*sponge_rate;
                        if(x[0] < 0.0)
                        {
                            const double erf_start_lo  = double(0.75)*(domain_xlo[0]-d_special_source_box_lo[0]) + d_special_source_box_lo[0]; //center of erf is 3/4 of the way into sponge
                            const double erf_offset_lo = double(-0.5) * erf((d_special_source_box_lo[0]-erf_start_lo)/(abs(erf_start_lo)*double(0.1))) + double(0.5); //value of erf at start of sponge
                            xi_b        = double(-0.5) * erf((x[0]-erf_start_lo)/(abs(erf_start_lo)*double(0.1))) + double(0.5) - erf_offset_lo; //subtract value of erf at start of sponge to start xi_b at zero
                            // xi_b        = 1.0;
                        }
                        else
                        {
                            const double erf_start_hi = double(0.75)*(domain_xhi[0]-d_special_source_box_hi[0]) + d_special_source_box_hi[0]; //center of erf is 3/4 of the way into sponge
                            const double erf_offset_hi = double(0.5) * erf((d_special_source_box_hi[0]-erf_start_hi)/(abs(erf_start_hi)*double(0.1))) + double(0.5); //value of erf at start of sponge
                            xi_b            = double(0.5) * erf((x[0]-erf_start_hi)/(abs(erf_start_hi)*double(0.1))) + double(0.5) - erf_offset_hi; //subtract value of erf at start of sponge to start xi_b at zero
                            // xi_b = 1.0;
                        }

                        const double rho_Y_0_p = rho_Y_0[idx_cons_var] - rho_Y_0_ref;
                        const double rho_Y_1_p = rho_Y_1[idx_cons_var] - rho_Y_1_ref;
                        const double rho_u_p   = rho_u[idx_cons_var]   - rho_u_ref;
                        const double rho_v_p   = rho_v[idx_cons_var]   - rho_v_ref;
                        const double rho_w_p   = rho_w[idx_cons_var]   - rho_w_ref;
                        const double E_p       = E[idx_cons_var]       - E_ref;
                    
                        S[0][idx_source] -= dt*sponge_rate_tot*xi_b*rho_Y_0_p;
                        S[1][idx_source] -= dt*sponge_rate_tot*xi_b*rho_Y_1_p;
                        S[2][idx_source] -= dt*sponge_rate_tot*xi_b*rho_u_p;
                        S[3][idx_source] -= dt*sponge_rate_tot*xi_b*rho_v_p;
                        S[4][idx_source] -= dt*sponge_rate_tot*xi_b*rho_w_p;
                        S[5][idx_source] -= dt*sponge_rate_tot*xi_b*E_p;
                    }
                }
            }
        }
    }
    else if (d_project_name == "3D smooth multi-mode Rayleigh-Taylor instability")
    {
        lambda               = lambda/4.0;
        eta_0                = lambda*0.04;
        const double shift   = 0.0; // location of interface.

        const double delta   = 0.04*lambda; // characteristic length of interface
        const int    waven   = 16;          // dominant wave number
        const double width   = lambda; // domain size in y direction

        // Discretize the domain in x-direction for the approximated integral.
        const int integral_N_x = 10000;
        
        std::vector<double>& integral_vector = d_special_source_vector; // Reference to the class member vector that stores the numerical integral.
        
        if (integral_vector.empty())
        {
            integral_vector.resize(integral_N_x + 3);
            std::ifstream f_in;
            std::string integral_filename = "integral.dat";
            f_in.open(integral_filename, std::ios::in | std::ios::binary);
            f_in.read((char*)&integral_vector[0], sizeof(double)*integral_vector.size());
            f_in.close();
        }
        
        const double x_domain_lo = integral_vector[integral_N_x + 0];
        const double x_domain_hi = integral_vector[integral_N_x + 1];
        const double dx_uniform  = integral_vector[integral_N_x + 2];
        
        const int patch_dims_0 = patch_dims[0];
        const int patch_dims_1 = patch_dims[1];
        const int patch_dims_2 = patch_dims[2];
        const int num_ghosts_source_0 = num_ghosts_source[0];
        const int num_ghosts_source_1 = num_ghosts_source[1];
        const int num_ghosts_source_2 = num_ghosts_source[2];
        const int ghostcell_dims_source_0 = ghostcell_dims_source[0];
        const int num_ghosts_cons_var_0 = num_ghosts_cons_var[0];
        const int num_ghosts_cons_var_1 = num_ghosts_cons_var[1];
        const int num_ghosts_cons_var_2 = num_ghosts_cons_var[2];
        const int ghostcell_dims_cons_var_0 = ghostcell_dims_cons_var[0];

        // Seed for "random" phase shifts.
        int random_seed = 0;
        if (d_source_terms_db->keyExists("random_seed"))
        {
            random_seed = d_source_terms_db->getInteger("random_seed");
            if (random_seed < 0 || random_seed > 15)
            {
                TBOX_ERROR(d_object_name << ": "
                    << "'random_seed' should be in between 0 - 15."
                    << std::endl);
            }
        }
        else
        {
            TBOX_WARNING(d_object_name << ": "
                << "'random_seed' not given for '2D smooth multi-mode Rayleigh-Taylor instability'!"
                << "  'random_seed' = 0 is assumed."
                << std::endl);
        }
        
        double rmod[9]; // random seed
        switch (random_seed)
        {
            case 0:
            {
                rmod[0] = 6.031966614958411e+000;
                rmod[1] = 1.273017034173460e+000;
                rmod[2] = 5.934447177754063e+000;
                rmod[3] = 3.101658133166612e+000;
                rmod[4] = 2.294026034817427e+000;
                rmod[5] = 4.916046917518752e+000;
                rmod[6] = 0.571212135466553e+000;
                rmod[7] = 4.966766749458944e+000;
                rmod[8] = 5.027899324302027e+000;
                
                break;
            }
            case 1:
            {
                rmod[0] = 2.620226532717789200e+000;
                rmod[1] = 4.525932273597345700e+000;
                rmod[2] = 7.186381718527406600e-004;
                rmod[3] = 1.899611578242180700e+000;
                rmod[4] = 9.220944569241362700e-001;
                rmod[5] = 5.801805019369201700e-001;
                rmod[6] = 1.170307423440345900e+000;
                rmod[7] = 2.171222082895173200e+000;
                rmod[8] = 2.492963564452900500e+000;
                
                break;
            }
            case 2:
            {
                rmod[0] = 2.739436763143839700e+00;
                rmod[1] = 1.628993188915385800e-01;
                rmod[2] = 3.453631204915429600e+00;
                rmod[3] = 2.735211261185420500e+00;
                rmod[4] = 2.641248797687487200e+00;
                rmod[5] = 2.075554893781340400e+00;
                rmod[6] = 1.285845290520944100e+00;
                rmod[7] = 3.890994236937394200e+00;
                rmod[8] = 1.882785842859457300e+00;
                
                break;
            }
            case 3:
            {
                rmod[0] = 3.460765288681905800e+00;
                rmod[1] = 4.449423994385291800e+00;
                rmod[2] = 1.827808381326725400e+00;
                rmod[3] = 3.209624503479690600e+00;
                rmod[4] = 5.610551183647944900e+00;
                rmod[5] = 5.631575567313184600e+00;
                rmod[6] = 7.890757775039627400e-01;
                rmod[7] = 1.302145406935464500e+00;
                rmod[8] = 3.233779755813989700e-01;
                
                break;
            }
            case 4:
            {
                rmod[0] = 6.076027676094973600e+00;
                rmod[1] = 3.438361627635736700e+00;
                rmod[2] = 6.111556079054740700e+00;
                rmod[3] = 4.491321348791744100e+00;
                rmod[4] = 4.383959499105254800e+00;
                rmod[5] = 1.357730343666468900e+00;
                rmod[6] = 6.134113310024843300e+00;
                rmod[7] = 3.914584796145817400e-02;
                rmod[8] = 1.589535062303236700e+00;
                
                break;
            }
            case 5:
            {
                rmod[0] = 1.394824230885255200e+00;
                rmod[1] = 5.470972432660288700e+00;
                rmod[2] = 1.298854759541258500e+00;
                rmod[3] = 5.771802559770447900e+00;
                rmod[4] = 3.068778005297785300e+00;
                rmod[5] = 3.843700051147186600e+00;
                rmod[6] = 4.812340990490530300e+00;
                rmod[7] = 3.257316284380881800e+00;
                rmod[8] = 1.864852550667249300e+00;
                
                break;
            }
            case 6:
            {
                rmod[0] = 5.610005784868826100e+00;
                rmod[1] = 2.085890634948696300e+00;
                rmod[2] = 5.159934759824945000e+00;
                rmod[3] = 2.619876261158565800e-01;
                rmod[4] = 6.764268695934091400e-01;
                rmod[5] = 3.738822386827532500e+00;
                rmod[6] = 3.328940665618581800e+00;
                rmod[7] = 2.631444681644834500e+00;
                rmod[8] = 2.107429670466887600e+00;
                
                break;
            }
            case 7:
            {
                rmod[0] = 4.794591226104558700e-01;
                rmod[1] = 4.900374296196336100e+00;
                rmod[2] = 2.754606441521316700e+00;
                rmod[3] = 4.545665775603436200e+00;
                rmod[4] = 6.144889332352788000e+00;
                rmod[5] = 3.383469340939719400e+00;
                rmod[6] = 3.148632734395143500e+00;
                rmod[7] = 4.527106224916906400e-01;
                rmod[8] = 1.686651855650350300e+00;
                
                break;
            }
            case 8:
            {
                rmod[0] = 5.487918790480180500e+00;
                rmod[1] = 6.085520462042457400e+00;
                rmod[2] = 5.461310364152818200e+00;
                rmod[3] = 3.335464681414797900e+00;
                rmod[4] = 1.462275210911384800e+00;
                rmod[5] = 7.162079955536335100e-02;
                rmod[6] = 2.704715354294335400e+00;
                rmod[7] = 2.528048153710494200e+00;
                rmod[8] = 3.284061815477384600e+00;
                
                break;
            }
            case 9:
            {
                rmod[0] = 6.518273126904997100e-02;
                rmod[1] = 3.153371063435702800e+00;
                rmod[2] = 3.115035471112504800e+00;
                rmod[3] = 8.408757300236915400e-01;
                rmod[4] = 8.929102841152991600e-01;
                rmod[5] = 1.373244659450403700e+00;
                rmod[6] = 2.629564450717737100e+00;
                rmod[7] = 1.558865616070162600e+00;
                rmod[8] = 5.281623651287690200e-01;
                
                break;
            }
            case 10:
            {
                rmod[0] = 4.846350532897925100e+00;
                rmod[1] = 1.303883433103263400e-01;
                rmod[2] = 3.981329279609052500e+00;
                rmod[3] = 4.704873552725635100e+00;
                rmod[4] = 3.132211935225629200e+00;
                rmod[5] = 1.412438980302679600e+00;
                rmod[6] = 1.244465681755566800e+00;
                rmod[7] = 4.778555396547323800e+00;
                rmod[8] = 1.062554723574571100e+00;
                
                break;
            }
            case 11:
            {
                rmod[0] = 1.132667860480351500e+00;
                rmod[1] = 1.223665511688170900e-01;
                rmod[2] = 2.910487839707776500e+00;
                rmod[3] = 4.554894212576069600e+00;
                rmod[4] = 2.640217114369509700e+00;
                rmod[5] = 3.050028410914633200e+00;
                rmod[6] = 8.030422644949872200e-02;
                rmod[7] = 3.062246122248716100e+00;
                rmod[8] = 5.917545720207830800e+00;
                
                break;
            }
            case 12:
            {
                rmod[0] = 9.686337061529999300e-01;
                rmod[1] = 4.649869379728302800e+00;
                rmod[2] = 1.654457034571007900e+00;
                rmod[3] = 3.353583514350031900e+00;
                rmod[4] = 9.157719014108255100e-02;
                rmod[5] = 5.772657702308401400e+00;
                rmod[6] = 5.659358337346415800e+00;
                rmod[7] = 2.099930230068142400e-01;
                rmod[8] = 6.012690009399070900e+00;
                
                break;
            }
            case 13:
            {
                rmod[0] = 4.886448359475572500e+00;
                rmod[1] = 1.492515503572874100e+00;
                rmod[2] = 5.179094765441459600e+00;
                rmod[3] = 6.067981171564244200e+00;
                rmod[4] = 6.111033028633724700e+00;
                rmod[5] = 2.849105648924096900e+00;
                rmod[6] = 3.826726653470131600e+00;
                rmod[7] = 4.872776801893367700e+00;
                rmod[8] = 4.031375540680533800e+00;
                
                break;
            }
            case 14:
            {
                rmod[0] = 3.229201266315314500e+00;
                rmod[1] = 4.857939295249377000e+00;
                rmod[2] = 5.469058445908473200e+00;
                rmod[3] = 5.056046877018226200e-02;
                rmod[4] = 1.946128216239969300e+00;
                rmod[5] = 6.016801745283995500e+00;
                rmod[6] = 3.224007387389883600e+00;
                rmod[7] = 1.999840021860475000e+00;
                rmod[8] = 3.387893124464984600e+00;
                
                break;
            }
            case 15:
            {
                rmod[0] = 5.333278883951943600e+00;
                rmod[1] = 1.124036246977920200e+00;
                rmod[2] = 3.415741493812250500e-01;
                rmod[3] = 2.271613052442061200e+00;
                rmod[4] = 1.730395068203421900e+00;
                rmod[5] = 3.330089625864812500e+00;
                rmod[6] = 1.922145236540217400e+00;
                rmod[7] = 1.913068819850232400e+00;
                rmod[8] = 7.020911446719929600e-01;
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name << ": "
                    << "'random_seed' should be in between 0 - 15."
                    << std::endl);
            }
        }
        

        for (int k = 0; k < patch_dims_2; k++)
        {
            for (int j = 0; j < patch_dims_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = 0; i < patch_dims_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_source = (i + num_ghosts_source[0]) +
                        (j + num_ghosts_source[1])*ghostcell_dims_source[0] + 
                        (k + num_ghosts_source[2])*ghostcell_dims_source[0]*ghostcell_dims_source[1];

                    const int idx_cons_var = (i + num_ghosts_cons_var[0]) +
                        (j + num_ghosts_cons_var[1])*ghostcell_dims_cons_var[0] +
                        (k + num_ghosts_cons_var[2])*ghostcell_dims_cons_var[0]*ghostcell_dims_cons_var[1];

                    // Compute the coordinates.
                    double x[3];
                    x[0] = patch_xlo[0] + (double(i) + double(1)/double(2))*dx[0];
                    x[1] = patch_xlo[1] + (double(j) + double(1)/double(2))*dx[1];
                    x[2] = patch_xlo[2] + (double(k) + double(1)/double(2))*dx[2];

                    double eta = 0.0;        
                    for (int m = waven - 4; m <= waven + 4; m++)
                    {
                        eta += eta_0/3.0*cos(2.0*M_PI*m/width*x[1] + rmod[m-waven+4])*cos(2.0*M_PI*m/width*x[2] + rmod[m-waven+4]);
                    }
                    
                    const double X_2_H  = 0.5*(1.0 + erf((x[0] - eta)/delta)); // mass fraction of second species (Y_2)
                    const double R_H    = R_1*(1.0 - X_2_H) + X_2_H*R_2;
                    
                    if (!(x[0] >= d_special_source_box_lo[0] && x[0] <= d_special_source_box_hi[0]))
                    {
                        double rho_Y_0_ref;
                        double rho_Y_1_ref;
                        double rho_u_ref;
                        double rho_v_ref;
                        double rho_w_ref;
                        double E_ref;
                        double sponge_rate_tot;
                        double xi_b;
                        
                        const int N_int = 100000; // number of numerical quadrature points
                        const double dx_p = (x[0] - shift)/(N_int - 1.0);
                    
                        double integral = 0.0;
                        double p_H = 0.0;
                        double rho_H = 0.0;
                        
                        
                        // Compute the integral with linear interpolation.
                        
                        const int idx_integral_vector_lo = int(std::floor((x[0] - x_domain_lo)/dx_uniform));
                        const int idx_integral_vector_hi = idx_integral_vector_lo + 1;
                        
                        const double x_integral_vector_lo = double(idx_integral_vector_lo)*dx_uniform + 0.5*dx_uniform + x_domain_lo;
                        const double x_integral_vector_hi = x_integral_vector_lo + dx_uniform;
                        
                        const double weight_lo = ( x_integral_vector_hi - x[0])/dx_uniform;
                        const double weight_hi = (-x_integral_vector_lo + x[0])/dx_uniform;
                        
                        integral = weight_lo*integral_vector[idx_integral_vector_lo] + weight_hi*integral_vector[idx_integral_vector_hi];
                        
                        //for (int ii = 0; ii < N_int; ii++)
                        //{
                        //    const double x_p = shift + ii*dx_p;  //Bug fixed 3.22.2023
                        //    integral += 1.0/(0.5*(R_2 - R_1)*erf((x_p - shift)/(delta)) + 0.5*(R_1 + R_2))*dx_p;
                        //}
                        
                        
                        p_H = p_i*exp(g/T_0*integral);
                        rho_H = p_H/(R_H*T_0);
                        
                        double rho_ref, p_ref;
                        rho_ref = rho_H;
                        p_ref   = p_H;
                        
                        rho_Y_0_ref = rho_ref*(1.0 - X_2_H);
                        rho_Y_1_ref = rho_ref*X_2_H;
                        
                        const double u_ref = 0.0;
                        const double v_ref = 0.0;
                        const double w_ref = 0.0;
                        
                        rho_u_ref = rho_ref*u_ref;
                        rho_v_ref = rho_ref*v_ref;
                        rho_w_ref = rho_ref*w_ref;
                        E_ref     = p_ref/(gamma - double(1)) + double(1)/double(2)*rho_ref*(u_ref*u_ref + v_ref*v_ref + w_ref*w_ref);
                        
                        sponge_rate_tot = sponge_rate*300.0; //300.0 is roughly (gamma*p_i/rho_i)^.05 OLD:(pow((gamma*p_ref/rho_ref),0.5))*sponge_rate;
                        if(x[0] < 0.0)
                        {
                            const double erf_start_lo  = double(0.75)*(domain_xlo[0]-d_special_source_box_lo[0]) + d_special_source_box_lo[0]; //center of erf is 3/4 of the way into sponge
                            const double erf_offset_lo = double(-0.5) * erf((d_special_source_box_lo[0]-erf_start_lo)/(abs(erf_start_lo)*double(0.1))) + double(0.5); //value of erf at start of sponge
                            xi_b        = double(-0.5) * erf((x[0]-erf_start_lo)/(abs(erf_start_lo)*double(0.1))) + double(0.5) - erf_offset_lo; //subtract value of erf at start of sponge to start xi_b at zero
                            // xi_b        = 1.0;
                        }
                        else
                        {
                            const double erf_start_hi = double(0.75)*(domain_xhi[0]-d_special_source_box_hi[0]) + d_special_source_box_hi[0]; //center of erf is 3/4 of the way into sponge
                            const double erf_offset_hi = double(0.5) * erf((d_special_source_box_hi[0]-erf_start_hi)/(abs(erf_start_hi)*double(0.1))) + double(0.5); //value of erf at start of sponge
                            xi_b            = double(0.5) * erf((x[0]-erf_start_hi)/(abs(erf_start_hi)*double(0.1))) + double(0.5) - erf_offset_hi; //subtract value of erf at start of sponge to start xi_b at zero
                            // xi_b = 1.0;
                        }

                        const double rho_Y_0_p = rho_Y_0[idx_cons_var] - rho_Y_0_ref;
                        const double rho_Y_1_p = rho_Y_1[idx_cons_var] - rho_Y_1_ref;
                        const double rho_u_p   = rho_u[idx_cons_var]   - rho_u_ref;
                        const double rho_v_p   = rho_v[idx_cons_var]   - rho_v_ref;
                        const double rho_w_p   = rho_w[idx_cons_var]   - rho_w_ref;
                        const double E_p       = E[idx_cons_var]       - E_ref;
                    
                        S[0][idx_source] -= dt*sponge_rate_tot*xi_b*rho_Y_0_p;
                        S[1][idx_source] -= dt*sponge_rate_tot*xi_b*rho_Y_1_p;
                        S[2][idx_source] -= dt*sponge_rate_tot*xi_b*rho_u_p;
                        S[3][idx_source] -= dt*sponge_rate_tot*xi_b*rho_v_p;
                        S[4][idx_source] -= dt*sponge_rate_tot*xi_b*rho_w_p;
                        S[5][idx_source] -= dt*sponge_rate_tot*xi_b*E_p;
                    }
                }
            }
        }
    }            
}


void
FlowModelSpecialSourceTerms::putToRestart(const HAMERS_SHARED_PTR<tbox::Database>& restart_source_terms_db)
{
    putToRestartBase(restart_source_terms_db);
    
    double random_seed = int(0);
    if (d_source_terms_db->keyExists("random_seed"))
    {
        random_seed = d_source_terms_db->getInteger("random_seed");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'random_seed' found in data for source terms."
            << std::endl);
    }
    restart_source_terms_db->putInteger("random_seed", random_seed);

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
    
    double sponge_nu_fact = double(0);
    if (d_source_terms_db->keyExists("sponge_nu_fact"))
    {
        sponge_nu_fact = d_source_terms_db->getDouble("sponge_nu_fact");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'sponge_nu_fact' found in data for source terms."
            << std::endl);
    }
    restart_source_terms_db->putDouble("sponge_nu_fact", sponge_nu_fact);
    
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
