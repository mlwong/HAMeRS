#include "flow/flow_models/FlowModelSpecialSourceTerms.hpp"

/*
 * Add the effects of the special source terms.
 */
void
FlowModelSpecialSourceTerms::computeSpecialSourceTermsOnPatch(
    HAMERS_SHARED_PTR<pdat::CellData<double> >& source,
    const hier::Patch& patch,
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& conservative_variables,
    const double time,
    const double dt,
    const int RK_step_number)
{
    // Follow Reckinger, Scott J., Daniel Livescu, and Oleg V. Vasilyev.
    // "Comprehensive numerical methodology for direct numerical simulations of compressible Rayleigh–Taylor instability."
    
    if ((d_project_name != "2D jet")) 
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Can only initialize data for 'project_name' = '2D jet' !\n"
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

    TBOX_ASSERT(d_source_terms_db->keyExists("U_jet"));
    
    double U_jet = double(0);
    if (d_source_terms_db->keyExists("U_jet"))
    {
        U_jet = d_source_terms_db->getDouble("U_jet");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'U_jet' found in data for source terms."
            << std::endl);
    }
    
    TBOX_ASSERT(d_source_terms_db->keyExists("theta_0"));
    
    double theta_0 = double(0);
    if (d_source_terms_db->keyExists("theta_0"))
    {
        theta_0 = d_source_terms_db->getDouble("theta_0");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'theta_0' found in data for source terms."
            << std::endl);
    }

    TBOX_ASSERT(d_source_terms_db->keyExists("D_jet"));
    
    double D_jet = double(0);
    if (d_source_terms_db->keyExists("D_jet"))
    {
        D_jet = d_source_terms_db->getDouble("D_jet");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'D_jet' found in data for source terms."
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
    
    const double W_0 = 0.0480; // molecular weight of heavier gas
    const double W_1 = 0.0160; // molecular weight of lighter gas
    
    const double p_ref = 100000.0; // interface pressure
    const double T_ref = 300.0;    // background temperature
    
    TBOX_ASSERT(d_source_terms_db != nullptr);
    // TBOX_ASSERT(d_source_terms_db->keyExists("has_gravity") || d_source_terms_db->keyExists("d_has_gravity"));
    
    // std::vector<double> gravity_vector;
    
    // if (d_source_terms_db->keyExists("has_gravity"))
    // {
    //     if (d_source_terms_db->keyExists("gravity"))
    //     {
    //         d_source_terms_db->getVector("gravity", gravity_vector);
    //     }
    //     else
    //     {
    //         TBOX_ERROR(d_object_name
    //             << ": "
    //             << "No key 'gravity' found in data for source terms."
    //             << std::endl);
    //     }
    // }
    // else if (d_source_terms_db->keyExists("d_has_gravity"))
    // {
    //     if (d_source_terms_db->keyExists("d_gravity"))
    //     {
    //         d_source_terms_db->getVector("d_gravity", gravity_vector);
    //     }
    //     else
    //     {
    //         TBOX_ERROR(d_object_name
    //             << ": "
    //             << "No key 'd_gravity' found in data for source terms."
    //             << std::endl);
    //     }
    // }
    
    // const double g = gravity_vector[0]; // gravity
    
    const double R_u = 8.31446261815324; // universal gas constant
    const double R_0 = R_u/W_0;          // gas constant of heavier gas
    const double R_1 = R_u/W_1;          // gas constant of lighter gas

    const double* const domain_xlo = d_grid_geometry->getXLower();
    
    if (d_project_name == "2D jet")
    {
        const double r_0 = D_jet/2.0;
        for (int j = 0; j < patch_dims[1]; j++)
        {
            for (int i = 0; i < patch_dims[0]; i++)
            {   
                // Compute the linear indices.
                const int idx_source = (i + num_ghosts_source[0]) +
                    (j + num_ghosts_source[1])*ghostcell_dims_source[0];
                
                const int idx_cons_var = (i + num_ghosts_cons_var[0]) +
                    (j + num_ghosts_cons_var[1])*ghostcell_dims_cons_var[0];
                
                // Compute the coordinates.

                double x[2];
                x[0] = patch_xlo[0] + (double(i) + double(1)/double(2))*dx[0];
                x[1] = patch_xlo[1] + (double(j) + double(1)/double(2))*dx[1];
                
                const double r = fabs(x[1]); //distance from jet center in y-direction

                // Check whether it is outside the special source box.
                if (x[0] <= d_special_source_box_lo[0])
                {                    
                    const double u_ref = U_jet*0.5*(1.0-tanh(r_0/(4.0*theta_0)*(r/r_0-r_0/r)));
                    const double v_ref = 0.0;
                    const double Z_ref = 0.5*(1.0-tanh(r_0/(4.0*theta_0)*(r/r_0-r_0/r)));
                    
                    const double rho0_ref = p_ref/(R_0*T_ref);
                    const double rho1_ref = p_ref/(R_1*T_ref);

                    const double rho_Y_0_ref = rho0_ref*Z_ref;
                    const double rho_Y_1_ref = rho1_ref*(1.0-Z_ref);
                    const double rho_ref     = rho_Y_0_ref + rho_Y_1_ref;

                    const double rho_u_ref = rho_ref * u_ref;
                    const double rho_v_ref = rho_ref * v_ref;
                    const double E_ref     = p_ref/(gamma - double(1)) + double(1)/double(2)*rho_ref*(u_ref*u_ref + v_ref*v_ref);

                    const double xi_b      = (1.0-(x[0]-domain_xlo[0])/(d_special_source_box_lo[0]-domain_xlo[0]))*sponge_rate; // mask value needs to be improved 

                    //sponge_rate_tot = (pow((p_ref/rho_ref),0.5))*sponge_rate;
                    //xi_b            = -(x[0])/(701.0-600.0); // mask value needs to be improved 
                    
                    const double rho_Y_0_p = rho_Y_0[idx_cons_var] - rho_Y_0_ref;
                    const double rho_Y_1_p = rho_Y_1[idx_cons_var] - rho_Y_1_ref;
                    const double rho_u_p   = rho_u[idx_cons_var]   - rho_u_ref;
                    const double rho_v_p   = rho_v[idx_cons_var]   - rho_v_ref;
                    const double E_p       = E[idx_cons_var]       - E_ref;
                    
                    S[0][idx_source] -= dt*xi_b*rho_Y_0_p;
                    S[1][idx_source] -= dt*xi_b*rho_Y_1_p;
                    S[2][idx_source] -= dt*xi_b*rho_u_p;
                    S[3][idx_source] -= dt*xi_b*rho_v_p;
                    S[4][idx_source] -= dt*xi_b*E_p;
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

    double U_jet = double(0);
    if (d_source_terms_db->keyExists("U_jet"))
    {
        U_jet = d_source_terms_db->getDouble("U_jet");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'U_jet' found in data for source terms."
            << std::endl);
    }
    
    restart_source_terms_db->putDouble("U_jet", U_jet);

    double theta_0 = double(0);
    if (d_source_terms_db->keyExists("theta_0"))
    {
        theta_0 = d_source_terms_db->getDouble("theta_0");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'theta_0' found in data for source terms."
            << std::endl);
    }
    
    restart_source_terms_db->putDouble("theta_0", theta_0);

    double D_jet = double(0);
    if (d_source_terms_db->keyExists("D_jet"))
    {
        D_jet = d_source_terms_db->getDouble("D_jet");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'D_jet' found in data for source terms."
            << std::endl);
    }
    
    restart_source_terms_db->putDouble("D_jet", D_jet);
}
