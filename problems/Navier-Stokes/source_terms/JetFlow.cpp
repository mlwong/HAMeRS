#include "flow/flow_models/FlowModelSpecialSourceTerms.hpp"

/*
 * Add the effects of the special source terms.
 */
void
FlowModelSpecialSourceTerms::computeSpecialSourceTermsOnPatch(
    HAMERS_SHARED_PTR<pdat::CellData<Real> >& source,
    const hier::Patch& patch,
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& conservative_variables,
    const std::unordered_map<std::string, Real>& monitoring_statistics_map,
    const double time,
    const double dt,
    const int RK_step_number)
{
    if ((d_project_name != "2D jet") && (d_project_name != "3D jet") ) 
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Can only initialize data for 'project_name' = '2D jet' or '3D jet' !\n"
            << "'project_name' = '"
            << d_project_name
            << "' is given."
            << std::endl);
    }
    
    if (d_dim != tbox::Dimension(2) && d_dim != tbox::Dimension(3) )
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Dimension of problem should be 2 oe 3!"
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
    
    Real sponge_rate = Real(0);
    if (d_source_terms_db->keyExists("sponge_rate"))
    {
        sponge_rate = d_source_terms_db->getReal("sponge_rate");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'sponge_rate' found in data for source terms."
            << std::endl);
    }
    
    TBOX_ASSERT(d_source_terms_db->keyExists("inflow_sponge_rate"));

    double inflow_sponge_rate = double(0);
    if (d_source_terms_db->keyExists("inflow_sponge_rate"))
    {
        inflow_sponge_rate = d_source_terms_db->getDouble("inflow_sponge_rate");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'inflow_sponge_rate' found in data for source terms."
            << std::endl);
    }

    TBOX_ASSERT(d_source_terms_db->keyExists("inflow_rate_power"));

    double inflow_rate_power = double(0);
    if (d_source_terms_db->keyExists("inflow_rate_power"))
    {
        inflow_rate_power = d_source_terms_db->getDouble("inflow_rate_power");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'inflow_rate_power' found in data for source terms."
            << std::endl);
    }

    
    TBOX_ASSERT(d_source_terms_db->keyExists("U_jet"));
    
    Real U_jet = Real(0);
    if (d_source_terms_db->keyExists("U_jet"))
    {
        U_jet = d_source_terms_db->getReal("U_jet");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'U_jet' found in data for source terms."
            << std::endl);
    }
    
    TBOX_ASSERT(d_source_terms_db->keyExists("theta_0"));
    
    Real theta_0 = Real(0);
    if (d_source_terms_db->keyExists("theta_0"))
    {
        theta_0 = d_source_terms_db->getReal("theta_0");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'theta_0' found in data for source terms."
            << std::endl);
    }

    TBOX_ASSERT(d_source_terms_db->keyExists("D_jet"));
    
    Real D_jet = Real(0);
    if (d_source_terms_db->keyExists("D_jet"))
    {
        D_jet = d_source_terms_db->getReal("D_jet");
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
    
    std::vector<Real*> S;
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
    
    HAMERS_SHARED_PTR<pdat::CellData<Real> > density         = conservative_variables[0];
    HAMERS_SHARED_PTR<pdat::CellData<Real> > momentum        = conservative_variables[1];
    HAMERS_SHARED_PTR<pdat::CellData<Real> > total_energy    = conservative_variables[2];
    
    Real* rho     = density->getPointer(0);
    Real* rho_u   = momentum->getPointer(0);
    Real* rho_v   = momentum->getPointer(1);
    Real* rho_w   = momentum->getPointer(2);
    Real* E       = total_energy->getPointer(0);
    
    const Real gamma = Real(7)/Real(5); // assume both gases have the same ratio of specific heat ratios
    
    const Real W_0 = Real(1.0000); // molecular weight of heavier gas
    //const Real W_1 = 0.0290; // molecular weight of lighter gas
    
    const Real p_ref = Real(100000.0); // interface pressure
    const Real T_ref = Real(300.0);    // background temperature
    
    TBOX_ASSERT(d_source_terms_db != nullptr);
    
    const Real R_u = Real(8.31446261815324); // universal gas constant
    const Real R_0 = R_u/W_0;                // gas constant
    
    const double* const domain_xlo = d_grid_geometry->getXLower();
    const double* const domain_xhi = d_grid_geometry->getXUpper();
    
    if (d_project_name == "2D jet")
    {
        const Real r_0 = D_jet/Real(2);
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
                Real x[2];
                x[0] = Real(patch_xlo[0]) + (Real(i) + Real(1)/Real(2))*Real(dx[0]);
                x[1] = Real(patch_xlo[1]) + (Real(j) + Real(1)/Real(2))*Real(dx[1]);
                
                const Real r = std::abs(x[1]); //distance from jet center in y-direction
                
                // Check whether it is outside the special source box.
                if (x[0] <= d_special_source_box_lo[0])
                {
                    const Real u_ref = U_jet*Real(0.5)*(Real(1)-std::tanh(r_0/(Real(4)*theta_0)*(r/r_0-r_0/r)));
                    const Real v_ref = Real(0);
                    
<<<<<<< HEAD
                    const double rho_ref = p_ref/(R_0*T_ref);

                    const double rho_u_ref = rho_ref * u_ref;
                    const double rho_v_ref = rho_ref * v_ref;
                    const double E_ref     = p_ref/(gamma - double(1)) + double(1)/double(2)*rho_ref*(u_ref*u_ref + v_ref*v_ref);

                    // 
                    // const double xi_b      = pow(fabs((x[0]-domain_xlo[0])/(d_special_source_box_lo[0]-domain_xlo[0])),inflow_rate_power)*inflow_sponge_rate; // mask value needs to be improved 
                    const double xi_b      = inflow_sponge_rate;                    
                    //

                    const double rho_p     = rho[idx_cons_var]     - rho_ref;
                    const double rho_u_p   = rho_u[idx_cons_var]   - rho_u_ref;
                    const double rho_v_p   = rho_v[idx_cons_var]   - rho_v_ref;
                    const double E_p       = E[idx_cons_var]       - E_ref;
                    

                    S[0][idx_source] -= dt*xi_b*rho_p;
                    S[1][idx_source] -= dt*xi_b*rho_u_p;
                    S[2][idx_source] -= dt*xi_b*rho_v_p;
                    S[3][idx_source] -= dt*xi_b*E_p;
=======
                    const Real rho_ref = p_ref/(R_0*T_ref);
                    
                    const Real rho_u_ref = rho_ref * u_ref;
                    const Real rho_v_ref = rho_ref * v_ref;
                    const Real E_ref     = p_ref/(gamma - Real(1)) + Real(1)/Real(2)*rho_ref*(u_ref*u_ref + v_ref*v_ref);
                    
                    const Real xi_b      = (Real(1)-(x[0]-Real(domain_xlo[0]))/(d_special_source_box_lo[0]-Real(domain_xlo[0])))*sponge_rate; // mask value needs to be improved 
                    
                    const Real rho_p     = rho[idx_cons_var]   - rho_ref;
                    const Real rho_u_p   = rho_u[idx_cons_var] - rho_u_ref;
                    const Real rho_v_p   = rho_v[idx_cons_var] - rho_v_ref;
                    const Real E_p       = E[idx_cons_var]     - E_ref;
                    
                    S[0][idx_source] -= Real(dt*double(xi_b*rho_p));
                    S[1][idx_source] -= Real(dt*double(xi_b*rho_u_p));
                    S[2][idx_source] -= Real(dt*double(xi_b*rho_v_p));
                    S[3][idx_source] -= Real(dt*double(xi_b*E_p));
>>>>>>> development
                }
                if (x[0] >= d_special_source_box_hi[0])
                {
                    const Real u_ref = Real(0);
                    const Real v_ref = Real(0);
                    
                    const Real rho_ref   = p_ref/(R_0*T_ref);
                    
                    const Real rho_u_ref = rho_ref * u_ref;
                    const Real rho_v_ref = rho_ref * v_ref;
                    const Real E_ref     = p_ref/(gamma - Real(1)) + Real(1)/Real(2)*rho_ref*(u_ref*u_ref + v_ref*v_ref);
                    
                    const Real xi_b      = (x[0]-d_special_source_box_hi[0])/(Real(domain_xhi[0])-d_special_source_box_hi[0])*sponge_rate/Real(1); // mask value needs to be improved 
                    
                    const Real rho_p     = rho[idx_cons_var] - rho_ref;
                    const Real rho_u_p   = rho_u[idx_cons_var] - rho_u_ref;
                    const Real rho_v_p   = rho_v[idx_cons_var] - rho_v_ref;
                    const Real E_p       = E[idx_cons_var]     - E_ref;
                    
                    S[0][idx_source] -= Real(dt*double(xi_b*rho_p));
                    S[1][idx_source] -= Real(dt*double(xi_b*rho_u_p));
                    S[2][idx_source] -= Real(dt*double(xi_b*rho_v_p));
                    S[3][idx_source] -= Real(dt*double(xi_b*E_p));
                }
                if (x[1] <= d_special_source_box_lo[1])
                {
                    const Real u_ref = Real(0);
                    const Real v_ref = Real(0);
                    
                    const Real rho_ref = p_ref/(R_0*T_ref);
                    
                    const Real rho_u_ref = rho_ref * u_ref;
                    const Real rho_v_ref = rho_ref * v_ref;
                    const Real E_ref     = p_ref/(gamma - Real(1)) + Real(1)/Real(2)*rho_ref*(u_ref*u_ref + v_ref*v_ref);
                    
                    const Real xi_b      = (Real(1)-(x[1]-Real(domain_xlo[1]))/(d_special_source_box_lo[1]-Real(domain_xlo[1])))*sponge_rate/Real(1); // mask value needs to be improved 
                    
                    const Real rho_p     = rho[idx_cons_var] - rho_ref;
                    const Real rho_u_p   = rho_u[idx_cons_var] - rho_u_ref;
                    const Real rho_v_p   = rho_v[idx_cons_var] - rho_v_ref;
                    const Real E_p       = E[idx_cons_var]     - E_ref;
                    
                    S[0][idx_source] -= Real(dt*double(xi_b*rho_p));
                    S[1][idx_source] -= Real(dt*double(xi_b*rho_u_p));
                    S[2][idx_source] -= Real(dt*double(xi_b*rho_v_p));
                    S[3][idx_source] -= Real(dt*double(xi_b*E_p));
                }
                if (x[1] >= d_special_source_box_hi[1])
                {
                    const Real u_ref = Real(0);
                    const Real v_ref = Real(0);
                    
                    const Real rho_ref = p_ref/(R_0*T_ref);
                    
                    const Real rho_u_ref = rho_ref * u_ref;
                    const Real rho_v_ref = rho_ref * v_ref;
                    const Real E_ref     = p_ref/(gamma - Real(1)) + Real(1)/Real(2)*rho_ref*(u_ref*u_ref + v_ref*v_ref);
                    
                    const Real xi_b      = (x[1]-d_special_source_box_hi[1])/(Real(domain_xhi[1])-d_special_source_box_hi[1])*sponge_rate/Real(1); // mask value needs to be improved
                    
                    const Real rho_p     = rho[idx_cons_var]   - rho_ref;
                    const Real rho_u_p   = rho_u[idx_cons_var] - rho_u_ref;
                    const Real rho_v_p   = rho_v[idx_cons_var] - rho_v_ref;
                    const Real E_p       = E[idx_cons_var]     - E_ref;
                    
                    S[0][idx_source] -= Real(dt*double(xi_b*rho_p));
                    S[1][idx_source] -= Real(dt*double(xi_b*rho_u_p));
                    S[2][idx_source] -= Real(dt*double(xi_b*rho_v_p));
                    S[3][idx_source] -= Real(dt*double(xi_b*E_p));
                }
            }
        }
    }
    else if (d_project_name == "3D jet")
    {
        const Real r_0 = D_jet/Real(2);
        
        for (int k = 0; k < patch_dims[2]; k++)
        {
            for (int j = 0; j < patch_dims[1]; j++)
            {
                for (int i = 0; i < patch_dims[0]; i++)
                {
                    // Compute the linear indices.
                    const int idx_source = (i + num_ghosts_source[0]) +
                        (j + num_ghosts_source[1])*ghostcell_dims_source[0] +
                        (k + num_ghosts_source[2])*ghostcell_dims_source[0]*ghostcell_dims_source[1];
                    
                    const int idx_cons_var = (i + num_ghosts_cons_var[0]) +
                        (j + num_ghosts_cons_var[1])*ghostcell_dims_cons_var[0] +
                        (k + num_ghosts_cons_var[2])*ghostcell_dims_cons_var[0]*ghostcell_dims_cons_var[1];
                    
                    // Compute the coordinates.
                    Real x[3];
                    x[0] = Real(patch_xlo[0]) + (Real(i) + Real(1)/Real(2))*Real(dx[0]);
                    x[1] = Real(patch_xlo[1]) + (Real(j) + Real(1)/Real(2))*Real(dx[1]);
                    x[2] = Real(patch_xlo[2]) + (Real(k) + Real(1)/Real(2))*Real(dx[2]);
                    
                    const Real r = std::sqrt(std::abs(x[1])*std::abs(x[1])+std::abs(x[2])*std::abs(x[2]));
                    
                    if (x[0] <= d_special_source_box_lo[0])
                    {
                        const Real u_ref = U_jet*Real(0.5)*(Real(1)-std::tanh(r_0/(Real(4)*theta_0)*(r/r_0-r_0/r)));
                        const Real v_ref = Real(0);
                        const Real w_ref = Real(0);
                        
                        const Real rho_ref = p_ref/(R_0*T_ref);
                        
                        const Real rho_u_ref = rho_ref * u_ref;
                        const Real rho_v_ref = rho_ref * v_ref;
                        const Real rho_w_ref = rho_ref * w_ref;
                        const Real E_ref     = p_ref/(gamma - Real(1)) + Real(1)/Real(2)*rho_ref*(u_ref*u_ref + v_ref*v_ref + w_ref*w_ref);
                        
                        const Real xi_b      = (Real(1)-(x[0]-Real(domain_xlo[0]))/(d_special_source_box_lo[0]-Real(domain_xlo[0])))*sponge_rate; // mask value needs to be improved
                        
                        const Real rho_p     = rho[idx_cons_var]   - rho_ref;
                        const Real rho_u_p   = rho_u[idx_cons_var] - rho_u_ref;
                        const Real rho_v_p   = rho_v[idx_cons_var] - rho_v_ref;
                        const Real rho_w_p   = rho_w[idx_cons_var] - rho_w_ref;
                        const Real E_p       = E[idx_cons_var]     - E_ref;
                        
                        S[0][idx_source] -= Real(dt*double(xi_b*rho_p));
                        S[1][idx_source] -= Real(dt*double(xi_b*rho_u_p));
                        S[2][idx_source] -= Real(dt*double(xi_b*rho_v_p));
                        S[3][idx_source] -= Real(dt*double(xi_b*rho_w_p));
                        S[4][idx_source] -= Real(dt*double(xi_b*E_p));
                    }
                    if (x[0] >= d_special_source_box_hi[0])
                    {
                        const Real u_ref = Real(0);
                        const Real v_ref = Real(0);
                        const Real w_ref = Real(0);
                        
                        const Real rho_ref   = p_ref/(R_0*T_ref);
                        
                        const Real rho_u_ref = rho_ref * u_ref;
                        const Real rho_v_ref = rho_ref * v_ref;
                        const Real rho_w_ref = rho_ref * w_ref;
                        const Real E_ref     = p_ref/(gamma - Real(1)) + Real(1)/Real(2)*rho_ref*(u_ref*u_ref + v_ref*v_ref + w_ref*w_ref);
                        
                        const Real xi_b      = (x[0]-d_special_source_box_hi[0])/(Real(domain_xhi[0])-d_special_source_box_hi[0])*sponge_rate/Real(1); // mask value needs to be improved
                        
                        const Real rho_p     = rho[idx_cons_var]   - rho_ref;
                        const Real rho_u_p   = rho_u[idx_cons_var] - rho_u_ref;
                        const Real rho_v_p   = rho_v[idx_cons_var] - rho_v_ref;
                        const Real rho_w_p   = rho_w[idx_cons_var] - rho_w_ref;
                        const Real E_p       = E[idx_cons_var]     - E_ref;
                        
                        S[0][idx_source] -= Real(dt*double(xi_b*rho_p));
                        S[1][idx_source] -= Real(dt*double(xi_b*rho_u_p));
                        S[2][idx_source] -= Real(dt*double(xi_b*rho_v_p));
                        S[3][idx_source] -= Real(dt*double(xi_b*rho_w_p));
                        S[4][idx_source] -= Real(dt*double(xi_b*E_p));
                    }
                    if (x[1] <= d_special_source_box_lo[1])
                    {
                        const Real u_ref = Real(0);
                        const Real v_ref = Real(0);
                        const Real w_ref = Real(0);
                        
                        const Real rho_ref = p_ref/(R_0*T_ref);
                        
                        const Real rho_u_ref = rho_ref * u_ref;
                        const Real rho_v_ref = rho_ref * v_ref;
                        const Real rho_w_ref = rho_ref * w_ref;
                        const Real E_ref     = p_ref/(gamma - Real(1)) + Real(1)/Real(2)*rho_ref*(u_ref*u_ref + v_ref*v_ref + w_ref*w_ref);
                        
                        const Real xi_b      = (Real(1)-(x[1]-Real(domain_xlo[1]))/(d_special_source_box_lo[1]-Real(domain_xlo[1])))*sponge_rate/Real(1); // mask value needs to be improved
                        
                        const Real rho_p     = rho[idx_cons_var]   - rho_ref;
                        const Real rho_u_p   = rho_u[idx_cons_var] - rho_u_ref;
                        const Real rho_v_p   = rho_v[idx_cons_var] - rho_v_ref;
                        const Real rho_w_p   = rho_w[idx_cons_var] - rho_w_ref;
                        const Real E_p       = E[idx_cons_var]     - E_ref;
                        
                        S[0][idx_source] -= Real(dt*double(xi_b*rho_p));
                        S[1][idx_source] -= Real(dt*double(xi_b*rho_u_p));
                        S[2][idx_source] -= Real(dt*double(xi_b*rho_v_p));
                        S[3][idx_source] -= Real(dt*double(xi_b*rho_w_p));
                        S[4][idx_source] -= Real(dt*double(xi_b*E_p));
                    }
                    if (x[1] >= d_special_source_box_hi[1])
                    {
                        const Real u_ref = Real(0);
                        const Real v_ref = Real(0);
                        const Real w_ref = Real(0);
                        
                        const Real rho_ref = p_ref/(R_0*T_ref);
                        
                        const Real rho_u_ref = rho_ref * u_ref;
                        const Real rho_v_ref = rho_ref * v_ref;
                        const Real rho_w_ref = rho_ref * w_ref;
                        const Real E_ref     = p_ref/(gamma - Real(1)) + Real(1)/Real(2)*rho_ref*(u_ref*u_ref + v_ref*v_ref + w_ref*w_ref);
                        
                        const Real xi_b      = (x[1]-d_special_source_box_hi[1])/(Real(domain_xhi[1])-d_special_source_box_hi[1])*sponge_rate/Real(1); // mask value needs to be improved
                        
                        const Real rho_p     = rho[idx_cons_var]   - rho_ref;
                        const Real rho_u_p   = rho_u[idx_cons_var] - rho_u_ref;
                        const Real rho_v_p   = rho_v[idx_cons_var] - rho_v_ref;
                        const Real rho_w_p   = rho_w[idx_cons_var] - rho_w_ref;
                        const Real E_p       = E[idx_cons_var]     - E_ref;
                        
                        S[0][idx_source] -= Real(dt*double(xi_b*rho_p));
                        S[1][idx_source] -= Real(dt*double(xi_b*rho_u_p));
                        S[2][idx_source] -= Real(dt*double(xi_b*rho_v_p));
                        S[3][idx_source] -= Real(dt*double(xi_b*rho_w_p));
                        S[4][idx_source] -= Real(dt*double(xi_b*E_p));
                    } 
                    if (x[2] <= d_special_source_box_lo[2])
                    {
                        const Real u_ref = Real(0);
                        const Real v_ref = Real(0);
                        const Real w_ref = Real(0);
                        
                        const Real rho_ref = p_ref/(R_0*T_ref);
                        
                        const Real rho_u_ref = rho_ref * u_ref;
                        const Real rho_v_ref = rho_ref * v_ref;
                        const Real rho_w_ref = rho_ref * w_ref;
                        const Real E_ref     = p_ref/(gamma - Real(1)) + Real(1)/Real(2)*rho_ref*(u_ref*u_ref + v_ref*v_ref +w_ref*w_ref);
                        
                        const Real xi_b      = (Real(1)-(x[2]-Real(domain_xlo[2]))/(d_special_source_box_lo[2]-Real(domain_xlo[2])))*sponge_rate/Real(1); // mask value needs to be improved
                        
                        const Real rho_p     = rho[idx_cons_var]   - rho_ref;
                        const Real rho_u_p   = rho_u[idx_cons_var] - rho_u_ref;
                        const Real rho_v_p   = rho_v[idx_cons_var] - rho_v_ref;
                        const Real rho_w_p   = rho_w[idx_cons_var] - rho_w_ref;
                        const Real E_p       = E[idx_cons_var]     - E_ref;
                        
                        S[0][idx_source] -= Real(dt*double(xi_b*rho_p));
                        S[1][idx_source] -= Real(dt*double(xi_b*rho_u_p));
                        S[2][idx_source] -= Real(dt*double(xi_b*rho_v_p));
                        S[3][idx_source] -= Real(dt*double(xi_b*rho_w_p));
                        S[4][idx_source] -= Real(dt*double(xi_b*E_p));
                    }
                    if (x[2] >= d_special_source_box_hi[2])
                    {
                        const Real u_ref = Real(0);
                        const Real v_ref = Real(0);
                        const Real w_ref = Real(0);
                        
                        const Real rho_ref = p_ref/(R_0*T_ref);
                        
                        const Real rho_u_ref = rho_ref * u_ref;
                        const Real rho_v_ref = rho_ref * v_ref;
                        const Real rho_w_ref = rho_ref * w_ref;
                        const Real E_ref     = p_ref/(gamma - Real(1)) + Real(1)/Real(2)*rho_ref*(u_ref*u_ref + v_ref*v_ref + w_ref*w_ref);
                        
                        const Real xi_b      = (x[2]-d_special_source_box_hi[2])/(Real(domain_xhi[2])-d_special_source_box_hi[2])*sponge_rate/Real(1); // mask value needs to be improved
                        
                        const Real rho_p     = rho[idx_cons_var]   - rho_ref;
                        const Real rho_u_p   = rho_u[idx_cons_var] - rho_u_ref;
                        const Real rho_v_p   = rho_v[idx_cons_var] - rho_v_ref;
                        const Real rho_w_p   = rho_w[idx_cons_var] - rho_w_ref;
                        const Real E_p       = E[idx_cons_var]     - E_ref;
                        
                        S[0][idx_source] -= Real(dt*double(xi_b*rho_p));
                        S[1][idx_source] -= Real(dt*double(xi_b*rho_u_p));
                        S[2][idx_source] -= Real(dt*double(xi_b*rho_v_p));
                        S[3][idx_source] -= Real(dt*double(xi_b*rho_w_p));
                        S[4][idx_source] -= Real(dt*double(xi_b*E_p));
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
    
    Real sponge_rate = Real(0);
    if (d_source_terms_db->keyExists("sponge_rate"))
    {
        sponge_rate = d_source_terms_db->getReal("sponge_rate");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'sponge_rate' found in data for source terms."
            << std::endl);
    }
    
<<<<<<< HEAD
    restart_source_terms_db->putDouble("sponge_rate", sponge_rate);
    
    double inflow_sponge_rate = double(0);
    if (d_source_terms_db->keyExists("inflow_sponge_rate"))
    {
        inflow_sponge_rate = d_source_terms_db->getDouble("inflow_sponge_rate");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'inflow_sponge_rate' found in data for source terms."
            << std::endl);
    }

    restart_source_terms_db->putDouble("inflow_sponge_rate", inflow_sponge_rate);

    double inflow_rate_power = double(0);
    if (d_source_terms_db->keyExists("inflow_rate_power"))
    {
        inflow_rate_power = d_source_terms_db->getDouble("inflow_rate_power");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'inflow_rate_power' found in data for source terms."
            << std::endl);
    }

    restart_source_terms_db->putDouble("inflow_rate_power", inflow_rate_power);


    double U_jet = double(0);
=======
    restart_source_terms_db->putReal("sponge_rate", sponge_rate);
    
    Real U_jet = Real(0);
>>>>>>> development
    if (d_source_terms_db->keyExists("U_jet"))
    {
        U_jet = d_source_terms_db->getReal("U_jet");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'U_jet' found in data for source terms."
            << std::endl);
    }
    
    restart_source_terms_db->putReal("U_jet", U_jet);
    
    Real theta_0 = Real(0);
    if (d_source_terms_db->keyExists("theta_0"))
    {
        theta_0 = d_source_terms_db->getReal("theta_0");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'theta_0' found in data for source terms."
            << std::endl);
    }
    
    restart_source_terms_db->putReal("theta_0", theta_0);
    
    Real D_jet = Real(0);
    if (d_source_terms_db->keyExists("D_jet"))
    {
        D_jet = d_source_terms_db->getReal("D_jet");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'D_jet' found in data for source terms."
            << std::endl);
    }
    
    restart_source_terms_db->putReal("D_jet", D_jet);
}
