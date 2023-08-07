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
    if ((d_project_name != "2D uniform flow") ) 
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
    
    TBOX_ASSERT(d_source_terms_db->keyExists("rho_inf"));
    
    Real rho_inf = Real(0);
    if (d_source_terms_db->keyExists("rho_inf"))
    {
        rho_inf = d_source_terms_db->getReal("rho_inf");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'rho_inf' found in data for source terms."
            << std::endl);
    }
    
    TBOX_ASSERT(d_source_terms_db->keyExists("u_inf"));
    
    Real u_inf = Real(0);
    if (d_source_terms_db->keyExists("u_inf"))
    {
        u_inf = d_source_terms_db->getReal("u_inf");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'u_inf' found in data for source terms."
            << std::endl);
    }
    
    TBOX_ASSERT(d_source_terms_db->keyExists("v_inf"));
    
    Real v_inf = Real(0);
    if (d_source_terms_db->keyExists("v_inf"))
    {
        v_inf = d_source_terms_db->getReal("v_inf");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'v_inf' found in data for source terms."
            << std::endl);
    }
    
    TBOX_ASSERT(d_source_terms_db->keyExists("p_inf"));
    
    Real p_inf = Real(0);
    if (d_source_terms_db->keyExists("p_inf"))
    {
        p_inf = d_source_terms_db->getReal("p_inf");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'p_inf' found in data for source terms."
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
    
    HAMERS_SHARED_PTR<pdat::CellData<Real> > density      = conservative_variables[0];
    HAMERS_SHARED_PTR<pdat::CellData<Real> > momentum     = conservative_variables[1];
    HAMERS_SHARED_PTR<pdat::CellData<Real> > total_energy = conservative_variables[2];
    
    Real* rho   = density->getPointer(0);
    Real* rho_u = momentum->getPointer(0);
    Real* rho_v = momentum->getPointer(1);
    Real* E     = total_energy->getPointer(0);
    
    const Real gamma = Real(7)/Real(5); // assume both gases have the same ratio of specific heat ratios
    
    TBOX_ASSERT(d_source_terms_db != nullptr);
    
    const double* const domain_xlo = d_grid_geometry->getXLower();
    const double* const domain_xhi = d_grid_geometry->getXUpper();
    
    if (d_project_name == "2D uniform flow")
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
                
                // Compute the coordinates.
                
                Real x[2];
                x[0] = patch_xlo[0] + (Real(i) + Real(1)/Real(2))*Real(dx[0]);
                x[1] = patch_xlo[1] + (Real(j) + Real(1)/Real(2))*Real(dx[1]);
                
                // Check whether it is outside the special source box.
                if (x[0] <= d_special_source_box_lo[0])
                {
                    const Real u_ref = u_inf;
                    const Real v_ref = v_inf;
                    
                    const Real rho_ref = rho_inf;
                    const Real p_ref   = p_inf;
                    
                    const Real rho_u_ref = rho_ref * u_ref;
                    const Real rho_v_ref = rho_ref * v_ref;
                    const Real E_ref     = p_ref/(gamma - Real(1)) + Real(1)/Real(2)*rho_ref*(u_ref*u_ref + v_ref*v_ref);
                    
                    Real xi_b = std::pow((Real(1) - (x[0] - domain_xlo[0])/(d_special_source_box_lo[0] - domain_xlo[0])), Real(3));
                    xi_b      *= sponge_rate; // mask value needs to be improved 
                    
                    const Real rho_p   = rho[idx_cons_var]   - rho_ref;
                    const Real rho_u_p = rho_u[idx_cons_var] - rho_u_ref;
                    const Real rho_v_p = rho_v[idx_cons_var] - rho_v_ref;
                    const Real E_p     = E[idx_cons_var]     - E_ref;
                    
                    S[0][idx_source] -= dt*xi_b*rho_p;
                    S[1][idx_source] -= dt*xi_b*rho_u_p;
                    S[2][idx_source] -= dt*xi_b*rho_v_p;
                    S[3][idx_source] -= dt*xi_b*E_p;
                }
                if (x[0] >= d_special_source_box_hi[0])
                {                    
                    const Real u_ref = u_inf;
                    const Real v_ref = v_inf;
                    
                    const Real rho_ref = rho_inf;
                    const Real p_ref   = p_inf;

                    const Real rho_u_ref = rho_ref * u_ref;
                    const Real rho_v_ref = rho_ref * v_ref;
                    const Real E_ref     = p_ref/(gamma - Real(1)) + Real(1)/Real(2)*rho_ref*(u_ref*u_ref + v_ref*v_ref);

                    Real xi_b      = std::pow((x[0]-d_special_source_box_hi[0])/(domain_xhi[0]-d_special_source_box_hi[0]), Real(3)); // mask value needs to be improved 
                    xi_b           = xi_b * sponge_rate;

                    const Real rho_p   = rho[idx_cons_var]   - rho_ref;
                    const Real rho_u_p = rho_u[idx_cons_var] - rho_u_ref;
                    const Real rho_v_p = rho_v[idx_cons_var] - rho_v_ref;
                    const Real E_p     = E[idx_cons_var]     - E_ref;
                    
                    S[0][idx_source] -= dt*xi_b*rho_p;
                    S[1][idx_source] -= dt*xi_b*rho_u_p;
                    S[2][idx_source] -= dt*xi_b*rho_v_p;
                    S[3][idx_source] -= dt*xi_b*E_p;
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
    
    restart_source_terms_db->putReal("sponge_rate", sponge_rate);

    Real rho_inf = Real(0);
    if (d_source_terms_db->keyExists("rho_inf"))
    {
        rho_inf = d_source_terms_db->getReal("rho_inf");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'rho_inf' found in data for source terms."
            << std::endl);
    }
    
    restart_source_terms_db->putReal("rho_inf", rho_inf);
    
    Real u_inf = Real(0);
    if (d_source_terms_db->keyExists("u_inf"))
    {
        u_inf = d_source_terms_db->getReal("u_inf");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'u_inf' found in data for source terms."
            << std::endl);
    }
    
    restart_source_terms_db->putReal("u_inf", u_inf);

    Real v_inf = Real(0);
    if (d_source_terms_db->keyExists("v_inf"))
    {
        v_inf = d_source_terms_db->getReal("v_inf");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'v_inf' found in data for source terms."
            << std::endl);
    }
    
    restart_source_terms_db->putReal("v_inf", v_inf);

    Real p_inf = Real(0);
    if (d_source_terms_db->keyExists("p_inf"))
    {
        p_inf = d_source_terms_db->getReal("p_inf");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'p_inf' found in data for source terms."
            << std::endl);
    }
    
    restart_source_terms_db->putReal("p_inf", p_inf);
}
