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
    // Follow Lundgren, T. S. "Linearly forced isotropic turbulence." Annual Research Briefs-2003 (2003): 461.
    
    if ((d_project_name != "3D HIT"))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Can only initialize data for 'project_name' = '3D HIT' !\n"
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
    
    TBOX_ASSERT(d_source_terms_db->keyExists("forcing_rate"));
    
    double eps_0 = double(0);
    if (d_source_terms_db->keyExists("forcing_rate"))
    {
        eps_0 = d_source_terms_db->getDouble("forcing_rate");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'forcing_rate' found in data for source terms."
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
    
    HAMERS_SHARED_PTR<pdat::CellData<double> > density      = conservative_variables[0];
    HAMERS_SHARED_PTR<pdat::CellData<double> > momentum     = conservative_variables[1];
    HAMERS_SHARED_PTR<pdat::CellData<double> > total_energy = conservative_variables[2];
    
    double* rho   = density->getPointer(0);
    double* rho_u = momentum->getPointer(0);
    double* rho_v = momentum->getPointer(1);
    double* rho_w = momentum->getPointer(2);
    double* E     = total_energy->getPointer(0);
    
    double gamma = double(5)/double(3);
    
    TBOX_ASSERT(d_source_terms_db != nullptr);
    
    const double TKE_avg = monitoring_statistics_map.at("KINETIC_ENERGY_AVG");
    const double Q_force = dt*eps_0/(TKE_avg*double(3));

    if (d_project_name == "3D HIT")
    {
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

                    S[1][idx_source] += Q_force*(rho_u[idx_cons_var]/rho[idx_cons_var]);
                    S[2][idx_source] += Q_force*(rho_v[idx_cons_var]/rho[idx_cons_var]);
                    S[3][idx_source] += Q_force*(rho_w[idx_cons_var]/rho[idx_cons_var]);
                }
            }
        }
    }
}

void
FlowModelSpecialSourceTerms::putToRestart(const HAMERS_SHARED_PTR<tbox::Database>& restart_source_terms_db)
{
    putToRestartBase(restart_source_terms_db);
    
    double eps_0 = double(0);
    if (d_source_terms_db->keyExists("forcing_rate"))
    {
        eps_0 = d_source_terms_db->getDouble("forcing_rate");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "No key 'forcing_rate' found in data for source terms."
            << std::endl);
    }
    
    restart_source_terms_db->putDouble("forcing_rate", eps_0);
}
