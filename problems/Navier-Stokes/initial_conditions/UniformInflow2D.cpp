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
    
    if (d_flow_model_type != FLOW_MODEL::SINGLE_SPECIES)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Flow model should be single-species!"
            << std::endl);
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
        const double* const dx = patch_geom->getDx();
        const double* const patch_xlo = patch_geom->getXLower();
        
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_cons_var = conservative_variables[0]->getGhostCellWidth();
        
        // Get the dimensions of the ghost cell boxes.
        const hier::Box ghost_box_cons_var = conservative_variables[0]->getGhostBox();
        const hier::IntVector ghostcell_dims_cons_var = ghost_box_cons_var.numberCells();
        
        /*
         * Initialize data for a 2D uniform flow.
         */
        
        HAMERS_SHARED_PTR<pdat::CellData<Real> > density      = conservative_variables[0];
        HAMERS_SHARED_PTR<pdat::CellData<Real> > momentum     = conservative_variables[1];
        HAMERS_SHARED_PTR<pdat::CellData<Real> > total_energy = conservative_variables[2];
        
        Real* rho   = density->getPointer(0);
        Real* rho_u = momentum->getPointer(0);
        Real* rho_v = momentum->getPointer(1);
        Real* E     = total_energy->getPointer(0);
        
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
        }
        
        for (int j = -num_ghosts_cons_var[1]; j < patch_dims[1] + num_ghosts_cons_var[1]; j++)
        {
            for (int i = -num_ghosts_cons_var[0]; i < patch_dims[0] + num_ghosts_cons_var[0]; i++)
            {
                // Compute index into linear data array.
                int idx_cell = (i + num_ghosts_cons_var[0]) +
                    (j + num_ghosts_cons_var[1])*ghostcell_dims_cons_var[0];
                
                Real x[2];
                x[0] = Real(patch_xlo[0]) + (Real(i) + half)*Real(dx[0]); // x-coordinates of the point.
                x[1] = Real(patch_xlo[1]) + (Real(j) + half)*Real(dx[1]); // y-coordinates of the point.
                
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
                
                if (x[0] < x_c)
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
                
                rho_u[idx_cell] = rho_inf*u_inf;
                rho_v[idx_cell] = rho_inf*v_inf;
                E[idx_cell]     = p_inf/(gamma - Real(1)) + half*rho_inf*(u_inf*u_inf + v_inf*v_inf);
            }
        }
    }
}
