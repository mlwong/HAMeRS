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
    
    if (d_project_name != "3D uniform flow")
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Can only initialize data for 'project_name' = '3D uniform flow'!\n"
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
    
    if (d_flow_model_type != FLOW_MODEL::SINGLE_SPECIES && d_flow_model_type != FLOW_MODEL::FIVE_EQN_ALLAIRE)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Flow model should be single-species or five-equation model by Allaire!"
            << std::endl);
    }
    
    if (d_flow_model_type == FLOW_MODEL::FIVE_EQN_ALLAIRE)
    {
        if (d_flow_model->getNumberOfSpecies() != 2)
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Number of species should be 2!"
                << std::endl);
        }
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
        const double* const dx = patch_geom->getDx();            //AFK
        const double* const patch_xlo = patch_geom->getXLower(); //AFK
        
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_cons_var = conservative_variables[0]->getGhostCellWidth();
        
        // Get the dimensions of the ghost cell boxes.
        const hier::Box ghost_box_cons_var = conservative_variables[0]->getGhostBox();
        const hier::IntVector ghostcell_dims_cons_var = ghost_box_cons_var.numberCells();
        
        /*
         * Initialize data for a 3D flow over a sphere problem.
         */
        
        if (d_flow_model_type == FLOW_MODEL::SINGLE_SPECIES)
        {
            HAMERS_SHARED_PTR<pdat::CellData<Real> > density      = conservative_variables[0];
            HAMERS_SHARED_PTR<pdat::CellData<Real> > momentum     = conservative_variables[1];
            HAMERS_SHARED_PTR<pdat::CellData<Real> > total_energy = conservative_variables[2];
            
            Real* rho   = density->getPointer(0);
            Real* rho_u = momentum->getPointer(0);
            Real* rho_v = momentum->getPointer(1);
            Real* rho_w = momentum->getPointer(2);
            Real* E     = total_energy->getPointer(0);
            
            // Initial conditions.
            Real rho_inf = Real(1);
            Real u_inf   = Real(1);
            Real v_inf   = Real(1);
            Real w_inf   = Real(1); 
            Real p_inf   = Real(1);
            Real gamma   = Real(1);
            Real x_c     = Real(1);
            Real y_c     = Real(1); 
            Real z_c     = Real(1); 
            Real D       = Real(1);
            Real r       = Real(1);
            Real theta   = Real(1);
            Real u_ic    = Real(1);
            Real v_ic    = Real(1);
            Real w_ic    = Real(1); 
            Real p_ic    = Real(1);
            Real u       = Real(1);
            Real v       = Real(1);
            Real w       = Real(1); 
            Real p       = Real(1);
            Real spongeR = Real(1);
            Real spongeL = Real(1);
            Real spongeB = Real(1);
            Real spongeT = Real(1);
            Real spongeF = Real(1); 
            Real spongeK = Real(1); 

            Real half = Real(1)/Real(2);
            
            if (d_initial_conditions_db != nullptr) 
            {
                TBOX_ASSERT(d_initial_conditions_db->keyExists("rho_inf"));
                TBOX_ASSERT(d_initial_conditions_db->keyExists("u_inf"));
                TBOX_ASSERT(d_initial_conditions_db->keyExists("v_inf"));
                TBOX_ASSERT(d_initial_conditions_db->keyExists("w_inf")); 
                TBOX_ASSERT(d_initial_conditions_db->keyExists("p_inf"));
                
                rho_inf = d_initial_conditions_db->getReal("rho_inf");
                u_inf   = d_initial_conditions_db->getReal("u_inf");
                v_inf   = d_initial_conditions_db->getReal("v_inf");
                w_inf   = d_initial_conditions_db->getReal("w_inf"); 
                p_inf   = d_initial_conditions_db->getReal("p_inf");
                gamma   = d_initial_conditions_db->getReal("gamma");
                x_c     = d_initial_conditions_db->getReal("x_c");
                y_c     = d_initial_conditions_db->getReal("y_c");
                z_c     = d_initial_conditions_db->getReal("z_c"); 
                D       = d_initial_conditions_db->getReal("D");
                spongeL = d_initial_conditions_db->getReal("spongeL");
                spongeR = d_initial_conditions_db->getReal("spongeR");
                spongeB = d_initial_conditions_db->getReal("spongeB"); 
                spongeT = d_initial_conditions_db->getReal("spongeT");
                spongeF = d_initial_conditions_db->getReal("spongeF");
                spongeK = d_initial_conditions_db->getReal("spongeK");
            }
            

            for (int k = -num_ghosts_cons_var[2]; k < patch_dims[2] + num_ghosts_cons_var[2]; k++)
            {
                for (int j = -num_ghosts_cons_var[1]; j < patch_dims[1] + num_ghosts_cons_var[1]; j++)
                {
                    for (int i = -num_ghosts_cons_var[0]; i < patch_dims[0] + num_ghosts_cons_var[0]; i++)
                    
                    {
                        // Compute index into linear data array.
                        int idx_cell = (i + num_ghosts_cons_var[0]) +
                            (j + num_ghosts_cons_var[1]) * ghostcell_dims_cons_var[0] + 
                            (k + num_ghosts_cons_var[2]) * ghostcell_dims_cons_var[1] * ghostcell_dims_cons_var[0];
                        
                        Real x[3];
                        x[0] = patch_xlo[0] + (Real(i) + half)*Real(dx[0]); // x coordinates of the point.
                        x[1] = patch_xlo[1] + (Real(j) + half)*Real(dx[1]); 
                        x[2] = patch_xlo[2] + (Real(k) + half)*Real(dx[2]); 

                        rho[idx_cell] = rho_inf;  
                        
                        // radial distance from the center
                        r       = std::pow((x[0] - x_c), 2) + std::pow((x[1] - y_c), 2) + std::pow(x[2] - z_c, 2);
                        r       = std::pow(r, half); 
                        
                        if (r < D/Real(4)) 
                        {
                
                            u_ic = 0;
                            v_ic = 0;
                            w_ic = 0;
                
                        }

                        else 
                        {
                
                            u_ic =  u_inf * (Real(1) - (Real(3) * std::pow(D, 3) * std::pow(x[0],2)) / (Real(16) * std::pow(r, 5)) + (std::pow(D, 3) / (Real(16)*std::pow(r, 3))));
                            v_ic = Real(0); //-u_inf*((Real(3) * std::pow(D,3) * x[0] * x[1])/(Real(16) * std::pow(r,5)));
                            w_ic = Real(0); //-u_inf*((Real(3) * std::pow(D,3) * x[0] * x[2])/(Real(16) * std::pow(r,5)));
                
                        }

                        p_ic = p_inf + half*rho_inf*(u_inf*u_inf - (u_ic*u_ic + v_ic*v_ic + w_ic*w_ic)); 
                        
                        if(x[0] < x_c)
                        {
                            u = half * (u_inf + u_ic) + half * (u_ic - u_inf)*erf((x[0]-spongeL)/(D/Real(5)));
                            v = half * (v_inf + v_ic) + half * (v_ic - v_inf)*erf((x[0]-spongeL)/(D/Real(5)));
                            w = half * (w_inf + w_ic) + half * (w_ic - w_inf)*erf((x[0]-spongeL)/(D/Real(5)));
                            p = half * (p_inf + p_ic) + half * (p_ic - p_inf)*erf((x[0]-spongeL)/(D/Real(5)));
                        }
                        else
                        {
                            u = half * (u_inf + u_ic) - half * (u_ic - u_inf)*erf((x[0]-spongeR)/(D/Real(5)));
                            v = half * (v_inf + v_ic) - half * (v_ic - v_inf)*erf((x[0]-spongeR)/(D/Real(5)));
                            w = half * (w_inf + w_ic) - half * (w_ic - w_inf)*erf((x[0]-spongeR)/(D/Real(5)));
                            p = half * (p_inf + p_ic) - half * (p_ic - p_inf)*erf((x[0]-spongeR)/(D/Real(5)));
                        }

                        if (x[1] < y_c) 
                        {
                            u = half * (u_inf + u_ic) + half * (u_ic - u_inf)*erf((x[1]-spongeB)/(D/Real(5)));
                            v = half * (v_inf + v_ic) + half * (v_ic - v_inf)*erf((x[1]-spongeB)/(D/Real(5)));
                            w = half * (w_inf + w_ic) + half * (w_ic - w_inf)*erf((x[1]-spongeB)/(D/Real(5)));
                            p = half * (p_inf + p_ic) + half * (p_ic - p_inf)*erf((x[1]-spongeB)/(D/Real(5)));
                        }
                        else
                        {
                            u = half * (u_inf + u_ic) - half * (u_ic - u_inf)*erf((x[1]-spongeT)/(D/Real(5)));
                            v = half * (v_inf + v_ic) - half * (v_ic - v_inf)*erf((x[1]-spongeT)/(D/Real(5)));
                            w = half * (w_inf + w_ic) - half * (w_ic - w_inf)*erf((x[1]-spongeT)/(D/Real(5)));
                            p = half * (p_inf + p_ic) - half * (p_ic - p_inf)*erf((x[1]-spongeT)/(D/Real(5)));

                        }
                        
                        if (x[2] < z_c) 
                        {
                            u = half * (u_inf + u_ic) + half * (u_ic - u_inf)*erf((x[2]-spongeK)/(D/Real(5)));
                            v = half * (v_inf + v_ic) + half * (v_ic - v_inf)*erf((x[2]-spongeK)/(D/Real(5)));
                            w = half * (w_inf + w_ic) + half * (w_ic - w_inf)*erf((x[2]-spongeK)/(D/Real(5)));
                            p = half * (p_inf + p_ic) + half * (p_ic - p_inf)*erf((x[2]-spongeK)/(D/Real(5)));
                        }
                        else

                        {
                            u = half * (u_inf + u_ic) - half * (u_ic - u_inf)*erf((x[2]-spongeF)/(D/Real(5)));
                            v = half * (v_inf + v_ic) - half * (v_ic - v_inf)*erf((x[2]-spongeF)/(D/Real(5)));
                            w = half * (w_inf + w_ic) - half * (w_ic - w_inf)*erf((x[2]-spongeF)/(D/Real(5)));
                            p = half * (p_inf + p_ic) - half * (p_ic - p_inf)*erf((x[2]-spongeF)/(D/Real(5)));
                        }

                        rho_u[idx_cell] = rho_inf*u;
                        rho_v[idx_cell] = rho_inf*v;
                        rho_w[idx_cell] = rho_inf*w;

                        E[idx_cell]     = (p / (gamma - Real(1)) + half * rho_inf * (u * u + v * v + w * w));
                    }
                } 
            }
        }
        else if (d_flow_model_type == FLOW_MODEL::FIVE_EQN_ALLAIRE) // 5 -eqn allaire.
        {
        /* It finds the total density/heat ratio, then loops over every cell to set condiitions and uses the 1-d array / linear index to access the values */
        
            HAMERS_SHARED_PTR<pdat::CellData<Real> > partial_density = conservative_variables[0]; // just data pointers 
            HAMERS_SHARED_PTR<pdat::CellData<Real> > momentum        = conservative_variables[1];
            HAMERS_SHARED_PTR<pdat::CellData<Real> > total_energy    = conservative_variables[2];
            HAMERS_SHARED_PTR<pdat::CellData<Real> > volume_fraction = conservative_variables[3];
            
            Real* Z_rho_1 = partial_density->getPointer(0);
            Real* Z_rho_2 = partial_density->getPointer(1);
            Real* rho_u   = momentum->getPointer(0);
            Real* rho_v   = momentum->getPointer(1);
            Real* E       = total_energy->getPointer(0);
            Real* Z_1     = volume_fraction->getPointer(0);
            Real* Z_2     = volume_fraction->getPointer(1);
            
            // Species 1.
            Real gamma_1 = Real(8)/Real(5); // 1.6
            
            // Species 2.
            Real gamma_2 = Real(7)/Real(5); // 1.4
            
            // Initial conditions.
            Real Z_rho_1_inf = Real(1);
            Real Z_rho_2_inf = Real(1);
            Real u_inf       = Real(1);
            Real v_inf       = Real(1);
            Real p_inf       = Real(1);
            Real Z_1_inf     = Real(1)/Real(2);
            Real Z_2_inf     = Real(1)/Real(2);
            
            if (d_initial_conditions_db != nullptr) // intial condiiton from database
            {
                TBOX_ASSERT(d_initial_conditions_db->keyExists("Z_rho_1_inf"));
                TBOX_ASSERT(d_initial_conditions_db->keyExists("Z_rho_2_inf"));
                TBOX_ASSERT(d_initial_conditions_db->keyExists("u_inf"));
                TBOX_ASSERT(d_initial_conditions_db->keyExists("v_inf"));
                TBOX_ASSERT(d_initial_conditions_db->keyExists("p_inf"));
                TBOX_ASSERT(d_initial_conditions_db->keyExists("Z_1_inf"));
                TBOX_ASSERT(d_initial_conditions_db->keyExists("Z_2_inf"));
                
                Z_rho_1_inf = d_initial_conditions_db->getReal("Z_rho_1_inf");
                Z_rho_2_inf = d_initial_conditions_db->getReal("Z_rho_2_inf");
                u_inf       = d_initial_conditions_db->getReal("u_inf");
                v_inf       = d_initial_conditions_db->getReal("v_inf");
                p_inf       = d_initial_conditions_db->getReal("p_inf");
                Z_1_inf     = d_initial_conditions_db->getReal("Z_1_inf");
                Z_2_inf     = d_initial_conditions_db->getReal("Z_2_inf");
            }
            
            const Real rho_inf = Z_rho_1_inf + Z_rho_2_inf; // inital density total 
            const Real gamma_m = Real(1)/(Z_1_inf/(gamma_1 - Real(1)) + Z_2_inf/(gamma_2 - Real(1))) + Real(1); // looks like heat ratio 
            
            for (int j = 0; j < patch_dims[1]; j++) // for every single cell 
            {
                for (int i = 0; i < patch_dims[0]; i++)
                {
                    // Compute index into linear data array.
                    int idx_cell = i + j*patch_dims[0];
                    
                    Z_rho_1[idx_cell] = Z_rho_1_inf;
                    Z_rho_2[idx_cell] = Z_rho_2_inf;
                    rho_u[idx_cell]   = rho_inf*u_inf; 
                    rho_v[idx_cell]   = rho_inf*v_inf;
                    E[idx_cell]       = p_inf/(gamma_m - Real(1)) + Real(1)/Real(2)*rho_inf*(u_inf*u_inf + v_inf*v_inf); // total E
                    Z_1[idx_cell]     = Z_1_inf;
                    Z_2[idx_cell]     = Z_2_inf;
                }
            }
        }
    }
}
